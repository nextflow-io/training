# भाग 1: बुनियादी संचालन चलाएं

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow Run प्रशिक्षण कोर्स के इस पहले भाग में, हम एक बहुत ही बुनियादी डोमेन-अज्ञेयवादी Hello World उदाहरण के साथ विषय में आसानी से प्रवेश करते हैं, जिसका उपयोग हम आवश्यक संचालन प्रदर्शित करने और संबंधित Nextflow कोड कॉम्पोनेंट्स को इंगित करने के लिए करेंगे।

??? info "Hello World उदाहरण क्या है?"

    "Hello World!" एक न्यूनतम उदाहरण है जिसका उद्देश्य प्रोग्रामिंग भाषा या software framework के बुनियादी syntax और संरचना को प्रदर्शित करना है।
    उदाहरण में आमतौर पर "Hello, World!" वाक्यांश को आउटपुट डिवाइस, जैसे console या terminal, पर प्रिंट करना, या इसे एक फ़ाइल में लिखना शामिल होता है।

---

## 1. Hello World सीधे चलाएं

आइए इस अवधारणा को एक सरल कमांड के साथ प्रदर्शित करें जिसे हम सीधे terminal में चलाते हैं, यह दिखाने के लिए कि यह क्या करता है इससे पहले कि हम इसे Nextflow में wrap करें।

!!! tip "सुझाव"

    याद रखो कि तुम अब `nextflow-run/` डायरेक्टरी के अंदर होने चाहिए जैसा कि [शुरू करना](00_orientation.md) पृष्ठ पर वर्णित है।

### 1.1. Terminal को hello कहलवाओ

अपने terminal में निम्नलिखित कमांड चलाओ।

```bash
echo 'Hello World!'
```

??? success "कमांड आउटपुट"

    ```console
    Hello World!
    ```

यह terminal में सीधे 'Hello World' टेक्स्ट आउटपुट करता है।

### 1.2. आउटपुट को एक फ़ाइल में लिखो

Pipelines चलाने में ज्यादातर फ़ाइलों से डेटा पढ़ना और परिणामों को अन्य फ़ाइलों में लिखना शामिल है, तो आइए उदाहरण को थोड़ा और प्रासंगिक बनाने के लिए कमांड को टेक्स्ट आउटपुट को एक फ़ाइल में लिखने के लिए संशोधित करें।

```bash
echo 'Hello World!' > output.txt
```

??? success "कमांड आउटपुट"

    ```console

    ```

यह terminal में कुछ भी आउटपुट नहीं करता है।

### 1.3. आउटपुट खोजो

'Hello World' टेक्स्ट अब हमारे द्वारा निर्दिष्ट आउटपुट फ़ाइल में होना चाहिए, जिसका नाम `output.txt` है।
तुम इसे file explorer में खोल सकते हो या कमांड लाइन से `cat` utility का उपयोग करके, उदाहरण के लिए।

??? abstract "फ़ाइल सामग्री"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

यह वही है जो हम अपने पहले Nextflow workflow के साथ replicate करने की कोशिश करने जा रहे हैं।

### सीख

अब तुम जानते हो कि terminal में एक सरल कमांड कैसे चलाया जाए जो कुछ टेक्स्ट आउटपुट करता है, और वैकल्पिक रूप से, इसे एक फ़ाइल में आउटपुट कैसे लिखवाया जाए।

### आगे क्या?

जानो कि एक Nextflow workflow चलाने में क्या लगता है जो समान परिणाम प्राप्त करता है।

---

## 2. Workflow चलाओ

हम तुम्हें `1-hello.nf` नाम की एक workflow script प्रदान करते हैं जो `--input` नाम के कमांड-लाइन आर्गुमेंट के माध्यम से एक इनपुट greeting लेती है और उस greeting को शामिल करते हुए एक टेक्स्ट फ़ाइल उत्पन्न करती है।

हम अभी कोड नहीं देखने जा रहे हैं; पहले आइए देखें कि इसे चलाना कैसा दिखता है।

### 2.1. Workflow लॉन्च करो और execution मॉनिटर करो

Terminal में, निम्नलिखित कमांड चलाओ:

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [a3/7be2fa] sayHello | 1 of 1 ✔
    ```

यदि तुम्हारा console आउटपुट कुछ इस तरह दिखता है, तो बधाई हो, तुमने अभी-अभी अपना पहला Nextflow workflow चलाया है!

यहां सबसे महत्वपूर्ण आउटपुट अंतिम पंक्ति है, जो ऊपर के आउटपुट में हाइलाइट किया गया है:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

यह हमें बताता है कि `sayHello` process एक बार सफलतापूर्वक execute हुआ (`1 of 1 ✔`)।

यह बढ़िया है, लेकिन तुम सोच रहे होगे: आउटपुट कहां है?

### 2.2. `results` डायरेक्टरी में आउटपुट फ़ाइल खोजो

यह workflow अपने आउटपुट को results डायरेक्टरी में publish करने के लिए कॉन्फ़िगर की गई है।
यदि तुम अपनी current डायरेक्टरी देखो, तुम देखोगे कि जब तुमने workflow चलाया, Nextflow ने `results` नाम की एक नई डायरेक्टरी बनाई, साथ ही उसके अंदर `1-hello` नाम की एक subdirectory, जिसमें `output.txt` नाम की एक फ़ाइल है।

```console title="results/"
results
└── 1-hello
    └── output.txt
```

फ़ाइल खोलो; सामग्री उस string से मेल खानी चाहिए जो तुमने कमांड लाइन पर निर्दिष्ट की थी।

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

बढ़िया, हमारे workflow ने वह किया जो उसे करना चाहिए था!

### 2.3. परिणामों को एक अलग डायरेक्टरी में save करो

डिफ़ॉल्ट रूप से, Nextflow pipeline outputs को तुम्हारे current path में `results` नाम की डायरेक्टरी में save करेगा।
यह बदलने के लिए कि तुम्हारी फ़ाइलें कहां publish होती हैं, `-output-dir` CLI flag का उपयोग करो (या संक्षेप में `-o`)

!!! danger

    ध्यान दो कि `--input` में दो hyphens हैं और `-output-dir` में एक!
    यह इसलिए है क्योंकि `--input` एक pipeline _parameter_ है और `-output-dir` एक core Nextflow CLI flag है।
    इन पर बाद में और जानकारी।

```bash
nextflow run 1-hello.nf --input 'Hello World!' -output-dir hello_results
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [hungry_celsius] DSL2 - revision: f048d6ea78

    executor >  local (1)
    [a3/1e1535] sayHello [100%] 1 of 1 ✔
    ```

तुम्हें देखना चाहिए कि तुम्हारे outputs अब `results` के बजाय `hello_results` नाम की डायरेक्टरी में publish हो रहे हैं:

```console title="hello_results/"
hello_results
└── 1-hello
    └── output.txt
```

इस डायरेक्टरी के अंदर की फ़ाइलें पहले जैसी ही हैं, बस top-level डायरेक्टरी अलग है।
हालांकि, दोनों मामलों में ध्यान रखो कि 'published' परिणाम वास्तव में उस आउटपुट की एक प्रति (या कुछ मामलों में एक symbolic link) है जो Nextflow ने workflow execute करते समय उत्पन्न किया था।

तो अब, हम hood के नीचे झांकने जा रहे हैं यह देखने के लिए कि Nextflow ने वास्तव में कहां काम execute किया।

!!! warning "चेतावनी"

    सभी workflows results डायरेक्टरी में outputs publish करने के लिए सेट अप नहीं होंगी, और/या डायरेक्टरी के नाम और संरचना भिन्न हो सकती है।
    इस सेक्शन में थोड़ा आगे, हम तुम्हें दिखाएंगे कि यह व्यवहार कहां निर्दिष्ट है यह कैसे पता करें।

### 2.4. `work/` डायरेक्टरी में मूल आउटपुट और लॉग खोजो

जब तुम एक workflow चलाते हो, Nextflow workflow में प्रत्येक process के हर एक invocation के लिए एक अलग 'task directory' बनाता है (=pipeline में प्रत्येक step)।
प्रत्येक के लिए, यह आवश्यक inputs को stage करेगा, प्रासंगिक instruction(s) execute करेगा और outputs और log फ़ाइलें उस एक डायरेक्टरी के भीतर लिखेगा, जिसे इसे unique बनाने के लिए hash का उपयोग करके स्वचालित रूप से नाम दिया जाता है।

ये सभी task directories तुम्हारी current डायरेक्टरी (जहां से तुम कमांड चला रहे हो) के भीतर `work` नाम की एक डायरेक्टरी के अंतर्गत रहेंगी।

यह भ्रमित करने वाला लग सकता है, तो आइए देखें कि व्यवहार में यह कैसा दिखता है।

पहले चलाए गए workflow के console आउटपुट पर वापस जाते हुए, हमारे पास यह पंक्ति थी:

```console
[a3/1e1535] sayHello [100%] 1 of 1 ✔
```

देखो कि पंक्ति `[a3/1e1535]` से कैसे शुरू होती है?
यह उस एक process call के लिए task directory path का एक छोटा रूप है, और तुम्हें बताता है कि `sayHello` process call का आउटपुट `work/` directory path के भीतर कहां खोजना है।

तुम निम्नलिखित कमांड टाइप करके पूर्ण path पा सकते हो (`a3/1e1535` को अपने terminal में जो दिखता है उससे बदलो) और path को autocomplete करने के लिए tab key दबाओ या asterisk जोड़ो:

```bash
ls work/a3/1e1535*
```

यह पूर्ण path directory path देना चाहिए: `work/a3/1e153543b0a7f9d2c4735ddb4ab231`

आइए देखें कि वहां क्या है।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    work
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "वही चीज़ नहीं दिख रही?"

    सटीक subdirectory नाम तुम्हारे system पर अलग होंगे।

    यदि तुम VSCode file explorer में task subdirectory की सामग्री browse करो, तुम्हें सभी फ़ाइलें तुरंत दिखेंगी।
    हालांकि, log फ़ाइलें terminal में invisible होने के लिए सेट हैं, इसलिए यदि तुम उन्हें देखने के लिए `ls` या `tree` का उपयोग करना चाहते हो, तुम्हें invisible फ़ाइलें प्रदर्शित करने के लिए प्रासंगिक विकल्प सेट करना होगा।

    ```bash
    tree -a work
    ```

`work/` में दो सेट directories हैं, जो हमारे द्वारा किए गए दो अलग-अलग pipeline runs से हैं।
प्रत्येक task execution को काम करने के लिए अपनी खुद की, isolated, डायरेक्टरी मिलती है।
इस मामले में pipeline ने दोनों बार एक ही काम किया, इसलिए प्रत्येक task डायरेक्टरी की सामग्री identical है।

तुम्हें तुरंत `output.txt` फ़ाइल पहचानना चाहिए, जो वास्तव में `sayHello` process का मूल आउटपुट है जो `results` डायरेक्टरी में publish हुआ था।
यदि तुम इसे खोलो, तुम्हें फिर से `Hello World!` greeting मिलेगी।

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/output.txt"
Hello World!
```

तो उन सभी अन्य फ़ाइलों के बारे में क्या?

ये helper और log फ़ाइलें हैं जो Nextflow ने task execution के भाग के रूप में लिखीं:

- **`.command.begin`**: Sentinel फ़ाइल जो task लॉन्च होते ही बनाई जाती है।
- **`.command.err`**: Process call द्वारा emit किए गए Error संदेश (`stderr`)
- **`.command.log`**: Process call द्वारा emit किया गया पूर्ण log आउटपुट
- **`.command.out`**: Process call द्वारा Regular आउटपुट (`stdout`)
- **`.command.run`**: Process call को execute करने के लिए Nextflow द्वारा चलाई गई पूर्ण script
- **`.command.sh`**: वह कमांड जो वास्तव में process call द्वारा चलाया गया
- **`.exitcode`**: कमांड से resulting exit code

`.command.sh` फ़ाइल विशेष रूप से उपयोगी है क्योंकि यह तुम्हें मुख्य कमांड दिखाती है जो Nextflow ने execute किया, सभी bookkeeping और task/environment setup को शामिल किए बिना।

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/.command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

तो यह पुष्टि करता है कि workflow ने वही कमांड compose किया जो हमने पहले सीधे कमांड-लाइन पर चलाया था।

जब कुछ गलत हो जाता है और तुम्हें troubleshoot करने की आवश्यकता होती है कि क्या हुआ, `command.sh` script को देखना उपयोगी हो सकता है ताकि यह जांचा जा सके कि Nextflow ने workflow instructions, variable interpolation आदि के आधार पर कौन सा कमांड compose किया।

### 2.5. अलग-अलग greetings के साथ workflow फिर से चलाओ

`--input` आर्गुमेंट के लिए अलग-अलग values के साथ workflow को कुछ बार फिर से चलाने की कोशिश करो, फिर task directories देखो।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    work/
    ├── 09
    │   └── 5ea8665939daf6f04724286c9b3c8a
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 92
    │   └── ceb95e05d87621c92a399da9bd2067
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 93
    │   └── 6708dbc20c7efdc6769cbe477061ec
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

तुम देखते हो कि प्रत्येक run के लिए output और log फ़ाइलों के पूर्ण सेट के साथ एक नई subdirectory बनाई गई है।

इसके विपरीत, यदि तुम `results` डायरेक्टरी देखो, वहां अभी भी केवल एक सेट परिणाम है, और output फ़ाइल की सामग्री जो तुमने आखिरी में चलाया उससे मेल खाती है।

??? abstract "डायरेक्टरी सामग्री"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

यह तुम्हें दिखाता है कि published results बाद के executions द्वारा overwrite हो जाएंगे, जबकि `work/` के अंतर्गत task directories संरक्षित रहती हैं।

### सीख

तुम जानते हो कि एक सरल Nextflow script कैसे चलाएं, इसके execution को कैसे monitor करें और इसके outputs कैसे खोजें।

### आगे क्या?

सीखो कि एक बुनियादी Nextflow script कैसे पढ़ें और पहचानें कि इसके components इसकी functionality से कैसे संबंधित हैं।

---

## 3. Hello World workflow starter script की जांच करो

जो हमने वहां किया वह मूल रूप से workflow script को एक black box की तरह treat करना था।
अब जब हमने देख लिया कि यह क्या करता है, आइए box खोलें और अंदर देखें।

यहां हमारा लक्ष्य Nextflow code का syntax याद करना नहीं है, बल्कि कुछ बुनियादी अंतर्ज्ञान बनाना है कि मुख्य components क्या हैं और वे कैसे organized हैं।

### 3.1. समग्र code संरचना की जांच करो

तुम्हें `1-hello.nf` script अपनी current डायरेक्टरी में मिलेगी, जो `nextflow-run` होनी चाहिए। इसे editor pane में खोलो।

??? full-code "पूर्ण code फ़ाइल"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

    /*
    * Pipeline parameters
    */
    params {
        input: String
    }

    workflow {

        main:
        // एक अभिवादन emit करें
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

एक Nextflow workflow script में आमतौर पर एक या अधिक **process** definitions, **workflow** स्वयं, और कुछ वैकल्पिक blocks जैसे **params** और **output** शामिल होते हैं।

प्रत्येक **process** वर्णन करता है कि pipeline में संबंधित step को क्या operation(s) पूरा करना चाहिए, जबकि **workflow** dataflow logic का वर्णन करता है जो विभिन्न steps को जोड़ता है।

आइए पहले **process** block पर करीब से नज़र डालें, फिर हम **workflow** block देखेंगे।

### 3.2. `process` definition

Code का पहला block एक [**process**](https://nextflow.io/docs/latest/process.html) का वर्णन करता है।
Process definition `process` keyword से शुरू होती है, उसके बाद process का नाम और अंत में curly braces द्वारा delimit किया गया process body।
Process body में एक script block होना चाहिए जो चलाने के लिए कमांड निर्दिष्ट करता है, जो कुछ भी हो सकता है जो तुम कमांड लाइन terminal में चला सको।

```groovy title="1-hello.nf" linenums="3"
/*
* Use echo to print a greeting to a file
*/
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

यहां हमारे पास `sayHello` नाम का एक **process** है जो `greeting` नाम का एक **input** variable लेता है और अपना **output** `output.txt` नाम की फ़ाइल में लिखता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

यह एक बहुत ही न्यूनतम process definition है जिसमें बस एक `input` definition, एक `output` definition और execute करने के लिए `script` है।

`input` definition में `val` qualifier शामिल है, जो Nextflow को बताता है कि किसी प्रकार का value expect करें (string, number, जो भी हो सकता है)।

`output` definition में `path` qualifier शामिल है, जो Nextflow को बताता है कि इसे path के रूप में handle किया जाना चाहिए (इसमें directory paths और files दोनों शामिल हैं)।

### 3.3. `workflow` definition

Code का दूसरा block [**workflow**](https://nextflow.io/docs/latest/workflow.html) स्वयं का वर्णन करता है।
Workflow definition `workflow` keyword से शुरू होती है, उसके बाद एक वैकल्पिक नाम, फिर curly braces द्वारा delimit किया गया workflow body।

यहां हमारे पास एक **workflow** है जिसमें एक `main:` block और एक `publish:` block है।
`main:` block workflow का मुख्य body है और `publish:` block उन outputs को सूचीबद्ध करता है जो `results` डायरेक्टरी में publish होने चाहिए।

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // एक अभिवादन emit करें
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

इस मामले में `main:` block में `sayHello` process का एक call है और इसे greeting के रूप में उपयोग करने के लिए `params.input` नाम का एक input देता है।

जैसा कि हम थोड़ी देर में और विस्तार से चर्चा करेंगे, `params.input` वह value रखता है जो हमने अपनी कमांड लाइन में `--input` parameter को दी थी।

`publish:` block `sayHello()` process call के output को सूचीबद्ध करता है, जिसे यह `sayHello.out` के रूप में संदर्भित करता है और `first_output` नाम देता है (यह कुछ भी हो सकता है जो workflow author चाहे)।

यह एक बहुत ही न्यूनतम **workflow** definition है।
एक real-world pipeline में, workflow में आमतौर पर **channels** द्वारा जुड़े **processes** के कई calls होते हैं, और variable inputs के लिए default values सेट अप हो सकते हैं।

हम कोर्स के भाग 2 में इसमें जाएंगे।
अभी के लिए, आइए इस पर करीब से नज़र डालें कि हमारा workflow inputs और outputs को कैसे handle कर रहा है।

### 3.4. कमांड-लाइन parameters की `params` system

`params.input` जो हम `sayHello()` process call को प्रदान करते हैं वह Nextflow code का एक अच्छा टुकड़ा है और इस पर एक अतिरिक्त मिनट खर्च करने लायक है।

जैसा कि ऊपर बताया गया है, इस तरह हम `--input` कमांड-लाइन parameter की value को `sayHello()` process call को pass करते हैं।
वास्तव में, बस `params.someParameterName` declare करना workflow को कमांड-लाइन से `--someParameterName` नाम का parameter देने के लिए पर्याप्त है।

यहां हमने उस parameter declaration को एक `params` block सेट अप करके formalize किया है जो उस प्रकार के input को निर्दिष्ट करता है जिसकी workflow expect करती है (Nextflow 25.10.2 और बाद में)।

```groovy title="1-hello.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String
}
```

समर्थित types में `String`, `Integer`, `Float`, `Boolean`, और `Path` शामिल हैं।
अधिक जानने के लिए, Nextflow reference documentation में [Workflow parameters](https://nextflow.io/docs/latest/config.html#workflow-parameters) देखें।

!!! tip "सुझाव"

    याद रखो कि `params` system का उपयोग करके declare किए गए _Workflow_ parameters कमांड लाइन पर हमेशा दो dashes (`--`) लेते हैं।
    यह उन्हें _Nextflow-level_ CLI flags से अलग करता है, जो केवल एक dash (`-`) लेते हैं।

### 3.5. `publish` directive

Workflow के दूसरे छोर पर, हमने पहले ही `publish:` block पर एक नज़र डाली है।
यह output handling system का एक आधा है; दूसरा आधा नीचे स्थित `output` block है।

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

यह निर्दिष्ट करता है कि `publish:` block में सूचीबद्ध `first_output` output को default `results` output डायरेक्टरी के अंतर्गत `1-hello` नाम की subdirectory में copy किया जाना चाहिए।

`mode 'copy'` पंक्ति system के default व्यवहार को override करती है, जो proper copy के बजाय `work/` डायरेक्टरी में मूल फ़ाइल का symbolic link (या symlink) बनाना है।

publishing व्यवहार को control करने के लिए यहां प्रदर्शित की तुलना में अधिक विकल्प हैं; हम बाद में कुछ cover करेंगे।
तुम यह भी देखोगे कि जब एक workflow कई outputs generate करता है, तो प्रत्येक को `output` block में इस तरह सूचीबद्ध किया जाता है।

अधिक जानने के लिए, Nextflow reference documentation में [Publishing outputs](https://nextflow.io/docs/latest/workflow.html#publishing-outputs) देखें।

??? info "`publishDir` का उपयोग करके outputs publish करने का पुराना syntax"

    हाल ही तक, outputs publish करने का स्थापित तरीका प्रत्येक individual process के स्तर पर `publishDir` directive का उपयोग करके करना था।

    तुम अभी भी पुरानी Nextflow pipelines और process modules में हर जगह यह code pattern पाओगे, इसलिए इसके बारे में जागरूक होना महत्वपूर्ण है।

    Workflow में `publish:` block और top level पर `output` block होने के बजाय, तुम `sayHello` process definition में एक `publishDir` पंक्ति देखोगे:

    ```groovy title="Syntax example" linenums="1" hl_lines="3"
    process sayHello {

        publishDir 'results/1-hello', mode: 'copy'

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }
    ```

    हालांकि, हम किसी भी नए काम में इसका उपयोग करने की अनुशंसा नहीं करते क्योंकि यह अंततः Nextflow भाषा के भविष्य के versions में disallowed हो जाएगा।

### सीख

अब तुम जानते हो कि एक सरल Nextflow workflow कैसे structured है, और बुनियादी components इसकी functionality से कैसे संबंधित हैं।

### आगे क्या?

सीखो कि अपने workflow executions को सुविधाजनक तरीके से कैसे manage करें।

---

## 4. Workflow executions manage करो

Workflows लॉन्च करना और outputs प्राप्त करना जानना बढ़िया है, लेकिन तुम जल्दी पाओगे कि workflow management के कुछ अन्य पहलू हैं जो तुम्हारी ज़िंदगी आसान बना देंगे।

यहां हम तुम्हें दिखाते हैं कि `resume` feature का लाभ कैसे उठाएं जब तुम्हें वही workflow फिर से लॉन्च करना हो, `nextflow log` के साथ execution logs का निरीक्षण कैसे करें, और `nextflow clean` के साथ पुरानी work directories कैसे delete करें।

### 4.1. `-resume` के साथ workflow फिर से लॉन्च करो

कभी-कभी, तुम एक pipeline फिर से चलाना चाहोगे जिसे तुमने पहले लॉन्च किया था बिना उस काम को दोहराए जो पहले ही सफलतापूर्वक पूरा हो चुका है।

Nextflow में `-resume` नाम का एक विकल्प है जो तुम्हें ऐसा करने की अनुमति देता है।
विशेष रूप से, इस mode में, कोई भी processes जो पहले से ही exact same code, settings और inputs के साथ run हो चुके हैं, skip कर दिए जाएंगे।
इसका मतलब है कि Nextflow केवल वे processes run करेगा जो तुमने पिछले run के बाद से add या modify किए हैं, या जिन्हें तुम नई settings या inputs प्रदान कर रहे हो।

ऐसा करने के दो मुख्य फायदे हैं:

- यदि तुम एक pipeline develop करने के बीच में हो, तुम अधिक तेज़ी से iterate कर सकते हो क्योंकि तुम्हें अपने changes test करने के लिए केवल वह process(es) run करने होंगे जिन पर तुम सक्रिय रूप से काम कर रहे हो।
- यदि तुम production में एक pipeline चला रहे हो और कुछ गलत हो जाता है, कई मामलों में तुम समस्या ठीक कर सकते हो और pipeline फिर से लॉन्च कर सकते हो, और यह failure के point से running resume करेगी, जो तुम्हारा बहुत समय और compute बचा सकती है।

इसका उपयोग करने के लिए, बस अपने कमांड में `-resume` जोड़ो और इसे run करो:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "कमांड आउटपुट"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

Console आउटपुट परिचित दिखना चाहिए, लेकिन पहले की तुलना में एक चीज़ थोड़ी अलग है।

`cached:` bit देखो जो process status पंक्ति (पंक्ति 5) में जोड़ी गई है, जिसका मतलब है कि Nextflow ने पहचान लिया है कि यह काम पहले ही कर चुका है और बस पिछले सफल run से परिणाम का पुन: उपयोग कर रहा है।

तुम यह भी देख सकते हो कि work subdirectory hash पिछले run जैसा ही है।
Nextflow literally तुम्हें पिछले execution की ओर इंगित कर रहा है और कह रहा है "मैंने वह वहां पहले ही कर दिया था।"

!!! tip "सुझाव"

    जब तुम `resume` के साथ pipeline फिर से चलाते हो, Nextflow पहले सफलतापूर्वक run हुए किसी भी executions द्वारा work directory के बाहर publish की गई किसी भी फ़ाइल को overwrite नहीं करता।

    अधिक जानने के लिए, Nextflow reference documentation में [Cache and resume](https://nextflow.io/docs/latest/cache-and-resume.html) देखें।

### 4.2. पिछले executions के log का निरीक्षण करो

जब भी तुम एक nextflow workflow लॉन्च करते हो, `history` नाम की एक log फ़ाइल में एक पंक्ति लिखी जाती है, जो current working directory में `.nextflow` नाम की एक hidden डायरेक्टरी के अंतर्गत होती है।

??? abstract "फ़ाइल सामग्री"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

यह फ़ाइल तुम्हें current working directory के भीतर से लॉन्च किए गए हर Nextflow run के लिए timestamp, run name, status, revision ID, session ID और पूर्ण कमांड लाइन देती है।

इस जानकारी तक पहुंचने का अधिक सुविधाजनक तरीका [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) कमांड का उपयोग करना है।

```bash
nextflow log
```

??? success "कमांड आउटपुट"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

यह log फ़ाइल की सामग्री को terminal में आउटपुट करेगा, एक header पंक्ति के साथ augmented।

तुम देखोगे कि session ID जब भी तुम एक नया `nextflow run` कमांड चलाते हो तब बदलता है, सिवाय इसके कि यदि तुम `-resume` विकल्प का उपयोग कर रहे हो।
उस मामले में, session ID वही रहता है।

Nextflow session ID का उपयोग `cache` डायरेक्टरी के अंतर्गत run caching जानकारी group करने के लिए करता है, जो `.nextflow` के अंतर्गत भी स्थित है।

### 4.3. पुरानी work directories delete करो

यदि तुम बहुत सारी pipelines चलाते हो, तुम कई subdirectories में बहुत सारी फ़ाइलें accumulate कर सकते हो।
चूंकि subdirectories randomly नाम दी जाती हैं, उनके नामों से यह बताना मुश्किल है कि कौन से पुराने हैं बनाम अधिक हाल के runs।

सौभाग्य से Nextflow में एक helpful कमांड शामिल है जिसे [`nextflow clean`](https://www.nextflow.io/docs/latest/reference/cli.html#clean) कहते हैं जो स्वचालित रूप से पिछले runs के लिए work subdirectories delete कर सकता है जिनकी तुम्हें अब परवाह नहीं है।

#### 4.3.1. Deletion criteria निर्धारित करो

यह निर्धारित करने के लिए कि क्या delete करना है, कई विकल्प हैं, जिन्हें तुम ऊपर linked documentation में explore कर सकते हो।
यहां हम तुम्हें एक उदाहरण दिखाते हैं जो एक दिए गए run से पहले के सभी runs की subdirectories delete करता है, जो इसके run name का उपयोग करके निर्दिष्ट किया गया है।

सबसे हाल का सफल run देखो जहां तुमने `-resume` का उपयोग नहीं किया; हमारे मामले में run name `backstabbing_swartz` था।

Run name machine-generated two-part string है जो `Launching (...)` console आउटपुट पंक्ति में square brackets में दिखाया जाता है।
तुम run को उसके timestamp और/या कमांड लाइन के आधार पर देखने के लिए Nextflow log का भी उपयोग कर सकते हो।

#### 4.3.2. Dry run करो

पहले हम dry run flag `-n` का उपयोग यह जांचने के लिए करते हैं कि कमांड को देखते हुए क्या delete होगा:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "कमांड आउटपुट"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

तुम्हारे आउटपुट में अलग task directory नाम होंगे और पंक्तियों की संख्या अलग हो सकती है, लेकिन यह उदाहरण के समान दिखना चाहिए।

यदि तुम कोई पंक्तियां आउटपुट नहीं देखते हो, तो या तो तुमने एक valid run name प्रदान नहीं किया या delete करने के लिए कोई पिछले runs नहीं हैं। उदाहरण कमांड में `backstabbing_swartz` को अपने log में जो भी संबंधित latest run name है उसमें बदलना सुनिश्चित करो।

#### 4.3.3. Deletion के साथ आगे बढ़ो

यदि आउटपुट अपेक्षित दिखता है और तुम deletion के साथ आगे बढ़ना चाहते हो, `-n` के बजाय `-f` flag के साथ कमांड फिर से चलाओ:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "कमांड आउटपुट"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

आउटपुट पहले जैसा ही होना चाहिए, लेकिन अब 'Would remove' के बजाय 'Removed' कह रहा है।
ध्यान दो कि यह two-character subdirectories (जैसे ऊपर `eb/`) को नहीं हटाता लेकिन यह उनकी सामग्री खाली कर देता है।

!!! warning "चेतावनी"

    पिछले runs की work subdirectories delete करने से उन्हें Nextflow के cache से हटा दिया जाता है और उन directories में stored कोई भी outputs delete हो जाते हैं।
    इसका मतलब है कि यह संबंधित processes को फिर से चलाए बिना execution resume करने की Nextflow की क्षमता को तोड़ देता है।

    तुम किसी भी outputs को save करने के लिए जिम्मेदार हो जिनकी तुम्हें परवाह है! यही मुख्य कारण है कि हम `publish` directive के लिए `symlink` mode के बजाय `copy` mode का उपयोग करना पसंद करते हैं।

### सीख

तुम जानते हो कि पहले से identical तरीके से run हुए steps को दोहराए बिना pipeline कैसे फिर से लॉन्च करें, execution log का निरीक्षण कैसे करें, और पुरानी work directories को साफ करने के लिए `nextflow clean` कमांड का उपयोग कैसे करें।

### आगे क्या?

थोड़ा ब्रेक लो! तुमने अभी-अभी Nextflow syntax और बुनियादी उपयोग निर्देशों के building blocks को absorb किया है।

इस प्रशिक्षण के अगले section में, हम Hello World pipeline के चार क्रमिक रूप से अधिक realistic versions देखने जा रहे हैं जो प्रदर्शित करेंगे कि Nextflow तुम्हें कई inputs को कुशलता से process करने, एक साथ जुड़े कई steps से बनी workflows चलाने, modular code components का लाभ उठाने, और अधिक reproducibility और portability के लिए containers का उपयोग करने की अनुमति कैसे देता है।

---

## Quiz

<quiz>
Console आउटपुट पंक्ति `[a3/7be2fa] SAYHELLO | 1 of 1 ✔` में, `[a3/7be2fa]` क्या दर्शाता है?
- [ ] Process version number
- [ ] एक unique run identifier
- [x] Task की work directory का truncated path
- [ ] Output फ़ाइल का checksum

और जानें: [2.4. `work/` डायरेक्टरी में मूल आउटपुट और लॉग खोजो](#24-work-डायरेक्टरी-में-मूल-आउटपुट-और-लॉग-खोजो)
</quiz>

<quiz>
Task directory में `.command.sh` फ़ाइल का उद्देश्य क्या है?
- [ ] यह task की configuration settings store करती है
- [x] यह process द्वारा execute किया गया actual कमांड दिखाती है
- [ ] इसमें failed tasks से error messages होते हैं
- [ ] यह task के लिए staged input फ़ाइलों को सूचीबद्ध करती है

और जानें: [2.4. `work/` डायरेक्टरी में मूल आउटपुट और लॉग खोजो](#24-work-डायरेक्टरी-में-मूल-आउटपुट-और-लॉग-खोजो)
</quiz>

<quiz>
जब तुम `-resume` के बिना workflow फिर से चलाते हो तो published results का क्या होता है?
- [ ] वे अलग timestamped directories में preserve होते हैं
- [x] वे नए execution द्वारा overwrite हो जाते हैं
- [ ] Nextflow overwriting रोकता है और fail होता है
- [ ] वे automatically backup हो जाते हैं

और जानें: [2.5. अलग-अलग greetings के साथ workflow फिर से चलाओ](#25-अलग-अलग-greetings-के-साथ-workflow-फिर-से-चलाओ)
</quiz>

<quiz>
यह console आउटपुट क्या indicate करता है?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] Task fail हुआ और skip किया गया
- [ ] Task queue में wait कर रहा है
- [x] Nextflow ने पिछले identical execution से results का पुन: उपयोग किया
- [ ] Task manually cancel किया गया

और जानें: [4.1. `-resume` के साथ workflow फिर से लॉन्च करो](#41--resume-के-साथ-workflow-फिर-से-लॉन्च-करो)
</quiz>

<quiz>
Nextflow `nextflow log` कमांड जो execution history प्रदर्शित करता है उसे कहां store करता है?
- [ ] results डायरेक्टरी में
- [ ] work डायरेक्टरी में
- [x] `.nextflow/history` फ़ाइल में
- [ ] `nextflow.config` में

और जानें: [4.2. पिछले executions के log का निरीक्षण करो](#42-पिछले-executions-के-log-का-निरीक्षण-करो)
</quiz>

<quiz>
Workflow फ़ाइल में `params` block का उद्देश्य क्या है?
- [ ] Process resource requirements define करना
- [ ] Executor configure करना
- [x] Workflow input parameters declare और type करना
- [ ] Output publishing options specify करना

और जानें: [3.4. कमांड-लाइन parameters की params system](#34-कमांड-लाइन-parameters-की-params-system)
</quiz>

<quiz>
Workflow के `output` block में, `mode 'copy'` क्या करता है?
- [ ] Work directory का backup बनाता है
- [x] Symbolic links के बजाय फ़ाइलों की पूर्ण copy बनाता है
- [ ] Workflow script को results में copy करता है
- [ ] Incremental file copying enable करता है

और जानें: [3.5. publish directive](#35-publish-directive)
</quiz>

<quiz>
वास्तव में फ़ाइलें delete करने से पहले `nextflow clean` कमांड के साथ कौन सा flag उपयोग करने की अनुशंसा की जाती है?
- [x] `-n` (dry run) यह preview करने के लिए कि क्या delete होगा
- [ ] `-v` (verbose) detailed आउटपुट देखने के लिए
- [ ] `-a` (all) सभी directories select करने के लिए
- [ ] `-q` (quiet) warnings suppress करने के लिए

और जानें: [4.3. पुरानी work directories delete करो](#43-पुरानी-work-directories-delete-करो)
</quiz>
