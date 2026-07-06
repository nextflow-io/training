# वर्कफ़्लो की समस्या निवारण

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

डीबगिंग एक महत्वपूर्ण कौशल है जो तुम्हें घंटों की निराशा से बचा सकती है और तुम्हें एक बेहतर Nextflow डेवलपर बनने में मदद कर सकती है। अपने करियर में, खासकर जब तुम शुरुआत कर रहे हो, वर्कफ़्लो बनाते और बनाए रखते समय तुम्हें बग्स का सामना करना पड़ेगा। व्यवस्थित डीबगिंग तरीके सीखने से तुम्हें समस्याओं को जल्दी पहचानने और हल करने में मदद मिलेगी।

### सीखने के लक्ष्य

इस साइड क्वेस्ट में, हम Nextflow वर्कफ़्लो के लिए **व्यवस्थित डीबगिंग तकनीकों** का पता लगाएंगे:

- **सिंटैक्स एरर डीबगिंग**: IDE फीचर्स और Nextflow एरर मैसेज का प्रभावी उपयोग
- **चैनल डीबगिंग**: डेटा फ्लो समस्याओं और चैनल संरचना की समस्याओं का निदान
- **प्रोसेस डीबगिंग**: एक्जीक्यूशन विफलताओं और रिसोर्स समस्याओं की जांच
- **बिल्ट-इन डीबगिंग टूल्स**: Nextflow के preview mode, stub running, और work directories का उपयोग
- **व्यवस्थित तरीके**: कुशल डीबगिंग के लिए चार-चरण की पद्धति

अंत में, तुम्हारे पास एक मजबूत डीबगिंग पद्धति होगी जो निराशाजनक एरर मैसेज को समाधान के स्पष्ट रोडमैप में बदल देगी।

### पूर्वापेक्षाएं

इस साइड क्वेस्ट को शुरू करने से पहले, तुम्हें:

- [Hello Nextflow](../../hello_nextflow/index.md) ट्यूटोरियल या समकक्ष शुरुआती कोर्स पूरा करना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (प्रोसेस, चैनल, ऑपरेटर) का उपयोग करने में सहज होना चाहिए।

**वैकल्पिक:** हम अनुशंसा करते हैं कि पहले [IDE Features for Nextflow Development](../dev_environment/index.md) साइड क्वेस्ट पूरा करो।
इसमें IDE फीचर्स का व्यापक कवरेज है जो डीबगिंग में सहायता करते हैं (सिंटैक्स हाइलाइटिंग, एरर डिटेक्शन, आदि), जिनका हम यहाँ भरपूर उपयोग करेंगे।

---

## 0. शुरू करना

#### ट्रेनिंग codespace खोलो

अगर तुमने अभी तक नहीं किया है, तो [Environment Setup](../../envsetup/index.md) में बताए अनुसार ट्रेनिंग वातावरण खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाओ

चलो उस डायरेक्टरी में जाते हैं जहाँ इस ट्यूटोरियल की फ़ाइलें हैं।

```bash
cd side-quests/debugging
```

तुम VSCode को इस डायरेक्टरी पर फोकस करने के लिए सेट कर सकते हो:

```bash
code .
```

#### सामग्री की समीक्षा करो

तुम्हें विभिन्न प्रकार के बग्स वाले उदाहरण वर्कफ़्लो का एक सेट मिलेगा जिनका उपयोग हम अभ्यास के लिए करेंगे:

??? abstract "डायरेक्टरी सामग्री"

    ```console
    .
    ├── bad_bash_var.nf
    ├── bad_channel_shape.nf
    ├── bad_channel_shape_viewed_debug.nf
    ├── bad_channel_shape_viewed.nf
    ├── bad_number_inputs.nf
    ├── badpractice_syntax.nf
    ├── bad_resources.nf
    ├── bad_syntax.nf
    ├── buggy_workflow.nf
    ├── data
    │   ├── sample_001.fastq.gz
    │   ├── sample_002.fastq.gz
    │   ├── sample_003.fastq.gz
    │   ├── sample_004.fastq.gz
    │   ├── sample_005.fastq.gz
    │   └── sample_data.csv
    ├── exhausted.nf
    ├── invalid_process.nf
    ├── missing_output.nf
    ├── missing_software.nf
    ├── missing_software_with_stub.nf
    ├── nextflow.config
    └── no_such_var.nf
    ```

ये फ़ाइलें वास्तविक दुनिया के विकास में आने वाले सामान्य डीबगिंग परिदृश्यों को दर्शाती हैं।

#### असाइनमेंट की समीक्षा करो

तुम्हारी चुनौती है कि प्रत्येक वर्कफ़्लो चलाओ, एरर(एरर्स) की पहचान करो, और उन्हें ठीक करो।

प्रत्येक बगी वर्कफ़्लो के लिए:

1. **वर्कफ़्लो चलाओ** और एरर देखो
2. **एरर मैसेज का विश्लेषण करो**: Nextflow तुम्हें क्या बता रहा है?
3. **दिए गए संकेतों का उपयोग करके कोड में समस्या ढूंढो**
4. **बग ठीक करो** और सत्यापित करो कि तुम्हारा समाधान काम करता है
5. **अगले सेक्शन पर जाने से पहले फ़ाइल रीसेट करो** (`git checkout <filename>` का उपयोग करो)

अभ्यास सरल सिंटैक्स एरर से लेकर अधिक सूक्ष्म रनटाइम समस्याओं तक बढ़ते हैं।
समाधान इनलाइन चर्चा किए गए हैं, लेकिन आगे पढ़ने से पहले प्रत्येक को खुद हल करने की कोशिश करो।

#### तैयारी की जांच सूची

क्या तुम शुरू करने के लिए तैयार हो?

- [ ] मैं इस कोर्स के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता/समझती हूँ
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी working directory उचित रूप से सेट की है
- [ ] मैं असाइनमेंट समझता/समझती हूँ

अगर तुम सभी बॉक्स चेक कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. सिंटैक्स एरर

सिंटैक्स एरर सबसे सामान्य प्रकार की एरर हैं जो तुम Nextflow कोड लिखते समय पाओगे। ये तब होती हैं जब कोड Nextflow DSL के अपेक्षित सिंटैक्स नियमों के अनुरूप नहीं होता। ये एरर तुम्हारे वर्कफ़्लो को बिल्कुल भी चलने से रोकती हैं, इसलिए यह सीखना महत्वपूर्ण है कि उन्हें जल्दी कैसे पहचाना और ठीक किया जाए।

### 1.1. गायब ब्रेसेज़

सबसे सामान्य सिंटैक्स एरर में से एक, और कभी-कभी डीबग करने के लिए अधिक जटिल एरर में से एक है **गायब या बेमेल ब्रैकेट**।

चलो एक व्यावहारिक उदाहरण से शुरू करते हैं।

#### पाइपलाइन चलाओ

```bash
nextflow run bad_syntax.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_syntax.nf` [stupefied_bhabha] DSL2 - revision: ca6327fad2

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

**सिंटैक्स एरर मैसेज के मुख्य तत्व:**

- **फ़ाइल और स्थान**: दिखाता है कि कौन सी फ़ाइल और लाइन/कॉलम में एरर है (`bad_syntax.nf:24:1`)
- **एरर विवरण**: बताता है कि पार्सर को क्या मिला जो उसे अपेक्षित नहीं था (`Unexpected input: '<EOF>'`)
- **EOF संकेतक**: `<EOF>` (End Of File) मैसेज इंगित करता है कि पार्सर फ़ाइल के अंत तक पहुंच गया जबकि अभी भी अधिक सामग्री की अपेक्षा थी - यह बंद न किए गए ब्रेसेज़ का एक क्लासिक संकेत है

#### कोड की जांच करो

अब, `bad_syntax.nf` की जांच करते हैं यह समझने के लिए कि एरर क्या कारण है:

```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
// प्रोसेस के लिए closing brace गायब है

workflow {

    // इनपुट चैनल बनाओ
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // इनपुट चैनल के साथ प्रोसेस को कॉल करो
    PROCESS_FILES(input_ch)
}
```

इस उदाहरण के लिए हमने तुम्हारे लिए एक कमेंट छोड़ा है जो दिखाता है कि एरर कहाँ है। Nextflow VSCode एक्सटेंशन भी तुम्हें कुछ संकेत दे रहा होगा कि क्या गलत हो सकता है, बेमेल ब्रेस को लाल रंग में दिखाकर और फ़ाइल के समय से पहले समाप्त होने को हाइलाइट करके:

![Bad syntax](../img/bad_syntax.png)

**ब्रैकेट एरर के लिए डीबगिंग रणनीति:**

1. VS Code के bracket matching का उपयोग करो (ब्रैकेट के बगल में कर्सर रखो)
2. ब्रैकेट-संबंधित मैसेज के लिए Problems panel जांचो
3. सुनिश्चित करो कि प्रत्येक opening `{` का एक corresponding closing `}` है

#### कोड ठीक करो

कमेंट को गायब closing brace से बदलो:

=== "बाद में"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }  // गायब closing brace जोड़ो

    workflow {

        // इनपुट चैनल बनाओ
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // इनपुट चैनल के साथ प्रोसेस को कॉल करो
        PROCESS_FILES(input_ch)
    }
    ```

=== "पहले"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    // प्रोसेस के लिए closing brace गायब है

    workflow {

        // इनपुट चैनल बनाओ
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // इनपुट चैनल के साथ प्रोसेस को कॉल करो
        PROCESS_FILES(input_ch)
    }
    ```

#### पाइपलाइन चलाओ

अब वर्कफ़्लो फिर से चलाओ यह पुष्टि करने के लिए कि यह काम करता है:

```bash
nextflow run bad_syntax.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_syntax.nf` [insane_faggin] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [48/cd7f54] PROCESS_FILES (1) | 3 of 3 ✔
    ```

### 1.2. गलत प्रोसेस कीवर्ड या निर्देशों का उपयोग

एक और सामान्य सिंटैक्स एरर है **अमान्य प्रोसेस परिभाषा**। यह तब हो सकता है जब तुम आवश्यक ब्लॉक परिभाषित करना भूल जाते हो या प्रोसेस परिभाषा में गलत निर्देशों का उपयोग करते हो।

#### पाइपलाइन चलाओ

```bash
nextflow run invalid_process.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `invalid_process.nf` [nasty_jepsen] DSL2 - revision: da9758d614

    Error invalid_process.nf:3:1: Invalid process definition -- check for missing or out-of-order section labels
    │   3 | process PROCESS_FILES {
    │     | ^^^^^^^^^^^^^^^^^^^^^^^
    │   4 |     inputs:
    │   5 |     val sample_name
    │   6 |
    ╰   7 |     output:

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### कोड की जांच करो

एरर "Invalid process definition" इंगित करती है और समस्या के आसपास का संदर्भ दिखाती है। लाइन 3-7 को देखते हुए, हम लाइन 4 पर `inputs:` देख सकते हैं, जो समस्या है। चलो `invalid_process.nf` की जांच करते हैं:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // एरर: 'inputs' नहीं, 'input' होना चाहिए
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // इनपुट चैनल बनाओ
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // इनपुट चैनल के साथ प्रोसेस को कॉल करो
    PROCESS_FILES(input_ch)
}
```

एरर संदर्भ में लाइन 4 को देखते हुए, हम समस्या देख सकते हैं: हम सही `input` निर्देश के बजाय `inputs` का उपयोग कर रहे हैं। Nextflow VSCode एक्सटेंशन भी इसे फ्लैग करेगा:

![Invalid process message](../img/invalid_process_message.png)

#### कोड ठीक करो

[दस्तावेज़ीकरण](https://www.nextflow.io/docs/latest/process.html#) का संदर्भ लेकर गलत कीवर्ड को सही से बदलो:

=== "बाद में"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // ठीक किया: 'inputs' को 'input' में बदला
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // इनपुट चैनल बनाओ
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // इनपुट चैनल के साथ प्रोसेस को कॉल करो
        PROCESS_FILES(input_ch)
    }
    ```

=== "पहले"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // एरर: 'inputs' नहीं, 'input' होना चाहिए
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // इनपुट चैनल बनाओ
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // इनपुट चैनल के साथ प्रोसेस को कॉल करो
        PROCESS_FILES(input_ch)
    }
    ```

#### पाइपलाइन चलाओ

अब वर्कफ़्लो फिर से चलाओ यह पुष्टि करने के लिए कि यह काम करता है:

```bash
nextflow run invalid_process.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `invalid_process.nf` [silly_fermi] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [b7/76cd9d] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.3. गलत वेरिएबल नामों का उपयोग

तुम्हारे script ब्लॉक में उपयोग किए जाने वाले वेरिएबल नाम वैध होने चाहिए, जो या तो इनपुट से या script से पहले डाले गए Groovy कोड से प्राप्त हों। लेकिन जब तुम पाइपलाइन विकास की शुरुआत में जटिलता से जूझ रहे होते हो, तो वेरिएबल नामकरण में गलतियाँ करना आसान है, और Nextflow तुम्हें जल्दी बता देगा।

#### पाइपलाइन चलाओ

```bash
nextflow run no_such_var.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `no_such_var.nf` [gloomy_meninsky] DSL2 - revision: 0c4d3bc28c

    Error no_such_var.nf:17:39: `undefined_var` is not defined
    │  17 |     echo "Using undefined variable: ${undefined_var}" >> ${output_pref
    ╰     |                                       ^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

एरर compile time पर पकड़ी जाती है और सीधे लाइन 17 पर अपरिभाषित वेरिएबल की ओर इशारा करती है, एक caret के साथ जो बिल्कुल वहाँ इंगित करता है जहाँ समस्या है।

#### कोड की जांच करो

चलो `no_such_var.nf` की जांच करते हैं:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // script से पहले Groovy कोड में वेरिएबल परिभाषित करो
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // एरर: undefined_var परिभाषित नहीं है
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

एरर मैसेज इंगित करता है कि वेरिएबल script टेम्पलेट में पहचाना नहीं गया है, और वहाँ तुम देख सकते हो - script ब्लॉक में `#!groovy ${undefined_var}` का उपयोग किया गया है, लेकिन कहीं और परिभाषित नहीं किया गया है।

#### कोड ठीक करो

अगर तुम्हें 'No such variable' एरर मिलती है, तो तुम इसे या तो वेरिएबल परिभाषित करके (इनपुट वेरिएबल नाम सही करके या script से पहले Groovy कोड संपादित करके), या script ब्लॉक से इसे हटाकर ठीक कर सकते हो अगर इसकी जरूरत नहीं है:

=== "बाद में"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // script से पहले Groovy कोड में वेरिएबल परिभाषित करो
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // undefined_var वाली लाइन हटाई
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "पहले"

    ```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // script से पहले Groovy कोड में वेरिएबल परिभाषित करो
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // एरर: undefined_var परिभाषित नहीं है
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### पाइपलाइन चलाओ

अब वर्कफ़्लो फिर से चलाओ यह पुष्टि करने के लिए कि यह काम करता है:

```bash
nextflow run no_such_var.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `no_such_var.nf` [suspicious_venter] DSL2 - revision: 6ba490f7c5

    executor >  local (3)
    [21/237300] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.4. Bash वेरिएबल का गलत उपयोग

Nextflow में शुरुआत करते समय, Nextflow (Groovy) और Bash वेरिएबल के बीच अंतर समझना मुश्किल हो सकता है। यह बुरे वेरिएबल एरर का एक और रूप उत्पन्न कर सकता है जो script ब्लॉक की Bash सामग्री में वेरिएबल का उपयोग करने की कोशिश करते समय प्रकट होता है।

#### पाइपलाइन चलाओ

```bash
nextflow run bad_bash_var.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_bash_var.nf` [infallible_mandelbrot] DSL2 - revision: 0853c11080

    Error bad_bash_var.nf:13:42: `prefix` is not defined
    │  13 |     echo "Processing ${sample_name}" > ${prefix}.txt
    ╰     |                                          ^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### कोड की जांच करो

एरर लाइन 13 की ओर इशारा करती है जहाँ `#!groovy ${prefix}` का उपयोग किया गया है। चलो `bad_bash_var.nf` की जांच करते हैं यह देखने के लिए कि समस्या क्या है:

```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    prefix="${sample_name}_output"
    echo "Processing ${sample_name}" > ${prefix}.txt  # एरर: ${prefix} Groovy सिंटैक्स है, Bash नहीं
    """
}
```

इस उदाहरण में, हम Bash में `prefix` वेरिएबल परिभाषित कर रहे हैं, लेकिन Nextflow प्रोसेस में `$` सिंटैक्स जिसका उपयोग हमने इसे संदर्भित करने के लिए किया (`#!groovy ${prefix}`) को Groovy वेरिएबल के रूप में व्याख्यायित किया जाता है, Bash के रूप में नहीं। वेरिएबल Groovy संदर्भ में मौजूद नहीं है, इसलिए हमें 'no such variable' एरर मिलती है।

#### कोड ठीक करो

अगर तुम Bash वेरिएबल का उपयोग करना चाहते हो, तो तुम्हें dollar sign को इस तरह escape करना होगा:

=== "बाद में"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > \${prefix}.txt  # ठीक किया: dollar sign escape किया
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "पहले"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > ${prefix}.txt  # एरर: ${prefix} Groovy सिंटैक्स है, Bash नहीं
        """
    }
    ```

यह Nextflow को इसे Bash वेरिएबल के रूप में व्याख्यायित करने के लिए कहता है।

#### पाइपलाइन चलाओ

अब वर्कफ़्लो फिर से चलाओ यह पुष्टि करने के लिए कि यह काम करता है:

```bash
nextflow run bad_bash_var.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_bash_var.nf` [naughty_franklin] DSL2 - revision: 58c1c83709

    executor >  local (3)
    [4e/560285] PROCESS_FILES (2) | 3 of 3 ✔
    ```

!!! tip "Groovy बनाम Bash वेरिएबल"

    सरल वेरिएबल मैनिपुलेशन जैसे string concatenation या prefix/suffix ऑपरेशन के लिए, script ब्लॉक में Bash वेरिएबल के बजाय script सेक्शन में Groovy वेरिएबल का उपयोग करना आमतौर पर अधिक पठनीय होता है:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    यह तरीका dollar signs को escape करने की जरूरत से बचाता है और कोड को पढ़ने और बनाए रखने में आसान बनाता है।

### 1.5. Workflow ब्लॉक के बाहर स्टेटमेंट

Nextflow VSCode एक्सटेंशन कोड संरचना की समस्याओं को हाइलाइट करता है जो एरर का कारण बनेंगी। एक सामान्य उदाहरण है `workflow {}` ब्लॉक के बाहर चैनल परिभाषित करना - यह अब एक सिंटैक्स एरर के रूप में लागू किया गया है।

#### पाइपलाइन चलाओ

```bash
nextflow run badpractice_syntax.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `badpractice_syntax.nf` [intergalactic_colden] DSL2 - revision: 5e4b291bde

    Error badpractice_syntax.nf:3:1: Statements cannot be mixed with script declarations -- move statements into a process or workflow
    │   3 | input_ch = channel.of('sample1', 'sample2', 'sample3')
    ╰     | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

एरर मैसेज स्पष्ट रूप से समस्या इंगित करता है: स्टेटमेंट (जैसे चैनल परिभाषाएं) को workflow या process ब्लॉक के बाहर script declarations के साथ मिश्रित नहीं किया जा सकता।

#### कोड की जांच करो

चलो `badpractice_syntax.nf` की जांच करते हैं यह देखने के लिए कि एरर क्या कारण है:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // एरर: workflow के बाहर चैनल परिभाषित

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // script से पहले Groovy कोड में वेरिएबल परिभाषित करो
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {
    PROCESS_FILES(input_ch)
}
```

VSCode एक्सटेंशन `input_ch` वेरिएबल को भी हाइलाइट करेगा जो workflow ब्लॉक के बाहर परिभाषित है:

![Non-lethal syntax error](../img/nonlethal.png)

#### कोड ठीक करो

चैनल परिभाषा को workflow ब्लॉक के अंदर ले जाओ:

=== "बाद में"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // script से पहले Groovy कोड में वेरिएबल परिभाषित करो
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // workflow ब्लॉक के अंदर ले जाया गया
        PROCESS_FILES(input_ch)
    }
    ```

=== "पहले"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // एरर: workflow के बाहर चैनल परिभाषित

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // script से पहले Groovy कोड में वेरिएबल परिभाषित करो
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        PROCESS_FILES(input_ch)
    }
    ```

#### पाइपलाइन चलाओ

यह पुष्टि करने के लिए कि फिक्स काम करता है, वर्कफ़्लो फिर से चलाओ:

```bash
nextflow run badpractice_syntax.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `badpractice_syntax.nf` [naughty_ochoa] DSL2 - revision: 5e4b291bde

    executor >  local (3)
    [6a/84a608] PROCESS_FILES (2) | 3 of 3 ✔
    ```

अपने इनपुट चैनल को workflow ब्लॉक के भीतर परिभाषित रखो, और सामान्य रूप से एक्सटेंशन द्वारा की गई किसी भी अन्य सिफारिश का पालन करो।

### सारांश

तुम Nextflow एरर मैसेज और IDE विज़ुअल संकेतकों का उपयोग करके सिंटैक्स एरर को व्यवस्थित रूप से पहचान और ठीक कर सकते हो। सामान्य सिंटैक्स एरर में गायब ब्रेसेज़, गलत प्रोसेस कीवर्ड, अपरिभाषित वेरिएबल, और Bash बनाम Nextflow वेरिएबल का अनुचित उपयोग शामिल हैं। VSCode एक्सटेंशन रनटाइम से पहले इनमें से कई को पकड़ने में मदद करता है। इन सिंटैक्स डीबगिंग कौशल के साथ, तुम सबसे सामान्य Nextflow सिंटैक्स एरर को जल्दी हल कर पाओगे और अधिक जटिल रनटाइम समस्याओं से निपटने के लिए आगे बढ़ पाओगे।

### आगे क्या है?

अधिक जटिल चैनल संरचना एरर को डीबग करना सीखो जो तब भी होती हैं जब सिंटैक्स सही हो।

---

## 2. चैनल संरचना एरर

चैनल संरचना एरर सिंटैक्स एरर से अधिक सूक्ष्म होती हैं क्योंकि कोड सिंटैक्टिकली सही है, लेकिन डेटा के आकार प्रोसेस की अपेक्षाओं से मेल नहीं खाते। Nextflow पाइपलाइन चलाने की कोशिश करेगा, लेकिन पा सकता है कि इनपुट की संख्या अपेक्षित से मेल नहीं खाती और विफल हो जाएगा। ये एरर आमतौर पर केवल रनटाइम पर दिखाई देती हैं और तुम्हारे वर्कफ़्लो से गुजरने वाले डेटा की समझ की आवश्यकता होती है।

!!! tip "`.view()` के साथ चैनल डीबग करना"

    इस सेक्शन में, याद रखो कि तुम अपने वर्कफ़्लो में किसी भी बिंदु पर चैनल सामग्री का निरीक्षण करने के लिए `.view()` ऑपरेटर का उपयोग कर सकते हो। यह चैनल संरचना समस्याओं को समझने के लिए सबसे शक्तिशाली डीबगिंग टूल में से एक है। हम इस तकनीक को सेक्शन 2.4 में विस्तार से देखेंगे, लेकिन उदाहरणों के माध्यम से काम करते समय इसका उपयोग करने के लिए स्वतंत्र महसूस करो।

    ```groovy
    my_channel.view()  // चैनल से क्या गुजर रहा है यह दिखाता है
    ```

### 2.1. इनपुट चैनल की गलत संख्या

यह एरर तब होती है जब तुम एक प्रोसेस की अपेक्षा से अलग संख्या में चैनल पास करते हो।

#### पाइपलाइन चलाओ

```bash
nextflow run bad_number_inputs.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_number_inputs.nf` [happy_swartz] DSL2 - revision: d83e58dcd3

    Error bad_number_inputs.nf:23:5: Incorrect number of call arguments, expected 1 but received 2
    │  23 |     PROCESS_FILES(samples_ch, files_ch)
    ╰     |     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### कोड की जांच करो

एरर मैसेज स्पष्ट रूप से बताता है कि call को 1 आर्गुमेंट की अपेक्षा थी लेकिन 2 मिले, और लाइन 23 की ओर इशारा करता है। चलो `bad_number_inputs.nf` की जांच करते हैं:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // प्रोसेस केवल 1 इनपुट की अपेक्षा करता है

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // दो अलग चैनल बनाओ
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // एरर: 2 चैनल पास किए जा रहे हैं लेकिन प्रोसेस केवल 1 की अपेक्षा करता है
    PROCESS_FILES(samples_ch, files_ch)
}
```

तुम्हें बेमेल `PROCESS_FILES` call दिखनी चाहिए, जो कई इनपुट चैनल प्रदान कर रही है जबकि प्रोसेस केवल एक परिभाषित करता है। VSCode एक्सटेंशन प्रोसेस call को लाल रंग में अंडरलाइन भी करेगा, और माउस ओवर करने पर एक diagnostic मैसेज प्रदान करेगा:

![Incorrect number of args message](../img/incorrect_num_args.png)

#### कोड ठीक करो

इस विशिष्ट उदाहरण के लिए, प्रोसेस एक single चैनल की अपेक्षा करता है और दूसरे चैनल की आवश्यकता नहीं है, इसलिए हम केवल `samples_ch` चैनल पास करके इसे ठीक कर सकते हैं:

=== "बाद में"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // प्रोसेस केवल 1 इनपुट की अपेक्षा करता है

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // दो अलग चैनल बनाओ
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ठीक किया: केवल वह चैनल पास करो जो प्रोसेस अपेक्षा करता है
        PROCESS_FILES(samples_ch)
    }
    ```

=== "पहले"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // प्रोसेस केवल 1 इनपुट की अपेक्षा करता है

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // दो अलग चैनल बनाओ
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // एरर: 2 चैनल पास किए जा रहे हैं लेकिन प्रोसेस केवल 1 की अपेक्षा करता है
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### पाइपलाइन चलाओ

```bash
nextflow run bad_number_inputs.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_number_inputs.nf` [big_euler] DSL2 - revision: e302bd87be

    executor >  local (3)
    [48/497f7b] PROCESS_FILES (3) | 3 of 3 ✔
    ```

इस उदाहरण की तुलना में अधिक सामान्य रूप से, तुम किसी प्रोसेस में अतिरिक्त इनपुट जोड़ सकते हो और workflow call को तदनुसार अपडेट करना भूल सकते हो, जो इस प्रकार की एरर का कारण बन सकता है। सौभाग्य से, यह समझने और ठीक करने में आसान एरर में से एक है, क्योंकि एरर मैसेज बेमेल के बारे में काफी स्पष्ट है।

### 2.2. चैनल समाप्ति (प्रोसेस अपेक्षा से कम बार चलता है)

कुछ चैनल संरचना एरर बहुत अधिक सूक्ष्म होती हैं और कोई एरर नहीं उत्पन्न करती हैं। इनमें से सबसे सामान्य शायद वह चुनौती है जो नए Nextflow उपयोगकर्ताओं को यह समझने में होती है कि queue channel समाप्त हो सकते हैं और आइटम खत्म हो सकते हैं, जिसका अर्थ है कि वर्कफ़्लो समय से पहले समाप्त हो जाता है।

#### पाइपलाइन चलाओ

```bash
nextflow run exhausted.nf
```

??? success "कमांड आउटपुट"

```console title="Exhausted channel output"
 N E X T F L O W   ~  version 25.10.4

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

यह वर्कफ़्लो बिना एरर के पूरा होता है, लेकिन यह केवल एक नमूना प्रोसेस करता है!

#### कोड की जांच करो

चलो `exhausted.nf` की जांच करते हैं यह देखने के लिए कि क्या यह सही है:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // script से पहले Groovy कोड में वेरिएबल परिभाषित करो
    output_prefix = "${reference}_${sample_name}"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {

    reference_ch = channel.of('baseline_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

प्रोसेस तीन बार के बजाय केवल एक बार चलता है क्योंकि `reference_ch` चैनल एक queue channel है जो पहले प्रोसेस एक्जीक्यूशन के बाद समाप्त हो जाता है। जब एक चैनल समाप्त हो जाता है, तो पूरी प्रोसेस रुक जाती है, भले ही अन्य चैनल में अभी भी आइटम हों।

यह एक सामान्य पैटर्न है जहाँ तुम्हारे पास एक single reference फ़ाइल है जिसे कई नमूनों में पुनः उपयोग करने की आवश्यकता है। समाधान है reference चैनल को एक value channel में बदलना जिसे अनिश्चित काल तक पुनः उपयोग किया जा सके।

#### कोड ठीक करो

इसे संबोधित करने के कुछ तरीके हैं जो इस बात पर निर्भर करते हैं कि कितनी फ़ाइलें प्रभावित हैं।

**विकल्प 1**: तुम्हारे पास एक single reference फ़ाइल है जिसे तुम बहुत अधिक पुनः उपयोग कर रहे हो। तुम बस एक value channel type बना सकते हो, जिसे बार-बार उपयोग किया जा सकता है। इसे करने के तीन तरीके हैं:

**1a** `channel.value()` का उपयोग करो:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Value channel को पुनः उपयोग किया जा सकता है
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** `first()` [ऑपरेटर](https://www.nextflow.io/docs/latest/reference/operator.html#first) का उपयोग करो:

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // value channel में बदलो
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** `collect()` [ऑपरेटर](https://www.nextflow.io/docs/latest/reference/operator.html#collect) का उपयोग करो:

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // value channel में बदलो
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**विकल्प 2**: अधिक जटिल परिदृश्यों में, शायद जहाँ sample channel में सभी नमूनों के लिए कई reference फ़ाइलें हैं, तुम `combine` ऑपरेटर का उपयोग करके एक नया चैनल बना सकते हो जो दोनों चैनल को tuples में जोड़ता है:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // cartesian product बनाता है

    PROCESS_FILES(combined_ch)
}
```

`.combine()` ऑपरेटर दोनों चैनल का cartesian product उत्पन्न करता है, इसलिए `reference_ch` में प्रत्येक आइटम `input_ch` में प्रत्येक आइटम के साथ जोड़ा जाएगा। यह प्रोसेस को reference का उपयोग करते हुए प्रत्येक नमूने के लिए चलने की अनुमति देता है।

इसके लिए प्रोसेस इनपुट को समायोजित करने की आवश्यकता है। हमारे उदाहरण में, प्रोसेस परिभाषा की शुरुआत को इस प्रकार समायोजित करने की आवश्यकता होगी:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

यह तरीका सभी स्थितियों में उपयुक्त नहीं हो सकता।

#### पाइपलाइन चलाओ

ऊपर दिए गए किसी एक फिक्स को आज़माओ और वर्कफ़्लो फिर से चलाओ:

```bash
nextflow run exhausted.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `exhausted.nf` [maniac_leavitt] DSL2 - revision: f372a56a7d

    executor >  local (3)
    [80/0779e9] PROCESS_FILES (3) | 3 of 3 ✔
    ```

अब तुम्हें केवल एक के बजाय तीनों नमूने प्रोसेस होते दिखने चाहिए।

### 2.3. गलत चैनल सामग्री संरचना

जब वर्कफ़्लो एक निश्चित स्तर की जटिलता तक पहुंचते हैं, तो प्रत्येक चैनल की आंतरिक संरचनाओं का ट्रैक रखना थोड़ा मुश्किल हो सकता है, और लोग आमतौर पर प्रोसेस की अपेक्षाओं और चैनल में वास्तव में क्या है के बीच बेमेल उत्पन्न करते हैं। यह पहले चर्चा की गई समस्या से अधिक सूक्ष्म है, जहाँ चैनल की संख्या गलत थी। इस मामले में, तुम्हारे पास इनपुट चैनल की सही संख्या हो सकती है, लेकिन उनमें से एक या अधिक की आंतरिक संरचना प्रोसेस की अपेक्षाओं से मेल नहीं खाती।

#### पाइपलाइन चलाओ

```bash
nextflow run bad_channel_shape.nf
```

??? failure "कमांड आउटपुट"

    ```console
    Launching `bad_channel_shape.nf` [hopeful_pare] DSL2 - revision: ffd66071a1

    executor >  local (3)
    executor >  local (3)
    [3f/c2dcb3] PROCESS_FILES (3) [  0%] 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `[sample1, file1.txt]_output.txt` expected by process `PROCESS_FILES (1)`


    Command executed:

      echo "Processing [sample1, file1.txt]" > [sample1, file1.txt]_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/d6/1fb69d1d93300bbc9d42f1875b981e

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### कोड की जांच करो

एरर मैसेज में square brackets यहाँ संकेत प्रदान करते हैं - प्रोसेस tuple को एक single value के रूप में मान रहा है, जो हम नहीं चाहते। चलो `bad_channel_shape.nf` की जांच करते हैं:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // single value की अपेक्षा करता है, tuple मिलता है

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // चैनल tuples emit करता है, लेकिन प्रोसेस single values की अपेक्षा करता है
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

तुम देख सकते हो कि हम tuples से बना एक चैनल उत्पन्न कर रहे हैं: `['sample1', 'file1.txt']`, लेकिन प्रोसेस एक single value, `val sample_name` की अपेक्षा करता है। एक्जीक्यूट किया गया कमांड दिखाता है कि प्रोसेस `[sample3, file3.txt]_output.txt` नाम की फ़ाइल बनाने की कोशिश कर रहा है, जो इच्छित आउटपुट नहीं है।

#### कोड ठीक करो

इसे ठीक करने के लिए, अगर प्रोसेस को दोनों इनपुट की आवश्यकता है तो हम प्रोसेस को tuple स्वीकार करने के लिए समायोजित कर सकते हैं:

=== "विकल्प 1: प्रोसेस में tuple स्वीकार करो"

    === "बाद में"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // ठीक किया: tuple स्वीकार करो

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // चैनल tuples emit करता है, लेकिन प्रोसेस single values की अपेक्षा करता है
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

    === "पहले"

        ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                val sample_name  // single value की अपेक्षा करता है, tuple मिलता है

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // चैनल tuples emit करता है, लेकिन प्रोसेस single values की अपेक्षा करता है
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "विकल्प 2: पहला तत्व निकालो"

    === "बाद में"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // चैनल tuples emit करता है, लेकिन प्रोसेस single values की अपेक्षा करता है
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // ठीक किया: पहला तत्व निकालो
        }
        ```

    === "पहले"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // चैनल tuples emit करता है, लेकिन प्रोसेस single values की अपेक्षा करता है
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### पाइपलाइन चलाओ

समाधानों में से एक चुनो और वर्कफ़्लो फिर से चलाओ:

```bash
nextflow run bad_channel_shape.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_channel_shape.nf` [clever_thompson] DSL2 - revision: 8cbcae3746

    executor >  local (3)
    [bb/80a958] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 2.4. चैनल डीबगिंग तकनीकें

#### चैनल निरीक्षण के लिए `.view()` का उपयोग

चैनल के लिए सबसे शक्तिशाली डीबगिंग टूल `.view()` ऑपरेटर है। `.view()` के साथ, तुम डीबगिंग में मदद के लिए सभी चरणों में अपने चैनल के आकार को समझ सकते हो।

#### पाइपलाइन चलाओ

इसे क्रिया में देखने के लिए `bad_channel_shape_viewed.nf` चलाओ:

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_channel_shape_viewed.nf` [maniac_poisson] DSL2 - revision: b4f24dc9da

    executor >  local (3)
    [c0/db76b3] PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

#### कोड की जांच करो

चलो `bad_channel_shape_viewed.nf` की जांच करते हैं यह देखने के लिए कि `.view()` का उपयोग कैसे किया जाता है:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // चैनल tuples emit करता है, लेकिन प्रोसेस single values की अपेक्षा करता है
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // डीबग: मूल चैनल सामग्री दिखाओ
    .map { tuple -> tuple[0] }        // रूपांतरण: पहला तत्व निकालो
    .view { "After mapping: $it" }    // डीबग: रूपांतरित चैनल सामग्री दिखाओ

    PROCESS_FILES(input_ch)
}
```

#### कोड ठीक करो

भविष्य में चैनल सामग्री को समझने के लिए अत्यधिक `.view()` ऑपरेशन का उपयोग करने से बचाने के लिए, कुछ कमेंट जोड़ना उचित है:

```groovy title="bad_channel_shape_viewed.nf (with comments)" linenums="16" hl_lines="8 9"
workflow {

    // चैनल tuples emit करता है, लेकिन प्रोसेस single values की अपेक्षा करता है
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

यह अधिक महत्वपूर्ण हो जाएगा जैसे-जैसे तुम्हारे वर्कफ़्लो जटिलता में बढ़ते हैं और चैनल संरचना अधिक अपारदर्शी हो जाती है।

#### पाइपलाइन चलाओ

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_channel_shape_viewed.nf` [marvelous_koch] DSL2 - revision: 03e79cdbad

    executor >  local (3)
    [ff/d67cec] PROCESS_FILES (2) | 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

### सारांश

वैध Nextflow सिंटैक्स के साथ कई चैनल संरचना एरर बनाई जा सकती हैं। तुम डेटा फ्लो को समझकर, निरीक्षण के लिए `.view()` ऑपरेटर का उपयोग करके, और एरर मैसेज पैटर्न को पहचानकर चैनल संरचना एरर को डीबग कर सकते हो जैसे square brackets जो अप्रत्याशित tuple संरचनाओं को इंगित करते हैं।

### आगे क्या है?

प्रोसेस परिभाषाओं द्वारा बनाई गई एरर के बारे में जानो।

---

## 3. प्रोसेस संरचना एरर

प्रोसेस से संबंधित अधिकांश एरर जो तुम्हें मिलेंगी, वे कमांड बनाने में की गई गलतियों या अंतर्निहित सॉफ़्टवेयर से संबंधित समस्याओं से संबंधित होंगी। फिर भी, ऊपर चैनल समस्याओं के समान, तुम प्रोसेस परिभाषा में ऐसी गलतियाँ कर सकते हो जो सिंटैक्स एरर के रूप में योग्य नहीं हैं, लेकिन रन टाइम पर एरर का कारण बनेंगी।

### 3.1. गायब आउटपुट फ़ाइलें

प्रोसेस लिखते समय एक सामान्य एरर है कुछ ऐसा करना जो प्रोसेस की अपेक्षाओं और जो उत्पन्न होता है के बीच बेमेल उत्पन्न करता है।

#### पाइपलाइन चलाओ

```bash
nextflow run missing_output.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `missing_output.nf` [zen_stone] DSL2 - revision: 37ff61f926

    executor >  local (3)
    executor >  local (3)
    [fd/2642e9] process > PROCESS_FILES (2) [ 66%] 2 of 3, failed: 2
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Missing output file(s) `sample3.txt` expected by process `PROCESS_FILES (3)`


    Command executed:

      echo "Processing sample3" > sample3_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/02/9604d49fb8200a74d737c72a6c98ed

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### कोड की जांच करो

एरर मैसेज इंगित करता है कि प्रोसेस `sample3.txt` नाम की आउटपुट फ़ाइल उत्पन्न करने की अपेक्षा करता था, लेकिन script वास्तव में `sample3_output.txt` बनाती है। चलो `missing_output.nf` में प्रोसेस परिभाषा की जांच करते हैं:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // अपेक्षित: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // बनाता है: sample3_output.txt
    """
}
```

तुम्हें देखना चाहिए कि `output:` ब्लॉक में आउटपुट फ़ाइल नाम और script में उपयोग किए गए नाम के बीच बेमेल है। यह बेमेल प्रोसेस को विफल करता है। अगर तुम्हें इस प्रकार की एरर मिलती है, तो वापस जाओ और जांचो कि तुम्हारी प्रोसेस परिभाषा और आउटपुट ब्लॉक के बीच आउटपुट मेल खाते हैं।

अगर समस्या अभी भी स्पष्ट नहीं है, तो वास्तव में बनाई गई आउटपुट फ़ाइलों की पहचान करने के लिए work directory की जांच करो:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

इस उदाहरण के लिए यह हमें उजागर करेगा कि जबकि हमारी अपेक्षित `sample3.txt` मौजूद नहीं थी, `sample3_output.txt` थी।

#### कोड ठीक करो

आउटपुट फ़ाइल नाम को सुसंगत बनाकर बेमेल ठीक करो:

=== "बाद में"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // ठीक किया: script आउटपुट से मेल खाओ

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }
    ```

=== "पहले"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}.txt"  // अपेक्षित: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // बनाता है: sample3_output.txt
        """
    }
    ```

#### पाइपलाइन चलाओ

```bash
nextflow run missing_output.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `missing_output.nf` [elated_hamilton] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [16/1c437c] PROCESS_FILES (3) | 3 of 3 ✔
    ```

### 3.2. गायब सॉफ़्टवेयर

एरर का एक और वर्ग सॉफ़्टवेयर प्रावधान में गलतियों के कारण होता है। `missing_software.nf` एक सिंटैक्टिकली वैध वर्कफ़्लो है, लेकिन यह उस `cowpy` कमांड के लिए कुछ बाहरी सॉफ़्टवेयर पर निर्भर करता है जिसका वह उपयोग करता है।

#### पाइपलाइन चलाओ

```bash
nextflow run missing_software.nf
```

??? failure "कमांड आउटपुट"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

प्रोसेस के पास उस कमांड तक पहुंच नहीं है जो हम निर्दिष्ट कर रहे हैं। कभी-कभी यह इसलिए होता है क्योंकि एक script वर्कफ़्लो `bin` डायरेक्टरी में मौजूद है, लेकिन executable नहीं बनाई गई है। अन्य बार यह इसलिए होता है क्योंकि सॉफ़्टवेयर उस कंटेनर या वातावरण में इंस्टॉल नहीं है जहाँ वर्कफ़्लो चल रहा है।

#### कोड की जांच करो

उस `127` exit code पर ध्यान दो - यह तुम्हें बिल्कुल समस्या बताता है। चलो `missing_software.nf` की जांच करते हैं:

```groovy title="missing_software.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

#### कोड ठीक करो

हम यहाँ थोड़े बेईमान रहे हैं, और वास्तव में कोड में कुछ भी गलत नहीं है। हमें बस प्रोसेस को इस तरह से चलाने के लिए आवश्यक कॉन्फ़िगरेशन निर्दिष्ट करने की आवश्यकता है कि उसके पास प्रश्न में कमांड तक पहुंच हो। इस मामले में प्रोसेस में एक container परिभाषा है, इसलिए हमें बस Docker enabled के साथ वर्कफ़्लो चलाना है।

#### पाइपलाइन चलाओ

हमने तुम्हारे लिए `nextflow.config` में एक Docker profile सेट किया है, इसलिए तुम वर्कफ़्लो इस तरह चला सकते हो:

```bash
nextflow run missing_software.nf -profile docker
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `missing_software.nf` [awesome_stonebraker] DSL2 - revision: 0296d12839

    executor >  local (3)
    [38/ab20d1] PROCESS_FILES (1) | 3 of 3 ✔
    ```

!!! note "नोट"

    Nextflow कंटेनर का उपयोग कैसे करता है इसके बारे में अधिक जानने के लिए, [Hello Nextflow](../../hello_nextflow/05_hello_containers.md) देखो

### 3.3. खराब रिसोर्स कॉन्फ़िगरेशन

प्रोडक्शन उपयोग में, तुम अपने प्रोसेस पर रिसोर्स कॉन्फ़िगर करोगे। उदाहरण के लिए `memory` तुम्हारे प्रोसेस के लिए उपलब्ध मेमोरी की अधिकतम मात्रा परिभाषित करता है, और अगर प्रोसेस उससे अधिक हो जाता है, तो तुम्हारा scheduler आमतौर पर प्रोसेस को kill कर देगा और `137` का exit code वापस करेगा। हम यहाँ इसे प्रदर्शित नहीं कर सकते क्योंकि हम `local` executor का उपयोग कर रहे हैं, लेकिन हम `time` के साथ कुछ इसी तरह दिखा सकते हैं।

#### पाइपलाइन चलाओ

`bad_resources.nf` में 1 millisecond के अवास्तविक time bound के साथ प्रोसेस कॉन्फ़िगरेशन है:

```bash
nextflow run bad_resources.nf -profile docker
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_resources.nf` [disturbed_elion] DSL2 - revision: 27d2066e86

    executor >  local (3)
    [c0/ded8e1] PROCESS_FILES (3) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (2)'

    Caused by:
      Process exceeded running time limit (1ms)

    Command executed:

      cowpy sample2 > sample2_output.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/53/f0a4cc56d6b3dc2a6754ff326f1349

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### कोड की जांच करो

चलो `bad_resources.nf` की जांच करते हैं:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // एरर: अवास्तविक time limit

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  # 1 सेकंड लेता है, लेकिन time limit 1ms है
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

हम जानते हैं कि प्रोसेस एक सेकंड से अधिक समय लेगा (हमने यह सुनिश्चित करने के लिए एक sleep जोड़ा है), लेकिन प्रोसेस 1 millisecond के बाद timeout होने के लिए सेट है। किसी ने अपने कॉन्फ़िगरेशन के साथ थोड़ा अवास्तविक काम किया है!

#### कोड ठीक करो

time limit को एक यथार्थवादी मान तक बढ़ाओ:

=== "बाद में"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // ठीक किया: यथार्थवादी time limit

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

=== "पहले"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '1 ms'  // एरर: अवास्तविक time limit

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  # 1 सेकंड लेता है, लेकिन time limit 1ms है
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### पाइपलाइन चलाओ

```bash
nextflow run bad_resources.nf -profile docker
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_resources.nf` [friendly_mcclintock] DSL2 - revision: 381567d2c1

    executor >  local (3)
    [c2/9b4c41] PROCESS_FILES (3) | 3 of 3 ✔
    ```

अगर तुम अपने एरर मैसेज ध्यान से पढ़ते हो तो इस तरह की विफलताएं तुम्हें लंबे समय तक परेशान नहीं करनी चाहिए। लेकिन सुनिश्चित करो कि तुम उन कमांड की रिसोर्स आवश्यकताओं को समझते हो जो तुम चला रहे हो ताकि तुम अपने resource निर्देशों को उचित रूप से कॉन्फ़िगर कर सको।

### 3.4. प्रोसेस डीबगिंग तकनीकें

जब प्रोसेस विफल होते हैं या अप्रत्याशित रूप से व्यवहार करते हैं, तो तुम्हें यह जांचने के लिए व्यवस्थित तकनीकों की आवश्यकता है कि क्या गलत हुआ। work directory में प्रोसेस एक्जीक्यूशन को डीबग करने के लिए आवश्यक सभी जानकारी होती है।

#### Work Directory निरीक्षण का उपयोग

प्रोसेस के लिए सबसे शक्तिशाली डीबगिंग टूल work directory की जांच करना है। जब कोई प्रोसेस विफल होता है, तो Nextflow उस विशिष्ट प्रोसेस एक्जीक्यूशन के लिए एक work directory बनाता है जिसमें यह समझने के लिए आवश्यक सभी फ़ाइलें होती हैं कि क्या हुआ।

#### पाइपलाइन चलाओ

work directory निरीक्षण प्रदर्शित करने के लिए पहले के `missing_output.nf` उदाहरण का उपयोग करते हैं (अगर जरूरत हो तो आउटपुट नामकरण बेमेल फिर से उत्पन्न करो):

```bash
nextflow run missing_output.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `missing_output.nf` [irreverent_payne] DSL2 - revision: 3d5117f7e2

    executor >  local (3)
    [5d/d544a4] PROCESS_FILES (2) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `sample1.txt` expected by process `PROCESS_FILES (1)`

    Command executed:

      echo "Processing sample1" > sample1_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/1e/2011154d0b0f001cd383d7364b5244

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Work Directory की जांच करो

जब तुम्हें यह एरर मिलती है, तो work directory में सभी डीबगिंग जानकारी होती है। एरर मैसेज से work directory path ढूंढो और इसकी सामग्री की जांच करो:

```bash
# एरर मैसेज से work directory ढूंढो
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

फिर तुम मुख्य फ़ाइलों की जांच कर सकते हो:

##### Command Script की जांच करो

`.command.sh` फ़ाइल दिखाती है कि वास्तव में कौन सा कमांड एक्जीक्यूट किया गया था:

```bash
# एक्जीक्यूट किया गया कमांड देखो
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

यह प्रकट करता है:

- **वेरिएबल प्रतिस्थापन**: क्या Nextflow वेरिएबल ठीक से expand किए गए थे
- **फ़ाइल पथ**: क्या इनपुट फ़ाइलें सही ढंग से स्थित थीं
- **कमांड संरचना**: क्या script सिंटैक्स सही है

देखने के लिए सामान्य समस्याएं:

- **गायब quotes**: spaces वाले वेरिएबल को उचित quoting की आवश्यकता है
- **गलत फ़ाइल पथ**: इनपुट फ़ाइलें जो मौजूद नहीं हैं या गलत स्थानों पर हैं
- **गलत वेरिएबल नाम**: वेरिएबल संदर्भों में टाइपो
- **गायब environment setup**: कमांड जो विशिष्ट वातावरण पर निर्भर करते हैं

##### एरर आउटपुट की जांच करो

`.command.err` फ़ाइल में वास्तविक एरर मैसेज होते हैं:

```bash
# एरर आउटपुट देखो
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

यह फ़ाइल दिखाएगी:

- **Exit codes**: 127 (command not found), 137 (killed), आदि
- **Permission errors**: फ़ाइल एक्सेस समस्याएं
- **Software errors**: एप्लिकेशन-विशिष्ट एरर मैसेज
- **Resource errors**: मेमोरी/time limit exceeded

##### Standard Output की जांच करो

`.command.out` फ़ाइल दिखाती है कि तुम्हारे कमांड ने क्या उत्पन्न किया:

```bash
# Standard output देखो
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

यह सत्यापित करने में मदद करता है:

- **अपेक्षित आउटपुट**: क्या कमांड ने सही परिणाम उत्पन्न किए
- **आंशिक एक्जीक्यूशन**: क्या कमांड शुरू हुआ लेकिन बीच में विफल हो गया
- **डीबग जानकारी**: तुम्हारी script से कोई diagnostic आउटपुट

##### Exit Code की जांच करो

`.exitcode` फ़ाइल में प्रोसेस के लिए exit code होता है:

```bash
# Exit code देखो
cat work/*/*/.exitcode
```

सामान्य exit codes और उनके अर्थ:

- **Exit code 127**: Command not found - सॉफ़्टवेयर इंस्टॉलेशन जांचो
- **Exit code 137**: Process killed - मेमोरी/time limits जांचो

##### फ़ाइल अस्तित्व की जांच करो

जब प्रोसेस गायब आउटपुट फ़ाइलों के कारण विफल होते हैं, तो जांचो कि वास्तव में कौन सी फ़ाइलें बनाई गई थीं:

```bash
# work directory में सभी फ़ाइलें सूचीबद्ध करो
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

यह पहचानने में मदद करता है:

- **फ़ाइल नामकरण बेमेल**: अपेक्षित से अलग नामों वाली आउटपुट फ़ाइलें
- **Permission समस्याएं**: फ़ाइलें जो नहीं बनाई जा सकीं
- **Path समस्याएं**: गलत डायरेक्टरी में बनाई गई फ़ाइलें

हमारे पहले के उदाहरण में, इसने हमें पुष्टि की कि जबकि हमारी अपेक्षित `sample3.txt` मौजूद नहीं थी, `sample3_output.txt` थी:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### सारांश

प्रोसेस डीबगिंग के लिए यह समझने के लिए work directories की जांच करना आवश्यक है कि क्या गलत हुआ। मुख्य फ़ाइलों में `.command.sh` (एक्जीक्यूट की गई script), `.command.err` (एरर मैसेज), और `.command.out` (standard output) शामिल हैं। Exit codes जैसे 127 (command not found) और 137 (process killed) विफलता के प्रकार के बारे में तत्काल diagnostic संकेत प्रदान करते हैं।

### आगे क्या है?

Nextflow के बिल्ट-इन डीबगिंग टूल और समस्या निवारण के व्यवस्थित तरीकों के बारे में जानो।

---

## 4. बिल्ट-इन डीबगिंग टूल और उन्नत तकनीकें

Nextflow वर्कफ़्लो एक्जीक्यूशन को डीबग और विश्लेषण करने के लिए कई शक्तिशाली बिल्ट-इन टूल प्रदान करता है। ये टूल तुम्हें यह समझने में मदद करते हैं कि क्या गलत हुआ, कहाँ गलत हुआ, और इसे कुशलतापूर्वक कैसे ठीक किया जाए।

### 4.1. रियल-टाइम प्रोसेस आउटपुट

कभी-कभी तुम्हें यह देखने की जरूरत होती है कि चल रहे प्रोसेस के अंदर क्या हो रहा है। तुम रियल-टाइम प्रोसेस आउटपुट enable कर सकते हो, जो तुम्हें दिखाता है कि प्रत्येक कार्य एक्जीक्यूट होते समय वास्तव में क्या कर रहा है।

#### पाइपलाइन चलाओ

हमारे पहले के उदाहरणों से `bad_channel_shape_viewed.nf` ने `.view()` का उपयोग करके चैनल सामग्री प्रिंट की, लेकिन हम प्रोसेस के भीतर से वेरिएबल echo करने के लिए `debug` निर्देश का भी उपयोग कर सकते हैं, जिसे हम `bad_channel_shape_viewed_debug.nf` में प्रदर्शित करते हैं। वर्कफ़्लो चलाओ:

```bash
nextflow run bad_channel_shape_viewed_debug.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_channel_shape_viewed_debug.nf` [agitated_crick] DSL2 - revision: ea3676d9ec

    executor >  local (3)
    [c6/2dac51] process > PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    Sample name inside process is sample2

    Sample name inside process is sample1

    Sample name inside process is sample3
    ```

#### कोड की जांच करो

चलो `bad_channel_shape_viewed_debug.nf` की जांच करते हैं यह देखने के लिए कि `debug` निर्देश कैसे काम करता है:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // रियल-टाइम आउटपुट enable करो

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Sample name inside process is ${sample_name}"
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}
```

`debug` निर्देश एक प्रोसेस के वातावरण को समझने का एक त्वरित और सुविधाजनक तरीका हो सकता है।

### 4.2. Preview Mode

कभी-कभी तुम किसी भी प्रोसेस के चलने से पहले समस्याओं को पकड़ना चाहते हो। Nextflow इस प्रकार की सक्रिय डीबगिंग के लिए एक flag प्रदान करता है: `-preview`।

#### पाइपलाइन चलाओ

Preview mode तुम्हें कमांड एक्जीक्यूट किए बिना वर्कफ़्लो लॉजिक का परीक्षण करने देता है। यह तुम्हारे वर्कफ़्लो की संरचना को जल्दी से जांचने और यह सुनिश्चित करने के लिए काफी उपयोगी हो सकता है कि प्रोसेस बिना किसी वास्तविक कमांड चलाए सही ढंग से जुड़े हैं।

!!! note "नोट"

    अगर तुमने पहले `bad_syntax.nf` ठीक किया था, तो इस कमांड को चलाने से पहले script ब्लॉक के बाद closing brace हटाकर सिंटैक्स एरर फिर से डालो।

यह कमांड चलाओ:

```bash
nextflow run bad_syntax.nf -preview
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `bad_syntax.nf` [magical_mercator] DSL2 - revision: 550b9a8873

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Preview mode किसी भी प्रोसेस को चलाए बिना सिंटैक्स एरर को जल्दी पकड़ने के लिए विशेष रूप से उपयोगी है। यह एक्जीक्यूशन से पहले वर्कफ़्लो संरचना और प्रोसेस कनेक्शन को validate करता है।

### 4.3. लॉजिक परीक्षण के लिए Stub Running

कभी-कभी एरर को डीबग करना मुश्किल होता है क्योंकि कमांड बहुत लंबे समय लेते हैं, विशेष सॉफ़्टवेयर की आवश्यकता होती है, या जटिल कारणों से विफल होते हैं। Stub running तुम्हें वास्तविक कमांड एक्जीक्यूट किए बिना वर्कफ़्लो लॉजिक का परीक्षण करने देता है।

#### पाइपलाइन चलाओ

जब तुम एक Nextflow प्रोसेस विकसित कर रहे हो, तो तुम `stub` निर्देश का उपयोग करके 'dummy' कमांड परिभाषित कर सकते हो जो वास्तविक कमांड चलाए बिना सही रूप के आउटपुट उत्पन्न करते हैं। यह तरीका विशेष रूप से मूल्यवान है जब तुम वास्तविक सॉफ़्टवेयर की जटिलताओं से निपटने से पहले यह सत्यापित करना चाहते हो कि तुम्हारा वर्कफ़्लो लॉजिक सही है।

उदाहरण के लिए, याद करो हमारा `missing_software.nf` पहले से? वह जहाँ हमारे पास गायब सॉफ़्टवेयर था जिसने वर्कफ़्लो को तब तक चलने से रोका जब तक हमने `-profile docker` नहीं जोड़ा? `missing_software_with_stub.nf` एक बहुत समान वर्कफ़्लो है। अगर हम इसे उसी तरह चलाते हैं, तो हम वही एरर उत्पन्न करेंगे:

```bash
nextflow run missing_software_with_stub.nf
```

??? failure "कमांड आउटपुट"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

हालांकि, यह वर्कफ़्लो `-stub-run` के साथ चलाने पर एरर उत्पन्न नहीं करेगा, यहाँ तक कि `docker` profile के बिना भी:

```bash
nextflow run missing_software_with_stub.nf -stub-run
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `missing_software_with_stub.nf` [astonishing_shockley] DSL2 - revision: f1f4f05d7d

    executor >  local (3)
    [b5/2517a3] PROCESS_FILES (3) | 3 of 3 ✔
    ```

#### कोड की जांच करो

चलो `missing_software_with_stub.nf` की जांच करते हैं:

```groovy title="missing_software.nf (with stub)" hl_lines="16-19" linenums="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """

    stub:
    """
    touch ${sample_name}_output.txt
    """
}
```

`missing_software.nf` की तुलना में, इस प्रोसेस में एक `stub:` निर्देश है जो एक कमांड निर्दिष्ट करता है जिसे `script:` में निर्दिष्ट के बजाय उपयोग किया जाएगा, इस स्थिति में कि Nextflow stub mode में चलाया जाए।

`touch` कमांड जिसका हम यहाँ उपयोग कर रहे हैं वह किसी भी सॉफ़्टवेयर या उचित इनपुट पर निर्भर नहीं करता, और सभी स्थितियों में चलेगा, जिससे हम प्रोसेस आंतरिक के बारे में चिंता किए बिना वर्कफ़्लो लॉजिक को डीबग कर सकते हैं।

**Stub running डीबग करने में मदद करता है:**

- चैनल संरचना और डेटा फ्लो
- प्रोसेस कनेक्शन और dependencies
- पैरामीटर प्रसार
- सॉफ़्टवेयर dependencies के बिना वर्कफ़्लो लॉजिक

### 4.4. व्यवस्थित डीबगिंग तरीका

अब जब तुमने व्यक्तिगत डीबगिंग तकनीकें सीखी हैं - trace files और work directories से लेकर preview mode, stub running, और resource monitoring तक - चलो उन्हें एक व्यवस्थित पद्धति में एक साथ जोड़ते हैं। एक संरचित तरीका होने से तुम जटिल एरर से अभिभूत होने से बचते हो और यह सुनिश्चित करते हो कि तुम महत्वपूर्ण संकेत नहीं चूकते।

यह पद्धति हमारे द्वारा कवर किए गए सभी टूल को एक कुशल वर्कफ़्लो में जोड़ती है:

**चार-चरण डीबगिंग विधि:**

**चरण 1: सिंटैक्स एरर समाधान (5 मिनट)**

1. VSCode या तुम्हारे IDE में लाल underlines की जांच करो
2. सिंटैक्स समस्याओं की पहचान करने के लिए `nextflow run workflow.nf -preview` चलाओ
3. सभी सिंटैक्स एरर ठीक करो (गायब ब्रेसेज़, trailing commas, आदि)
4. आगे बढ़ने से पहले सुनिश्चित करो कि वर्कफ़्लो सफलतापूर्वक parse होता है

**चरण 2: त्वरित मूल्यांकन (5 मिनट)**

1. रनटाइम एरर मैसेज ध्यान से पढ़ो
2. जांचो कि यह रनटाइम, लॉजिक, या resource एरर है
3. बुनियादी वर्कफ़्लो लॉजिक का परीक्षण करने के लिए preview mode का उपयोग करो

**चरण 3: विस्तृत जांच (15-30 मिनट)**

1. विफल कार्य की work directory ढूंढो
2. log files की जांच करो
3. चैनल निरीक्षण के लिए `.view()` ऑपरेटर जोड़ो
4. एक्जीक्यूशन के बिना वर्कफ़्लो लॉजिक का परीक्षण करने के लिए `-stub-run` का उपयोग करो

**चरण 4: ठीक करो और validate करो (15 मिनट)**

1. न्यूनतम लक्षित फिक्स करो
2. resume के साथ परीक्षण करो: `nextflow run workflow.nf -resume`
3. पूर्ण वर्कफ़्लो एक्जीक्यूशन सत्यापित करो

!!! tip "कुशल डीबगिंग के लिए Resume का उपयोग"

    एक बार जब तुमने समस्या की पहचान कर ली, तो तुम्हें अपने वर्कफ़्लो के सफल हिस्सों को फिर से चलाने में समय बर्बाद किए बिना अपने फिक्स का परीक्षण करने का एक कुशल तरीका चाहिए। Nextflow की `-resume` कार्यक्षमता डीबगिंग के लिए अमूल्य है।

    तुम `-resume` से परिचित होगे अगर तुमने [Hello Nextflow](../../hello_nextflow/index.md) के माध्यम से काम किया है, और यह महत्वपूर्ण है कि तुम डीबगिंग करते समय इसका अच्छा उपयोग करो ताकि तुम्हारी समस्या प्रोसेस से पहले के प्रोसेस चलने का इंतजार करते हुए समय बर्बाद न हो।

    **Resume डीबगिंग रणनीति:**

    1. विफलता तक वर्कफ़्लो चलाओ
    2. विफल कार्य के लिए work directory की जांच करो
    3. विशिष्ट समस्या ठीक करो
    4. केवल फिक्स का परीक्षण करने के लिए resume करो
    5. वर्कफ़्लो पूरा होने तक दोहराओ

#### डीबगिंग कॉन्फ़िगरेशन Profile

इस व्यवस्थित तरीके को और भी कुशल बनाने के लिए, तुम एक समर्पित डीबगिंग कॉन्फ़िगरेशन बना सकते हो जो स्वचालित रूप से तुम्हारी जरूरत के सभी टूल enable करता है:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // डीबगिंग के लिए conservative resources
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

फिर तुम इस profile enabled के साथ पाइपलाइन चला सकते हो:

```bash
nextflow run workflow.nf -profile debug
```

यह profile रियल-टाइम आउटपुट enable करता है, work directories को preserve करता है, और आसान डीबगिंग के लिए parallelization को सीमित करता है।

### 4.5. व्यावहारिक डीबगिंग अभ्यास

अब व्यवस्थित डीबगिंग तरीके कोव्यवहार में लाने का समय है। वर्कफ़्लो `buggy_workflow.nf` में कई सामान्य एरर हैं जो वास्तविक दुनिया के विकास में आने वाली समस्याओं के प्रकारों का प्रतिनिधित्व करती हैं।

!!! exercise "अभ्यास"

    `buggy_workflow.nf` में सभी एरर की पहचान करने और उन्हें ठीक करने के लिए व्यवस्थित डीबगिंग तरीके का उपयोग करो। यह वर्कफ़्लो एक CSV फ़ाइल से नमूना डेटा प्रोसेस करने की कोशिश करता है लेकिन इसमें सामान्य डीबगिंग परिदृश्यों का प्रतिनिधित्व करने वाले कई जानबूझकर बग हैं।

    पहली एरर देखने के लिए वर्कफ़्लो चलाकर शुरू करो:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "कमांड आउटपुट"

        ```console
        N E X T F L O W   ~  version 25.10.4

        Launching `buggy_workflow.nf` [wise_ramanujan] DSL2 - revision: d51a8e83fd

        ERROR ~ Range [11, 12) out of bounds for length 11

         -- Check '.nextflow.log' file for details
        ```

        यह cryptic एरर `params{}` ब्लॉक में लाइन 11-12 के आसपास एक parsing समस्या इंगित करती है। v2 parser संरचनात्मक समस्याओं को जल्दी पकड़ता है।

    तुमने जो चार-चरण डीबगिंग विधि सीखी है उसे लागू करो:

    **चरण 1: सिंटैक्स एरर समाधान**
    - VSCode या तुम्हारे IDE में लाल underlines की जांच करो
    - सिंटैक्स समस्याओं की पहचान करने के लिए `nextflow run workflow.nf -preview` चलाओ
    - सभी सिंटैक्स एरर ठीक करो (गायब ब्रेसेज़, trailing commas, आदि)
    - आगे बढ़ने से पहले सुनिश्चित करो कि वर्कफ़्लो सफलतापूर्वक parse होता है

    **चरण 2: त्वरित मूल्यांकन**
    - रनटाइम एरर मैसेज ध्यान से पढ़ो
    - पहचानो कि एरर रनटाइम, लॉजिक, या resource-संबंधित हैं
    - बुनियादी वर्कफ़्लो लॉजिक का परीक्षण करने के लिए `-preview` mode का उपयोग करो

    **चरण 3: विस्तृत जांच**
    - विफल कार्यों के लिए work directories की जांच करो
    - चैनल निरीक्षण के लिए `.view()` ऑपरेटर जोड़ो
    - work directories में log files जांचो
    - एक्जीक्यूशन के बिना वर्कफ़्लो लॉजिक का परीक्षण करने के लिए `-stub-run` का उपयोग करो

    **चरण 4: ठीक करो और validate करो**
    - लक्षित फिक्स करो
    - फिक्स को कुशलतापूर्वक परीक्षण करने के लिए `-resume` का उपयोग करो
    - पूर्ण वर्कफ़्लो एक्जीक्यूशन सत्यापित करो

    **तुम्हारे पास उपलब्ध डीबगिंग टूल:**
    ```bash
    # सिंटैक्स जांच के लिए preview mode
    nextflow run buggy_workflow.nf -preview

    # विस्तृत आउटपुट के लिए debug profile
    nextflow run buggy_workflow.nf -profile debug

    # लॉजिक परीक्षण के लिए stub running
    nextflow run buggy_workflow.nf -stub-run

    # फिक्स के बाद resume
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution "समाधान"
        `buggy_workflow.nf` में 9 या 10 अलग-अलग एरर हैं (गिनने के तरीके पर निर्भर करते हुए) जो सभी प्रमुख डीबगिंग श्रेणियों को कवर करती हैं। यहाँ प्रत्येक एरर और उसे कैसे ठीक करें का व्यवस्थित विवरण है।

        चलो उन सिंटैक्स एरर से शुरू करते हैं:

        **एरर 1: सिंटैक्स एरर - Trailing Comma**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // एरर: trailing comma
        ```
        **फिक्स:** trailing comma हटाओ
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **एरर 2: सिंटैक्स एरर - गायब Closing Brace**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // एरर: processFiles प्रोसेस के लिए closing brace गायब है
        ```
        **फिक्स:** गायब closing brace जोड़ो
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // गायब closing brace जोड़ो
        ```

        **एरर 3: वेरिएबल नाम एरर**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // एरर: sample_id होना चाहिए
        cat ${input_file} > ${sample}_result.txt  // एरर: sample_id होना चाहिए
        ```
        **फिक्स:** सही इनपुट वेरिएबल नाम का उपयोग करो
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **एरर 4: अपरिभाषित वेरिएबल एरर**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // एरर: sample_ids अपरिभाषित है
        ```
        **फिक्स:** सही चैनल का उपयोग करो और sample IDs निकालो
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        इस बिंदु पर वर्कफ़्लो चलेगा, लेकिन हमें अभी भी एरर मिलेंगी (जैसे `processFiles` में `Path value cannot be null`), जो खराब चैनल संरचना के कारण हैं।

        **एरर 5: चैनल संरचना एरर - गलत Map आउटपुट**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // एरर: processFiles tuple की अपेक्षा करता है
        ```
        **फिक्स:** वह tuple संरचना वापस करो जो processFiles अपेक्षा करता है
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        लेकिन यह ऊपर `heavyProcess()` चलाने के लिए हमारे फिक्स को तोड़ देगा, इसलिए हमें उस प्रोसेस को केवल sample IDs पास करने के लिए map का उपयोग करना होगा:

        **एरर 6: heavyProcess के लिए खराब चैनल संरचना**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // एरर: input_ch में अब प्रति emission 2 तत्व हैं - heavyProcess को केवल 1 (पहला) चाहिए
        ```
        **फिक्स:** सही चैनल का उपयोग करो और sample IDs निकालो
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        अब हम थोड़ा आगे बढ़ते हैं लेकिन `No such variable: i` के बारे में एरर मिलती है, क्योंकि हमने Bash वेरिएबल को escape नहीं किया।

        **एरर 7: Bash वेरिएबल Escaping एरर**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // एरर: $i escape नहीं किया
        ```
        **फिक्स:** bash वेरिएबल escape करो
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        अब हमें `Process exceeded running time limit (1ms)` मिलती है, इसलिए हम संबंधित प्रोसेस के लिए run time limit ठीक करते हैं:

        **एरर 8: Resource कॉन्फ़िगरेशन एरर**
        ```groovy linenums="36"
        time '1 ms'  // एरर: अवास्तविक time limit
        ```
        **फिक्स:** एक यथार्थवादी time limit तक बढ़ाओ
        ```groovy linenums="36"
        time '100 s'
        ```

        अगला हमारे पास हल करने के लिए एक `Missing output file(s)` एरर है:

        **एरर 9: आउटपुट फ़ाइल नाम बेमेल**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // एरर: गलत फ़ाइलनाम, output declaration से मेल खाना चाहिए
        ```
        **फिक्स:** output declaration से मेल खाओ
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        पहले दो प्रोसेस चले, लेकिन तीसरा नहीं।

        **एरर 10: आउटपुट फ़ाइल नाम बेमेल**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // एरर: pwd से इनपुट लेने की कोशिश कर रहा है न कि किसी प्रोसेस से
        handleFiles(file_ch)
        ```
        **फिक्स:** पिछले प्रोसेस से आउटपुट लो
        ```groovy linenums="88"
        file_ch = handleFiles(heavy_ch)
        ```

        इसके साथ, पूरा वर्कफ़्लो चलना चाहिए।

        **पूर्ण सही वर्कफ़्लो:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * डीबगिंग अभ्यास के लिए बगी वर्कफ़्लो
        * इस वर्कफ़्लो में सीखने के उद्देश्यों के लिए कई जानबूझकर बग हैं
        */

        params{
            // गायब validation के साथ पैरामीटर
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * इनपुट/आउटपुट बेमेल वाला प्रोसेस
        */
        process processFiles {

            input:
                tuple val(sample_id), path(input_file)

            output:
                path "${sample_id}_result.txt"

            script:
            """
            echo "Processing: ${sample_id}"
            cat ${input_file} > ${sample_id}_result.txt
            """
        }

        /*
        * resource समस्याओं वाला प्रोसेस
        */
        process heavyProcess {

            time '100 s'

            input:
                val sample_id

            output:
                path "${sample_id}_heavy.txt"

            script:
            """
            # भारी computation simulate करो
            for i in {1..1000000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * फ़ाइल handling समस्याओं वाला प्रोसेस
        */
        process handleFiles {

            input:
                path input_file

            output:
                path "processed_${input_file}"

            script:
            """
            if [ -f "${input_file}" ]; then
                cp ${input_file} processed_${input_file}
            fi
            """
        }

        /*
        * चैनल समस्याओं वाला मुख्य वर्कफ़्लो
        */
        workflow {
            main:
            // गलत उपयोग वाला चैनल
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map{it[0]})

            file_ch = handleFiles(heavy_ch)

            publish:
            processed = processed_ch
            heavy = heavy_ch
            files = file_ch
        }

        output {
            processed {
                path 'processed'
            }
            heavy {
                path 'heavy'
            }
            files {
                path 'files'
            }
        }
        ```

**कवर की गई एरर श्रेणियां:**

- **सिंटैक्स एरर**: गायब ब्रेसेज़, trailing commas, अपरिभाषित वेरिएबल
- **चैनल संरचना एरर**: गलत डेटा आकार, अपरिभाषित चैनल
- **प्रोसेस एरर**: आउटपुट फ़ाइल बेमेल, वेरिएबल escaping
- **Resource एरर**: अवास्तविक time limits

**मुख्य डीबगिंग सबक:**

1. **एरर मैसेज ध्यान से पढ़ो** - वे अक्सर सीधे समस्या की ओर इशारा करते हैं
2. **व्यवस्थित तरीकों का उपयोग करो** - एक बार में एक एरर ठीक करो और `-resume` के साथ परीक्षण करो
3. **डेटा फ्लो समझो** - चैनल संरचना एरर अक्सर सबसे सूक्ष्म होती हैं
4. **work directories जांचो** - जब प्रोसेस विफल होते हैं, तो logs तुम्हें बिल्कुल बताते हैं कि क्या गलत हुआ

---

## सारांश

इस साइड क्वेस्ट में, तुमने Nextflow वर्कफ़्लो को डीबग करने के लिए व्यवस्थित तकनीकों का एक सेट सीखा है।
अपने काम में इन तकनीकों को लागू करने से तुम अपने कंप्यूटर से लड़ने में कम समय बिताओगे, समस्याओं को तेज़ी से हल करोगे और भविष्य की समस्याओं से खुद को बचाओगे।

### मुख्य पैटर्न

**1. सिंटैक्स एरर की पहचान और उन्हें कैसे ठीक करें:**

- Nextflow एरर मैसेज की व्याख्या करना और समस्याओं का पता लगाना
- सामान्य सिंटैक्स एरर: गायब ब्रेसेज़, गलत कीवर्ड, अपरिभाषित वेरिएबल
- Nextflow (Groovy) और Bash वेरिएबल के बीच अंतर करना
- प्रारंभिक एरर डिटेक्शन के लिए VS Code एक्सटेंशन फीचर्स का उपयोग करना

```groovy
// गायब brace - IDE में लाल underlines देखो
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- गायब!

// गलत कीवर्ड
inputs:  // 'input:' होना चाहिए

// अपरिभाषित वेरिएबल - Bash वेरिएबल के लिए backslash से escape करो
echo "${undefined_var}"      // Nextflow वेरिएबल (परिभाषित न होने पर एरर)
echo "\${bash_var}"          // Bash वेरिएबल (escaped)
```

**2. चैनल संरचना समस्याओं को कैसे डीबग करें:**

- चैनल cardinality और exhaustion समस्याओं को समझना
- चैनल सामग्री संरचना बेमेल को डीबग करना
- चैनल निरीक्षण के लिए `.view()` ऑपरेटर का उपयोग करना
- आउटपुट में square brackets जैसे एरर पैटर्न को पहचानना

```groovy
// चैनल सामग्री निरीक्षण करो
my_channel.view { "Content: $it" }

// queue को value channel में बदलो (exhaustion रोकता है)
reference_ch = channel.value('ref.fa')
// या
reference_ch = channel.of('ref.fa').first()
```

**3. प्रोसेस एक्जीक्यूशन समस्याओं का निवारण कैसे करें:**

- गायब आउटपुट फ़ाइल एरर का निदान करना
- exit codes समझना (गायब सॉफ़्टवेयर के लिए 127, मेमोरी समस्याओं के लिए 137)
- work directories और command files की जांच करना
- resources को उचित रूप से कॉन्फ़िगर करना

```bash
# वास्तव में क्या एक्जीक्यूट किया गया था जांचो
cat work/ab/cdef12/.command.sh

# एरर आउटपुट जांचो
cat work/ab/cdef12/.command.err

# Exit code 127 = command not found
# Exit code 137 = killed (memory/time limit)
```

**4. Nextflow के बिल्ट-इन डीबगिंग टूल का उपयोग कैसे करें:**

- preview mode और रियल-टाइम डीबगिंग का लाभ उठाना
- लॉजिक परीक्षण के लिए stub running लागू करना
- कुशल डीबगिंग cycles के लिए resume लागू करना
- चार-चरण व्यवस्थित डीबगिंग पद्धति का पालन करना

!!! tip "त्वरित डीबगिंग संदर्भ"

    **सिंटैक्स एरर?** → VSCode चेतावनियाँ जांचो, `nextflow run workflow.nf -preview` चलाओ

    **चैनल समस्याएं?** → सामग्री निरीक्षण के लिए `.view()` का उपयोग करो: `my_channel.view()`

    **प्रोसेस विफलताएं?** → work directory फ़ाइलें जांचो:

    - `.command.sh` - एक्जीक्यूट की गई script
    - `.command.err` - एरर मैसेज
    - `.exitcode` - exit status (127 = command not found, 137 = killed)

    **रहस्यमय व्यवहार?** → वर्कफ़्लो लॉजिक का परीक्षण करने के लिए `-stub-run` के साथ चलाओ

    **फिक्स किए?** → परीक्षण में समय बचाने के लिए `-resume` का उपयोग करो: `nextflow run workflow.nf -resume`

---

### अतिरिक्त संसाधन

- [Nextflow troubleshooting guide](https://www.nextflow.io/docs/latest/troubleshooting.html): आधिकारिक troubleshooting दस्तावेज़ीकरण
- [Understanding Nextflow channels](https://www.nextflow.io/docs/latest/channel.html): चैनल प्रकारों और व्यवहार में गहरी जानकारी
- [Process directives reference](https://www.nextflow.io/docs/latest/process.html#directives): सभी उपलब्ध प्रोसेस कॉन्फ़िगरेशन विकल्प
- [nf-test](https://www.nf-test.com/): Nextflow पाइपलाइन के लिए परीक्षण framework
- [Nextflow Slack community](https://www.nextflow.io/slack-invite.html): समुदाय से मदद लो

प्रोडक्शन वर्कफ़्लो के लिए, विचार करो:

- बड़े पैमाने पर monitoring और डीबगिंग के लिए [Seqera Platform](https://seqera.io/platform/) सेट करना
- reproducible सॉफ़्टवेयर वातावरण के लिए [Wave containers](https://seqera.io/wave/) का उपयोग करना

**याद रखो:** प्रभावी डीबगिंग एक कौशल है जो अभ्यास के साथ बेहतर होता है। यहाँ तुमने जो व्यवस्थित पद्धति और व्यापक toolkit हासिल की है, वह तुम्हारी Nextflow विकास यात्रा में तुम्हारी अच्छी सेवा करेगी।

---

## आगे क्या है?

[साइड क्वेस्ट के मेनू](../index.md) पर वापस जाओ या सूची में अगले विषय पर जाने के लिए पृष्ठ के नीचे दाईं ओर बटन पर क्लिक करो।
