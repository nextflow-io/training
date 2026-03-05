# वर्कफ़्लो की डिबगिंग

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

डिबगिंग एक महत्वपूर्ण कौशल है जो आपको घंटों की निराशा से बचा सकता है और आपको एक अधिक प्रभावी Nextflow डेवलपर बनने में मदद कर सकता है। अपने करियर के दौरान, खासकर जब आप शुरुआत कर रहे हों, आपको अपने वर्कफ़्लो को बनाते और बनाए रखते समय बग्स का सामना करना पड़ेगा। व्यवस्थित डिबगिंग दृष्टिकोण सीखने से आपको समस्याओं की पहचान करने और उन्हें जल्दी हल करने में मदद मिलेगी।

### सीखने के लक्ष्य

इस साइड क्वेस्ट में, हम Nextflow वर्कफ़्लो के लिए **व्यवस्थित डिबगिंग तकनीकों** का अन्वेषण करेंगे:

- **सिंटैक्स एरर डिबगिंग**: IDE फीचर्स और Nextflow एरर मैसेज का प्रभावी ढंग से उपयोग करना
- **Channel डिबगिंग**: डेटा फ़्लो समस्याओं और channel संरचना समस्याओं का निदान करना
- **Process डिबगिंग**: निष्पादन विफलताओं और संसाधन समस्याओं की जांच करना
- **बिल्ट-इन डिबगिंग टूल्स**: Nextflow के preview mode, stub running, और work directories का लाभ उठाना
- **व्यवस्थित दृष्टिकोण**: कुशल डिबगिंग के लिए चार-चरण पद्धति

अंत तक, आपके पास एक मजबूत डिबगिंग पद्धति होगी जो निराशाजनक एरर मैसेज को समाधान के लिए स्पष्ट रोडमैप में बदल देती है।

### पूर्वापेक्षाएँ

इस साइड क्वेस्ट को शुरू करने से पहले, आपको:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष शुरुआती पाठ्यक्रम पूरा किया होना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्र (processes, channels, operators) का उपयोग करने में सहज होना चाहिए

**वैकल्पिक:** हम अनुशंसा करते हैं कि पहले [IDE Features for Nextflow Development](./ide_features.md) साइड क्वेस्ट को पूरा करें।
यह IDE फीचर्स की व्यापक कवरेज प्रदान करता है जो डिबगिंग का समर्थन करते हैं (सिंटैक्स हाइलाइटिंग, एरर डिटेक्शन, आदि), जिन्हें हम यहाँ भारी मात्रा में उपयोग करेंगे।

---

## 0. शुरुआत करें

#### प्रशिक्षण codespace खोलें

यदि आपने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित अनुसार प्रशिक्षण वातावरण खोलना सुनिश्चित करें।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाएं

आइए उस डायरेक्टरी में चलें जहाँ इस ट्यूटोरियल के लिए फ़ाइलें स्थित हैं।

```bash
cd side-quests/debugging
```

आप VSCode को इस डायरेक्टरी पर फोकस करने के लिए सेट कर सकते हैं:

```bash
code .
```

#### सामग्री की समीक्षा करें

आपको विभिन्न प्रकार के बग्स के साथ उदाहरण वर्कफ़्लो का एक सेट मिलेगा जिन्हें हम अभ्यास के लिए उपयोग करेंगे:

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

ये फ़ाइलें सामान्य डिबगिंग परिदृश्यों का प्रतिनिधित्व करती हैं जिनका सामना आप वास्तविक दुनिया के विकास में करेंगे।

#### असाइनमेंट की समीक्षा करें

आपकी चुनौती प्रत्येक वर्कफ़्लो को चलाना, एरर(एरर्स) की पहचान करना और उन्हें ठीक करना है।

प्रत्येक बगी वर्कफ़्लो के लिए:

1. **वर्कफ़्लो चलाएं** और एरर का अवलोकन करें
2. **एरर मैसेज का विश्लेषण करें**: Nextflow आपको क्या बता रहा है?
3. **कोड में समस्या का पता लगाएं** प्रदान किए गए सुरागों का उपयोग करके
4. **बग को ठीक करें** और सत्यापित करें कि आपका समाधान काम करता है
5. अगले अनुभाग पर जाने से पहले **फ़ाइल को रीसेट करें** (`git checkout <filename>` का उपयोग करें)

अभ्यास सरल सिंटैक्स एरर से अधिक सूक्ष्म रनटाइम समस्याओं तक प्रगति करते हैं।
समाधान इनलाइन चर्चा किए गए हैं, लेकिन आगे पढ़ने से पहले प्रत्येक को स्वयं हल करने का प्रयास करें।

#### तैयारी चेकलिस्ट

क्या आपको लगता है कि आप शुरू करने के लिए तैयार हैं?

- [ ] मैं इस पाठ्यक्रम के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता हूं
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी वर्किंग डायरेक्टरी उचित रूप से सेट की है
- [ ] मैं असाइनमेंट को समझता हूं

यदि आप सभी बॉक्स को चेक कर सकते हैं, तो आप जाने के लिए तैयार हैं।

---

## 1. सिंटैक्स एरर

सिंटैक्स एरर सबसे सामान्य प्रकार के एरर हैं जिनका सामना आप Nextflow कोड लिखते समय करेंगे। ये तब होते हैं जब कोड Nextflow DSL के अपेक्षित सिंटैक्स नियमों के अनुरूप नहीं होता है। ये एरर आपके वर्कफ़्लो को बिल्कुल चलने से रोकते हैं, इसलिए उन्हें जल्दी पहचानना और ठीक करना सीखना महत्वपूर्ण है।

### 1.1. गुम ब्रेसिज़

सबसे सामान्य सिंटैक्स एरर में से एक, और कभी-कभी डिबग करने के लिए अधिक जटिल में से एक **गुम या बेमेल ब्रैकेट** है।

आइए एक व्यावहारिक उदाहरण से शुरू करें।

#### पाइपलाइन चलाएं

```bash
nextflow run bad_syntax.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [stupefied_bhabha] DSL2 - revision: ca6327fad2

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

**सिंटैक्स एरर मैसेज के मुख्य तत्व:**

- **फ़ाइल और स्थान**: दिखाता है कि किस फ़ाइल और लाइन/कॉलम में एरर है (`bad_syntax.nf:24:1`)
- **एरर विवरण**: बताता है कि पार्सर को क्या मिला जो उसे उम्मीद नहीं थी (`Unexpected input: '<EOF>'`)
- **EOF संकेतक**: `<EOF>` (End Of File) मैसेज इंगित करता है कि पार्सर अभी भी अधिक सामग्री की अपेक्षा करते हुए फ़ाइल के अंत तक पहुंच गया - बंद न हुए ब्रेसिज़ का एक क्लासिक संकेत

#### कोड की जांच करें

अब, आइए `bad_syntax.nf` की जांच करें यह समझने के लिए कि एरर का कारण क्या है:

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
// process के लिए closing brace गुम है

workflow {

    // इनपुट channel बनाएं
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // इनपुट channel के साथ process को call करें
    PROCESS_FILES(input_ch)
}
```

इस उदाहरण के उद्देश्य के लिए हमने आपको दिखाने के लिए एक टिप्पणी छोड़ी है कि एरर कहाँ है। Nextflow VSCode एक्सटेंशन को भी आपको कुछ संकेत देना चाहिए कि क्या गलत हो सकता है, बेमेल ब्रेस को लाल रंग में रखते हुए और फ़ाइल के समय से पहले अंत को हाइलाइट करते हुए:

![Bad syntax](img/bad_syntax.png)

**ब्रैकेट एरर के लिए डिबगिंग रणनीति:**

1. VS Code की ब्रैकेट मैचिंग का उपयोग करें (कर्सर को ब्रैकेट के बगल में रखें)
2. ब्रैकेट-संबंधित मैसेज के लिए Problems पैनल की जांच करें
3. सुनिश्चित करें कि प्रत्येक opening `{` में एक संबंधित closing `}` है

#### कोड को ठीक करें

गुम closing brace के साथ टिप्पणी को बदलें:

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
    }  // गुम closing brace जोड़ें

    workflow {

        // इनपुट channel बनाएं
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // इनपुट channel के साथ process को call करें
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
    // process के लिए closing brace गुम है

    workflow {

        // इनपुट channel बनाएं
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // इनपुट channel के साथ process को call करें
        PROCESS_FILES(input_ch)
    }
    ```

#### पाइपलाइन चलाएं

अब वर्कफ़्लो को फिर से चलाएं यह पुष्टि करने के लिए कि यह काम करता है:

```bash
nextflow run bad_syntax.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [insane_faggin] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [48/cd7f54] PROCESS_FILES (1) | 3 of 3 ✔
    ```

### 1.2. गलत process कीवर्ड या निर्देशों का उपयोग करना

एक अन्य सामान्य सिंटैक्स एरर एक **अमान्य process परिभाषा** है। यह तब हो सकता है जब आप आवश्यक ब्लॉक परिभाषित करना भूल जाते हैं या process परिभाषा में गलत निर्देशों का उपयोग करते हैं।

#### पाइपलाइन चलाएं

```bash
nextflow run invalid_process.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

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

#### कोड की जांच करें

एरर एक "Invalid process definition" को इंगित करता है और समस्या के आसपास संदर्भ दिखाता है। लाइन 3-7 को देखते हुए, हम लाइन 4 पर `inputs:` देख सकते हैं, जो समस्या है। आइए `invalid_process.nf` की जांच करें:

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

    // इनपुट channel बनाएं
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // इनपुट channel के साथ process को call करें
    PROCESS_FILES(input_ch)
}
```

एरर संदर्भ में लाइन 4 को देखते हुए, हम समस्या को देख सकते हैं: हम सही `input` निर्देश के बजाय `inputs` का उपयोग कर रहे हैं। Nextflow VSCode एक्सटेंशन भी इसे फ्लैग करेगा:

![Invalid process message](img/invalid_process_message.png)

#### कोड को ठीक करें

[डॉक्यूमेंटेशन](https://www.nextflow.io/docs/latest/process.html#) का संदर्भ देकर गलत कीवर्ड को सही के साथ बदलें:

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

        // इनपुट channel बनाएं
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // इनपुट channel के साथ process को call करें
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

        // इनपुट channel बनाएं
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // इनपुट channel के साथ process को call करें
        PROCESS_FILES(input_ch)
    }
    ```

#### पाइपलाइन चलाएं

अब वर्कफ़्लो को फिर से चलाएं यह पुष्टि करने के लिए कि यह काम करता है:

```bash
nextflow run invalid_process.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [silly_fermi] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [b7/76cd9d] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.3. खराब वेरिएबल नामों का उपयोग करना

आपके स्क्रिप्ट ब्लॉक में उपयोग किए जाने वाले वेरिएबल नाम मान्य होने चाहिए, या तो इनपुट से या स्क्रिप्ट से पहले डाले गए groovy कोड से प्राप्त होने चाहिए। लेकिन जब आप पाइपलाइन विकास की शुरुआत में जटिलता से जूझ रहे हों, तो वेरिएबल नामकरण में गलतियां करना आसान है, और Nextflow आपको जल्दी से बता देगा।

#### पाइपलाइन चलाएं

```bash
nextflow run no_such_var.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [gloomy_meninsky] DSL2 - revision: 0c4d3bc28c

    Error no_such_var.nf:17:39: `undefined_var` is not defined
    │  17 |     echo "Using undefined variable: ${undefined_var}" >> ${output_pref
    ╰     |                                       ^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

एरर को कंपाइल समय पर पकड़ा जाता है और सीधे लाइन 17 पर अपरिभाषित वेरिएबल की ओर इशारा करता है, जिसमें एक कैरेट समस्या कहाँ है यह बिल्कुल इंगित करता है।

#### कोड की जांच करें

आइए `no_such_var.nf` की जांच करें:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // script से पहले Groovy code में variables परिभाषित करें
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

एरर मैसेज इंगित करता है कि वेरिएबल को स्क्रिप्ट टेम्पलेट में पहचाना नहीं गया है, और वहाँ आप जाते हैं- आपको स्क्रिप्ट ब्लॉक में उपयोग किया गया `${undefined_var}` देखने में सक्षम होना चाहिए, लेकिन कहीं और परिभाषित नहीं किया गया है।

#### कोड को ठीक करें

यदि आपको 'No such variable' एरर मिलती है, तो आप वेरिएबल को परिभाषित करके (इनपुट वेरिएबल नामों को सही करके या स्क्रिप्ट से पहले groovy कोड को संपादित करके), या यदि इसकी आवश्यकता नहीं है तो इसे स्क्रिप्ट ब्लॉक से हटाकर इसे ठीक कर सकते हैं:

=== "बाद में"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // script से पहले Groovy code में variables परिभाषित करें
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // undefined_var वाली लाइन हटाई गई
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
        // script से पहले Groovy code में variables परिभाषित करें
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

#### पाइपलाइन चलाएं

अब वर्कफ़्लो को फिर से चलाएं यह पुष्टि करने के लिए कि यह काम करता है:

```bash
nextflow run no_such_var.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [suspicious_venter] DSL2 - revision: 6ba490f7c5

    executor >  local (3)
    [21/237300] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.4. Bash वेरिएबल का खराब उपयोग

Nextflow में शुरुआत करते समय, Nextflow (Groovy) और Bash वेरिएबल के बीच अंतर को समझना मुश्किल हो सकता है। यह एक अन्य रूप की खराब वेरिएबल एरर उत्पन्न कर सकता है जो स्क्रिप्ट ब्लॉक की Bash सामग्री में वेरिएबल का उपयोग करने का प्रयास करते समय दिखाई देती है।

#### पाइपलाइन चलाएं

```bash
nextflow run bad_bash_var.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [infallible_mandelbrot] DSL2 - revision: 0853c11080

    Error bad_bash_var.nf:13:42: `prefix` is not defined
    │  13 |     echo "Processing ${sample_name}" > ${prefix}.txt
    ╰     |                                          ^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### कोड की जांच करें

एरर लाइन 13 की ओर इशारा करता है जहाँ `${prefix}` का उपयोग किया गया है। आइए `bad_bash_var.nf` की जांच करें यह देखने के लिए कि समस्या का कारण क्या है:

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
    echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
    """
}
```

इस उदाहरण में, हम Bash में `prefix` वेरिएबल को परिभाषित कर रहे हैं, लेकिन एक Nextflow process में हमने इसे संदर्भित करने के लिए उपयोग किया गया `$` सिंटैक्स (`${prefix}`) को Groovy वेरिएबल के रूप में व्याख्या किया जाता है, न कि Bash। वेरिएबल Groovy संदर्भ में मौजूद नहीं है, इसलिए हमें 'no such variable' एरर मिलता है।

#### कोड को ठीक करें

यदि आप एक Bash वेरिएबल का उपयोग करना चाहते हैं, तो आपको डॉलर चिह्न को इस तरह एस्केप करना होगा:

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
        echo "Processing ${sample_name}" > \${prefix}.txt  # Fixed: Escaped the dollar sign
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
        echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
        """
    }
    ```

यह Nextflow को बताता है कि इसे Bash वेरिएबल के रूप में व्याख्या करें।

#### पाइपलाइन चलाएं

अब वर्कफ़्लो को फिर से चलाएं यह पुष्टि करने के लिए कि यह काम करता है:

```bash
nextflow run bad_bash_var.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [naughty_franklin] DSL2 - revision: 58c1c83709

    executor >  local (3)
    [4e/560285] PROCESS_FILES (2) | 3 of 3 ✔
    ```

!!! tip "Groovy बनाम Bash वेरिएबल"

    स्ट्रिंग संयोजन या prefix/suffix ऑपरेशन जैसे सरल वेरिएबल मैनिपुलेशन के लिए, स्क्रिप्ट ब्लॉक में Bash वेरिएबल के बजाय script सेक्शन में Groovy वेरिएबल का उपयोग करना आमतौर पर अधिक पठनीय होता है:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    यह दृष्टिकोण डॉलर चिह्न को एस्केप करने की आवश्यकता से बचता है और कोड को पढ़ने और बनाए रखने में आसान बनाता है।

### 1.5. Workflow ब्लॉक के बाहर स्टेटमेंट

Nextflow VSCode एक्सटेंशन कोड संरचना के साथ समस्याओं को हाइलाइट करता है जो एरर का कारण बनेंगे। एक सामान्य उदाहरण `workflow {}` ब्लॉक के बाहर channels को परिभाषित करना है - यह अब एक सिंटैक्स एरर के रूप में लागू किया गया है।

#### पाइपलाइन चलाएं

```bash
nextflow run badpractice_syntax.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [intergalactic_colden] DSL2 - revision: 5e4b291bde

    Error badpractice_syntax.nf:3:1: Statements cannot be mixed with script declarations -- move statements into a process or workflow
    │   3 | input_ch = channel.of('sample1', 'sample2', 'sample3')
    ╰     | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

एरर मैसेज स्पष्ट रूप से समस्या को इंगित करता है: स्टेटमेंट (जैसे channel परिभाषाएं) को workflow या process ब्लॉक के बाहर स्क्रिप्ट घोषणाओं के साथ मिश्रित नहीं किया जा सकता है।

#### कोड की जांच करें

आइए `badpractice_syntax.nf` की जांच करें यह देखने के लिए कि एरर का कारण क्या है:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // एरर: Channel workflow के बाहर परिभाषित है

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // script से पहले Groovy code में variables परिभाषित करें
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

VSCode एक्सटेंशन `input_ch` वेरिएबल को workflow ब्लॉक के बाहर परिभाषित होने के रूप में भी हाइलाइट करेगा:

![Non-lethal syntax error](img/nonlethal.png)

#### कोड को ठीक करें

channel परिभाषा को workflow ब्लॉक के अंदर ले जाएं:

=== "बाद में"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // script से पहले Groovy code में variables परिभाषित करें
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
    #!/usr्bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // एरर: Channel workflow के बाहर परिभाषित है

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // script से पहले Groovy code में variables परिभाषित करें
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

#### पाइपलाइन चलाएं

फिक्स काम करती है यह पुष्टि करने के लिए वर्कफ़्लो को फिर से चलाएं:

```bash
nextflow run badpractice_syntax.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [naughty_ochoa] DSL2 - revision: 5e4b291bde

    executor >  local (3)
    [6a/84a608] PROCESS_FILES (2) | 3 of 3 ✔
    ```

अपने इनपुट channels को workflow ब्लॉक के भीतर परिभाषित रखें, और सामान्य रूप से एक्सटेंशन द्वारा की गई किसी भी अन्य सिफारिश का पालन करें।

### निष्कर्ष

आप Nextflow एरर मैसेज और IDE दृश्य संकेतकों का उपयोग करके सिंटैक्स एरर को व्यवस्थित रूप से पहचान और ठीक कर सकते हैं। सामान्य सिंटैक्स एरर में गुम ब्रेसिज़, गलत process कीवर्ड, अपरिभाषित वेरिएबल, और Bash बनाम Nextflow वेरिएबल का अनुचित उपयोग शामिल हैं। VSCode एक्सटेंशन रनटाइम से पहले इनमें से कई को पकड़ने में मदद करता है। आपके टूलकिट में इन सिंटैक्स डिबगिंग कौशल के साथ, आप सबसे सामान्य Nextflow सिंटैक्स एरर को जल्दी हल करने में सक्षम होंगे और अधिक जटिल रनटाइम समस्याओं से निपटने के लिए आगे बढ़ेंगे।

### आगे क्या है?

अधिक जटिल channel संरचना एरर को डिबग करना सीखें जो तब भी होते हैं जब सिंटैक्स सही है।

---

## 2. Channel संरचना एरर

Channel संरचना एरर सिंटैक्स एरर से अधिक सूक्ष्म हैं क्योंकि कोड सिंटैक्टिकली सही है, लेकिन डेटा आकार उस से मेल नहीं खाते जो processes अपेक्षा करते हैं। Nextflow पाइपलाइन चलाने का प्रयास करेगा, लेकिन यह पा सकता है कि इनपुट की संख्या उसकी अपेक्षा से मेल नहीं खाती और विफल हो जाती है। ये एरर आम तौर पर केवल रनटाइम पर दिखाई देते हैं और आपके वर्कफ़्लो के माध्यम से प्रवाहित होने वाले डेटा की समझ की आवश्यकता होती है।

!!! tip "`.view()` के साथ Channels की डिबगिंग"

    इस सेक्शन के दौरान, याद रखें कि आप अपने वर्कफ़्लो में किसी भी बिंदु पर channel सामग्री का निरीक्षण करने के लिए `.view()` operator का उपयोग कर सकते हैं। यह channel संरचना समस्याओं को समझने के लिए सबसे शक्तिशाली डिबगिंग टूल में से एक है। हम सेक्शन 2.4 में इस तकनीक का विस्तार से पता लगाएंगे, लेकिन उदाहरणों के माध्यम से काम करते समय इसका उपयोग करने के लिए स्वतंत्र महसूस करें।

    ```groovy
    my_channel.view()  // channel में क्या बह रहा है दिखाता है
    ```

### 2.1. गलत संख्या में इनपुट Channels

यह एरर तब होती है जब आप एक process की अपेक्षा से अलग संख्या में channels पास करते हैं।

#### पाइपलाइन चलाएं

```bash
nextflow run bad_number_inputs.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [happy_swartz] DSL2 - revision: d83e58dcd3

    Error bad_number_inputs.nf:23:5: Incorrect number of call arguments, expected 1 but received 2
    │  23 |     PROCESS_FILES(samples_ch, files_ch)
    ╰     |     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### कोड की जांच करें

एरर मैसेज स्पष्ट रूप से कहता है कि कॉल को 1 argument की उम्मीद थी लेकिन 2 प्राप्त हुए, और लाइन 23 की ओर इशारा करता है। आइए `bad_number_inputs.nf` की जांच करें:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Process केवल 1 इनपुट की अपेक्षा करती है

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // दो अलग-अलग channels बनाएं
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: 2 channels पास कर रहे हैं लेकिन process केवल 1 की अपेक्षा करता है
    PROCESS_FILES(samples_ch, files_ch)
}
```

आपको बेमेल `PROCESS_FILES` कॉल देखनी चाहिए, जब process केवल एक को परिभाषित करती है तो कई इनपुट channels की आपूर्ति करना। VSCode एक्सटेंशन भी process कॉल को लाल रंग में रेखांकित करेगा, और जब आप माउस ओवर करते हैं तो एक diagnostic मैसेज प्रदान करेगा:

![Incorrect number of args message](img/incorrect_num_args.png)

#### कोड को ठीक करें

इस विशिष्ट उदाहरण के लिए, process एक single channel की उम्मीद करती है और दूसरे channel की आवश्यकता नहीं है, इसलिए हम केवल `samples_ch` channel को पास करके इसे ठीक कर सकते हैं:

=== "बाद में"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process केवल 1 इनपुट की अपेक्षा करती है

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // दो अलग-अलग channels बनाएं
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ठीक किया: केवल वह channel पास करें जो process अपेक्षा करता है
        PROCESS_FILES(samples_ch)
    }
    ```

=== "पहले"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process केवल 1 इनपुट की अपेक्षा करती है

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // दो अलग-अलग channels बनाएं
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERROR: 2 channels पास कर रहे हैं लेकिन process केवल 1 की अपेक्षा करता है
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### पाइपलाइन चलाएं

```bash
nextflow run bad_number_inputs.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [big_euler] DSL2 - revision: e302bd87be

    executor >  local (3)
    [48/497f7b] PROCESS_FILES (3) | 3 of 3 ✔
    ```

इस उदाहरण की तुलना में अधिक सामान्यतः, आप एक process में अतिरिक्त इनपुट जोड़ सकते हैं और तदनुसार workflow कॉल को अपडेट करना भूल सकते हैं, जो इस प्रकार की एरर का कारण बन सकता है। सौभाग्य से, यह समझने और ठीक करने में आसान एरर में से एक है, क्योंकि एरर मैसेज बेमेल के बारे में काफी स्पष्ट है।

### 2.2. Channel समाप्ति (Process अपेक्षा से कम बार चलती है)

कुछ channel संरचना एरर बहुत अधिक सूक्ष्म हैं और बिल्कुल भी एरर उत्पन्न नहीं करते हैं। शायद इनमें से सबसे सामान्य एक चुनौती को दर्शाता है जिसका सामना नए Nextflow उपयोगकर्ताओं को यह समझने में करना पड़ता है कि queue channels समाप्त हो सकते हैं और आइटम खत्म हो सकते हैं, जिसका अर्थ है कि वर्कफ़्लो समय से पहले समाप्त हो जाता है।

#### पाइपलाइन चलाएं

```bash
nextflow run exhausted.nf
```

??? success "कमांड आउटपुट"

```console title="Exhausted channel output"
 N E X T F L O W   ~  version 25.10.2

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

यह वर्कफ़्लो बिना एरर के पूरा होता है, लेकिन यह केवल एक single नमूने को प्रोसेस करता है!

#### कोड की जांच करें

आइए `exhausted.nf` की जांच करें यह देखने के लिए कि क्या यह सही है:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // script से पहले Groovy code में variables परिभाषित करें
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

Process केवल एक बार चलती है तीन बार के बजाय क्योंकि `reference_ch` channel एक queue channel है जो पहले process निष्पादन के बाद समाप्त हो जाती है। जब एक channel समाप्त हो जाती है, तो पूरी process रुक जाती है, भले ही अन्य channels में अभी भी आइटम हों।

यह एक सामान्य पैटर्न है जहां आपके पास एक single संदर्भ फ़ाइल है जिसे कई नमूनों में पुन: उपयोग करने की आवश्यकता है। समाधान संदर्भ channel को एक value channel में परिवर्तित करना है जिसे अनिश्चित काल तक पुन: उपयोग किया जा सकता है।

#### कोड को ठीक करें

इसे संबोधित करने के कुछ तरीके हैं जो इस बात पर निर्भर करते हैं कि कितनी फ़ाइलें प्रभावित हैं।

**विकल्प 1**: आपके पास एक single संदर्भ फ़ाइल है जिसे आप बहुत पुन: उपयोग कर रहे हैं। आप बस एक value channel प्रकार बना सकते हैं, जिसका उपयोग बार-बार किया जा सकता है। ऐसा करने के तीन तरीके हैं:

**1a** `channel.value()` का उपयोग करें:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Value channel को पुन: उपयोग किया जा सकता है
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** `first()` [operator](https://www.nextflow.io/docs/latest/reference/operator.html#first) का उपयोग करें:

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Value channel में बदलें
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** `collect()` [operator](https://www.nextflow.io/docs/latest/reference/operator.html#collect) का उपयोग करें:

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Value channel में बदलें
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**विकल्प 2**: अधिक जटिल परिदृश्यों में, शायद जहां आपके पास sample channel में सभी नमूनों के लिए कई संदर्भ फ़ाइलें हैं, आप `combine` operator का उपयोग कर सकते हैं एक नया channel बनाने के लिए जो दोनों channels को tuples में संयोजित करता है:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // कार्टेशियन उत्पाद बनाता है

    PROCESS_FILES(combined_ch)
}
```

`.combine()` operator दोनों channels का एक कार्टेशियन उत्पाद उत्पन्न करता है, इसलिए `reference_ch` में प्रत्येक आइटम को `input_ch` में प्रत्येक आइटम के साथ जोड़ा जाएगा। यह process को संदर्भ का उपयोग करते हुए प्रत्येक नमूने के लिए चलाने की अनुमति देता है।

इसके लिए process इनपुट को समायोजित करने की आवश्यकता है। हमारे उदाहरण में, process परिभाषा की शुरुआत को निम्नानुसार समायोजित करने की आवश्यकता होगी:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

यह दृष्टिकोण सभी स्थितियों में उपयुक्त नहीं हो सकता है।

#### पाइपलाइन चलाएं

ऊपर दिए गए फिक्स में से एक को आजमाएं और वर्कफ़्लो को फिर से चलाएं:

```bash
nextflow run exhausted.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `exhausted.nf` [maniac_leavitt] DSL2 - revision: f372a56a7d

    executor >  local (3)
    [80/0779e9] PROCESS_FILES (3) | 3 of 3 ✔
    ```

आपको अब केवल एक के बजाय तीनों नमूनों को प्रोसेस किए जाते हुए देखना चाहिए।

### 2.3. गलत Channel सामग्री संरचना

जब वर्कफ़्लो जटिलता के एक निश्चित स्तर तक पहुंचते हैं, तो प्रत्येक channel की आंतरिक संरचनाओं पर नज़र रखना थोड़ा मुश्किल हो सकता है, और लोग आमतौर पर process की अपेक्षा और channel में वास्तव में क्या है के बीच बेमेल उत्पन्न करते हैं। यह उस समस्या से अधिक सूक्ष्म है जिस पर हमने पहले चर्चा की, जहां channels की संख्या गलत थी। इस मामले में, आपके पास सही संख्या में इनपुट channels हो सकते हैं, लेकिन उन channels में से एक या अधिक की आंतरिक संरचना उस से मेल नहीं खाती जो process अपेक्षा करती है।

#### पाइपलाइन चलाएं

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

#### कोड की जांच करें

एरर मैसेज में वर्ग कोष्ठक यहाँ सुराग प्रदान करते हैं - process tuple को एकल मान के रूप में मान रहा है, जो हम नहीं चाहते। आइए `bad_channel_shape.nf` की जांच करें:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Expects single value, gets tuple

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Channel tuples emit करता है, लेकिन process single values की अपेक्षा करता है
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

आप देख सकते हैं कि हम tuples से बना एक channel उत्पन्न कर रहे हैं: `['sample1', 'file1.txt']`, लेकिन process एकल मान की अपेक्षा करता है, `val sample_name`। निष्पादित कमांड दिखाता है कि process `[sample3, file3.txt]_output.txt` नाम की एक फ़ाइल बनाने की कोशिश कर रहा है, जो अभीष्ट आउटपुट नहीं है।

#### कोड को ठीक करें

इसे ठीक करने के लिए, यदि process को दोनों इनपुट की आवश्यकता है तो हम process को एक tuple स्वीकार करने के लिए समायोजित कर सकते हैं:

=== "विकल्प 1: Process में tuple स्वीकार करें"

    === "बाद में"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Fixed: Accept tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel tuples emit करता है, लेकिन process single values की अपेक्षा करता है
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
                val sample_name  // Expects single value, gets tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel tuples emit करता है, लेकिन process single values की अपेक्षा करता है
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "विकल्प 2: पहला तत्व निकालें"

    === "बाद में"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel tuples emit करता है, लेकिन process single values की अपेक्षा करता है
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Fixed: Extract first element
        }
        ```

    === "पहले"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel tuples emit करता है, लेकिन process single values की अपेक्षा करता है
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### पाइपलाइन चलाएं

कोई एक समाधान चुनें और वर्कफ़्लो को फिर से चलाएं:

```bash
nextflow run bad_channel_shape.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape.nf` [clever_thompson] DSL2 - revision: 8cbcae3746

    executor >  local (3)
    [bb/80a958] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 2.4. Channel डिबगिंग तकनीकें

#### Channel निरीक्षण के लिए `.view()` का उपयोग

Channels के लिए सबसे शक्तिशाली डिबगिंग टूल `.view()` operator है। `.view()` के साथ, आप डिबगिंग में सहायता के लिए सभी चरणों में अपने channels के आकार को समझ सकते हैं।

#### पाइपलाइन चलाएं

इसे क्रिया में देखने के लिए `bad_channel_shape_viewed.nf` चलाएं:

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

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

#### कोड की जांच करें

आइए `bad_channel_shape_viewed.nf` की जांच करें यह देखने के लिए कि `.view()` का उपयोग कैसे किया जाता है:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Channel tuples emit करता है, लेकिन process single values की अपेक्षा करता है
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debug: Show original channel content
    .map { tuple -> tuple[0] }        // Transform: Extract first element
    .view { "After mapping: $it" }    // Debug: Show transformed channel content

    PROCESS_FILES(input_ch)
}
```

#### कोड को ठीक करें

भविष्य में channel सामग्री को समझने के लिए `.view()` ऑपरेशन का अत्यधिक उपयोग करने से बचने के लिए, कुछ टिप्पणियाँ जोड़ना उचित है:

```groovy title="bad_channel_shape_viewed.nf (टिप्पणियों के साथ)" linenums="16" hl_lines="8 9"
workflow {

    // Channel tuples emit करता है, लेकिन process single values की अपेक्षा करता है
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

यह आपके वर्कफ़्लो के जटिलता में बढ़ने और channel संरचना के अधिक अपारदर्शी होने पर अधिक महत्वपूर्ण हो जाएगा।

#### पाइपलाइन चलाएं

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

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

### निष्कर्ष

कई channel संरचना एरर मान्य Nextflow सिंटैक्स के साथ बनाई जा सकती हैं। आप डेटा फ्लो को समझकर, निरीक्षण के लिए `.view()` operators का उपयोग करके, और अप्रत्याशित tuple संरचनाओं को इंगित करने वाले वर्ग कोष्ठक जैसे एरर मैसेज पैटर्न को पहचानकर channel संरचना एरर को डिबग कर सकते हैं।

### आगे क्या है?

Process परिभाषाओं द्वारा बनाई गई एरर के बारे में जानें।

---

## 3. Process संरचना एरर

Processes से संबंधित अधिकांश एरर जो आपको मिलेंगे वे कमांड बनाने में आपकी गलतियों से, या अंतर्निहित सॉफ़्टवेयर से संबंधित समस्याओं से संबंधित होंगी। कहा जाता है, ऊपर दी गई channel समस्याओं के समान, आप process परिभाषा में ऐसी गलतियाँ कर सकते हैं जो सिंटैक्स एरर के रूप में योग्य नहीं हैं, लेकिन जो रन टाइम पर एरर का कारण बनेंगी।

### 3.1. गुम आउटपुट फ़ाइलें

Processes लिखते समय एक सामान्य एरर ऐसा कुछ करना है जो process की अपेक्षा और जो उत्पन्न होता है उसके बीच बेमेल पैदा करता है।

#### पाइपलाइन चलाएं

```bash
nextflow run missing_output.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

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

#### कोड की जांच करें

एरर मैसेज इंगित करता है कि process ने `sample3.txt` नाम की एक आउटपुट फ़ाइल उत्पन्न करने की अपेक्षा की थी, लेकिन स्क्रिप्ट वास्तव में `sample3_output.txt` बनाती है। आइए `missing_output.nf` में process परिभाषा की जांच करें:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Expects: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
    """
}
```

आपको देखना चाहिए कि `output:` ब्लॉक में आउटपुट फ़ाइल नाम और स्क्रिप्ट में उपयोग किए गए नाम के बीच बेमेल है। यह बेमेल process को विफल करता है। यदि आपको इस प्रकार की एरर मिलती है, तो वापस जाएं और जांचें कि आउटपुट आपकी process परिभाषा और आपके output ब्लॉक के बीच मेल खाते हैं।

यदि समस्या अभी भी स्पष्ट नहीं है, तो वास्तव में बनाई गई आउटपुट फ़ाइलों की पहचान करने के लिए work डायरेक्टरी की जांच करें:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

इस उदाहरण के लिए यह हमें हाइलाइट करेगा कि एक `_output` suffix आउटपुट फ़ाइल नाम में शामिल किया जा रहा है, जो हमारी `output:` परिभाषा के विपरीत है।

#### कोड को ठीक करें

आउटपुट फ़ाइलनाम को सुसंगत बनाकर बेमेल को ठीक करें:

=== "बाद में"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Fixed: Match the script output

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
        path "${sample_name}.txt"  // Expects: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
        """
    }
    ```

#### पाइपलाइन चलाएं

```bash
nextflow run missing_output.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [elated_hamilton] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [16/1c437c] PROCESS_FILES (3) | 3 of 3 ✔
    ```

### 3.2. गुम सॉफ़्टवेयर

एरर का एक अन्य वर्ग सॉफ़्टवेयर प्रावधान में गलतियों के कारण होता है। `missing_software.nf` एक सिंटैक्टिकली मान्य वर्कफ़्लो है, लेकिन यह `cowpy` कमांड प्रदान करने के लिए कुछ बाहरी सॉफ़्टवेयर पर निर्भर करता है।

#### पाइपलाइन चलाएं

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

Process के पास हमारे द्वारा निर्दिष्ट कमांड तक पहुंच नहीं है। कभी-कभी ऐसा इसलिए होता है क्योंकि एक स्क्रिप्ट वर्कफ़्लो `bin` डायरेक्टरी में मौजूद है, लेकिन इसे executable नहीं बनाया गया है। अन्य समय में यह इसलिए होता है क्योंकि सॉफ़्टवेयर container या उस वातावरण में स्थापित नहीं है जहां वर्कफ़्लो चल रहा है।

#### कोड की जांच करें

उस `127` एग्जिट कोड पर ध्यान दें - यह आपको बिल्कुल समस्या बताता है। आइए `missing_software.nf` की जांच करें:

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

#### कोड को ठीक करें

हम यहाँ थोड़ा बेईमान रहे हैं, और वास्तव में कोड में कुछ भी गलत नहीं है। हमें बस process को इस तरह चलाने के लिए आवश्यक कॉन्फ़िगरेशन निर्दिष्ट करने की आवश्यकता है कि उसके पास प्रश्न में कमांड तक पहुंच हो। इस मामले में process में एक container परिभाषा है, इसलिए हमें बस Docker सक्षम करके वर्कफ़्लो चलाना है।

#### पाइपलाइन चलाएं

हमने आपके लिए `nextflow.config` में एक Docker प्रोफ़ाइल सेट किया है, इसलिए आप वर्कफ़्लो को इसके साथ चला सकते हैं:

```bash
nextflow run missing_software.nf -profile docker
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software.nf` [awesome_stonebraker] DSL2 - revision: 0296d12839

    executor >  local (3)
    [38/ab20d1] PROCESS_FILES (1) | 3 of 3 ✔
    ```

!!! note

    Nextflow containers कैसे उपयोग करता है इसके बारे में अधिक जानने के लिए, [Hello Nextflow](../hello_nextflow/05_hello_containers.md) देखें

### 3.3. खराब संसाधन कॉन्फ़िगरेशन

उत्पादन उपयोग में, आप अपनी processes पर संसाधनों को कॉन्फ़िगर करेंगे। उदाहरण के लिए `memory` आपकी process के लिए उपलब्ध अधिकतम मेमोरी मात्रा को परिभाषित करता है, और यदि process उससे अधिक हो जाती है, तो आपका scheduler आमतौर पर process को समाप्त कर देगा और `137` का एग्जिट कोड लौटाएगा। हम यहाँ इसका प्रदर्शन नहीं कर सकते क्योंकि हम `local` executor का उपयोग कर रहे हैं, लेकिन हम `time` के साथ कुछ ऐसा ही दिखा सकते हैं।

#### पाइपलाइन चलाएं

`bad_resources.nf` में 1 मिलीसेकंड की अवास्तविक सीमा के साथ process कॉन्फ़िगरेशन है:

```bash
nextflow run bad_resources.nf -profile docker
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

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

#### कोड की जांच करें

आइए `bad_resources.nf` की जांच करें:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // ERROR: Unrealistic time limit

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // Takes 1 second, but time limit is 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

हम जानते हैं कि process एक सेकंड से अधिक समय लेगा (हमने इसे सुनिश्चित करने के लिए एक sleep जोड़ा है), लेकिन process 1 मिलीसेकंड के बाद समय समाप्त करने के लिए सेट है। किसी ने अपने कॉन्फ़िगरेशन के साथ थोड़ा अवास्तविक रहा है!

#### कोड को ठीक करें

समय सीमा को एक वास्तविक मान तक बढ़ाएं:

=== "बाद में"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // Fixed: Realistic time limit

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

        time '1 ms'  // ERROR: Unrealistic time limit

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // Takes 1 second, but time limit is 1ms
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### पाइपलाइन चलाएं

```bash
nextflow run bad_resources.nf -profile docker
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [friendly_mcclintock] DSL2 - revision: 381567d2c1

    executor >  local (3)
    [c2/9b4c41] PROCESS_FILES (3) | 3 of 3 ✔
    ```

यदि आप अपने एरर मैसेज पढ़ने का ध्यान रखें तो इस प्रकार की विफलताओं से आप बहुत लंबे समय तक भ्रमित नहीं रहेंगे। लेकिन सुनिश्चित करें कि आप चलाए जा रहे कमांड की संसाधन आवश्यकताओं को समझते हैं ताकि आप अपने संसाधन निर्देशों को उचित रूप से कॉन्फ़िगर कर सकें।

### 3.4. Process डिबगिंग तकनीकें

जब processes विफल होते हैं या अप्रत्याशित रूप से व्यवहार करते हैं, तो आपको यह जांचने के लिए व्यवस्थित तकनीकों की आवश्यकता होती है कि क्या गलत हुआ। Work डायरेक्टरी में process निष्पादन को डिबग करने के लिए आवश्यक सभी जानकारी होती है।

#### Work डायरेक्टरी निरीक्षण का उपयोग

Processes के लिए सबसे शक्तिशाली डिबगिंग टूल work डायरेक्टरी की जांच करना है। जब एक process विफल होती है, तो Nextflow उस विशिष्ट process निष्पादन के लिए एक work डायरेक्टरी बनाता है जिसमें क्या हुआ यह समझने के लिए आवश्यक सभी फ़ाइलें होती हैं।

#### पाइपलाइन चलाएं

Work डायरेक्टरी निरीक्षण प्रदर्शित करने के लिए पहले के `missing_output.nf` उदाहरण का उपयोग करें (यदि आवश्यक हो तो आउटपुट नामकरण बेमेल को फिर से उत्पन्न करें):

```bash
nextflow run missing_output.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

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

#### Work डायरेक्टरी की जांच करें

जब आपको यह एरर मिलती है, तो work डायरेक्टरी में सभी डिबगिंग जानकारी होती है। एरर मैसेज से work डायरेक्टरी पथ खोजें और इसकी सामग्री की जांच करें:

```bash
# एरर मैसेज से work डायरेक्टरी खोजें
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

फिर आप मुख्य फ़ाइलों की जांच कर सकते हैं:

##### कमांड स्क्रिप्ट की जांच करें

`.command.sh` फ़ाइल बिल्कुल दिखाती है कि कौन सा कमांड निष्पादित किया गया था:

```bash
# निष्पादित कमांड देखें
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

यह प्रकट करता है:

- **वेरिएबल प्रतिस्थापन**: क्या Nextflow वेरिएबल ठीक से विस्तारित किए गए थे
- **फ़ाइल पथ**: क्या इनपुट फ़ाइलें सही ढंग से स्थित थीं
- **कमांड संरचना**: क्या स्क्रिप्ट सिंटैक्स सही है

देखने योग्य सामान्य समस्याएं:

- **गुम कोट्स**: स्पेस वाले वेरिएबल को उचित कोटिंग की आवश्यकता है
- **गलत फ़ाइल पथ**: इनपुट फ़ाइलें जो मौजूद नहीं हैं या गलत स्थानों पर हैं
- **गलत वेरिएबल नाम**: वेरिएबल संदर्भों में टाइपो
- **गुम पर्यावरण सेटअप**: कमांड जो विशिष्ट वातावरणों पर निर्भर करते हैं

##### एरर आउटपुट की जांच करें

`.command.err` फ़ाइल में वास्तविक एरर मैसेज होते हैं:

```bash
# एरर आउटपुट देखें
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

यह फ़ाइल दिखाएगी:

- **एग्जिट कोड**: 127 (कमांड नहीं मिला), 137 (समाप्त किया गया), आदि
- **अनुमति एरर**: फ़ाइल एक्सेस समस्याएं
- **सॉफ़्टवेयर एरर**: एप्लिकेशन-विशिष्ट एरर मैसेज
- **संसाधन एरर**: मेमोरी/समय सीमा पार हो गई

##### स्टैंडर्ड आउटपुट की जांच करें

`.command.out` फ़ाइल दिखाती है कि आपके कमांड ने क्या उत्पन्न किया:

```bash
# स्टैंडर्ड आउटपुट देखें
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

यह सत्यापित करने में मदद करता है:

- **अपेक्षित आउटपुट**: क्या कमांड ने सही परिणाम उत्पन्न किए
- **आंशिक निष्पादन**: क्या कमांड शुरू हुआ लेकिन बीच में विफल हो गया
- **डिबग जानकारी**: आपकी स्क्रिप्ट से कोई भी निदान आउटपुट

##### एग्जिट कोड की जांच करें

`.exitcode` फ़ाइल में process का एग्जिट कोड होता है:

```bash
# एग्जिट कोड देखें
cat work/*/*/.exitcode
```

सामान्य एग्जिट कोड और उनके अर्थ:

- **एग्जिट कोड 127**: कमांड नहीं मिला - सॉफ़्टवेयर इंस्टॉलेशन की जांच करें
- **एग्जिट कोड 137**: Process समाप्त किया गया - मेमोरी/समय सीमा की जांच करें

##### फ़ाइल अस्तित्व की जांच करें

जब गुम आउटपुट फ़ाइलों के कारण processes विफल होती हैं, तो जांचें कि वास्तव में कौन सी फ़ाइलें बनाई गई थीं:

```bash
# work डायरेक्टरी में सभी फ़ाइलों की सूची
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

यह पहचानने में मदद करता है:

- **फ़ाइल नामकरण बेमेल**: अपेक्षित से अलग नामों वाली आउटपुट फ़ाइलें
- **अनुमति समस्याएं**: फ़ाइलें जो बनाई नहीं जा सकीं
- **पथ समस्याएं**: गलत डायरेक्टरी में बनाई गई फ़ाइलें

पहले के हमारे उदाहरण में, इसने हमें पुष्टि की कि जबकि हमारा अपेक्षित `sample3.txt` मौजूद नहीं था, `sample3_output.txt` था:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### निष्कर्ष

Process डिबगिंग के लिए work डायरेक्टरी की जांच करना आवश्यक है यह समझने के लिए कि क्या गलत हुआ। मुख्य फ़ाइलों में `.command.sh` (निष्पादित स्क्रिप्ट), `.command.err` (एरर मैसेज), और `.command.out` (स्टैंडर्ड आउटपुट) शामिल हैं। एग्जिट कोड जैसे 127 (कमांड नहीं मिला) और 137 (process समाप्त किया गया) विफलता के प्रकार के बारे में तत्काल निदान सुराग प्रदान करते हैं।

### आगे क्या है?

Nextflow के बिल्ट-इन डिबगिंग टूल्स और समस्या निवारण के व्यवस्थित दृष्टिकोणों के बारे में जानें।

---

## 4. बिल्ट-इन डिबगिंग टूल्स और उन्नत तकनीकें

Nextflow वर्कफ़्लो निष्पादन को डिबग और विश्लेषण करने के लिए कई शक्तिशाली बिल्ट-इन टूल्स प्रदान करता है। ये टूल्स आपको समझने में मदद करते हैं कि क्या गलत हुआ, कहाँ गलत हुआ, और इसे कुशलतापूर्वक कैसे ठीक किया जाए।

### 4.1. रियल-टाइम Process आउटपुट

कभी-कभी आपको यह देखने की आवश्यकता होती है कि चल रही processes के अंदर क्या हो रहा है। आप रियल-टाइम process आउटपुट सक्षम कर सकते हैं, जो आपको बिल्कुल दिखाता है कि प्रत्येक कार्य निष्पादित होते समय क्या कर रहा है।

#### पाइपलाइन चलाएं

हमारे पहले के उदाहरणों से `bad_channel_shape_viewed.nf` ने `.view()` का उपयोग करके channel सामग्री प्रिंट की, लेकिन हम process के भीतर से ही वेरिएबल को echo करने के लिए `debug` निर्देश का भी उपयोग कर सकते हैं, जिसे हम `bad_channel_shape_viewed_debug.nf` में प्रदर्शित करते हैं। वर्कफ़्लो चलाएं:

```bash
nextflow run bad_channel_shape_viewed_debug.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

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

#### कोड की जांच करें

आइए `bad_channel_shape_viewed_debug.nf` की जांच करें यह देखने के लिए कि `debug` निर्देश कैसे काम करता है:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Enable real-time output

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

`debug` निर्देश एक process के वातावरण को समझने का एक त्वरित और सुविधाजनक तरीका हो सकता है।

### 4.2. Preview Mode

कभी-कभी आप किसी भी process के चलने से पहले समस्याओं को पकड़ना चाहते हैं। Nextflow इस प्रकार की सक्रिय डिबगिंग के लिए एक flag प्रदान करता है: `-preview`।

#### पाइपलाइन चलाएं

Preview mode आपको बिना कमांड निष्पादित किए वर्कफ़्लो logic का परीक्षण करने देता है। यह आपके वर्कफ़्लो की संरचना की जल्दी जांच करने और यह सुनिश्चित करने के लिए काफी उपयोगी हो सकता है कि processes बिना कोई वास्तविक कमांड चलाए सही ढंग से जुड़े हुए हैं।

!!! note

    यदि आपने पहले `bad_syntax.nf` को ठीक किया था, तो इस कमांड को चलाने से पहले script ब्लॉक के बाद closing brace हटाकर सिंटैक्स एरर को फिर से शुरू करें।

यह कमांड चलाएं:

```bash
nextflow run bad_syntax.nf -preview
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [magical_mercator] DSL2 - revision: 550b9a8873

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Preview mode किसी भी process को चलाए बिना सिंटैक्स एरर को जल्दी पकड़ने के लिए विशेष रूप से उपयोगी है। यह निष्पादन से पहले वर्कफ़्लो संरचना और process कनेक्शन को मान्य करता है।

### 4.3. Logic Testing के लिए Stub Running

कभी-कभी एरर को डिबग करना मुश्किल होता है क्योंकि कमांड बहुत लंबे समय तक चलते हैं, विशेष सॉफ़्टवेयर की आवश्यकता होती है, या जटिल कारणों से विफल होते हैं। Stub running आपको वास्तविक कमांड निष्पादित किए बिना वर्कफ़्लो logic का परीक्षण करने देता है।

#### पाइपलाइन चलाएं

जब आप एक Nextflow process विकसित कर रहे हों, तो आप `stub` निर्देश का उपयोग करके 'डमी' कमांड परिभाषित कर सकते हैं जो वास्तविक कमांड चलाए बिना सही रूप के आउटपुट उत्पन्न करते हैं। यह दृष्टिकोण विशेष रूप से तब मूल्यवान होता है जब आप वास्तविक सॉफ़्टवेयर की जटिलताओं से निपटने से पहले यह सत्यापित करना चाहते हैं कि आपकी वर्कफ़्लो logic सही है।

उदाहरण के लिए, पहले के हमारे `missing_software.nf` को याद करें? जहां हमारे पास गुम सॉफ़्टवेयर था जिसने `-profile docker` जोड़ने तक वर्कफ़्लो को चलने से रोक दिया? `missing_software_with_stub.nf` एक बहुत ही समान वर्कफ़्लो है। यदि हम इसे उसी तरह चलाते हैं, तो हम वही एरर उत्पन्न करेंगे:

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

हालांकि, यदि हम इसे `-stub-run` के साथ चलाते हैं, तो `docker` प्रोफ़ाइल के बिना भी यह वर्कफ़्लो एरर उत्पन्न नहीं करेगा:

```bash
nextflow run missing_software_with_stub.nf -stub-run
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software_with_stub.nf` [astonishing_shockley] DSL2 - revision: f1f4f05d7d

    executor >  local (3)
    [b5/2517a3] PROCESS_FILES (3) | 3 of 3 ✔
    ```

#### कोड की जांच करें

आइए `missing_software_with_stub.nf` की जांच करें:

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

`missing_software.nf` के सापेक्ष, इस process में एक `stub:` निर्देश है जो `script:` में निर्दिष्ट कमांड के बजाय उपयोग किए जाने वाले कमांड को निर्दिष्ट करता है, जब Nextflow stub mode में चलाया जाता है।

हम यहाँ जिस `touch` कमांड का उपयोग कर रहे हैं वह किसी भी सॉफ़्टवेयर या उचित इनपुट पर निर्भर नहीं करता है, और सभी स्थितियों में चलेगा, जो हमें process के आंतरिक भागों की चिंता किए बिना वर्कफ़्लो logic को डिबग करने की अनुमति देता है।

**Stub running डिबग करने में मदद करता है:**

- Channel संरचना और डेटा फ्लो
- Process कनेक्शन और निर्भरताएं
- पैरामीटर प्रसार
- सॉफ़्टवेयर निर्भरताओं के बिना वर्कफ़्लो logic

### 4.4. व्यवस्थित डिबगिंग दृष्टिकोण

अब जब आपने व्यक्तिगत डिबगिंग तकनीकें सीख ली हैं - ट्रेस फ़ाइलों और work डायरेक्टरी से लेकर preview mode, stub running, और संसाधन निगरानी तक - आइए उन्हें एक व्यवस्थित पद्धति में एक साथ बांधें। एक संरचित दृष्टिकोण होने से आप जटिल एरर से अभिभूत होने से बचते हैं और सुनिश्चित करते हैं कि आप महत्वपूर्ण सुरागों को न चूकें।

यह पद्धति हमारे द्वारा कवर किए गए सभी टूल्स को एक कुशल वर्कफ़्लो में जोड़ती है:

**चार-चरण डिबगिंग विधि:**

**चरण 1: सिंटैक्स एरर समाधान (5 मिनट)**

1. VSCode या अपने IDE में लाल रेखांकन की जांच करें
2. सिंटैक्स समस्याओं की पहचान के लिए `nextflow run workflow.nf -preview` चलाएं
3. सभी सिंटैक्स एरर ठीक करें (गुम ब्रेसिज़, trailing कॉमा, आदि)
4. आगे बढ़ने से पहले सुनिश्चित करें कि वर्कफ़्लो सफलतापूर्वक पार्स होता है

**चरण 2: त्वरित आकलन (5 मिनट)**

1. रनटाइम एरर मैसेज को ध्यान से पढ़ें
2. जांचें कि यह रनटाइम, logic, या संसाधन एरर है
3. बुनियादी वर्कफ़्लो logic का परीक्षण करने के लिए preview mode का उपयोग करें

**चरण 3: विस्तृत जांच (15-30 मिनट)**

1. विफल कार्य की work डायरेक्टरी खोजें
2. लॉग फ़ाइलों की जांच करें
3. Channels का निरीक्षण करने के लिए `.view()` operators जोड़ें
4. निष्पादन के बिना वर्कफ़्लो logic का परीक्षण करने के लिए `-stub-run` का उपयोग करें

**चरण 4: ठीक करें और मान्य करें (15 मिनट)**

1. न्यूनतम लक्षित सुधार करें
2. resume के साथ परीक्षण करें: `nextflow run workflow.nf -resume`
3. पूर्ण वर्कफ़्लो निष्पादन सत्यापित करें

!!! tip "कुशल डिबगिंग के लिए Resume का उपयोग"

    एक बार जब आपने एक समस्या की पहचान कर ली है, तो आपको अपने वर्कफ़्लो के सफल भागों को फिर से चलाने में समय बर्बाद किए बिना अपने सुधारों का परीक्षण करने के लिए एक कुशल तरीके की आवश्यकता है। Nextflow की `-resume` कार्यक्षमता डिबगिंग के लिए अमूल्य है।

    यदि आपने [Hello Nextflow](../hello_nextflow/) के माध्यम से काम किया है तो आपने `-resume` का सामना किया होगा, और यह महत्वपूर्ण है कि आप डिबगिंग करते समय इसका अच्छा उपयोग करें ताकि आपकी समस्या process से पहले की processes चलने के दौरान प्रतीक्षा से बच सकें।

    **Resume डिबगिंग रणनीति:**

    1. विफलता तक वर्कफ़्लो चलाएं
    2. विफल कार्य के लिए work डायरेक्टरी की जांच करें
    3. विशिष्ट समस्या ठीक करें
    4. केवल सुधार का परीक्षण करने के लिए resume करें
    5. वर्कफ़्लो पूरा होने तक दोहराएं

#### डिबगिंग कॉन्फ़िगरेशन प्रोफ़ाइल

इस व्यवस्थित दृष्टिकोण को और भी अधिक कुशल बनाने के लिए, आप एक समर्पित डिबगिंग कॉन्फ़िगरेशन बना सकते हैं जो स्वचालित रूप से आपको आवश्यक सभी टूल्स सक्षम करता है:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Conservative resources for debugging
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

फिर आप इस प्रोफ़ाइल को सक्षम करके पाइपलाइन चला सकते हैं:

```bash
nextflow run workflow.nf -profile debug
```

यह प्रोफ़ाइल रियल-टाइम आउटपुट सक्षम करती है, work डायरेक्टरी को संरक्षित करती है, और आसान डिबगिंग के लिए समानांतरता को सीमित करती है।

### 4.5. व्यावहारिक डिबगिंग अभ्यास

अब व्यवस्थित डिबगिंग दृष्टिकोण को अभ्यास में लाने का समय है। वर्कफ़्लो `buggy_workflow.nf` में कई सामान्य एरर हैं जो वास्तविक दुनिया के विकास में आपको मिलने वाली समस्याओं के प्रकारों का प्रतिनिधित्व करती हैं।

!!! exercise

    `buggy_workflow.nf` में सभी एरर की पहचान करने और उन्हें ठीक करने के लिए व्यवस्थित डिबगिंग दृष्टिकोण का उपयोग करें। यह वर्कफ़्लो CSV फ़ाइल से नमूना डेटा को प्रोसेस करने का प्रयास करता है लेकिन इसमें सामान्य डिबगिंग परिदृश्यों का प्रतिनिधित्व करने वाली कई जानबूझकर बग हैं।

    पहली एरर देखने के लिए वर्कफ़्लो चलाकर शुरू करें:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "कमांड आउटपुट"

        ```console
        N E X T F L O W   ~  version 25.10.2

        Launching `buggy_workflow.nf` [wise_ramanujan] DSL2 - revision: d51a8e83fd

        ERROR ~ Range [11, 12) out of bounds for length 11

         -- Check '.nextflow.log' file for details
        ```

        यह गूढ़ एरर `params{}` ब्लॉक में लाइन 11-12 के आसपास एक पार्सिंग समस्या को इंगित करती है। v2 पार्सर संरचनात्मक समस्याओं को जल्दी पकड़ता है।

    आपने जो चार-चरण डिबगिंग विधि सीखी है उसे लागू करें:

    **चरण 1: सिंटैक्स एरर समाधान**
    - VSCode या अपने IDE में लाल रेखांकन की जांच करें
    - सिंटैक्स समस्याओं की पहचान के लिए `nextflow run workflow.nf -preview` चलाएं
    - सभी सिंटैक्स एरर ठीक करें (गुम ब्रेसिज़, trailing कॉमा, आदि)
    - आगे बढ़ने से पहले सुनिश्चित करें कि वर्कफ़्लो सफलतापूर्वक पार्स होता है

    **चरण 2: त्वरित आकलन**
    - रनटाइम एरर मैसेज को ध्यान से पढ़ें
    - पहचानें कि एरर रनटाइम, logic, या संसाधन-संबंधित हैं
    - बुनियादी वर्कफ़्लो logic का परीक्षण करने के लिए `-preview` mode का उपयोग करें

    **चरण 3: विस्तृत जांच**
    - विफल कार्यों के लिए work डायरेक्टरी की जांच करें
    - Channels का निरीक्षण करने के लिए `.view()` operators जोड़ें
    - Work डायरेक्टरी में लॉग फ़ाइलों की जांच करें
    - निष्पादन के बिना वर्कफ़्लो logic का परीक्षण करने के लिए `-stub-run` का उपयोग करें

    **चरण 4: ठीक करें और मान्य करें**
    - लक्षित सुधार करें
    - सुधारों को कुशलतापूर्वक परीक्षण करने के लिए `-resume` का उपयोग करें
    - पूर्ण वर्कफ़्लो निष्पादन सत्यापित करें

    **आपके पास उपलब्ध डिबगिंग टूल्स:**
    ```bash
    # सिंटैक्स जांच के लिए Preview mode
    nextflow run buggy_workflow.nf -preview

    # विस्तृत आउटपुट के लिए Debug profile
    nextflow run buggy_workflow.nf -profile debug

    # Logic testing के लिए Stub running
    nextflow run buggy_workflow.nf -stub-run

    # सुधार के बाद Resume
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution
        `buggy_workflow.nf` में 9 या 10 अलग-अलग एरर हैं (आप कैसे गिनते हैं इसके आधार पर) जो सभी प्रमुख डिबगिंग श्रेणियों को कवर करती हैं। यहाँ प्रत्येक एरर और इसे कैसे ठीक करें का एक व्यवस्थित विश्लेषण है

        सिंटैक्स एरर से शुरू करते हैं:

        **एरर 1: सिंटैक्स एरर - Trailing कॉमा**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // ERROR: Trailing comma
        ```
        **फिक्स:** Trailing कॉमा हटाएं
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **एरर 2: सिंटैक्स एरर - गुम Closing Brace**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // ERROR: Missing closing brace for processFiles process
        ```
        **फिक्स:** गुम closing brace जोड़ें
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Add missing closing brace
        ```

        **एरर 3: वेरिएबल नाम एरर**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // ERROR: should be sample_id
        cat ${input_file} > ${sample}_result.txt  // ERROR: should be sample_id
        ```
        **फिक्स:** सही इनपुट वेरिएबल नाम का उपयोग करें
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **एरर 4: अपरिभाषित वेरिएबल एरर**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // ERROR: sample_ids undefined
        ```
        **फिक्स:** सही channel का उपयोग करें और sample IDs निकालें
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        इस बिंदु पर वर्कफ़्लो चलेगा, लेकिन हमें अभी भी एरर मिलेंगी (जैसे `processFiles` में `Path value cannot be null`), जो खराब channel संरचना के कारण है।

        **एरर 5: Channel संरचना एरर - गलत Map आउटपुट**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // ERROR: processFiles expects tuple
        ```
        **फिक्स:** processFiles की अपेक्षित tuple संरचना लौटाएं
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        लेकिन यह ऊपर `heavyProcess()` चलाने के हमारे फिक्स को तोड़ देगा, इसलिए हमें उस process को केवल sample IDs पास करने के लिए map का उपयोग करना होगा:

        **एरर 6: heavyProcess के लिए खराब channel संरचना**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // ERROR: input_ch now has 2 elements per emission- heavyProcess only needs 1 (the first)
        ```
        **फिक्स:** सही channel का उपयोग करें और sample IDs निकालें
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        अब हम थोड़ा आगे बढ़ते हैं लेकिन `No such variable: i` के बारे में एरर प्राप्त करते हैं, क्योंकि हमने एक Bash वेरिएबल को escape नहीं किया।

        **एरर 7: Bash वेरिएबल Escaping एरर**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // ERROR: $i not escaped
        ```
        **फिक्स:** bash वेरिएबल को escape करें
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        अब हमें `Process exceeded running time limit (1ms)` मिलता है, इसलिए हम संबंधित process के लिए रन टाइम लिमिट ठीक करते हैं:

        **एरर 8: संसाधन कॉन्फ़िगरेशन एरर**
        ```groovy linenums="36"
        time '1 ms'  // ERROR: Unrealistic time limit
        ```
        **फिक्स:** एक वास्तविक समय सीमा तक बढ़ाएं
        ```groovy linenums="36"
        time '100 s'
        ```

        अगला हमारे पास हल करने के लिए एक `Missing output file(s)` एरर है:

        **एरर 9: आउटपुट फ़ाइल नाम बेमेल**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // ERROR: Wrong filename, should match output declaration
        ```
        **फिक्स:** आउटपुट घोषणा से मिलान करें
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        पहली दो processes चलीं, लेकिन तीसरी नहीं।

        **एरर 10: आउटपुट फ़ाइल नाम बेमेल**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Error: attempting to take input from the pwd rather than a process
        handleFiles(file_ch)
        ```
        **फिक्स:** पिछली process से आउटपुट लें
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        इसके साथ, पूरा वर्कफ़्लो चलना चाहिए।

        **पूर्ण सही वर्कफ़्लो:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Buggy workflow for debugging exercises
        * This workflow contains several intentional bugs for learning purposes
        */

        params{
            // Parameters with missing validation
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Process with input/output mismatch
        */
        process processFiles {
            publishDir "${params.output}/processed", mode: 'copy'

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
        * Process with resource issues
        */
        process heavyProcess {
            publishDir "${params.output}/heavy", mode: 'copy'

            time '100 s'

            input:
                val sample_id

            output:
                path "${sample_id}_heavy.txt"

            script:
            """
            # Simulate heavy computation
            for i in {1..1000000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * Process with file handling issues
        */
        process handleFiles {
            publishDir "${params.output}/files", mode: 'copy'

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
        * Main workflow with channel issues
        */
        workflow {

            // Channel with incorrect usage
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map{it[0]})

            handleFiles(heavyProcess.out)
        }
        ```

**कवर की गई एरर श्रेणियां:**

- **सिंटैक्स एरर**: गुम ब्रेसिज़, trailing कॉमा, अपरिभाषित वेरिएबल
- **Channel संरचना एरर**: गलत डेटा आकार, अपरिभाषित channels
- **Process एरर**: आउटपुट फ़ाइल बेमेल, वेरिएबल escaping
- **संसाधन एरर**: अवास्तविक समय सीमाएं

**मुख्य डिबगिंग पाठ:**

1. **एरर मैसेज ध्यान से पढ़ें** - वे अक्सर सीधे समस्या की ओर इंगित करते हैं
2. **व्यवस्थित दृष्टिकोण का उपयोग करें** - एक बार में एक एरर ठीक करें और `-resume` के साथ परीक्षण करें
3. **डेटा फ्लो को समझें** - channel संरचना एरर अक्सर सबसे सूक्ष्म होती हैं
4. **Work डायरेक्टरी की जांच करें** - जब processes विफल होते हैं, तो लॉग आपको बिल्कुल बताते हैं कि क्या गलत हुआ

---

## सारांश

इस साइड क्वेस्ट में, आपने Nextflow वर्कफ़्लो को डिबग करने के लिए व्यवस्थित तकनीकों का एक सेट सीखा है।
इन तकनीकों को अपने काम में लागू करने से आप अपने कंप्यूटर से लड़ने में कम समय बिताने, समस्याओं को तेज़ी से हल करने और भविष्य की समस्याओं से खुद को बचाने में सक्षम होंगे।

### मुख्य पैटर्न

**1. सिंटैक्स एरर की पहचान और ठीक कैसे करें**:

- Nextflow एरर मैसेज की व्याख्या और समस्याओं का पता लगाना
- सामान्य सिंटैक्स एरर: गुम ब्रेसिज़, गलत कीवर्ड, अपरिभाषित वेरिएबल
- Nextflow (Groovy) और Bash वेरिएबल के बीच अंतर करना
- प्रारंभिक एरर डिटेक्शन के लिए VS Code एक्सटेंशन फीचर्स का उपयोग करना

```groovy
// गुम brace - IDE में लाल रेखांकन देखें
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- गुम है!

// गलत कीवर्ड
inputs:  // 'input:' होना चाहिए

// अपरिभाषित वेरिएबल - Bash वेरिएबल के लिए बैकस्लैश से escape करें
echo "${undefined_var}"      // Nextflow वेरिएबल (एरर यदि परिभाषित नहीं)
echo "\${bash_var}"          // Bash वेरिएबल (escaped)
```

**2. Channel संरचना समस्याओं को डिबग कैसे करें**:

- Channel cardinality और exhaustion समस्याओं को समझना
- Channel सामग्री संरचना बेमेल को डिबग करना
- Channel निरीक्षण के लिए `.view()` operators का उपयोग करना
- आउटपुट में वर्ग कोष्ठक जैसे एरर पैटर्न को पहचानना

```groovy
// Channel सामग्री का निरीक्षण करें
my_channel.view { "Content: $it" }

// Queue को value channel में बदलें (exhaustion को रोकता है)
reference_ch = channel.value('ref.fa')
// या
reference_ch = channel.of('ref.fa').first()
```

**3. Process निष्पादन समस्याओं का समाधान कैसे करें**:

- गुम आउटपुट फ़ाइल एरर का निदान करना
- एग्जिट कोड को समझना (गुम सॉफ़्टवेयर के लिए 127, मेमोरी समस्याओं के लिए 137)
- Work डायरेक्टरी और कमांड फ़ाइलों की जांच करना
- संसाधनों को उचित रूप से कॉन्फ़िगर करना

```bash
# जांचें कि वास्तव में क्या निष्पादित किया गया
cat work/ab/cdef12/.command.sh

# एरर आउटपुट जांचें
cat work/ab/cdef12/.command.err

# एग्जिट कोड 127 = कमांड नहीं मिला
# एग्जिट कोड 137 = समाप्त किया गया (मेमोरी/समय सीमा)
```

**4. Nextflow के बिल्ट-इन डिबगिंग टूल्स का उपयोग कैसे करें**:

- Preview mode और रियल-टाइम डिबगिंग का लाभ उठाना
- Logic testing के लिए stub running को लागू करना
- कुशल डिबगिंग चक्रों के लिए resume लागू करना
- चार-चरण व्यवस्थित डिबगिंग पद्धति का पालन करना

!!! tip "त्वरित डिबगिंग संदर्भ"

    **सिंटैक्स एरर?** → VSCode चेतावनियां जांचें, `nextflow run workflow.nf -preview` चलाएं

    **Channel समस्याएं?** → सामग्री का निरीक्षण करने के लिए `.view()` का उपयोग करें: `my_channel.view()`

    **Process विफलताएं?** → Work डायरेक्टरी फ़ाइलों की जांच करें:

    - `.command.sh` - निष्पादित स्क्रिप्ट
    - `.command.err` - एरर मैसेज
    - `.exitcode` - एग्जिट स्टेटस (127 = कमांड नहीं मिला, 137 = समाप्त किया गया)

    **रहस्यमय व्यवहार?** → वर्कफ़्लो logic का परीक्षण करने के लिए `-stub-run` के साथ चलाएं

    **सुधार किए?** → परीक्षण में समय बचाने के लिए `-resume` का उपयोग करें: `nextflow run workflow.nf -resume`

---

### अतिरिक्त संसाधन

- [Nextflow समस्या निवारण गाइड](https://www.nextflow.io/docs/latest/troubleshooting.html): आधिकारिक समस्या निवारण दस्तावेज़
- [Nextflow channels को समझना](https://www.nextflow.io/docs/latest/channel.html): Channel प्रकारों और व्यवहार में गहराई से गोता लगाएं
- [Process directives संदर्भ](https://www.nextflow.io/docs/latest/process.html#directives): सभी उपलब्ध process कॉन्फ़िगरेशन विकल्प
- [nf-test](https://www.nf-test.com/): Nextflow pipelines के लिए परीक्षण फ्रेमवर्क
- [Nextflow Slack समुदाय](https://www.nextflow.io/slack-invite.html): समुदाय से सहायता प्राप्त करें

उत्पादन वर्कफ़्लो के लिए, विचार करें:

- स्केल पर निगरानी और डिबगिंग के लिए [Seqera Platform](https://seqera.io/platform/) सेट अप करना
- पुनरुत्पादनीय सॉफ़्टवेयर वातावरण के लिए [Wave containers](https://seqera.io/wave/) का उपयोग करना

**याद रखें:** प्रभावी डिबगिंग एक कौशल है जो अभ्यास के साथ सुधरती है। यहाँ आपने जो व्यवस्थित पद्धति और व्यापक टूलकिट प्राप्त किया है, वह आपके पूरे Nextflow विकास यात्रा में आपकी अच्छी सेवा करेगा।

---

## आगे क्या है?

[Side Quests के मेनू](./index.md) पर लौटें या अगले विषय पर जाने के लिए पृष्ठ के निचले दाएं कोने में बटन पर क्लिक करें।
