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

    -- Check '.nextflow.log' file for
