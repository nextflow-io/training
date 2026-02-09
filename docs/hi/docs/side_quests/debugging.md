# वर्कफ़्लो डिबगिंग

डिबगिंग एक महत्वपूर्ण कौशल है जो तुम्हें घंटों की निराशा से बचा सकता है और तुम्हें एक अधिक प्रभावी Nextflow डेवलपर बनने में मदद कर सकता है। अपने करियर के दौरान, खासकर जब तुम शुरुआत कर रहे हो, तुम्हें अपने वर्कफ़्लो बनाते और बनाए रखते समय बग का सामना करना पड़ेगा। व्यवस्थित डिबगिंग दृष्टिकोण सीखने से तुम्हें समस्याओं की पहचान करने और उन्हें जल्दी हल करने में मदद मिलेगी।

### सीखने के लक्ष्य

इस साइड क्वेस्ट में, हम Nextflow वर्कफ़्लो के लिए **व्यवस्थित डिबगिंग तकनीकों** का पता लगाएंगे:

- **सिंटैक्स एरर डिबगिंग**: IDE फीचर्स और Nextflow एरर मैसेज का प्रभावी ढंग से उपयोग करना
- **चैनल डिबगिंग**: डेटा फ़्लो समस्याओं और चैनल संरचना समस्याओं का निदान करना
- **प्रोसेस डिबगिंग**: एक्ज़ीक्यूशन विफलताओं और रिसोर्स समस्याओं की जांच करना
- **बिल्ट-इन डिबगिंग टूल्स**: Nextflow के प्रीव्यू मोड, स्टब रनिंग, और वर्क डायरेक्टरी का लाभ उठाना
- **व्यवस्थित दृष्टिकोण**: कुशल डिबगिंग के लिए चार-चरण पद्धति

अंत तक, तुम्हारे पास एक मजबूत डिबगिंग पद्धति होगी जो निराशाजनक एरर मैसेज को समाधान के लिए स्पष्ट रोडमैप में बदल देती है।

### पूर्वापेक्षाएँ

इस साइड क्वेस्ट को शुरू करने से पहले, तुम्हें:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष शुरुआती कोर्स पूरा कर लेना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (प्रोसेस, चैनल, ऑपरेटर) का उपयोग करने में सहज होना चाहिए

**वैकल्पिक:** हम पहले [IDE Features for Nextflow Development](./ide_features.md) साइड क्वेस्ट पूरा करने की सलाह देते हैं।
यह IDE फीचर्स का व्यापक कवरेज प्रदान करता है जो डिबगिंग का समर्थन करते हैं (सिंटैक्स हाइलाइटिंग, एरर डिटेक्शन, आदि), जिनका हम यहाँ भारी उपयोग करेंगे।

---

## 0. शुरू करना

#### ट्रेनिंग कोडस्पेस खोलें

अगर तुमने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित अनुसार ट्रेनिंग एनवायरनमेंट खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाएं

चलो उस डायरेक्टरी में चलते हैं जहाँ इस ट्यूटोरियल के लिए फ़ाइलें स्थित हैं।

```bash
cd side-quests/debugging
```

तुम VSCode को इस डायरेक्टरी पर फोकस करने के लिए सेट कर सकते हो:

```bash
code .
```

#### सामग्री की समीक्षा करें

तुम्हें विभिन्न प्रकार के बग वाले उदाहरण वर्कफ़्लो का एक सेट मिलेगा जिनका हम अभ्यास के लिए उपयोग करेंगे:

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

ये फ़ाइलें सामान्य डिबगिंग परिदृश्यों का प्रतिनिधित्व करती हैं जिनका तुम वास्तविक दुनिया के विकास में सामना करोगे।

#### असाइनमेंट की समीक्षा करें

तुम्हारी चुनौती है प्रत्येक वर्कफ़्लो को चलाना, एरर की पहचान करना, और उन्हें ठीक करना।

प्रत्येक बगी वर्कफ़्लो के लिए:

1. **वर्कफ़्लो चलाओ** और एरर देखो
2. **एरर मैसेज का विश्लेषण करो**: Nextflow तुम्हें क्या बता रहा है?
3. **कोड में समस्या का पता लगाओ** दिए गए सुरागों का उपयोग करके
4. **बग को ठीक करो** और सत्यापित करो कि तुम्हारा समाधान काम करता है
5. अगले सेक्शन पर जाने से पहले **फ़ाइल को रीसेट करो** (`git checkout <filename>` का उपयोग करो)

अभ्यास सरल सिंटैक्स एरर से अधिक सूक्ष्म रनटाइम समस्याओं तक प्रगति करते हैं।
समाधान इनलाइन चर्चा किए गए हैं, लेकिन आगे पढ़ने से पहले प्रत्येक को स्वयं हल करने का प्रयास करो।

#### तैयारी चेकलिस्ट

लगता है कि तुम गोता लगाने के लिए तैयार हो?

- [ ] मैं इस कोर्स के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता हूँ
- [ ] मेरा कोडस्पेस चालू है और चल रहा है
- [ ] मैंने अपनी वर्किंग डायरेक्टरी उचित रूप से सेट कर ली है
- [ ] मैं असाइनमेंट को समझता हूँ

अगर तुम सभी बॉक्स चेक कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. सिंटैक्स एरर

सिंटैक्स एरर सबसे आम प्रकार की एरर हैं जिनका तुम Nextflow कोड लिखते समय सामना करोगे। ये तब होती हैं जब कोड Nextflow DSL के अपेक्षित सिंटैक्स नियमों के अनुरूप नहीं होता है। ये एरर तुम्हारे वर्कफ़्लो को बिल्कुल भी चलने से रोकती हैं, इसलिए उन्हें जल्दी पहचानना और ठीक करना सीखना महत्वपूर्ण है।

### 1.1. लापता ब्रेसेस

सबसे आम सिंटैक्स एरर में से एक, और कभी-कभी डिबग करने के लिए अधिक जटिल में से एक है **लापता या बेमेल ब्रैकेट**।

चलो एक व्यावहारिक उदाहरण से शुरू करते हैं।

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
- **EOF संकेतक**: `<EOF>` (End Of File) मैसेज इंगित करता है कि पार्सर फ़ाइल के अंत तक पहुंच गया जबकि अभी भी अधिक सामग्री की उम्मीद थी - अनक्लोज़्ड ब्रेसेस का एक क्लासिक संकेत

#### कोड जांचें

अब, चलो `bad_syntax.nf` की जांच करते हैं यह समझने के लिए कि एरर का कारण क्या है:

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
// Missing closing brace for the process

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
```

इस उदाहरण के उद्देश्य के लिए हमने तुम्हारे लिए एक कमेंट छोड़ा है जो दिखाता है कि एरर कहाँ है। Nextflow VSCode एक्सटेंशन को भी तुम्हें कुछ संकेत देने चाहिए कि क्या गलत हो सकता है, बेमेल ब्रेस को लाल में डालते हुए और फ़ाइल के समय से पहले अंत को हाइलाइट करते हुए:

![Bad syntax](img/bad_syntax.png)

**ब्रैकेट एरर के लिए डिबगिंग रणनीति:**

1. VS Code की ब्रैकेट मैचिंग का उपयोग करें (कर्सर को ब्रैकेट के बगल में रखें)
2. ब्रैकेट-संबंधित मैसेज के लिए Problems पैनल जांचें
3. सुनिश्चित करें कि प्रत्येक ओपनिंग `{` का एक संबंधित क्लोजिंग `}` है

#### कोड ठीक करें

कमेंट को लापता क्लोजिंग ब्रेस से बदलें:

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
    }  // Add the missing closing brace

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
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
    // Missing closing brace for the process

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

#### पाइपलाइन चलाएं

अब वर्कफ़्लो को फिर से चलाओ यह पुष्टि करने के लिए कि यह काम करता है:

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

### 1.2. गलत प्रोसेस कीवर्ड या निर्देशों का उपयोग करना

एक और आम सिंटैक्स एरर एक **अमान्य प्रोसेस परिभाषा** है। यह तब हो सकता है जब तुम आवश्यक ब्लॉक को परिभाषित करना भूल जाते हो या प्रोसेस परिभाषा में गलत निर्देशों का उपयोग करते हो।

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

#### कोड जांचें

एरर "Invalid process definition" इंगित करता है और समस्या के आसपास के संदर्भ को दिखाता है। लाइन 3-7 को देखते हुए, हम लाइन 4 पर `inputs:` देख सकते हैं, जो समस्या है। चलो `invalid_process.nf` की जांच करते हैं:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // ERROR: Should be 'input' not 'inputs'
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
```

एरर संदर्भ में लाइन 4 को देखते हुए, हम समस्या को पहचान सकते हैं: हम सही `input` निर्देश के बजाय `inputs` का उपयोग कर रहे हैं। Nextflow VSCode एक्सटेंशन भी इसे फ्लैग करेगा:

![Invalid process message](img/invalid_process_message.png)

#### कोड ठीक करें

[डॉक्यूमेंटेशन](https://www.nextflow.io/docs/latest/process.html#) का संदर्भ देकर गलत कीवर्ड को सही से बदलें:

=== "बाद में"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Fixed: Changed 'inputs' to 'input'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

=== "पहले"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // ERROR: Should be 'input' not 'inputs'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

#### पाइपलाइन चलाएं

अब वर्कफ़्लो को फिर से चलाओ यह पुष्टि करने के लिए कि यह काम करता है:

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

तुम्हारे स्क्रिप्ट ब्लॉक में उपयोग किए जाने वाले वेरिएबल नाम मान्य होने चाहिए, या तो इनपुट से प्राप्त या स्क्रिप्ट से पहले डाले गए groovy कोड से। लेकिन जब तुम पाइपलाइन विकास की शुरुआत में जटिलता से जूझ रहे होते हो, तो वेरिएबल नामकरण में गलतियाँ करना आसान है, और Nextflow तुम्हें जल्दी बता देगा।

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

एरर कंपाइल टाइम पर पकड़ी जाती है और सीधे लाइन 17 पर अपरिभाषित वेरिएबल की ओर इशारा करती है, एक कैरेट के साथ जो बिल्कुल इंगित करता है कि समस्या कहाँ है।

#### कोड जांचें

चलो `no_such_var.nf` की जांच करते हैं:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

एरर मैसेज इंगित करता है कि वेरिएबल स्क्रिप्ट टेम्पलेट में पहचाना नहीं गया है, और वहाँ तुम जाओ- तुम्हें स्क्रिप्ट ब्लॉक में उपयोग किया गया `${undefined_var}` देखने में सक्षम होना चाहिए, लेकिन कहीं और परिभाषित नहीं।

#### कोड ठीक करें

अगर तुम्हें 'No such variable' एरर मिलती है, तो तुम इसे या तो वेरिएबल को परिभाषित करके (इनपुट वेरिएबल नामों को सही करके या स्क्रिप्ट से पहले groovy कोड को संपादित करके), या स्क्रिप्ट ब्लॉक से इसे हटाकर ठीक कर सकते हो अगर इसकी आवश्यकता नहीं है:

=== "बाद में"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Removed the line with undefined_var
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
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### पाइपलाइन चलाएं

अब वर्कफ़्लो को फिर से चलाओ यह पुष्टि करने के लिए कि यह काम करता है:

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

Nextflow में शुरुआत करते समय, Nextflow (Groovy) और Bash वेरिएबल के बीच अंतर को समझना मुश्किल हो सकता है। यह खराब वेरिएबल एरर का एक और रूप उत्पन्न कर सकता है जो स्क्रिप्ट ब्लॉक की Bash सामग्री में वेरिएबल का उपयोग करने की कोशिश करते समय दिखाई देता है।

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

#### कोड जांचें

एरर लाइन 13 की ओर इशारा करती है जहाँ `${prefix}` का उपयोग किया गया है। चलो `bad_bash_var.nf` की जांच करते हैं यह देखने के लिए कि समस्या का कारण क्या है:

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
    echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} Groovy सिंटैक्स है, Bash नहीं
    """
}
```

इस उदाहरण में, हम Bash में `prefix` वेरिएबल को परिभाषित कर रहे हैं, लेकिन एक Nextflow प्रोसेस में `$` सिंटैक्स जिसका हमने इसे संदर्भित करने के लिए उपयोग किया (`${prefix}`) को Groovy वेरिएबल के रूप में व्याख्यायित किया जाता है, Bash नहीं। वेरिएबल Groovy संदर्भ में मौजूद नहीं है, इसलिए हमें 'no such variable' एरर मिलती है।

#### कोड ठीक करें

अगर तुम Bash वेरिएबल का उपयोग करना चाहते हो, तो तुम्हें डॉलर साइन को इस तरह एस्केप करना होगा:

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
        echo "Processing ${sample_name}" > \${prefix}.txt  # ठीक किया: डॉलर साइन को एस्केप किया
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
        echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} Groovy सिंटैक्स है, Bash नहीं
        """
    }
    ```

यह Nextflow को बताता है कि इसे Bash वेरिएबल के रूप में व्याख्यायित करे।

#### पाइपलाइन चलाएं

अब वर्कफ़्लो को फिर से चलाओ यह पुष्टि करने के लिए कि यह काम करता है:

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

    स्ट्रिंग कॉन्कैटेनेशन या प्रीफिक्स/सफिक्स ऑपरेशन जैसे सरल वेरिएबल मैनिपुलेशन के लिए, आमतौर पर स्क्रिप्ट ब्लॉक में Bash वेरिएबल के बजाय स्क्रिप्ट सेक्शन में Groovy वेरिएबल का उपयोग करना अधिक पठनीय होता है:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    यह दृष्टिकोण डॉलर साइन को एस्केप करने की आवश्यकता से बचता है और कोड को पढ़ने और बनाए रखने में आसान बनाता है।

### 1.5. Workflow ब्लॉक के बाहर स्टेटमेंट

Nextflow VSCode एक्सटेंशन कोड संरचना के साथ समस्याओं को हाइलाइट करता है जो एरर का कारण बनेंगी। एक सामान्य उदाहरण `workflow {}` ब्लॉक के बाहर चैनल को परिभाषित करना है - यह अब एक सिंटैक्स एरर के रूप में लागू किया गया है।

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

एरर मैसेज स्पष्ट रूप से समस्या को इंगित करता है: स्टेटमेंट (जैसे चैनल परिभाषाएँ) को workflow या process ब्लॉक के बाहर स्क्रिप्ट घोषणाओं के साथ मिश्रित नहीं किया जा सकता है।

#### कोड जांचें

चलो `badpractice_syntax.nf` की जांच करते हैं यह देखने के लिए कि एरर का कारण क्या है:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Channel defined outside workflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
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

#### कोड ठीक करें

चैनल परिभाषा को workflow ब्लॉक के अंदर ले जाएं:

=== "बाद में"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Moved inside workflow block
        PROCESS_FILES(input_ch)
    }
    ```

=== "पहले"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Channel defined outside workflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
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

फिक्स काम करती है यह पुष्टि करने के लिए वर्कफ़्लो को फिर से चलाओ:

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

अपने इनपुट चैनल को workflow ब्लॉक के भीतर परिभाषित रखो, और सामान्य रूप से एक्सटेंशन द्वारा दी गई किसी भी अन्य सिफारिशों का पालन करो।

### सारांश

तुम Nextflow एरर मैसेज और IDE विज़ुअल संकेतकों का उपयोग करके व्यवस्थित रूप से सिंटैक्स एरर की पहचान और उन्हें ठीक कर सकते हो। सामान्य सिंटैक्स एरर में लापता ब्रेसेस, गलत प्रोसेस कीवर्ड, अपरिभाषित वेरिएबल, और Bash बनाम Nextflow वेरिएबल का अनुचित उपयोग शामिल है। VSCode एक्सटेंशन रनटाइम से पहले इनमें से कई को पकड़ने में मदद करता है। अपने टूलकिट में इन सिंटैक्स डिबगिंग कौशल के साथ, तुम सबसे आम Nextflow सिंटैक्स एरर को जल्दी से हल करने और अधिक जटिल रनटाइम समस्याओं से निपटने के लिए आगे बढ़ने में सक्षम होगे।

### आगे क्या है?

अधिक जटिल चैनल संरचना एरर को डिबग करना सीखो जो तब भी होती हैं जब सिंटैक्स सही हो।

---

## 2. चैनल संरचना एरर

चैनल संरचना एरर सिंटैक्स एरर से अधिक सूक्ष्म हैं क्योंकि कोड सिंटैक्टिक रूप से सही है, लेकिन डेटा आकार प्रोसेस की अपेक्षाओं से मेल नहीं खाते हैं। Nextflow पाइपलाइन को चलाने की कोशिश करेगा, लेकिन यह पा सकता है कि इनपुट की संख्या उसकी अपेक्षा से मेल नहीं खाती है और विफल हो जाता है। ये एरर आमतौर पर केवल रनटाइम पर दिखाई देती हैं और तुम्हारे वर्कफ़्लो के माध्यम से बहने वाले डेटा की समझ की आवश्यकता होती है।

!!! tip "`.view()` के साथ चैनल डिबग करना"

    इस पूरे सेक्शन में, याद रखो कि तुम अपने वर्कफ़्लो में किसी भी बिंदु पर चैनल सामग्री का निरीक्षण करने के लिए `.view()` ऑपरेटर का उपयोग कर सकते हो। यह चैनल संरचना समस्याओं को समझने के लिए सबसे शक्तिशाली डिबगिंग टूल में से एक है। हम सेक्शन 2.4 में इस तकनीक का विस्तार से पता लगाएंगे, लेकिन उदाहरणों के माध्यम से काम करते समय इसका उपयोग करने के लिए स्वतंत्र महसूस करो।

    ```groovy
    my_channel.view()  // Shows what's flowing through the channel
    ```

### 2.1. गलत संख्या में इनपुट चैनल

यह एरर तब होती है जब तुम एक प्रोसेस की अपेक्षा से अलग संख्या में चैनल पास करते हो।

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

#### कोड जांचें

एरर मैसेज स्पष्ट रूप से बताता है कि कॉल को 1 आर्गुमेंट की उम्मीद थी लेकिन 2 प्राप्त हुए, और लाइन 23 की ओर इशारा करता है। चलो `bad_number_inputs.nf` की जांच करते हैं:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Process expects only 1 input

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create two separate channels
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: Passing 2 channels but process expects only 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

तुम्हें बेमेल `PROCESS_FILES` कॉल देखनी चाहिए, जो कई इनपुट चैनल की आपूर्ति करती है जब प्रोसेस केवल एक को परिभाषित करता है। VSCode एक्सटेंशन भी प्रोसेस कॉल को लाल में रेखांकित करेगा, और जब तुम माउस ओवर करोगे तो एक डायग्नोस्टिक मैसेज प्रदान करेगा:

![Incorrect number of args message](img/incorrect_num_args.png)

#### कोड ठीक करें

इस विशिष्ट उदाहरण के लिए, प्रोसेस एक एकल चैनल की उम्मीद करता है और दूसरे चैनल की आवश्यकता नहीं है, इसलिए हम केवल `samples_ch` चैनल पास करके इसे ठीक कर सकते हैं:

=== "बाद में"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Fixed: Pass only the channel the process expects
        PROCESS_FILES(samples_ch)
    }
    ```

=== "पहले"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERROR: Passing 2 channels but process expects only 1
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

इस उदाहरण की तुलना में अधिक सामान्यतः, तुम एक प्रोसेस में अतिरिक्त इनपुट जोड़ सकते हो और तदनुसार वर्कफ़्लो कॉल को अपडेट करना भूल सकते हो, जो इस प्रकार की एरर का कारण बन सकता है। सौभाग्य से, यह समझने और ठीक करने में आसान एरर में से एक है, क्योंकि एरर मैसेज बेमेल के बारे में काफी स्पष्ट है।

### 2.2. चैनल एग्ज़ॉशन (प्रोसेस अपेक्षा से कम बार चलता है)

कुछ चैनल संरचना एरर बहुत अधिक सूक्ष्म हैं और बिल्कुल भी कोई एरर उत्पन्न नहीं करती हैं। शायद इनमें से सबसे आम एक चुनौती को दर्शाता है जिसका नए Nextflow उपयोगकर्ता सामना करते हैं यह समझने में कि queue चैनल समाप्त हो सकते हैं और आइटम खत्म हो सकते हैं, जिसका अर्थ है कि वर्कफ़्लो समय से पहले समाप्त हो जाता है।

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

यह वर्कफ़्लो बिना एरर के पूरा होता है, लेकिन यह केवल एक एकल नमूने को प्रोसेस करता है!

#### कोड जांचें

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
    // Define variables in Groovy code before the script
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

प्रोसेस तीन बार के बजाय केवल एक बार चलता है क्योंकि `reference_ch` चैनल एक queue चैनल है जो पहले प्रोसेस एक्ज़ीक्यूशन के बाद समाप्त हो जाता है। जब एक चैनल समाप्त हो जाता है, तो पूरी प्रोसेस रुक जाती है, भले ही अन्य चैनल में अभी भी आइटम हों।

यह एक सामान्य पैटर्न है जहाँ तुम्हारे पास एक एकल संदर्भ फ़ाइल है जिसे कई नमूनों में पुन: उपयोग करने की आवश्यकता है। समाधान संदर्भ चैनल को एक value चैनल में बदलना है जिसे अनिश्चित काल तक पुन: उपयोग किया जा सकता है।

#### कोड ठीक करें

इसे संबोधित करने के कुछ तरीके हैं जो इस बात पर निर्भर करते हैं कि कितनी फ़ाइलें प्रभावित हैं।

**विकल्प 1**: तुम्हारे पास एक एकल संदर्भ फ़ाइल है जिसे तुम बहुत पुन: उपयोग कर रहे हो। तुम बस एक value चैनल प्रकार बना सकते हो, जिसे बार-बार उपयोग किया जा सकता है। इसे करने के तीन तरीके हैं:

**1a** `channel.value()` का उपयोग करें:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Value channel can be reused
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** `first()` [ऑपरेटर](https://www.nextflow.io/docs/latest/reference/operator.html#first) का उपयोग करें:

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** `collect()` [ऑपरेटर](https://www.nextflow.io/docs/latest/reference/operator.html#collect) का उपयोग करें:

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**विकल्प 2**: अधिक जटिल परिदृश्यों में, शायद जहाँ तुम्हारे पास नमूना चैनल में सभी नमूनों के लिए कई संदर्भ फ़ाइलें हैं, तुम `combine` ऑपरेटर का उपयोग कर सकते हो एक नया चैनल बनाने के लिए जो दो चैनलों को टपल में जोड़ता है:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Creates cartesian product

    PROCESS_FILES(combined_ch)
}
```

`.combine()` ऑपरेटर दो चैनलों का कार्टेशियन उत्पाद उत्पन्न करता है, इसलिए `reference_ch` में प्रत्येक आइटम `input_ch` में प्रत्येक आइटम के साथ जोड़ा जाएगा। यह प्रोसेस को प्रत्येक नमूने के लिए चलने की अनुमति देता है जबकि अभी भी संदर्भ का उपयोग करता है।

इसके लिए प्रोसेस इनपुट को समायोजित करने की आवश्यकता है। हमारे उदाहरण में, प्रोसेस परिभाषा की शुरुआत को निम्नानुसार समायोजित करने की आवश्यकता होगी:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

यह दृष्टिकोण सभी स्थितियों में उपयुक्त नहीं हो सकता है।

#### पाइपलाइन चलाएं

ऊपर दिए गए फिक्स में से एक को आज़माओ और वर्कफ़्लो को फिर से चलाओ:

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

अब तुम्हें केवल एक के बजाय सभी तीन नमूनों को प्रोसेस होते देखना चाहिए।

### 2.3. गलत चैनल सामग्री संरचना

जब वर्कफ़्लो जटिलता के एक निश्चित स्तर तक पहुंचते हैं, तो प्रत्येक चैनल की आंतरिक संरचनाओं का ट्रैक रखना थोड़ा मुश्किल हो सकता है, और लोग आमतौर पर प्रोसेस की अपेक्षा और चैनल में वास्तव में क्या है के बीच बेमेल उत्पन्न करते हैं। यह उस समस्या से अधिक सूक्ष्म है जिसकी हमने पहले चर्चा की थी, जहाँ चैनलों की संख्या गलत थी। इस मामले में, तुम्हारे पास सही संख्या में इनपुट चैनल हो सकते हैं, लेकिन उन चैनलों में से एक या अधिक की आंतरिक संरचना प्रोसेस की अपेक्षा से मेल नहीं खाती है।

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

#### कोड जांचें

एरर मैसेज में स्क्वायर ब्रैकेट यहाँ सुराग प्रदान करते हैं - प्रोसेस टपल को एक एकल मान के रूप में मान रहा है, जो हम नहीं चाहते हैं। चलो `bad_channel_shape.nf` की जांच करते हैं:

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

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

तुम देख सकते हो कि हम टपल से बना एक चैनल उत्पन्न कर रहे हैं: `['sample1', 'file1.txt']`, लेकिन प्रोसेस एक एकल मान की उम्मीद करता है, `val sample_name`। निष्पादित कमांड दिखाता है कि प्रोसेस `[sample3, file3.txt]_output.txt` नामक फ़ाइल बनाने की कोशिश कर रहा है, जो इच्छित आउटपुट नहीं है।

#### कोड ठीक करें

इसे ठीक करने के लिए, अगर प्रोसेस को दोनों इनपुट की आवश्यकता है तो हम प्रोसेस को टपल स्वीकार करने के लिए समायोजित कर सकते हैं:

=== "विकल्प 1: प्रोसेस में टपल स्वीकार करें"

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

            // Channel emits tuples, but process expects single values
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

            // Channel emits tuples, but process expects single values
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

            // Channel emits tuples, but process expects single values
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

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### पाइपलाइन चलाएं

समाधानों में से एक चुनो और वर्कफ़्लो को फिर से चलाओ:

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

### 2.4. चैनल डिबगिंग तकनीकें

#### चैनल निरीक्षण के लिए `.view()` का उपयोग करना

चैनलों के लिए सबसे शक्तिशाली डिबगिंग टूल `.view()` ऑपरेटर है। `.view()` के साथ, तुम डिबगिंग में मदद के लिए सभी चरणों में अपने चैनलों के आकार को समझ सकते हो।

#### पाइपलाइन चलाएं

इसे क्रिया में देखने के लिए `bad_channel_shape_viewed.nf` चलाओ:

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

#### कोड जांचें

चलो `bad_channel_shape_viewed.nf` की जांच करते हैं यह देखने के लिए कि `.view()` का उपयोग कैसे किया जाता है:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Channel emits tuples, but process expects single values
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

#### कोड ठीक करें

भविष्य में चैनल सामग्री को समझने के लिए `.view()` ऑपरेशन का अत्यधिक उपयोग करने से बचाने के लिए, मदद के लिए कुछ कमेंट जोड़ना उचित है:

```groovy title="bad_channel_shape_viewed.nf (with comments)" linenums="16" hl_lines="8 9"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

यह अधिक महत्वपूर्ण हो जाएगा क्योंकि तुम्हारे वर्कफ़्लो जटिलता में बढ़ते हैं और चैनल संरचना अधिक अपारदर्शी हो जाती है।

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

### सारांश

कई चैनल संरचना एरर मान्य Nextflow सिंटैक्स के साथ बनाई जा सकती हैं। तुम डेटा फ़्लो को समझकर, निरीक्षण के लिए `.view()` ऑपरेटर का उपयोग करके, और अप्रत्याशित टपल संरचनाओं को इंगित करने वाले स्क्वायर ब्रैकेट जैसे एरर मैसेज पैटर्न को पहचानकर चैनल संरचना एरर को डिबग कर सकते हो।

### आगे क्या है?

प्रोसेस परिभाषाओं द्वारा बनाई गई एरर के बारे में जानें।

---

## 3. प्रोसेस संरचना एरर

प्रोसेस से संबंधित अधिकांश एरर जो तुम्हें मिलेंगी वे कमांड बनाने में की गई गलतियों से संबंधित होंगी, या अंतर्निहित सॉफ़्टवेयर से संबंधित समस्याओं से। उस ने कहा, ऊपर चैनल समस्याओं के समान, तुम प्रोसेस परिभाषा में गलतियाँ कर सकते हो जो सिंटैक्स एरर के रूप में योग्य नहीं हैं, लेकिन जो रन टाइम पर एरर का कारण बनेंगी।

### 3.1. लापता आउटपुट फ़ाइलें

प्रोसेस लिखते समय एक सामान्य एरर कुछ ऐसा करना है जो प्रोसेस की अपेक्षा और जो उत्पन्न होता है के बीच बेमेल उत्पन्न करता है।

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

#### कोड जांचें

एरर मैसेज इंगित करता है कि प्रोसेस को `sample3.txt` नामक आउटपुट फ़ाइल उत्पन्न करने की उम्मीद थी, लेकिन स्क्रिप्ट वास्तव में `sample3_output.txt` बनाती है। चलो `missing_output.nf` में प्रोसेस परिभाषा की जांच करते हैं:

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

तुम्हें देखना चाहिए कि `output:` ब्लॉक में आउटपुट फ़ाइल नाम और स्क्रिप्ट में उपयोग किए गए के बीच बेमेल है। यह बेमेल प्रोसेस को विफल करने का कारण बनता है। अगर तुम्हें इस तरह की एरर मिलती है, तो वापस जाओ और जांचें कि आउटपुट तुम्हारी प्रोसेस परिभाषा और तुम्हारे आउटपुट ब्लॉक के बीच मेल खाते हैं।

अगर समस्या अभी भी स्पष्ट नहीं है, तो वास्तविक आउटपुट फ़ाइलों की पहचान करने के लिए वर्क डायरेक्टरी की जांच करो:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

इस उदाहरण के लिए यह हमें हाइलाइट करेगा कि आउटपुट फ़ाइल नाम में एक `_output` सफिक्स शामिल किया जा रहा है, हमारी `output:` परिभाषा के विपरीत।

#### कोड ठीक करें

आउटपुट फ़ाइलनाम को सुसंगत बनाकर बेमेल को ठीक करो:

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

### 3.2. लापता सॉफ़्टवेयर

एरर का एक और वर्ग सॉफ़्टवेयर प्रावधान में गलतियों के कारण होता है। `missing_software.nf` एक सिंटैक्टिक रूप से मान्य वर्कफ़्लो है, लेकिन यह कुछ बाहरी सॉफ़्टवेयर पर निर्भर करता है जो `cowpy` कमांड प्रदान करता है जिसका यह उपयोग करता है।

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

प्रोसेस के पास हम जो कमांड निर्दिष्ट कर रहे हैं उसकी पहुंच नहीं है। कभी-कभी यह इसलिए होता है क्योंकि एक स्क्रिप्ट वर्कफ़्लो `bin` डायरेक्टरी में मौजूद है, लेकिन निष्पादन योग्य नहीं बनाई गई है। अन्य समय यह इसलिए होता है क्योंकि सॉफ़्टवेयर कंटेनर या एनवायरनमेंट में इंस्टॉल नहीं है जहाँ वर्कफ़्लो चल रहा है।

#### कोड जांचें

उस `127` एक्ज़िट कोड पर ध्यान दो - यह तुम्हें बिल्कुल समस्या बताता है। चलो `missing_software.nf` की जांच करते हैं:

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

#### कोड ठीक करें

हम यहाँ थोड़े बेईमान रहे हैं, और वास्तव में कोड में कुछ भी गलत नहीं है। हमें बस प्रोसेस को इस तरह से चलाने के लिए आवश्यक कॉन्फ़िगरेशन निर्दिष्ट करने की आवश्यकता है कि इसके पास प्रश्न में कमांड तक पहुंच हो। इस मामले में प्रोसेस में एक कंटेनर परिभाषा है, इसलिए हमें बस Docker सक्षम के साथ वर्कफ़्लो चलाने की आवश्यकता है।

#### पाइपलाइन चलाएं

हमने तुम्हारे लिए `nextflow.config` में एक Docker प्रोफ़ाइल सेट अप की है, इसलिए तुम वर्कफ़्लो को इसके साथ चला सकते हो:

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

    Nextflow कंटेनर का उपयोग कैसे करता है इसके बारे में अधिक जानने के लिए, [Hello Nextflow](../hello_nextflow/05_hello_containers.md) देखें

### 3.3. खराब रिसोर्स कॉन्फ़िगरेशन

प्रोडक्शन उपयोग में, तुम अपनी प्रोसेस पर रिसोर्स कॉन्फ़िगर कर रहे होगे। उदाहरण के लिए `memory` तुम्हारी प्रोसेस के लिए उपलब्ध मेमोरी की अधिकतम मात्रा को परिभाषित करता है, और अगर प्रोसेस उससे अधिक हो जाती है, तो तुम्हारा शेड्यूलर आमतौर पर प्रोसेस को मार देगा और `137` का एक्ज़िट कोड लौटाएगा। हम यहाँ इसे प्रदर्शित नहीं कर सकते क्योंकि हम `local` एक्ज़ीक्यूटर का उपयोग कर रहे हैं, लेकिन हम `time` के साथ कुछ समान दिखा सकते हैं।

#### पाइपलाइन चलाएं

`bad_resources.nf` में 1 मिलीसेकंड के समय की अवास्तविक सीमा के साथ प्रोसेस कॉन्फ़िगरेशन है:

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

#### कोड जांचें

चलो `bad_resources.nf` की जांच करते हैं:

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

हम जानते हैं कि प्रोसेस एक सेकंड से अधिक समय लेगी (हमने यह सुनिश्चित करने के लिए वहाँ एक sleep जोड़ा है), लेकिन प्रोसेस 1 मिलीसेकंड के बाद टाइम आउट होने के लिए सेट है। किसी ने अपने कॉन्फ़िगरेशन के साथ थोड़ा अवास्तविक रहा है!

#### कोड ठीक करें

समय सीमा को एक यथार्थवादी मान तक बढ़ाओ:

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

अगर तुम सुनिश्चित करते हो कि तुम अपने एरर मैसेज पढ़ते हो तो इस तरह की विफलताएं तुम्हें बहुत लंबे समय तक परेशान नहीं करनी चाहिए। लेकिन सुनिश्चित करो कि तुम उन कमांड की रिसोर्स आवश्यकताओं को समझते हो जो तुम चला रहे हो ताकि तुम अपने रिसोर्स निर्देशों को उचित रूप से कॉन्फ़िगर कर सको।

### 3.4. प्रोसेस डिबगिंग तकनीकें

जब प्रोसेस विफल होती हैं या अप्रत्याशित रूप से व्यवहार करती हैं, तो तुम्हें यह जांचने के लिए व्यवस्थित तकनीकों की आवश्यकता होती है कि क्या गलत हुआ। वर्क डायरेक्टरी में वह सभी जानकारी होती है जो तुम्हें प्रोसेस एक्ज़ीक्यूशन को डिबग करने के लिए चाहिए।

#### वर्क डायरेक्टरी निरीक्षण का उपयोग करना

प्रोसेस के लिए सबसे शक्तिशाली डिबगिंग टूल वर्क डायरेक्टरी की जांच करना है। जब एक प्रोसेस विफल होती है, तो Nextflow उस विशिष्ट प्रोसेस एक्ज़ीक्यूशन के लिए एक वर्क डायरेक्टरी बनाता है जिसमें वे सभी फ़ाइलें होती हैं जो यह समझने के लिए आवश्यक हैं कि क्या हुआ।

#### पाइपलाइन चलाएं

चलो वर्क डायरेक्टरी निरीक्षण को प्रदर्शित करने के लिए पहले के `missing_output.nf` उदाहरण का उपयोग करते हैं (अगर तुम्हें आवश्यकता हो तो आउटपुट नामकरण बेमेल को फिर से उत्पन्न करो):

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

#### वर्क डायरेक्टरी जांचें

जब तुम्हें यह एरर मिलती है, तो वर्क डायरेक्टरी में सभी डिबगिंग जानकारी होती है। एरर मैसेज से वर्क डायरेक्टरी पथ खोजो और इसकी सामग्री की जांच करो:

```bash
# एरर मैसेज से वर्क डायरेक्टरी खोजें
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

तुम फिर मुख्य फ़ाइलों की जांच कर सकते हो:

##### कमांड स्क्रिप्ट जांचें

`.command.sh` फ़ाइल बिल्कुल दिखाती है कि कौन सी कमांड निष्पादित की गई थी:

```bash
# निष्पादित कमांड देखें
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

यह प्रकट करता है:

- **वेरिएबल प्रतिस्थापन**: क्या Nextflow वेरिएबल ठीक से विस्तारित किए गए थे
- **फ़ाइल पथ**: क्या इनपुट फ़ाइलें सही ढंग से स्थित थीं
- **कमांड संरचना**: क्या स्क्रिप्ट सिंटैक्स सही है

देखने के लिए सामान्य समस्याएं:

- **लापता उद्धरण**: स्पेस वाले वेरिएबल को उचित उद्धरण की आवश्यकता होती है
- **गलत फ़ाइल पथ**: इनपुट फ़ाइलें जो मौजूद नहीं हैं या गलत स्थानों में हैं
- **गलत वेरिएबल नाम**: वेरिएबल संदर्भों में टाइपो
- **लापता एनवायरनमेंट सेटअप**: कमांड जो विशिष्ट एनवायरनमेंट पर निर्भर करती हैं

##### एरर आउटपुट जांचें

`.command.err` फ़ाइल में वास्तविक एरर मैसेज होते हैं:

```bash
# एरर आउटपुट देखें
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

यह फ़ाइल दिखाएगी:

- **एक्ज़िट कोड**: 127 (कमांड नहीं मिली), 137 (मारा गया), आदि।
- **अनुमति एरर**: फ़ाइल एक्सेस समस्याएं
- **सॉफ़्टवेयर एरर**: एप्लिकेशन-विशिष्ट एरर मैसेज
- **रिसोर्स एरर**: मेमोरी/समय सीमा पार हो गई

##### स्टैंडर्ड आउटपुट जांचें

`.command.out` फ़ाइल दिखाती है कि तुम्हारी कमांड ने क्या उत्पन्न किया:

```bash
# स्टैंडर्ड आउटपुट देखें
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

यह सत्यापित करने में मदद करता है:

- **अपेक्षित आउटपुट**: क्या कमांड ने सही परिणाम उत्पन्न किए
- **आंशिक एक्ज़ीक्यूशन**: क्या कमांड शुरू हुई लेकिन बीच में विफल हो गई
- **डिबग जानकारी**: तुम्हारी स्क्रिप्ट से कोई डायग्नोस्टिक आउटपुट

##### एक्ज़िट कोड जांचें

`.exitcode` फ़ाइल में प्रोसेस के लिए एक्ज़िट कोड होता है:

```bash
# एक्ज़िट कोड देखें
cat work/*/*/.exitcode
```

सामान्य एक्ज़िट कोड और उनके अर्थ:

- **एक्ज़िट कोड 127**: कमांड नहीं मिली - सॉफ़्टवेयर इंस्टॉलेशन जांचें
- **एक्ज़िट कोड 137**: प्रोसेस मारी गई - मेमोरी/समय सीमा जांचें

##### फ़ाइल अस्तित्व जांचें

जब प्रोसेस लापता आउटपुट फ़ाइलों के कारण विफल होती हैं, तो जांचें कि वास्तव में कौन सी फ़ाइलें बनाई गईं:

```bash
# वर्क डायरेक्टरी में सभी फ़ाइलें सूचीबद्ध करें
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

यह पहचानने में मदद करता है:

- **फ़ाइल नामकरण बेमेल**: अपेक्षित से अलग नामों वाली आउटपुट फ़ाइलें
- **अनुमति समस्याएं**: फ़ाइलें जो नहीं बनाई जा सकीं
- **पथ समस्याएं**: गलत डायरेक्टरी में बनाई गई फ़ाइलें

हमारे पहले के उदाहरण में, इसने हमें पुष्टि की कि जबकि हमारी अपेक्षित `sample3.txt` मौजूद नहीं थी, `sample3_output.txt` थी:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### सारांश

प्रोसेस डिबगिंग के लिए यह समझने के लिए वर्क डायरेक्टरी की जांच करने की आवश्यकता होती है कि क्या गलत हुआ। मुख्य फ़ाइलों में `.command.sh` (निष्पादित स्क्रिप्ट), `.command.err` (एरर मैसेज), और `.command.out` (स्टैंडर्ड आउटपुट) शामिल हैं। 127 (कमांड नहीं मिली) और 137 (प्रोसेस मारी गई) जैसे एक्ज़िट कोड विफलता के प्रकार के बारे में तत्काल डायग्नोस्टिक सुराग प्रदान करते हैं।

### आगे क्या है?

Nextflow के बिल्ट-इन डिबगिंग टूल और समस्या निवारण के लिए व्यवस्थित दृष्टिकोण के बारे में जानें।

---

## 4. बिल्ट-इन डिबगिंग टूल और उन्नत तकनीकें

Nextflow डिबगिंग और वर्कफ़्लो एक्ज़ीक्यूशन का विश्लेषण करने के लिए कई शक्तिशाली बिल्ट-इन टूल प्रदान करता है। ये टूल तुम्हें यह समझने में मदद करते हैं कि क्या गलत हुआ, यह कहाँ गलत हुआ, और इसे कुशलता से कैसे ठीक किया जाए।

### 4.1. रियल-टाइम प्रोसेस आउटपुट

कभी-कभी तुम्हें यह देखने की आवश्यकता होती है कि चल रही प्रोसेस के अंदर क्या हो रहा है। तुम रियल-टाइम प्रोसेस आउटपुट सक्षम कर सकते हो, जो तुम्हें बिल्कुल दिखाता है कि प्रत्येक कार्य क्या कर रहा है जैसे यह निष्पादित होता है।

#### पाइपलाइन चलाएं

हमारे पहले के उदाहरणों से `bad_channel_shape_viewed.nf` ने `.view()` का उपयोग करके चैनल सामग्री प्रिंट की, लेकिन हम प्रोसेस के भीतर से ही वेरिएबल को इको करने के लिए `debug` निर्देश का भी उपयोग कर सकते हैं, जिसे हम `bad_channel_shape_viewed_debug.nf` में प्रदर्शित करते हैं। वर्कफ़्लो चलाओ:

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

#### कोड जांचें

चलो `bad_channel_shape_viewed_debug.nf` की जांच करते हैं यह देखने के लिए कि `debug` निर्देश कैसे काम करता है:

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

`debug` निर्देश एक प्रोसेस के एनवायरनमेंट को समझने का एक त्वरित और सुविधाजनक तरीका हो सकता है।

### 4.2. प्रीव्यू मोड

कभी-कभी तुम किसी भी प्रोसेस के चलने से पहले समस्याओं को पकड़ना चाहते हो। Nextflow इस तरह की सक्रिय डिबगिंग के लिए एक फ्लैग प्रदान करता है: `-preview`।

#### पाइपलाइन चलाएं

प्रीव्यू मोड तुम्हें कमांड निष्पादित किए बिना वर्कफ़्लो लॉजिक का परीक्षण करने देता है। यह तुम्हारे वर्कफ़्लो की संरचना की जल्दी जांच करने और यह सुनिश्चित करने के लिए काफी उपयोगी हो सकता है कि प्रोसेस बिना किसी वास्तविक कमांड को चलाए सही ढंग से जुड़े हुए हैं।

!!! note

    अगर तुमने पहले `bad_syntax.nf` को ठीक किया था, तो इस कमांड को चलाने से पहले स्क्रिप्ट ब्लॉक के बाद क्लोजिंग ब्रेस को हटाकर सिंटैक्स एरर को फिर से पेश करो।

यह कमांड चलाओ:

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

प्रीव्यू मोड किसी भी प्रोसेस को चलाए बिना सिंटैक्स एरर को जल्दी पकड़ने के लिए विशेष रूप से उपयोगी है। यह एक्ज़ीक्यूशन से पहले वर्कफ़्लो संरचना और प्रोसेस कनेक्शन को मान्य करता है।

### 4.3. लॉजिक टेस्टिंग के लिए स्टब रनिंग

कभी-कभी एरर को डिबग करना मुश्किल होता है क्योंकि कमांड बहुत लंबे समय तक चलती हैं, विशेष सॉफ़्टवेयर की आवश्यकता होती है, या जटिल कारणों से विफल होती हैं। स्टब रनिंग तुम्हें वास्तविक कमांड निष्पादित किए बिना वर्कफ़्लो लॉजिक का परीक्षण करने देती है।

#### पाइपलाइन चलाएं

जब तुम एक Nextflow प्रोसेस विकसित कर रहे होते हो, तो तुम 'डमी' कमांड को परिभाषित करने के लिए `stub` निर्देश का उपयोग कर सकते हो जो वास्तविक कमांड चलाए बिना सही रूप के आउटपुट उत्पन्न करते हैं। यह दृष्टिकोण विशेष रूप से मूल्यवान है जब तुम यह सत्यापित करना चाहते हो कि तुम्हारा वर्कफ़्लो लॉजिक वास्तविक सॉफ़्टवेयर की जटिलताओं से निपटने से पहले सही है।

उदाहरण के लिए, पहले से हमारे `missing_software.nf` को याद करो? वह जहाँ हमारे पास लापता सॉफ़्टवेयर था जिसने वर्कफ़्लो को चलने से रोका जब तक हमने `-profile docker` नहीं जोड़ा? `missing_software_with_stub.nf` एक बहुत समान वर्कफ़्लो है। अगर हम इसे उसी तरह चलाते हैं, तो हम वही एरर उत्पन्न करेंगे:

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

हालांकि, यह वर्कफ़्लो `-stub-run` के साथ चलाने पर एरर उत्पन्न नहीं करेगा, यहां तक कि `docker` प्रोफ़ाइल के बिना भी:

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

#### कोड जांचें

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

`missing_software.nf` के सापेक्ष, इस प्रोसेस में एक `stub:` निर्देश है जो `script:` में निर्दिष्ट के बजाय उपयोग की जाने वाली कमांड निर्दिष्ट करता है, इस घटना में कि Nextflow स्टब मोड में चलाया जाता है।

`touch` कमांड जिसका हम यहाँ उपयोग कर रहे हैं किसी भी सॉफ़्टवेयर या उपयुक्त इनपुट पर निर्भर नहीं करता है, और सभी स्थितियों में चलेगा, जिससे हमें प्रोसेस आंतरिक के बारे में चिंता किए बिना वर्कफ़्लो लॉजिक को डिबग करने की अनुमति मिलती है।

**स्टब रनिंग डिबग करने में मदद करती है:**

- चैनल संरचना और डेटा फ़्लो
- प्रोसेस कनेक्शन और निर्भरताएं
- पैरामीटर प्रचार
- सॉफ़्टवेयर निर्भरताओं के बिना वर्कफ़्लो लॉजिक

### 4.4. व्यवस्थित डिबगिंग दृष्टिकोण

अब जब तुमने व्यक्तिगत डिबगिंग तकनीकें सीख ली हैं - ट्रेस फ़ाइलों और वर्क डायरेक्टरी से लेकर प्रीव्यू मोड, स्टब रनिंग, और रिसोर्स मॉनिटरिंग तक - चलो उन्हें एक व्यवस्थित पद्धति में एक साथ बांधते हैं। एक संरचित दृष्टिकोण होने से तुम्हें जटिल एरर से अभिभूत होने से रोकता है और सुनिश्चित करता है कि```markdown
 तुम महत्वपूर्ण सुरागों को नहीं चूकते।

यह पद्धति हमने कवर किए गए सभी टूल को एक कुशल वर्कफ़्लो में जोड़ती है:

**चार-चरण डिबगिंग विधि:**

**चरण 1: सिंटैक्स एरर समाधान (5 मिनट)**

1. VSCode या तुम्हारे IDE में लाल अंडरलाइन जांचें
2. सिंटैक्स समस्याओं की पहचान करने के लिए `nextflow run workflow.nf -preview` चलाओ
3. सभी सिंटैक्स एरर ठीक करो (लापता ब्रेसेस, ट्रेलिंग कॉमा, आदि)
4. आगे बढ़ने से पहले सुनिश्चित करो कि वर्कफ़्लो सफलतापूर्वक पार्स होता है

**चरण 2: त्वरित मूल्यांकन (5 मिनट)**

1. रनटाइम एरर मैसेज को ध्यान से पढ़ो
2. जांचें कि यह रनटाइम, लॉजिक, या रिसोर्स एरर है
3. बुनियादी वर्कफ़्लो लॉजिक का परीक्षण करने के लिए प्रीव्यू मोड का उपयोग करो

**चरण 3: विस्तृत जांच (15-30 मिनट)**

1. विफल कार्य की वर्क डायरेक्टरी खोजो
2. लॉग फ़ाइलों की जांच करो
3. चैनलों का निरीक्षण करने के लिए `.view()` ऑपरेटर जोड़ो
4. एक्ज़ीक्यूशन के बिना वर्कफ़्लो लॉजिक का परीक्षण करने के लिए `-stub-run` का उपयोग करो

**चरण 4: ठीक करो और मान्य करो (15 मिनट)**

1. न्यूनतम लक्षित फिक्स करो
2. resume के साथ परीक्षण करो: `nextflow run workflow.nf -resume`
3. पूर्ण वर्कफ़्लो एक्ज़ीक्यूशन सत्यापित करो

!!! tip "कुशल डिबगिंग के लिए Resume का उपयोग करना"

    एक बार जब तुमने एक समस्या की पहचान कर ली है, तो तुम्हें अपने वर्कफ़्लो के सफल भागों को फिर से चलाने में समय बर्बाद किए बिना अपने फिक्स का परीक्षण करने का एक कुशल तरीका चाहिए। Nextflow की `-resume` कार्यक्षमता डिबगिंग के लिए अमूल्य है।

    अगर तुमने [Hello Nextflow](../hello_nextflow/) के माध्यम से काम किया है तो तुम `-resume` का सामना कर चुके होगे, और यह महत्वपूर्ण है कि तुम डिबगिंग करते समय अपने आप को समय बचाने के लिए इसका अच्छा उपयोग करो जबकि तुम्हारी समस्या प्रोसेस से पहले की प्रोसेस चलती हैं।

    **Resume डिबगिंग रणनीति:**

    1. विफलता तक वर्कफ़्लो चलाओ
    2. विफल कार्य के लिए वर्क डायरेक्टरी की जांच करो
    3. विशिष्ट समस्या को ठीक करो
    4. केवल फिक्स का परीक्षण करने के लिए resume करो
    5. वर्कफ़्लो पूर्ण होने तक दोहराओ

#### डिबगिंग कॉन्फ़िगरेशन प्रोफ़ाइल

इस व्यवस्थित दृष्टिकोण को और भी अधिक कुशल बनाने के लिए, तुम एक समर्पित डिबगिंग कॉन्फ़िगरेशन बना सकते हो जो स्वचालित रूप से सभी टूल को सक्षम करता है जिनकी तुम्हें आवश्यकता है:

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

फिर तुम इस प्रोफ़ाइल को सक्षम करके पाइपलाइन चला सकते हो:

```bash
nextflow run workflow.nf -profile debug
```

यह प्रोफ़ाइल रियल-टाइम आउटपुट सक्षम करती है, वर्क डायरेक्टरी को संरक्षित करती है, और आसान डिबगिंग के लिए समानांतरण को सीमित करती है।

### 4.5. व्यावहारिक डिबगिंग अभ्यास

अब व्यवस्थित डिबगिंग दृष्टिकोण को व्यवहार में लाने का समय है। वर्कफ़्लो `buggy_workflow.nf` में कई सामान्य एरर हैं जो उन प्रकार की समस्याओं का प्रतिनिधित्व करती हैं जिनका तुम वास्तविक दुनिया के विकास में सामना करोगे।

!!! exercise

    `buggy_workflow.nf` में सभी एरर की पहचान करने और उन्हें ठीक करने के लिए व्यवस्थित डिबगिंग दृष्टिकोण का उपयोग करो। यह वर्कफ़्लो CSV फ़ाइल से नमूना डेटा को प्रोसेस करने का प्रयास करता है लेकिन इसमें सामान्य डिबगिंग परिदृश्यों का प्रतिनिधित्व करने वाले कई जानबूझकर बग हैं।

    पहली एरर देखने के लिए वर्कफ़्लो चलाकर शुरू करो:

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

        यह रहस्यमय एरर `params{}` ब्लॉक में लाइन 11-12 के आसपास एक पार्सिंग समस्या को इंगित करती है। v2 पार्सर संरचनात्मक समस्याओं को जल्दी पकड़ता है।

    तुमने सीखी चार-चरण डिबगिंग विधि लागू करो:

    **चरण 1: सिंटैक्स एरर समाधान**
    - VSCode या तुम्हारे IDE में लाल अंडरलाइन जांचें
    - सिंटैक्स समस्याओं की पहचान करने के लिए `nextflow run workflow.nf -preview` चलाओ
    - सभी सिंटैक्स एरर ठीक करो (लापता ब्रेसेस, ट्रेलिंग कॉमा, आदि)
    - आगे बढ़ने से पहले सुनिश्चित करो कि वर्कफ़्लो सफलतापूर्वक पार्स होता है

    **चरण 2: त्वरित मूल्यांकन**
    - रनटाइम एरर मैसेज को ध्यान से पढ़ो
    - पहचानो कि एरर रनटाइम, लॉजिक, या रिसोर्स-संबंधित हैं या नहीं
    - बुनियादी वर्कफ़्लो लॉजिक का परीक्षण करने के लिए `-preview` मोड का उपयोग करो

    **चरण 3: विस्तृत जांच**
    - विफल कार्यों के लिए वर्क डायरेक्टरी की जांच करो
    - चैनलों का निरीक्षण करने के लिए `.view()` ऑपरेटर जोड़ो
    - वर्क डायरेक्टरी में लॉग फ़ाइलों की जांच करो
    - वर्कफ़्लो लॉजिक का एक्ज़ीक्यूशन के बिना परीक्षण करने के लिए `-stub-run` का उपयोग करो

    **चरण 4: ठीक करो और मान्य करो**
    - लक्षित फिक्स करो
    - फिक्स को कुशलता से परीक्षण करने के लिए `-resume` का उपयोग करो
    - पूर्ण वर्कफ़्लो एक्ज़ीक्यूशन सत्यापित करो

    **तुम्हारे निपटान में डिबगिंग टूल:**
    ```bash
    # सिंटैक्स जांच के लिए प्रीव्यू मोड
    nextflow run buggy_workflow.nf -preview

    # विस्तृत आउटपुट के लिए डिबग प्रोफ़ाइल
    nextflow run buggy_workflow.nf -profile debug

    # लॉजिक टेस्टिंग के लिए स्टब रनिंग
    nextflow run buggy_workflow.nf -stub-run

    # फिक्स के बाद Resume
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution
        `buggy_workflow.nf` में 9 या 10 अलग एरर हैं (तुम कैसे गिनते हो इस पर निर्भर करते हुए) जो सभी प्रमुख डिबगिंग श्रेणियों को कवर करती हैं। यहाँ प्रत्येक एरर और इसे कैसे ठीक करें का एक व्यवस्थित विवरण है

        चलो उन सिंटैक्स एरर से शुरू करते हैं:

        **एरर 1: सिंटैक्स एरर - ट्रेलिंग कॉमा**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // ERROR: Trailing comma
        ```
        **फिक्स:** ट्रेलिंग कॉमा हटाओ
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **एरर 2: सिंटैक्स एरर - लापता क्लोजिंग ब्रेस**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // ERROR: Missing closing brace for processFiles process
        ```
        **फिक्स:** लापता क्लोजिंग ब्रेस जोड़ो
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
        **फिक्स:** सही इनपुट वेरिएबल नाम का उपयोग करो
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **एरर 4: अपरिभाषित वेरिएबल एरर**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // ERROR: sample_ids undefined
        ```
        **फिक्स:** सही चैनल का उपयोग करो और नमूना ID निकालो
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        इस बिंदु पर वर्कफ़्लो चलेगा, लेकिन हमें अभी भी एरर मिल रही होंगी (जैसे `processFiles` में `Path value cannot be null`), खराब चैनल संरचना के कारण।

        **एरर 5: चैनल संरचना एरर - गलत Map आउटपुट**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // ERROR: processFiles expects tuple
        ```
        **फिक्स:** टपल संरचना लौटाओ जो processFiles की उम्मीद करता है
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        लेकिन यह ऊपर `heavyProcess()` चलाने के लिए हमारे फिक्स को तोड़ देगा, इसलिए हमें केवल नमूना ID को उस प्रोसेस में पास करने के लिए एक map का उपयोग करने की आवश्यकता होगी:

        **एरर 6: heavyProcess के लिए खराब चैनल संरचना**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // ERROR: input_ch now has 2 elements per emission- heavyProcess only needs 1 (the first)
        ```
        **फिक्स:** सही चैनल का उपयोग करो और नमूना ID निकालो
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        अब हम थोड़ा आगे बढ़ते हैं लेकिन `No such variable: i` के बारे में एक एरर प्राप्त करते हैं, क्योंकि हमने Bash वेरिएबल को एस्केप नहीं किया।

        **एरर 7: Bash वेरिएबल एस्केपिंग एरर**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // ERROR: $i not escaped
        ```
        **फिक्स:** bash वेरिएबल को एस्केप करो
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        अब हमें `Process exceeded running time limit (1ms)` मिलता है, इसलिए हम संबंधित प्रोसेस के लिए रन टाइम लिमिट को ठीक करते हैं:

        **एरर 8: रिसोर्स कॉन्फ़िगरेशन एरर**
        ```groovy linenums="36"
        time '1 ms'  // ERROR: Unrealistic time limit
        ```
        **फिक्स:** यथार्थवादी समय सीमा तक बढ़ाओ
        ```groovy linenums="36"
        time '100 s'
        ```

        अगला हमारे पास हल करने के लिए `Missing output file(s)` एरर है:

        **एरर 9: आउटपुट फ़ाइल नाम बेमेल**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // ERROR: Wrong filename, should match output declaration
        ```
        **फिक्स:** आउटपुट घोषणा से मेल खाएं
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        पहली दो प्रोसेस चलीं, लेकिन तीसरी नहीं।

        **एरर 10: आउटपुट फ़ाइल नाम बेमेल**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Error: attempting to take input from the pwd rather than a process
        handleFiles(file_ch)
        ```
        **फिक्स:** पिछली प्रोसेस से आउटपुट लो
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        इसके साथ, पूरा वर्कफ़्लो चलना चाहिए।

        **पूर्ण सही वर्कफ़्लो:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * डिबगिंग अभ्यास के लिए बगी वर्कफ़्लो
        * यह वर्कफ़्लो सीखने के उद्देश्यों के लिए कई जानबूझकर बग शामिल करता है
        */

        params{
            // Parameters with missing validation
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * इनपुट/आउटपुट बेमेल के साथ प्रोसेस
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
        * रिसोर्स समस्याओं के साथ प्रोसेस
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
            # भारी गणना का अनुकरण करें
            for i in {1..1000000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * फ़ाइल हैंडलिंग समस्याओं के साथ प्रोसेस
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
        * चैनल समस्याओं के साथ मुख्य वर्कफ़्लो
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

- **सिंटैक्स एरर**: लापता ब्रेसेस, ट्रेलिंग कॉमा, अपरिभाषित वेरिएबल
- **चैनल संरचना एरर**: गलत डेटा आकार, अपरिभाषित चैनल
- **प्रोसेस एरर**: आउटपुट फ़ाइल बेमेल, वेरिएबल एस्केपिंग
- **रिसोर्स एरर**: अवास्तविक समय सीमा

**मुख्य डिबगिंग सबक:**

1. **एरर मैसेज को ध्यान से पढ़ो** - वे अक्सर सीधे समस्या की ओर इशारा करते हैं
2. **व्यवस्थित दृष्टिकोण का उपयोग करो** - एक समय में एक एरर ठीक करो और `-resume` के साथ परीक्षण करो
3. **डेटा फ़्लो को समझो** - चैनल संरचना एरर अक्सर सबसे सूक्ष्म होती हैं
4. **वर्क डायरेक्टरी जांचें** - जब प्रोसेस विफल होती हैं, तो लॉग तुम्हें बिल्कुल बताते हैं कि क्या गलत हुआ

---

## सारांश

इस साइड क्वेस्ट में, तुमने Nextflow वर्कफ़्लो को डिबग करने के लिए व्यवस्थित तकनीकों का एक सेट सीखा है।
अपने स्वयं के काम में इन तकनीकों को लागू करने से तुम्हें अपने कंप्यूटर से लड़ने में कम समय बिताने, समस्याओं को तेजी से हल करने और भविष्य की समस्याओं से खुद को बचाने में सक्षम होगे।

### मुख्य पैटर्न

**1. सिंटैक्स एरर की पहचान और उन्हें कैसे ठीक करें**:

- Nextflow एरर मैसेज की व्याख्या करना और समस्याओं का पता लगाना
- सामान्य सिंटैक्स एरर: लापता ब्रेसेस, गलत कीवर्ड, अपरिभाषित वेरिएबल
- Nextflow (Groovy) और Bash वेरिएबल के बीच अंतर करना
- प्रारंभिक एरर डिटेक्शन के लिए VS Code एक्सटेंशन फीचर्स का उपयोग करना

```groovy
// Missing brace - look for red underlines in IDE
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- missing!

// Wrong keyword
inputs:  // Should be 'input:'

// Undefined variable - escape with backslash for Bash variables
echo "${undefined_var}"      // Nextflow variable (error if not defined)
echo "\${bash_var}"          // Bash variable (escaped)
```

**2. चैनल संरचना समस्याओं को कैसे डिबग करें**:

- चैनल कार्डिनैलिटी और एग्ज़ॉशन समस्याओं को समझना
- चैनल सामग्री संरचना बेमेल को डिबग करना
- चैनल निरीक्षण के लिए `.view()` ऑपरेटर का उपयोग करना
- आउटपुट में स्क्वायर ब्रैकेट जैसे एरर मैसेज पैटर्न को पहचानना जो अप्रत्याशित टपल संरचनाओं को इंगित करते हैं

```groovy
// Inspect channel content
my_channel.view { "Content: $it" }

// Convert queue to value channel (prevents exhaustion)
reference_ch = channel.value('ref.fa')
// or
reference_ch = channel.of('ref.fa').first()
```

**3. प्रोसेस एक्ज़ीक्यूशन समस्याओं को कैसे ट्रबलशूट करें**:

- लापता आउटपुट फ़ाइल एरर का निदान करना
- एक्ज़िट कोड को समझना (लापता सॉफ़्टवेयर के लिए 127, मेमोरी समस्याओं के लिए 137)
- वर्क डायरेक्टरी और कमांड फ़ाइलों की जांच करना
- रिसोर्स को उचित रूप से कॉन्फ़िगर करना

```bash
# जांचें कि वास्तव में क्या निष्पादित किया गया था
cat work/ab/cdef12/.command.sh

# एरर आउटपुट जांचें
cat work/ab/cdef12/.command.err

# एक्ज़िट कोड 127 = कमांड नहीं मिली
# एक्ज़िट कोड 137 = मारा गया (मेमोरी/समय सीमा)
```

**4. Nextflow के बिल्ट-इन डिबगिंग टूल का उपयोग कैसे करें**:

- प्रीव्यू मोड और रियल-टाइम डिबगिंग का लाभ उठाना
- लॉजिक टेस्टिंग के लिए स्टब रनिंग लागू करना
- कुशल डिबगिंग चक्रों के लिए resume लागू करना
- चार-चरण व्यवस्थित डिबगिंग पद्धति का पालन करना

!!! tip "त्वरित डिबगिंग संदर्भ"

    **सिंटैक्स एरर?** → VSCode चेतावनियां जांचें, `nextflow run workflow.nf -preview` चलाओ

    **चैनल समस्याएं?** → सामग्री का निरीक्षण करने के लिए `.view()` का उपयोग करो: `my_channel.view()`

    **प्रोसेस विफलताएं?** → वर्क डायरेक्टरी फ़ाइलें जांचें:

    - `.command.sh` - निष्पादित स्क्रिप्ट
    - `.command.err` - एरर मैसेज
    - `.exitcode` - एक्ज़िट स्टेटस (127 = कमांड नहीं मिली, 137 = मारा गया)

    **रहस्यमय व्यवहार?** → वर्कफ़्लो लॉजिक का परीक्षण करने के लिए `-stub-run` के साथ चलाओ

    **फिक्स किए?** → परीक्षण समय बचाने के लिए `-resume` का उपयोग करो: `nextflow run workflow.nf -resume`

---

### अतिरिक्त संसाधन

- [Nextflow troubleshooting guide](https://www.nextflow.io/docs/latest/troubleshooting.html): आधिकारिक समस्या निवारण डॉक्यूमेंटेशन
- [Understanding Nextflow channels](https://www.nextflow.io/docs/latest/channel.html): चैनल प्रकारों और व्यवहार में गहरी गोता
- [Process directives reference](https://www.nextflow.io/docs/latest/process.html#directives): सभी उपलब्ध प्रोसेस कॉन्फ़िगरेशन विकल्प
- [nf-test](https://www.nf-test.com/): Nextflow पाइपलाइन के लिए टेस्टिंग फ्रेमवर्क
- [Nextflow Slack community](https://www.nextflow.io/slack-invite.html): समुदाय से मदद प्राप्त करो

प्रोडक्शन वर्कफ़्लो के लिए, विचार करो:

- स्केल पर मॉनिटरिंग और डिबगिंग के लिए [Seqera Platform](https://seqera.io/platform/) सेट अप करना
- पुनरुत्पादनीय सॉफ़्टवेयर एनवायरनमेंट के लिए [Wave containers](https://seqera.io/wave/) का उपयोग करना

**याद रखो:** प्रभावी डिबगिंग एक कौशल है जो अभ्यास के साथ सुधरता है। व्यवस्थित पद्धति और व्यापक टूलकिट जो तुमने यहाँ प्राप्त की है, तुम्हारी पूरी Nextflow विकास यात्रा के दौरान तुम्हारी अच्छी सेवा करेगी।

---

## आगे क्या है?

[Side Quests के मेनू](./index.md) पर वापस जाओ या सूची में अगले विषय पर जाने के लिए पृष्ठ के निचले दाएं कोने में बटन पर क्लिक करो।
