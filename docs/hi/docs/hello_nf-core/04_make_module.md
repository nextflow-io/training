# भाग 4: nf-core मॉड्यूल बनाएं

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core प्रशिक्षण पाठ्यक्रम के इस चौथे भाग में, हम आपको दिखाएंगे कि मुख्य परंपराओं को लागू करके एक nf-core मॉड्यूल कैसे बनाया जाए जो मॉड्यूल को पोर्टेबल और मेंटेनेबल बनाती हैं।

nf-core प्रोजेक्ट एक कमांड (`nf-core modules create`) प्रदान करता है जो स्वचालित रूप से उचित रूप से संरचित मॉड्यूल टेम्पलेट जेनरेट करता है, जैसा कि हमने भाग 2 में workflow के लिए उपयोग किया था।
हालाँकि, शिक्षण उद्देश्यों के लिए, हम मैन्युअल रूप से शुरू करने जा रहे हैं: आपके `core-hello` पाइपलाइन में स्थानीय `cowpy` मॉड्यूल को चरण-दर-चरण nf-core-शैली के मॉड्यूल में बदलना।
इसके बाद, हम आपको दिखाएंगे कि भविष्य में अधिक कुशलता से काम करने के लिए टेम्पलेट-आधारित मॉड्यूल निर्माण का उपयोग कैसे करें।

??? info "इस सेक्शन से कैसे शुरू करें"

    यह सेक्शन मानता है कि आपने [भाग 3: nf-core मॉड्यूल का उपयोग करें](./03_use_module.md) पूरा कर लिया है और अपनी पाइपलाइन में `CAT_CAT` मॉड्यूल को एकीकृत कर लिया है।

    यदि आपने भाग 3 पूरा नहीं किया है या इस भाग के लिए नया शुरू करना चाहते हैं, तो आप अपने शुरुआती बिंदु के रूप में `core-hello-part3` समाधान का उपयोग कर सकते हैं।
    `hello-nf-core/` डायरेक्टरी के अंदर से ये कमांड चलाएं:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    यह आपको एक पाइपलाइन देता है जिसमें `CAT_CAT` मॉड्यूल पहले से एकीकृत है।
    आप निम्नलिखित कमांड चलाकर परीक्षण कर सकते हैं कि यह सफलतापूर्वक चलता है:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. `cowpy` को nf-core मॉड्यूल में बदलें

इस सेक्शन में, हम आपकी `core-hello` पाइपलाइन में स्थानीय `cowpy` मॉड्यूल पर nf-core परंपराओं को लागू करेंगे, इसे एक ऐसे मॉड्यूल में बदलेंगे जो nf-core समुदाय के मानकों का पालन करता है।

यह `cowpy` प्रोसेस मॉड्यूल के लिए वर्तमान कोड है:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// cowpy के साथ ASCII art जनरेट करें (https://github.com/jeffbuttars/cowpy)
process cowpy {

    publishDir 'results', mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """
}
```

हम निम्नलिखित nf-core परंपराओं को क्रमिक रूप से लागू करेंगे:

1. **प्रोसेस नाम को `COWPY` में अपरकेस करें** ताकि परंपरा का पालन हो सके।
2. **`COWPY` को मेटाडेटा टपल का उपयोग करने के लिए अपडेट करें** ताकि workflow के माध्यम से नमूना मेटाडेटा प्रसारित हो सके।
3. **`ext.args` के साथ टूल आर्गुमेंट कॉन्फ़िगरेशन को केंद्रीकृत करें** ताकि मॉड्यूल की बहुमुखी प्रतिभा बढ़े जबकि इंटरफ़ेस न्यूनतम रहे।
4. **`ext.prefix` के साथ आउटपुट नामकरण को मानकीकृत करें** ताकि स्थिरता को बढ़ावा मिले।
5. **प्रकाशन कॉन्फ़िगरेशन को केंद्रीकृत करें** ताकि स्थिरता को बढ़ावा मिले।

प्रत्येक चरण के बाद, हम यह परीक्षण करने के लिए पाइपलाइन चलाएंगे कि सब कुछ अपेक्षित रूप से काम कर रहा है।

!!! warning "वर्किंग डायरेक्टरी"

    सुनिश्चित करें कि आप इस सेक्शन में सभी फ़ाइल संपादनों और कमांड निष्पादन के लिए `core-hello` डायरेक्टरी (आपकी पाइपलाइन रूट) में हैं।

    ```bash
    cd core-hello
    ```

### 1.1. प्रोसेस नाम को अपरकेस करें

यह पूरी तरह से एक शैलीगत परंपरा है (कोई तकनीकी औचित्य नहीं है) लेकिन चूंकि यह nf-core मॉड्यूल के लिए मानक है, तो चलिए इसका पालन करते हैं।

हमें तीन सेट परिवर्तन करने की आवश्यकता है:

1. मॉड्यूल में प्रोसेस नाम अपडेट करें
2. workflow हेडर में मॉड्यूल इम्पोर्ट स्टेटमेंट अपडेट करें
3. workflow बॉडी में प्रोसेस कॉल और emit घोषणा अपडेट करें

चलिए शुरू करते हैं!

#### 1.1.1. मॉड्यूल में प्रोसेस नाम अपडेट करें

`cowpy.nf` मॉड्यूल फ़ाइल (`core-hello/modules/local/` के अंतर्गत) खोलें और प्रोसेस नाम को अपरकेस में संशोधित करें:

=== "बाद में"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // cowpy के साथ ASCII art जनरेट करें (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "पहले"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // cowpy के साथ ASCII art जनरेट करें (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

इस मामले में अपरकेसिंग पूरी तरह से सीधी है।

यदि प्रोसेस नाम कई शब्दों से बना होता, उदाहरण के लिए यदि हमारे पास मूल रूप से camel case में MyCowpyTool नामक एक प्रोसेस होता, तो nf-core परंपरा उन्हें अलग करने के लिए अंडरस्कोर का उपयोग करना होगी, जिससे MY_COWPY_TOOL मिलेगा।

#### 1.1.2. मॉड्यूल इम्पोर्ट स्टेटमेंट अपडेट करें

प्रोसेस नाम केस-सेंसिटिव हैं, इसलिए अब जब हमने प्रोसेस नाम बदल दिया है, हमें `hello.nf` के workflow हेडर में मॉड्यूल इम्पोर्ट स्टेटमेंट को तदनुसार अपडेट करने की आवश्यकता है:

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { cowpy                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

हम प्रोसेस कॉल को अपडेट करने से बचने के लिए इम्पोर्ट स्टेटमेंट में एक उपनाम का उपयोग कर सकते हैं, लेकिन इससे अपरकेसिंग परंपरा को अपनाने का उद्देश्य कुछ हद तक विफल हो जाएगा।

#### 1.1.3. प्रोसेस कॉल और emit घोषणा अपडेट करें

तो अब चलिए `hello.nf` के workflow ब्लॉक में प्रोसेस के दोनों संदर्भों को अपडेट करते हैं:

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
    COWPY(CAT_CAT.out.file_out)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    cowpy_hellos   = COWPY.out.cowpy_output
    versions       = ch_versions
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
    cowpy(CAT_CAT.out.file_out)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    cowpy_hellos   = cowpy.out.cowpy_output
    versions       = ch_versions
    ```

**दोनों** परिवर्तन करना सुनिश्चित करें, अन्यथा जब आप इसे चलाएंगे तो आपको एक त्रुटि मिलेगी।

#### 1.1.4. इसे परीक्षण करने के लिए पाइपलाइन चलाएं

चलिए इन परिवर्तनों के बाद यह परीक्षण करने के लिए workflow चलाते हैं कि सब कुछ सही तरीके से काम कर रहा है।

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [elegant_plateau] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2026-01-06_04-51-29

    Core Nextflow options
      runName                   : elegant_plateau
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [7b/66ceb5] CORE_HELLO:HELLO:sayHello (3)       | 3 of 3 ✔
    [8e/1bafb9] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [bb/203575] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [39/715489] CORE_HELLO:HELLO:COWPY              | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

बढ़िया, यह काम करता है! अब चलिए अधिक महत्वपूर्ण परिवर्तन करने की ओर बढ़ते हैं।

### 1.2. `COWPY` को मेटाडेटा टपल का उपयोग करने के लिए अपडेट करें

`core-hello` पाइपलाइन के वर्तमान संस्करण में, हम `COWPY` को पास करने के लिए `CAT_CAT` के आउटपुट टपल से फ़ाइल निकाल रहे हैं, जैसा कि नीचे दिए गए आरेख के ऊपरी भाग में दिखाया गया है।

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

यह बेहतर होगा कि `COWPY` सीधे मेटाडेटा टपल स्वीकार करे, जिससे मेटाडेटा workflow के माध्यम से प्रवाहित हो सके, जैसा कि आरेख के निचले भाग में दिखाया गया है।

इसके लिए, हमें निम्नलिखित परिवर्तन करने की आवश्यकता होगी:

1. इनपुट और आउटपुट परिभाषाओं को अपडेट करें
2. workflow में प्रोसेस कॉल को अपडेट करें
3. workflow में emit ब्लॉक को अपडेट करें

एक बार जब हम यह सब कर लेते हैं, तो हम यह परीक्षण करने के लिए पाइपलाइन चलाएंगे कि सब कुछ अभी भी पहले की तरह काम करता है।

#### 1.2.1. इनपुट और आउटपुट परिभाषाओं को अपडेट करें

`cowpy.nf` मॉड्यूल फ़ाइल पर वापस जाएं और इसे नीचे दिखाए अनुसार मेटाडेटा टपल स्वीकार करने के लिए संशोधित करें।

=== "बाद में"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output
    ```

=== "पहले"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"
    ```

जैसा कि आप देख सकते हैं, हमने **मुख्य इनपुट** और **आउटपुट** दोनों को एक टपल में बदल दिया जो भाग 3 में प्रस्तुत `tuple val(meta), path(input_file)` पैटर्न का पालन करता है।
आउटपुट के लिए, हमने आउटपुट चैनल को एक वर्णनात्मक नाम देने के लिए `emit: cowpy_output` जोड़ने का अवसर भी लिया।

अब जब हमने प्रोसेस की अपेक्षा बदल दी है, हमें प्रोसेस कॉल में जो हम प्रदान करते हैं उसे तदनुसार अपडेट करने की आवश्यकता है।

#### 1.2.2. workflow में प्रोसेस कॉल को अपडेट करें

अच्छी खबर यह है कि यह परिवर्तन प्रोसेस कॉल को सरल बना देगा।
अब जब `CAT_CAT` का आउटपुट और `COWPY` का इनपुट एक ही 'आकार' के हैं, यानी दोनों में `tuple val(meta), path(input_file)` संरचना है, तो हम उन्हें सीधे कनेक्ट कर सकते हैं बजाय `CAT_CAT` प्रोसेस के आउटपुट से फ़ाइल को स्पष्ट रूप से निकालने के।

`hello.nf` workflow फ़ाइल (`core-hello/workflows/` के अंतर्गत) खोलें और नीचे दिखाए अनुसार `COWPY` के कॉल को अपडेट करें।

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // tuple से file निकालें क्योंकि cowpy अभी metadata का उपयोग नहीं करता
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
        COWPY(ch_for_cowpy, params.character)
    ```

अब हम `CAT_CAT.out.file_out` पर सीधे `COWPY` को कॉल करते हैं।

परिणामस्वरूप, हमें अब `ch_for_cowpy` चैनल बनाने की आवश्यकता नहीं है, इसलिए उस लाइन (और इसकी टिप्पणी लाइन) को पूरी तरह से हटाया जा सकता है।

#### 1.2.3. workflow में emit ब्लॉक को अपडेट करें

चूंकि `COWPY` अब एक नामांकित आउटपुट, `cowpy_output`, उत्सर्जित करता है, हम `hello.nf` workflow के `emit:` ब्लॉक को उसका उपयोग करने के लिए अपडेट कर सकते हैं।

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out.cowpy_output
        versions       = ch_versions
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out
        versions       = ch_versions
    ```

तकनीकी रूप से यह आवश्यक नहीं है, लेकिन जब भी संभव हो नामांकित आउटपुट का संदर्भ देना अच्छा अभ्यास है।

#### 1.2.4. इसे परीक्षण करने के लिए पाइपलाइन चलाएं

चलिए इन परिवर्तनों के बाद यह परीक्षण करने के लिए workflow चलाते हैं कि सब कुछ सही तरीके से काम कर रहा है।

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [modest_saha] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-16-55

    Core Nextflow options
      runName                   : modest_saha
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [a8/447993] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [00/1fc59c] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [57/ac800d] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [b7/092f2b] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

पाइपलाइन सफलतापूर्वक चलनी चाहिए, मेटाडेटा अब `CAT_CAT` से `COWPY` के माध्यम से प्रवाहित हो रहा है।

यह पूरा करता है कि हमें `COWPY` को मेटाडेटा टपल संभालने के लिए क्या करने की आवश्यकता थी।
अब, चलिए देखते हैं कि हम nf-core मॉड्यूल पैटर्न का लाभ उठाने के लिए और क्या कर सकते हैं।

### 1.3. `ext.args` के साथ टूल आर्गुमेंट कॉन्फ़िगरेशन को केंद्रीकृत करें

अपनी वर्तमान स्थिति में, `COWPY` प्रोसेस `character` पैरामीटर के लिए एक मान प्राप्त करने की अपेक्षा करता है।
परिणामस्वरूप, हमें हर बार प्रोसेस को कॉल करते समय एक मान प्रदान करना होगा, भले ही हम टूल द्वारा सेट किए गए डिफ़ॉल्ट से खुश हों।
`COWPY` के लिए यह स्वीकार्य रूप से एक बड़ी समस्या नहीं है, लेकिन कई वैकल्पिक पैरामीटर वाले टूल के लिए, यह काफी बोझिल हो सकता है।

nf-core प्रोजेक्ट [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) नामक एक Nextflow फीचर का उपयोग करने की सिफारिश करता है ताकि कॉन्फ़िगरेशन फ़ाइलों के माध्यम से टूल आर्गुमेंट को अधिक सुविधाजनक रूप से प्रबंधित किया जा सके।

हर टूल विकल्प के लिए प्रोसेस इनपुट घोषित करने के बजाय, आप मॉड्यूल को अपनी कमांड लाइन के निर्माण में `ext.args` का संदर्भ देने के लिए लिखते हैं।
फिर यह केवल `modules.config` फ़ाइल में `ext.args` वेरिएबल को सेट करने की बात है, जो सभी मॉड्यूल के लिए कॉन्फ़िगरेशन विवरण को समेकित करती है, जिसमें आप उपयोग करना चाहते हैं आर्गुमेंट और मान रखते हैं।
Nextflow रनटाइम पर उन आर्गुमेंट को उनके मानों के साथ टूल कमांड लाइन में जोड़ देगा।

चलिए इस दृष्टिकोण को `COWPY` मॉड्यूल पर लागू करते हैं।
हमें निम्नलिखित परिवर्तन करने की आवश्यकता होगी:

1. `COWPY` मॉड्यूल को अपडेट करें
2. `modules.config` फ़ाइल में `ext.args` को कॉन्फ़िगर करें
3. `hello.nf` workflow को अपडेट करें

एक बार जब हम यह सब कर लेते हैं, तो हम यह परीक्षण करने के लिए पाइपलाइन चलाएंगे कि सब कुछ अभी भी पहले की तरह काम करता है।

#### 1.3.1. `COWPY` मॉड्यूल को अपडेट करें

चलिए इसे करते हैं।
`cowpy.nf` मॉड्यूल फ़ाइल (`core-hello/modules/local/` के अंतर्गत) खोलें और इसे नीचे दिखाए अनुसार `ext.args` का संदर्भ देने के लिए संशोधित करें।

=== "बाद में"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // cowpy के साथ ASCII art जनरेट करें (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

=== "पहले"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 20"
    #!/usr/bin/env nextflow

    // cowpy के साथ ASCII art जनरेट करें (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

आप देख सकते हैं कि हमने तीन परिवर्तन किए।

1. **`input:` ब्लॉक में, हमने `val character` इनपुट को हटा दिया।**
   आगे बढ़ते हुए, हम नीचे वर्णित `ext.args` कॉन्फ़िगरेशन के माध्यम से उस आर्गुमेंट को प्रदान करेंगे।

2. **`script:` ब्लॉक में, हमने लाइन `def args = task.ext.args ?: ''` जोड़ी।**
   वह लाइन `args` वेरिएबल का मान निर्धारित करने के लिए `?:` ऑपरेटर का उपयोग करती है: यदि यह खाली नहीं है तो `task.ext.args` की सामग्री, या यदि यह है तो एक खाली स्ट्रिंग।
   ध्यान दें कि जबकि हम आम तौर पर `ext.args` का संदर्भ देते हैं, इस कोड को मॉड्यूल-स्तरीय `ext.args` कॉन्फ़िगरेशन को बाहर निकालने के लिए `task.ext.args` का संदर्भ देना होगा।

3. **कमांड लाइन में, हमने `-c "$character"` को `$args` से बदल दिया।**
   यह वह जगह है जहाँ Nextflow `modules.config` फ़ाइल में `ext.args` में सेट किए गए किसी भी टूल आर्गुमेंट को इंजेक्ट करेगा।

परिणामस्वरूप, मॉड्यूल इंटरफ़ेस अब सरल है: यह केवल आवश्यक मेटाडेटा और फ़ाइल इनपुट की अपेक्षा करता है।

!!! note

    `?:` ऑपरेटर को अक्सर 'Elvis ऑपरेटर' कहा जाता है क्योंकि यह बगल से Elvis Presley के चेहरे की तरह दिखता है, `?` कैरेक्टर उनके बालों में लहर का प्रतीक है।

#### 1.3.2. `modules.config` फ़ाइल में `ext.args` को कॉन्फ़िगर करें

अब जब हमने मॉड्यूल से `character` घोषणा निकाल दी है, तो हमें इसे `modules.config` कॉन्फ़िगरेशन फ़ाइल में `ext.args` में जोड़ना होगा।

विशेष रूप से, हम `process {}` ब्लॉक में कोड का यह छोटा हिस्सा जोड़ने जा रहे हैं:

```groovy title="जोड़ने के लिए कोड"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

`withName:` सिंटैक्स इस कॉन्फ़िगरेशन को केवल `COWPY` प्रोसेस को असाइन करता है, और `ext.args = { "-c ${params.character}" }` बस एक स्ट्रिंग बनाता है जिसमें `character` पैरामीटर का मान शामिल होगा।
घुंघराले ब्रेसिज़ के उपयोग पर ध्यान दें, जो Nextflow को रनटाइम पर पैरामीटर का मान मूल्यांकन करने के लिए कहते हैं।

समझ में आया? चलिए इसे जोड़ते हैं।

`conf/modules.config` खोलें और नीचे दिखाए अनुसार `process {}` ब्लॉक के अंदर कॉन्फ़िगरेशन कोड जोड़ें।

=== "बाद में"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]

        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "पहले"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    ```

उम्मीद है कि आप कल्पना कर सकते हैं कि एक पाइपलाइन में सभी मॉड्यूल की अपनी `ext.args` इस फ़ाइल में निर्दिष्ट हों, निम्नलिखित लाभों के साथ:

- **मॉड्यूल इंटरफ़ेस सरल रहता है** - यह केवल आवश्यक मेटाडेटा और फ़ाइल इनपुट स्वीकार करता है
- **पाइपलाइन अभी भी `params.character` को उजागर करती है** - अंतिम उपयोगकर्ता अभी भी इसे पहले की तरह कॉन्फ़िगर कर सकते हैं
- **मॉड्यूल अब पोर्टेबल है** - इसे किसी विशिष्ट पैरामीटर नाम की अपेक्षा के बिना अन्य पाइपलाइनों में पुन: उपयोग किया जा सकता है
- कॉन्फ़िगरेशन `modules.config` में **केंद्रीकृत** है, workflow लॉजिक को साफ रखते हुए

`modules.config` फ़ाइल का उपयोग उस स्थान के रूप में करके जहाँ सभी पाइपलाइनें प्रति-मॉड्यूल कॉन्फ़िगरेशन को केंद्रीकृत करती हैं, हम अपने मॉड्यूल को विभिन्न पाइपलाइनों में अधिक पुन: उपयोग करने योग्य बनाते हैं।

#### 1.3.3. `hello.nf` workflow को अपडेट करें

चूंकि `COWPY` मॉड्यूल को अब इनपुट के रूप में `character` पैरामीटर की आवश्यकता नहीं है, हमें तदनुसार workflow कॉल को अपडेट करने की आवश्यकता है।

`hello.nf` workflow फ़ाइल (`core-hello/workflows/` के अंतर्गत) खोलें और नीचे दिखाए अनुसार `COWPY` के कॉल को अपडेट करें।

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
        COWPY(CAT_CAT.out.file_out)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

workflow कोड अब साफ है: हमें प्रोसेस को सीधे `params.character` पास करने की आवश्यकता नहीं है।
मॉड्यूल इंटरफ़ेस को न्यूनतम रखा गया है, जिससे यह अधिक पोर्टेबल हो जाता है, जबकि पाइपलाइन अभी भी कॉन्फ़िगरेशन के माध्यम से स्पष्ट विकल्प प्रदान करती है।

#### 1.3.4. इसे परीक्षण करने के लिए पाइपलाइन चलाएं

चलिए परीक्षण करते हैं कि workflow अभी भी अपेक्षित रूप से काम करता है, एक अलग कैरेक्टर निर्दिष्ट करके यह सत्यापित करें कि `ext.args` कॉन्फ़िगरेशन काम कर रहा है।

`kosh` का उपयोग करके यह कमांड चलाएं, जो अधिक... रहस्यमय विकल्पों में से एक है:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [exotic_planck] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-23-13

    Core Nextflow options
      runName                   : exotic_planck
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [13/9e3c0e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [e2/5b0ee5] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b6/4fb569] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [38/eb29ea] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

यह पहले की तरह सफलतापूर्वक चलना चाहिए।

चलिए सत्यापित करें कि `ext.args` कॉन्फ़िगरेशन ने काम किया आउटपुट की जाँच करके।
फ़ाइल ब्राउज़र में आउटपुट खोजें या आउटपुट फ़ाइल को देखने के लिए task hash (ऊपर के उदाहरण में `38/eb29ea` भाग) का उपयोग करें:

```bash
cat work/38/eb29ea*/cowpy-test.txt
```

??? success "कमांड आउटपुट"

    ```console
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
        \
        \
          \
      ___       _____     ___
    /   \     /    /|   /   \
    |     |   /    / |  |     |
    |     |  /____/  |  |     |
    |     |  |    |  |  |     |
    |     |  | {} | /   |     |
    |     |  |____|/    |     |
    |     |    |==|     |     |
    |      \___________/      |
    |                         |
    |                         |
    ```

आपको `kosh` कैरेक्टर के साथ ASCII art प्रदर्शित दिखनी चाहिए, जो पुष्टि करती है कि `ext.args` कॉन्फ़िगरेशन ने काम किया!

??? info "(वैकल्पिक) कमांड फ़ाइल का निरीक्षण करें"

    यदि आप यह देखना चाहते हैं कि कॉन्फ़िगरेशन को वास्तव में कैसे लागू किया गया था, तो आप `.command.sh` फ़ाइल का निरीक्षण कर सकते हैं:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    आप `-c kosh` आर्गुमेंट के साथ `cowpy` कमांड देखेंगे:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    यह दर्शाता है कि `.command.sh` फ़ाइल `ext.args` कॉन्फ़िगरेशन के आधार पर सही ढंग से जेनरेट की गई थी।

सोचने के लिए एक क्षण लें कि हमने यहाँ क्या हासिल किया।
यह दृष्टिकोण मॉड्यूल इंटरफ़ेस को आवश्यक डेटा (फ़ाइलें, मेटाडेटा, और कोई भी अनिवार्य प्रति-नमूना पैरामीटर) पर केंद्रित रखता है, जबकि टूल के व्यवहार को नियंत्रित करने वाले विकल्पों को कॉन्फ़िगरेशन के माध्यम से अलग से संभाला जाता है।

यह `cowpy` जैसे सरल टूल के लिए अनावश्यक लग सकता है, लेकिन यह डेटा विश्लेषण टूल के लिए बड़ा अंतर ला सकता है जिनमें बहुत सारे वैकल्पिक आर्गुमेंट हैं।

इस दृष्टिकोण के लाभों को सारांशित करने के लिए:

- **साफ इंटरफ़ेस**: मॉड्यूल आवश्यक डेटा इनपुट (मेटाडेटा और फ़ाइलें) पर केंद्रित है
- **लचीलापन**: उपयोगकर्ता कॉन्फ़िगरेशन के माध्यम से टूल आर्गुमेंट निर्दिष्ट कर सकते हैं, नमूना-विशिष्ट मानों सहित
- **स्थिरता**: सभी nf-core मॉड्यूल इस पैटर्न का पालन करते हैं
- **पोर्टेबिलिटी**: मॉड्यूल को हार्डकोडेड टूल विकल्पों के बिना पुन: उपयोग किया जा सकता है
- **कोई workflow परिवर्तन नहीं**: टूल विकल्पों को जोड़ने या बदलने के लिए workflow कोड को अपडेट करने की आवश्यकता नहीं है

!!! note

    `ext.args` सिस्टम में शक्तिशाली अतिरिक्त क्षमताएँ हैं जो यहाँ कवर नहीं की गई हैं, जिसमें मेटाडेटा के आधार पर आर्गुमेंट मानों को गतिशील रूप से स्विच करना शामिल है। अधिक विवरण के लिए [nf-core मॉड्यूल विशिष्टताएँ](https://nf-co.re/docs/guidelines/components/modules) देखें।

### 1.4. `ext.prefix` के साथ आउटपुट नामकरण को मानकीकृत करें

अब जब हमने `COWPY` प्रोसेस को metamap तक पहुँच दी है, तो हम एक और उपयोगी nf-core पैटर्न का लाभ उठाना शुरू कर सकते हैं: मेटाडेटा के आधार पर आउटपुट फ़ाइलों का नामकरण।

यहाँ हम `ext.prefix` नामक एक Nextflow फीचर का उपयोग करने जा रहे हैं जो हमें `meta.id` (metamap में शामिल पहचानकर्ता) का उपयोग करके मॉड्यूलों में आउटपुट फ़ाइल नामकरण को मानकीकृत करने की अनुमति देगा, जबकि अभी भी वांछित होने पर मॉड्यूल को व्यक्तिगत रूप से कॉन्फ़िगर करने में सक्षम होंगे।

यह `ext.args` के साथ हमने जो किया उसके समान होगा, कुछ अंतरों के साथ जिन्हें हम आगे बढ़ते हुए विस्तार से बताएंगे।

चलिए इस दृष्टिकोण को `COWPY` मॉड्यूल पर लागू करते हैं।
हमें निम्नलिखित परिवर्तन करने की आवश्यकता होगी:

1. `COWPY` मॉड्यूल को अपडेट करें
2. `modules.config` फ़ाइल में `ext.prefix` को कॉन्फ़िगर करें

(workflow में कोई परिवर्तन आवश्यक नहीं।)

एक बार जब हम यह कर लेते हैं, तो हम यह परीक्षण करने के लिए पाइपलाइन चलाएंगे कि सब कुछ अभी भी पहले की तरह काम करता है।

#### 1.4.1. `COWPY` मॉड्यूल को अपडेट करें

`cowpy.nf` मॉड्यूल फ़ाइल (`core-hello/modules/local/` के अंतर्गत) खोलें और इसे नीचे दिखाए अनुसार `ext.prefix` का संदर्भ देने के लिए संशोधित करें।

=== "बाद में"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 6 8"
        output:
            tuple val(meta), path("${prefix}.txt"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat $input_file | cowpy $args > ${prefix}.txt
        """
    }
    ```

=== "पहले"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 7"
        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

आप देख सकते हैं कि हमने तीन परिवर्तन किए।

1. **`script:` ब्लॉक में, हमने लाइन `prefix = task.ext.prefix ?: "${meta.id}"` जोड़ी।**
   वह लाइन `prefix` वेरिएबल का मान निर्धारित करने के लिए `?:` ऑपरेटर का उपयोग करती है: यदि यह खाली नहीं है तो `task.ext.prefix` की सामग्री, या यदि यह है तो metamap से पहचानकर्ता (`meta.id`)।
   ध्यान दें कि जबकि हम आम तौर पर `ext.prefix` का संदर्भ देते हैं, इस कोड को मॉड्यूल-स्तरीय `ext.prefix` कॉन्फ़िगरेशन को बाहर निकालने के लिए `task.ext.prefix` का संदर्भ देना होगा।

2. **कमांड लाइन में, हमने `cowpy-${input_file}` को `${prefix}.txt` से बदल दिया।**
   यह वह जगह है जहाँ Nextflow ऊपर की लाइन द्वारा निर्धारित `prefix` का मान इंजेक्ट करेगा।

3. **`output:` ब्लॉक में, हमने `path("cowpy-${input_file}")` को `path("${prefix}.txt")` से बदल दिया।**
   यह केवल यह दोहराता है कि कमांड लाइन में लिखे अनुसार फ़ाइल पथ क्या होगा।

परिणामस्वरूप, आउटपुट फ़ाइल नाम अब उचित फ़ाइल फ़ॉर्मेट एक्सटेंशन के साथ संयुक्त समझदार डिफ़ॉल्ट (metamap से पहचानकर्ता) का उपयोग करके बनाया गया है।

#### 1.4.2. `modules.config` फ़ाइल में `ext.prefix` को कॉन्फ़िगर करें

इस मामले में समझदार डिफ़ॉल्ट हमारे स्वाद के लिए पर्याप्त रूप से अभिव्यंजक नहीं है; हम एक कस्टम नामकरण पैटर्न का उपयोग करना चाहते हैं जिसमें टूल नाम शामिल हो, `cowpy-<id>.txt`, जैसा कि हमारे पास पहले था।

हम `modules.config` में `ext.prefix` को कॉन्फ़िगर करके ऐसा करेंगे, जैसा कि हमने `ext.args` के साथ `character` पैरामीटर के लिए किया था, सिवाय इस बार `withName: 'COWPY' {}` ब्लॉक पहले से मौजूद है, और हमें बस निम्नलिखित लाइन जोड़ने की आवश्यकता है:

```groovy title="जोड़ने के लिए कोड"
ext.prefix = { "cowpy-${meta.id}" }
```

यह वह स्ट्रिंग तैयार करेगा जो हम चाहते हैं।
ध्यान दें कि एक बार फिर हम घुंघराले ब्रेसिज़ का उपयोग करते हैं, इस बार Nextflow को रनटाइम पर `meta.id` के मान का मूल्यांकन करने के लिए कहने के लिए।

चलिए इसे जोड़ते हैं।

`conf/modules.config` खोलें और नीचे दिखाए अनुसार `process {}` ब्लॉक के अंदर कॉन्फ़िगरेशन कोड जोड़ें।

=== "बाद में"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "पहले"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    ```

यदि आप सोच रहे हैं, तो `ext.prefix` closure को सही मेटाडेटा तक पहुँच है क्योंकि कॉन्फ़िगरेशन को प्रोसेस निष्पादन के संदर्भ में मूल्यांकन किया जाता है, जहाँ मेटाडेटा उपलब्ध है।

#### 1.4.3. इसे परीक्षण करने के लिए पाइपलाइन चलाएं

चलिए परीक्षण करते हैं कि workflow अभी भी अपेक्षित रूप से काम करता है।

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [admiring_turing] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-29-02

    Core Nextflow options
      runName                   : admiring_turing
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [b2/e08524] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [13/88939f] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [23/4554e1] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [a3/c6cbe9] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

results डायरेक्टरी में आउटपुट देखें।
आपको पहले की तरह ही नामकरण के साथ cowpy आउटपुट फ़ाइल दिखनी चाहिए: `cowpy-test.txt`, डिफ़ॉल्ट batch नाम के आधार पर।

??? abstract "डायरेक्टरी सामग्री"

    ```console hl_lines="3"
    results
    ├── Bonjour-output.txt
    ├── cowpy-test.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

यह संतुष्ट करने के लिए स्वतंत्र महसूस करें कि आप मॉड्यूल या workflow कोड में कोई परिवर्तन किए बिना नामकरण पैटर्न को बदल सकते हैं, `conf/modules.config` में `ext.prefix` कॉन्फ़िगरेशन को बदलकर।

वैकल्पिक रूप से, आप यह भी परीक्षण कर सकते हैं कि वह भाग अभी भी तुरंत अनुकूलन योग्य है, कमांड लाइन पर एक अलग `--batch` पैरामीटर निर्दिष्ट करके इसे फिर से चलाकर।

यह दर्शाता है कि कैसे `ext.prefix` आपको अपनी पसंदीदा नामकरण परंपरा को बनाए रखने की अनुमति देता है जबकि मॉड्यूल इंटरफ़ेस को लचीला रखता है।

इस दृष्टिकोण के लाभों को सारांशित करने के लिए:

- **मानकीकृत नामकरण**: आउटपुट फ़ाइलें आमतौर पर मेटाडेटा से नमूना ID का उपयोग करके नामित की जाती हैं
- **कॉन्फ़िगर करने योग्य**: उपयोगकर्ता आवश्यकता पड़ने पर डिफ़ॉल्ट नामकरण को ओवरराइड कर सकते हैं
- **स्थिर**: सभी nf-core मॉड्यूल इस पैटर्न का पालन करते हैं
- **अनुमानित**: यह जानना आसान है कि आउटपुट फ़ाइलें क्या कहलाएंगी

काफी अच्छा, है ना?
खैर, हमें nf-core दिशानिर्देशों के अनुरूप अपने मॉड्यूल को सुधारने के लिए एक और महत्वपूर्ण परिवर्तन करने की आवश्यकता है।

### 1.5. प्रकाशन कॉन्फ़िगरेशन को केंद्रीकृत करें

आपने देखा होगा कि हम दो अलग-अलग डायरेक्टरी में आउटपुट प्रकाशित कर रहे हैं:

- **`results`** — मूल आउटपुट डायरेक्टरी जिसका हम शुरुआत से अपने स्थानीय मॉड्यूल के लिए उपयोग कर रहे हैं, प्रति-मॉड्यूल `publishDir` निर्देशों का उपयोग करके व्यक्तिगत रूप से सेट की गई;
- **`core-hello-results`** — कमांड लाइन पर `--outdir` के साथ सेट की गई आउटपुट डायरेक्टरी, जो nf-core लॉग और `CAT_CAT` द्वारा प्रकाशित परिणाम प्राप्त कर रही है।

यह गड़बड़ और उप-इष्टतम है; सब कुछ के लिए एक स्थान होना बेहतर होगा।
बेशक, हम अपने प्रत्येक स्थानीय मॉड्यूल में जा सकते हैं और `publishDir` निर्देश को मैन्युअल रूप से `core-hello-results` डायरेक्टरी का उपयोग करने के लिए अपडेट कर सकते हैं, लेकिन अगली बार जब हम आउटपुट डायरेक्टरी बदलने का निर्णय लेते हैं तो क्या होगा?

व्यक्तिगत मॉड्यूल को प्रकाशन निर्णय लेने देना स्पष्ट रूप से सही तरीका नहीं है, विशेष रूप से एक ऐसी दुनिया में जहाँ एक ही मॉड्यूल का उपयोग बहुत सारी विभिन्न पाइपलाइनों में किया जा सकता है, उन लोगों द्वारा जिनकी अलग-अलग आवश्यकताएँ या प्राथमिकताएँ हैं।
हम workflow कॉन्फ़िगरेशन के स्तर पर यह नियंत्रित करने में सक्षम होना चाहते हैं कि आउटपुट कहाँ प्रकाशित होते हैं।

"अरे," आप कह सकते हैं, "`CAT_CAT` अपने आउटपुट को `--outdir` में भेज रहा है। शायद हमें इसके `publishDir` निर्देश की नकल करनी चाहिए?"

हाँ, यह एक बढ़िया विचार है।

सिवाय इसके कि इसमें `publishDir` निर्देश नहीं है। (आगे बढ़ें, मॉड्यूल कोड देखें।)

ऐसा इसलिए है क्योंकि nf-core पाइपलाइनें व्यक्तिगत मॉड्यूल में बजाय `conf/modules.config` में `publishDir` को कॉन्फ़िगर करके workflow स्तर पर नियंत्रण को केंद्रीकृत करती हैं।
विशेष रूप से, nf-core टेम्पलेट एक डिफ़ॉल्ट `publishDir` निर्देश (एक पूर्वनिर्धारित डायरेक्टरी संरचना के साथ) घोषित करता है जो सभी मॉड्यूल पर लागू होता है जब तक कि एक ओवरराइडिंग निर्देश प्रदान नहीं किया जाता।

क्या यह अद्भुत नहीं लगता? क्या यह हो सकता है कि इस डिफ़ॉल्ट निर्देश का लाभ उठाने के लिए, हमें बस अपने स्थानीय मॉड्यूल से वर्तमान `publishDir` निर्देश को हटाने की आवश्यकता है?

चलिए `COWPY` पर इसे आज़माते हैं यह देखने के लिए कि क्या होता है, फिर हम यह समझने के लिए डिफ़ॉल्ट कॉन्फ़िगरेशन के लिए कोड देखेंगे कि यह कैसे काम करता है।

अंत में, हम प्रदर्शित करेंगे कि वांछित होने पर डिफ़ॉल्ट व्यवहार को कैसे ओवरराइड किया जाए।

#### 1.5.1. `COWPY` से `publishDir` निर्देश हटाएं

चलिए यह करते हैं।
`cowpy.nf` मॉड्यूल फ़ाइल (`core-hello/modules/local/` के अंतर्गत) खोलें और नीचे दिखाए अनुसार `publishDir` निर्देश हटाएं।

=== "बाद में"

    ```groovy title="core-hello/modules/local/cowpy.nf (excerpt)" linenums="1"
    #!/usr/bin/env nextflow

    // cowpy के साथ ASCII art जनरेट करें (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "पहले"

    ```groovy title="core-hello/modules/local/cowpy.nf (excerpt)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // cowpy के साथ ASCII art जनरेट करें (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

बस इतना ही!

#### 1.5.2. इसे परीक्षण करने के लिए पाइपलाइन चलाएं

चलिए देखते हैं कि अगर हम अब पाइपलाइन चलाते हैं तो क्या होता है।

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [silly_caravaggio] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-35-56

    Core Nextflow options
      runName                   : silly_caravaggio
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [db/39978e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [b5/bf6a8d] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b7/c61842] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [46/5839d6] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

अपनी वर्तमान वर्किंग डायरेक्टरी देखें।
अब `core-hello-results` में `COWPY` मॉड्यूल के आउटपुट भी शामिल हैं।

??? abstract "डायरेक्टरी सामग्री"

    ```console hl_lines="4-5"
    core-hello-results/
    ├── cat
    │   └── test.txt
    ├── cowpy
    │   └── cowpy-test.txt
    └── pipeline_info
        ├── execution_report_2025-12-27_06-16-55.html
        ├── execution_report_2025-12-27_06-23-13.html
        ├── execution_report_2025-12-27_06-29-02.html
        ├── execution_report_2025-12-27_06-35-56.html
        ├── execution_timeline_2025-12-27_06-16-55.html
        ├── execution_timeline_2025-12-27_06-23-13.html
        ├── execution_timeline_2025-12-27_06-29-02.html
        ├── execution_timeline_2025-12-27_06-35-56.html
        ├── execution_trace_2025-12-27_06-16-55.txt
        ├── execution_trace_2025-12-27_06-23-13.txt
        ├── execution_trace_2025-12-27_06-29-02.txt
        ├── execution_trace_2025-12-27_06-35-56.txt
        ├── hello_software_versions.yml
        ├── params_2025-12-27_06-17-00.json
        ├── params_2025-12-27_06-23-17.json
        ├── params_2025-12-27_06-29-07.json
        ├── params_2025-12-27_06-36-01.json
        ├── pipeline_dag_2025-12-27_06-16-55.html
        ├── pipeline_dag_2025-12-27_06-23-13.html
        ├── pipeline_dag_2025-12-27_06-29-02.html
        └── pipeline_dag_2025-12-27_06-35-56.html
    ```

आप देख सकते हैं कि Nextflow ने workflow और मॉड्यूल के नामों के आधार पर डायरेक्टरी की यह पदानुक्रम बनाई।

जिम्मेदार कोड `conf/modules.config` फ़ाइल में रहता है।
यह डिफ़ॉल्ट `publishDir` कॉन्फ़िगरेशन है जो nf-core टेम्पलेट का हिस्सा है और सभी प्रोसेस पर लागू होता है:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

यह जटिल लग सकता है, तो चलिए तीन घटकों में से प्रत्येक को देखते हैं:

- **`path:`** प्रोसेस नाम के आधार पर आउटपुट डायरेक्टरी निर्धारित करता है।
  `task.process` में निहित प्रोसेस का पूरा नाम workflow और मॉड्यूल इम्पोर्ट की पदानुक्रम शामिल करता है (जैसे `CORE_HELLO:HELLO:CAT_CAT`)।
  `tokenize` ऑपरेशन उस पदानुक्रम को हटाकर केवल प्रोसेस नाम प्राप्त करते हैं, फिर किसी भी अंडरस्कोर से पहले पहला भाग लेते हैं (यदि लागू हो), और इसे लोअरकेस में बदल देते हैं।
  यह वह है जो निर्धारित करता है कि `CAT_CAT` के परिणाम `${params.outdir}/cat/` में प्रकाशित होते हैं।
- **`mode:`** नियंत्रित करता है कि फ़ाइलें कैसे प्रकाशित की जाती हैं (copy, symlink, आदि)।
  यह `params.publish_dir_mode` पैरामीटर के माध्यम से कॉन्फ़िगर करने योग्य है।
- **`saveAs:`** फ़िल्टर करता है कि कौन सी फ़ाइलें प्रकाशित की जाएं।
  यह उदाहरण `versions.yml` फ़ाइलों को उनके लिए `null` लौटाकर बाहर करता है, उन्हें प्रकाशित होने से रोकता है।

यह आउटपुट को व्यवस्थित करने के लिए एक सुसंगत तर्क प्रदान करता है।

जब पाइपलाइन में सभी मॉड्यूल इस परंपरा को अपनाते हैं तो आउटपुट और भी बेहतर दिखता है, इसलिए अपनी पाइपलाइन में अन्य मॉड्यूल से `publishDir` निर्देशों को हटाने के लिए स्वतंत्र महसूस करें।
यह डिफ़ॉल्ट उन मॉड्यूल पर भी लागू होगा जिन्हें हमने nf-core दिशानिर्देशों का पालन करने के लिए स्पष्ट रूप से संशोधित नहीं किया।

उस ने कहा, आप तय कर सकते हैं कि आप अपने इनपुट को अलग तरीके से व्यवस्थित करना चाहते हैं, और अच्छी खबर यह है कि ऐसा करना आसान है।

#### 1.5.3. डिफ़ॉल्ट को ओवरराइड करें

डिफ़ॉल्ट `publishDir` निर्देश को ओवरराइड करने के लिए, आप बस `conf/modules.config` फ़ाइल में अपने स्वयं के निर्देश जोड़ सकते हैं।

उदाहरण के लिए, आप `withName:` सिलेक्टर का उपयोग करके एक एकल प्रोसेस के लिए डिफ़ॉल्ट को ओवरराइड कर सकते हैं, जैसा कि इस उदाहरण में है जहाँ हम 'COWPY' प्रोसेस के लिए एक कस्टम `publishDir` निर्देश जोड़ते हैं।

```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    withName: 'COWPY' {
        ext.args = { "-c ${params.character}" }
        publishDir = [
            path: 'my_custom_results'
        ]
    }
}
```

हम वास्तव में वह परिवर्तन नहीं करने जा रहे हैं, लेकिन इसके साथ खेलने और देखने के लिए स्वतंत्र महसूस करें कि आप क्या तर्क लागू कर सकते हैं।

बात यह है कि यह सिस्टम आपको दोनों दुनिया का सर्वश्रेष्ठ देता है: डिफ़ॉल्ट रूप से स्थिरता और मांग पर कॉन्फ़िगरेशन को अनुकूलित करने की लचीलापन।

सारांशित करने के लिए, आपको मिलता है:

- **सत्य का एकल स्रोत**: सभी प्रकाशन कॉन्फ़िगरेशन `modules.config` में रहता है
- **उपयोगी डिफ़ॉल्ट**: प्रोसेस प्रति-मॉड्यूल कॉन्फ़िगरेशन के बिना बॉक्स से बाहर काम करते हैं
- **आसान अनुकूलन**: मॉड्यूल कोड में नहीं, कॉन्फ़िग में प्रकाशन व्यवहार को ओवरराइड करें
- **पोर्टेबल मॉड्यूल**: मॉड्यूल आउटपुट स्थानों को हार्डकोड नहीं करते

यह nf-core मॉड्यूल फीचर्स के सेट को पूरा करता है जिन्हें आपको निश्चित रूप से उपयोग करना सीखना चाहिए, लेकिन अन्य भी हैं जिनके बारे में आप [nf-core मॉड्यूल विशिष्टताओं](https://nf-co.re/docs/guidelines/components/modules) में पढ़ सकते हैं।

### सारांश

अब आप जानते हैं कि nf-core परंपराओं का पालन करने के लिए स्थानीय मॉड्यूल को कैसे अनुकूलित किया जाए:

- अपने मॉड्यूल को मेटाडेटा टपल स्वीकार करने और प्रसारित करने के लिए डिज़ाइन करें;
- मॉड्यूल इंटरफ़ेस को न्यूनतम और पोर्टेबल रखने के लिए `ext.args` का उपयोग करें;
- कॉन्फ़िगर करने योग्य, मानकीकृत आउटपुट फ़ाइल नामकरण के लिए `ext.prefix` का उपयोग करें;
- एक सुसंगत परिणाम डायरेक्टरी संरचना के लिए डिफ़ॉल्ट केंद्रीकृत `publishDir` निर्देश को अपनाएं।

### आगे क्या है?

जानें कि आसान तरीके से मॉड्यूल बनाने के लिए nf-core के बिल्ट-इन टेम्पलेट-आधारित टूल का उपयोग कैसे करें।

---

## 2. nf-core टूलिंग के साथ एक मॉड्यूल बनाएं

अब जब आपने उन्हें मैन्युअल रूप से लागू करके nf-core मॉड्यूल पैटर्न सीख लिए हैं, तो चलिए देखते हैं कि आप व्यवहार में मॉड्यूल कैसे बनाएंगे।

### 2.1. टेम्पलेट से एक मॉड्यूल scaffold जेनरेट करें

पाइपलाइन बनाने के लिए जो मौजूद है उसके समान, nf-core प्रोजेक्ट एक टेम्पलेट के आधार पर उचित रूप से संरचित मॉड्यूल जेनरेट करने के लिए टूलिंग प्रदान करता है, जिसमें ये सभी पैटर्न शुरुआत से ही बनाए गए हैं।

#### 2.1.1. मॉड्यूल निर्माण कमांड चलाएं

`nf-core modules create` कमांड एक मॉड्यूल टेम्पलेट जेनरेट करता है जो पहले से ही उन सभी परंपराओं का पालन करता है जो आपने सीखी हैं।

चलिए एक न्यूनतम टेम्पलेट के साथ `COWPY` मॉड्यूल का एक नया संस्करण बनाते हैं यह कमांड चलाकर:

```bash
nf-core modules create --empty-template COWPY
```

`--empty-template` फ्लैग अतिरिक्त कोड के बिना एक साफ स्टार्टर टेम्पलेट बनाता है, जिससे आवश्यक संरचना को देखना आसान हो जाता है।

कमांड इंटरैक्टिव रूप से चलता है, आपको सेटअप के माध्यम से मार्गदर्शन करता है।
यह स्वचालित रूप से Bioconda और bio.tools जैसे पैकेज रिपॉजिटरी से टूल जानकारी खोजता है ताकि मेटाडेटा को पूर्व-भरा जा सके।

आपको कई कॉन्फ़िगरेशन विकल्पों के लिए प्रॉम्प्ट किया जाएगा:

- **लेखक जानकारी**: एट्रिब्यूशन के लिए आपका GitHub उपयोगकर्ता नाम
- **संसाधन लेबल**: कम्प्यूटेशनल आवश्यकताओं का एक पूर्वनिर्धारित सेट।
  nf-core प्रोजेक्ट हल्के टूल के लिए `process_single` और मांग वाले टूल के लिए `process_high` जैसे मानक लेबल प्रदान करता है।
  ये लेबल विभिन्न निष्पादन वातावरणों में संसाधन आवंटन को प्रबंधित करने में मदद करते हैं।
- **मेटाडेटा आवश्यकता**: क्या मॉड्यूल को `meta` map के माध्यम से नमूना-विशिष्ट जानकारी की आवश्यकता है (आमतौर पर डेटा प्रोसेसिंग मॉड्यूल के लिए हाँ)।

टूल पैकेज जानकारी खोजने और संरचना सेट करने की जटिलता को संभालता है, जिससे आप टूल के विशिष्ट तर्क को लागू करने पर ध्यान केंद्रित कर सकते हैं।

#### 2.1.2. मॉड्यूल scaffold की जाँच करें

टूल `modules/local/` में एक पूर्ण मॉड्यूल संरचना बनाता है (या `modules/nf-core/` यदि आप nf-core/modules रिपॉजिटरी में हैं):

??? abstract "डायरेक्टरी सामग्री"

    ```console
    modules/local/cowpy
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
    ```

प्रत्येक फ़ाइल एक विशिष्ट उद्देश्य पूरा करती है:

- **`main.nf`**: सभी nf-core पैटर्न के साथ प्रोसेस परिभाषा बनाई गई
- **`meta.yml`**: इनपुट, आउटपुट और टूल का वर्णन करने वाला मॉड्यूल दस्तावेज़ीकरण
- **`environment.yml`**: निर्भरताओं के लिए Conda वातावरण विशिष्टता
- **`tests/main.nf.test`**: मॉड्यूल के काम को मान्य करने के लिए nf-test टेस्ट केस

!!! tip "परीक्षण के बारे में अधिक जानें"

    जेनरेट की गई टेस्ट फ़ाइल nf-test का उपयोग करती है, जो Nextflow पाइपलाइनों और मॉड्यूल के लिए एक परीक्षण फ्रेमवर्क है। इन परीक्षणों को लिखने और चलाने का तरीका जानने के लिए, [nf-test साइड क्वेस्ट](../side_quests/nf-test.md) देखें।

जेनरेट की गई `main.nf` में वे सभी पैटर्न शामिल हैं जो आपने अभी सीखे हैं, साथ ही कुछ अतिरिक्त फीचर:

```groovy title="modules/local/cowpy/main.nf" hl_lines="11 21 22"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(input)        // पैटर्न 1: मेटाडेटा टपल ✓

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                // पैटर्न 2: ext.args ✓
    def prefix = task.ext.prefix ?: "${meta.id}"  // पैटर्न 3: ext.prefix ✓

    """
    // यहाँ अपना टूल कमांड जोड़ें

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
}
```

ध्यान दें कि ऊपर आपने मैन्युअल रूप से लागू किए गए सभी पैटर्न पहले से ही वहाँ हैं!

टेम्पलेट में कई अतिरिक्त nf-core परंपराएँ भी शामिल हैं।
इनमें से कुछ बॉक्स से बाहर काम करते हैं, जबकि अन्य प्लेसहोल्डर हैं जिन्हें हमें भरने की आवश्यकता होगी, जैसा कि नीचे वर्णित है।

**फीचर जो जैसे हैं वैसे काम करते हैं:**

- **`tag "$meta.id"`**: आसान ट्रैकिंग के लिए लॉग में प्रोसेस नामों में नमूना ID जोड़ता है
- **`label 'process_single'`**: CPU/मेमोरी आवश्यकताओं को कॉन्फ़िगर करने के लिए संसाधन लेबल
- **`when:` ब्लॉक**: `task.ext.when` कॉन्फ़िगरेशन के माध्यम से सशर्त निष्पादन की अनुमति देता है

ये फीचर पहले से ही कार्यात्मक हैं और मॉड्यूल को अधिक मेंटेनेबल बनाते हैं।

**प्लेसहोल्डर जिन्हें हम नीचे अनुकूलित करेंगे:**

- **`input:` और `output:` ब्लॉक**: सामान्य घोषणाएँ जिन्हें हम अपने टूल से मेल खाने के लिए अपडेट करेंगे
- **`script:` ब्लॉक**: एक टिप्पणी शामिल है जहाँ हम `cowpy` कमांड जोड़ेंगे
- **`stub:` ब्लॉक**: टेम्पलेट जिसे हम सही आउटपुट उत्पन्न करने के लिए अपडेट करेंगे
- **कंटेनर और वातावरण**: प्लेसहोल्डर जिन्हें हम पैकेज जानकारी से भरेंगे

अगले सेक्शन इन अनुकूलन को पूरा करने के माध्यम से चलते हैं।

### 2.2. कंटेनर और conda वातावरण सेट करें

nf-core दिशानिर्देशों की आवश्यकता है कि हम मॉड्यूल के हिस्से के रूप में एक कंटेनर और एक Conda वातावरण दोनों निर्दिष्ट करें।

#### 2.2.1. कंटेनर

कंटेनर के लिए, आप किसी भी Conda पैकेज से स्वचालित रूप से एक कंटेनर बनाने के लिए [Seqera Containers](https://seqera.io/containers/) का उपयोग कर सकते हैं, जिसमें conda-forge पैकेज शामिल हैं।
इस मामले में हम पहले की तरह ही पूर्व-निर्मित कंटेनर का उपयोग कर रहे हैं।

डिफ़ॉल्ट कोड Docker और Singularity के बीच टॉगल करने की पेशकश करता है, लेकिन हम उस लाइन को सरल बनाने जा रहे हैं और बस Docker कंटेनर निर्दिष्ट करेंगे जो हमें ऊपर Seqera Containers से मिला।

=== "बाद में"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273"
```

=== "पहले"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"
```

#### 2.2.2. Conda वातावरण

Conda वातावरण के लिए, मॉड्यूल कोड `conda "${moduleDir}/environment.yml"` निर्दिष्ट करता है जिसका अर्थ है कि इसे `environment.yml` फ़ाइल में कॉन्फ़िगर किया जाना चाहिए।

मॉड्यूल निर्माण टूल ने हमें चेतावनी दी कि यह Bioconda (बायोइन्फॉर्मेटिक्स टूल के लिए प्राथमिक चैनल) में `cowpy` पैकेज नहीं ढूंढ सका।
हालाँकि, `cowpy` conda-forge में उपलब्ध है, इसलिए आप `environment.yml` को इस तरह पूरा कर सकते हैं:

=== "बाद में"

    ```yaml title="modules/local/cowpy/environment.yml"  linenums="1" hl_lines="1 3 5"
    name: COWPY
    channels:
      - conda-forge
    dependencies:
      - cowpy=1.1.5
    ```

=== "पहले"

    ```yaml title="modules/local/cowpy/environment.yml" linenums="1"
    ---
    # yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
    channels:
      - conda-forge
      - bioconda
    dependencies:
      # TODO nf-core: List required Conda package(s).
      #               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
      #               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
      - "YOUR-TOOL-HERE"
    ```

nf-core में सबमिशन के लिए, हमें डिफ़ॉल्ट का अधिक बारीकी से पालन करना होगा, लेकिन अपने स्वयं के उपयोग के लिए हम इस तरह से कोड को सरल बना सकते हैं।

!!! tip "Bioconda बनाम conda-forge पैकेज"

    - **Bioconda पैकेज**: स्वचालित रूप से BioContainers बनाए जाते हैं, उपयोग के लिए तैयार कंटेनर प्रदान करते हैं
    - **conda-forge पैकेज**: Conda रेसिपी से मांग पर कंटेनर बनाने के लिए Seqera Containers का उपयोग कर सकते हैं

    अधिकांश बायोइन्फॉर्मेटिक्स टूल Bioconda में हैं, लेकिन conda-forge टूल के लिए, Seqera Containers कंटेनरीकरण के लिए एक आसान समाधान प्रदान करता है।

### 2.3. `COWPY` तर्क में प्लग करें

अब चलिए उन कोड तत्वों को अपडेट करते हैं जो `COWPY` प्रोसेस क्या करता है के लिए विशिष्ट हैं: इनपुट और आउटपुट, और script ब्लॉक।

#### 2.3.1. इनपुट और आउटपुट

जेनरेट किया गया टेम्पलेट सामान्य इनपुट और आउटपुट घोषणाएँ शामिल करता है जिन्हें आपको अपने विशिष्ट टूल के लिए अनुकूलित करने की आवश्यकता होगी।
सेक्शन 1 से हमारे मैन्युअल `COWPY` मॉड्यूल को देखते हुए, हम उसे एक गाइड के रूप में उपयोग कर सकते हैं।

इनपुट और आउटपुट ब्लॉक अपडेट करें:

=== "बाद में"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: cowpy_output
    path "versions.yml"           , emit: versions
    ```

=== "पहले"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions
    ```

यह निर्दिष्ट करता है:

- इनपुट फ़ाइल पैरामीटर नाम (`input_file` सामान्य `input` के बजाय)
- कॉन्फ़िगर करने योग्य prefix पैटर्न का उपयोग करके आउटपुट फ़ाइल नाम (`${prefix}.txt` वाइल्डकार्ड `*` के बजाय)
- एक वर्णनात्मक emit नाम (`cowpy_output` सामान्य `output` के बजाय)

यदि आप सिंटैक्स को मान्य करने के लिए Nextflow भाषा सर्वर का उपयोग कर रहे हैं, तो `${prefix}` भाग को इस स्तर पर एक त्रुटि के रूप में फ्लैग किया जाएगा क्योंकि हमने इसे अभी तक script ब्लॉक में नहीं जोड़ा है।
चलिए अब उस पर पहुँचते हैं।

#### 2.3.2. script ब्लॉक

टेम्पलेट script ब्लॉक में एक टिप्पणी प्लेसहोल्डर प्रदान करता है जहाँ आपको वास्तविक टूल कमांड जोड़नी चाहिए।

हमारे द्वारा पहले मैन्युअल रूप से लिखे गए मॉड्यूल के आधार पर, हमें निम्नलिखित संपादन करने चाहिए:

=== "बाद में"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="3 6"
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat $input_file | cowpy $args > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "पहले"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="6"
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    // यहाँ अपना टूल कमांड जोड़ें

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

मुख्य परिवर्तन:

- `def prefix` को केवल `prefix` में बदलें (`def` के बिना) ताकि इसे आउटपुट ब्लॉक में सुलभ बनाया जा सके
- टिप्पणी को वास्तविक `cowpy` कमांड से बदलें जो `$args` और `${prefix}.txt` दोनों का उपयोग करता है

ध्यान दें कि यदि हमने पहले से ही `COWPY` प्रोसेस के लिए `ext.args` और `ext.prefix` कॉन्फ़िगरेशन को `modules.config` फ़ाइल में जोड़ने का काम नहीं किया होता, तो हमें अब ऐसा करने की आवश्यकता होगी।

#### 2.3.3. stub ब्लॉक को लागू करना

Nextflow संदर्भ में, एक [stub](https://www.nextflow.io/docs/latest/process.html#stub) ब्लॉक आपको वास्तविक कमांड को निष्पादित किए बिना पाइपलाइन के तर्क के तेज़ प्रोटोटाइपिंग और परीक्षण के लिए उपयोग की जाने वाली एक हल्की, डमी स्क्रिप्ट परिभाषित करने की अनुमति देता है।

<!-- TODO (future) This is super glossed over but should really be explained or at least link out to an explanation about stubs (the reference doc isn't terribly helpful either). Right now this is likely to be mostly meaningless to anyone who doesn't already know about stubs. -->

यदि यह रहस्यमय लगता है तो बहुत चिंता न करें; हम इसे पूर्णता के लिए शामिल करते हैं लेकिन यदि आप इससे निपटना नहीं चाहते हैं तो आप stub सेक्शन को भी हटा सकते हैं, क्योंकि यह पूरी तरह से वैकल्पिक है।

=== "बाद में"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "पहले"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

मुख्य परिवर्तन:

- script ब्लॉक से मेल खाने के लिए `def prefix` को केवल `prefix` में बदलें
- `echo $args` लाइन हटाएं (जो केवल टेम्पलेट प्लेसहोल्डर कोड था)
- stub एक खाली `${prefix}.txt` फ़ाइल बनाता है जो script ब्लॉक द्वारा उत्पादित की जाती है

यह आपको वास्तविक टूल के चलने की प्रतीक्षा किए बिना workflow तर्क और फ़ाइल हैंडलिंग का परीक्षण करने की अनुमति देता है।

एक बार जब आप वातावरण सेटअप (सेक्शन 2.2), इनपुट/आउटपुट (सेक्शन 2.3.1), script ब्लॉक (सेक्शन 2.3.2), और stub ब्लॉक (सेक्शन 2.3.3) पूरा कर लेते हैं, तो मॉड्यूल परीक्षण के लिए तैयार है!

### 2.4. नए `COWPY` मॉड्यूल में स्वैप करें और पाइपलाइन चलाएं

इस नए संस्करण के `COWPY` मॉड्यूल को आज़माने के लिए हमें बस `hello.nf` workflow फ़ाइल में इम्पोर्ट स्टेटमेंट को नई फ़ाइल की ओर इंगित करने के लिए स्विच करना है।

=== "बाद में"

    ```groovy title="workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "पहले"

    ```groovy title="modules/local/cowpy/main.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

चलिए इसे परीक्षण करने के लिए पाइपलाइन चलाते हैं।

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "कमांड आउटपुट"

    ```console hl_lines="33"
      N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [prickly_neumann] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_08-23-51

    Core Nextflow options
      runName                   : prickly_neumann
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [e9/008ede] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [f0/d70cfe] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [be/0ecc58] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [11/8e082f] CORE_HELLO:HELLO:COWPY (test)       [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

यह पहले की तरह ही परिणाम उत्पन्न करता है।

### सारांश

अब आप जानते हैं कि सब कुछ शुरुआत से लिखने के बजाय टेम्पलेट का उपयोग करके कुशलता से मॉड्यूल बनाने के लिए बिल्ट-इन nf-core टूलिंग का उपयोग कैसे करें।

### आगे क्या है?

जानें कि nf-core में मॉड्यूल योगदान करने के क्या लाभ हैं और मुख्य चरण और आवश्यकताएँ क्या शामिल हैं।

---

## 3. nf-core में मॉड्यूल वापस योगदान करना

[nf-core/modules](https://github.com/nf-core/modules) रिपॉजिटरी अच्छी तरह से परीक्षण किए गए, मानकीकृत मॉड्यूल के योगदान का स्वागत करती है।

### 3.1. योगदान क्यों करें?

अपने मॉड्यूल को nf-core में योगदान करना:

- आपके टूल को [nf-co.re/modules](https://nf-co.re/modules) पर मॉड्यूल कैटलॉग के माध्यम से पूरे nf-core समुदाय के लिए उपलब्ध कराता है
- निरंतर समुदाय रखरखाव और सुधार सुनिश्चित करता है
- कोड समीक्षा और स्वचालित परीक्षण के माध्यम से गुणवत्ता आश्वासन प्रदान करता है
- आपके काम को दृश्यता और मान्यता देता है

### 3.2. योगदानकर्ता की चेकलिस्ट

nf-core में एक मॉड्यूल योगदान करने के लिए, आपको निम्नलिखित चरणों से गुजरना होगा:

1. जाँचें कि क्या यह पहले से [nf-co.re/modules](https://nf-co.re/modules) पर मौजूद है
2. [nf-core/modules](https://github.com/nf-core/modules) रिपॉजिटरी को fork करें
3. टेम्पलेट जेनरेट करने के लिए `nf-core modules create` का उपयोग करें
4. मॉड्यूल तर्क और परीक्षण भरें
5. `nf-core modules test tool/subtool` के साथ परीक्षण करें
6. `nf-core modules lint tool/subtool` के साथ lint करें
7. एक pull request सबमिट करें

विस्तृत निर्देशों के लिए, [nf-core components tutorial](https://nf-co.re/docs/tutorials/nf-core_components/components) देखें।

### 3.3. संसाधन

- **Components tutorial**: [मॉड्यूल बनाने और योगदान करने के लिए पूर्ण गाइड](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Module specifications**: [तकनीकी आवश्यकताएँ और दिशानिर्देश](https://nf-co.re/docs/guidelines/components/modules)
- **Community support**: [nf-core Slack](https://nf-co.re/join) - `#modules` चैनल में शामिल हों

### सारांश

अब आप जानते हैं कि nf-core मॉड्यूल कैसे बनाएं! आपने चार मुख्य पैटर्न सीखे जो मॉड्यूल को पोर्टेबल और मेंटेनेबल बनाते हैं:

- **मेटाडेटा टपल** workflow के माध्यम से मेटाडेटा प्रसारित करते हैं
- **`ext.args`** कॉन्फ़िगरेशन के माध्यम से वैकल्पिक आर्गुमेंट को संभालकर मॉड्यूल इंटरफ़ेस को सरल बनाता है
- **`ext.prefix`** आउटपुट फ़ाइल नामकरण को मानकीकृत करता है
- **केंद्रीकृत प्रकाशन** मॉड्यूल में हार्डकोड किए जाने के बजाय `modules.config` में कॉन्फ़िगर किए गए `publishDir` के माध्यम से

`COWPY` को चरण-दर-चरण बदलकर, आपने इन पैटर्न की गहरी समझ विकसित की, जिससे आप nf-core मॉड्यूल के साथ काम करने, डिबग करने और बनाने के लिए सुसज्जित हैं।
व्यवहार में, आप शुरुआत से ही इन पैटर्न के साथ उचित रूप से संरचित मॉड्यूल जेनरेट करने के लिए `nf-core modules create` का उपयोग करेंगे।

अंत में, आपने सीखा कि nf-core समुदाय में मॉड्यूल कैसे योगदान करें, टूल को दुनिया भर के शोधकर्ताओं के लिए उपलब्ध कराते हुए निरंतर समुदाय रखरखाव से लाभान्वित होते हुए।

### आगे क्या है?

जब आप तैयार हों, तो [भाग 5: इनपुट सत्यापन](./05_input_validation.md) पर जारी रखें ताकि आप अपनी पाइपलाइन में स्कीमा-आधारित इनपुट सत्यापन जोड़ना सीख सकें।
