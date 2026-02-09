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
    // cowpy के साथ ASCII art जेनरेट करें (https://github.com/jeffbuttars/cowpy)
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
        // cowpy के साथ अभिवादनों का ASCII art generate करें
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // tuple से file extract करें क्योंकि cowpy अभी metadata का उपयोग नहीं करता
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
        // cowpy के साथ अभिवादनों का ASCII art generate करें
        COWPY(CAT_CAT.out.file_out)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // cowpy के साथ अभिवादनों का ASCII art generate करें
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

यदि आप सोच रहे हैं, तो `ext.prefix` closure को सही metadata तक पहुंच है क्योंकि configuration का मूल्यांकन process execution के संदर्भ में किया जाता है, जहां metadata उपलब्ध है।

#### 1.4.3. इसे टेस्ट करने के लिए pipeline चलाएं

आइए जांचें कि workflow अभी भी अपेक्षित रूप से काम करता है।

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

results directory में आउटपुट पर नज़र डालें।
आपको cowpy आउटपुट फ़ाइल पहले जैसी नामकरण के साथ दिखनी चाहिए: `cowpy-test.txt`, डिफ़ॉल्ट बैच नाम पर आधारित।

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

`conf/modules.config` में `ext.prefix` configuration को बदलकर संतुष्ट होने के लिए स्वतंत्र महसूस करें कि आप module या workflow code में कोई बदलाव किए बिना naming pattern बदल सकते हैं।

वैकल्पिक रूप से, आप command line पर एक अलग `--batch` parameter निर्दिष्ट करके इसे फिर से चलाने की कोशिश कर सकते हैं ताकि आप संतुष्ट हो सकें कि वह हिस्सा अभी भी तुरंत अनुकूलन योग्य है।

यह दर्शाता है कि `ext.prefix` आपको module interface को लचीला रखते हुए अपनी पसंदीदा naming convention बनाए रखने की अनुमति देता है।

इस दृष्टिकोण के लाभों का सारांश:

- **मानकीकृत नामकरण**: आउटपुट फ़ाइलें आमतौर पर metadata से sample IDs का उपयोग करके नामित होती हैं
- **कॉन्फ़िगर करने योग्य**: उपयोगकर्ता आवश्यकतानुसार डिफ़ॉल्ट नामकरण को ओवरराइड कर सकते हैं
- **सुसंगत**: सभी nf-core module इस pattern का पालन करते हैं
- **पूर्वानुमान योग्य**: यह जानना आसान है कि आउटपुट फ़ाइलें क्या कहलाएंगी

बहुत अच्छा, है ना?
खैर, nf-core दिशानिर्देशों के अनुसार हमारे module को बेहतर बनाने के लिए एक और महत्वपूर्ण बदलाव करना बाकी है।

### 1.5. Publishing configuration को केंद्रीकृत करें

आपने देखा होगा कि हम दो अलग-अलग directories में outputs publish कर रहे हैं:

- **`results`** — मूल output directory जिसे हम शुरू से अपने local modules के लिए उपयोग कर रहे हैं, प्रति-module `publishDir` directives का उपयोग करके अलग-अलग सेट किया गया;
- **`core-hello-results`** — command line पर `--outdir` के साथ सेट की गई output directory, जो nf-core logs और `CAT_CAT` द्वारा published results प्राप्त कर रही है।

यह गड़बड़ और अपर्याप्त है; सब कुछ के लिए एक स्थान होना बेहतर होगा।
बेशक, हम अपने प्रत्येक local module में जाकर `publishDir` directive को manually अपडेट कर सकते हैं ताकि `core-hello-results` directory का उपयोग किया जा सके, लेकिन अगली बार जब हम output directory बदलने का फैसला करें तो क्या होगा?

व्यक्तिगत modules को publishing निर्णय लेने देना स्पष्ट रूप से सही तरीका नहीं है, विशेष रूप से ऐसी दुनिया में जहां एक ही module कई अलग-अलग pipelines में उपयोग किया जा सकता है, ऐसे लोगों द्वारा जिनकी अलग-अलग ज़रूरतें या प्राथमिकताएं हैं।
हम workflow configuration के स्तर पर यह नियंत्रित करने में सक्षम होना चाहते हैं कि outputs कहाँ publish हों।

"अरे," आप कह सकते हैं, "`CAT_CAT` अपने outputs `--outdir` में भेज रहा है। शायद हमें इसकी `publishDir` directive कॉपी करनी चाहिए?"

हाँ, यह एक बढ़िया विचार है।

सिवाय इसके कि इसमें `publishDir` directive नहीं है। (आगे बढ़ें, module code देखें।)

ऐसा इसलिए है क्योंकि nf-core pipelines workflow स्तर पर नियंत्रण को केंद्रीकृत करती हैं, `conf/modules.config` में `publishDir` कॉन्फ़िगर करके, न कि individual modules में।
विशेष रूप से, nf-core template एक default `publishDir` directive (एक पूर्वनिर्धारित directory structure के साथ) declare करता है जो सभी modules पर लागू होता है जब तक कि एक overriding directive प्रदान नहीं किया जाता।

क्या यह बहुत बढ़िया नहीं लगता? क्या यह हो सकता है कि इस default directive का लाभ उठाने के लिए, हमें बस अपने local modules से वर्तमान `publishDir` directive हटाने की ज़रूरत है?

आइए `COWPY` पर यह आज़माएं और देखें क्या होता है, फिर हम default configuration के code को देखेंगे ताकि समझ सकें कि यह कैसे काम करता है।

अंत में, हम दिखाएंगे कि यदि वांछित हो तो default behavior को कैसे override करें।

#### 1.5.1. `COWPY` से `publishDir` directive हटाएं

चलिए यह करते हैं।
`cowpy.nf` module file (`core-hello/modules/local/` के अंतर्गत) खोलें और नीचे दिखाए अनुसार `publishDir` directive हटाएं।

=== "बाद में"

    ```groovy title="core-hello/modules/local/cowpy.nf (excerpt)" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "पहले"

    ```groovy title="core-hello/modules/local/cowpy.nf (excerpt)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

बस इतना ही!

#### 1.5.2. इसे टेस्ट करने के लिए pipeline चलाएं

आइए देखें कि अब pipeline चलाने पर क्या होता है।

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

अपनी वर्तमान working directory पर नज़र डालें।
अब `core-hello-results` में `COWPY` module के outputs भी शामिल हैं।

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

आप देख सकते हैं कि Nextflow ने workflow और module के नामों के आधार पर directories की यह hierarchy बनाई।

ज़िम्मेदार code `conf/modules.config` फ़ाइल में है।
यह default `publishDir` configuration है जो nf-core template का हिस्सा है और सभी processes पर लागू होता है:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

यह जटिल लग सकता है, तो आइए तीन components में से प्रत्येक को देखें:

- **`path:`** process नाम के आधार पर output directory निर्धारित करता है।
  `task.process` में निहित एक process का पूरा नाम workflow और module imports की hierarchy शामिल करता है (जैसे `CORE_HELLO:HELLO:CAT_CAT`)।
  `tokenize` operations उस hierarchy को हटाकर बस process नाम प्राप्त करते हैं, फिर किसी भी underscore से पहले का पहला भाग लेते हैं (यदि लागू हो), और lowercase में बदलते हैं।
  यही निर्धारित करता है कि `CAT_CAT` के results `${params.outdir}/cat/` में publish होते हैं।
- **`mode:`** यह नियंत्रित करता है कि files कैसे publish होती हैं (copy, symlink, आदि)।
  यह `params.publish_dir_mode` parameter के माध्यम से कॉन्फ़िगर करने योग्य है।
- **`saveAs:`** यह फ़िल्टर करता है कि कौन सी files publish करनी हैं।
  यह उदाहरण `versions.yml` files को उनके लिए `null` return करके बाहर करता है, जिससे उन्हें publish होने से रोका जाता है।

यह outputs को व्यवस्थित करने के लिए एक सुसंगत logic प्रदान करता है।

जब pipeline में सभी module इस convention को अपनाते हैं तो output और भी बेहतर दिखता है, इसलिए अपनी pipeline के अन्य modules से `publishDir` directives हटाने के लिए स्वतंत्र महसूस करें।
यह default उन modules पर भी लागू होगा जिन्हें हमने nf-core दिशानिर्देशों का पालन करने के लिए स्पष्ट रूप से संशोधित नहीं किया।

उस ने कहा, आप तय कर सकते हैं कि आप अपने inputs को अलग तरीके से व्यवस्थित करना चाहते हैं, और अच्छी बात यह है कि ऐसा करना आसान है।

#### 1.5.3. Default को override करें

Default `publishDir` directive को override करने के लिए, आप बस `conf/modules.config` फ़ाइल में अपनी directives जोड़ सकते हैं।

उदाहरण के लिए, आप `withName:` selector का उपयोग करके एक single process के लिए default को override कर सकते हैं, जैसा कि इस उदाहरण में जहां हम 'COWPY' process के लिए एक custom `publishDir` directive जोड़ते हैं।

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

हम वास्तव में यह बदलाव नहीं करने जा रहे हैं, लेकिन इसके साथ खेलने और देखने के लिए स्वतंत्र महसूस करें कि आप कौन सी logic implement कर सकते हैं।

बात यह है कि यह system आपको दोनों दुनिया का सबसे अच्छा देता है: default रूप से consistency और मांग पर configuration को customize करने की flexibility।

सारांश में, आपको मिलता है:

- **Single Source of Truth**: सभी publishing configuration `modules.config` में रहती है
- **उपयोगी default**: Processes प्रति-module configuration के बिना तुरंत काम करते हैं
- **आसान customization**: Publishing behavior को config में override करें, module code में नहीं
- **Portable modules**: Modules में output locations hardcode नहीं होते

यह nf-core module features का वह सेट पूरा करता है जो आपको निश्चित रूप से सीखना चाहिए, लेकिन और भी हैं जिनके बारे में आप [nf-core module specifications](https://nf-co.re/docs/guidelines/components/modules) में पढ़ सकते हैं।

### सारांश

अब आप जानते हैं कि nf-core conventions का पालन करने के लिए local modules को कैसे अनुकूलित करें:

- अपने modules को metadata tuples स्वीकार करने और propagate करने के लिए डिज़ाइन करें;
- Module interfaces को minimal और portable रखने के लिए `ext.args` का उपयोग करें;
- कॉन्फ़िगर करने योग्य, standardized output file naming के लिए `ext.prefix` का उपयोग करें;
- एक consistent results directory structure के लिए default centralized `publishDir` directive अपनाएं।

### आगे क्या है?

जानें कि modules को आसान तरीके से बनाने के लिए nf-core के built-in template-based tools का उपयोग कैसे करें।

---

## 2. nf-core tooling के साथ एक module बनाएं

अब जब आपने nf-core module patterns को manually लागू करके सीख लिया है, तो आइए देखें कि आप practice में modules कैसे बनाएंगे।

### 2.1. एक template से module scaffold generate करें

Pipelines बनाने के लिए जो मौजूद है उसके समान, nf-core project एक template के आधार पर सही ढंग से structured modules generate करने के लिए tooling प्रदान करता है, जिसमें ये सभी patterns शुरू से ही built in होते हैं।

#### 2.1.1. Module creation command चलाएं

`nf-core modules create` command एक module template generate करता है जो पहले से ही आपके द्वारा सीखी गई सभी conventions का पालन करता है।

इस command को चलाकर एक minimal template के साथ `COWPY` module का नया version बनाएं:

```bash
nf-core modules create --empty-template COWPY
```

`--empty-template` flag बिना extra code के एक clean starter template बनाता है, जिससे essential structure को देखना आसान हो जाता है।

Command interactively चलता है, setup के माध्यम से आपका मार्गदर्शन करता है।
यह metadata को pre-populate करने के लिए Bioconda और bio.tools जैसे package repositories से tool information को automatically look up करता है।

आपको कई configuration options के लिए prompt किया जाएगा:

- **Author information**: Attribution के लिए आपका GitHub username
- **Resource label**: Computational requirements का एक predefined set।
  nf-core project lightweight tools के लिए `process_single` और demanding tools के लिए `process_high` जैसे standard labels प्रदान करता है।
  ये labels विभिन्न execution environments में resource allocation प्रबंधित करने में मदद करते हैं।
- **Metadata requirement**: क्या module को `meta` map के माध्यम से sample-specific information की आवश्यकता है (data processing modules के लिए आमतौर पर हाँ)।

Tool package information खोजने और structure सेट करने की complexity को संभालता है, जिससे आप tool की specific logic implement करने पर ध्यान केंद्रित कर सकते हैं।

#### 2.1.2. Module scaffold की जांच करें

Tool `modules/local/` में (या `modules/nf-core/` में यदि आप nf-core/modules repository में हैं) एक complete module structure बनाता है:

??? abstract "डायरेक्टरी सामग्री"

    ```console
    modules/local/cowpy
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
    ```

प्रत्येक file एक विशिष्ट उद्देश्य पूरा करती है:

- **`main.nf`**: सभी nf-core patterns built in के साथ process definition
- **`meta.yml`**: Inputs, outputs और tool का वर्णन करने वाली module documentation
- **`environment.yml`**: Dependencies के लिए Conda environment specification
- **`tests/main.nf.test`**: Module के काम करने को validate करने के लिए nf-test test cases

!!! tip "Testing के बारे में और जानें"

    Generated test file nf-test का उपयोग करती है, जो Nextflow pipelines और modules के लिए एक testing framework है। इन tests को लिखना और चलाना सीखने के लिए, [nf-test side quest](../side_quests/nf-test.md) देखें।

Generated `main.nf` में वे सभी patterns शामिल हैं जो आपने अभी सीखे, साथ ही कुछ अतिरिक्त features:

```groovy title="modules/local/cowpy/main.nf" hl_lines="11 21 22"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(input)        // Pattern 1: Metadata tuples ✓

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                // Pattern 2: ext.args ✓
    def prefix = task.ext.prefix ?: "${meta.id}"  // Pattern 3: ext.prefix ✓

    """
    // Add your tool command here

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

ध्यान दें कि वे सभी patterns जो आपने ऊपर manually लागू किए थे, पहले से ही मौजूद हैं!

Template में कई अतिरिक्त nf-core conventions भी शामिल हैं।
इनमें से कुछ तुरंत काम करते हैं, जबकि अन्य placeholders हैं जिन्हें हमें भरना होगा, जैसा कि नीचे वर्णित है।

**Features जो तुरंत काम करते हैं:**

- **`tag "$meta.id"`**: Logs में process names में sample ID जोड़ता है ताकि tracking आसान हो
- **`label 'process_single'`**: CPU/memory requirements configure करने के लिए resource label
- **`when:` block**: `task.ext.when` configuration के माध्यम से conditional execution की अनुमति देता है

ये features पहले से functional हैं और modules को अधिक maintainable बनाते हैं।

**Placeholders जिन्हें हम नीचे customize करेंगे:**

- **`input:` और `output:` blocks**: Generic declarations जिन्हें हम अपने tool से match करने के लिए update करेंगे
- **`script:` block**: एक comment शामिल है जहां हम `cowpy` command जोड़ेंगे
- **`stub:` block**: Template जिसे हम correct outputs produce करने के लिए update करेंगे
- **Container और environment**: Placeholders जिन्हें हम package information से भरेंगे

अगले sections इन customizations को पूरा करने के बारे में बताते हैं।

### 2.2. Container और Conda environment सेट करें

nf-core guidelines के अनुसार हमें module के हिस्से के रूप में container और Conda environment दोनों specify करने होंगे।

#### 2.2.1. Container

Container के लिए, आप किसी भी Conda package से automatically container बनाने के लिए [Seqera Containers](https://seqera.io/containers/) का उपयोग कर सकते हैं, जिसमें conda-forge packages भी शामिल हैं।
इस मामले में हम पहले जैसा ही prebuilt container उपयोग कर रहे हैं।

Default code Docker और Singularity के बीच toggle करने की पेशकश करता है, लेकिन हम उस line को simplify करेंगे और बस ऊपर Seqera Containers से प्राप्त Docker container specify करेंगे।

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

#### 2.2.2. Conda environment

Conda environment के लिए, module code `conda "${moduleDir}/environment.yml"` specify करता है जिसका अर्थ है कि इसे `environment.yml` file में configure किया जाना चाहिए।

Module creation tool ने हमें चेतावनी दी कि यह Bioconda (bioinformatics tools के लिए primary channel) में `cowpy` package नहीं ढूंढ सका।
हालांकि, `cowpy` conda-forge में उपलब्ध है, इसलिए आप `environment.yml` को इस प्रकार पूरा कर सकते हैं:

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

nf-core में submission के लिए, हमें defaults का अधिक बारीकी से पालन करना होगा, लेकिन अपने स्वयं के उपयोग के लिए हम code को इस तरह simplify कर सकते हैं।

!!! tip "Bioconda vs conda-forge packages"

    - **Bioconda packages**: Automatically BioContainers built मिलते हैं, ready-to-use containers प्रदान करते हैं
    - **conda-forge packages**: Conda recipe से on-demand containers build करने के लिए Seqera Containers का उपयोग कर सकते हैं

    अधिकांश bioinformatics tools Bioconda में हैं, लेकिन conda-forge tools के लिए, Seqera Containers containerization के लिए एक आसान समाधान प्रदान करता है।

### 2.3. `COWPY` logic को जोड़ें

अब उन code elements को update करें जो `COWPY` process क्या करता है उसके लिए specific हैं: inputs और outputs, और script block।

#### 2.3.1. Inputs और outputs

Generated template में generic input और output declarations शामिल हैं जिन्हें आपको अपने specific tool के लिए customize करना होगा।
Section 1 से हमारे manual `COWPY` module को देखते हुए, हम उसे guide के रूप में उपयोग कर सकते हैं।

Input और output blocks को update करें:

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

यह specify करता है:

- Input file parameter name (`input_file` generic `input` के बजाय)
- Configurable prefix pattern का उपयोग करके output filename (`${prefix}.txt` wildcard `*` के बजाय)
- एक descriptive emit name (`cowpy_output` generic `output` के बजाय)

यदि आप syntax validate करने के लिए Nextflow language server का उपयोग कर रहे हैं, तो `${prefix}` part इस stage पर error के रूप में flag किया जाएगा क्योंकि हमने इसे अभी तक script block में नहीं जोड़ा है।
अब उस पर आते हैं।

#### 2.3.2. Script block

Template script block में एक comment placeholder प्रदान करता है जहां आपको actual tool command जोड़ना चाहिए।

पहले manually लिखे गए module के आधार पर, हमें निम्नलिखित edits करने चाहिए:

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
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

मुख्य बदलाव:

- `def prefix` को सिर्फ `prefix` में बदलें (बिना `def`) ताकि इसे output block में accessible बनाया जा सके
- Comment को actual `cowpy` command से replace करें जो `$args` और `${prefix}.txt` दोनों का उपयोग करता है

ध्यान दें कि यदि हमने पहले से `COWPY` process के लिए `ext.args` और `ext.prefix` configuration को `modules.config` file में जोड़ने का काम नहीं किया होता, तो हमें इसे अभी करना होता।

#### 2.3.3. Stub block implement करना

Nextflow context में, एक [stub](https://www.nextflow.io/docs/latest/process.html#stub) block आपको एक lightweight, dummy script define करने की अनुमति देता है जिसका उपयोग actual command execute किए बिना pipeline logic के rapid prototyping और testing के लिए किया जाता है।

<!-- TODO (future) This is super glossed over but should really be explained or at least link out to an explanation about stubs (the reference doc isn't terribly helpful either). Right now this is likely to be mostly meaningless to anyone who doesn't already know about stubs. -->

यदि यह रहस्यमय लगता है तो बहुत चिंता न करें; हम इसे completeness के लिए शामिल करते हैं लेकिन यदि आप इससे निपटना नहीं चाहते तो stub section को delete भी कर सकते हैं, क्योंकि यह पूरी तरह से optional है।

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

मुख्य बदलाव:

- Script block से match करने के लिए `def prefix` को सिर्फ `prefix` में बदलें
- `echo $args` line हटाएं (जो सिर्फ template placeholder code था)
- Stub एक खाली `${prefix}.txt` file बनाता है जो script block के output से match करता है

यह आपको actual tool चलने का इंतज़ार किए बिना workflow logic और file handling test करने की अनुमति देता है।

एक बार जब आपने environment setup (section 2.2), inputs/outputs (section 2.3.1), script block (section 2.3.2), और stub block (section 2.3.3) पूरा कर लिया, तो module test के लिए तैयार है!

### 2.4. नया `COWPY` module लगाएं और pipeline चलाएं

`COWPY` module के इस नए version को आज़माने के लिए हमें बस `hello.nf` workflow file में import statement को नई file की ओर point करने के लिए switch करना होगा।

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

इसे test करने के लिए pipeline चलाएं।

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

यह पहले जैसे ही results produce करता है।

### सारांश

अब आप जानते हैं कि सब कुछ scratch से लिखने के बजाय templates का उपयोग करके efficiently modules बनाने के लिए built-in nf-core tooling का उपयोग कैसे करें।

### आगे क्या है?

जानें कि nf-core में modules contribute करने के क्या लाभ हैं और इसमें शामिल मुख्य कदम और आवश्यकताएं क्या हैं।

---

## 3. nf-core में modules वापस contribute करना

[nf-core/modules](https://github.com/nf-core/modules) repository अच्छी तरह से tested, standardized modules के contributions का स्वागत करता है।

### 3.1. क्यों contribute करें?

nf-core में अपने modules contribute करने से:

- आपके tools [nf-co.re/modules](https://nf-co.re/modules) पर modules catalog के माध्यम से पूरी nf-core community के लिए उपलब्ध होते हैं
- निरंतर community maintenance और improvements सुनिश्चित होती है
- Code review और automated testing के माध्यम से quality assurance मिलती है
- आपके काम को visibility और recognition मिलती है

### 3.2. Contributor's checklist

nf-core में एक module contribute करने के लिए, आपको निम्नलिखित steps से गुज़रना होगा:

1. जांचें कि क्या यह पहले से [nf-co.re/modules](https://nf-co.re/modules) पर मौजूद है
2. [nf-core/modules](https://github.com/nf-core/modules) repository को fork करें
3. Template generate करने के लिए `nf-core modules create` का उपयोग करें
4. Module logic और tests भरें
5. `nf-core modules test tool/subtool` से test करें
6. `nf-core modules lint tool/subtool` से lint करें
7. Pull request submit करें

विस्तृत निर्देशों के लिए, [nf-core components tutorial](https://nf-co.re/docs/tutorials/nf-core_components/components) देखें।

### 3.3. संसाधन

- **Components tutorial**: [Modules बनाने और contribute करने के लिए पूर्ण guide](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Module specifications**: [तकनीकी आवश्यकताएं और guidelines](https://nf-co.re/docs/guidelines/components/modules)
- **Community support**: [nf-core Slack](https://nf-co.re/join) - `#modules` channel join करें

### सारांश

अब आप जानते हैं कि nf-core modules कैसे बनाएं! आपने चार key patterns सीखे जो modules को portable और maintainable बनाते हैं:

- **Metadata tuples** workflow के माध्यम से metadata propagate करते हैं
- **`ext.args`** configuration के माध्यम से optional arguments handle करके module interfaces को simplify करता है
- **`ext.prefix`** output file naming को standardize करता है
- **Centralized publishing** `publishDir` के माध्यम से `modules.config` में configured, modules में hardcode किए जाने के बजाय

`COWPY` को step-by-step transform करके, आपने इन patterns की गहरी समझ विकसित की है, जो आपको nf-core modules के साथ काम करने, debug करने और बनाने के लिए तैयार करती है।
Practice में, आप शुरू से ही इन patterns के साथ correctly structured modules generate करने के लिए `nf-core modules create` का उपयोग करेंगे।

अंत में, आपने सीखा कि nf-core community में modules कैसे contribute करें, दुनिया भर के researchers के लिए tools उपलब्ध कराते हुए निरंतर community maintenance से लाभ उठाएं।

### आगे क्या है?

जब आप तैयार हों, तो अपने pipeline में schema-based input validation जोड़ने के लिए [Part 5: Input validation](./05_input_validation.md) पर जारी रखें।
