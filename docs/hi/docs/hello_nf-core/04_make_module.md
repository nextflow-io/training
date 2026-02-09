# भाग 4: एक nf-core मॉड्यूल बनाओ

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [और जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core प्रशिक्षण पाठ्यक्रम के इस चौथे भाग में, हम तुम्हें दिखाएंगे कि कैसे एक nf-core मॉड्यूल बनाया जाए, उन मुख्य परंपराओं को लागू करके जो मॉड्यूल को पोर्टेबल और मेंटेनेबल बनाती हैं।

nf-core प्रोजेक्ट एक कमांड (`nf-core modules create`) प्रदान करता है जो स्वचालित रूप से सही तरीके से संरचित मॉड्यूल टेम्पलेट जेनरेट करता है, जैसा कि हमने भाग 2 में workflow के लिए उपयोग किया था।
हालांकि, शिक्षण उद्देश्यों के लिए, हम मैन्युअली शुरू करने जा रहे हैं: तुम्हारी `core-hello` पाइपलाइन में लोकल `cowpy` मॉड्यूल को nf-core-स्टाइल मॉड्यूल में चरण-दर-चरण बदलना।
उसके बाद, हम तुम्हें दिखाएंगे कि भविष्य में अधिक कुशलता से काम करने के लिए टेम्पलेट-आधारित मॉड्यूल निर्माण का उपयोग कैसे करें।

??? info "इस सेक्शन से कैसे शुरू करें"

    यह सेक्शन मानता है कि तुमने [भाग 3: एक nf-core मॉड्यूल का उपयोग करो](./03_use_module.md) पूरा कर लिया है और `CAT_CAT` मॉड्यूल को अपनी पाइपलाइन में एकीकृत कर लिया है।

    यदि तुमने भाग 3 पूरा नहीं किया है या इस भाग के लिए नए सिरे से शुरू करना चाहते हो, तो तुम `core-hello-part3` समाधान को अपने शुरुआती बिंदु के रूप में उपयोग कर सकते हो।
    `hello-nf-core/` डायरेक्टरी के अंदर से ये कमांड चलाओ:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    यह तुम्हें एक पाइपलाइन देता है जिसमें `CAT_CAT` मॉड्यूल पहले से एकीकृत है।
    तुम निम्नलिखित कमांड चलाकर परीक्षण कर सकते हो कि यह सफलतापूर्वक चलता है:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. `cowpy` को एक nf-core मॉड्यूल में बदलो

इस सेक्शन में, हम अपनी `core-hello` पाइपलाइन में लोकल `cowpy` मॉड्यूल पर nf-core परंपराओं को लागू करेंगे, इसे एक ऐसे मॉड्यूल में बदलते हुए जो nf-core कम्युनिटी मानकों का पालन करता है।

यह `cowpy` प्रोसेस मॉड्यूल के लिए वर्तमान कोड है:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
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

1. **प्रोसेस नाम को `COWPY` में अपरकेस करो** परंपरा का पालन करने के लिए।
2. **`COWPY` को मेटाडेटा टपल का उपयोग करने के लिए अपडेट करो** ताकि workflow के माध्यम से सैंपल मेटाडेटा को प्रचारित किया जा सके।
3. **`ext.args` के साथ टूल आर्गुमेंट कॉन्फ़िगरेशन को केंद्रीकृत करो** ताकि इंटरफ़ेस को न्यूनतम रखते हुए मॉड्यूल की बहुमुखी प्रतिभा बढ़ाई जा सके।
4. **`ext.prefix` के साथ आउटपुट नामकरण को मानकीकृत करो** ताकि स्थिरता को बढ़ावा मिले।
5. **पब्लिशिंग कॉन्फ़िगरेशन को केंद्रीकृत करो** ताकि स्थिरता को बढ़ावा मिले।

प्रत्येक चरण के बाद, हम पाइपलाइन चलाएंगे यह परीक्षण करने के लिए कि सब कुछ अपेक्षित रूप से काम कर रहा है।

!!! warning "वर्किंग डायरेक्टरी"

    सुनिश्चित करो कि तुम इस सेक्शन में सभी फ़ाइल संपादन और कमांड निष्पादन के लिए `core-hello` डायरेक्टरी (तुम्हारी पाइपलाइन रूट) में हो।

    ```bash
    cd core-hello
    ```

### 1.1. प्रोसेस नाम को अपरकेस करो

यह विशुद्ध रूप से एक स्टाइलिस्टिक परंपरा है (कोई तकनीकी औचित्य नहीं है) लेकिन चूंकि यह nf-core मॉड्यूल के लिए मानक है, तो चलो इसका पालन करते हैं।

हमें तीन सेट परिवर्तन करने होंगे:

1. मॉड्यूल में प्रोसेस नाम अपडेट करो
2. workflow हेडर में मॉड्यूल इंपोर्ट स्टेटमेंट अपडेट करो
3. workflow बॉडी में प्रोसेस कॉल और emit डिक्लेरेशन अपडेट करो

चलो शुरू करते हैं!

#### 1.1.1. मॉड्यूल में प्रोसेस नाम अपडेट करो

`cowpy.nf` मॉड्यूल फ़ाइल (`core-hello/modules/local/` के अंदर) खोलो और प्रोसेस नाम को अपरकेस में संशोधित करो:

=== "बाद में"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "पहले"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

इस मामले में अपरकेसिंग पूरी तरह से सीधी है।

यदि प्रोसेस नाम कई शब्दों से बना होता, उदाहरण के लिए यदि हमारे पास MyCowpyTool नाम की एक प्रोसेस होती जो मूल रूप से camel case में होती, तो nf-core परंपरा उन्हें अलग करने के लिए अंडरस्कोर का उपयोग करना होगा, जिससे MY_COWPY_TOOL मिलेगा।

#### 1.1.2. मॉड्यूल इंपोर्ट स्टेटमेंट अपडेट करो

प्रोसेस नाम केस-सेंसिटिव हैं, इसलिए अब जब हमने प्रोसेस नाम बदल दिया है, तो हमें `hello.nf` के workflow हेडर में मॉड्यूल इंपोर्ट स्टेटमेंट को तदनुसार अपडेट करना होगा:

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

हम इंपोर्ट स्टेटमेंट में एक alias का उपयोग कर सकते हैं ताकि प्रोसेस के कॉल को अपडेट करने से बचा जा सके, लेकिन यह अपरकेसिंग परंपरा को अपनाने के उद्देश्य को कुछ हद तक विफल कर देगा।

#### 1.1.3. प्रोसेस कॉल और emit डिक्लेरेशन अपडेट करो

तो अब चलो `hello.nf` के workflow ब्लॉक में प्रोसेस के दो संदर्भों को अपडेट करते हैं:

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // generate ASCII art of the greetings with cowpy
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
    // generate ASCII art of the greetings with cowpy
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

**दोनों** परिवर्तन करना सुनिश्चित करो, अन्यथा जब तुम इसे चलाओगे तो तुम्हें एक एरर मिलेगी।

#### 1.1.4. इसे परीक्षण करने के लिए पाइपलाइन चलाओ

चलो इन परिवर्तनों के बाद यह परीक्षण करने के लिए workflow चलाते हैं कि सब कुछ सही तरीके से काम कर रहा है।

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

बढ़िया, यह काम करता है! अब चलो अधिक महत्वपूर्ण परिवर्तन करने की ओर बढ़ते हैं।

### 1.2. `COWPY` को मेटाडेटा टपल का उपयोग करने के लिए अपडेट करो

`core-hello` पाइपलाइन के वर्तमान संस्करण में, हम `CAT_CAT` के आउटपुट टपल से फ़ाइल निकाल रहे हैं ताकि `COWPY` को पास किया जा सके, जैसा कि नीचे डायग्राम के ऊपरी आधे हिस्से में दिखाया गया है।

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

यह बेहतर होगा कि `COWPY` सीधे मेटाडेटा टपल स्वीकार करे, जिससे मेटाडेटा workflow के माध्यम से प्रवाहित हो सके, जैसा कि डायग्राम के निचले आधे हिस्से में दिखाया गया है।

इस उद्देश्य के लिए, हमें निम्नलिखित परिवर्तन करने होंगे:

1. इनपुट और आउटपुट डेफिनिशन अपडेट करो
2. workflow में प्रोसेस कॉल अपडेट करो
3. workflow में emit ब्लॉक अपडेट करो

एक बार जब हम यह सब कर लेंगे, तो हम पाइपलाइन चलाएंगे यह परीक्षण करने के लिए कि सब कुछ पहले की तरह काम कर रहा है।

#### 1.2.1. इनपुट और आउटपुट डेफिनिशन अपडेट करो

`cowpy.nf` मॉड्यूल फ़ाइल पर वापस जाओ और इसे नीचे दिखाए अनुसार मेटाडेटा टपल स्वीकार करने के लिए संशोधित करो।

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

जैसा कि तुम देख सकते हो, हमने **मुख्य इनपुट** और **आउटपुट** दोनों को एक टपल में बदल दिया जो `tuple val(meta), path(input_file)` पैटर्न का पालन करता है जो इस प्रशिक्षण के भाग 3 में पेश किया गया था।
आउटपुट के लिए, हमने आउटपुट चैनल को एक वर्णनात्मक नाम देने के लिए `emit: cowpy_output` भी जोड़ा।

अब जब हमने बदल दिया है कि प्रोसेस क्या अपेक्षा करती है, तो हमें प्रोसेस कॉल में जो हम प्रदान करते हैं उसे अपडेट करना होगा।

#### 1.2.2. workflow में प्रोसेस कॉल अपडेट करो

अच्छी खबर यह है कि यह परिवर्तन प्रोसेस कॉल को सरल बना देगा।
अब जब `CAT_CAT` का आउटपुट और `COWPY` का इनपुट एक ही 'आकार' के हैं, यानी दोनों में `tuple val(meta), path(input_file)` संरचना है, तो हम उन्हें सीधे कनेक्ट कर सकते हैं बजाय `CAT_CAT` प्रोसेस के आउटपुट से फ़ाइल को स्पष्ट रूप से निकालने के।

`hello.nf` workflow फ़ाइल (`core-hello/workflows/` के अंदर) खोलो और नीचे दिखाए अनुसार `COWPY` के कॉल को अपडेट करो।

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generate ASCII art of the greetings with cowpy
        COWPY(ch_for_cowpy, params.character)
    ```

अब हम `COWPY` को सीधे `CAT_CAT.out.file_out` पर कॉल करते हैं।

परिणामस्वरूप, हमें अब `ch_for_cowpy` चैनल बनाने की आवश्यकता नहीं है, इसलिए उस लाइन (और उसकी कमेंट लाइन) को पूरी तरह से हटाया जा सकता है।

#### 1.2.3. workflow में emit ब्लॉक अपडेट करो

चूंकि `COWPY` अब एक नामित आउटपुट, `cowpy_output`, emit करता है, हम `hello.nf` workflow के `emit:` ब्लॉक को उसका उपयोग करने के लिए अपडेट कर सकते हैं।

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

यह तकनीकी रूप से आवश्यक नहीं है, लेकिन जब भी संभव हो नामित आउटपुट का संदर्भ देना अच्छा अभ्यास है।

#### 1.2.4. इसे परीक्षण करने के लिए पाइपलाइन चलाओ

चलो इन परिवर्तनों के बाद यह परीक्षण करने के लिए workflow चलाते हैं कि सब कुछ सही तरीके से काम कर रहा है।

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

यह पूरा करता है कि हमें `COWPY` को मेटाडेटा टपल संभालने के लिए क्या करना था।
अब, चलो देखते हैं कि हम nf-core मॉड्यूल पैटर्न का लाभ उठाने के लिए और क्या कर सकते हैं।

### 1.3. `ext.args` के साथ टूल आर्गुमेंट कॉन्फ़िगरेशन को केंद्रीकृत करो

अपनी वर्तमान स्थिति में, `COWPY` प्रोसेस `character` पैरामीटर के लिए एक मान प्राप्त करने की अपेक्षा करती है।
परिणामस्वरूप, हमें हर बार प्रोसेस को कॉल करते समय एक मान प्रदान करना होता है, भले ही हम टूल द्वारा सेट किए गए डिफ़ॉल्ट से खुश हों।
`COWPY` के लिए यह स्वीकार्य रूप से बड़ी समस्या नहीं है, लेकिन कई वैकल्पिक पैरामीटर वाले टूल के लिए, यह काफी बोझिल हो सकता है।

nf-core प्रोजेक्ट [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) नामक एक Nextflow फीचर का उपयोग करने की सिफारिश करता है ताकि कॉन्फ़िगरेशन फ़ाइलों के माध्यम से टूल आर्गुमेंट को अधिक सुविधाजनक तरीके से प्रबंधित किया जा सके।

हर टूल ऑप्शन के लिए प्रोसेस इनपुट घोषित करने के बजाय, तुम मॉड्यूल को अपनी कमांड लाइन के निर्माण में `ext.args` का संदर्भ देने के लिए लिखते हो।
फिर यह सिर्फ `modules.config` फ़ाइल में `ext.args` वेरिएबल को सेट करने की बात है, जो सभी मॉड्यूल के लिए कॉन्फ़िगरेशन विवरण को समेकित करती है, ताकि तुम जो आर्गुमेंट और मान उपयोग करना चाहते हो उन्हें रखा जा सके।
Nextflow रनटाइम पर उन आर्गुमेंट को उनके मानों के साथ टूल कमांड लाइन में जोड़ देगा।

चलो इस दृष्टिकोण को `COWPY` मॉड्यूल पर लागू करते हैं।
हमें निम्नलिखित परिवर्तन करने होंगे:

1. `COWPY` मॉड्यूल अपडेट करो
2. `modules.config` फ़ाइल में `ext.args` कॉन्फ़िगर करो
3. `hello.nf` workflow अपडेट करो

एक बार जब हम यह सब कर लेंगे, तो हम पाइपलाइन चलाएंगे यह परीक्षण करने के लिए कि सब कुछ पहले की तरह काम कर रहा है।

#### 1.3.1. `COWPY` मॉड्यूल अपडेट करो

चलो करते हैं।
`cowpy.nf` मॉड्यूल फ़ाइल (`core-hello/modules/local/` के अंदर) खोलो और इसे नीचे दिखाए अनुसार `ext.args` का संदर्भ देने के लिए संशोधित करो।

=== "बाद में"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
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

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
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

तुम देख सकते हो कि हमने तीन परिवर्तन किए।

1. **`input:` ब्लॉक में, हमने `val character` इनपुट हटा दिया।**
   आगे बढ़ते हुए, हम उस आर्गुमेंट को `ext.args` कॉन्फ़िगरेशन के माध्यम से प्रदान करेंगे जैसा कि नीचे वर्णित है।

2. **`script:` ब्लॉक में, हमने लाइन `def args = task.ext.args ?: ''` जोड़ी।**
   वह लाइन `?:` ऑपरेटर का उपयोग करती है `args` वेरिएबल के मान को निर्धारित करने के लिए: `task.ext.args` की सामग्री यदि यह खाली नहीं है, या एक खाली स्ट्रिंग यदि यह है।
   ध्यान दें कि जबकि हम आम तौर पर `ext.args` का संदर्भ देते हैं, इस कोड को मॉड्यूल-स्तरीय `ext.args` कॉन्फ़िगरेशन को निकालने के लिए `task.ext.args` का संदर्भ देना चाहिए।

3. **कमांड लाइन में, हमने `-c "$character"` को `$args` से बदल दिया।**
   यह वह जगह है जहां Nextflow `modules.config` फ़ाइल में `ext.args` में सेट किए गए किसी भी टूल आर्गुमेंट को इंजेक्ट करेगा।

परिणामस्वरूप, मॉड्यूल इंटरफ़ेस अब सरल है: यह केवल आवश्यक मेटाडेटा और फ़ाइल इनपुट की अपेक्षा करता है।

!!! note

    `?:` ऑपरेटर को अक्सर 'Elvis ऑपरेटर' कहा जाता है क्योंकि यह एक साइडवेज Elvis Presley चेहरे की तरह दिखता है, `?` कैरेक्टर उनके बालों में लहर का प्रतीक है।

#### 1.3.2. `modules.config` फ़ाइल में `ext.args` कॉन्फ़िगर करो

अब जब हमने मॉड्यूल से `character` डिक्लेरेशन निकाल दिया है, तो हमें इसे `modules.config` कॉन्फ़िगरेशन फ़ाइल में `ext.args` में जोड़ना होगा।

विशेष रूप से, हम `process {}` ब्लॉक में कोड का यह छोटा टुकड़ा जोड़ने जा रहे हैं:

```groovy title="जोड़ने के लिए कोड"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

`withName:` सिंटैक्स इस कॉन्फ़िगरेशन को केवल `COWPY` प्रोसेस को असाइन करता है, और `ext.args = { "-c ${params.character}" }` बस एक स्ट्रिंग बनाता है जिसमें `character` पैरामीटर का मान शामिल होगा।
कर्ली ब्रेसेस के उपयोग पर ध्यान दें, जो Nextflow को रनटाइम पर पैरामीटर के मान का मूल्यांकन करने के लिए कहते हैं।

समझ में आया? चलो इसे जोड़ते हैं।

`conf/modules.config` खोलो और नीचे दिखाए अनुसार `process {}` ब्लॉक के अंदर कॉन्फ़िगरेशन कोड जोड़ो।

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

उम्मीद है कि तुम कल्पना कर सकते हो कि एक पाइपलाइन में सभी मॉड्यूल के लिए उनके `ext.args` इस फ़ाइल में निर्दिष्ट हों, निम्नलिखित लाभों के साथ:

- **मॉड्यूल इंटरफ़ेस सरल रहता है** - यह केवल आवश्यक मेटाडेटा और फ़ाइल इनपुट स्वीकार करता है
- **पाइपलाइन अभी भी `params.character` को एक्सपोज़ करती है** - अंतिम उपयोगकर्ता अभी भी इसे पहले की तरह कॉन्फ़िगर कर सकते हैं
- **मॉड्यूल अब पोर्टेबल है** - इसे अन्य पाइपलाइन में किसी विशिष्ट पैरामीटर नाम की अपेक्षा के बिना पुन: उपयोग किया जा सकता है
- कॉन्फ़िगरेशन `modules.config` में **केंद्रीकृत** है, workflow लॉजिक को साफ रखते हुए

`modules.config` फ़ाइल को उस स्थान के रूप में उपयोग करके जहां सभी पाइपलाइन प्रति-मॉड्यूल कॉन्फ़िगरेशन को केंद्रीकृत करती हैं, हम अपने मॉड्यूल को विभिन्न पाइपलाइन में अधिक पुन: प्रयोज्य बनाते हैं।

#### 1.3.3. `hello.nf` workflow अपडेट करो

चूंकि `COWPY` मॉड्यूल को अब इनपुट के रूप में `character` पैरामीटर की आवश्यकता नहीं है, हमें तदनुसार workflow कॉल को अपडेट करना होगा।

`hello.nf` workflow फ़ाइल (`core-hello/workflows/` के अंदर) खोलो और नीचे दिखाए अनुसार `COWPY` के कॉल को अपडेट करो।

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

workflow कोड अब साफ है: हमें प्रोसेस को सीधे `params.character` पास करने की आवश्यकता नहीं है।
मॉड्यूल इंटरफ़ेस को न्यूनतम रखा गया है, जिससे यह अधिक पोर्टेबल हो जाता है, जबकि पाइपलाइन अभी भी कॉन्फ़िगरेशन के माध्यम से स्पष्ट विकल्प प्रदान करती है।

#### 1.3.4. इसे परीक्षण करने के लिए पाइपलाइन चलाओ

चलो परीक्षण करते हैं कि workflow अभी भी अपेक्षित रूप से काम करता है, एक अलग कैरेक्टर निर्दिष्ट करके यह सत्यापित करने के लिए कि `ext.args` कॉन्फ़िगरेशन काम कर रहा है।

`kosh` का उपयोग करके यह कमांड चलाओ, जो अधिक... रहस्यमय विकल्पों में से एक है:

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

चलो सत्यापित करते हैं कि `ext.args` कॉन्फ़िगरेशन ने काम किया आउटपुट की जांच करके।
फ़ाइल ब्राउज़र में आउटपुट खोजो या आउटपुट फ़ाइल को देखने के लिए टास्क हैश (ऊपर के उदाहरण में `38/eb29ea` भाग) का उपयोग करो:

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

तुम्हें `kosh` कैरेक्टर के साथ प्रदर्शित ASCII आर्ट दिखाई देनी चाहिए, जो पुष्टि करती है कि `ext.args` कॉन्फ़िगरेशन ने काम किया!

??? info "(वैकल्पिक) कमांड फ़ाइल का निरीक्षण करो"

    यदि तुम देखना चाहते हो कि कॉन्फ़िगरेशन कैसे लागू किया गया था, तो तुम `.command.sh` फ़ाइल का निरीक्षण कर सकते हो:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    तुम्हें `-c kosh` आर्गुमेंट के साथ `cowpy` कमांड दिखाई देगी:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    यह दिखाता है कि `.command.sh` फ़ाइल `ext.args` कॉन्फ़िगरेशन के आधार पर सही तरीके से जेनरेट की गई थी।

एक पल के लिए सोचो कि हमने यहां क्या हासिल किया।
यह दृष्टिकोण मॉड्यूल इंटरफ़ेस को आवश्यक डेटा (फ़ाइलें, मेटाडेटा, और कोई भी अनिवार्य प्रति-सैंपल पैरामीटर) पर केंद्रित रखता है, जबकि टूल के व्यवहार को नियंत्रित करने वाले विकल्पों को कॉन्फ़िगरेशन के माध्यम से अलग से संभाला जाता है।

यह `cowpy` जैसे सरल टूल के लिए अनावश्यक लग सकता है, लेकिन यह डेटा विश्लेषण टूल के लिए बड़ा अंतर ला सकता है जिनमें बहुत सारे वैकल्पिक आर्गुमेंट होते हैं।

इस दृष्टिकोण के लाभों को संक्षेप में प्रस्तुत करने के लिए:

- **साफ इंटरफ़ेस**: मॉड्यूल आवश्यक डेटा इनपुट (मेटाडेटा और फ़ाइलें) पर केंद्रित है
- **लचीलापन**: उपयोगकर्ता कॉन्फ़िगरेशन के माध्यम से टूल आर्गुमेंट निर्दिष्ट कर सकते हैं, जिसमें सैंपल-विशिष्ट मान शामिल हैं
- **स्थिरता**: सभी nf-core मॉड्यूल इस पैटर्न का पालन करते हैं
- **पोर्टेबिलिटी**: मॉड्यूल को हार्डकोडेड टूल विकल्पों के बिना पुन: उपयोग किया जा सकता है
- **कोई workflow परिवर्तन नहीं**: टूल विकल्पों को जोड़ने या बदलने के लिए workflow कोड को अपडेट करने की आवश्यकता नहीं है

!!! note

    `ext.args` सिस्टम में शक्तिशाली अतिरिक्त क्षमताएं हैं जो यहां कवर नहीं की गई हैं, जिसमें मेटाडेटा के आधार पर आर्गुमेंट मानों को गतिशील रूप से स्विच करना शामिल है। अधिक विवरण के लिए [nf-core मॉड्यूल विनिर्देश](https://nf-co.re/docs/guidelines/components/modules) देखें।

### 1.4. `ext.prefix` के साथ आउटपुट नामकरण को मानकीकृत करो

अब जब हमने `COWPY` प्रोसेस को metamap तक पहुंच दी है, तो हम एक और उपयोगी nf-core पैटर्न का लाभ उठाना शुरू कर सकते हैं: मेटाडेटा के आधार पर आउटपुट फ़ाइलों का नामकरण।

यहां हम `ext.prefix` नामक एक Nextflow फीचर का उपयोग करने जा रहे हैं जो हमें `meta.id` (metamap में शामिल आइडेंटिफ़ायर) का उपयोग करके मॉड्यूल में आउटपुट फ़ाइल नामकरण को मानकीकृत करने की अनुमति देगा, जबकि यदि वांछित हो तो अभी भी मॉड्यूल को व्यक्तिगत रूप से कॉन्फ़िगर करने में सक्षम होंगे।

यह `ext.args` के साथ हमने जो किया उसके समान होगा, कुछ अंतरों के साथ जिन्हें हम आगे बढ़ते हुए विस्तार से बताएंगे।

चलो इस दृष्टिकोण को `COWPY` मॉड्यूल पर लागू करते हैं।
हमें निम्नलिखित परिवर्तन करने होंगे:

1. `COWPY` मॉड्यूल अपडेट करो
2. `modules.config` फ़ाइल में `ext.prefix` कॉन्फ़िगर करो

(workflow में कोई परिवर्तन की आवश्यकता नहीं है।)

एक बार जब हम यह कर लेंगे, तो हम पाइपलाइन चलाएंगे यह परीक्षण करने के लिए कि सब कुछ पहले की तरह काम कर रहा है।

#### 1.4.1. `COWPY` मॉड्यूल अपडेट करो

`cowpy.nf` मॉड्यूल फ़ाइल (`core-hello/modules/local/` के अंदर) खोलो और इसे नीचे दिखाए अनुसार `ext.prefix` का संदर्भ देने के लिए संशोधित करो।

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

तुम देख सकते हो कि हमने तीन परिवर्तन किए।

1. **`script:` ब्लॉक में, हमने लाइन `prefix = task.ext.prefix ?: "${meta.id}"` जोड़ी।**
   वह लाइन `?:` ऑपरेटर का उपयोग करती है `prefix` वेरिएबल के मान को निर्धारित करने के लिए: `task.ext.prefix` की सामग्री यदि यह खाली नहीं है, या metamap से आइडेंटिफ़ायर (`meta.id`) यदि यह है।
   ध्यान दें कि जबकि हम आम तौर पर `ext.prefix` का संदर्भ देते हैं, इस कोड को मॉड्यूल-स्तरीय `ext.prefix` कॉन्फ़िगरेशन को निकालने के लिए `task.ext.prefix` का संदर्भ देना चाहिए।

2. **कमांड लाइन में, हमने `cowpy-${input_file}` को `${prefix}.txt` से बदल दिया।**
   यह वह जगह है जहां Nextflow ऊपर की लाइन द्वारा निर्धारित `prefix` के मान को इंजेक्ट करेगा।

3. **`output:` ब्लॉक में, हमने `path("cowpy-${input_file}")` को `path("${prefix}.txt")` से बदल दिया।**
   यह बस दोहराता है कि कमांड लाइन में लिखे अनुसार फ़ाइल पथ क्या होगा।

परिणामस्वरूप, आउटपुट फ़ाइल नाम अब एक समझदार डिफ़ॉल्ट (metamap से आइडेंटिफ़ायर) का उपयोग करके उपयुक्त फ़ाइल फ़ॉर्मेट एक्सटेंशन के साथ बनाया गया है।

#### 1.4.2. `modules.config` फ़ाइल में `ext.prefix` कॉन्फ़िगर करो

इस मामले में समझदार डिफ़ॉल्ट हमारे स्वाद के लिए पर्याप्त अभिव्यंजक नहीं है; हम एक कस्टम नामकरण पैटर्न का उपयोग करना चाहते हैं जिसमें टूल नाम शामिल हो, `cowpy-<id>.txt`, जैसा कि हमारे पास पहले था।

हम `modules.config` में `ext.prefix` को कॉन्फ़िगर करके ऐसा करेंगे, जैसा कि हमने `ext.args` के साथ `character` पैरामीटर के लिए किया था, सिवाय इस बार `withName: 'COWPY' {}` ब्लॉक पहले से मौजूद है, और हमें बस निम्नलिखित लाइन जोड़नी है:

```groovy title="जोड़ने के लिए कोड"
ext.prefix = { "cowpy-${meta.id}" }
```

यह वह स्ट्रिंग बनाएगा जो हम चाहते हैं।
ध्यान दें कि एक बार फिर हम कर्ली ब्रेसेस का उपयोग करते हैं, इस बार Nextflow को रनटाइम पर `meta.id` के मान का मूल्यांकन करने के लिए कहने के लिए।

चलो इसे जोड़ते हैं।

`conf/modules.config` खोलो और नीचे दिखाए अनुसार `process {}` ब्लॉक के अंदर कॉन्फ़िगरेशन कोड जोड़ो।

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

यदि तुम सोच रहे हो, तो `ext.prefix` closure को मेटाडेटा के सही टुकड़े तक पहुंच है क्योंकि कॉन्फ़िगरेशन का मूल्यांकन प्रोसेस निष्पादन के संदर्भ में किया जाता है, जहां मेटाडेटा उपलब्ध है।

#### 1.4.3. इसे परीक्षण करने के लिए पाइपलाइन चलाओ

चलो परीक्षण करते हैं कि workflow अभी भी अपेक्षित रूप से काम करता है।

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

results डायरेक्टरी में आउटपुट पर एक नज़र डालो।
तुम्हें पहले की तरह ही नामकरण के साथ cowpy आउटपुट फ़ाइल दिखाई देनी चाहिए: `cowpy-test.txt`, डिफ़ॉल्ट बैच नाम के आधार पर।

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

अपने आप को संतुष्ट करने के लिए `conf/modules.config` में `ext.prefix` कॉन्फ़िगरेशन को बदलने के लिए स्वतंत्र महसूस करो कि तुम मॉड्यूल या workflow कोड में कोई परिवर्तन किए बिना नामकरण पैटर्न को बदल सकते हो।

वैकल्पिक रूप से, तुम कमांड लाइन पर निर्दिष्ट एक अलग `--batch` पैरामीटर के साथ इसे फिर से चलाने की कोशिश कर सकते हो ताकि अपने आप को संतुष्ट कर सको कि वह हिस्सा अभी भी तुरंत अनुकूलन योग्य है।

यह दर्शाता है कि कैसे `ext.prefix` तुम्हें अपनी पसंदीदा नामकरण परंपरा को बनाए रखने की अनुमति देता है जबकि मॉड्यूल इंटरफ़ेस को लचीला रखता है।

इस दृष्टिकोण के लाभों को संक्षेप में प्रस्तुत करने के लिए:

- **मानकीकृत नामकरण**: आउटपुट फ़ाइलें आमतौर पर मेटाडेटा से सैंपल ID का उपयोग करके नामित की जाती हैं
- **कॉन्फ़िगर करने योग्य**: यदि आवश्यक हो तो उपयोगकर्ता डिफ़ॉल्ट नामकरण को ओवरराइड कर सकते हैं
- **स्थिर**: सभी nf-core मॉड्यूल इस पैटर्न का पालन करते हैं
- **अनुमानित**: यह जानना आसान है कि आउटपुट फ़ाइलें क्या कहलाएंगी

काफी अच्छा, है ना?
खैर, हमारे मॉड्यूल को nf-core दिशानिर्देशों के अनुरूप बनाने के लिए एक और महत्वपूर्ण परिवर्तन करना है।

### 1.5. पब्लिशिंग कॉन्फ़िगरेशन को केंद्रीकृत करो

तुमने देखा होगा कि हम दो अलग-अलग डायरेक्टरी में आउटपुट प्रकाशित कर रहे हैं:

- **`results`** — मूल आउटपुट डायरेक्टरी जिसे हम शुरुआत से अपने लोकल मॉड्यूल के लिए उपयोग कर रहे हैं, प्रति-मॉड्यूल `publishDir` निर्देशों का उपयोग करके व्यक्तिगत रूप से सेट की गई;
- **`core-hello-results`** — कमांड लाइन पर `--outdir` के साथ सेट की गई आउटपुट डायरेक्टरी, जो nf-core लॉग और `CAT_CAT` द्वारा प्रकाशित परिणाम प्राप्त कर रही है।

यह गड़बड़ और उप-इष्टतम है; सब कुछ के लिए एक स्थान होना बेहतर होगा।
बेशक, हम अपने प्रत्येक लोकल मॉड्यूल में जा सकते हैं और `core-hello-results` डायरेक्टरी का उपयोग करने के लिए `publishDir` निर्देश को मैन्युअली अपडेट कर सकते हैं, लेकिन अगली बार जब हम आउटपुट डायरेक्टरी बदलने का निर्णय लेते हैं तो क्या होगा?

व्यक्तिगत मॉड्यूल को प्रकाशन निर्णय लेने देना स्पष्ट रूप से जाने का तरीका नहीं है, विशेष रूप से एक ऐसी दुनिया में जहां एक ही मॉड्यूल का उपयोग बहुत सारी अलग-अलग पाइपलाइन में किया जा सकता है, उन लोगों द्वारा जिनकी अलग-अलग आवश्यकताएं या प्राथमिकताएं हैं।
हम workflow कॉन्फ़िगरेशन के स्तर पर यह नियंत्रित करने में सक्षम होना चाहते हैं कि आउटपुट कहां प्रकाशित होते हैं।

"अरे," तुम कह सकते हो, "`CAT_CAT` अपने आउटपुट को `--outdir` में भेज रहा है। शायद हमें इसके `publishDir` निर्देश की नकल करनी चाहिए?"

हां, यह एक बढ़िया विचार है।

सिवाय इसके कि इसमें `publishDir` निर्देश नहीं है। (आगे बढ़ो, मॉड्यूल कोड देखो।)

ऐसा इसलिए है क्योंकि nf-core पाइपलाइन व्यक्तिगत मॉड्यूल में `publishDir` को कॉन्फ़िगर करने के बजाय `conf/modules.config` में workflow स्तर पर नियंत्रण को केंद्रीकृत करती हैं।
विशेष रूप से, nf-core टेम्पलेट एक डिफ़ॉल्ट `publishDir` निर्देश (एक पूर्वनिर्धारित डायरेक्टरी संरचना के साथ) घोषित करता है जो सभी मॉड्यूल पर लागू होता है जब तक कि एक ओवरराइडिंग निर्देश प्रदान नहीं किया जाता है।

क्या यह अद्भुत नहीं लगता? क्या यह हो सकता है कि इस डिफ़ॉल्ट निर्देश का लाभ उठाने के लिए, हमें बस अपने लोकल मॉड्यूल से वर्तमान `publishDir` निर्देश को हटाना है?

चलो यह देखने के लिए `COWPY` पर कोशिश करते हैं कि क्या होता है, फिर हम यह समझने के लिए डिफ़ॉल्ट कॉन्फ़िगरेशन के लिए कोड देखेंगे कि यह कैसे काम करता है।

अंत में, हम प्रदर्शित करेंगे कि यदि वांछित हो तो डिफ़ॉल्ट व्यवहार को कैसे ओवरराइड किया जाए।

#### 1.5.1. `COWPY` से `publishDir` निर्देश हटाओ

चलो यह करते हैं।
`cowpy.nf` मॉड्यूल फ़ाइल (`core-hello/modules/local/` के अंदर) खोलो और नीचे दिखाए अनुसार `publishDir` निर्देश हटाओ।

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

#### 1.5.2. इसे परीक्षण करने के लिए पाइपलाइन चलाओ

चलो देखते हैं कि अगर हम अब पाइपलाइन चलाते हैं तो क्या होता है।

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

अपनी वर्तमान वर्किंग डायरेक्टरी पर एक नज़र डालो।
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

तुम देख सकते हो कि Nextflow ने workflow और मॉड्यूल के नामों के आधार पर डायरेक्टरी की यह पदानुक्रम बनाई।

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

यह जटिल लग सकता है, तो चलो तीन घटकों में से प्रत्येक को देखते हैं:

- **`path:`** प्रोसेस नाम के आधार पर आउटपुट डायरेक्टरी निर्धारित करता है।
  `task.process` में निहित एक प्रोसेस का पूरा नाम workflow और मॉड्यूल इंपोर्ट की पदानुक्रम शामिल करता है (जैसे `CORE_HELLO:HELLO:CAT_CAT`)।
  `tokenize` ऑपरेशन उस पदानुक्रम को हटाकर केवल प्रोसेस नाम प्राप्त करते हैं, फिर किसी भी अंडरस्कोर से पहले पहला भाग लेते हैं (यदि लागू हो), और इसे लोअरकेस में बदल देते हैं।
  यह वह है जो निर्धारित करता है कि `CAT_CAT` के परिणाम `${params.outdir}/cat/` में प्रकाशित होते हैं।
- **`mode:`** नियंत्रित करता है कि फ़ाइलें कैसे प्रकाशित की जाती हैं (copy, symlink, आदि)।
  यह `params.publish_dir_mode` पैरामीटर के माध्यम से कॉन्फ़िगर करने योग्य है।
- **`saveAs:`** फ़िल्टर करता है कि कौन सी फ़ाइलें प्रकाशित की जाएं।
  यह उदाहरण उनके लिए `null` लौटाकर `versions.yml` फ़ाइलों को बाहर करता है, उन्हें प्रकाशित होने से रोकता है।

यह आउटपुट को व्यवस्थित करने के लिए एक सुसंगत तर्क प्रदान करता है।

जब एक पाइपलाइन में सभी मॉड्यूल इस परंपरा को अपनाते हैं तो आउटपुट और भी बेहतर दिखता है, इसलिए अपनी पाइपलाइन में अन्य मॉड्यूल से `publishDir` निर्देशों को हटाने के लिए स्वतंत्र महसूस करो।
यह डिफ़ॉल्ट उन मॉड्यूल पर भी लागू होगा जिन्हें हमने nf-core दिशानिर्देशों का पालन करने के लिए स्पष्ट रूप से संशोधित नहीं किया।

उस ने कहा, तुम निर्णय ले सकते हो कि तुम अपने इनपुट को अलग तरीके से व्यवस्थित करना चाहते हो, और अच्छी खबर यह है कि ऐसा करना आसान है।

#### 1.5.3. डिफ़ॉल्ट को ओवरराइड करो

डिफ़ॉल्ट `publishDir` निर्देश को ओवरराइड करने के लिए, तुम बस `conf/modules.config` फ़ाइल में अपने स्वयं के निर्देश जोड़ सकते हो।

उदाहरण के लिए, तुम `withName:` सिलेक्टर का उपयोग करके एक एकल प्रोसेस के लिए डिफ़ॉल्ट को ओवरराइड कर सकते हो, जैसा कि इस उदाहरण में है जहां हम 'COWPY' प्रोसेस के लिए एक कस्टम `publishDir` निर्देश जोड़ते हैं।

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

हम वास्तव में वह परिवर्तन नहीं करने जा रहे हैं, लेकिन इसके साथ खेलने और देखने के लिए स्वतंत्र महसूस करो कि तुम कौन सा तर्क लागू कर सकते हो।

बात यह है कि यह सिस्टम तुम्हें दोनों दुनिया का सर्वश्रेष्ठ देता है: डिफ़ॉल्ट रूप से स्थिरता और मांग पर कॉन्फ़िगरेशन को अनुकूलित करने की लचीलापन।

संक्षेप में, तुम्हें मिलता है:

- **सत्य का एकल स्रोत**: सभी प्रकाशन कॉन्फ़िगरेशन `modules.config` में रहता है
- **उपयोगी डिफ़ॉल्ट**: प्रोसेस प्रति-मॉड्यूल कॉन्फ़िगरेशन के बिना out-of-the-box काम करती हैं
- **आसान अनुकूलन**: मॉड्यूल कोड में नहीं, config में प्रकाशन व्यवहार को ओवरराइड करो
- **पोर्टेबल मॉड्यूल**: मॉड्यूल आउटपुट स्थानों को हार्डकोड नहीं करते हैं

यह nf-core मॉड्यूल फीचर्स के सेट को पूरा करता है जिन्हें तुम्हें बिल्कुल उपयोग करना सीखना चाहिए, लेकिन अन्य भी हैं जिनके बारे में तुम [nf-core मॉड्यूल विनिर्देश](https://nf-co.re/docs/guidelines/components/modules) में पढ़ सकते हो।

### सारांश

अब तुम जानते हो कि nf-core परंपराओं का पालन करने के लिए लोकल मॉड्यूल को कैसे अनुकूलित किया जाए:

- अपने मॉड्यूल को मेटाडेटा टपल स्वीकार करने और प्रचारित करने के लिए डिज़ाइन करो;
- मॉड्यूल इंटरफ़ेस को न्यूनतम और पोर्टेबल रखने के लिए `ext.args` का उपयोग करो;
- कॉन्फ़िगर करने योग्य, मानकीकृत आउटपुट फ़ाइल नामकरण के लिए `ext.prefix` का उपयोग करो;
- एक सुसंगत results डायरेक्टरी संरचना के लिए डिफ़ॉल्ट केंद्रीकृत `publishDir` निर्देश को अपनाओ।

### आगे क्या है?

जानो कि आसान तरीके से मॉड्यूल बनाने के लिए nf-core के बिल्ट-इन टेम्पलेट-आधारित टूल का उपयोग कैसे करें।

---

## 2. nf-core टूलिंग के साथ एक मॉड्यूल बनाओ

अब जब तुमने मैन्युअली उन्हें लागू करके nf-core मॉड्यूल पैटर्न सीख लिए हैं, तो चलो देखते हैं कि तुम व्यवहार में मॉड्यूल कैसे बनाओगे।

### 2.1. एक टेम्पलेट से एक मॉड्यूल scaffold जेनरेट करो

पाइपलाइन बनाने के लिए जो मौजूद है उसके समान, nf-core प्रोजेक्ट एक टेम्पलेट के आधार पर सही तरीके से संरचित मॉड्यूल जेनरेट करने के लिए टूलिंग प्रदान करता है, जिसमें ये सभी पैटर्न शुरुआत से ही बनाए गए हैं।

#### 2.1.1. मॉड्यूल निर्माण कमांड चलाओ

`nf-core modules create` कमांड एक मॉड्यूल टेम्पलेट जेनरेट करता है जो पहले से ही उन सभी परंपराओं का पालन करता है जो तुमने सीखी हैं।

चलो एक न्यूनतम टेम्पलेट के साथ `COWPY` मॉड्यूल का एक नया संस्करण बनाते हैं यह कमांड चलाकर:

```bash
nf-core modules create --empty-template COWPY
```

`--empty-template` फ्लैग अतिरिक्त कोड के बिना एक साफ स्टार्टर टेम्पलेट बनाता है, जिससे आवश्यक संरचना को देखना आसान हो जाता है।

कमांड इंटरैक्टिव रूप से चलता है, तुम्हें सेटअप के माध्यम से मार्गदर्शन करता है।
यह स्वचालित रूप से Bioconda और bio.tools जैसे पैकेज रिपॉजिटरी से टूल जानकारी खोजता है ताकि मेटाडेटा को पूर्व-भरा जा सके।

तुम्हें कई कॉन्फ़िगरेशन विकल्पों के लिए प्रॉम्प्ट किया जाएगा:

- **लेखक जानकारी**: एट्रिब्यूशन के लिए तुम्हारा GitHub यूज़रनेम
- **रिसोर्स लेबल**: कम्प्यूटेशनल आवश्यकताओं का एक पूर्वनिर्धारित सेट।
  nf-core प्रोजेक्ट हल्के टूल के लिए `process_single` और मांग वाले टूल के लिए `process_high` जैसे मानक लेबल प्रदान करता है।
  ये लेबल विभिन्न निष्पादन वातावरण में रिसोर्स आवंटन को प्रबंधित करने में मदद करते हैं।
- **मेटाडेटा आवश्यकता**: क्या मॉड्यूल को `meta` map के माध्यम से सैंपल-विशिष्ट जानकारी की आवश्यकता है (आमतौर पर डेटा प्रोसेसिंग मॉड्यूल के लिए हां)।

टूल पैकेज जानकारी खोजने और संरचना सेट करने की जटिलता को संभालता है, जिससे तुम टूल के विशिष्ट तर्क को लागू करने पर ध्यान केंद्रित कर सकते हो।

#### 2.1.2. मॉड्यूल scaffold की जांच करो

टूल `modules/local/` में एक पूर्ण मॉड्यूल संरचना बनाता है (या `modules/nf-core/` यदि तुम nf-core/modules रिपॉजिटरी में हो):

??? abstract "डायरेक्टरी सामग्री"

    ```console
    modules/local/cowpy
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
    ```

प्रत्येक फ़ाइल एक विशिष्ट उद्देश्य की पूर्ति करती है:

- **`main.nf`**: सभी nf-core पैटर्न के साथ प्रोसेस डेफिनिशन बिल्ट-इन
- **`meta.yml`**: इनपुट, आउटपुट और टूल का वर्णन करने वाला मॉड्यूल डॉक्यूमेंटेशन
- **`environment.yml`**: निर्भरताओं के लिए Conda वातावरण विनिर्देश
- **`tests/main.nf.test`**: मॉड्यूल के काम करने को मान्य करने के लिए nf-test टेस्ट केस

!!! tip "परीक्षण के बारे में और जानें"

    जेनरेट की गई टेस्ट फ़ाइल nf-test का उपयोग करती है, जो Nextflow पाइपलाइन और मॉड्यूल के लिए एक परीक्षण फ्रेमवर्क है। इन परीक्षणों को कैसे लिखें और चलाएं यह जानने के लिए, [nf-test साइड क्वेस्ट](../side_quests/nf-test.md) देखें।

जेनरेट किया गया `main.nf` उन सभी पैटर्न को शामिल करता है जो तुमने अभी सीखे, साथ ही कुछ अतिरिक्त फीचर:

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

ध्यान दें कि ऊपर तुमने मैन्युअली लागू किए गए सभी पैटर्न पहले से ही वहां हैं!

टेम्पलेट में कई अतिरिक्त nf-core परंपराएं भी शामिल हैं।
इनमें से कुछ out of the box काम करती हैं, जबकि अन्य प्लेसहोल्डर हैं जिन्हें हमें भरना होगा, जैसा कि नीचे वर्णित है।

**फीचर जो as-is काम करते हैं:**

- **`tag "$meta.id"`**: आसान ट्रैकिंग के लिए लॉग में प्रोसेस नामों में सैंपल ID जोड़ता है
- **`label 'process_single'`**: CPU/मेमोरी आवश्यकताओं को कॉन्फ़िगर करने के लिए रिसोर्स लेबल
- **`when:` ब्लॉक**: `task.ext.when` कॉन्फ़िगरेशन के माध्यम से सशर्त निष्पादन की अनुमति देता है

ये फीचर पहले से ही कार्यात्मक हैं और मॉड्यूल को अधिक मेंटेनेबल बनाते हैं।

**प्लेसहोल्डर जिन्हें हम नीचे अनुकूलित करेंगे:**

- **`input:` और `output:` ब्लॉक**: सामान्य डिक्लेरेशन जिन्हें हम अपने टूल से मेल खाने के लिए अपडेट करेंगे
- **`script:` ब्लॉक**: एक कमेंट शामिल है जहां हम `cowpy` कमांड जोड़ेंगे
- **`stub:` ब्लॉक**: टेम्पलेट जिसे हम सही आउटपुट उत्पन्न करने के लिए अपडेट करेंगे
- **Container और environment**: प्लेसहोल्डर जिन्हें हम पैकेज जानकारी से भरेंगे

अगले सेक्शन इन अनुकूलन को पूरा करने के माध्यम से चलते हैं।

### 2.2. container और conda environment सेट अप करो

nf-core दिशानिर्देशों की आवश्यकता है कि हम मॉड्यूल के हिस्से के रूप में एक container और एक Conda environment दोनों निर्दिष्ट करें।

#### 2.2.1. Container

container के लिए, तुम किसी भी Conda पैकेज से स्वचालित रूप से एक container बनाने के लिए [Seqera Containers](https://seqera.io/containers/) का उपयोग कर सकते हो, जिसमें conda-forge पैकेज शामिल हैं।
इस मामले में हम पहले की तरह ही prebuilt container का उपयोग कर रहे हैं।

डिफ़ॉल्ट कोड Docker और Singularity के बीच टॉगल करने की पेशकश करता है, लेकिन हम उस लाइन को सरल बनाने जा रहे हैं और बस Docker container निर्दिष्ट करने जा रहे हैं जो हमें ऊपर Seqera Containers से मिला।

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

Conda environment के लिए, मॉड्यूल कोड `conda "${moduleDir}/environment.yml"` निर्दिष्ट करता है जिसका अर्थ है कि इसे `environment.yml` फ़ाइल में कॉन्फ़िगर किया जाना चाहिए।

मॉड्यूल निर्माण टूल ने हमें चेतावनी दी कि वह Bioconda (बायोइन्फॉर्मेटिक्स टूल के लिए प्राथमिक चैनल) में `cowpy` पैकेज नहीं ढूंढ सका।
हालांकि, `cowpy` conda-forge में उपलब्ध है, इसलिए तुम `environment.yml` को इस तरह पूरा कर सकते हो:

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

    - **Bioconda पैकेज**: स्वचालित रूप से BioContainers बनाए जाते हैं, उपयोग के लिए तैयार container प्रदान करते हैं
    - **conda-forge पैकेज**: Conda रेसिपी से मांग पर container बनाने के लिए Seqera Containers का उपयोग कर सकते हैं

    अधिकांश बायोइन्फॉर्मेटिक्स टूल Bioconda में हैं, लेकिन conda-forge टूल के लिए, Seqera Containers कंटेनराइज़ेशन के लिए एक आसान समाधान प्रदान करता है।

### 2.3. `COWPY` लॉजिक प्लग इन करो

अब चलो उन कोड तत्वों को अपडेट करते हैं जो `COWPY` प्रोसेस क्या करती है उसके लिए विशिष्ट हैं: इनपुट और आउटपुट, और script ब्लॉक।

#### 2.3.1. इनपुट और आउटपुट

जेनरेट किया गया टेम्पलेट सामान्य इनपुट और आउटपुट डिक्लेरेशन शामिल करता है जिन्हें तुम्हें अपने विशिष्ट टूल के लिए अनुकूलित करना होगा।
सेक्शन 1 से हमारे मैनुअल `COWPY` मॉड्यूल को देखते हुए, हम उसे एक गाइड के रूप में उपयोग कर सकते हैं।

इनपुट और आउटपुट ब्लॉक अपडेट करो:

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

- इनपुट फ़ाइल पैरामीटर नाम (सामान्य `input` के बजाय `input_file`)
- कॉन्फ़िगर करने योग्य prefix पैटर्न का उपयोग करके आउटपुट फ़ाइलनाम (वाइल्डकार्ड `*` के बजाय `${prefix}.txt`)
- एक वर्णनात्मक emit नाम (सामान्य `output` के बजाय `cowpy_output`)

यदि तुम सिंटैक्स को मान्य करने के लिए Nextflow language server का उपयोग कर रहे हो, तो `${prefix}` भाग को इस स्तर पर एक एरर के रूप में फ्लैग किया जाएगा क्योंकि हमने इसे अभी तक script ब्लॉक में नहीं जोड़ा है।
चलो अब उस पर पहुंचते हैं।

#### 2.3.2. script ब्लॉक

टेम्पलेट script ब्लॉक में एक कमेंट प्लेसहोल्डर प्रदान करता है जहां तुम्हें वास्तविक टूल कमांड जोड़नी चाहिए।

पहले हमने मैन्युअली लिखे गए मॉड्यूल के आधार पर, हमें निम्नलिखित संपादन करने चाहिए:

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

मुख्य परिवर्तन:

- `def prefix` को सिर्फ `prefix` में बदलो (`def` के बिना) ताकि इसे आउटपुट ब्लॉक में सुलभ बनाया जा सके
- कमेंट को वास्तविक `cowpy` कमांड से बदलो जो `$args` और `${prefix}.txt` दोनों का उपयोग करता है

ध्यान दें कि यदि हमने पहले से ही `modules.config` फ़ाइल में `COWPY` प्रोसेस के लिए `ext.args` और `ext.prefix` कॉन्फ़िगरेशन जोड़ने का काम नहीं किया होता, तो हमें अब ऐसा करना होगा।

#### 2.3.3. stub ब्लॉक को लागू करना

Nextflow संदर्भ में, एक [stub](https://www.nextflow.io/docs/latest/process.html#stub) ब्लॉक तुम्हें वास्तविक कमांड को निष्पादित किए बिना पाइपलाइन के तर्क के तेजी से प्रोटोटाइपिंग और परीक्षण के लिए उपयोग की जाने वाली एक हल्की, डमी स्क्रिप्ट को परिभाषित करने की अनुमति देता है।

<!-- TODO (future) This is super glossed over but should really be explained or at least link out to an explanation about stubs (the reference doc isn't terribly helpful either). Right now this is likely to be mostly meaningless to anyone who doesn't already know about stubs. -->

बहुत चिंता मत करो यदि यह रहस्यमय लगता है; हम इसे पूर्णता के लिए शामिल करते हैं लेकिन तुम stub सेक्शन को भी हटा सकते हो यदि तुम इससे निपटना नहीं चाहते हो, क्योंकि यह पूरी तरह से वैकल्पिक है।

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

- script ब्लॉक से मेल खाने के लिए `def prefix` को सिर्फ `prefix` में बदलो
- `echo $args` लाइन हटाओ (जो सिर्फ टेम्पलेट प्लेसहोल्डर कोड था)
- stub एक खाली `${prefix}.txt` फ़ाइल बनाता है जो script ब्लॉक द्वारा उत्पादित की गई से मेल खाती है

यह तुम्हें वास्तविक टूल के चलने की प्रतीक्षा किए बिना workflow तर्क और फ़ाइल हैंडलिंग का परीक्षण करने की अनुमति देता है।

एक बार जब तुमने environment सेटअप (सेक्शन 2.2), इनपुट/आउटपुट (सेक्शन 2.3.1), script ब्लॉक (सेक्शन 2.3.2), और stub ब्लॉक (सेक्शन 2.3.3) पूरा कर लिया है, तो मॉड्यूल परीक्षण के लिए तैयार है!

### 2.4. नए `COWPY` मॉड्यूल में स्वैप करो और पाइपलाइन चलाओ

हमें इस नए संस्करण के `COWPY` मॉड्यूल को आज़माने के लिए बस `hello.nf` workflow फ़ाइल में इंपोर्ट स्टेटमेंट को नई फ़ाइल की ओर इशारा करने के लिए स्विच करना है।

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

चलो इसे परीक्षण करने के लिए पाइपलाइन चलाते हैं।

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

अब तुम जानते हो कि सब कुछ स्क्रैच से लिखने के बजाय टेम्पलेट का उपयोग करके कुशलता से मॉड्यूल बनाने के लिए बिल्ट-इन nf-core टूलिंग का उपयोग कैसे करें।

### आगे क्या है?

जानो कि nf-core में मॉड्यूल योगदान करने के क्या लाभ हैं और मुख्य चरण और आवश्यकताएं क्या शामिल हैं।

---

## 3. nf-core में मॉड्यूल वापस योगदान करना

[nf-core/modules](https://github.com/nf-core/modules) रिपॉजिटरी अच्छी तरह से परीक्षण किए गए, मानकीकृत मॉड्यूल के योगदान का स्वागत करती है।

### 3.1. योगदान क्यों करें?

अपने मॉड्यूल को nf-core में योगदान करना:

- [nf-co.re/modules](https://nf-co.re/modules) पर मॉड्यूल कैटलॉग के माध्यम से पूरे nf-core कम्युनिटी के लिए तुम्हारे टूल उपलब्ध कराता है
- चल रही कम्युनिटी रखरखाव और सुधार सुनिश्चित करता है
- कोड रिव्यू और स्वचालित परीक्षण के माध्यम से गुणवत्ता आश्वासन प्रदान करता है
- तुम्हारे काम को दृश्यता और मान्यता देता है

### 3.2. योगदानकर्ता की चेकलिस्ट

nf-core```markdown
में एक मॉड्यूल योगदान करने के लिए, तुम्हें निम्नलिखित चरणों से गुजरना होगा:

1. जांचो कि क्या यह पहले से [nf-co.re/modules](https://nf-co.re/modules) पर मौजूद है
2. [nf-core/modules](https://github.com/nf-core/modules) रिपॉजिटरी को फोर्क करो
3. टेम्पलेट जेनरेट करने के लिए `nf-core modules create` का उपयोग करो
4. मॉड्यूल लॉजिक और टेस्ट भरो
5. `nf-core modules test tool/subtool` के साथ परीक्षण करो
6. `nf-core modules lint tool/subtool` के साथ lint करो
7. एक pull request सबमिट करो

विस्तृत निर्देशों के लिए, [nf-core components ट्यूटोरियल](https://nf-co.re/docs/tutorials/nf-core_components/components) देखें।

### 3.3. रिसोर्स

- **Components ट्यूटोरियल**: [मॉड्यूल बनाने और योगदान करने के लिए पूर्ण गाइड](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **मॉड्यूल विनिर्देश**: [तकनीकी आवश्यकताएं और दिशानिर्देश](https://nf-co.re/docs/guidelines/components/modules)
- **कम्युनिटी सपोर्ट**: [nf-core Slack](https://nf-co.re/join) - `#modules` चैनल में शामिल हों

### सारांश

अब तुम जानते हो कि nf-core मॉड्यूल कैसे बनाए जाते हैं! तुमने चार मुख्य पैटर्न सीखे जो मॉड्यूल को पोर्टेबल और मेंटेनेबल बनाते हैं:

- **मेटाडेटा टपल** workflow के माध्यम से मेटाडेटा को प्रचारित करते हैं
- **`ext.args`** कॉन्फ़िगरेशन के माध्यम से वैकल्पिक आर्गुमेंट को संभालकर मॉड्यूल इंटरफ़ेस को सरल बनाता है
- **`ext.prefix`** आउटपुट फ़ाइल नामकरण को मानकीकृत करता है
- **केंद्रीकृत प्रकाशन** मॉड्यूल में हार्डकोड करने के बजाय `modules.config` में कॉन्फ़िगर किए गए `publishDir` के माध्यम से

`COWPY` को चरण-दर-चरण बदलकर, तुमने इन पैटर्न की गहरी समझ विकसित की, जिससे तुम nf-core मॉड्यूल के साथ काम करने, डिबग करने और बनाने के लिए सुसज्जित हो।
व्यवहार में, तुम शुरुआत से ही इन पैटर्न के साथ बने सही तरीके से संरचित मॉड्यूल जेनरेट करने के लिए `nf-core modules create` का उपयोग करोगे।

अंत में, तुमने सीखा कि nf-core कम्युनिटी में मॉड्यूल कैसे योगदान करें, टूल को दुनिया भर के शोधकर्ताओं के लिए उपलब्ध कराते हुए चल रही कम्युनिटी रखरखाव से लाभान्वित होते हुए।

### आगे क्या है?

जब तुम तैयार हो, तो [भाग 5: इनपुट सत्यापन](./05_input_validation.md) पर जारी रखो ताकि सीख सको कि अपनी पाइपलाइन में स्कीमा-आधारित इनपुट सत्यापन कैसे जोड़ें।
```
