# फ़ाइल इनपुट प्रोसेसिंग

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

वैज्ञानिक विश्लेषण workflows में अक्सर बड़ी संख्या में फ़ाइलों को प्रोसेस करना शामिल होता है।
Nextflow फ़ाइलों को कुशलता से संभालने के लिए शक्तिशाली टूल प्रदान करता है, जो तुम्हें न्यूनतम कोड के साथ अपने डेटा को व्यवस्थित और प्रोसेस करने में मदद करता है।

### सीखने के लक्ष्य

इस side quest में, हम जानेंगे कि Nextflow फ़ाइलों को कैसे संभालता है, बुनियादी फ़ाइल ऑपरेशंस से लेकर फ़ाइल संग्रहों के साथ काम करने की अधिक उन्नत तकनीकों तक।
तुम सीखोगे कि फ़ाइलनामों से मेटाडेटा कैसे निकालना है, जो वैज्ञानिक विश्लेषण pipelines में एक सामान्य आवश्यकता है।

इस side quest के अंत तक, तुम निम्नलिखित करने में सक्षम होंगे:

- Nextflow की `file()` मेथड का उपयोग करके फ़ाइल path strings से Path ऑब्जेक्ट्स बनाना
- फ़ाइल विशेषताओं जैसे name, extension, और parent डायरेक्टरी तक पहुँचना
- URIs का उपयोग करके लोकल और रिमोट फ़ाइलों को पारदर्शी रूप से संभालना
- `channel.fromPath()` और `channel.fromFilePairs()` के साथ फ़ाइल हैंडलिंग को स्वचालित करने के लिए channels का उपयोग करना
- string manipulation का उपयोग करके फ़ाइलनामों से मेटाडेटा निकालना और संरचित करना
- पैटर्न मैचिंग और glob expressions का उपयोग करके संबंधित फ़ाइलों को ग्रुप करना
- उचित इनपुट हैंडलिंग के साथ फ़ाइल ऑपरेशंस को Nextflow processes में इंटीग्रेट करना
- मेटाडेटा-संचालित डायरेक्टरी संरचनाओं का उपयोग करके process आउटपुट को व्यवस्थित करना

ये कौशल तुम्हें ऐसे workflows बनाने में मदद करेंगे जो विभिन्न प्रकार के फ़ाइल इनपुट को बहुत लचीलेपन के साथ संभाल सकते हैं।

### पूर्वापेक्षाएँ

इस side quest को शुरू करने से पहले, तुम्हें:

- [Hello Nextflow](../../hello_nextflow/) ट्यूटोरियल या समकक्ष बिगिनर कोर्स पूरा कर लेना चाहिए।
- बुनियादी Nextflow अवधारणाओं और मैकेनिज्म (processes, channels, operators) का उपयोग करने में सहज होना चाहिए

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. शुरू करें

#### प्रशिक्षण codespace खोलें

यदि तुमने अभी तक नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित अनुसार प्रशिक्षण वातावरण खोलना सुनिश्चित करें।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाएं

चलो उस डायरेक्टरी में जाते हैं जहाँ इस ट्यूटोरियल की फ़ाइलें स्थित हैं।

```bash
cd side-quests/working_with_files
```

तुम VSCode को इस डायरेक्टरी पर फोकस करने के लिए सेट कर सकते हो:

```bash
code .
```

#### सामग्री की समीक्षा करें

तुम्हें `main.nf` नामक एक सरल workflow फ़ाइल, दो मॉड्यूल फ़ाइलों वाली `modules` डायरेक्टरी, और कुछ उदाहरण डेटा फ़ाइलों वाली `data` डायरेक्टरी मिलेगी।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    .
    ├── data
    │   ├── patientA_rep1_normal_R1_001.fastq.gz
    │   ├── patientA_rep1_normal_R2_001.fastq.gz
    │   ├── patientA_rep1_tumor_R1_001.fastq.gz
    │   ├── patientA_rep1_tumor_R2_001.fastq.gz
    │   ├── patientA_rep2_normal_R1_001.fastq.gz
    │   ├── patientA_rep2_normal_R2_001.fastq.gz
    │   ├── patientA_rep2_tumor_R1_001.fastq.gz
    │   ├── patientA_rep2_tumor_R2_001.fastq.gz
    │   ├── patientB_rep1_normal_R1_001.fastq.gz
    │   ├── patientB_rep1_normal_R2_001.fastq.gz
    │   ├── patientB_rep1_tumor_R1_001.fastq.gz
    │   ├── patientB_rep1_tumor_R2_001.fastq.gz
    │   ├── patientC_rep1_normal_R1_001.fastq.gz
    │   ├── patientC_rep1_normal_R2_001.fastq.gz
    │   ├── patientC_rep1_tumor_R1_001.fastq.gz
    │   └── patientC_rep1_tumor_R2_001.fastq.gz
    ├── main.nf
    └── modules
        ├── analyze_reads.nf
        └── count_lines.nf
    ```

इस डायरेक्टरी में तीन मरीजों (A, B, C) से paired-end सीक्वेंसिंग डेटा है।

प्रत्येक मरीज के लिए, हमारे पास नमूने हैं जो `tumor` प्रकार के हैं (आमतौर पर ट्यूमर बायोप्सी से उत्पन्न) या `normal` (स्वस्थ ऊतक या रक्त से लिए गए)।
यदि तुम कैंसर विश्लेषण से परिचित नहीं हो, बस जान लो कि यह एक प्रयोगात्मक मॉडल से मेल खाता है जो कंट्रास्टिव विश्लेषण करने के लिए paired tumor/normal नमूनों का उपयोग करता है।

विशेष रूप से मरीज A के लिए, हमारे पास तकनीकी प्रतिकृतियों (दोहराव) के दो सेट हैं।

सीक्वेंसिंग डेटा फ़ाइलों को एक सामान्य `_R1_` और `_R2_` कन्वेंशन के साथ नामित किया गया है जो 'forward reads' और 'reverse reads' के रूप में जानी जाती हैं।

_अगर तुम इस प्रयोगात्मक डिज़ाइन से परिचित नहीं हो तो चिंता मत करो, यह इस ट्यूटोरियल को समझने के लिए महत्वपूर्ण नहीं है।_

#### असाइनमेंट की समीक्षा करें

तुम्हारी चुनौती एक Nextflow workflow लिखना है जो:

1. Nextflow की फ़ाइल हैंडलिंग मेथड्स का उपयोग करके इनपुट फ़ाइलों को **लोड** करेगा
2. फ़ाइलनाम संरचना से मेटाडेटा (patient ID, replicate, sample type) **निकालेगा**
3. `channel.fromFilePairs()` का उपयोग करके paired फ़ाइलों (R1/R2) को एक साथ **ग्रुप** करेगा
4. प्रदान किए गए विश्लेषण मॉड्यूल के साथ फ़ाइलों को **प्रोसेस** करेगा
5. निकाले गए मेटाडेटा के आधार पर आउटपुट को डायरेक्टरी संरचना में **व्यवस्थित** करेगा

#### तैयारी चेकलिस्ट

क्या तुम शुरू करने के लिए तैयार हो?

- [ ] मैं इस कोर्स का लक्ष्य और इसकी पूर्वापेक्षाएँ समझता/समझती हूँ
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी working डायरेक्टरी उचित रूप से सेट कर ली है
- [ ] मैं असाइनमेंट समझता/समझती हूँ

यदि तुम सभी बॉक्स चेक कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. बुनियादी फ़ाइल ऑपरेशंस

### 1.1. `.class` के साथ ऑब्जेक्ट के प्रकार की पहचान करें

workflow फ़ाइल `main.nf` पर एक नज़र डालो:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // एक string path से Path ऑब्जेक्ट बनाएं
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

यह एक मिनी-workflow है (बिना किसी processes के) जो अपने workflow में एकल फ़ाइल path को संदर्भित करता है, फिर इसे कंसोल पर प्रिंट करता है, इसकी class के साथ।

??? info "`.class` क्या है?"

    Nextflow में, `.class` हमें बताता है कि हम किस प्रकार के ऑब्जेक्ट के साथ काम कर रहे हैं। यह पूछने जैसा है "यह किस तरह की चीज़ है?" यह पता लगाने के लिए कि यह एक string है, एक number है, एक फ़ाइल है, या कुछ और।
    यह हमें अगले अनुभागों में एक plain string और एक Path ऑब्जेक्ट के बीच अंतर को स्पष्ट करने में मदद करेगा।

चलो workflow चलाते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

जैसा कि तुम देख सकते हो, Nextflow ने string path को ठीक वैसे ही प्रिंट किया जैसा हमने लिखा था।

यह सिर्फ टेक्स्ट आउटपुट है; Nextflow ने अभी तक इसके साथ कुछ विशेष नहीं किया है।
हमने यह भी पुष्टि कर ली है कि जहाँ तक Nextflow का संबंध है, यह केवल एक string है (class `java.lang.String` की)।
यह समझ में आता है, क्योंकि हमने अभी तक Nextflow को नहीं बताया है कि यह एक फ़ाइल से मेल खाती है।

### 1.2. file() के साथ Path ऑब्जेक्ट बनाएं

हम path strings से [Path ऑब्जेक्ट्स](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) बनाकर Nextflow को बता सकते हैं कि फ़ाइलों को कैसे संभालना है।

हमारे workflow में, हम `file()` मेथड का उपयोग करके string path `data/patientA_rep1_normal_R1_001.fastq.gz` को Path ऑब्जेक्ट में बदल सकते हैं, जो फ़ाइल प्रॉपर्टीज़ और ऑपरेशंस तक पहुँच प्रदान करती है।

`main.nf` को एडिट करें और string को `file()` के साथ wrap करें जैसा कि नीचे दिखाया गया है:

=== "बाद"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // एक string path से Path ऑब्जेक्ट बनाएं
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // एक string path से Path ऑब्जेक्ट बनाएं
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

अब workflow फिर से चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

इस बार, तुम देखते हो कि हमारे द्वारा इनपुट के रूप में प्रदान किए गए relative path के बजाय full absolute path दिख रहा है।

Nextflow ने हमारी string को Path ऑब्जेक्ट में बदल दिया है और इसे सिस्टम पर वास्तविक फ़ाइल स्थान पर resolve कर दिया है।
फ़ाइल path अब absolute होगा, जैसे `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`।

यह भी ध्यान दो कि Path ऑब्जेक्ट class `sun.nio.fs.UnixPath` है: यह Nextflow का लोकल फ़ाइलों का प्रतिनिधित्व करने का तरीका है।
जैसा कि हम बाद में देखेंगे, रिमोट फ़ाइलों की अलग-अलग class names होंगी (जैसे HTTP फ़ाइलों के लिए `nextflow.file.http.XPath`), लेकिन वे सभी ठीक उसी तरह काम करती हैं और तुम्हारे workflows में समान रूप से उपयोग की जा सकती हैं।

!!! tip "सुझाव"

    **मुख्य अंतर:**

    - **Path string**: बस टेक्स्ट जिसे Nextflow अक्षरों के रूप में मानता है
    - **Path ऑब्जेक्ट**: एक स्मार्ट फ़ाइल रेफरेंस जिसके साथ Nextflow काम कर सकता है

    इसे इस तरह सोचो: एक path string कागज पर पता लिखने जैसा है, जबकि Path ऑब्जेक्ट GPS डिवाइस में लोड किए गए पते जैसा है जो जानता है कि वहाँ कैसे नेविगेट करना है और यात्रा के बारे में विवरण बता सकता है।

### 1.3. फ़ाइल विशेषताओं तक पहुँचें

यह क्यों सहायक है? खैर, अब जब Nextflow समझ गया है कि `myFile` एक Path ऑब्जेक्ट है न कि सिर्फ एक string, हम Path ऑब्जेक्ट की विभिन्न विशेषताओं तक पहुँच सकते हैं।

चलो अपने workflow को अपडेट करते हैं ताकि built-in फ़ाइल विशेषताएँ प्रिंट हों:

=== "बाद"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // एक string path से Path ऑब्जेक्ट बनाएं
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ प्रिंट करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // एक string path से Path ऑब्जेक्ट बनाएं
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

तुम ऊपर कंसोल पर प्रिंट की गई विभिन्न फ़ाइल विशेषताएँ देख सकते हो।

### 1.4. फ़ाइल को process में फीड करें

strings और Path ऑब्जेक्ट्स के बीच का अंतर तब महत्वपूर्ण हो जाता है जब तुम processes के साथ वास्तविक workflows बनाना शुरू करते हो।
अब तक हमने सत्यापित किया है कि Nextflow अब हमारी इनपुट फ़ाइल को फ़ाइल के रूप में मान रहा है, लेकिन देखते हैं कि क्या हम वास्तव में उस फ़ाइल पर एक process में कुछ चला सकते हैं।

#### 1.4.1. process को import करें और कोड की जाँच करें

हम तुम्हें `COUNT_LINES` नामक एक पूर्व-लिखित process मॉड्यूल प्रदान करते हैं जो एक फ़ाइल इनपुट लेता है और गिनता है कि इसमें कितनी लाइनें हैं।

workflow में process का उपयोग करने के लिए, तुम्हें workflow block से पहले एक include स्टेटमेंट जोड़ना होगा:

=== "बाद"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

तुम मॉड्यूल फ़ाइल को खोलकर इसके कोड की जाँच कर सकते हो:

```groovy title="modules/count_lines.nf" linenums="1"
#!/usr/bin/env nextflow

process COUNT_LINES {
    debug true

    input:
    path input_file

    script:
    """
    set -o pipefail
    echo "Processing file: $input_file"
    gzip -dc $input_file | wc -l
    """
}
```

जैसा कि तुम देख सकते हो, यह एक काफी सीधी छोटी script है जो फ़ाइल को unzip करती है और गिनती करती है कि इसमें कितनी लाइनें हैं।

??? info "`debug true` क्या करता है?"

    process definition में `debug true` निर्देश Nextflow को तुम्हारी script से आउटपुट (जैसे line count "40") सीधे execution log में प्रिंट करने का कारण बनता है।
    इसके बिना, तुम केवल process execution status देखोगे लेकिन तुम्हारी script से वास्तविक आउटपुट नहीं।

    Nextflow processes को debug करने के बारे में अधिक जानकारी के लिए, [Debugging Nextflow Workflows](debugging.md) side quest देखें।

#### 1.4.2. `COUNT_LINES` को कॉल जोड़ें

अब जब process workflow के लिए उपलब्ध है, हम इनपुट फ़ाइल पर इसे चलाने के लिए `COUNT_LINES` process को कॉल जोड़ सकते हैं।

workflow में निम्नलिखित एडिट करें:

=== "बाद"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // एक string path से Path ऑब्जेक्ट बनाएं
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ प्रिंट करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // फ़ाइल में लाइनें गिनें
        COUNT_LINES(myFile)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // एक string path से Path ऑब्जेक्ट बनाएं
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ प्रिंट करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

और अब workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_hypatia] DSL2 - revision: 281d13c414

    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    executor >  local (1)
    [e9/341c05] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

यह दिखाता है कि हम एक process के अंदर फ़ाइल पर उचित रूप से ऑपरेट करने में सक्षम हैं।

विशेष रूप से, Nextflow ने निम्नलिखित ऑपरेशंस सफलतापूर्वक किए:

- फ़ाइल को working डायरेक्टरी में stage किया
- .gz फ़ाइल को decompress किया
- लाइनें गिनीं (इस मामले में 40 लाइनें)
- बिना error के पूरा हुआ

इस सुगम ऑपरेशन की कुंजी यह है कि हम Nextflow को स्पष्ट रूप से बता रहे हैं कि हमारा इनपुट एक फ़ाइल है और इसे इस प्रकार माना जाना चाहिए।

### 1.5. बुनियादी फ़ाइल इनपुट errors का समस्या निवारण

यह अक्सर Nextflow के नए उपयोगकर्ताओं को परेशान करता है, तो चलो कुछ मिनट लेते हैं यह देखने के लिए कि जब तुम इसे गलत करते हो तो क्या होता है।

दो मुख्य स्थान हैं जहाँ तुम फ़ाइल हैंडलिंग गलत कर सकते हो: workflow के स्तर पर, और process के स्तर पर।

#### 1.5.1. Workflow-स्तर error

देखते हैं क्या होता है अगर हम workflow block में इनपुट निर्दिष्ट करते समय फ़ाइल को string के रूप में मानते हैं।

workflow में निम्नलिखित एडिट करें, path-विशिष्ट print statements को comment out करना सुनिश्चित करें:

=== "बाद"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // एक string path से Path ऑब्जेक्ट बनाएं
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // फ़ाइल विशेषताएँ प्रिंट करें
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // फ़ाइल में लाइनें गिनें
        COUNT_LINES(myFile)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // एक string path से Path ऑब्जेक्ट बनाएं
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ प्रिंट करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // फ़ाइल में लाइनें गिनें
        COUNT_LINES(myFile)
    ```

और अब workflow चलाओ:

```bash
nextflow run main.nf
```

??? failure "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_goodall] DSL2 - revision: ae50609b20

    [-        ] COUNT_LINES -
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'



    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

यह महत्वपूर्ण हिस्सा है:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

जब तुम `path` इनपुट निर्दिष्ट करते हो, Nextflow सत्यापित करता है कि तुम वास्तविक फ़ाइल references पास कर रहे हो, न कि सिर्फ strings।
यह error तुम्हें बता रहा है कि `'data/patientA_rep1_normal_R1_001.fastq.gz'` एक valid path value नहीं है क्योंकि यह एक string है, Path ऑब्जेक्ट नहीं।

Nextflow ने तुरंत समस्या का पता लगाया और process शुरू करने से पहले ही रुक गया।

#### 1.5.2. Process-स्तर error

दूसरी जगह जहाँ हम निर्दिष्ट करना भूल सकते हैं कि हम चाहते हैं कि Nextflow इनपुट को फ़ाइल के रूप में माने, वह process definition में है।

!!! warning "चेतावनी"

    इस टेस्ट के सही ढंग से काम करने के लिए, workflow को उसकी टूटी हुई अवस्था में रखो (`file()` के बजाय plain string का उपयोग करके)।
    जब process में `val` के साथ जोड़ा जाता है, तो यह नीचे दिखाया गया error उत्पन्न करता है।

मॉड्यूल में निम्नलिखित एडिट करें:

=== "बाद"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "पहले"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

और अब workflow फिर से चलाओ:

```bash
nextflow run main.nf
```

??? failure "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_golick] DSL2 - revision: ae50609b20

    executor >  local (1)
    [b3/b3023c] COUNT_LINES [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Process `COUNT_LINES` terminated with an error exit status (1)


    Command executed:

      set -o pipefail
      echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
      gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l

    Command exit status:
      1

    Command output:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      0

    Command error:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
      0

    Work dir:
      /workspaces/training/side-quests/working_with_files/work/b3/b3023cb2ccb986851301d8e369e79f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

यह error के बारे में बहुत सारे विवरण दिखाता है क्योंकि process debugging जानकारी आउटपुट करने के लिए सेट है, जैसा कि ऊपर नोट किया गया है।

ये सबसे प्रासंगिक अनुभाग हैं:

```console
Command executed:

  set -o pipefail
  echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
  gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l
```

```console
Command error:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
  0
```

यह कहता है कि सिस्टम फ़ाइल नहीं ढूंढ सका; हालाँकि अगर तुम path देखो, उस नाम की एक फ़ाइल उस स्थान पर है।

जब हमने यह चलाया, Nextflow ने string value को script में pass किया, लेकिन इसने वास्तविक फ़ाइल को working डायरेक्टरी में _stage_ नहीं किया।
इसलिए process ने relative string, `data/patientA_rep1_normal_R1_001.fastq.gz` का उपयोग करने की कोशिश की, लेकिन वह फ़ाइल process working डायरेक्टरी के भीतर मौजूद नहीं है।

इन दोनों उदाहरणों को मिलाकर, तुम देख सकते हो कि Nextflow को बताना कितना महत्वपूर्ण है कि क्या किसी इनपुट को फ़ाइल के रूप में संभाला जाना चाहिए।

!!! note "नोट"

    अगले अनुभाग पर जाने से पहले दोनों जानबूझकर errors को ठीक करना सुनिश्चित करो।

### सीख

- Path strings बनाम Path ऑब्जेक्ट्स: Strings बस टेक्स्ट हैं, Path ऑब्जेक्ट्स स्मार्ट फ़ाइल references हैं
- `file()` मेथड एक string path को Path ऑब्जेक्ट में बदलती है जिसके साथ Nextflow काम कर सकता है
- तुम [file attributes का उपयोग करके](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes) फ़ाइल प्रॉपर्टीज़ जैसे `name`, `simpleName`, `extension`, और `parent` तक पहुँच सकते हो
- strings के बजाय Path ऑब्जेक्ट्स का उपयोग Nextflow को तुम्हारे workflow में फ़ाइलों को ठीक से प्रबंधित करने की अनुमति देता है
- Process इनपुट परिणाम: उचित फ़ाइल हैंडलिंग के लिए strings नहीं, Path ऑब्जेक्ट्स की आवश्यकता होती है, यह सुनिश्चित करने के लिए कि फ़ाइलें processes द्वारा उपयोग के लिए सही ढंग से staged और accessible हैं।

---

## 2. रिमोट फ़ाइलों का उपयोग

Nextflow की प्रमुख विशेषताओं में से एक लोकल फ़ाइलों (उसी मशीन पर) से इंटरनेट पर accessible रिमोट फ़ाइलों पर seamlessly switch करने की क्षमता है।

यदि तुम इसे सही ढंग से कर रहे हो, तो विभिन्न स्थानों से आने वाली फ़ाइलों को accommodate करने के लिए तुम्हें अपने workflow के लॉजिक को कभी बदलने की आवश्यकता नहीं होनी चाहिए।
रिमोट फ़ाइल का उपयोग करने के लिए तुम्हें बस इतना करना है कि फ़ाइल path में उचित prefix निर्दिष्ट करो जब तुम इसे workflow को supply कर रहे हो।

उदाहरण के लिए, `/path/to/data` में कोई prefix नहीं है, जो दर्शाता है कि यह एक 'normal' लोकल फ़ाइल path है, जबकि `s3://path/to/data` में `s3://` prefix शामिल है, जो दर्शाता है कि यह Amazon के S3 ऑब्जेक्ट स्टोरेज में स्थित है।

कई अलग-अलग protocols समर्थित हैं:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

इनमें से किसी का भी उपयोग करने के लिए, बस string में relevant prefix निर्दिष्ट करो, जिसे तकनीकी रूप से फ़ाइल path के बजाय Uniform Resource Identifier (URI) कहा जाता है।
Nextflow authentication और फ़ाइलों को सही जगह पर staging को संभालेगा, downloading या uploading और अन्य सभी फ़ाइल ऑपरेशंस जो तुम expect करोगे।

इस सिस्टम की मुख्य ताकत यह है कि यह हमें किसी भी pipeline लॉजिक को बदले बिना environments के बीच switch करने में सक्षम बनाता है।
उदाहरण के लिए, तुम एक छोटे, लोकल टेस्ट सेट के साथ develop कर सकते हो और फिर बस URI बदलकर रिमोट स्टोरेज में स्थित full-scale टेस्ट सेट पर switch कर सकते हो।

### 2.1. इंटरनेट से फ़ाइल का उपयोग करें

चलो इसे टेस्ट करते हैं अपने workflow को दिए जा रहे लोकल path को एक HTTPS path से बदलकर जो Github में stored उसी डेटा की एक कॉपी को point करता है।

!!! warning "चेतावनी"

    यह तभी काम करेगा जब तुम्हारे पास active इंटरनेट कनेक्शन है।

`main.nf` फिर से खोलो और इनपुट path को निम्नानुसार बदलो:

=== "बाद"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // इंटरनेट से रिमोट फ़ाइल का उपयोग
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ प्रिंट करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // एक string path से Path ऑब्जेक्ट बनाएं
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ प्रिंट करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

चलो workflow चलाते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    File object class: class nextflow.file.http.XPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /nextflow-io/training/master/side-quests/working_with_files/data
    executor >  local (1)
    [8a/2ab7ca] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

यह काम करता है! तुम देख सकते हो कि बहुत कम बदला है।

कंसोल आउटपुट में एक अंतर यह है कि path object class अब `nextflow.file.http.XPath` है, जबकि लोकल path के लिए class `sun.nio.fs.UnixPath` थी।
तुम्हें इन classes को याद रखने की आवश्यकता नहीं है; हम इसका उल्लेख बस यह प्रदर्शित करने के लिए करते हैं कि Nextflow विभिन्न स्थानों को उचित रूप से पहचानता और संभालता है।

पर्दे के पीछे, Nextflow ने फ़ाइल को work डायरेक्टरी के भीतर स्थित staging डायरेक्टरी में download किया।
उस staged फ़ाइल को फिर लोकल फ़ाइल के रूप में माना जा सकता है और relevant process डायरेक्टरी में symlinked किया जा सकता है।

तुम process के hash value पर स्थित working डायरेक्टरी की सामग्री देखकर यह सत्यापित कर सकते हो कि ऐसा हुआ।

??? abstract "Work डायरेक्टरी सामग्री"

    यदि process hash `8a/2ab7ca` था, तुम work डायरेक्टरी explore कर सकते हो:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    Symlink रिमोट फ़ाइल की staged copy को point करता है जिसे Nextflow ने automatically download किया।

ध्यान दो कि बड़ी फ़ाइलों के लिए, downloading step लोकल फ़ाइलों पर चलने की तुलना में कुछ अतिरिक्त समय लेगा।
हालाँकि, Nextflow check करता है कि क्या उसके पास पहले से staged copy है ताकि अनावश्यक downloads से बचा जा सके।
इसलिए यदि तुम उसी फ़ाइल पर फिर से run करते हो और staged फ़ाइल को delete नहीं किया है, Nextflow staged copy का उपयोग करेगा।

यह दिखाता है कि Nextflow का उपयोग करके लोकल और रिमोट डेटा के बीच switch करना कितना आसान है, जो Nextflow की एक प्रमुख विशेषता है।

!!! note "नोट"

    इस सिद्धांत का एक महत्वपूर्ण अपवाद यह है कि तुम HTTPS के साथ glob patterns या डायरेक्टरी paths का उपयोग नहीं कर सकते क्योंकि HTTPS multiple फ़ाइलों को list नहीं कर सकता, इसलिए तुम्हें exact फ़ाइल URLs निर्दिष्ट करने होंगे।
    हालाँकि, अन्य storage protocols जैसे blob storage (`s3://`, `az://`, `gs://`) globs और डायरेक्टरी paths दोनों का उपयोग कर सकते हैं।

    यहाँ बताया गया है कि तुम cloud storage के साथ glob patterns कैसे उपयोग कर सकते हो:

    ```groovy title="Cloud storage examples (इस वातावरण में runnable नहीं)"
    // S3 with glob patterns - multiple फ़ाइलों को match करेगा
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage with glob patterns
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage with glob patterns
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    हम तुम्हें अगले अनुभाग में practice में globs के साथ काम करना दिखाएंगे।

### 2.2. लोकल फ़ाइल पर वापस switch करें

हम इस side quest के बाकी हिस्से के लिए अपनी लोकल example फ़ाइलों का उपयोग करने जा रहे हैं, तो चलो workflow इनपुट को वापस original फ़ाइल पर switch करते हैं:

=== "बाद"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // एक string path से Path ऑब्जेक्ट बनाएं
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ प्रिंट करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // एक string path से Path ऑब्जेक्ट बनाएं
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ प्रिंट करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### सीख

- रिमोट डेटा को URI (HTTP, FTP, S3, Azure, Google Cloud) का उपयोग करके access किया जाता है
- Nextflow automatically डेटा को download और सही जगह पर stage करेगा, जब तक ये paths processes को feed किए जा रहे हैं
- रिमोट फ़ाइलों को download या upload करने के लिए लॉजिक मत लिखो!
- लोकल और रिमोट फ़ाइलें अलग-अलग object types produce करती हैं लेकिन समान रूप से काम करती हैं
- **महत्वपूर्ण**: HTTP/HTTPS केवल single फ़ाइलों के साथ काम करते हैं (glob patterns नहीं)
- Cloud storage (S3, Azure, GCS) single फ़ाइलों और glob patterns दोनों को support करता है
- तुम code लॉजिक बदले बिना लोकल और रिमोट डेटा sources के बीच seamlessly switch कर सकते हो (जब तक protocol तुम्हारे required ऑपरेशंस को support करता है)

---

## 3. `fromPath()` channel factory का उपयोग

अब तक हम एक समय में एक फ़ाइल के साथ काम कर रहे थे, लेकिन Nextflow में, हम आमतौर पर प्रोसेस करने के लिए multiple इनपुट फ़ाइलों के साथ एक इनपुट channel बनाना चाहेंगे।

ऐसा करने का एक naive तरीका `file()` मेथड को [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) के साथ combine करना होगा जैसे:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

यह काम करता है, लेकिन यह अजीब है।

!!! tip "सुझाव"

    कब `file()` बनाम `channel.fromPath()` का उपयोग करें

    - `file()` का उपयोग करो जब तुम्हें direct manipulation के लिए single Path ऑब्जेक्ट की आवश्यकता हो (check करना कि फ़ाइल exists है, इसकी attributes पढ़ना, या single process invocation को pass करना)
    - `channel.fromPath()` का उपयोग करो जब तुम्हें ऐसा channel चाहिए जो multiple फ़ाइलें hold कर सके, विशेष रूप से glob patterns के साथ, या जब फ़ाइलें multiple processes से flow करेंगी

यहाँ [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) आता है: एक सुविधाजनक channel factory जो एक या अधिक static फ़ाइल strings के साथ-साथ glob patterns से channel generate करने के लिए हमें आवश्यक सभी functionality bundle करती है।

### 3.1. channel factory जोड़ें

चलो अपने workflow को अपडेट करते हैं `channel.fromPath` का उपयोग करने के लिए।

=== "बाद"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // फ़ाइल विशेषताएँ प्रिंट करें
        /* इन्हें अभी comment out करो, हम इन पर वापस आएंगे!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // फ़ाइल में लाइनें गिनें
        // COUNT_LINES(myFile)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // एक string path से Path ऑब्जेक्ट बनाएं
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ प्रिंट करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // फ़ाइल में लाइनें गिनें
        COUNT_LINES(myFile)
    ```

हमने अभी के लिए attributes print करने वाले code को भी comment out कर दिया है, और इसके बजाय बस filename print करने के लिए `.view` statement जोड़ा है।

workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

जैसा कि तुम देख सकते हो, फ़ाइल path channel में `Path` type ऑब्जेक्ट के रूप में load हो रहा है।
यह उसी के समान है जो `file()` ने किया होता, सिवाय इसके कि अब हमारे पास एक channel है जिसमें हम चाहें तो और फ़ाइलें load कर सकते हैं।

`channel.fromPath()` का उपयोग फ़ाइलों की सूची द्वारा populated नया channel बनाने का एक सुविधाजनक तरीका है।

### 3.2. channel में फ़ाइलों की attributes देखें

Channel factory का उपयोग करने के हमारे पहले प्रयास में, हमने code को सरल बनाया और बस फ़ाइल का नाम print किया।

चलो पूर्ण फ़ाइल attributes print करने पर वापस जाते हैं:

=== "बाद"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // फ़ाइल में लाइनें गिनें
        COUNT_LINES(ch_files)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // फ़ाइल में लाइनें गिनें
        // COUNT_LINES(ch_files)
    ```

हम `COUNT_LINES` process call को भी re-enable कर रहे हैं यह verify करने के लिए कि फ़ाइल processing अभी भी हमारे channel-based approach के साथ सही ढंग से काम करती है।

workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [furious_swanson] DSL2 - revision: c35c34950d

    executor >  local (1)
    [9d/6701a6] COUNT_LINES (1) [100%] 1 of 1 ✔
    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

और वहाँ तुम हो, पहले जैसे ही परिणाम लेकिन अब हमारे पास फ़ाइल एक channel में है, इसलिए हम और जोड़ सकते हैं।

### 3.3. multiple फ़ाइलों को match करने के लिए glob का उपयोग

channel में और फ़ाइलें load करने के कई तरीके हैं।
यहाँ हम तुम्हें glob patterns का उपयोग करना दिखाने जा रहे हैं, जो wildcard characters के आधार पर फ़ाइल और डायरेक्टरी names को match और retrieve करने का एक सुविधाजनक तरीका है।
इन patterns को match करने की प्रक्रिया को "globbing" या "filename expansion" कहा जाता है।

!!! note "नोट"

    जैसा कि पहले noted किया गया था, Nextflow अधिकांश मामलों में इनपुट और आउटपुट फ़ाइलों को manage करने के लिए globbing को support करता है, HTTPS filepaths के साथ छोड़कर क्योंकि HTTPS multiple फ़ाइलों को list नहीं कर सकता।

मान लो हम किसी दिए गए patient, `patientA` से जुड़ी फ़ाइलों की pair में दोनों फ़ाइलें retrieve करना चाहते हैं:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

चूँकि filenames के बीच एकमात्र अंतर replicate number है, _यानी_ `R` के बाद का number, हम wildcard character `*` का उपयोग number के स्थान पर निम्नानुसार कर सकते हैं:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

यह वह glob pattern है जिसकी हमें आवश्यकता है।

अब हमें बस channel factory में फ़ाइल path को उस glob pattern का उपयोग करने के लिए update करना है जैसा कि निम्नानुसार है:

=== "बाद"

    ```groovy title="main.nf" linenums="7"
      // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7"
      // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow automatically पहचान लेगा कि यह एक glob pattern है और इसे appropriately handle करेगा।

इसे test करने के लिए workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

    executor >  local (2)
    [3c/a65de5] COUNT_LINES (2) [100%] 2 of 2 ✔
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R2_001.fastq.gz
    Simple name: patientA_rep1_normal_R2_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

जैसा कि तुम देख सकते हो, अब हमारे channel में दो Path ऑब्जेक्ट्स हैं, जो दिखाता है कि Nextflow ने filename expansion सही ढंग से की है, और दोनों फ़ाइलों को expected के अनुसार load और process किया है।

इस मेथड का उपयोग करके, हम बस glob pattern बदलकर जितनी चाहें उतनी या कम फ़ाइलें retrieve कर सकते हैं। यदि हम इसे अधिक generous बनाते हैं, उदाहरण के लिए filenames के सभी variable parts को `*` से replace करके (_जैसे_ `data/patient*_rep*_*_R*_001.fastq.gz`) हम `data` डायरेक्टरी में सभी example फ़ाइलें grab कर सकते हैं।

### सीख

- `channel.fromPath()` एक pattern से match करने वाली फ़ाइलों के साथ channel बनाता है
- प्रत्येक फ़ाइल channel में एक separate element के रूप में emit होती है
- हम multiple फ़ाइलों को match करने के लिए glob pattern का उपयोग कर सकते हैं
- फ़ाइलें automatically full attributes के साथ Path ऑब्जेक्ट्स में convert हो जाती हैं
- `.view()` मेथड channel contents का inspection करने की अनुमति देती है

---

## 4. फ़ाइलनामों से बुनियादी मेटाडेटा निकालना

अधिकांश वैज्ञानिक domains में, डेटा वाली फ़ाइलों के नामों में मेटाडेटा encoded होना बहुत आम है।
उदाहरण के लिए, bioinformatics में, सीक्वेंसिंग डेटा वाली फ़ाइलों को अक्सर इस तरह नाम दिया जाता है जो नमूने, condition, replicate, और read number के बारे में जानकारी encode करता है।

यदि filenames एक consistent convention के अनुसार constructed हैं, तो तुम उस मेटाडेटा को standardized manner में extract कर सकते हो और अपने विश्लेषण के दौरान इसका उपयोग कर सकते हो।
यह एक बड़ा 'अगर' है, बेशक, और जब भी तुम filename structure पर rely करते हो तो तुम्हें बहुत सावधान रहना चाहिए; लेकिन वास्तविकता यह है कि यह approach बहुत व्यापक रूप से उपयोग किया जाता है, तो चलो देखते हैं कि Nextflow में यह कैसे किया जाता है।

हमारे example डेटा के मामले में, हम जानते हैं कि filenames में consistently structured मेटाडेटा शामिल है।
उदाहरण के लिए, filename `patientA_rep1_normal_R2_001` निम्नलिखित encode करता है:

- patient ID: `patientA`
- replicate ID: `rep1`
- sample type: `normal` (`tumor` के विपरीत)
- read set: `R1` (`R2` के विपरीत)

हम अपने workflow को तीन चरणों में इस जानकारी को retrieve करने के लिए modify करने जा रहे हैं:

1. फ़ाइल का `simpleName` retrieve करें, जिसमें मेटाडेटा शामिल है
2. `tokenize()` नामक मेथड का उपयोग करके मेटाडेटा को separate करें
3. मेटाडेटा को organize करने के लिए map का उपयोग करें

!!! warning "चेतावनी"

    तुम्हें कभी भी filenames में sensitive जानकारी encode नहीं करनी चाहिए, जैसे patient names या अन्य identifying characteristics, क्योंकि यह patient privacy या अन्य relevant security restrictions को compromise कर सकता है।

### 4.1. `simpleName` retrieve करें

`simpleName` एक फ़ाइल attribute है जो path और extension से stripped filename से correspond करता है।

workflow में निम्नलिखित edits करें:

=== "बाद"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

यह `simpleName` retrieve करता है और इसे `map()` ऑपरेशन का उपयोग करके full file ऑब्जेक्ट के साथ associate करता है।

यह काम करता है इसे test करने के लिए workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_mahavira] DSL2 - revision: ae8edc4e48

    executor >  local (2)
    [e9/55774b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [patientA_rep1_normal_R2_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [patientA_rep1_normal_R1_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

channel में प्रत्येक element अब एक tuple है जिसमें `simpleName` और original file ऑब्जेक्ट है।

### 4.2. `simplename` से मेटाडेटा extract करें

इस point पर, हमें जो मेटाडेटा चाहिए वह `simplename` में embedded है, लेकिन हम individual items को directly access नहीं कर सकते।
इसलिए हमें `simplename` को इसके components में split करना होगा।
सौभाग्य से, वे components original filename में बस underscores द्वारा separated हैं, इसलिए हम `tokenize()` नामक एक common Nextflow मेथड apply कर सकते हैं जो इस कार्य के लिए perfect है।

workflow में निम्नलिखित edits करें:

=== "बाद"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

`tokenize()` मेथड `simpleName` string को जहाँ भी underscores पाती है वहाँ split करेगी, और substrings वाली list return करेगी।

workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [gigantic_gauss] DSL2 - revision: a39baabb57

    executor >  local (2)
    [e7/da2f4b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [[patientA, rep1, normal, R2, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[patientA, rep1, normal, R1, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

अब हमारे channel में प्रत्येक element के लिए tuple में मेटाडेटा की list (_जैसे_ `[patientA, rep1, normal, R1, 001]`) और original file ऑब्जेक्ट है।

बहुत बढ़िया!
हमने अपनी patient जानकारी को single string से strings की list में break down किया है।
अब हम patient जानकारी के प्रत्येक भाग को separately handle कर सकते हैं।

### 4.3. मेटाडेटा को organize करने के लिए map का उपयोग करें

हमारा मेटाडेटा इस समय बस एक flat list है।
इसे use करना आसान है लेकिन पढ़ना मुश्किल है।

```console
[patientA, rep1, normal, R1, 001]
```

index 3 पर item क्या है? क्या तुम मेटाडेटा structure की original explanation को refer किए बिना बता सकते हो?

यह एक key-value store का उपयोग करने का एक बढ़िया अवसर है, जहाँ प्रत्येक item में keys और उनके associated values का set होता है, इसलिए तुम corresponding value प्राप्त करने के लिए आसानी से प्रत्येक key को refer कर सकते हो।

हमारे example में, इसका मतलब है इस organization से जाना:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

इस तक:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

Nextflow में, इसे [map](https://nextflow.io/docs/latest/script.html#maps) कहा जाता है।

चलो अब अपनी flat list को map में convert करते हैं।
workflow में निम्नलिखित edits करें:

=== "बाद"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
              [
                id: patient,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum.replace('R', ''),
              ],
              myFile
            ]
        }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

यहाँ मुख्य changes हैं:

- **Destructuring assignment**: `def (patient, replicate, type, readNum) = ...` tokenized values को एक line में named variables में extract करता है
- **Map literal syntax**: `[id: patient, replicate: ...]` एक map बनाता है जहाँ प्रत्येक key (जैसे `id`) एक value (जैसे `patient`) से associated है
- **Nested structure**: Outer list `[..., myFile]` मेटाडेटा map को original file ऑब्जेक्ट के साथ pair करती है

हमने `replace()` नामक string replacement मेथड का उपयोग करके कुछ मेटाडेटा strings को भी simplify किया ताकि कुछ unnecessary characters remove हों (_जैसे_ `replicate.replace('rep', '')` replicate IDs से केवल number रखने के लिए)।

चलो workflow फिर से चलाते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

    executor >  local (2)
    [1b/e7fb27] COUNT_LINES (1) [100%] 2 of 2 ✔
    [[id:patientA, replicate:1, type:normal, readNum:2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[id:patientA, replicate:1, type:normal, readNum:1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

अब मेटाडेटा neatly labeled है (_जैसे_ `[id:patientA, replicate:1, type:normal, readNum:2]`) इसलिए यह बताना बहुत आसान है कि क्या है।

workflow में मेटाडेटा के elements का actually उपयोग करना भी बहुत आसान होगा, और हमारे code को पढ़ने में आसान और अधिक maintainable बनाएगा।

### सीख

- हम Nextflow में filenames को एक full programming language की power के साथ handle कर सकते हैं
- हम relevant जानकारी extract करने के लिए filenames को strings के रूप में treat कर सकते हैं
- `tokenize()` और `replace()` जैसी मेथड्स का उपयोग हमें filename में strings को manipulate करने की अनुमति देता है
- `.map()` ऑपरेशन structure को preserve करते हुए channel elements को transform करता है
- Structured मेटाडेटा (maps) positional lists की तुलना में code को अधिक readable और maintainable बनाता है

अगला, हम देखेंगे कि paired डेटा फ़ाइलों को कैसे handle करना है।

---

## 5. paired डेटा फ़ाइलों को handle करना

कई experimental designs paired डेटा फ़ाइलें produce करते हैं जो explicitly paired तरीके से handle किए जाने से लाभान्वित होती हैं।
उदाहरण के लिए, bioinformatics में, सीक्वेंसिंग डेटा अक्सर paired reads के रूप में generate होता है, जिसका मतलब sequence strings जो DNA के same fragment से originate होती हैं (अक्सर 'forward' और 'reverse' कहा जाता है क्योंकि वे opposite ends से read की जाती हैं)।

यही हमारे example डेटा का मामला है, जहाँ R1 और R2 reads के दो sets को refer करते हैं।

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow इस तरह की paired फ़ाइलों के साथ काम करने के लिए एक specialized channel factory प्रदान करता है जिसे `channel.fromFilePairs()` कहा जाता है, जो automatically shared naming pattern के आधार पर फ़ाइलों को group करती है। यह तुम्हें कम प्रयास के साथ paired फ़ाइलों को अधिक tightly associate करने की अनुमति देता है।

हम इसका लाभ उठाने के लिए अपने workflow को modify करने जा रहे हैं।
इसमें दो steps लगेंगे:

1. channel factory को `channel.fromFilePairs()` में switch करें
2. मेटाडेटा extract और map करें

### 5.1. channel factory को `channel.fromFilePairs()` में switch करें

`channel.fromFilePairs` का उपयोग करने के लिए, हमें वह pattern specify करना होगा जिसका Nextflow को pair में दो members identify करने के लिए उपयोग करना चाहिए।

हमारे example डेटा पर वापस जाते हुए, हम naming pattern को निम्नानुसार formalize कर सकते हैं:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

यह उस glob pattern के समान है जो हमने पहले उपयोग किया था, सिवाय इसके कि यह specifically उन substrings को enumerate करता है (या तो `1` या `2` R के ठीक बाद आ रहा है) जो pair के दो members को identify करते हैं।

चलो `main.nf` workflow को accordingly update करते हैं:

=== "बाद"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // channel.fromFilePairs के साथ फ़ाइलें लोड करें
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* अभी के लिए mapping comment out करो, हम इस पर वापस आएंगे!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // channel.fromFilePairs के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        .view()
    ```

हमने channel factory switch की है और file matching pattern adapt किया है, और जब हम इस पर थे, हमने map ऑपरेशन को comment out कर दिया।
हम इसे बाद में वापस add करेंगे, कुछ modifications के साथ।

इसे test करने के लिए workflow चलाओ:

```bash
nextflow run main.nf
```

??? failure "कमांड आउटपुट"

    ```console hl_lines="7-8"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_koch] DSL2 - revision: 44fdf66105

    [-        ] COUNT_LINES -
    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ERROR ~ Error executing process > 'COUNT_LINES (1)'

    Caused by:
      Not a valid path value: 'patientA_rep1_normal_R'



    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

अरे, इस बार run fail हो गया!

error message का relevant हिस्सा यहाँ है:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

ऐसा इसलिए है क्योंकि हमने channel factory बदल दी है।
अब तक, original input channel में केवल file paths थे।
हम जो सभी मेटाडेटा manipulation कर रहे थे उसने actually channel contents को affect नहीं किया।

अब जब हम `.fromFilePairs` channel factory का उपयोग कर रहे हैं, resulting channel की contents different हैं।
हम केवल एक channel element देखते हैं, जो दो items वाले tuple से composed है: दो फ़ाइलों द्वारा shared `simpleName` का हिस्सा, जो एक identifier के रूप में serve करता है, और दो file ऑब्जेक्ट्स वाला tuple, format `id, [ file1, file2 ]` में।

यह बढ़िया है, क्योंकि Nextflow ने shared prefix examine करके और इसे patient identifier के रूप में उपयोग करके patient name extract करने का कठिन काम किया है।

हालाँकि, यह हमारे current workflow को break करता है।
यदि हम अभी भी `COUNT_LINES` को process बदले बिना same तरीके से run करना चाहते थे, तो हमें file paths extract करने के लिए mapping ऑपरेशन apply करना होगा।
लेकिन हम ऐसा नहीं करने जा रहे हैं, क्योंकि हमारा ultimate goal एक different process, `ANALYZE_READS` का उपयोग करना है, जो file pairs को appropriately handle करता है।

तो चलो बस `COUNT_LINES` की call को comment out (या delete) करते हैं और आगे बढ़ते हैं।

=== "बाद"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // फ़ाइल में लाइनें गिनें
        // COUNT_LINES(ch_files)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // फ़ाइल में लाइनें गिनें
        COUNT_LINES(ch_files)
    ```

तुम `COUNT_LINES` include statement को भी comment out या delete कर सकते हो, लेकिन उसका कोई functional effect नहीं होगा।

अब चलो workflow फिर से चलाते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

इस बार workflow succeed होता है!

हालाँकि, हमें अभी भी `id` field से बाकी मेटाडेटा निकालना है।

### 5.2. file pairs से मेटाडेटा extract और organize करें

पहले का हमारा `map` ऑपरेशन काम नहीं करेगा क्योंकि यह data structure से match नहीं करता, लेकिन हम इसे काम करने के लिए modify कर सकते हैं।

हमारे पास पहले से ही actual patient identifier का access उस string में है जो `fromFilePairs()` ने identifier के रूप में उपयोग किया, इसलिए हम इसका उपयोग मेटाडेटा extract करने के लिए कर सकते हैं बिना Path ऑब्जेक्ट से `simpleName` प्राप्त किए जैसा हमने पहले किया था।

workflow में map ऑपरेशन को uncomment करो और निम्नलिखित edits करो:

=== "बाद"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // channel.fromFilePairs के साथ फ़ाइलें लोड करें
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id, files ->
            def (sample, replicate, type) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                files
            ]
        }
        .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // channel.fromFilePairs के साथ फ़ाइलें लोड करें
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* अभी के लिए mapping comment out करो, हम इस पर वापस आएंगे!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

इस बार map बस `myFile` के बजाय `id, files` से शुरू होता है, और `tokenize()` `myFile.simpleName` के बजाय `id` पर apply होता है।

यह भी notice करो कि हमने `tokenize()` line से `readNum` हटा दिया है; कोई भी substrings जिन्हें हम specifically name नहीं करते (left से शुरू करके) silently drop हो जाएंगी।
हम यह इसलिए कर सकते हैं क्योंकि paired फ़ाइलें अब tightly associated हैं, इसलिए हमें metadata map में `readNum` की अब आवश्यकता नहीं है।

चलो workflow चलाते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

और वहाँ है: हमारे पास metadata map (`[id:patientA, replicate:1, type:normal]`) output tuple की first position में है, इसके बाद paired फ़ाइलों का tuple है, जैसा intended था।

बेशक, यह केवल उस specific pair of फ़ाइलों को pick up और process करेगा।
यदि तुम multiple pairs को process करने के साथ experiment करना चाहते हो, तुम input pattern में wildcards add करने की कोशिश कर सकते हो और देख सकते हो क्या होता है।
उदाहरण के लिए, `data/patientA_rep1_*_R{1,2}_001.fastq.gz` try करो।

### सीख

- [`channel.fromFilePairs()` automatically related फ़ाइलों को find और pair करती है](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- यह तुम्हारी pipeline में paired-end reads को handle करना simplify करता है
- Paired फ़ाइलों को `[id, [file1, file2]]` tuples के रूप में group किया जा सकता है
- मेटाडेटा extraction individual फ़ाइलों के बजाय paired file ID से की जा सकती है

---

## 6. processes में file ऑपरेशंस का उपयोग

अब चलो यह सब एक simple process में एक साथ रखते हैं यह reinforce करने के लिए कि Nextflow process के अंदर file ऑपरेशंस कैसे उपयोग करें।

हम तुम्हें `ANALYZE_READS` नामक एक pre-written process मॉड्यूल प्रदान करते हैं जो मेटाडेटा और इनपुट फ़ाइलों की pair का tuple लेता है और उनका analyze करता है।
हम imagine कर सकते हैं कि यह sequence alignment, या variant calling या कोई अन्य step कर रहा है जो इस data type के लिए sense बनाता है।

चलो शुरू करते हैं।

### 6.1. process को import करें और code examine करें

workflow में इस process का उपयोग करने के लिए, हमें बस workflow block से पहले एक module include statement add करना होगा।

workflow में निम्नलिखित edit करो:

=== "बाद"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

तुम मॉड्यूल फ़ाइल को खोलकर इसके code examine कर सकते हो:

```groovy title="modules/analyze_reads.nf - process example" linenums="1"
#!/usr/bin/env nextflow

process ANALYZE_READS {
    tag { meta.id }

    publishDir { "results/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta.id), path("${meta.id}_stats.txt")

    script:
    """
    echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Read 1: ${files[0]}" >> ${meta.id}_stats.txt
    echo "Read 2: ${files[1]}" >> ${meta.id}_stats.txt
    echo "Read 1 size: \$(gunzip -dc ${files[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    echo "Read 2 size: \$(gunzip -dc ${files[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    """
}
```

!!! note "नोट"

    `tag` और `publishDir` निर्देश string interpolation (`"${...}"`) के बजाय closure syntax (`{ ... }`) का उपयोग करते हैं।
    ऐसा इसलिए है क्योंकि ये निर्देश input variables (`meta`) को reference करते हैं जो runtime तक available नहीं हैं।
    Closure syntax evaluation को तब तक defer करता है जब तक process actually run नहीं होता।

!!! note "नोट"

    हम अपने metadata map को convention से `meta` कह रहे हैं।
    Meta maps में deeper dive के लिए, [Metadata and meta maps](./metadata.md) side quest देखें।

### 6.2. workflow में process को call करें

अब जब process workflow के लिए available है, हम इसे run करने के लिए `ANALYZE_READS` process को call add कर सकते हैं।

इसे हमारे example डेटा पर run करने के लिए, हमें दो चीजें करनी होंगी:

1. Remapped channel को name दें
2. Process को call add करें

#### 6.2.1. Remapped input channel को name दें

हमने पहले mapping manipulations को directly input channel पर apply किया था।
Remapped contents को `ANALYZE_READS` process को feed करने के लिए (और इसे इस तरह से करने के लिए जो clear और पढ़ने में easy हो) हम `ch_samples` नाम का एक नया channel बनाना चाहते हैं।

हम यह [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) ऑपरेटर का उपयोग करके कर सकते हैं।

Main workflow में, `.view()` ऑपरेटर को `.set { ch_samples }` से replace करो, और एक line add करो testing कि हम channel को name से refer कर सकते हैं।

=== "बाद"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // channel.fromFilePairs के साथ फ़ाइलें लोड करें
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
            .set { ch_samples }

        // Temporary: ch_samples में peek करें
        ch_samples.view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // channel.fromFilePairs के साथ फ़ाइलें लोड करें
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
        .view()
    }
    ```

चलो यह run करते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

यह confirm करता है कि हम अब channel को name से refer कर सकते हैं।

#### 6.2.2. data पर process call करें

अब चलो actually `ch_samples` channel पर `ANALYZE_READS` process call करते हैं।

Main workflow में, निम्नलिखित code changes करो:

=== "बाद"

    ```groovy title="main.nf" linenums="23"
        // analysis चलाएं
        ANALYZE_READS(ch_samples)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="23"
        // Temporary: ch_samples में peek करें
        ch_samples.view()
    ```

चलो यह run करते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

यह process अपने outputs को `results` डायरेक्टरी में publish करने के लिए set up है, तो वहाँ एक नज़र डालो।

??? abstract "डायरेक्टरी और फ़ाइल contents"

    ```console
    results
    └── patientA
        └── patientA_stats.txt
    ```

    ```txt title="patientA_stats.txt"
    Sample metadata: patientA
    Replicate: 1
    Type: normal
    Read 1: patientA_rep1_normal_R1_001.fastq.gz
    Read 2: patientA_rep1_normal_R2_001.fastq.gz
    Read 1 size: 10 reads
    Read 2 size: 10 reads
    ```

Process ने हमारे inputs लिए और patient metadata वाली एक new फ़ाइल बनाई, जैसा designed था।
शानदार!

### 6.3. कई और patients include करें

बेशक, यह single patient के लिए बस single pair of फ़ाइलों को process कर रहा है, जो exactly उस तरह का high throughput नहीं है जिसकी तुम Nextflow के साथ उम्मीद कर रहे हो।
तुम शायद एक समय में बहुत अधिक डेटा process करना चाहोगे।

याद रखो `channel.fromPath()` input के रूप में एक _glob_ accept करता है, जिसका मतलब है कि यह pattern से match करने वाली किसी भी संख्या में फ़ाइलें accept कर सकता है।
इसलिए यदि हम सभी patients को include करना चाहते हैं, हम बस अधिक patients include करने के लिए input string modify कर सकते हैं, जैसा कि पहले passing में note किया गया था।

मान लो हम जितना possible हो उतना greedy होना चाहते हैं।
workflow में निम्नलिखित edits करो:

=== "बाद"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // channel.fromFilePairs के साथ फ़ाइलें लोड करें
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // channel.fromFilePairs के साथ फ़ाइलें लोड करें
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Pipeline फिर से run करो:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Results डायरेक्टरी में अब सभी available डेटा के लिए results होने चाहिए।

??? abstract "डायरेक्टरी contents"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

सफलता! हमने सभी patients को एक बार में analyze कर दिया! सही?

शायद नहीं।
यदि तुम अधिक closely देखो, हमारे पास एक problem है: patientA के लिए हमारे पास दो replicates हैं, लेकिन केवल एक output फ़ाइल!
हम हर बार output फ़ाइल को overwrite कर रहे हैं।

### 6.4. published फ़ाइलों को unique बनाएं

चूँकि हमारे पास patient metadata तक access है, हम इसका उपयोग published फ़ाइलों को unique बनाने के लिए कर सकते हैं differentiating metadata को include करके, या तो directory structure में या filenames में ही।

workflow में निम्नलिखित change करो:

=== "बाद"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "पहले"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

यहाँ हम sample types और replicates के account के लिए additional directory levels उपयोग करने का option दिखाते हैं, लेकिन तुम filename level पर भी इसे करने के साथ experiment कर सकते हो।

अब pipeline एक और बार run करो, लेकिन पहले results डायरेक्टरी remove करना सुनिश्चित करो ताकि तुम्हें clean workspace मिले:

```bash
rm -r results
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

अब results डायरेक्टरी check करो:

??? abstract "डायरेक्टरी contents"

    ```console
    results/
    ├── normal
    │   ├── patientA
    │   │   ├── 1
    │   │   │   └── patientA_stats.txt
    │   │   └── 2
    │   │       └── patientA_stats.txt
    │   ├── patientB
    │   │   └── 1
    │   │       └── patientB_stats.txt
    │   └── patientC
    │       └── 1
    │           └── patientC_stats.txt
    └── tumor
        ├── patientA
        │   ├── 1
        │   │   └── patientA_stats.txt
        │   └── 2
        │       └── patientA_stats.txt
        ├── patientB
        │   └── 1
        │       └── patientB_stats.txt
        └── patientC
            └── 1
                └── patientC_stats.txt
    ```

और वहाँ है, हमारा सारा मेटाडेटा, neatly organized। यह सफलता है!

एक बार जब तुम्हारा मेटाडेटा इस तरह map में loaded है तो और भी बहुत कुछ कर सकते हो:

1. Patient attributes के आधार पर organized output directories बनाना
2. Patient properties के आधार पर processes में decisions लेना
3. Metadata values के आधार पर data को split, join, और recombine करना

मेटाडेटा को explicit और data से attached रखने का यह pattern (filenames में encoded करने के बजाय) Nextflow में एक core best practice है जो robust, maintainable analysis workflows बनाने में enable करता है।
तुम इसके बारे में [Metadata and meta maps](./metadata.md) side quest में और जान सकते हो।

### सीख

- `publishDir` निर्देश metadata values के आधार पर outputs organize कर सकता है
- Tuples में metadata results का structured organization enable करता है
- यह approach clear data provenance के साथ maintainable workflows बनाता है
- Processes मेटाडेटा और फ़ाइलों के tuples को input के रूप में ले सकते हैं
- `tag` निर्देश execution logs में process identification प्रदान करता है
- Workflow structure channel creation को process execution से separate करती है

---

## सारांश

इस side quest में, तुमने सीखा कि Nextflow में फ़ाइलों के साथ कैसे काम करना है, basic ऑपरेशंस से लेकर फ़ाइलों के collections को handle करने की अधिक advanced techniques तक।

अपने own काम में इन techniques को apply करना तुम्हें अधिक efficient और maintainable workflows बनाने में enable करेगा, खासकर जब complex naming conventions के साथ बड़ी संख्या में फ़ाइलों के साथ काम कर रहे हो।

### मुख्य पैटर्न

1.  **Basic File ऑपरेशंस:** हमने `file()` के साथ Path ऑब्जेक्ट्स बनाए और name, extension, और parent directory जैसी फ़ाइल attributes access कीं, strings और Path ऑब्जेक्ट्स के बीच अंतर सीखा।

    - `file()` के साथ Path ऑब्जेक्ट बनाएं

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - फ़ाइल attributes प्राप्त करें

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Remote Files का उपयोग**: हमने सीखा कि URIs का उपयोग करके local और remote फ़ाइलों के बीच transparently switch कैसे करना है, Nextflow की विभिन्न sources से फ़ाइलों को workflow logic बदले बिना handle करने की क्षमता demonstrate करते हुए।

    - Local फ़ाइल

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - FTP

    ```groovy
    myFile = file('ftp://path/to/file.txt')
    ```

    - HTTPS

    ```groovy
    myFile = file('https://path/to/file.txt')
    ```

    - Amazon S3

    ```groovy
    myFile = file('s3://path/to/file.txt')
    ```

    - Azure Blob Storage

    ```groovy
    myFile = file('az://path/to/file.txt')
    ```

    - Google Cloud Storage

    ```groovy
    myFile = file('gs://path/to/file.txt')
    ```

3.  **`fromPath()` channel factory का उपयोग करके फ़ाइलें load करना:** हमने `channel.fromPath()` के साथ file patterns से channels बनाए और object types सहित उनकी file attributes देखीं।

    - File pattern से channel बनाएं

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - फ़ाइल attributes प्राप्त करें

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Filenames से Patient Metadata Extract करना:** हमने filenames से metadata extract और structure करने के लिए `tokenize()` और `replace()` का उपयोग किया, उन्हें organized maps में convert करते हुए।

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **channel.fromFilePairs के साथ Simplify करना:** हमने related फ़ाइलों को automatically pair करने और paired file IDs से metadata extract करने के लिए `channel.fromFilePairs()` का उपयोग किया।

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Processes में File ऑपरेशंस का उपयोग:** हमने proper input handling के साथ file ऑपरेशंस को Nextflow processes में integrate किया, metadata के आधार पर outputs organize करने के लिए `publishDir` का उपयोग करते हुए।

    - Process inputs के साथ meta map associate करें

    ```groovy
    ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_files.map { id,  files ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
             files
        ]
    }
        .set { ch_samples }

    ANALYZE_READS(ch_samples)
    ```

    - Metadata के आधार पर outputs organize करें

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### अतिरिक्त संसाधन

- [Nextflow Documentation: Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## आगे क्या?

[Side Quests के menu](./index.md) पर वापस जाओ या सूची में अगले topic पर जाने के लिए page के bottom right में button पर click करो।
