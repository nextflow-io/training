# फ़ाइल इनपुट प्रोसेसिंग

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

वैज्ञानिक विश्लेषण वर्कफ़्लो में अक्सर बड़ी संख्या में फ़ाइलों को प्रोसेस करना होता है।
Nextflow फ़ाइलों को कुशलतापूर्वक संभालने के लिए शक्तिशाली टूल प्रदान करता है, जो तुम्हें न्यूनतम कोड के साथ अपने डेटा को व्यवस्थित और प्रोसेस करने में मदद करता है।

### सीखने के लक्ष्य

इस side quest में, हम देखेंगे कि Nextflow फ़ाइलों को कैसे संभालता है, बुनियादी फ़ाइल ऑपरेशन से लेकर फ़ाइल संग्रह के साथ काम करने की अधिक उन्नत तकनीकों तक।
तुम सीखोगे कि फ़ाइलनामों से मेटाडेटा कैसे निकाला जाए, जो वैज्ञानिक विश्लेषण पाइपलाइन में एक सामान्य आवश्यकता है।

इस side quest के अंत तक, तुम सक्षम होगे:

- Nextflow के `file()` मेथड का उपयोग करके फ़ाइल पाथ strings से Path objects बनाना
- name, extension, और parent directory जैसी फ़ाइल विशेषताओं तक पहुँचना
- URIs का उपयोग करके local और remote फ़ाइलों को पारदर्शी रूप से संभालना
- `channel.fromPath()` और `channel.fromFilePairs()` के साथ फ़ाइल हैंडलिंग को स्वचालित करने के लिए channels का उपयोग करना
- string manipulation का उपयोग करके फ़ाइलनामों से मेटाडेटा निकालना और संरचित करना
- pattern matching और glob expressions का उपयोग करके संबंधित फ़ाइलों को ग्रुप करना
- उचित इनपुट हैंडलिंग के साथ Nextflow processes में फ़ाइल ऑपरेशन को एकीकृत करना
- मेटाडेटा-आधारित डायरेक्टरी संरचनाओं का उपयोग करके process आउटपुट को व्यवस्थित करना

ये कौशल तुम्हें ऐसे वर्कफ़्लो बनाने में मदद करेंगे जो विभिन्न प्रकार के फ़ाइल इनपुट को बड़ी लचीलेपन के साथ संभाल सकते हैं।

### पूर्वापेक्षाएँ

इस side quest को शुरू करने से पहले, तुम्हें:

- [Hello Nextflow](../../hello_nextflow/) ट्यूटोरियल या समकक्ष शुरुआती कोर्स पूरा करना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (processes, channels, operators) का उपयोग करने में सहज होना चाहिए।

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. शुरू करना

#### Training codespace खोलो

अगर तुमने अभी तक नहीं किया है, तो [Environment Setup](../envsetup/index.md) में बताए अनुसार training वातावरण खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाओ

चलो उस डायरेक्टरी में जाते हैं जहाँ इस ट्यूटोरियल की फ़ाइलें हैं।

```bash
cd side-quests/working_with_files
```

तुम VSCode को इस डायरेक्टरी पर फ़ोकस करने के लिए सेट कर सकते हो:

```bash
code .
```

#### सामग्री की समीक्षा करो

तुम्हें `main.nf` नामक एक सरल वर्कफ़्लो फ़ाइल, दो module फ़ाइलें वाली एक `modules` डायरेक्टरी, और कुछ उदाहरण डेटा फ़ाइलें वाली एक `data` डायरेक्टरी मिलेगी।

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

इस डायरेक्टरी में तीन मरीज़ों (A, B, C) के paired-end sequencing डेटा हैं।

प्रत्येक मरीज़ के लिए, हमारे पास `tumor` प्रकार के नमूने हैं (आमतौर पर tumor biopsies से) या `normal` (स्वस्थ ऊतक या रक्त से लिए गए)।
अगर तुम cancer analysis से परिचित नहीं हो, तो बस यह जानो कि यह एक प्रयोगात्मक मॉडल से मेल खाता है जो contrastive analyses करने के लिए paired tumor/normal नमूनों का उपयोग करता है।

विशेष रूप से patient A के लिए, हमारे पास दो sets के technical replicates (दोहराव) हैं।

Sequencing डेटा फ़ाइलों को 'forward reads' और 'reverse reads' के लिए एक सामान्य `_R1_` और `_R2_` convention के साथ नाम दिया गया है।

_चिंता मत करो अगर तुम इस experimental design से परिचित नहीं हो, यह इस ट्यूटोरियल को समझने के लिए महत्वपूर्ण नहीं है।_

#### असाइनमेंट की समीक्षा करो

तुम्हारी चुनौती एक Nextflow वर्कफ़्लो लिखना है जो:

1. Nextflow के फ़ाइल हैंडलिंग मेथड का उपयोग करके इनपुट फ़ाइलें **लोड** करे
2. फ़ाइलनाम संरचना से मेटाडेटा (patient ID, replicate, sample type) **निकाले**
3. `channel.fromFilePairs()` का उपयोग करके paired फ़ाइलें (R1/R2) **ग्रुप** करे
4. एक प्रदान किए गए analysis module के साथ फ़ाइलों को **प्रोसेस** करे
5. निकाले गए मेटाडेटा के आधार पर डायरेक्टरी संरचना में आउटपुट **व्यवस्थित** करे

#### तैयारी की जाँच सूची

क्या तुम शुरू करने के लिए तैयार हो?

- [ ] मैं इस कोर्स के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता/समझती हूँ
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी working directory उचित रूप से सेट की है
- [ ] मैं असाइनमेंट को समझता/समझती हूँ

अगर तुम सभी बॉक्स चेक कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. बुनियादी फ़ाइल ऑपरेशन

### 1.1. `.class` के साथ किसी object के प्रकार की पहचान करो

वर्कफ़्लो फ़ाइल `main.nf` देखो:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // एक string path से Path object बनाएँ
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

यह एक mini-workflow है (बिना किसी process के) जो अपने workflow में एक single फ़ाइल पाथ को refer करता है, फिर उसे console पर उसके class के साथ print करता है।

??? info "`.class` क्या है?"

    Nextflow में, `.class` हमें बताता है कि हम किस प्रकार के object के साथ काम कर रहे हैं। यह पूछने जैसा है "यह किस तरह की चीज़ है?" यह जानने के लिए कि यह एक string है, एक number है, एक फ़ाइल है, या कुछ और।
    यह हमें अगले sections में एक plain string और एक Path object के बीच का अंतर दिखाने में मदद करेगा।

चलो वर्कफ़्लो चलाते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

जैसा तुम देख सकते हो, Nextflow ने string path को ठीक वैसे ही print किया जैसा हमने लिखा था।

यह सिर्फ text आउटपुट है; Nextflow ने अभी तक इसके साथ कुछ खास नहीं किया है।
हमने यह भी पुष्टि की है कि जहाँ तक Nextflow का संबंध है, यह केवल एक string है (class `java.lang.String` की)।
यह समझ में आता है, क्योंकि हमने अभी तक Nextflow को नहीं बताया है कि यह एक फ़ाइल से मेल खाता है।

### 1.2. `file()` के साथ Path object बनाओ

हम Nextflow को फ़ाइलों को संभालने का तरीका बता सकते हैं, path strings से [Path objects](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) बनाकर।

हमारे वर्कफ़्लो में, हम `file()` मेथड का उपयोग करके string path `data/patientA_rep1_normal_R1_001.fastq.gz` को एक Path object में बदल सकते हैं, जो फ़ाइल properties और operations तक पहुँच प्रदान करता है।

`main.nf` को edit करो और string को `file()` से इस प्रकार wrap करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // एक string path से Path object बनाएँ
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // एक string path से Path object बनाएँ
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

अब वर्कफ़्लो फिर से चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

इस बार, तुम हमारे द्वारा इनपुट के रूप में दिए गए relative path के बजाय पूरा absolute path देखते हो।

Nextflow ने हमारी string को एक Path object में बदल दिया है और इसे सिस्टम पर actual फ़ाइल location पर resolve कर दिया है।
फ़ाइल path अब absolute होगा, जैसे `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`।

यह भी ध्यान दो कि Path object class `sun.nio.fs.UnixPath` है: यह Nextflow का local फ़ाइलों को represent करने का तरीका है।
जैसा कि हम बाद में देखेंगे, remote फ़ाइलों के अलग class names होंगे (जैसे HTTP फ़ाइलों के लिए `nextflow.file.http.XPath`), लेकिन वे सभी बिल्कुल एक ही तरह से काम करते हैं और तुम्हारे वर्कफ़्लो में समान रूप से उपयोग किए जा सकते हैं।

!!! tip "सुझाव"

    **मुख्य अंतर:**

    - **Path string**: सिर्फ text जिसे Nextflow characters के रूप में treat करता है
    - **Path object**: एक smart फ़ाइल reference जिसके साथ Nextflow काम कर सकता है

    इसे इस तरह सोचो: एक path string कागज़ पर पता लिखने जैसा है, जबकि एक Path object GPS device में पता लोड करने जैसा है जो वहाँ navigate करना जानता है और यात्रा के बारे में विवरण बता सकता है।

### 1.3. फ़ाइल विशेषताओं तक पहुँचो

यह क्यों उपयोगी है? अब जब Nextflow समझता है कि `myFile` एक Path object है न कि सिर्फ एक string, हम Path object की विभिन्न विशेषताओं तक पहुँच सकते हैं।

चलो अपने वर्कफ़्लो को built-in फ़ाइल विशेषताओं को print करने के लिए update करते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // एक string path से Path object बनाएँ
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ print करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // एक string path से Path object बनाएँ
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

वर्कफ़्लो चलाओ:

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

तुम ऊपर console पर print की गई विभिन्न फ़ाइल विशेषताएँ देखते हो।

### 1.4. फ़ाइल को एक process में फ़ीड करो

strings और Path objects के बीच का अंतर तब महत्वपूर्ण हो जाता है जब तुम processes के साथ actual वर्कफ़्लो बनाना शुरू करते हो।
अब तक हमने verify किया है कि Nextflow अब हमारी इनपुट फ़ाइल को एक फ़ाइल के रूप में treat कर रहा है, लेकिन देखते हैं कि क्या हम वास्तव में किसी process में उस फ़ाइल पर कुछ चला सकते हैं।

#### 1.4.1. Process को import करो और code की जाँच करो

हम तुम्हें `COUNT_LINES` नामक एक pre-written process module प्रदान करते हैं जो एक फ़ाइल इनपुट लेता है और उसमें कितनी lines हैं यह गिनता है।

वर्कफ़्लो में process का उपयोग करने के लिए, तुम्हें बस workflow block से पहले एक include statement जोड़नी होगी:

=== "बाद में"

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

तुम module फ़ाइल खोलकर उसके code की जाँच कर सकते हो:

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

जैसा तुम देख सकते हो, यह एक काफी सीधा-सादा script है जो फ़ाइल को unzip करता है और उसमें कितनी lines हैं यह गिनता है।

??? info "`debug true` क्या करता है?"

    process definition में `debug true` directive Nextflow को तुम्हारे script के आउटपुट (जैसे line count "40") को सीधे execution log में print करने का कारण बनता है।
    इसके बिना, तुम केवल process execution status देखोगे लेकिन तुम्हारे script का actual आउटपुट नहीं।

    Nextflow processes को debug करने के बारे में अधिक जानकारी के लिए, [Debugging Nextflow Workflows](debugging.md) side quest देखो।

#### 1.4.2. `COUNT_LINES` को call करो

अब जब process वर्कफ़्लो के लिए उपलब्ध है, हम इनपुट फ़ाइल पर इसे चलाने के लिए `COUNT_LINES` process को call कर सकते हैं।

वर्कफ़्लो में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // एक string path से Path object बनाएँ
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ print करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // फ़ाइल में lines गिनें
        COUNT_LINES(myFile)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // एक string path से Path object बनाएँ
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ print करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

और अब वर्कफ़्लो चलाओ:

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

यह दिखाता है कि हम एक process के अंदर फ़ाइल पर उचित रूप से काम करने में सक्षम हैं।

विशेष रूप से, Nextflow ने निम्नलिखित ऑपरेशन सफलतापूर्वक किए:

- फ़ाइल को working directory में stage किया
- .gz फ़ाइल को decompress किया
- Lines गिनीं (इस मामले में 40 lines)
- बिना किसी error के पूरा किया

इस smooth operation की कुंजी यह है कि हम Nextflow को स्पष्ट रूप से बता रहे हैं कि हमारा इनपुट एक फ़ाइल है और इसे ऐसे ही treat किया जाना चाहिए।

### 1.5. बुनियादी फ़ाइल इनपुट errors को troubleshoot करो

यह अक्सर Nextflow के नए users को परेशान करता है, इसलिए चलो कुछ मिनट देखते हैं कि जब तुम इसे गलत करते हो तो क्या होता है।

दो मुख्य जगहें हैं जहाँ तुम फ़ाइल हैंडलिंग गलत कर सकते हो: वर्कफ़्लो के स्तर पर, और process के स्तर पर।

#### 1.5.1. Workflow-level error

देखते हैं कि क्या होता है अगर हम workflow block में इनपुट specify करते समय फ़ाइल को string के रूप में treat करने पर वापस जाते हैं।

वर्कफ़्लो में निम्नलिखित बदलाव करो, path-specific print statements को comment out करना सुनिश्चित करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // एक string path से Path object बनाएँ
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // फ़ाइल विशेषताएँ print करें
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // फ़ाइल में lines गिनें
        COUNT_LINES(myFile)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // एक string path से Path object बनाएँ
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ print करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // फ़ाइल में lines गिनें
        COUNT_LINES(myFile)
    ```

और अब वर्कफ़्लो चलाओ:

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

जब तुम `path` इनपुट specify करते हो, Nextflow validate करता है कि तुम actual फ़ाइल references pass कर रहे हो, न कि सिर्फ strings।
यह error तुम्हें बता रही है कि `'data/patientA_rep1_normal_R1_001.fastq.gz'` एक valid path value नहीं है क्योंकि यह एक string है, Path object नहीं।

Nextflow ने तुरंत समस्या detect की और process शुरू करने से पहले ही रुक गया।

#### 1.5.2. Process-level error

दूसरी जगह जहाँ हम यह specify करना भूल सकते हैं कि हम Nextflow को इनपुट को फ़ाइल के रूप में treat करना चाहते हैं, वह process definition में है।

!!! warning "1.5.1 से workflow error को रखो"

    इस test के सही ढंग से काम करने के लिए, workflow को उसकी broken state में रखो (plain string का उपयोग करते हुए `file()` के बजाय)।
    Process में `val` के साथ combine होने पर, यह नीचे दिखाई गई error produce करता है।

Module में निम्नलिखित बदलाव करो:

=== "बाद में"

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

और अब वर्कफ़्लो फिर से चलाओ:

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

यह error के बारे में बहुत सारे विवरण दिखाता है क्योंकि process ऊपर बताए अनुसार debugging जानकारी output करने के लिए सेट है।

ये सबसे relevant sections हैं:

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

यह कहता है कि सिस्टम फ़ाइल नहीं ढूंढ सका; हालाँकि अगर तुम path देखो, तो उस location पर उस नाम की एक फ़ाइल है।

जब हमने यह चलाया, Nextflow ने string value को script में pass किया, लेकिन उसने working directory में actual फ़ाइल को _stage_ नहीं किया।
इसलिए process ने relative string, `data/patientA_rep1_normal_R1_001.fastq.gz`, का उपयोग करने की कोशिश की, लेकिन वह फ़ाइल process working directory के अंदर मौजूद नहीं है।

इन दोनों उदाहरणों को मिलाकर, ये तुम्हें दिखाते हैं कि Nextflow को यह बताना कितना महत्वपूर्ण है कि एक इनपुट को फ़ाइल के रूप में handle किया जाना चाहिए।

!!! note "नोट"

    अगले section में जाने से पहले दोनों जानबूझकर की गई errors को ठीक करना सुनिश्चित करो।

### सारांश

- Path strings बनाम Path objects: Strings सिर्फ text हैं, Path objects smart फ़ाइल references हैं
- `file()` मेथड एक string path को एक Path object में convert करता है जिसके साथ Nextflow काम कर सकता है
- तुम `name`, `simpleName`, `extension`, और `parent` जैसी फ़ाइल properties तक [फ़ाइल विशेषताओं का उपयोग करके](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes) पहुँच सकते हो
- strings के बजाय Path objects का उपयोग करने से Nextflow तुम्हारे वर्कफ़्लो में फ़ाइलों को ठीक से manage कर सकता है
- Process Input Outcomes: उचित फ़ाइल हैंडलिंग के लिए Path objects की आवश्यकता होती है, strings की नहीं, यह सुनिश्चित करने के लिए कि फ़ाइलें processes द्वारा उपयोग के लिए सही ढंग से staged और accessible हों।

---

## 2. Remote फ़ाइलों का उपयोग करना

Nextflow की एक प्रमुख विशेषता local फ़ाइलों (उसी machine पर) से internet पर accessible remote फ़ाइलों के बीच seamlessly switch करने की क्षमता है।

अगर तुम इसे सही तरीके से कर रहे हो, तो तुम्हें अलग-अलग locations से आने वाली फ़ाइलों को accommodate करने के लिए अपने वर्कफ़्लो के logic को कभी बदलने की ज़रूरत नहीं होनी चाहिए।
Remote फ़ाइल का उपयोग करने के लिए तुम्हें बस इतना करना है कि जब तुम इसे वर्कफ़्लो में supply कर रहे हो तो file path में उचित prefix specify करो।

उदाहरण के लिए, `/path/to/data` में कोई prefix नहीं है, जो दर्शाता है कि यह एक 'normal' local file path है, जबकि `s3://path/to/data` में `s3://` prefix शामिल है, जो दर्शाता है कि यह Amazon के S3 object storage में स्थित है।

कई अलग-अलग protocols supported हैं:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

इनमें से किसी का भी उपयोग करने के लिए, बस string में relevant prefix specify करो, जिसे तकनीकी रूप से file path के बजाय Uniform Resource Identifier (URI) कहा जाता है।
Nextflow authentication और फ़ाइलों को सही जगह staging, downloading या uploading और अन्य सभी file operations को handle करेगा जिनकी तुम उम्मीद करोगे।

इस system की मुख्य ताकत यह है कि यह हमें किसी भी pipeline logic को बदले बिना environments के बीच switch करने में सक्षम बनाता है।
उदाहरण के लिए, तुम URI बदलकर remote storage में स्थित full-scale test set पर switch करने से पहले एक छोटे, local test set के साथ develop कर सकते हो।

### 2.1. Internet से एक फ़ाइल का उपयोग करो

चलो इसे test करते हैं, हमारे वर्कफ़्लो को दिए जा रहे local path को एक HTTPS path से बदलकर जो Github में stored उसी डेटा की एक copy की ओर point करता है।

!!! warning "चेतावनी"

    यह केवल तभी काम करेगा जब तुम्हारे पास active internet connection हो।

`main.nf` फिर से खोलो और input path को इस प्रकार बदलो:

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Internet से एक remote फ़ाइल का उपयोग करना
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ print करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // एक string path से Path object बनाएँ
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ print करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

चलो वर्कफ़्लो चलाते हैं:

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

Console आउटपुट में एक अंतर यह है कि path object class अब `nextflow.file.http.XPath` है, जबकि local path के लिए class `sun.nio.fs.UnixPath` था।
तुम्हें इन classes को याद रखने की ज़रूरत नहीं है; हम इसका उल्लेख केवल यह demonstrate करने के लिए करते हैं कि Nextflow विभिन्न locations को identify और handle करता है।

पर्दे के पीछे, Nextflow ने फ़ाइल को work directory के अंदर स्थित एक staging directory में download किया।
उस staged फ़ाइल को फिर एक local फ़ाइल के रूप में treat किया जा सकता है और relevant process directory में symlink किया जा सकता है।

तुम verify कर सकते हो कि process के hash value पर स्थित working directory की सामग्री देखकर ऐसा हुआ।

??? abstract "Work डायरेक्टरी सामग्री"

    अगर process hash `8a/2ab7ca` था, तो तुम work directory explore कर सकते हो:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    Symlink उस remote फ़ाइल की staged copy की ओर point करता है जिसे Nextflow ने automatically download किया।

ध्यान दो कि बड़ी फ़ाइलों के लिए, downloading step local फ़ाइलों पर चलाने की तुलना में कुछ extra समय लेगी।
हालाँकि, Nextflow जाँचता है कि क्या उसके पास पहले से एक staged copy है ताकि अनावश्यक downloads से बचा जा सके।
इसलिए अगर तुम उसी फ़ाइल पर फिर से चलाते हो और staged फ़ाइल को delete नहीं किया है, तो Nextflow staged copy का उपयोग करेगा।

यह दिखाता है कि Nextflow का उपयोग करके local और remote डेटा के बीच switch करना कितना आसान है, जो Nextflow की एक प्रमुख विशेषता है।

!!! note "नोट"

    इस सिद्धांत का एक महत्वपूर्ण अपवाद यह है कि तुम HTTPS के साथ glob patterns या directory paths का उपयोग नहीं कर सकते क्योंकि HTTPS कई फ़ाइलों को list नहीं कर सकता, इसलिए तुम्हें exact file URLs specify करने होंगे।
    हालाँकि, अन्य storage protocols जैसे blob storage (`s3://`, `az://`, `gs://`) दोनों globs और directory paths का उपयोग कर सकते हैं।

    यहाँ बताया गया है कि तुम cloud storage के साथ glob patterns का उपयोग कैसे कर सकते हो:

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // S3 with glob patterns - would match multiple files
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage with glob patterns
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage with glob patterns
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    हम अगले section में व्यवहार में globs के साथ काम करने का तरीका दिखाएंगे।

### 2.2. Local फ़ाइल पर वापस जाओ

हम इस side quest के बाकी हिस्से के लिए अपनी local example फ़ाइलों का उपयोग करने वापस जाएंगे, इसलिए चलो वर्कफ़्लो इनपुट को original फ़ाइल पर वापस switch करते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // एक string path से Path object बनाएँ
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ print करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // एक string path से Path object बनाएँ
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ print करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### सारांश

- Remote डेटा को URI (HTTP, FTP, S3, Azure, Google Cloud) का उपयोग करके access किया जाता है
- Nextflow automatically डेटा को सही जगह download और stage करेगा, जब तक ये paths processes को feed की जा रही हों
- Remote फ़ाइलों को download या upload करने के लिए logic मत लिखो!
- Local और remote फ़ाइलें अलग-अलग object types produce करती हैं लेकिन समान रूप से काम करती हैं
- **महत्वपूर्ण**: HTTP/HTTPS केवल single फ़ाइलों के साथ काम करता है (कोई glob patterns नहीं)
- Cloud storage (S3, Azure, GCS) single फ़ाइलों और glob patterns दोनों को support करता है
- तुम code logic बदले बिना local और remote data sources के बीच seamlessly switch कर सकते हो (जब तक protocol तुम्हारे required operations को support करता हो)

---

## 3. `fromPath()` channel factory का उपयोग करना

अब तक हम एक समय में एक फ़ाइल के साथ काम कर रहे थे, लेकिन Nextflow में, हम आमतौर पर process करने के लिए कई इनपुट फ़ाइलों के साथ एक input channel बनाना चाहेंगे।

ऐसा करने का एक naive तरीका `file()` मेथड को [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) के साथ इस तरह combine करना होगा:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

यह काम करता है, लेकिन यह clunky है।

!!! tip "`file()` बनाम `channel.fromPath()` का उपयोग कब करें"

    - `file()` का उपयोग तब करो जब तुम्हें direct manipulation के लिए एक single Path object चाहिए (जाँचना कि फ़ाइल exist करती है, उसकी attributes पढ़ना, या एक single process invocation को pass करना)
    - `channel.fromPath()` का उपयोग तब करो जब तुम्हें एक channel चाहिए जो कई फ़ाइलें hold कर सके, विशेष रूप से glob patterns के साथ, या जब फ़ाइलें कई processes से flow करेंगी

यहीं पर [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) आता है: एक convenient channel factory जो एक या अधिक static file strings के साथ-साथ glob patterns से channel generate करने के लिए हमें जो functionality चाहिए वह सब bundle करता है।

### 3.1. Channel factory जोड़ो

चलो अपने वर्कफ़्लो को `channel.fromPath` का उपयोग करने के लिए update करते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // फ़ाइल विशेषताएँ print करें
        /* अभी के लिए इन्हें comment out करें, हम बाद में वापस आएंगे!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // फ़ाइल में lines गिनें
        // COUNT_LINES(myFile)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // एक string path से Path object बनाएँ
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // फ़ाइल विशेषताएँ print करें
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // फ़ाइल में lines गिनें
        COUNT_LINES(myFile)
    ```

हमने अभी के लिए attributes print करने वाले code को भी comment out कर दिया है, और इसके बजाय सिर्फ filename print करने के लिए एक `.view` statement जोड़ा है।

वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

जैसा तुम देख सकते हो, file path channel में एक `Path` type object के रूप में load हो रहा है।
यह `file()` ने जो किया होता उसके समान है, सिवाय इसके कि अब हमारे पास एक channel है जिसमें हम चाहें तो और फ़ाइलें load कर सकते हैं।

`channel.fromPath()` का उपयोग करना फ़ाइलों की list से populated एक नया channel बनाने का एक convenient तरीका है।

### 3.2. Channel में फ़ाइलों की विशेषताएँ देखो

Channel factory का उपयोग करने के हमारे पहले प्रयास में, हमने code को simplify किया और सिर्फ file name print किया।

चलो full file attributes print करने पर वापस जाते हैं:

=== "बाद में"

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

        // फ़ाइल में lines गिनें
        COUNT_LINES(ch_files)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // फ़ाइल में lines गिनें
        // COUNT_LINES(ch_files)
    ```

हम `COUNT_LINES` process call को भी re-enable कर रहे हैं यह verify करने के लिए कि file processing हमारे channel-based approach के साथ सही ढंग से काम करती है।

वर्कफ़्लो चलाओ:

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

और यह रहा, पहले जैसे ही results लेकिन अब हमारे पास फ़ाइल एक channel में है, इसलिए हम और जोड़ सकते हैं।

### 3.3. कई फ़ाइलों को match करने के लिए glob का उपयोग करो

Channel में और फ़ाइलें load करने के कई तरीके हैं।
यहाँ हम तुम्हें glob patterns का उपयोग करने का तरीका दिखाएंगे, जो wildcard characters के आधार पर फ़ाइल और directory names को match और retrieve करने का एक convenient तरीका है।
इन patterns को match करने की प्रक्रिया को "globbing" या "filename expansion" कहा जाता है।

!!! note "नोट"

    जैसा कि पहले बताया गया था, Nextflow अधिकांश मामलों में input और output फ़ाइलों को manage करने के लिए globbing को support करता है, HTTPS filepaths को छोड़कर क्योंकि HTTPS कई फ़ाइलों को list नहीं कर सकता।

मान लो हम किसी दिए गए patient, `patientA` से जुड़ी फ़ाइलों की एक pair में दोनों फ़ाइलें retrieve करना चाहते हैं:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

चूँकि filenames के बीच एकमात्र अंतर replicate number है, _यानी_ `R` के बाद की संख्या, हम wildcard character `*` का उपयोग संख्या के स्थान पर इस प्रकार कर सकते हैं:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

यही वह glob pattern है जिसकी हमें ज़रूरत है।

अब हमें बस channel factory में file path को उस glob pattern का उपयोग करने के लिए update करना है:

=== "बाद में"

    ```groovy title="main.nf" linenums="7"
      // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7"
      // channel.fromPath के साथ फ़ाइलें लोड करें
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow automatically पहचान लेगा कि यह एक glob pattern है और इसे उचित रूप से handle करेगा।

इसे test करने के लिए वर्कफ़्लो चलाओ:

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

जैसा तुम देख सकते हो, अब हमारे channel में दो Path objects हैं, जो दिखाता है कि Nextflow ने filename expansion सही ढंग से की है, और दोनों फ़ाइलों को load और process किया है।

इस method का उपयोग करके, हम glob pattern बदलकर जितनी चाहें उतनी या कम फ़ाइलें retrieve कर सकते हैं। अगर हम इसे अधिक generous बनाते, उदाहरण के लिए filenames के सभी variable parts को `*` से replace करके (_जैसे_ `data/patient*_rep*_*_R*_001.fastq.gz`) तो हम `data` directory में सभी example फ़ाइलें grab कर सकते थे।

### सारांश

- `channel.fromPath()` एक pattern से match होने वाली फ़ाइलों के साथ एक channel बनाता है
- प्रत्येक फ़ाइल channel में एक अलग element के रूप में emit होती है
- हम कई फ़ाइलों को match करने के लिए glob pattern का उपयोग कर सकते हैं
- फ़ाइलें automatically full attributes के साथ Path objects में convert हो जाती हैं
- `.view()` method channel contents की inspection की अनुमति देता है

---

## 4. Filenames से बुनियादी मेटाडेटा निकालना

अधिकांश वैज्ञानिक domains में, उन फ़ाइलों के नामों में मेटाडेटा encode होना बहुत सामान्य है जिनमें डेटा होता है।
उदाहरण के लिए, bioinformatics में, sequencing डेटा वाली फ़ाइलों को अक्सर इस तरह नाम दिया जाता है जो sample, condition, replicate, और read number के बारे में जानकारी encode करता है।

अगर filenames एक consistent convention के अनुसार बनाए गए हैं, तो तुम उस मेटाडेटा को एक standardized तरीके से निकाल सकते हो और अपने analysis के दौरान इसका उपयोग कर सकते हो।
यह एक बड़ा 'if' है, निश्चित रूप से, और जब भी तुम filename structure पर rely करते हो तो तुम्हें बहुत सावधान रहना चाहिए; लेकिन वास्तविकता यह है कि यह approach बहुत व्यापक रूप से उपयोग की जाती है, इसलिए चलो देखते हैं कि Nextflow में यह कैसे किया जाता है।

हमारे example data के मामले में, हम जानते हैं कि filenames में consistently structured मेटाडेटा शामिल है।
उदाहरण के लिए, filename `patientA_rep1_normal_R2_001` निम्नलिखित encode करता है:

- patient ID: `patientA`
- replicate ID: `rep1`
- sample type: `normal` (`tumor` के विपरीत)
- read set: `R1` (`R2` के विपरीत)

हम अपने वर्कफ़्लो को तीन steps में इस जानकारी को retrieve करने के लिए modify करेंगे:

1. फ़ाइल का `simpleName` retrieve करो, जिसमें मेटाडेटा शामिल है
2. `tokenize()` नामक method का उपयोग करके मेटाडेटा को अलग करो
3. मेटाडेटा को organize करने के लिए map का उपयोग करो

!!! warning "चेतावनी"

    तुम्हें filenames में कभी भी sensitive जानकारी encode नहीं करनी चाहिए, जैसे patient names या अन्य identifying characteristics, क्योंकि यह patient privacy या अन्य relevant security restrictions को compromise कर सकता है।

### 4.1. `simpleName` retrieve करो

`simpleName` एक file attribute है जो path और extension से stripped filename से मेल खाता है।

वर्कफ़्लो में निम्नलिखित बदलाव करो:

=== "बाद में"

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

यह `simpleName` को retrieve करता है और इसे `map()` operation का उपयोग करके full file object के साथ associate करता है।

इसे test करने के लिए वर्कफ़्लो चलाओ:

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

Channel में प्रत्येक element अब `simpleName` और original file object वाला एक tuple है।

### 4.2. `simplename` से मेटाडेटा निकालो

इस point पर, हम जो मेटाडेटा चाहते हैं वह `simplename` में embedded है, लेकिन हम individual items को directly access नहीं कर सकते।
इसलिए हमें `simplename` को उसके components में split करना होगा।
सौभाग्य से, वे components original filename में underscores द्वारा simply separated हैं, इसलिए हम `tokenize()` नामक एक common Nextflow method apply कर सकते हैं जो इस कार्य के लिए perfect है।

वर्कफ़्लो में निम्नलिखित बदलाव करो:

=== "बाद में"

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

`tokenize()` method `simpleName` string को जहाँ भी underscores मिलेंगे वहाँ split करेगा, और substrings वाली एक list return करेगा।

वर्कफ़्लो चलाओ:

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

अब हमारे channel में प्रत्येक element के लिए tuple में मेटाडेटा की list (_जैसे_ `[patientA, rep1, normal, R1, 001]`) और original file object शामिल है।

बहुत अच्छा!
हमने अपनी patient information को एक single string से strings की list में break down किया है।
अब हम patient information के प्रत्येक भाग को अलग से handle कर सकते हैं।

### 4.3. मेटाडेटा को organize करने के लिए map का उपयोग करो

हमारा मेटाडेटा अभी सिर्फ एक flat list है।
इसका उपयोग करना काफी आसान है लेकिन पढ़ना मुश्किल है।

```console
[patientA, rep1, normal, R1, 001]
```

Index 3 पर item क्या है? क्या तुम मेटाडेटा structure की original explanation पर वापस जाए बिना बता सकते हो?

यह एक key-value store का उपयोग करने का एक अच्छा अवसर है, जहाँ प्रत्येक item में keys का एक set और उनके associated values होते हैं, इसलिए तुम corresponding value प्राप्त करने के लिए आसानी से प्रत्येक key को refer कर सकते हो।

हमारे उदाहरण में, इसका मतलब है इस organization से जाना:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

इस पर:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

Nextflow में, इसे [map](https://nextflow.io/docs/latest/script.html#maps) कहा जाता है।

चलो अब अपनी flat list को एक map में convert करते हैं।
वर्कफ़्लो में निम्नलिखित बदलाव करो:

=== "बाद में"

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

यहाँ मुख्य बदलाव हैं:

- **Destructuring assignment**: `def (patient, replicate, type, readNum) = ...` tokenized values को एक line में named variables में extract करता है
- **Map literal syntax**: `[id: patient, replicate: ...]` एक map बनाता है जहाँ प्रत्येक key (जैसे `id`) एक value (जैसे `patient`) के साथ associated है
- **Nested structure**: बाहरी list `[..., myFile]` metadata map को original file object के साथ pair करती है

हमने `replace()` नामक एक string replacement method का उपयोग करके कुछ metadata strings को भी simplify किया ताकि कुछ unnecessary characters को remove किया जा सके (_जैसे_ `replicate.replace('rep', '')` replicate IDs से केवल number रखने के लिए)।

चलो वर्कफ़्लो फिर से चलाते हैं:

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

अब मेटाडेटा clearly labeled है (_जैसे_ `[id:patientA, replicate:1, type:normal, readNum:2]`) इसलिए यह बताना बहुत आसान है कि क्या है।

वर्कफ़्लो में मेटाडेटा के elements का वास्तव में उपयोग करना भी बहुत आसान होगा, और हमारे code को पढ़ना और maintain करना आसान बनाएगा।

### सारांश

- हम Nextflow में filenames को एक full programming language की शक्ति के साथ handle कर सकते हैं
- हम relevant information निकालने के लिए filenames को strings के रूप में treat कर सकते हैं
- `tokenize()` और `replace()` जैसे methods का उपयोग हमें filename में strings को manipulate करने की अनुमति देता है
- `.map()` operation structure को preserve करते हुए channel elements को transform करता है
- Structured मेटाडेटा (maps) code को positional lists की तुलना में अधिक readable और maintainable बनाता है

आगे, हम paired data files को handle करने का तरीका देखेंगे।

---

## 5. Paired data files को handle करना

कई experimental designs paired data files produce करते हैं जिन्हें explicitly paired तरीके से handle करने से फायदा होता है।
उदाहरण के लिए, bioinformatics में, sequencing डेटा अक्सर paired reads के रूप में generate होता है, यानी sequence strings जो DNA के एक ही fragment से originate होती हैं (अक्सर 'forward' और 'reverse' कहलाती हैं क्योंकि वे opposite ends से read की जाती हैं)।

यही हमारे example data का मामला है, जहाँ R1 और R2 reads के दो sets को refer करते हैं।

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow इस तरह की paired files के साथ काम करने के लिए `channel.fromFilePairs()` नामक एक specialized channel factory प्रदान करता है, जो automatically shared naming pattern के आधार पर फ़ाइलों को group करता है। यह तुम्हें कम effort के साथ paired files को अधिक tightly associate करने की अनुमति देता है।

हम अपने वर्कफ़्लो को इसका फायदा उठाने के लिए modify करेंगे।
इसमें दो steps लगेंगे:

1. Channel factory को `channel.fromFilePairs()` पर switch करो
2. मेटाडेटा extract और map करो

### 5.1. Channel factory को `channel.fromFilePairs()` पर switch करो

`channel.fromFilePairs` का उपयोग करने के लिए, हमें वह pattern specify करना होगा जिसका उपयोग Nextflow को एक pair में दो members की पहचान करने के लिए करना चाहिए।

हमारे example data पर वापस जाते हुए, हम naming pattern को इस प्रकार formalize कर सकते हैं:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

यह हमारे पहले उपयोग किए गए glob pattern के समान है, सिवाय इसके कि यह specifically उन substrings को enumerate करता है (R के ठीक बाद आने वाला `1` या `2`) जो pair के दो members की पहचान करते हैं।

चलो वर्कफ़्लो `main.nf` को accordingly update करते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // channel.fromFilePairs के साथ फ़ाइलें लोड करें
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* अभी के लिए mapping को comment out करें, हम बाद में वापस आएंगे!
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

हमने channel factory switch किया और file matching pattern को adapt किया, और साथ ही map operation को comment out किया।
हम इसे बाद में कुछ modifications के साथ वापस जोड़ेंगे।

इसे test करने के लिए वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? failure "कमांड आउटपुट"

    ```console hl_lines="7-8"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_koch] DSL2 - revision: 44fdf66105

    [-        ] COUNT_LINES -
    [-        ] COUNT_LINES -
    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ERROR ~ Error executing process > 'COUNT_LINES (1)'

    Caused by:
      Not a valid path value: 'patientA_rep1_normal_R'



    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

अरे नहीं, इस बार run fail हो गया!

Error message का relevant हिस्सा यहाँ है:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

यह इसलिए है क्योंकि हमने channel factory बदल दी है।
अब तक, original input channel में केवल file paths थे।
हम जो सभी metadata manipulation कर रहे थे उसने वास्तव में channel contents को affect नहीं किया।

अब जब हम `.fromFilePairs` channel factory का उपयोग कर रहे हैं, resulting channel की contents अलग हैं।
हम केवल एक channel element देखते हैं, जो दो items वाले tuple से composed है: दोनों फ़ाइलों द्वारा shared `simpleName` का हिस्सा, जो identifier के रूप में serve करता है, और दो file objects वाला एक tuple, format `id, [ file1, file2 ]` में।

यह बहुत अच्छा है, क्योंकि Nextflow ने shared prefix की जाँच करके और इसे patient identifier के रूप में उपयोग करके patient name extract करने का कठिन काम किया है।

हालाँकि, यह हमारे current वर्कफ़्लो को break करता है।
अगर हम process को बदले बिना `COUNT_LINES` को उसी तरह चलाना चाहते, तो हमें file paths extract करने के लिए एक mapping operation apply करनी होगी।
लेकिन हम ऐसा नहीं करेंगे, क्योंकि हमारा ultimate goal एक अलग process, `ANALYZE_READS`, का उपयोग करना है जो file pairs को उचित रूप से handle करता है।

इसलिए चलो बस `COUNT_LINES` को call करना comment out (या delete) करते हैं और आगे बढ़ते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // फ़ाइल में lines गिनें
        // COUNT_LINES(ch_files)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // फ़ाइल में lines गिनें
        COUNT_LINES(ch_files)
    ```

तुम `COUNT_LINES` include statement को भी comment out या delete कर सकते हो, लेकिन इसका कोई functional effect नहीं होगा।

अब चलो वर्कफ़्लो फिर से चलाते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

बढ़िया, इस बार वर्कफ़्लो सफल हुआ!

हालाँकि, हमें अभी भी `id` field से बाकी मेटाडेटा निकालना है।

### 5.2. File pairs से मेटाडेटा extract और organize करो

हमारी पहले की `map` operation काम नहीं करेगी क्योंकि यह data structure से match नहीं करती, लेकिन हम इसे काम करने के लिए modify कर सकते हैं।

हमारे पास पहले से ही उस string में actual patient identifier तक पहुँच है जिसे `fromFilePairs()` ने identifier के रूप में उपयोग किया, इसलिए हम पहले की तरह Path object से `simpleName` प्राप्त किए बिना मेटाडेटा extract करने के लिए इसका उपयोग कर सकते हैं।

वर्कफ़्लो में map operation को uncomment करो और निम्नलिखित बदलाव करो:

=== "बाद में"

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
        /* अभी के लिए mapping को comment out करें, हम बाद में वापस आएंगे!
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

इस बार map `id, files` से शुरू होता है न कि सिर्फ `myFile` से, और `tokenize()` को `myFile.simpleName` के बजाय `id` पर apply किया जाता है।

यह भी ध्यान दो कि हमने `tokenize()` line से `readNum` को drop किया है; कोई भी substrings जिन्हें हम specifically name नहीं करते (बाईं ओर से शुरू करके) silently drop हो जाएंगे।
हम यह कर सकते हैं क्योंकि paired files अब tightly associated हैं, इसलिए हमें metadata map में `readNum` की अब ज़रूरत नहीं है।

चलो वर्कफ़्लो चलाते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

और यह रहा: हमारे पास output tuple की पहली position में metadata map (`[id:patientA, replicate:1, type:normal]`) है, उसके बाद paired files का tuple है, जैसा intended था।

बेशक, यह केवल उस specific pair of files को pick up और process करेगा।
अगर तुम कई pairs को process करने के साथ experiment करना चाहते हो, तो input pattern में wildcards जोड़ने की कोशिश करो और देखो क्या होता है।
उदाहरण के लिए, `data/patientA_rep1_*_R{1,2}_001.fastq.gz` try करो।

### सारांश

- [`channel.fromFilePairs()` automatically related files को ढूंढता और pair करता है](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- यह तुम्हारी pipeline में paired-end reads को handle करना simplify करता है
- Paired files को `[id, [file1, file2]]` tuples के रूप में group किया जा सकता है
- मेटाडेटा extraction individual files के बजाय paired file ID से की जा सकती है

---

## 6. Processes में file operations का उपयोग करना

अब चलो इन सबको एक simple process में एक साथ रखते हैं ताकि Nextflow process के अंदर file operations का उपयोग करने के तरीके को reinforce किया जा सके।

हम तुम्हें `ANALYZE_READS` नामक एक pre-written process module प्रदान करते हैं जो metadata का एक tuple और input files की एक pair लेता है और उनका analysis करता है।
हम imagine कर सकते हैं कि यह sequence alignment, या variant calling या इस data type के लिए कोई अन्य step कर रहा है।

चलो शुरू करते हैं।

### 6.1. Process को import करो और code की जाँच करो

इस process को वर्कफ़्लो में उपयोग करने के लिए, हमें बस workflow block से पहले एक module include statement जोड़नी होगी।

वर्कफ़्लो में निम्नलिखित बदलाव करो:

=== "बाद में"

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

तुम module फ़ाइल खोलकर उसके code की जाँच कर सकते हो:

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

    `tag` और `publishDir` directives string interpolation (`"${...}"`) के बजाय closure syntax (`{ ... }`) का उपयोग करते हैं।
    यह इसलिए है क्योंकि ये directives input variables (`meta`) को reference करते हैं जो runtime तक available नहीं होते।
    Closure syntax evaluation को तब तक defer करता है जब तक process actually run नहीं होता।

!!! note "नोट"

    हम convention के अनुसार अपने metadata map को `meta` कह रहे हैं।
    Meta maps में deeper dive के लिए, [Metadata and meta maps](../metadata/) side quest देखो।

### 6.2. Workflow में process को call करो

अब जब process वर्कफ़्लो के लिए उपलब्ध है, हम इसे चलाने के लिए `ANALYZE_READS` process को call कर सकते हैं।

हमारे example data पर इसे चलाने के लिए, हमें दो काम करने होंगे:

1. Remapped channel को एक नाम दो
2. Process को call करो

#### 6.2.1. Remapped input channel को नाम दो

हमने पहले mapping manipulations को directly input channel पर apply किया था।
`ANALYZE_READS` process को remapped contents feed करने के लिए (और इसे clear और easy to read तरीके से करने के लिए) हम `ch_samples` नामक एक नया channel बनाना चाहते हैं।

हम [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) operator का उपयोग करके ऐसा कर सकते हैं।

Main workflow में, `.view()` operator को `.set { ch_samples }` से replace करो, और एक line जोड़ो जो test करे कि हम channel को नाम से refer कर सकते हैं।

=== "बाद में"

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

        // अस्थायी: ch_samples में झाँकें
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

चलो यह चलाते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

यह confirm करता है कि अब हम channel को नाम से refer कर सकते हैं।

#### 6.2.2. Data पर process को call करो

अब चलो `ch_samples` channel पर `ANALYZE_READS` process को actually call करते हैं।

Main workflow में, निम्नलिखित code changes करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="23"
        // Analysis चलाएँ
        ANALYZE_READS(ch_samples)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="23"
        // अस्थायी: ch_samples में झाँकें
        ch_samples.view()
    ```

चलो यह चलाते हैं:

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

यह process अपने outputs को एक `results` directory में publish करने के लिए set है, इसलिए वहाँ देखो।

??? abstract "डायरेक्टरी और फ़ाइल सामग्री"

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

Process ने हमारे inputs लिए और designed के अनुसार patient metadata वाली एक नई फ़ाइल बनाई।
शानदार!

### 6.3. और अधिक patients शामिल करो

बेशक, यह एक single patient के लिए files की एक single pair को process कर रहा है, जो Nextflow के साथ तुम्हें मिलने वाले high throughput का प्रकार नहीं है।
तुम शायद एक साथ बहुत अधिक डेटा process करना चाहोगे।

याद करो `channel.fromPath()` input के रूप में एक _glob_ accept करता है, जिसका मतलब है कि यह pattern से match होने वाली किसी भी संख्या में फ़ाइलें accept कर सकता है।
इसलिए अगर हम सभी patients को include करना चाहते हैं, तो हम बस input string को modify कर सकते हैं ताकि अधिक patients शामिल हों, जैसा कि पहले passing में बताया गया था।

चलो pretend करते हैं कि हम जितना possible हो उतना greedy होना चाहते हैं।
वर्कफ़्लो में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // channel.fromFilePairs के साथ फ़ाइलें लोड करें
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // channel.fromFilePairs के साथ फ़ाइलें लोड करें
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Pipeline फिर से चलाओ:

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

Results directory में अब सभी available data के results होने चाहिए।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

सफलता! हमने एक ही बार में सभी patients का analysis किया! है ना?

शायद नहीं।
अगर तुम अधिक ध्यान से देखो, तो हमारे पास एक समस्या है: patientA के लिए हमारे पास दो replicates हैं, लेकिन केवल एक output फ़ाइल है!
हम हर बार output फ़ाइल को overwrite कर रहे हैं।

### 6.4. Published files को unique बनाओ

चूँकि हमारे पास patient metadata तक पहुँच है, हम इसका उपयोग published files को unique बनाने के लिए कर सकते हैं, या तो directory structure में या filenames में differentiating metadata शामिल करके।

Workflow में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "पहले"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

यहाँ हम sample types और replicates के लिए account करने के लिए additional directory levels का उपयोग करने का option दिखाते हैं, लेकिन तुम filename level पर भी ऐसा करने के साथ experiment कर सकते हो।

अब pipeline एक और बार चलाओ, लेकिन पहले results directory को remove करना सुनिश्चित करो ताकि तुम्हें एक clean workspace मिले:

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

अब results directory जाँचो:

??? abstract "डायरेक्टरी सामग्री"

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

और यह रहा, हमारा सारा मेटाडेटा, नeatly organized। यह सफलता है!

एक बार जब तुम्हारा मेटाडेटा इस तरह एक map में load हो जाता है तो तुम बहुत कुछ और कर सकते हो:

1. Patient attributes के आधार पर organized output directories बनाओ
2. Patient properties के आधार पर processes में decisions लो
3. Metadata values के आधार पर data को split, join, और recombine करो

मेटाडेटा को explicit और data के साथ attached रखने का यह pattern (filenames में encode करने के बजाय) Nextflow में एक core best practice है जो robust, maintainable analysis workflows बनाने में सक्षम बनाता है।
तुम इसके बारे में [Metadata and meta maps](../metadata/) side quest में अधिक जान सकते हो।

### सारांश

- `publishDir` directive metadata values के आधार पर outputs को organize कर सकता है
- Tuples में metadata results के structured organization को enable करता है
- यह approach clear data provenance के साथ maintainable workflows बनाता है
- Processes metadata और files के tuples को input के रूप में ले सकते हैं
- `tag` directive execution logs में process identification प्रदान करता है
- Workflow structure channel creation को process execution से अलग करती है

---

## सारांश

इस side quest में, तुमने Nextflow में फ़ाइलों के साथ काम करना सीखा, बुनियादी operations से लेकर files के collections को handle करने की अधिक advanced techniques तक।

अपने काम में इन techniques को apply करने से तुम्हें अधिक efficient और maintainable workflows बनाने में मदद मिलेगी, विशेष रूप से जब complex naming conventions वाली बड़ी संख्या में फ़ाइलों के साथ काम करते हो।

### मुख्य patterns

1.  **बुनियादी File Operations:** हमने `file()` के साथ Path objects बनाए और name, extension, और parent directory जैसी file attributes तक पहुँचे, strings और Path objects के बीच का अंतर सीखा।

    - `file()` के साथ Path object बनाओ

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - File attributes प्राप्त करो

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Remote Files का उपयोग करना**: हमने सीखा कि URIs का उपयोग करके local और remote files के बीच transparently कैसे switch करें, Nextflow की workflow logic बदले बिना विभिन्न sources से files handle करने की क्षमता demonstrate करते हुए।

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

3.  **`fromPath()` channel factory का उपयोग करके फ़ाइलें लोड करना:** हमने `channel.fromPath()` के साथ file patterns से channels बनाए और उनकी file attributes देखीं, object types सहित।

    - File pattern से channel बनाओ

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - File attributes प्राप्त करो

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Filenames से Patient Metadata निकालना:** हमने filenames से metadata extract और structure करने के लिए `tokenize()` और `replace()` का उपयोग किया, उन्हें organized maps में convert किया।

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **`channel.fromFilePairs` के साथ Simplify करना:** हमने `channel.fromFilePairs()` का उपयोग करके automatically related files को pair किया और paired file IDs से metadata extract किया।

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Processes में File Operations का उपयोग करना:** हमने proper input handling के साथ Nextflow processes में file operations को integrate किया, metadata के आधार पर outputs organize करने के लिए `publishDir` का उपयोग किया।

    - Process inputs के साथ meta map associate करो

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

    - Metadata के आधार पर outputs organize करो

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### अतिरिक्त संसाधन

- [Nextflow Documentation: Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## आगे क्या है?

[Side Quests के menu](../) पर वापस जाओ या list में अगले topic पर जाने के लिए page के नीचे दाईं ओर button click करो।
