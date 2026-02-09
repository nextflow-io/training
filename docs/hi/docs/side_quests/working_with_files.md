# फ़ाइल इनपुट प्रोसेसिंग

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [और जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

वैज्ञानिक विश्लेषण workflows में अक्सर बड़ी संख्या में फ़ाइलों को प्रोसेस करना शामिल होता है।
Nextflow फ़ाइलों को कुशलता से संभालने के लिए शक्तिशाली टूल प्रदान करता है, जो तुम्हें न्यूनतम कोड के साथ अपने डेटा को व्यवस्थित और प्रोसेस करने में मदद करता है।

### सीखने के लक्ष्य

इस side quest में, हम यह जानेंगे कि Nextflow फ़ाइलों को कैसे संभालता है, बुनियादी फ़ाइल ऑपरेशन से लेकर फ़ाइल संग्रह के साथ काम करने की अधिक उन्नत तकनीकों तक।
तुम सीखोगे कि फ़ाइल नामों से metadata कैसे निकाला जाए, जो वैज्ञानिक विश्लेषण pipelines में एक सामान्य आवश्यकता है।

इस side quest के अंत तक, तुम सक्षम हो जाओगे:

- Nextflow की `file()` method का उपयोग करके फ़ाइल path strings से Path objects बनाना
- फ़ाइल attributes जैसे name, extension, और parent directory तक पहुँचना
- URIs का उपयोग करके local और remote फ़ाइलों को पारदर्शी रूप से संभालना
- `channel.fromPath()` और `channel.fromFilePairs()` के साथ फ़ाइल हैंडलिंग को स्वचालित करने के लिए channels का उपयोग करना
- string manipulation का उपयोग करके फ़ाइल नामों से metadata निकालना और संरचित करना
- pattern matching और glob expressions का उपयोग करके संबंधित फ़ाइलों को group करना
- उचित input हैंडलिंग के साथ Nextflow processes में फ़ाइल ऑपरेशन को integrate करना
- metadata-driven directory structures का उपयोग करके process outputs को व्यवस्थित करना

ये कौशल तुम्हें ऐसे workflows बनाने में मदद करेंगे जो विभिन्न प्रकार के फ़ाइल inputs को बेहतरीन लचीलेपन के साथ संभाल सकते हैं।

### पूर्वापेक्षाएँ

इस side quest को शुरू करने से पहले, तुम्हें चाहिए:

- [Hello Nextflow](../../hello_nextflow/) tutorial या समकक्ष beginner's course पूरा किया हो।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (processes, channels, operators) का उपयोग करने में सहज हो

---

## 0. शुरू करना

#### training codespace खोलें

यदि तुमने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित अनुसार training environment खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### project directory में जाओ

चलो उस directory में चलते हैं जहाँ इस tutorial के लिए फ़ाइलें स्थित हैं।

```bash
cd side-quests/working_with_files
```

तुम VSCode को इस directory पर focus करने के लिए सेट कर सकते हो:

```bash
code .
```

#### सामग्री की समीक्षा करो

तुम्हें `main.nf` नामक एक सरल workflow फ़ाइल, दो module फ़ाइलों वाली एक `modules` directory, और कुछ उदाहरण डेटा फ़ाइलों वाली एक `data` directory मिलेगी।

??? abstract "Directory contents"

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

इस directory में तीन patients (A, B, C) से paired-end sequencing डेटा है।

प्रत्येक patient के लिए, हमारे पास `tumor` (आमतौर पर tumor biopsies से उत्पन्न) या `normal` (स्वस्थ ऊतक या रक्त से लिए गए) प्रकार के नमूने हैं।
यदि तुम cancer विश्लेषण से परिचित नहीं हो, तो बस इतना जान लो कि यह एक प्रयोगात्मक मॉडल से मेल खाता है जो contrastive विश्लेषण करने के लिए paired tumor/normal नमूनों का उपयोग करता है।

विशेष रूप से patient A के लिए, हमारे पास दो सेट technical replicates (repeats) हैं।

sequencing डेटा फ़ाइलों का नाम 'forward reads' और 'reverse reads' के रूप में जाने जाने वाले के लिए एक विशिष्ट `_R1_` और `_R2_` convention के साथ रखा गया है।

_यदि तुम इस प्रयोगात्मक डिज़ाइन से परिचित नहीं हो तो चिंता मत करो, यह इस tutorial को समझने के लिए महत्वपूर्ण नहीं है।_

#### assignment की समीक्षा करो

तुम्हारी चुनौती एक Nextflow workflow लिखना है जो:

1. **Load** Nextflow की फ़ाइल हैंडलिंग methods का उपयोग करके input फ़ाइलें
2. **Extract** फ़ाइल नाम संरचना से metadata (patient ID, replicate, sample type)
3. **Group** `channel.fromFilePairs()` का उपयोग करके paired फ़ाइलें (R1/R2) एक साथ
4. **Process** प्रदान किए गए विश्लेषण module के साथ फ़ाइलें
5. **Organize** निकाले गए metadata के आधार पर एक directory संरचना में outputs

#### तैयारी checklist

क्या तुम dive in करने के लिए तैयार हो?

- [ ] मैं इस course के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता हूँ
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी working directory उचित रूप से सेट कर ली है
- [ ] मैं assignment को समझता हूँ

यदि तुम सभी boxes को check कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. बुनियादी फ़ाइल ऑपरेशन

### 1.1. `.class` के साथ किसी object के प्रकार की पहचान करो

workflow फ़ाइल `main.nf` पर एक नज़र डालो:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Create a Path object from a string path
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

यह एक mini-workflow है (बिना किसी processes के) जो अपने workflow में एक single फ़ाइल path को refer करता है, फिर इसे console पर print करता है, इसकी class के साथ।

??? info "`.class` क्या है?"

    Nextflow में, `.class` हमें बताता है कि हम किस प्रकार के object के साथ काम कर रहे हैं। यह पूछने जैसा है "यह किस प्रकार की चीज़ है?" यह पता लगाने के लिए कि यह एक string है, एक number है, एक फ़ाइल है, या कुछ और।
    यह हमें अगले sections में एक plain string और एक Path object के बीच अंतर को स्पष्ट करने में मदद करेगा।

चलो workflow चलाते हैं:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

जैसा कि तुम देख सकते हो, Nextflow ने string path को बिल्कुल वैसे ही print किया जैसे हमने इसे लिखा था।

यह सिर्फ text output है; Nextflow ने अभी तक इसके साथ कुछ खास नहीं किया है।
हमने यह भी confirm किया है कि जहाँ तक Nextflow का संबंध है, यह केवल एक string है (class `java.lang.String` की)।
यह समझ में आता है, क्योंकि हमने अभी तक Nextflow को नहीं बताया है कि यह एक फ़ाइल से मेल खाता है।

### 1.2. file() के साथ एक Path object बनाओ

हम Nextflow को बता सकते हैं कि path strings से [Path objects](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) बनाकर फ़ाइलों को कैसे संभालना है।

हमारे workflow में, हम `file()` method का उपयोग करके string path `data/patientA_rep1_normal_R1_001.fastq.gz` को एक Path object में convert कर सकते हैं, जो फ़ाइल properties और operations तक पहुँच प्रदान करता है।

`main.nf` को निम्नानुसार `file()` के साथ string को wrap करने के लिए edit करो:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

अब workflow को फिर से चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

इस बार, तुम हमारे द्वारा input के रूप में प्रदान किए गए relative path के बजाय पूर्ण absolute path देखते हो।

Nextflow ने हमारी string को एक Path object में convert कर दिया है और इसे system पर वास्तविक फ़ाइल स्थान पर resolve कर दिया है।
फ़ाइल path अब absolute होगा, जैसे `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`।

यह भी ध्यान दो कि Path object class `sun.nio.fs.UnixPath` है: यह local फ़ाइलों को represent करने का Nextflow का तरीका है।
जैसा कि हम बाद में देखेंगे, remote फ़ाइलों के अलग-अलग class नाम होंगे (जैसे HTTP फ़ाइलों के लिए `nextflow.file.http.XPath`), लेकिन वे सभी बिल्कुल उसी तरह काम करती हैं और तुम्हारे workflows में समान रूप से उपयोग की जा सकती हैं।

!!! tip

    **मुख्य अंतर:**

    - **Path string**: सिर्फ text जिसे Nextflow characters के रूप में treat करता है
    - **Path object**: एक smart फ़ाइल reference जिसके साथ Nextflow काम कर सकता है

    इसे इस तरह सोचो: एक path string कागज पर एक address लिखने जैसा है, जबकि एक Path object एक GPS device में लोड किए गए address की तरह है जो वहाँ navigate करना जानता है और तुम्हें यात्रा के बारे में विवरण बता सकता है।

### 1.3. फ़ाइल attributes तक पहुँचो

यह क्यों मददगार है? खैर, अब जब Nextflow समझता है कि `myFile` एक Path object है न कि सिर्फ एक string, हम Path object के विभिन्न attributes तक पहुँच सकते हैं।

चलो अपने workflow को built-in फ़ाइल attributes को print करने के लिए update करते हैं:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

तुम ऊपर console पर print किए गए विभिन्न फ़ाइल attributes देखते हो।

### 1.4. फ़ाइल को एक process में feed करो

strings और Path objects के बीच का अंतर तब महत्वपूर्ण हो जाता है जब तुम processes के साथ वास्तविक workflows बनाना शुरू करते हो।
अब तक हमने verify किया है कि Nextflow अब हमारी input फ़ाइल को एक फ़ाइल के रूप में treat कर रहा है, लेकिन चलो देखते हैं कि क्या हम वास्तव में उस फ़ाइल पर एक process में कुछ चला सकते हैं।

#### 1.4.1. process को import करो और code की जाँच करो

हम तुम्हें `COUNT_LINES` नामक एक pre-written process module प्रदान करते हैं जो एक फ़ाइल input लेता है और गिनता है कि इसमें कितनी lines हैं।

workflow में process का उपयोग करने के लिए, तुम्हें बस workflow block से पहले एक include statement जोड़ना होगा:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

तुम इसके code की जाँच करने के लिए module फ़ाइल खोल सकते हो:

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

जैसा कि तुम देख सकते हो, यह एक काफी सीधी छोटी script है जो फ़ाइल को unzip करती है और गिनती करती है कि इसमें कितनी lines हैं।

??? info "`debug true` क्या करता है?"

    process definition में `debug true` directive Nextflow को तुम्हारी script से output (जैसे line count "40") को सीधे execution log में print करने का कारण बनता है।
    इसके बिना, तुम केवल process execution status देखोगे लेकिन तुम्हारी script से वास्तविक output नहीं।

    Nextflow processes को debug करने के बारे में अधिक जानकारी के लिए, [Debugging Nextflow Workflows](debugging.md) side quest देखो।

#### 1.4.2. `COUNT_LINES` में एक call जोड़ो

अब जब process workflow के लिए उपलब्ध है, हम input फ़ाइल पर इसे चलाने के लिए `COUNT_LINES` process में एक call जोड़ सकते हैं।

workflow में निम्नलिखित edits करो:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
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

??? success "Command output"

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

यह दिखाता है कि हम एक process के अंदर फ़ाइल पर उचित रूप से operate करने में सक्षम हैं।

विशेष रूप से, Nextflow ने निम्नलिखित ऑपरेशन सफलतापूर्वक किए:

- फ़ाइल को working directory में staged किया
- .gz फ़ाइल को decompress किया
- lines को गिना (इस मामले में 40 lines)
- बिना error के पूरा किया

इस smooth operation की कुंजी यह है कि हम स्पष्ट रूप से Nextflow को बता रहे हैं कि हमारा input एक फ़ाइल है और इसे इस तरह treat किया जाना चाहिए।

### 1.5. बुनियादी फ़ाइल input errors को troubleshoot करो

यह अक्सर Nextflow में नए लोगों को परेशान करता है, तो चलो कुछ मिनट लेते हैं यह देखने के लिए कि जब तुम इसे गलत करते हो तो क्या होता है।

दो मुख्य स्थान हैं जहाँ तुम फ़ाइल हैंडलिंग को गलत कर सकते हो: workflow के स्तर पर, और process के स्तर पर।

#### 1.5.1. Workflow-level error

चलो देखते हैं कि क्या होता है यदि हम workflow block में input specify करते समय फ़ाइल को एक string के रूप में treat करने के लिए revert करते हैं।

workflow में निम्नलिखित edits करो, path-specific print statements को comment out करना सुनिश्चित करो:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Print file attributes
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

और अब workflow चलाओ:

```bash
nextflow run main.nf
```

??? failure "Command output"

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

यह महत्वपूर्ण bit है:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

जब तुम एक `path` input specify करते हो, Nextflow validate करता है कि तुम वास्तविक फ़ाइल references pass कर रहे हो, सिर्फ strings नहीं।
यह error तुम्हें बता रहा है कि `'data/patientA_rep1_normal_R1_001.fastq.gz'` एक valid path value नहीं है क्योंकि यह एक string है, Path object नहीं।

Nextflow ने तुरंत समस्या का पता लगाया और process शुरू करने से पहले ही रुक गया।

#### 1.5.2. Process-level error

दूसरी जगह जहाँ हम यह specify करना भूल सकते हैं कि हम चाहते हैं कि Nextflow input को एक फ़ाइल के रूप में treat करे, वह process definition में है।

!!! warning "1.5.1 से workflow error रखो"

    इस test के सही तरीके से काम करने के लिए, workflow को इसकी broken state में रखो (`file()` के बजाय एक plain string का उपयोग करते हुए)।
    जब process में `val` के साथ combined किया जाता है, तो यह नीचे दिखाई गई error उत्पन्न करता है।

module में निम्नलिखित edit करो:

=== "After"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "Before"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

और अब workflow को फिर से चलाओ:

```bash
nextflow run main.nf
```

??? failure "Command output"

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

यह error के बारे में बहुत सारे विवरण दिखाता है क्योंकि process debugging जानकारी output करने के लिए set है, जैसा कि ऊपर noted है।

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

यह कहता है कि system फ़ाइल नहीं ढूंढ सका; हालाँकि यदि तुम path को देखो, तो उस स्थान पर उस नाम की एक फ़ाइल है।

जब हमने इसे चलाया, तो Nextflow ने string value को script में pass कर दिया, लेकिन इसने वास्तविक फ़ाइल को working directory में _stage_ नहीं किया।
तो process ने relative string, `data/patientA_rep1_normal_R1_001.fastq.gz`, का उपयोग करने की कोशिश की, लेकिन वह फ़ाइल process working directory के भीतर मौजूद नहीं है।

एक साथ लिए गए, ये दो उदाहरण तुम्हें दिखाते हैं कि यदि कोई input एक फ़ाइल के रूप में संभाला जाना चाहिए तो Nextflow को बताना कितना महत्वपूर्ण है।

!!! note

    अगले section में जारी रखने से पहले दोनों intentional errors को वापस जाकर ठीक करना सुनिश्चित करो।

### सारांश

- Path strings बनाम Path objects: Strings सिर्फ text हैं, Path objects smart फ़ाइल references हैं
- `file()` method एक string path को एक Path object में convert करती है जिसके साथ Nextflow काम कर सकता है
- तुम `name`, `simpleName`, `extension`, और `parent` जैसी फ़ाइल properties तक पहुँच सकते हो [फ़ाइल attributes का उपयोग करके](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- strings के बजाय Path objects का उपयोग करना Nextflow को तुम्हारे workflow में फ़ाइलों को ठीक से manage करने की अनुमति देता है
- Process Input Outcomes: उचित फ़ाइल हैंडलिंग के लिए Path objects की आवश्यकता होती है, strings की नहीं, यह सुनिश्चित करने के लिए कि फ़ाइलें सही तरीके से staged और processes द्वारा उपयोग के लिए accessible हैं।

---

## 2. remote फ़ाइलों का उपयोग करना

Nextflow की एक मुख्य विशेषता local फ़ाइलों (उसी machine पर) से internet पर accessible remote फ़ाइलों में seamlessly switch करने की क्षमता है।

यदि तुम इसे सही तरीके से कर रहे हो, तो तुम्हें विभिन्न स्थानों से आने वाली फ़ाइलों को accommodate करने के लिए अपने workflow के logic को कभी भी बदलने की आवश्यकता नहीं होनी चाहिए।
तुम्हें बस इतना करना है कि जब तुम इसे workflow को supply कर रहे हो तो फ़ाइल path में उपयुक्त prefix specify करो।

उदाहरण के लिए, `/path/to/data` में कोई prefix नहीं है, यह दर्शाता है कि यह एक 'normal' local फ़ाइल path है, जबकि `s3://path/to/data` में `s3://` prefix शामिल है, यह दर्शाता है कि यह Amazon के S3 object storage में स्थित है।

कई अलग-अलग protocols समर्थित हैं:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

इनमें से किसी का भी उपयोग करने के लिए, बस string में relevant prefix specify करो, जिसे फिर technically एक Uniform Resource Identifier (URI) कहा जाता है न कि फ़ाइल path।
Nextflow authentication और फ़ाइलों को सही जगह पर staging, downloading या uploading और अन्य सभी फ़ाइल ऑपरेशन जो तुम expect करोगे, को संभालेगा।

इस system की मुख्य ताकत यह है कि यह हमें किसी भी pipeline logic को बदले बिना environments के बीच switch करने में सक्षम बनाता है।
उदाहरण के लिए, तुम एक छोटे, local test set के साथ develop कर सकते हो इससे पहले कि तुम remote storage में स्थित एक full-scale test set पर switch करो बस URI को बदलकर।

### 2.1. internet से एक फ़ाइल का उपयोग करो

चलो इसे test करते हैं हमारे workflow को प्रदान कर रहे local path को एक HTTPS path के साथ switch करके जो Github में stored उसी डेटा की एक copy की ओर इशारा करता है।

!!! warning

    यह केवल तभी काम करेगा जब तुम्हारे पास एक active internet connection हो।

`main.nf` को फिर से खोलो और input path को निम्नानुसार बदलो:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Using a remote file from the internet
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
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

??? success "Command output"

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

console output में एक अंतर यह है कि path object class अब `nextflow.file.http.XPath` है, जबकि local path के लिए class `sun.nio.fs.UnixPath` थी।
तुम्हें इन classes को याद रखने की आवश्यकता नहीं है; हम बस यह demonstrate करने के लिए इसका उल्लेख करते हैं कि Nextflow विभिन्न स्थानों की पहचान करता है और उन्हें उचित रूप से संभालता है।

पर्दे के पीछे, Nextflow ने फ़ाइल को work directory के भीतर स्थित एक staging directory में download किया।
उस staged फ़ाइल को फिर एक local फ़ाइल के रूप में treat किया जा सकता है और relevant process directory में symlinked किया जा सकता है।

तुम यह verify कर सकते हो कि यह यहाँ हुआ process के hash value पर स्थित working directory की contents को देखकर।

??? abstract "Work directory contents"

    यदि process hash `8a/2ab7ca` था, तो तुम work directory को explore कर सकते हो:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    symlink remote फ़ाइल की एक staged copy की ओर इशारा करता है जिसे Nextflow ने automatically download किया।

ध्यान दो कि बड़ी फ़ाइलों के लिए, downloading step local फ़ाइलों पर चलने की तुलना में कुछ अतिरिक्त समय लेगा।
हालाँकि, Nextflow check करता है कि क्या उसके पास पहले से एक staged copy है ताकि अनावश्यक downloads से बचा जा सके।
तो यदि तुम उसी फ़ाइल पर फिर से चलाते हो और staged फ़ाइल को delete नहीं किया है, तो Nextflow staged copy का उपयोग करेगा।

यह दिखाता है कि Nextflow का उपयोग करके local और remote डेटा के बीच switch करना कितना आसान है, जो Nextflow की एक मुख्य विशेषता है।

!!! note

    इस सिद्धांत का एक महत्वपूर्ण अपवाद यह है कि तुम HTTPS के साथ glob patterns या directory paths का उपयोग नहीं कर सकते क्योंकि HTTPS कई फ़ाइलों को list नहीं कर सकता, इसलिए तुम्हें exact फ़ाइल URLs specify करने होंगे।
    हालाँकि, अन्य storage protocols जैसे blob storage (`s3://`, `az://`, `gs://`) globs और directory paths दोनों का उपयोग कर सकते हैं।

    यहाँ बताया गया है कि तुम cloud storage के साथ glob patterns का उपयोग कैसे कर सकते हो:

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // S3 with glob patterns - would match multiple files
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage with glob patterns
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage with glob patterns
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    हम तुम्हें अगले section में व्यावहारिक रूप से globs के साथ काम करने का तरीका दिखाएंगे।

### 2.2. local फ़ाइल पर वापस switch करो

हम इस side quest के बाकी हिस्सों के लिए अपनी local उदाहरण फ़ाइलों का उपयोग करने जा रहे हैं, तो चलो workflow input को मूल फ़ाइल पर वापस switch करते हैं:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### सारांश

- Remote डेटा को एक URI (HTTP, FTP, S3, Azure, Google Cloud) का उपयोग करके access किया जाता है
- Nextflow automatically डेटा को download और सही जगह पर stage करेगा, जब तक ये paths processes को feed किए जा रहे हैं
- Remote फ़ाइलों को download या upload करने के लिए logic मत लिखो!
- Local और remote फ़ाइलें अलग-अलग object types उत्पन्न करती हैं लेकिन समान रूप से काम करती हैं
- **महत्वपूर्ण**: HTTP/HTTPS केवल single फ़ाइलों के साथ काम करता है (कोई glob patterns नहीं)
- Cloud storage (S3, Azure, GCS) single फ़ाइलों और glob patterns दोनों को support करता है
- तुम code logic को बदले बिना local और remote डेटा sources के बीच seamlessly switch कर सकते हो (जब तक protocol तुम्हारे आवश्यक operations को support करता है)

---

## 3. `fromPath()` channel factory का उपयोग करना

अब तक हम एक बार में एक single फ़ाइल के साथ काम कर रहे हैं, लेकिन Nextflow में, हम आमतौर पर process करने के लिए कई input फ़ाइलों के साथ एक input channel बनाना चाहेंगे।

ऐसा करने का एक naive तरीका `file()` method को [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) के साथ combine करना होगा जैसे:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

यह काम करता है, लेकिन यह clunky है।

!!! tip "`file()` बनाम `channel.fromPath()` का उपयोग कब करें"

    - `file()` का उपयोग करो जब तुम्हें direct manipulation के लिए एक single Path object की आवश्यकता हो (यह check करना कि कोई फ़ाइल मौजूद है, इसकी attributes को पढ़ना, या एक single process invocation में pass करना)
    - `channel.fromPath()` का उपयोग करो जब तुम्हें एक channel की आवश्यकता हो जो कई फ़ाइलों को hold कर सके, विशेष रूप से glob patterns के साथ, या जब फ़ाइलें कई processes के माध्यम से flow होंगी

यहीं पर [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) आता है: एक सुविधाजनक channel factory जो एक या अधिक static फ़ाइल strings के साथ-साथ glob patterns से एक channel generate करने के लिए हमें आवश्यक सभी functionality को bundle करता है।

### 3.1. channel factory जोड़ो

चलो अपने workflow को `channel.fromPath` का उपयोग करने के लिए update करते हैं।

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Print file attributes
        /* Comment these out for now, we'll come back to them!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        // COUNT_LINES(myFile)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

हमने अभी के लिए attributes को print करने वाले code को भी comment out कर दिया है, और बस filename को print करने के लिए एक `.view` statement जोड़ा है।

workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

जैसा कि तुम देख सकते हो, फ़ाइल path को channel में एक `Path` type object के रूप में load किया जा रहा है।
यह `file()` ने जो किया होगा उसके समान है, सिवाय इसके कि अब हमारे पास एक channel है जिसमें हम चाहें तो और फ़ाइलें load कर सकते हैं।

`channel.fromPath()` का उपयोग करना फ़ाइलों की एक list द्वारा populated एक नया channel बनाने का एक सुविधाजनक तरीका है।

### 3.2. channel में फ़ाइलों की attributes देखो

channel factory का उपयोग करने के हमारे पहले pass में, हमने code को सरल बनाया और बस फ़ाइल नाम को print किया।

चलो पूर्ण फ़ाइल attributes को print करने के लिए वापस जाते हैं:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

हम यह verify करने के लिए `COUNT_LINES` process call को भी फिर से enable कर रहे हैं कि हमारे channel-based approach के साथ फ़ाइल processing अभी भी सही तरीके से काम करती है।

workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

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

और वहाँ तुम हो, पहले जैसे ही परिणाम लेकिन अब हमारे पास channel में फ़ाइल है, तो हम और जोड़ सकते हैं।

### 3.3. कई फ़ाइलों को match करने के लिए एक glob का उपयोग करना

कई तरीके हैं जिनसे हम channel में और फ़ाइलें load कर सकते हैं।
यहाँ हम तुम्हें दिखाने जा रहे हैं कि glob patterns का उपयोग कैसे करें, जो wildcard characters के आधार पर फ़ाइल और directory नामों को match और retrieve करने का एक सुविधाजनक तरीका है।
इन patterns को match करने की प्रक्रिया को "globbing" या "filename expansion" कहा जाता है।

!!! note

    जैसा कि पहले noted किया गया है, Nextflow अधिकांश मामलों में input और output फ़ाइलों को manage करने के लिए globbing को support करता है, सिवाय HTTPS filepaths के साथ क्योंकि HTTPS कई फ़ाइलों को list नहीं कर सकता।

मान लो हम एक दिए गए patient, `patientA` से जुड़ी फ़ाइलों की एक pair में दोनों फ़ाइलों को retrieve करना चाहते हैं:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

चूंकि filenames के बीच एकमात्र अंतर replicate number है, _यानी_ `R` के बाद की संख्या, हम निम्नानुसार संख्या के लिए खड़े होने के लिए wildcard character `*` का उपयोग कर सकते हैं:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

यही glob pattern है जिसकी हमें आवश्यकता है।

अब हमें बस इतना करना है कि channel factory में फ़ाइल path को उस glob pattern का उपयोग करने के लिए update करें जैसे:

=== "After"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow automatically पहचान लेगा कि यह एक glob pattern है और इसे उचित रूप से संभालेगा।

इसे test करने के लिए workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

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

जैसा कि तुम देख सकते हो, अब हमारे channel में दो Path objects हैं, जो दिखाता है कि Nextflow ने filename expansion सही तरीके से किया है, और दोनों फ़ाइलों को expected के रूप में load और process किया है।

इस method का उपयोग करके, हम बस glob pattern को बदलकर जितनी चाहें उतनी या कम फ़ाइलें retrieve कर सकते हैं। यदि हमने इसे अधिक generous बनाया, उदाहरण के लिए filenames के सभी variable parts को `*` से replace करके (_जैसे_ `data/patient*_rep*_*_R*_001.fastq.gz`) हम `data` directory में सभी उदाहरण फ़ाइलों को grab कर सकते हैं।

### सारांश

- `channel.fromPath()` एक pattern से matching फ़ाइलों के साथ एक channel बनाता है
- प्रत्येक फ़ाइल को channel में एक अलग element के रूप में emit किया जाता है
- हम कई फ़ाइलों को match करने के लिए एक glob pattern का उपयोग कर सकते हैं
- फ़ाइलें automatically पूर्ण attributes के साथ Path objects में convert हो जाती हैं
- `.view()` method channel contents के inspection की अनुमति देता है

---

## 4. filenames से बुनियादी metadata निकालना

अधिकांश वैज्ञानिक domains में, डेटा वाली फ़ाइलों के नामों में metadata encoded होना बहुत आम है।
उदाहरण के लिए, bioinformatics में, sequencing डेटा वाली फ़ाइलों का नाम अक्सर इस तरह से रखा जाता है जो sample, condition, replicate, और read number के बारे में जानकारी encode करता है।

यदि filenames एक consistent convention के अनुसार constructed हैं, तो तुम उस metadata को एक standardized तरीके से extract कर सकते हो और अपने विश्लेषण के दौरान इसका उपयोग कर सकते हो।
यह एक बड़ा 'if' है, निश्चित रूप से, और तुम्हें बहुत सावधान रहना चाहिए जब भी तुम filename structure पर भरोसा करते हो; लेकिन वास्तविकता यह है कि यह approach बहुत व्यापक रूप से उपयोग किया जाता है, तो चलो देखते हैं कि यह Nextflow में कैसे किया जाता है।

हमारे उदाहरण डेटा के मामले में, हम जानते हैं कि filenames में consistently structured metadata शामिल है।
उदाहरण के लिए, filename `patientA_rep1_normal_R2_001` निम्नलिखित को encode करता है:

- patient ID: `patientA`
- replicate ID: `rep1`
- sample type: `normal` (`tumor` के विपरीत)
- read set: `R1` (`R2` के विपरीत)

हम तीन steps में इस जानकारी को retrieve करने के लिए अपने workflow को modify करने जा रहे हैं:

1. फ़ाइल की `simpleName` को retrieve करो, जिसमें metadata शामिल है
2. `tokenize()` नामक एक method का उपयोग करके metadata को अलग करो
3. metadata को organize करने के लिए एक map का उपयोग करो

!!! warning

    तुम्हें कभी भी filenames में sensitive जानकारी encode नहीं करनी चाहिए, जैसे patient नाम या अन्य identifying characteristics, क्योंकि यह patient privacy या अन्य relevant security restrictions को compromise कर सकता है।

### 4.1. `simpleName` को retrieve करो

`simpleName` एक फ़ाइल attribute है जो filename से मेल खाता है जिसे इसके path और extension से stripped किया गया है।

workflow में निम्नलिखित edits करो:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

यह `simpleName` को retrieve करता है और इसे एक `map()` operation का उपयोग करके पूर्ण फ़ाइल object के साथ associate करता है।

यह test करने के लिए workflow चलाओ कि यह काम करता है:

```bash
nextflow run main.nf
```

??? success "Command output"

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

channel में प्रत्येक element अब एक tuple है जिसमें `simpleName` और मूल फ़ाइल object है।

### 4.2. `simplename` से metadata निकालो

इस बिंदु पर, हम जो metadata चाहते हैं वह `simplename` में embedded है, लेकिन हम सीधे individual items तक नहीं पहुँच सकते।
तो हमें `simplename` को इसके components में split करने की आवश्यकता है।
सौभाग्य से, वे components मूल filename में बस underscores द्वारा अलग किए गए हैं, तो हम `tokenize()` नामक एक सामान्य Nextflow method को apply कर सकते हैं जो इस कार्य के लिए perfect है।

workflow में निम्नलिखित edits करो:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

`tokenize()` method `simpleName` string को जहाँ भी underscores मिलते हैं वहाँ split करेगा, और substrings वाली एक list return करेगा।

workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

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

अब हमारे channel में प्रत्येक element के लिए tuple में metadata की list (_जैसे_ `[patientA, rep1, normal, R1, 001]`) और मूल फ़ाइल object है।

यह बहुत अच्छा है!
हमने अपनी patient जानकारी को एक single string से strings की एक list में तोड़ दिया है।
अब हम patient जानकारी के प्रत्येक भाग को अलग से संभाल सकते हैं।

### 4.3. metadata को organize करने के लिए एक map का उपयोग करो

हमारा metadata इस समय सिर्फ एक flat list है।
इसका उपयोग करना काफी आसान है लेकिन पढ़ना मुश्किल है।

```console
[patientA, rep1, normal, R1, 001]
```

index 3 पर item क्या है? क्या तुम metadata structure की मूल व्याख्या को वापस refer किए बिना बता सकते हो?

यह एक key-value store का उपयोग करने का एक शानदार अवसर है, जहाँ प्रत्येक item में keys का एक set और उनके associated values होते हैं, तो तुम आसानी से corresponding value प्राप्त करने के लिए प्रत्येक key को refer कर सकते हो।

हमारे उदाहरण में, इसका मतलब है इस organization से जाना:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

इस एक में:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

Nextflow में, इसे [map](https://nextflow.io/docs/latest/script.html#maps) कहा जाता है।

चलो अब अपनी flat list को एक map में convert करते हैं।
workflow में निम्नलिखित edits करो:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Load files with channel.fromPath
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

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

यहाँ मुख्य परिवर्तन हैं:

- **Destructuring assignment**: `def (patient, replicate, type, readNum) = ...` tokenized values को एक line में named variables में extract करता है
- **Map literal syntax**: `[id: patient, replicate: ...]` एक map बनाता है जहाँ प्रत्येक key (जैसे `id`) एक value (जैसे `patient`) से associated है
- **Nested structure**: outer list `[..., myFile]` metadata map को मूल फ़ाइल object के साथ pair करती है

हमने कुछ metadata strings को `replace()` नामक एक string replacement method का उपयोग करके भी सरल बनाया ताकि कुछ characters को हटाया जा सके जो अनावश्यक हैं (_जैसे_ `replicate.replace('rep', '')` replicate IDs से केवल संख्या रखने के लिए)।

चलो workflow को फिर से चलाते हैं:

```bash
nextflow run main.nf
```

??? success "Command output"

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

अब metadata साफ-सुथरे ढंग से labeled है (_जैसे_ `[id:patientA, replicate:1, type:normal, readNum:2]`) तो यह बताना बहुत आसान है कि क्या क्या है।

workflow में metadata के elements का वास्तव में उपयोग करना भी बहुत आसान होगा, और हमारे code को पढ़ने में आसान और अधिक maintainable बना देगा।

### सारांश

- हम Nextflow में एक पूर्ण programming language की शक्ति के साथ filenames को संभाल सकते हैं
- हम relevant जानकारी निकालने के लिए filenames को strings के रूप में treat कर सकते हैं
- `tokenize()` और `replace()` जैसी methods का उपयोग हमें filename में strings को manipulate करने की अनुमति देता है
- `.map()` operation structure को preserve करते हुए channel elements को transform करता है
- Structured metadata (maps) positional lists की तुलना में code को अधिक readable और maintainable बनाता है

अगला, हम देखेंगे कि paired डेटा फ़ाइलों को कैसे संभालें।

---

## 5. paired डेटा फ़ाइलों को संभालना

कई प्रयोगात्मक डिज़ाइन paired डेटा फ़ाइलें उत्पन्न करते हैं जो स्पष्ट रूप से paired तरीके से संभाले जाने से लाभान्वित होती हैं।
उदाहरण के लिए, bioinformatics में, sequencing डेटा अक्सर paired reads के रूप में उत्पन्न होता है, जिसका अर्थ है sequence strings जो DNA के एक ही fragment से उत्पन्न होती हैं (अक्सर 'forward' और 'reverse' कहा जाता है क्योंकि उन्हें विपरीत ends से पढ़ा जाता है)।

यह हमारे उदाहरण डेटा का मामला है, जहाँ R1 और R2 reads के दो sets को refer करते हैं।

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow इस तरह की paired फ़ाइलों के साथ काम करने के लिए `channel.fromFilePairs()` नामक एक specialized channel factory प्रदान करता है, जो automatically एक shared naming pattern के आधार पर फ़ाइलों को group करता है। यह तुम्हें कम प्रयास के साथ paired फ़ाइलों को अधिक कसकर associate करने की अनुमति देता है।

हम इसका लाभ उठाने के लिए अपने workflow को modify करने जा रहे हैं।
इसमें दो steps लगने वाले हैं:

1. channel factory को `channel.fromFilePairs()` में switch करो
2. metadata को extract और map करो

### 5.1. channel factory को `channel.fromFilePairs()` में switch करो

`channel.fromFilePairs` का उपयोग करने के लिए, हमें वह pattern specify करना होगा जिसका उपयोग Nextflow को एक pair में दो members की पहचान करने के लिए करना चाहिए।

हमारे उदाहरण डेटा पर वापस जाते हुए, हम naming pattern को निम्नानुसार formalize कर सकते हैं:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

यह हमारे द्वारा पहले उपयोग किए गए glob pattern के समान है, सिवाय इसके कि यह विशेष रूप से substrings (या तो `1` या `2` R के ठीक बाद आने वाला) को enumerate करता है जो pair के दो members की पहचान करते हैं।

चलो workflow `main.nf` को तदनुसार update करते हैं:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comment out the mapping for now, we'll come back to it!
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

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
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

हमने channel factory को switch किया है और फ़ाइल matching pattern को adapt किया है, और जब हम इस पर थे, हमने map operation को comment out कर दिया।
हम इसे बाद में वापस जोड़ेंगे, कुछ modifications के साथ।

इसे test करने के लिए workflow चलाओ:

```bash
nextflow run main.nf
```

??? failure "Command output"

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

उह-ओह, इस बार run fail हो गया!

error message का relevant bit यहाँ है:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

ऐसा इसलिए है क्योंकि हमने channel factory को बदल दिया है।
अब तक, मूल input channel में केवल फ़ाइल paths थे।
हम जो सभी metadata manipulation कर रहे थे वह वास्तव में channel contents को प्रभावित नहीं करता था।

अब जब हम `.fromFilePairs` channel factory का उपयोग कर रहे हैं, तो resulting channel की contents अलग हैं।
हम केवल एक channel element देखते हैं, जो दो items वाले एक tuple से composed है: दो फ़ाइलों द्वारा shared `simpleName` का भाग, जो एक identifier के रूप में कार्य करता है, और दो फ़ाइल objects वाला एक tuple, format `id, [ file1, file2 ]` में।

यह बहुत अच्छा है, क्योंकि Nextflow ने shared prefix की जाँच करके और इसे एक patient identifier के रूप में उपयोग करके patient नाम निकालने का कठिन काम किया है।

हालाँकि, यह हमारे वर्तमान workflow को तोड़ता है।
यदि हम अभी भी process को बदले बिना `COUNT_LINES` को उसी तरह चलाना चाहते हैं, तो हमें फ़ाइल paths को extract करने के लिए एक mapping operation apply करना होगा।
लेकिन हम ऐसा नहीं करने जा रहे हैं, क्योंकि हमारा ultimate लक्ष्य एक अलग process, `ANALYZE_READS`, का उपयोग करना है, जो फ़ाइल pairs को उचित रूप से संभालता है।

तो चलो बस `COUNT_LINES` की call को comment out (या delete) करते हैं और आगे बढ़ते हैं।

=== "After"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

तुम `COUNT_LINES` include statement को भी comment out या delete कर सकते हो, लेकिन इसका कोई functional effect नहीं होगा।

अब चलो workflow को फिर से चलाते हैं:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

यय, इस बार workflow succeed होता है!

हालाँकि, हमें अभी भी `id` field से बाकी metadata निकालने की आवश्यकता है।

### 5.2. फ़ाइल pairs से metadata को extract और organize करो

पहले से हमारा `map` operation काम नहीं करेगा क्योंकि यह डेटा structure से match नहीं करता है, लेकिन हम इसे काम करने के लिए modify कर सकते हैं।

हमारे पास पहले से ही string में वास्तविक patient identifier तक पहुँच है जिसे `fromFilePairs()` ने एक identifier के रूप में उपयोग किया, तो हम Path object से `simpleName` प्राप्त किए बिना metadata को extract करने के लिए इसका उपयोग कर सकते हैं जैसा कि हमने पहले किया था।

workflow में map operation को uncomment करो और निम्नलिखित edits करो:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Load files with channel.fromFilePairs
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

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comment out the mapping for now, we'll come back to it!
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

इस बार map सिर्फ `myFile` के बजाय `id, files` से शुरू होता है, और `tokenize()` को `myFile.simpleName` के बजाय `id` पर apply किया जाता है।

यह भी ध्यान दो कि हमने `tokenize()` line से `readNum` को drop कर दिया है; कोई भी substrings जिन्हें हम विशेष रूप से name नहीं करते (बाएं से शुरू करते हुए) चुपचाप drop हो जाएंगे।
हम ऐसा कर सकते हैं क्योंकि paired फ़ाइलें अब कसकर associated हैं, तो हमें अब metadata map में `readNum` की आवश्यकता नहीं है।

चलो workflow चलाते हैं:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

और वहाँ यह है: हमारे पास output tuple की पहली position में metadata map (`[id:patientA, replicate:1, type:normal]`) है, इसके बाद paired फ़ाइलों का tuple है, जैसा कि intended था।

निश्चित रूप से, यह केवल फ़ाइलों की उस specific pair को pick up और process करेगा।
यदि तुम कई pairs को process करने के साथ experiment करना चाहते हो, तो तुम input pattern में wildcards जोड़ने की कोशिश कर सकते हो और देख सकते हो कि क्या होता है।
उदाहरण के लिए, `data/patientA_rep1_*_R{1,2}_001.fastq.gz` का उपयोग करने की कोशिश करो

### सारांश

- [`channel.fromFilePairs()` automatically संबंधित फ़ाइलों को ढूंढता और pair करता है](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- यह तुम्हारी pipeline में paired-end reads को संभालने को सरल बनाता है
- Paired फ़ाइलों को `[id, [file1, file2]]` tuples के रूप में grouped किया जा सकता है
- Metadata extraction individual फ़ाइलों के बजाय paired फ़ाइल ID से किया जा सकता है

---

## 6. processes में फ़ाइल operations का उपयोग करना

अब चलो यह सब एक साथ एक सरल process में डालते हैं ताकि यह reinforce किया जा सके कि Nextflow process के अंदर फ़ाइल operations का उपयोग कैसे करें।

हम तुम्हें `ANALYZE_READS` नामक एक pre-written process module प्रदान करते हैं जो metadata और input फ़ाइलों की एक pair का एक tuple लेता है और उनका विश्लेषण करता है।
हम कल्पना कर सकते हैं कि यह sequence alignment, या variant calling या इस डेटा प्रकार के लिए समझ में आने वाला कोई अन्य step कर रहा है।

चलो शुरू करते हैं।

### 6.1. process को import करो और code की जाँच करो

workflow में इस process का उपयोग करने के लिए, हमें बस workflow block से पहले एक module include statement जोड़ना होगा।

workflow में निम्नलिखित edit करो:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

तुम इसके code की जाँच करने के लिए module फ़ाइल खोल सकते हो:

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

!!! note

    `tag` और `publishDir` directives string interpolation (`"${...}"`) के बजाय closure syntax (`{ ... }`) का उपयोग करते हैं।
    ऐसा इसलिए है क्योंकि ये directives input variables (`meta`) को reference करते हैं जो runtime तक उपलब्ध नहीं हैं।
    closure syntax evaluation को तब तक defer करता है जब तक process वास्तव में चलता नहीं है।

!!! note

    हम convention द्वारा अपने metadata map को `meta` कह रहे हैं।
    meta maps में गहराई से dive के लिए, [Metadata and meta maps](./metadata.md) side quest देखो।

### 6.2. workflow में process को call करो

अब जब process workflow के लिए उपलब्ध है, हम इसे चलाने के लिए `ANALYZE_READS` process में एक call जोड़ सकते हैं।

इसे हमारे उदाहरण डेटा पर चलाने के लिए, हमें दो चीजें करनी होंगी:

1. remapped channel को एक नाम दो
2. process में एक call जोड़ो

#### 6.2.1. remapped input channel को नाम दो

हमने पहले mapping manipulations को सीधे input channel पर apply किया था।
remapped contents को `ANALYZE_READS` process को feed करने के लिए (और ऐसा करने के लिए जो clear और पढ़ने में आसान हो) हम `ch_samples` नामक एक नया channel बनाना चाहते हैं।

हम [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) operator का उपयोग करके ऐसा कर सकते हैं।

मुख्य workflow में, `.view()` operator को `.set { ch_samples }` से replace करो, और एक line जोड़ो जो test करे कि हम channel को नाम से refer कर सकते हैं।

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Load files with channel.fromFilePairs
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

        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Load files with channel.fromFilePairs
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

चलो इसे चलाते हैं:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

यह confirm करता है कि अब हम channel को नाम से refer कर सकते हैं।

#### 6.2.2. डेटा पर process को call करो

अब चलो वास्तव में `ch_samples` channel पर `ANALYZE_READS` process को call करते हैं।

मुख्य workflow में, निम्नलिखित code changes करो:

=== "After"

    ```groovy title="main.nf" linenums="23"
        // Run the analysis
        ANALYZE_READS(ch_samples)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="23"
        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

चलो इसे चलाते हैं:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

यह process अपने outputs को एक `results` directory में publish करने के लिए set up है, तो वहाँ एक नज़र डालो।

??? abstract "Directory and file contents"

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

process ने हमारे inputs लिए और patient metadata वाली एक नई फ़ाइल बनाई, जैसा कि designed था।
शानदार!

### 6.3. कई और patients को शामिल करो

निश्चित रूप से, यह सिर्फ एक single patient के लिए फ़ाइलों की एक single pair को process कर रहा है, जो बिल्कुल वह high throughput नहीं है जिसकी तुम Nextflow के साथ उम्मीद कर रहे हो।
तुम शायद एक बार में बहुत अधिक डेटा process करना चाहोगे।

याद रखो `channel.fromPath()` input के रूप में एक _glob_ स्वीकार करता है, जिसका अर्थ है कि यह pattern से matching किसी भी संख्या में फ़ाइलों को स्वीकार कर सकता है।
इसलिए यदि हम सभी patients को शामिल करना चाहते हैं, तो हम बस अधिक patients को शामिल करने के लिए input string को modify कर सकते हैं, जैसा कि पहले passing में noted किया गया था।

चलो मान लेते हैं कि हम जितना संभव हो उतना greedy होना चाहते हैं।
workflow में निम्नलिखित edits करो:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

pipeline को फिर से चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

results directory में अब सभी उपलब्ध डेटा के लिए परिणाम होने चाहिए।

??? abstract "Directory contents"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

सफलता! हमने एक बार में सभी patients का विश्लेषण किया है! सही?

शायद नहीं।
यदि तुम अधिक ध्यान से देखो, तो हमारे पास एक समस्या है: हमारे पास patientA के लिए दो replicates हैं, लेकिन केवल एक output फ़ाइल!
हम हर बार output फ़ाइल को overwrite कर रहे हैं।

### 6.4. published फ़ाइलों को unique बनाओ

चूंकि हमारे पास patient metadata तक पहुँच है, हम इसका उपयोग differentiating metadata को शामिल करके published फ़ाइलों को unique बनाने के लिए कर सकते हैं, या तो directory structure में या filenames में ही।

workflow में निम्नलिखित change करो:

=== "After"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Before"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

यहाँ हम sample types और replicates के लिए account करने के लिए अतिरिक्त directory levels का उपयोग करने का विकल्प दिखाते हैं, लेकिन तुम filename level पर भी ऐसा करने के साथ experiment कर सकते हो।

अब pipeline को एक बार और चलाओ, लेकिन अपने आप को एक clean workspace देने के लिए पहले results directory को हटाना सुनिश्चित करो:

```bash
rm -r results
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

अब results directory को check करो:

??? abstract "Directory contents"

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

और वहाँ यह है, हमारा सभी metadata, साफ-सुथरे ढंग से organized। यह सफलता है!

एक बार जब तुम्हारे पास तुम्हारा metadata इस तरह एक map में loaded हो जाता है तो तुम बहुत कुछ कर सकते हो:

1. patient attributes के आधार पर organized output directories बनाओ
2. patient properties के आधार पर processes में निर्णय लो
3. metadata values के आधार पर डेटा को split, join, और recombine करो

metadata को explicit रखने और डेटा के साथ attached रखने का यह pattern (filenames में encoded के बजाय) Nextflow में एक core best practice है जो robust, maintainable विश्लेषण workflows बनाने में सक्षम बनाता है।
तुम [Metadata and meta maps](./metadata.md) side quest में इसके बारे में अधिक जान सकते हो।

### सारांश

- `publishDir` directive metadata values के आधार पर outputs को organize कर सकता है
- tuples में metadata structured परिणामों के organization को सक्षम बनाता है
- यह approach clear डेटा provenance के साथ maintainable workflows बनाता है
- Processes input के रूप में metadata और फ़ाइलों के tuples ले सकते हैं
- `tag` directive execution logs में process identification प्रदान करता है
- Workflow structure channel creation को process execution से अलग करता है

---

## सारांश

इस side quest में, तुमने सीखा है कि Nextflow में फ़ाइलों के साथ कैसे काम करें, बुनियादी operations से लेकर फ़ाइलों के संग्रह को संभालने की अधिक उन्नत तकनीकों तक।

अपने स्वयं के काम में इन तकनीकों को apply करना तुम्हें अधिक कुशल और maintainable workflows बनाने में सक्षम बनाएगा, विशेष रूप से जब जटिल naming conventions के साथ बड़ी संख्या में फ़ाइलों के साथ काम कर रहे हो।

### मुख्य patterns

1.  **बुनियादी फ़ाइल ऑपरेशन:** हमने `file()` के साथ Path objects बनाए और name, extension, और parent directory जैसी फ़ाइल attributes तक पहुँचे, strings और Path objects के बीच अंतर सीखा।

    - `file()` के साथ एक Path object बनाओ

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - फ़ाइल attributes प्राप्त करो

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **remote फ़ाइलों का उपयोग करना**: हमने सीखा कि URIs का उपयोग करके local और remote फ़ाइलों के बीच पारदर्शी रूप से कैसे switch करें, Nextflow की workflow logic को बदले बिना विभिन्न sources से फ़ाइलों को संभालने की क्षमता को demonstrate करते हुए।

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

3.  **`fromPath()` channel factory का उपयोग करके फ़ाइलें load करना:** हमने `channel.fromPath()` के साथ फ़ाइल patterns से channels बनाए और उनकी फ़ाइल attributes देखीं, object types सहित।

    - एक फ़ाइल pattern से एक channel बनाओ

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - फ़ाइल attributes प्राप्त करो

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **filenames से Patient Metadata निकालना:** हमने filenames से metadata को extract और structure करने के लिए `tokenize()` और `replace()` का उपयोग किया, उन्हें organized maps में convert करते हुए।

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **channel.fromFilePairs के साथ सरलीकरण:** हमने automatically संबंधित फ़ाइलों को pair करने और paired फ़ाइल IDs से metadata extract करने के लिए `channel.fromFilePairs()` का उपयोग किया।

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **processes में फ़ाइल ऑपरेशन का उपयोग करना:** हमने उचित input हैंडलिंग के साथ Nextflow processes में फ़ाइल operations को integrate किया, metadata के आधार पर outputs को organize करने के लिए `publishDir` का उपयोग करते हुए।

    - process inputs के साथ एक meta map को associate करो

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

    - metadata के आधार पर outputs को organize करो

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### अतिरिक्त संसाधन

- [Nextflow Documentation: Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## आगे क्या है?

[Side Quests के menu](./index.md) पर वापस जाओ या list में अगले topic पर जाने के लिए page के निचले दाएं कोने में button पर click करो।
