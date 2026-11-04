# भाग 2: पाइपलाइन को कॉन्फ़िगर करना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[भाग 1](./01_run_nextflow.md) में, तुमने एक पूरी मल्टी-स्टेप पाइपलाइन चलाई जो कंटेनर का उपयोग करके कई इनपुट को समानांतर में प्रोसेस करती है।
अब हम देखेंगे कि `nextflow.config` का उपयोग करके पाइपलाइन के व्यवहार को कैसे कॉन्फ़िगर किया जाए: पहले उस कॉन्फ़िगरेशन फ़ाइल की जांच करके जो हमने तुम्हें पहले से दी है, फिर कॉन्फ़िगरेशन देने के कुछ और तरीके खोजकर, और अंत में यह नियंत्रित करके कि आउटपुट कैसे और कहाँ प्रकाशित होते हैं।

---

## 1. मुख्य कॉन्फ़िगरेशन फ़ाइल की जांच करना

Nextflow स्वचालित रूप से working directory से `nextflow.config` उठाता है और इसकी सेटिंग्स हर रन पर लागू करता है।

हम तुम्हें एक कॉन्फ़िगरेशन फ़ाइल देते हैं जो चार क्षेत्रों को कवर करती है: सॉफ़्टवेयर पैकेजिंग, प्रोसेस सेटिंग्स, पाइपलाइन पैरामीटर, और execution प्रोफ़ाइल।

??? full-code "nextflow.config"

    ```groovy title="nextflow.config" linenums="1"
    /*
     * सॉफ़्टवेयर पैकेजिंग
     */
    docker.enabled = true

    /*
     * प्रोसेस सेटिंग्स
     */
    process {
        cpus = 1
        memory = 1.GB
    }

    /*
     * पाइपलाइन पैरामीटर
     */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * प्रोफ़ाइल
     */
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

आओ हर एक को देखें, फिर प्रोफ़ाइल को उपयोग में लाएं और पाइपलाइन को उसके साथ चलाएं।

!!! note "नोट"

    यह कॉन्फ़िग एक ही मशीन पर लोकल execution को कवर करता है।
    Nextflow HPC schedulers (SLURM, PBS, LSF) और cloud executors (AWS Batch, Google Cloud Batch, Azure Batch) को भी सपोर्ट करता है, जो सभी उसी `nextflow.config` मेकेनिज्म के ज़रिए कॉन्फ़िगर होते हैं।
    इन विकल्पों की पूरी जानकारी के लिए [Configure Execution](../config_exec/01_packaging_and_execution.md) कोर्स में [भाग 1: अपने compute वातावरण के अनुसार ढालना](../config_exec/index.md) देखो।

### 1.1. सॉफ़्टवेयर पैकेजिंग

सॉफ़्टवेयर पैकेजिंग वह तरीका है जिससे Nextflow तुम्हारे प्रोसेस को ज़रूरी टूल्स देता है, चाहे वह कंटेनर इमेज हो, Conda वातावरण हो, या कुछ और।

```groovy title="nextflow.config" linenums="1"
/*
 * सॉफ़्टवेयर पैकेजिंग
 */
docker.enabled = true
```

यह लाइन हर प्रोसेस के लिए Docker को सक्षम करती है।
जो भी प्रोसेस `container` निर्देश घोषित करता है, वह निर्दिष्ट इमेज के अंदर चलता है।

### 1.2. प्रोसेस सेटिंग्स

याद रखो कि एक प्रोसेस तुम्हारी पाइपलाइन का एक ही कदम होता है, जैसे `sayHello` या `cowpy`।
Nextflow तुम्हें यह कॉन्फ़िगर करने देता है कि हर प्रोसेस वास्तव में कैसे चलता है: उसे कितना CPU और memory मिलता है, वह कौन सा कंटेनर या Conda वातावरण उपयोग करता है, और भी बहुत कुछ।

```groovy title="nextflow.config" linenums="6"
/*
 * प्रोसेस सेटिंग्स
 */
process {
    cpus = 1
    memory = 1.GB
}
```

यह हर प्रोसेस को एक ही CPU और 1 GB memory तक सीमित करता है।

Nextflow तुम्हें अलग-अलग नामित प्रोसेस या प्रोसेस के समूहों के लिए अलग-अलग मान भी सेट करने देता है; तुम यह [Configure Execution](../config_exec/02_resources_and_retries.md#12-set-resource-allocations-for-a-specific-process) कोर्स के [भाग 2: compute संसाधन और विफलताओं का प्रबंधन](../config_exec/index.md) में सीखोगे।

### 1.3. पाइपलाइन पैरामीटर

पैरामीटर पाइपलाइन के command-line इनपुट हैं, वही `--input`, `--batch` और `--character` फ्लैग जो तुम पहले से सीधे command line पर सेट करते आए हो।
यहाँ उनके लिए डिफ़ॉल्ट सेट करने का मतलब है कि तुम्हें हर बार उन्हें टाइप नहीं करना पड़ता, हालांकि जैसा तुम इस भाग में बाद में देखोगे, उन्हें देने के कुछ और तरीके भी हैं।

```groovy title="nextflow.config" linenums="14"
/*
 * पाइपलाइन पैरामीटर
 */
params {
    input = 'data/greetings.csv'
    batch = 'batch'
    character = 'turkey'
}
```

ये डिफ़ॉल्ट तब काम आते हैं जब command line पर कोई पैरामीटर नहीं दिया जाता, इसलिए बिना किसी फ्लैग के `nextflow run main.nf` चलाना भी काम करता है।

### 1.4. प्रोफ़ाइल

प्रोफ़ाइल तुम्हें एक नाम के तहत सेटिंग्स का एक सेट बंडल करने देते हैं, ताकि तुम हर बार हाथ से मान बदलने की बजाय एक फ्लैग से पूरी कॉन्फ़िगरेशन बदल सको।

```groovy title="nextflow.config" linenums="23"
/*
 * प्रोफ़ाइल
 */
profiles {
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'tux'
    }
    conda {
        docker.enabled = false
        conda.enabled = true
    }
}
```

`test` प्रोफ़ाइल तीन पैरामीटर को ओवरराइड करके पाइपलाइन को एक छोटे, सुपरिभाषित इनपुट सेट के साथ चलाता है; हर nf-core पाइपलाइन त्वरित सत्यापन के लिए इनमें से एक के साथ आती है, और यह एक ऐसा तरीका है जिसे तुम्हारी अपनी पाइपलाइन में भी अपनाना चाहिए।

`conda` प्रोफ़ाइल सॉफ़्टवेयर पैकेजिंग को Docker से Conda में बदल देता है।

तुम command line पर `-profile <name>` पास करके एक प्रोफ़ाइल सक्रिय करते हो।

आओ `test` प्रोफ़ाइल को उपयोग में लाएं।

```bash
nextflow run main.nf -profile test
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [reverent_heisenberg] revision: ce74f81996

    executor >  local (8)
    [3d/8a12c7] sayHello (3)       | 3 of 3 ✔
    [e7/e0934f] convertToUpper (2) | 3 of 3 ✔
    [1d/616569] collectGreetings   | 1 of 1 ✔
    [44/7d46cf] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - test/intermediates/Bonjour-output.txt
        - test/intermediates/Hello-output.txt
        - test/intermediates/Hola-output.txt

      uppercased:
        - test/intermediates/UPPER-Bonjour-output.txt
        - test/intermediates/UPPER-Hello-output.txt
        - test/intermediates/UPPER-Hola-output.txt

      collected: test/intermediates/COLLECTED-test-output.txt

      batch_report: test/test-report.txt

      cowpy_art: test/cowpy-COLLECTED-test-output.txt
    ```

पाइपलाइन `batch = 'test'` और `character = 'tux'` के साथ चलती है।
`results/test/` देखो: batch का नाम अब डायरेक्टरी पाथ का हिस्सा है, और ASCII आर्ट में turkey की जगह tux penguin है।

!!! note "नोट"

    तुम एक साथ कई प्रोफ़ाइल सक्रिय कर सकते हो, और कुछ भी चलाने से पहले पूरी तरह से resolved परिणाम देखने के लिए `nextflow config -profile <name>,<name>` उपयोग कर सकते हो।
    प्रोफ़ाइल को मिलाना, और Nextflow उनके बीच टकराव को कैसे सुलझाता है, इसे [Configure Execution](../config_exec/03_profiles.md) कोर्स के [भाग 3: कॉन्फ़िगरेशन बदलने के लिए प्रोफ़ाइल का उपयोग](../config_exec/index.md) में विस्तार से कवर किया गया है।

### सारांश

तुम जानते हो कि `nextflow.config` फ़ाइल के सबसे सामान्य तत्व क्या करते हैं, और एक प्रोफ़ाइल को कैसे सक्रिय करते हैं।

### आगे क्या है?

मुख्य `nextflow.config` फ़ाइल को बदले बिना कॉन्फ़िगरेशन मान देने के कुछ और तरीके सीखो, जो अलग-अलग रन को कॉन्फ़िगर करने और किसी के साथ सेटिंग्स का एक सटीक सेट साझा करने के लिए उपयोगी हैं।

---

## 2. अतिरिक्त फ़ाइलों के ज़रिए कॉन्फ़िगरेशन देना

`nextflow.config` में डिफ़ॉल्ट सेट करना उन मानों के लिए अच्छा काम करता है जो शायद ही कभी बदलते हैं।
Nextflow तुम्हें दो और लक्षित तरीके भी देता है: एक रन-विशिष्ट कॉन्फ़िगरेशन फ़ाइल जो execution को किसी विशेष वातावरण के अनुसार ढालती है, और एक पैरामीटर फ़ाइल जो किसी सहयोगी के साथ इनपुट मानों का एक सटीक सेट साझा करने के लिए है।

### 2.1. रन-विशिष्ट कॉन्फ़िगरेशन फ़ाइल का उपयोग करना

मान लो तुम पाइपलाइन को एक ऐसी मशीन पर ले जा रहे हो जिसमें Docker नहीं है, और तुम हर प्रोसेस को काम करने के लिए ज़्यादा जगह देना चाहते हो।
सिर्फ उन ओवरराइड के साथ एक नई कॉन्फ़िगरेशन फ़ाइल बनाओ जिनकी तुम्हें ज़रूरत है:

```groovy title="custom.config" linenums="1"
process {
    cpus = 2
    memory = 2.GB
}

docker.enabled = false
conda.enabled = true
```

इसे अपनी मुख्य पाइपलाइन के साथ `-c` से पास करो:

```bash
nextflow run main.nf -c custom.config
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [exotic_cray] revision: ce74f81996

    executor >  local (8)
    [77/e02315] sayHello (1)       | 3 of 3 ✔
    [a6/ccf44b] convertToUpper (3) | 3 of 3 ✔
    [56/fd1296] collectGreetings   | 1 of 1 ✔
    [2c/205a94] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow `custom.config` को पाइपलाइन के अपने `nextflow.config` के ऊपर मर्ज करता है, इसलिए अब हर प्रोसेस को डिफ़ॉल्ट की जगह 2 CPU और 2 GB memory मिलती है, और Docker की जगह Conda के ज़रिए चलता है।
`cowpy` एकमात्र प्रोसेस है जिसके पास अपने कंटेनर के साथ-साथ Conda पैकेज भी घोषित है, इसलिए यही वह प्रोसेस है जिसके लिए तुम Nextflow को वास्तव में एक वातावरण बनाते देखोगे:

```console
Creating env using conda: conda-forge::cowpy==1.1.5 [cache /path/to/work/conda/env-898314d566668b6587ad714ae06b8520]
```

एक छोटी फ़ाइल जो केवल resource allocation और पैकेजिंग को ओवरराइड करती है, बिना पाइपलाइन पैरामीटर को छुए, यही वह तरीका है जिसकी nf-core पाइपलाइन institutional configs से उम्मीद करती हैं।
वास्तविक दुनिया के उदाहरणों के लिए [nf-core/configs](https://github.com/nf-core/configs) रिपॉजिटरी देखो।

यह तुम्हें अपनी सामान्य कॉन्फ़िगरेशन को छुए बिना पाइपलाइन को नए वातावरण के अनुसार ढालने का एक अस्थायी तरीका देता है।

### 2.2. पैरामीटर फ़ाइल का उपयोग करना

मान लो इसके बजाय तुम्हें किसी सहयोगी के साथ रन पैरामीटर का एक सटीक सेट साझा करना है, या किसी प्रकाशन के लिए उन्हें रिकॉर्ड करना है।

Nextflow तुम्हें YAML या JSON फॉर्मेट में [parameter files](https://nextflow.io/docs/latest/config.html#parameter-file) देने की सुविधा देता है, जो मानों का एक सटीक, reproducible सेट वितरित करने का एक सरल तरीका है।

`test-params.yaml` नाम की एक पैरामीटर फ़ाइल तुम्हारी working directory में पहले से मौजूद है:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

इस सिंटैक्स में `nextflow.config` में उपयोग होने वाले equal signs (`=`) की जगह colons (`:`) का उपयोग होता है, क्योंकि यह फ़ाइल Groovy की बजाय plain YAML है।

!!! info "जानकारी"

    एक JSON संस्करण, `test-params.json`, भी दिया गया है। इसे अपने आप आज़माने के लिए स्वतंत्र हो; इसे पास करने का सिंटैक्स बिल्कुल एक जैसा है।

फ़ाइल को `-params-file` से पास करो:

```bash
nextflow run main.nf -params-file test-params.yaml
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sharp_faraday] revision: ce74f81996

    executor >  local (8)
    [1c/9ff63e] sayHello (1)       | 3 of 3 ✔
    [3b/bb5691] convertToUpper (2) | 3 of 3 ✔
    [cd/2c1f6e] collectGreetings   | 1 of 1 ✔
    [89/c333bc] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - yaml/intermediates/Bonjour-output.txt
        - yaml/intermediates/Hello-output.txt
        - yaml/intermediates/Hola-output.txt

      uppercased:
        - yaml/intermediates/UPPER-Bonjour-output.txt
        - yaml/intermediates/UPPER-Hello-output.txt
        - yaml/intermediates/UPPER-Hola-output.txt

      collected: yaml/intermediates/COLLECTED-yaml-output.txt

      batch_report: yaml/yaml-report.txt

      cowpy_art: yaml/cowpy-COLLECTED-yaml-output.txt
    ```

??? abstract "फ़ाइल सामग्री"

    ```console title="results/yaml/cowpy-COLLECTED-yaml-output.txt"
     _________
    / BONJOUR \
    | HOLA    |
    \ HELLO   /
     ---------
    \                             .       .
     \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
       \                 |    \  \ - ~ ~ - /  /    |
             _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
           ---------   \__/                        \__/
          .'  O    \     /               /       \  "
         (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

एक पैरामीटर फ़ाइल विशेष रूप से तब मूल्यवान होती है जब पाइपलाइन में कुछ से ज़्यादा पैरामीटर हों: यह तुम्हें उन सभी को एक साथ देने देती है, बिना किसी लंबे-चौड़े command line या workflow स्क्रिप्ट में किसी बदलाव के, और इसे अपने परिणामों के साथ वितरित करना आसान है।

### सारांश

तुम कॉन्फ़िगरेशन देने के दो और तरीके जानते हो: एक रन-विशिष्ट कॉन्फ़िगरेशन फ़ाइल जो execution को नए वातावरण के अनुसार ढालती है, और एक पैरामीटर फ़ाइल जो सटीक, reproducible इनपुट मान साझा करने के लिए है।

### आगे क्या है?

यह सीखो कि तुम्हारी पाइपलाइन के आउटपुट कैसे और कहाँ प्रकाशित होते हैं, इसे कैसे नियंत्रित करें।

---

## 3. पाइपलाइन आउटपुट का प्रबंधन

एक पाइपलाइन लेखक तय करता है कि आउटपुट कोड में कैसे व्यवस्थित होते हैं, लेकिन यह नियंत्रित करने के लिए कि वे कहाँ जाते हैं या कैसे जाते हैं, तुम्हें उस कोड को छूने की ज़रूरत नहीं है।
Nextflow इसके लिए config-स्तरीय तरीके देता है: एक बेस आउटपुट डायरेक्टरी सेट करो, और चुनो कि फ़ाइलें कॉपी होती हैं या symlink।

### 3.1. आउटपुट डायरेक्टरी को कस्टमाइज़ करना

डिफ़ॉल्ट रूप से, Nextflow आउटपुट `results/` के अंतर्गत प्रकाशित करता है।
इसे `-output-dir` (या इसके संक्षिप्त रूप `-o`) से कहीं और इंगित करो:

```bash
nextflow run main.nf -output-dir outputs
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [serene_kimura] revision: ce74f81996

    executor >  local (8)
    [31/df5c15] sayHello (3)       | 3 of 3 ✔
    [08/8bb2e5] convertToUpper (1) | 3 of 3 ✔
    [e5/5814da] collectGreetings   | 1 of 1 ✔
    [cf/8ab8c6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/outputs

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    outputs/batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
    ```

आउटपुट अब बिल्ट-इन `results/batch/` डिफ़ॉल्ट की जगह `outputs/batch/` के अंतर्गत आते हैं।
पाइपलाइन का अपना कोड अभी भी उस बेस डायरेक्टरी के भीतर की संरचना तय करता है, जैसे `batch/` और `intermediates/` सबडायरेक्टरी; `-output-dir` केवल यह नियंत्रित करता है कि वह संरचना कहाँ से शुरू होती है।

`-output-dir` वास्तव में `outputDir` कॉन्फ़िगरेशन विकल्प के लिए एक command-line शॉर्टकट है, इसलिए यह कहीं भी जा सकता है जहाँ कॉन्फ़िगरेशन जा सकती है: सीधे `nextflow.config` में, किसी प्रोफ़ाइल के अंदर, या `-c` ओवरले फ़ाइल में जैसी तुमने इस भाग में पहले उपयोग की थी।
उदाहरण के लिए, यह snippet command line पर पास करने की बजाय सीधे `nextflow.config` में रखी गई वही सेटिंग दिखाता है:

```groovy title="nextflow.config"
outputDir = 'outputs'
```

इस जैसे कॉन्फ़िगरेशन विकल्प के लिए सभी संभावित स्थानों की पूरी सूची के लिए Nextflow संदर्भ में [Configuration file](https://nextflow.io/docs/latest/config.html) देखो।

### 3.2. आउटपुट कैसे प्रकाशित होते हैं, यह चुनना

डिफ़ॉल्ट रूप से, Nextflow आउटपुट को symlinks के रूप में प्रकाशित करता है जो `work/` के अंतर्गत आउटपुट के स्थानों की ओर इंगित करते हैं, न कि वास्तविक कॉपी:

```console
$ ls -l results/batch/intermediates/Hello-output.txt
lrwxr-xr-x  ...  Hello-output.txt -> /workspaces/training/nextflow-run/work/b7/b8c4d1.../Hello-output.txt
```

पाइपलाइन लेखक workflow कोड में प्रत्येक अलग-अलग प्रोसेस के लिए 'publish mode' को `'copy'` या `'move'` पर सेट कर सकते हैं।
वे आमतौर पर यह पाइपलाइन के अंतिम आउटपुट के लिए करते हैं, जबकि intermediate फ़ाइलों के लिए डिफ़ॉल्ट `'symlink'` व्यवहार छोड़ देते हैं जिन्हें पूरी पाइपलाइन चलने के बाद हटाया जा सकता है।

यह डिस्क पर डेटा को डुप्लिकेट करने से बचाता है, लेकिन इसका मतलब है कि तुम link को तोड़े बिना `work/` के अंतर्गत task डायरेक्टरी नहीं हटा सकते, जिससे `-resume` उपयोग करने की क्षमता खो जाती है।
अगर तुम चाहते हो कि सभी आउटपुट फ़ाइलें ठीक से कॉपी हों, तो अपनी पाइपलाइन कॉन्फ़िगरेशन में [`workflow.output.mode`](https://nextflow.io/docs/latest/reference/config.html#workflow) को `'copy'` पर सेट करो। (`-output-dir` के विपरीत, इसके लिए कोई command-line फ्लैग नहीं है; यह केवल config में है।)

इसे `nextflow.config` में सेट करके देखो:

=== "बाद में"

    ```groovy title="nextflow.config" hl_lines="7"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    workflow.output.mode = 'copy'
    ```

=== "पहले"

    ```groovy title="nextflow.config"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

फिर पाइपलाइन चलाओ, batch का नाम बदलकर ताकि तुम आउटपुट में अंतर देख सको:

```bash
nextflow run main.nf --batch withmode
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [angry_noether] revision: ce74f81996

    executor >  local (8)
    [41/2b478d] sayHello (3)       | 3 of 3 ✔
    [bf/dd2840] convertToUpper (2) | 3 of 3 ✔
    [ea/364e97] collectGreetings   | 1 of 1 ✔
    [10/76fe7b] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - withmode/intermediates/Bonjour-output.txt
        - withmode/intermediates/Hello-output.txt
        - withmode/intermediates/Hola-output.txt

      uppercased:
        - withmode/intermediates/UPPER-Bonjour-output.txt
        - withmode/intermediates/UPPER-Hello-output.txt
        - withmode/intermediates/UPPER-Hola-output.txt

      collected: withmode/intermediates/COLLECTED-withmode-output.txt

      batch_report: withmode/withmode-report.txt

      cowpy_art: withmode/cowpy-COLLECTED-withmode-output.txt
    ```

पहले की तरह आउटपुट फ़ाइलों में से एक देखो:

```console
$ ls -l results/withmode/intermediates/Hello-output.txt
-rw-r--r--  ...  Hello-output.txt
```

अब यह एक वास्तविक, स्वतंत्र फ़ाइल है जो `work/` साफ होने पर भी उपलब्ध रहेगी।

!!! warning "चेतावनी"

    `workflow.output.mode` सेटिंग केवल उन आउटपुट के लिए डिफ़ॉल्ट भरती है जिनका पाइपलाइन कोड में पहले से कोई mode सेट नहीं है।
    यह उस mode को ओवरराइड नहीं कर सकती जिसे लेखक ने hardcode किया है, चाहे तुम इसे कुछ भी सेट करो।

### सारांश

तुम जानते हो कि बेस आउटपुट डायरेक्टरी को कैसे कस्टमाइज़ करें और कॉपी किए गए और symlinked आउटपुट के बीच कैसे चुनें, दोनों बिना पाइपलाइन के कोड को छुए।

### आगे क्या है?

[भाग 3](./03_manage_executions.md) पर जाओ, जहाँ तुम पिछले रन के इतिहास की जांच करना, execution रिपोर्ट बनाना, और पुरानी work डायरेक्टरी साफ करना सीखोगे।

---

## सारांश

इस भाग में तुमने सीखा:

- `nextflow.config` और प्रोफ़ाइल का उपयोग करके पाइपलाइन व्यवहार को कॉन्फ़िगर करना
- रन-विशिष्ट कॉन्फ़िगरेशन फ़ाइल या पैरामीटर फ़ाइल के ज़रिए कॉन्फ़िगरेशन देना
- आउटपुट डायरेक्टरी को कस्टमाइज़ करना और कॉपी किए गए और symlinked आउटपुट के बीच चुनना
