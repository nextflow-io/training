# भाग 3: कॉन्फ़िगरेशन बदलने के लिए profiles का उपयोग करें

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[भाग 1](./01_packaging_and_execution.md) और [भाग 2](./02_resources_and_retries.md) में, तुमने कुछ कॉन्फ़िगरेशन विकल्प जमा किए: सॉफ़्टवेयर पैकेजिंग, execution प्लेटफ़ॉर्म, और resource आवंटन।
व्यवहार में, तुम अक्सर इन विकल्पों के पूरे सेट के बीच स्विच करना चाहोगे, यह इस बात पर निर्भर करता है कि तुम कहाँ चला रहे हो — उदाहरण के लिए, development के लिए laptop और production के लिए HPC cluster।

Nextflow तुम्हें किसी भी संख्या में [profiles](https://nextflow.io/docs/latest/config.html#profiles) सेट करने देता है जो अलग-अलग कॉन्फ़िगरेशन का वर्णन करते हैं, और runtime पर एक single flag के साथ एक (या कई) को चुनने देता है।

तुम पहले से एक का उपयोग कर चुके हो: [Nextflow Run](../nextflow_run/index.md) का `test` profile इनपुट पैरामीटर को एक छोटे, सुपरिभाषित सेट से override करता है।
अब तुम अपने खुद के infrastructure profiles बनाओगे और उन्हें इसके साथ मिलाओगे।

---

## 1. अलग-अलग वातावरण के लिए profiles बनाएं

### 1.1. Profiles सेट करें

`nextflow.config` में दो profiles जोड़ो: एक Docker के साथ सामान्य laptop पर चलाने के लिए, और एक Slurm scheduler और Conda के साथ university HPC cluster के लिए।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="35" hl_lines="10-19"
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
        my_laptop {
            process.executor = 'local'
            docker.enabled = true
        }
        univ_hpc {
            process.executor = 'slurm'
            conda.enabled = true
            process.resourceLimits = [
                memory: 750.GB,
                cpus: 200,
                time: 30.d
            ]
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="35"
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

`univ_hpc` profile resource limits भी सेट करता है, क्योंकि shared HPC infrastructure पर यह आमतौर पर आवश्यक होता है।

### 1.2. Profile के साथ वर्कफ़्लो चलाएं

Runtime पर `-profile` के साथ एक profile चुनो।

```bash
nextflow run main.nf -profile my_laptop
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nice_kimura] revision: c3c85dec78

    executor >  local (8)
    [df/b86b2f] sayHello (2)       | 3 of 3 ✔
    [b2/e2cb9c] convertToUpper (3) | 3 of 3 ✔
    [af/2ecff7] collectGreetings   | 1 of 1 ✔
    [c6/713c96] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Hola-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

!!! warning "चेतावनी"

    `univ_hpc` profile training वातावरण में नहीं चलेगा, क्योंकि वहाँ कोई Slurm scheduler उपलब्ध नहीं है।

अगर तुम्हें अन्य settings मिलती हैं जो हमेशा एक साथ होती हैं, तो उन्हें संबंधित profile में जोड़ो।
तुम किसी भी अन्य combination को group करने के लिए अतिरिक्त profiles भी बना सकते हो।

### 1.3. कई profiles के साथ चलाएं

Profiles परस्पर अनन्य नहीं हैं।
तुम `-profile <profile1>,<profile2>` के साथ एक साथ कई को activate कर सकते हो।
`my_laptop` को उस `test` profile के साथ मिलाओ जो तुम Nextflow Run से पहले से जानते हो।

```bash
nextflow run main.nf -profile my_laptop,test
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_volhard] revision: c3c85dec78

    executor >  local (8)
    [22/81de4f] sayHello (3)       | 3 of 3 ✔
    [93/e16b04] convertToUpper (1) | 3 of 3 ✔
    [ca/3e6a91] collectGreetings   | 1 of 1 ✔
    [16/0cc7e3] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-test-output.txt

      batch_report: full_pipeline/test-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-test-output.txt
    ```

अलग-अलग फ़ाइल नाम सही तरीके से `test` profile से `batch = 'test'` को उठाते हैं (`COLLECTED-test-output.txt`, और इसी तरह)।

अगर तुम ऐसे profiles को मिलाते हो जो एक ही option सेट करते हैं, तो Nextflow उस conflict को उस value से resolve करता है जिसे वह अंत में पढ़ता है, यानी जो फ़ाइल में बाद में आती है।
अगर conflicting settings पूरी तरह से अलग-अलग कॉन्फ़िगरेशन sources से आती हैं, तो मानक [precedence का क्रम](https://www.nextflow.io/docs/latest/config.html) लागू होता है।

### सारांश

तुम जानते हो कि infrastructure-specific कॉन्फ़िगरेशन को bundle करने वाले profiles कैसे define करें, runtime पर `-profile` के साथ एक को कैसे चुनें, एक single run में कई profiles को कैसे मिलाएं, और जब एक से अधिक profile एक ही option सेट करती है तो Nextflow conflicts को कैसे resolve करता है।

### आगे क्या है?

जानो कि कुछ भी चलाने से पहले पूरी तरह से resolved कॉन्फ़िगरेशन को कैसे inspect करें।

---

## 2. Resolved कॉन्फ़िगरेशन inspect करें

तुमने पहले से [Nextflow Run](../nextflow_run/02_configure_pipeline.md) में `nextflow config -profile test` का उपयोग किया है यह जाँचने के लिए कि एक single profile किसमें resolve होती है।
यह कमांड विशेष रूप से तब उपयोगी हो जाती है जब तुम कई profiles को मिला रहे हो: जैसा कि तुमने अभी देखा, जब दो profiles एक ही option सेट करती हैं, तो यह हाथ से पता लगाना मुश्किल हो सकता है कि कौन सी value वास्तव में जीतती है।
`nextflow config` कमांड pipeline चलाए बिना यह सब तुम्हारे लिए resolve कर देता है।

### 2.1. Default कॉन्फ़िगरेशन resolve करें

```bash
nextflow config
```

??? success "कमांड आउटपुट"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'batch'
       character = 'turkey'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
    }
    ```

यह बिल्कुल वही है जो लागू होगा अगर तुम pipeline को बिना किसी extra flags के चलाओ।

### 2.2. Profiles activate करके कॉन्फ़िगरेशन resolve करें

वही profiles जोड़ो जो तुम actual run के लिए उपयोग करोगे।

```bash
nextflow config -profile my_laptop,test
```

??? success "कमांड आउटपुट"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'test'
       character = 'tux'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
       executor = 'local'
    }
    ```

दोनों की तुलना करने से पुष्टि होती है कि क्या बदला: `params.batch`, `params.character`, और `process.executor` सभी `my_laptop,test` profiles को reflect करते हैं।
यह कॉन्फ़िगरेशन की कई layers वाले pipelines के लिए विशेष रूप से मूल्यवान हो जाता है, जहाँ resolved settings को हाथ से निकालना थकाऊ और error-prone होगा।

### सारांश

तुम जानते हो कि कुछ भी चलाने से पहले, profiles के किसी भी combination के लिए पूरी तरह से resolved कॉन्फ़िगरेशन inspect करने के लिए `nextflow config` का उपयोग कैसे करें।

### आगे क्या है?

तुमने Nextflow pipelines को configure करने की मूल बातें cover कर ली हैं।
आगे कहाँ जाना है, इसके लिए [Course summary](next_steps.md) देखो।

---

## सारांश

इस भाग में तुमने सीखा:

- Infrastructure-specific कॉन्फ़िगरेशन को bundle करने वाले profiles define करना
- एक single run में कई profiles को मिलाना, और यह समझना कि उनके बीच conflicts कैसे resolve होते हैं
- Fully resolved कॉन्फ़िगरेशन inspect करने के लिए `nextflow config` का उपयोग करना
