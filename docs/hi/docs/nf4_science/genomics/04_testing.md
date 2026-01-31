# भाग 4: टेस्ट जोड़ना

इस कोर्स के पहले भाग में, आपने एक variant calling pipeline बनाई जो पूरी तरह से linear थी और प्रत्येक नमूने के डेटा को दूसरों से स्वतंत्र रूप से प्रोसेस करती थी।

दूसरे भाग में, हमने आपको दिखाया कि GATK के साथ joint variant calling को लागू करने के लिए channels और channel operators का उपयोग कैसे करें।

तीसरे भाग में, हमने pipeline को modularize किया।

प्रशिक्षण के इस भाग में, हम आपको दिखाने जा रहे हैं कि [**nf-test**](https://www.nf-test.com/) का उपयोग कैसे करें, एक testing framework जो Nextflow के साथ अच्छी तरह integrate होता है और आपकी pipeline में module-level और workflow-level दोनों प्रकार के टेस्ट जोड़ना सरल बनाता है। प्रशिक्षण के इस भाग का पालन करने के लिए, आपको भाग 1, भाग 2, और भाग 3 के साथ-साथ [nf-test side quest](../../side_quests/nf-test.md) पूरा करना चाहिए, जो nf-test की मूल बातें और testing क्यों महत्वपूर्ण है, इस पर आधारित है।

---

## 0. Warmup

!!! note

    सुनिश्चित करें कि आप सही working डायरेक्टरी में हैं:
    `cd /workspaces/training/nf4-science/genomics`

यदि आपने इस प्रशिक्षण कोर्स के पिछले भागों को पूरा किया है, तो आपके पास उचित modules डायरेक्टरी संरचना के साथ genomics pipeline का एक working संस्करण होना चाहिए।

??? abstract "डायरेक्टरी की सामग्री"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf
    ```

यह modules डायरेक्टरी `solutions` डायरेक्टरी में पाई जा सकती है यदि आपको इसकी आवश्यकता हो।

हम भाग 3 में जैसी ही workflow के साथ शुरू करने जा रहे हैं, जिसे हमने आपके लिए `genomics-4.nf` फ़ाइल में प्रदान किया है। [nf-test side quest](../../side_quests/nf-test.md) की तरह ही, हम इस pipeline में तीन processes के लिए कुछ अलग-अलग प्रकार के टेस्ट जोड़ने जा रहे हैं, साथ ही एक workflow-level टेस्ट भी।

### 0.1. जांचें कि workflow चलती है

टेस्ट जोड़ना शुरू करने से पहले, सुनिश्चित करें कि workflow अपेक्षा के अनुसार चलती है।

```bash
nextflow run genomics-4.nf -resume
```

यदि आप इस प्रशिक्षण कोर्स को शुरू से कर रहे हैं तो यह अब तक बहुत परिचित होना चाहिए।

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-4.nf` [gloomy_poincare] DSL2 - revision: 43203316e0

    executor >  local (7)
    [18/89dfa4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [30/b2522b] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    [a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

पहले की तरह, अब आपकी project डायरेक्टरी के अंदर एक `work` डायरेक्टरी और एक `results_genomics` डायरेक्टरी होगी। हम वास्तव में इन परिणामों का उपयोग बाद में अपनी testing में करेंगे। लेकिन अब से हम pipeline को टेस्ट करने के लिए `nf-test` package का उपयोग करने जा रहे हैं।

### 0.2. `nf-test` को initialize करें

[nf-test side quest](../../side_quests/nf-test.md) की तरह, हमें `nf-test` package को initialize करना होगा।

```bash
nf-test init
```

??? success "कमांड आउटपुट"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

    Project configured. Configuration is stored in nf-test.config
    ```

??? abstract "nf-test.config की सामग्री"

    ```groovy title="nf-test.config"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

यह एक configuration फ़ाइल stub युक्त `tests` डायरेक्टरी भी बनाता है।

### निष्कर्ष

अब हम अपनी genomics pipeline के लिए टेस्ट लिखने के लिए तैयार हैं।

### आगे क्या है?

मूल टेस्ट लिखें जो मूल्यांकन करते हैं कि process calls सफल रहे और सही आउटपुट उत्पन्न हुए।

---

## 1. सफलता और मिलान करते आउटपुट के लिए एक process को टेस्ट करें

हम `SAMTOOLS_INDEX` process को टेस्ट करके शुरू करेंगे, जो कुशल random access को सक्षम करने के लिए BAM फ़ाइलों के लिए इंडेक्स फ़ाइलें बनाती है। यह एक अच्छा पहला टेस्ट केस है क्योंकि:

1. इसमें एक एकल, सुपरिभाषित इनपुट (एक BAM फ़ाइल) है
2. यह एक अनुमानित आउटपुट (एक BAI इंडेक्स फ़ाइल) उत्पन्न करता है
3. समान इनपुट के लिए आउटपुट समान होना चाहिए

### 1.1. एक टेस्ट फ़ाइल stub उत्पन्न करें

सबसे पहले, एक टेस्ट फ़ाइल stub उत्पन्न करें:

```bash
nf-test generate process modules/samtools/index/main.nf
```

??? success "कमांड आउटपुट"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/samtools/index/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/samtools/index/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

यह `main.nf` के समान डायरेक्टरी में एक फ़ाइल बनाता है।
आप file explorer में डायरेक्टरी पर जा सकते हैं और फ़ाइल खोल सकते हैं, जिसमें निम्नलिखित कोड होना चाहिए:

```groovy title="tests/modules/samtools/index/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

शुरुआती assertions [nf-test side quest](../../side_quests/nf-test.md) से परिचित होने चाहिए:

- `assert process.success` कहता है कि हम उम्मीद करते हैं कि process सफलतापूर्वक चले और बिना किसी विफलता के पूरा हो।
- `snapshot(process.out).match()` कहता है कि हम उम्मीद करते हैं कि रन का परिणाम पिछले रन में प्राप्त परिणाम (यदि लागू हो) के समान हो।
  हम इस पर बाद में अधिक विस्तार से चर्चा करते हैं।

इसे एक शुरुआती बिंदु के रूप में उपयोग करते हुए, हमें samtools index process के लिए सही टेस्ट इनपुट जोड़ने होंगे, और यदि लागू हो तो कोई भी पैरामीटर।

### 1.2. टेस्ट फ़ाइल को स्थानांतरित करें और script path को अपडेट करें

टेस्ट को भरने के लिए काम करने से पहले, हमें फ़ाइल को उसके निश्चित स्थान पर ले जाना होगा। हमने प्रत्येक मॉड्यूल के लिए एक डायरेक्टरी जोड़ने का एक कारण यह है कि अब हम प्रत्येक मॉड्यूल की `main.nf` फ़ाइल के साथ co-located `tests` डायरेक्टरी में टेस्ट ship कर सकते हैं। वह डायरेक्टरी बनाएं और टेस्ट फ़ाइल को वहाँ ले जाएं।

```bash
mkdir -p modules/samtools/index/tests
mv tests/modules/samtools/index/main.nf.test modules/samtools/index/tests/
```

अब हम टेस्ट फ़ाइल के `script` सेक्शन को एक relative path में सरल बना सकते हैं:

=== "बाद में"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "../main.nf"
    process "SAMTOOLS_INDEX"
    ```

=== "पहले"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"
    ```

यह टेस्ट को बताता है कि मॉड्यूल की `main.nf` फ़ाइल कहाँ खोजनी है, पूरा path निर्दिष्ट किए बिना।

### 1.3. SAMTOOLS_INDEX के लिए टेस्ट इनपुट प्रदान करें

stub फ़ाइल में एक placeholder शामिल है जिसे हमें `samtools index` के इनपुट के लिए उपयुक्त वास्तविक टेस्ट इनपुट से बदलना होगा। उपयुक्त इनपुट एक BAM फ़ाइल है, जो हमारे पास `data/bam` डायरेक्टरी में उपलब्ध है।

=== "बाद में"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        input[0] = file("${projectDir}/data/bam/reads_son.bam")
        """
    }
    ```

=== "पहले"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        // define inputs of the process here. Example:
        // input[0] = file("test-file.txt")
        """
    }
    ```

### 1.4. कार्यक्षमता के आधार पर टेस्ट का नाम दें

जैसा कि हमने पहले सीखा, टेस्ट का नाम बदलना अच्छा अभ्यास है ताकि यह टेस्ट के संदर्भ में समझ में आए।

=== "बाद में"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should index reads_son.bam correctly") {
    ```

    यह एक arbitrary string लेता है, इसलिए हम जो चाहें डाल सकते हैं।
    यहां हम फ़ाइल नाम और उसके प्रारूप को संदर्भित करना चुनते हैं।

=== "पहले"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should run without failures") {
    ```

### 1.5. टेस्ट चलाएं और आउटपुट की जांच करें

टेस्ट चलाएं:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "कमांड आउटपुट"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.717s)
      Snapshots:
        1 created [Should index reads_son.bam correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 7.727s
    ```

जैसा कि हमने पहले सीखा, इसने process की सफलता के बारे में मूल assertion को सत्यापित किया और process के आउटपुट के आधार पर एक snapshot फ़ाइल बनाई। हम `tests/modules/samtools/index/tests/main.nf.test.snap` फ़ाइल में snapshot फ़ाइल की सामग्री देख सकते हैं:

```json title="modules/samtools/index/tests/main.nf.test.snap" linenums="1"
{
  "Should index reads_son.bam correctly": {
    "content": [
      {
        "0": [
          [
            "reads_son.bam:md5,af5956d9388ba017944bef276b71d809",
            "reads_son.bam.bai:md5,a2ca7b84998218ee77eff14af8eb8ca2"
          ]
        ]
      }
    ],
    "meta": {
      "nf-test": "0.9.3",
      "nextflow": "25.10.2"
    },
    "timestamp": "2026-01-27T15:09:48.394063389"
  }
}
```

हम टेस्ट को फिर से चला सकते हैं और देख सकते हैं कि यह पास हो जाता है, क्योंकि आउटपुट snapshot के समान है:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "कमांड आउटपुट"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.938s)


    SUCCESS: Executed 1 tests in 7.987s
    ```

### 1.6. `SAMTOOLS_INDEX` में अधिक टेस्ट जोड़ें

कभी-कभी यह सुनिश्चित करने के लिए कि हम विभिन्न संभावित मुद्दों के लिए टेस्ट कर रहे हैं, विभिन्न इनपुट फ़ाइलों की एक श्रृंखला को टेस्ट करना उपयोगी होता है। हमारे टेस्ट डेटा से trio में mother और father की BAM फ़ाइलों के लिए टेस्ट जोड़ें।

```groovy
    test("Should index reads_mother.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_mother.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

    test("Should index reads_father.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_father.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

फिर आप टेस्ट को फिर से चला सकते हैं:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "कमांड आउटपुट"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.185s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (6.576s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (6.31s)
      Snapshots:
        2 created [Should index reads_father.bam correctly, Should index reads_mother.bam correctly]


    Snapshot Summary:
      2 created

    SUCCESS: Executed 3 tests in 20.117s
    ```

चेतावनी पर ध्यान दें, जो `--update-snapshot` पैरामीटर के प्रभाव को संदर्भित करती है।

!!! note

    यहां हम टेस्ट डेटा का उपयोग कर रहे हैं जिसे हमने पहले pipeline के वैज्ञानिक आउटपुट का प्रदर्शन करने के लिए उपयोग किया था।
    यदि हम इन टेस्टों को production वातावरण में संचालित करने की योजना बना रहे होते, तो हमने testing उद्देश्यों के लिए छोटे इनपुट उत्पन्न किए होते।

    सामान्य तौर पर, process functionality का मूल्यांकन करने के लिए आवश्यक और पर्याप्त डेटा के सबसे छोटे टुकड़ों का उपयोग करके unit tests को यथासंभव हल्का रखना महत्वपूर्ण है, अन्यथा कुल runtime गंभीरता से जोड़ सकता है।
    एक टेस्ट suite जो नियमित रूप से चलने में बहुत लंबा समय लेता है, एक ऐसा टेस्ट suite है जिसे सुविधा के हित में छोड़े जाने की संभावना है।

### निष्कर्ष

आपने एक genomics process के लिए अपना पहला मॉड्यूल टेस्ट लिखा है, यह सत्यापित करते हुए कि `SAMTOOLS_INDEX` विभिन्न BAM फ़ाइलों के लिए सही तरीके से इंडेक्स फ़ाइलें बनाता है। टेस्ट suite सुनिश्चित करता है कि:

1. Process सफलतापूर्वक चलता है
2. इंडेक्स फ़ाइलें बनाई गई हैं
3. आउटपुट runs में consistent हैं
4. Process सभी नमूना BAM फ़ाइलों के लिए काम करता है

### आगे क्या है?

हमारी genomics workflow में अन्य processes के लिए टेस्ट लिखना सीखें, chained processes को संभालने के लिए setup method का उपयोग करते हुए। हम यह भी मूल्यांकन करेंगे कि आउटपुट, विशेष रूप से हमारी VCF फ़ाइलें, अपेक्षित variant calls शामिल हैं या नहीं।

---

## 2. एक chained process में टेस्ट जोड़ें और सामग्री के लिए टेस्ट करें

`GATK_HAPLOTYPECALLER` को टेस्ट करने के लिए, हमें process को इनपुट के रूप में `SAMTOOLS_INDEX` आउटपुट प्रदान करना होगा। हम `SAMTOOLS_INDEX` को चलाकर, इसके आउटपुट को पुनर्प्राप्त करके, और उन्हें workflow के लिए टेस्ट डेटा के साथ संग्रहीत करके ऐसा कर सकते हैं। वास्तव में यह एक polished pipeline के लिए अनुशंसित दृष्टिकोण है, लेकिन nf-test `setup` method का उपयोग करके एक वैकल्पिक दृष्टिकोण प्रदान करता है।

setup method के साथ, हम टेस्ट setup के हिस्से के रूप में `SAMTOOLS_INDEX` process को ट्रिगर कर सकते हैं, और फिर इसके आउटपुट को `GATK_HAPLOTYPECALLER` के लिए इनपुट के रूप में उपयोग कर सकते हैं। इसकी एक लागत है: हमें हर बार `GATK_HAPLOTYPECALLER` के लिए टेस्ट चलाते समय `SAMTOOLS_INDEX` process चलानी होगी। हालांकि, शायद हम अभी भी workflow विकसित कर रहे हैं और टेस्ट डेटा pre-generate नहीं करना चाहते जिसे हमें बाद में बदलना पड़ सकता है। `SAMTOOLS_INDEX` process भी बहुत तेज़ है, इसलिए शायद इसके आउटपुट को pre-generate और स्टोर करने के लाभ नगण्य हैं। यहां बताया गया है कि setup method कैसे काम करता है।

### 2.1. टेस्ट फ़ाइल उत्पन्न करें और रखें

पहले की तरह, सबसे पहले हम फ़ाइल stub उत्पन्न करते हैं:

```bash
nf-test generate process modules/gatk/haplotypecaller/main.nf
```

??? success "कमांड आउटपुट"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/haplotypecaller/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/haplotypecaller/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

यह निम्नलिखित टेस्ट stub उत्पन्न करता है:

```groovy title="tests/modules/gatk/haplotypecaller/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_HAPLOTYPECALLER"
    script "modules/gatk/haplotypecaller/main.nf"
    process "GATK_HAPLOTYPECALLER"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 2.2. टेस्ट फ़ाइल को स्थानांतरित करें और script path को अपडेट करें

हम मॉड्यूल की `main.nf` फ़ाइल के साथ co-located टेस्ट फ़ाइल के लिए एक डायरेक्टरी बनाते हैं:

```bash
mkdir -p modules/gatk/haplotypecaller/tests
```

और टेस्ट stub फ़ाइल को वहाँ स्थानांतरित करें:

```bash
mv tests/modules/gatk/haplotypecaller/main.nf.test modules/gatk/haplotypecaller/tests/
```

अंत में, script path को अपडेट करना न भूलें:

=== "बाद में"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "../main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

=== "पहले"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "modules/gatk/haplotypecaller/main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

### 2.3. setup method का उपयोग करके इनपुट प्रदान करें

हम `when` block से पहले एक `setup` block डालते हैं, जहां हम हमारी मूल इनपुट फ़ाइलों में से एक पर `SAMTOOLS_INDEX` process का एक रन ट्रिगर कर सकते हैं। साथ ही, पहले की तरह टेस्ट के नाम को कुछ सार्थक में बदलना याद रखें।

=== "बाद में"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7" hl_lines="1-12"
        test("Should call son's haplotype correctly") {

            setup {
                run("SAMTOOLS_INDEX") {
                    script "../../../samtools/index/main.nf"
                    process {
                        """
                        input[0] =  file("${projectDir}/data/bam/reads_son.bam")
                        """
                    }
                }
            }

            when {
    ```

=== "पहले"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7"  hl_lines="1"
    test("Should run without failures") {

        when {
    ```

फिर हम उस process के आउटपुट को `when` block में संदर्भित कर सकते हैं जहां हम टेस्ट इनपुट निर्दिष्ट करते हैं:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="20"
        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }
```

वह परिवर्तन करें और टेस्ट को फिर से चलाएं:

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "कमांड आउटपुट"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Snapshots:
        1 created [Should call son's haplotype correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 40.555s
    ```

यह पहले की तरह एक snapshot फ़ाइल भी उत्पन्न करता है।

### 2.4. फिर से चलाएं और विफलता देखें

दिलचस्प बात यह है कि यदि आप बिल्कुल वही कमांड फिर से चलाते हैं, तो इस बार टेस्ट विफल हो जाएगा।

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? failure "कमांड आउटपुट"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' FAILED (40.123s)

      java.lang.RuntimeException: Different Snapshot:
      [                                                                                           [
          {                                                                                           {
              "0": [                                                                                      "0": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ],                                                                                          ],
              "1": [                                                                                      "1": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "idx": [                                                                                    "idx": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "vcf": [                                                                                    "vcf": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ]                                                                                           ]
          }                                                                                           }
      ]                                                                                           ]

      Nextflow stdout:

      Nextflow stderr:


        Obsolete snapshots can only be checked if all tests of a file are executed successful.


    FAILURE: Executed 1 tests in 40.156s (1 failed)
    ```

त्रुटि संदेश आपको बताता है कि दो runs के लिए snapshots के बीच अंतर थे; विशेष रूप से, VCF फ़ाइलों के लिए md5sum values अलग हैं।

क्यों? एक लंबी कहानी को छोटा करने के लिए, HaplotypeCaller tool में VCF header में एक timestamp शामिल है जो हर बार अलग है (परिभाषा के अनुसार)।
परिणामस्वरूप, हम केवल यह उम्मीद नहीं कर सकते कि फ़ाइलों में समान md5sums होंगे भले ही उनमें variant calls के संदर्भ में समान सामग्री हो।

हम इससे कैसे निपटते हैं?

### 2.5. एक विशिष्ट variant की जांच करने के लिए content assertion method का उपयोग करें

समस्या को हल करने का एक तरीका [एक अलग प्रकार के assertion](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions) का उपयोग करना है।
इस मामले में, हम identity को assert करने के बजाय विशिष्ट सामग्री की जांच करने जा रहे हैं।
अधिक सटीक रूप से, हम tool को VCF फ़ाइल की lines पढ़ने और विशिष्ट lines के अस्तित्व की जांच करने देंगे।

व्यवहार में, हम `then` block में दूसरे assertion को निम्नानुसार बदल देते हैं:

=== "बाद में"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3 4"
            then {
                assert process.success
                assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_son')
                assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3282	GT:DP:GQ:MIN_DP:PL	0/0:25:72:24:0,72,719')
            }
    ```

=== "पहले"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3"
    then {
        assert process.success
        assert snapshot(process.out).match()
    }
    ```

यहां हम VCF आउटपुट फ़ाइल की पूरी सामग्री को पढ़ रहे हैं और content match की खोज कर रहे हैं, जो एक छोटी टेस्ट फ़ाइल पर करना ठीक है, लेकिन आप किसी बड़ी फ़ाइल पर ऐसा नहीं करना चाहेंगे।
आप इसके बजाय विशिष्ट lines पढ़ना चुन सकते हैं।

इस दृष्टिकोण के लिए अधिक सावधानी से चुनने की आवश्यकता है कि हम 'signal' के रूप में किसका उपयोग करना चाहते हैं।
सकारात्मक पक्ष में, इसका उपयोग बड़ी सटीकता के साथ यह टेस्ट करने के लिए किया जा सकता है कि क्या एक analysis tool लगातार 'कठिन' सुविधाओं (जैसे rare variants) की पहचान कर सकता है क्योंकि यह आगे के विकास से गुजरता है।

### 2.6. फिर से चलाएं और सफलता देखें

एक बार जब हमने इस तरह से टेस्ट को संशोधित कर लिया है, तो हम टेस्ट को कई बार चला सकते हैं, और यह लगातार पास होगा।

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "कमांड आउटपुट"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)


    SUCCESS: Executed 1 tests in 40.555s
    ```

### 2.7. अधिक टेस्ट जोड़ें

mother और father नमूनों के लिए समान टेस्ट जोड़ें:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="43"
    test("Should call mother's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_mother.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3278	GT:DP:GQ:MIN_DP:PL	0/0:38:99:37:0,102,1530')
        }
    }

    test("Should call father's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_father.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3281	GT:DP:GQ:MIN_DP:PL	0/0:44:99:42:0,120,1800')
        }
    }
```

### 2.8. टेस्ट कमांड चलाएं

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "कमांड आउटपुट"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (41.47s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (45.556s)


    SUCCESS: Executed 3 tests in 127.586s
    ```

यह pipeline में इस दूसरे चरण के लिए मूल टेस्ट योजना को पूरा करता है। तीसरे और अंतिम module-level टेस्ट पर आगे बढ़ें!

### निष्कर्ष

आपने सीखा कि कैसे:

1. उन processes को टेस्ट करें जो अन्य processes के आउटपुट पर निर्भर हैं
2. VCF आउटपुट फ़ाइलों में विशिष्ट genomic variants सत्यापित करें
3. विशिष्ट सामग्री की जांच करके non-deterministic आउटपुट को संभालें
4. कई नमूनों में variant calling टेस्ट करें

### आगे क्या है?

joint genotyping चरण के लिए pre-generated टेस्ट डेटा का उपयोग करने वाले टेस्ट लिखना सीखें।

---

## 3. Pre-generated टेस्ट डेटा का उपयोग करें

joint genotyping चरण के लिए, हम एक अलग दृष्टिकोण का उपयोग करेंगे - pre-generated टेस्ट डेटा का उपयोग करना। यह अक्सर इसके लिए बेहतर होता है:

1. कई dependencies के साथ जटिल processes
2. Processes जो चलने में लंबा समय लेती हैं
3. Processes जो एक स्थिर, production pipeline का हिस्सा हैं

### 3.1. टेस्ट डेटा उत्पन्न करें

इस खंड की शुरुआत में हमने जो परिणाम उत्पन्न किए थे, उनका निरीक्षण करें:

```bash
tree results_genomics/
```

```console title="परिणाम डायरेक्टरी की सामग्री"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam -> /workspaces/training/nf4-science/genomics/work/42/a3bf19dbfaf1f3672b16a5d5e6a8be/reads_father.bam
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/cf/289c2d264f496d60a69e3e9ba6463e/reads_father.bam.bai
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/work/af/f31a6ade82cc0cf853c4f61c8bc473/reads_mother.bam
    ├── reads_mother.bam.bai -> /workspaces/training/nf4-science/genomics/work/18/89dfa40a3def17e45421e54431a126/reads_mother.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/work/9f/9615dd553d6f13d8bec4f006ac395f/reads_son.bam
    └── reads_son.bam.bai -> /workspaces/training/nf4-science/genomics/work/4d/cb384a97db5687cc9daab002017c7c/reads_son.bam.bai

2 directories, 14 files
```

joint genotyping चरण को haplotype caller चरणों द्वारा उत्पादित VCF फ़ाइलों को indices के साथ इनपुट के रूप में चाहिए। तो चलिए हमारे पास जो परिणाम हैं उन्हें `jointgenotyping` मॉड्यूल की टेस्ट डायरेक्टरी में कॉपी करें।

```bash
mkdir -p modules/gatk/jointgenotyping/tests/inputs/
cp results_genomics/gvcf/*.g.vcf results_genomics/gvcf/*.g.vcf.idx modules/gatk/jointgenotyping/tests/inputs/
```

अब हम इन फ़ाइलों को joint genotyping चरण के लिए हम जो टेस्ट लिखने जा रहे हैं उसके इनपुट के रूप में उपयोग कर सकते हैं।

### 3.2. टेस्ट फ़ाइल stub उत्पन्न करें

पहले की तरह, सबसे पहले हम फ़ाइल stub उत्पन्न करते हैं:

```bash
nf-test generate process modules/gatk/jointgenotyping/main.nf
```

??? success "कमांड आउटपुट"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/jointgenotyping/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/jointgenotyping/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

यह निम्नलिखित टेस्ट stub उत्पन्न करता है:

```groovy title="tests/modules/gatk/jointgenotyping/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 3.3. टेस्ट फ़ाइल को स्थानांतरित करें और script path को अपडेट करें

इस बार हमारे पास पहले से ही मॉड्यूल की `main.nf` फ़ाइल के साथ co-located टेस्ट के लिए एक डायरेक्टरी है, इसलिए हम टेस्ट stub फ़ाइल को वहाँ स्थानांतरित कर सकते हैं:

```bash
mv tests/modules/gatk/jointgenotyping/main.nf.test modules/gatk/jointgenotyping/tests/
```

और script path को अपडेट करना न भूलें:

=== "बाद में"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "../main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

=== "पहले"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

### 3.4. इनपुट प्रदान करें

process input definitions के आधार पर इनपुट भरें और टेस्ट का नाम तदनुसार बदलें:

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="7"
    test("Should call trio's joint genotype correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf")
                ]
                input[1] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf.idx")
                ]
                input[2] = file("${projectDir}/data/ref/intervals.bed")
                input[3] = "family_trio"
                input[4] = file("${projectDir}/data/ref/ref.fasta")
                input[5] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[6] = file("${projectDir}/data/ref/ref.dict")
                """
            }
        }
```

### 3.5. Content assertions का उपयोग करें

joint genotyping चरण का आउटपुट एक और VCF फ़ाइल है, इसलिए हम फिर से content assertion का उपयोग करने जा रहे हैं।

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="25"
    then {
        assert process.success
        assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son')
        assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0')
    }
```

आउटपुट फ़ाइल में एक विशिष्ट variant की सामग्री की जांच करके, यह टेस्ट सत्यापित करता है कि:

1. joint genotyping process सफलतापूर्वक चलता है
2. आउटपुट VCF में सही क्रम में सभी तीन नमूने शामिल हैं
3. एक विशिष्ट variant को सही तरीके से call किया गया है:
   - प्रत्येक नमूने के लिए सटीक genotypes (father के लिए 0/1, mother और son के लिए 1/1)
   - सही read depths और genotype qualities
   - Population-level आंकड़े जैसे allele frequency (AF=0.833)

हमने पूरी फ़ाइल को snapshot नहीं किया है, लेकिन एक विशिष्ट variant की जांच करके, हम आश्वस्त हो सकते हैं कि joint genotyping process अपेक्षा के अनुसार काम कर रही है।

### 3.6. टेस्ट चलाएं

```bash
nf-test test modules/gatk/jointgenotyping/tests/main.nf.test
```

??? success "कमांड आउटपुट"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (53.827s)


    SUCCESS: Executed 1 tests in 53.837s
    ```

टेस्ट पास हो जाता है, यह सत्यापित करते हुए कि हमारी joint genotyping process सही तरीके से:

1. व्यक्तिगत नमूना VCFs को जोड़ती है
2. joint variant calling करती है
3. runs में consistent genotype calls के साथ एक multi-sample VCF उत्पन्न करती है

### निष्कर्ष

आप जानते हैं कि कैसे:

- टेस्ट के लिए इनपुट के रूप में पहले से उत्पन्न परिणामों का उपयोग करें
- Pre-generated टेस्ट डेटा का उपयोग करके टेस्ट लिखें

### आगे क्या है?

पूरी variant calling pipeline end-to-end सही तरीके से काम करती है यह सत्यापित करने के लिए एक workflow-level टेस्ट जोड़ें।

---

## 4. एक workflow-level टेस्ट जोड़ें

अब हम पूरी variant calling pipeline को टेस्ट करेंगे, BAM फ़ाइलों से लेकर joint genotypes तक। यह सत्यापित करता है कि:

1. सभी processes एक साथ सही तरीके से काम करती हैं
2. डेटा चरणों के बीच उचित रूप से flow करता है
3. अंतिम variant calls consistent हैं

### 4.1. workflow टेस्ट उत्पन्न करें

पूरी pipeline के लिए एक टेस्ट फ़ाइल उत्
