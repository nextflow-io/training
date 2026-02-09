# भाग 3: एक nf-core मॉड्यूल का उपयोग करें

Hello nf-core प्रशिक्षण कोर्स के इस तीसरे भाग में, हम तुम्हें दिखाएंगे कि अपनी पाइपलाइन में मौजूदा nf-core मॉड्यूल को कैसे खोजें, इंस्टॉल करें और उपयोग करें।

nf-core के साथ काम करने के बड़े फायदों में से एक है [nf-core/modules](https://github.com/nf-core/modules) रिपॉजिटरी से पहले से बने, टेस्ट किए गए मॉड्यूल्स का लाभ उठाना।
हर प्रोसेस को शुरू से लिखने के बजाय, तुम कम्युनिटी द्वारा मेंटेन किए गए मॉड्यूल्स को इंस्टॉल और उपयोग कर सकते हो जो बेस्ट प्रैक्टिसेज़ को फॉलो करते हैं।

यह कैसे काम करता है यह दिखाने के लिए, हम `core-hello` पाइपलाइन में कस्टम `collectGreetings` मॉड्यूल को nf-core/modules से `cat/cat` मॉड्यूल से बदलेंगे।

??? info "इस सेक्शन से कैसे शुरू करें"

    कोर्स का यह सेक्शन मानता है कि तुमने [भाग 2: nf-core के लिए Hello को फिर से लिखें](./02_rewrite_hello.md) पूरा कर लिया है और तुम्हारे पास एक काम करती हुई `core-hello` पाइपलाइन है।

    अगर तुमने भाग 2 पूरा नहीं किया या इस भाग के लिए नए सिरे से शुरू करना चाहते हो, तो तुम `core-hello-part2` सॉल्यूशन को अपने शुरुआती बिंदु के रूप में उपयोग कर सकते हो।
    `hello-nf-core/` डायरेक्टरी के अंदर से यह कमांड चलाओ:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    यह तुम्हें मॉड्यूल्स जोड़ने के लिए तैयार एक पूरी तरह से फंक्शनल nf-core पाइपलाइन देता है।
    तुम यह टेस्ट कर सकते हो कि यह सफलतापूर्वक चलती है निम्नलिखित कमांड चलाकर:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. एक उपयुक्त nf-core मॉड्यूल खोजें और इंस्टॉल करें

पहले, आओ सीखें कि मौजूदा nf-core मॉड्यूल कैसे खोजें और इसे अपनी पाइपलाइन में कैसे इंस्टॉल करें।

हम `collectGreetings` प्रोसेस को बदलने का लक्ष्य रखेंगे, जो कई ग्रीटिंग फ़ाइलों को एक में जोड़ने के लिए Unix `cat` कमांड का उपयोग करती है।
फ़ाइलों को जोड़ना एक बहुत सामान्य ऑपरेशन है, इसलिए यह मानना उचित है कि nf-core में पहले से ही इस उद्देश्य के लिए डिज़ाइन किया गया एक मॉड्यूल हो सकता है।

चलो शुरू करते हैं।

### 1.1. nf-core वेबसाइट पर उपलब्ध मॉड्यूल्स ब्राउज़ करें

nf-core प्रोजेक्ट [https://nf-co.re/modules](https://nf-co.re/modules) पर मॉड्यूल्स का एक केंद्रीकृत कैटलॉग मेंटेन करता है।

अपने वेब ब्राउज़र में मॉड्यूल्स पेज पर जाओ और 'concatenate' खोजने के लिए सर्च बार का उपयोग करो।

![module search results](./img/module-search-results.png)

जैसा कि तुम देख सकते हो, काफी कुछ रिज़ल्ट्स हैं, उनमें से कई मॉड्यूल्स बहुत विशिष्ट प्रकार की फ़ाइलों को जोड़ने के लिए डिज़ाइन किए गए हैं।
उनमें से, तुम्हें `cat_cat` नाम का एक दिखना चाहिए जो सामान्य-उद्देश्य वाला है।

!!! note "मॉड्यूल नामकरण परंपरा"

    अंडरस्कोर (`_`) का उपयोग मॉड्यूल नामों में स्लैश (`/`) कैरेक्टर के स्थान पर किया जाता है।

    nf-core मॉड्यूल्स `software/command` नामकरण परंपरा का पालन करते हैं जब कोई टूल कई कमांड्स प्रदान करता है, जैसे `samtools/view` (samtools पैकेज, view कमांड) या `gatk/haplotypecaller` (GATK पैकेज, HaplotypeCaller कमांड)।
    जो टूल्स केवल एक मुख्य कमांड प्रदान करते हैं, उनके मॉड्यूल्स `fastqc` या `multiqc` जैसे एकल स्तर का उपयोग करते हैं।

मॉड्यूल डॉक्यूमेंटेशन देखने के लिए `cat_cat` मॉड्यूल बॉक्स पर क्लिक करो।

मॉड्यूल पेज दिखाता है:

- एक संक्षिप्त विवरण: "A module for concatenation of gzipped or uncompressed files"
- इंस्टॉलेशन कमांड: `nf-core modules install cat/cat`
- इनपुट और आउटपुट चैनल संरचना
- उपलब्ध पैरामीटर्स

### 1.2. कमांड लाइन से उपलब्ध मॉड्यूल्स की सूची बनाएं

वैकल्पिक रूप से, तुम nf-core टूल्स का उपयोग करके सीधे कमांड लाइन से भी मॉड्यूल्स खोज सकते हो।

```bash
nf-core modules list remote
```

यह nf-core/modules रिपॉजिटरी में सभी उपलब्ध मॉड्यूल्स की सूची प्रदर्शित करेगा, हालांकि यह थोड़ा कम सुविधाजनक है अगर तुम पहले से ही उस मॉड्यूल का नाम नहीं जानते जिसे तुम खोज रहे हो।
हालांकि, अगर तुम जानते हो, तो तुम विशिष्ट मॉड्यूल्स खोजने के लिए सूची को `grep` में पाइप कर सकते हो:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "कमांड आउटपुट"

    ```console
    │ cat/cat
    ```

बस ध्यान रखो कि `grep` दृष्टिकोण केवल उन रिज़ल्ट्स को निकालेगा जिनके नाम में सर्च टर्म है, जो `cat_cat` के लिए काम नहीं करेगा।

### 1.3. मॉड्यूल के बारे में विस्तृत जानकारी प्राप्त करें

कमांड लाइन से किसी विशिष्ट मॉड्यूल के बारे में विस्तृत जानकारी देखने के लिए, `info` कमांड का उपयोग करो:

```bash
nf-core modules info cat/cat
```

यह मॉड्यूल के बारे में डॉक्यूमेंटेशन प्रदर्शित करता है, जिसमें इसके इनपुट्स, आउटपुट्स और बुनियादी उपयोग जानकारी शामिल है।

??? success "कमांड आउटपुट"

    ```console

                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    ╭─ Module: cat/cat  ─────────────────────────────────────────────────╮
    │ 🌐 Repository: https://github.com/nf-core/modules.git              │
    │ 🔧 Tools: cat                                                      │
    │ 📖 Description: A module for concatenation of gzipped or           │
    │ uncompressed files                                                 │
    ╰────────────────────────────────────────────────────────────────────╯
                      ╷                                          ╷
    📥 Inputs        │Description                               │Pattern
    ╺━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
    input[0]         │                                          │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      meta  (map)     │Groovy Map containing sample information  │
                      │e.g. [ id:'test', single_end:false ]      │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      files_in  (file)│List of compressed / uncompressed files   │      *
                      ╵                                          ╵
                          ╷                                 ╷
    📥 Outputs           │Description                      │     Pattern
    ╺━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
    file_out             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      meta  (map)         │Groovy Map containing sample     │
                          │information                      │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      ${prefix}  (file)   │Concatenated file. Will be       │ ${file_out}
                          │gzipped if file_out ends with    │
                          │".gz"                            │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
    versions             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      versions.yml  (file)│File containing software versions│versions.yml
                          ╵                                 ╵

    💻  Installation command: nf-core modules install cat/cat

    ```

यह वही जानकारी है जो तुम वेबसाइट पर पा सकते हो।

### 1.4. cat/cat मॉड्यूल इंस्टॉल करें

अब जब हमने वह मॉड्यूल ढूंढ लिया है जो हम चाहते हैं, हमें इसे अपनी पाइपलाइन के सोर्स कोड में जोड़ना होगा।

अच्छी खबर यह है कि nf-core प्रोजेक्ट में कुछ टूलिंग शामिल है जो इस भाग को आसान बनाती है।
विशेष रूप से, `nf-core modules install` कमांड कोड को रिट्रीव करने और इसे एक ही स्टेप में अपने प्रोजेक्ट के लिए उपलब्ध कराने को स्वचालित करना संभव बनाता है।

अपनी पाइपलाइन डायरेक्टरी में जाओ और इंस्टॉलेशन कमांड चलाओ:

```bash
cd core-hello
nf-core modules install cat/cat
```

टूल पहले तुम्हें रिपॉजिटरी टाइप निर्दिष्ट करने के लिए प्रॉम्प्ट कर सकता है।
(अगर नहीं, तो "अंत में, टूल मॉड्यूल इंस्टॉल करने के लिए आगे बढ़ेगा" तक स्किप करो।)

??? success "कमांड आउटपुट"

    ```console

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.4.1 - https://nf-co.re


    WARNING  'repository_type' not defined in .nf-core.yml
    ? Is this repository a pipeline or a modules repository? (Use arrow keys)
    » Pipeline
      Modules repository
    ```

अगर ऐसा है, तो डिफ़ॉल्ट रिस्पॉन्स (`Pipeline`) को स्वीकार करने के लिए एंटर दबाओ और जारी रखो।

टूल फिर भविष्य में इस प्रॉम्प्ट से बचने के लिए तुम्हारे प्रोजेक्ट की कॉन्फ़िगरेशन को संशोधित करने की पेशकश करेगा।

??? success "कमांड आउटपुट"

    ```console
        INFO     To avoid this prompt in the future, add the 'repository_type' key to your .nf-core.yml file.
        ? Would you like me to add this config now? [y/n] (y):
    ```

इस सुविधाजनक टूलिंग का लाभ उठाना अच्छा रहेगा!
डिफ़ॉल्ट रिस्पॉन्स (yes) को स्वीकार करने के लिए एंटर दबाओ।

अंत में, टूल मॉड्यूल इंस्टॉल करने के लिए आगे बढ़ेगा।

??? success "कमांड आउटपुट"

    ```console
    INFO Config added to '.nf-core.yml'
    INFO Reinstalling modules found in 'modules.json' but missing from directory:
    INFO Installing 'cat/cat'
    INFO Use the following statement to include this module:

        include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
    ```

कमांड स्वचालित रूप से:

- मॉड्यूल फ़ाइलों को `modules/nf-core/cat/cat/` में डाउनलोड करता है
- इंस्टॉल किए गए मॉड्यूल को ट्रैक करने के लिए `modules.json` को अपडेट करता है
- तुम्हें अपने workflow में उपयोग करने के लिए सही `include` स्टेटमेंट प्रदान करता है

!!! tip

    मॉड्यूल इंस्टॉलेशन कमांड चलाने से पहले हमेशा सुनिश्चित करो कि तुम्हारी वर्तमान वर्किंग डायरेक्टरी तुम्हारे पाइपलाइन प्रोजेक्ट की रूट है।

आओ जांचें कि मॉड्यूल सही तरीके से इंस्टॉल हुआ:

```bash
tree -L 4 modules
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    modules
    ├── local
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nf-core
        └── cat
            └── cat
                ├── environment.yml
                ├── main.nf
                ├── meta.yml
                └── tests

    5 directories, 7 files
    ```

तुम स्थानीय रूप से इंस्टॉल किए गए मॉड्यूल्स की सूची बनाने के लिए nf-core यूटिलिटी से पूछकर भी इंस्टॉलेशन को वेरिफाई कर सकते हो:

```bash
nf-core modules list local
```

??? success "कमांड आउटपुट"

    ```console
    INFO     Repository type: pipeline
    INFO     Modules installed in '.':

    ┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
    ┃ Module Name ┃ Repository      ┃ Version SHA ┃ Message                                ┃ Date       ┃
    ┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
    │ cat/cat     │ nf-core/modules │ 41dfa3f     │ update meta.yml of all modules (#8747) │ 2025-07-07 │
    └─────────────┴─────────────────┴─────────────┴────────────────────────────────────────┴────────────┘
    ```

यह पुष्टि करता है कि `cat/cat` मॉड्यूल अब तुम्हारे प्रोजेक्ट के सोर्स कोड का हिस्सा है।

हालांकि, नए मॉड्यूल का वास्तव में उपयोग करने के लिए, हमें इसे अपनी पाइपलाइन में इम्पोर्ट करना होगा।

### 1.5. मॉड्यूल इम्पोर्ट्स को अपडेट करें

आओ `workflows/hello.nf` workflow के इम्पोर्ट्स सेक्शन में `collectGreetings` मॉड्यूल के लिए `include` स्टेटमेंट को `CAT_CAT` के लिए वाले से बदलें।

याद दिलाने के लिए, मॉड्यूल इंस्टॉल टूल ने हमें उपयोग करने के लिए सटीक स्टेटमेंट दिया था:

```groovy title="Import statement produced by install command"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

ध्यान दो कि nf-core परंपरा मॉड्यूल नामों को इम्पोर्ट करते समय अपरकेस का उपयोग करना है।

[core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) खोलो और निम्नलिखित प्रतिस्थापन करो:

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
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    include { cowpy                  } from '../modules/local/cowpy.nf'
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
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

ध्यान दो कि nf-core मॉड्यूल का पाथ लोकल मॉड्यूल्स से कैसे अलग है:

- **nf-core मॉड्यूल**: `'../modules/nf-core/cat/cat/main'` (`main.nf` को रेफरेंस करता है)
- **लोकल मॉड्यूल**: `'../modules/local/collectGreetings.nf'` (सिंगल फ़ाइल रेफरेंस)

मॉड्यूल अब workflow के लिए उपलब्ध है, इसलिए हमें बस `collectGreetings` की कॉल को `CAT_CAT` का उपयोग करने के लिए स्वैप करना है। सही?

इतनी जल्दी नहीं।

इस बिंदु पर, तुम कोड एडिट करना शुरू करने के लिए ललचा सकते हो, लेकिन यह ध्यान से जांचने लायक है कि नया मॉड्यूल क्या अपेक्षा करता है और क्या उत्पन्न करता है।

हम इसे एक अलग सेक्शन के रूप में निपटाने जा रहे हैं क्योंकि इसमें एक नया मैकेनिज्म शामिल है जिसे हमने अभी तक कवर नहीं किया है: मेटाडेटा मैप्स।

!!! note

    तुम वैकल्पिक रूप से `collectGreetings.nf` फ़ाइल को डिलीट कर सकते हो:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    हालांकि, तुम इसे लोकल और nf-core मॉड्यूल्स के बीच अंतर को समझने के लिए एक रेफरेंस के रूप में रखना चाह सकते हो।

### सारांश

तुम जानते हो कि nf-core मॉड्यूल कैसे खोजें और इसे अपने प्रोजेक्ट के लिए उपलब्ध कैसे कराएं।

### आगे क्या है?

आकलन करो कि एक नया मॉड्यूल क्या आवश्यकता रखता है और इसे पाइपलाइन में एकीकृत करने के लिए आवश्यक किसी भी महत्वपूर्ण परिवर्तन की पहचान करो।

---

## 2. नए मॉड्यूल की आवश्यकताओं का आकलन करें

विशेष रूप से, हमें मॉड्यूल के **इंटरफ़ेस** की जांच करनी होगी, यानी इसकी इनपुट और आउटपुट परिभाषाएं, और इसकी तुलना उस मॉड्यूल के इंटरफ़ेस से करनी होगी जिसे हम बदलना चाहते हैं।
यह हमें यह निर्धारित करने की अनुमति देगा कि क्या हम नए मॉड्यूल को सिर्फ एक ड्रॉप-इन रिप्लेसमेंट के रूप में ट्रीट कर सकते हैं या क्या हमें कुछ वायरिंग को अनुकूलित करने की आवश्यकता होगी।

आदर्श रूप से यह कुछ ऐसा है जो तुम्हें मॉड्यूल इंस्टॉल करने से _पहले_ करना चाहिए, लेकिन अरे, देर से बेहतर कभी नहीं।
(जो भी हो, उन मॉड्यूल्स से छुटकारा पाने के लिए एक `uninstall` कमांड है जिन्हें तुम अब नहीं चाहते।)

!!! note

    CAT_CAT प्रोसेस में विभिन्न कम्प्रेशन टाइप्स, फ़ाइल एक्सटेंशन्स आदि की काफी चतुर हैंडलिंग शामिल है जो सख्ती से उस चीज़ के लिए प्रासंगिक नहीं हैं जो हम तुम्हें यहां दिखाने की कोशिश कर रहे हैं, इसलिए हम इसमें से अधिकांश को अनदेखा करेंगे और केवल उन हिस्सों पर ध्यान केंद्रित करेंगे जो महत्वपूर्ण हैं।

### 2.1. दोनों मॉड्यूल्स के इंटरफ़ेस की तुलना करें

याद दिलाने के लिए, हमारे `collectGreetings` मॉड्यूल के इंटरफ़ेस की तरह दिखता है:

```groovy title="modules/local/collectGreetings.nf (excerpt)" linenums="1" hl_lines="6-7 10"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
```

`collectGreetings` मॉड्यूल दो इनपुट्स लेता है:

- `input_files` में प्रोसेस करने के लिए एक या अधिक इनपुट फ़ाइलें होती हैं;
- `batch_name` एक वैल्यू है जिसका उपयोग हम आउटपुट फ़ाइल को रन-विशिष्ट नाम असाइन करने के लिए करते हैं, जो मेटाडेटा का एक रूप है।

पूर्णता पर, `collectGreetings` एक सिंगल फ़ाइल पाथ आउटपुट करता है, जो `outfile` टैग के साथ emit किया जाता है।

तुलना में, `cat/cat` मॉड्यूल का इंटरफ़ेस अधिक जटिल है:

```groovy title="modules/nf-core/cat/cat/main.nf (excerpt)" linenums="1" hl_lines="11 14"
process CAT_CAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"               , emit: versions
```

CAT_CAT मॉड्यूल एक सिंगल इनपुट लेता है, लेकिन वह इनपुट एक टपल है जिसमें दो चीजें होती हैं:

- `meta` एक संरचना है जिसमें मेटाडेटा होता है, जिसे मेटामैप कहा जाता है;
- `files_in` में प्रोसेस करने के लिए एक या अधिक इनपुट फ़ाइलें होती हैं, `collectGreetings` के `input_files` के बराबर।

पूर्णता पर, CAT_CAT अपने आउटपुट्स को दो भागों में डिलीवर करता है:

- एक और टपल जिसमें मेटामैप और कॉन्कैटेनेटेड आउटपुट फ़ाइल होती है, `file_out` टैग के साथ emit किया जाता है;
- एक `versions.yml` फ़ाइल जो उपयोग किए गए सॉफ़्टवेयर वर्जन के बारे में जानकारी कैप्चर करती है, `versions` टैग के साथ emit की जाती है।

यह भी ध्यान दो कि डिफ़ॉल्ट रूप से, आउटपुट फ़ाइल को एक आइडेंटिफ़ायर के आधार पर नाम दिया जाएगा जो मेटाडेटा का हिस्सा है (कोड यहां नहीं दिखाया गया)।

यह सिर्फ कोड को देखकर ट्रैक करने के लिए बहुत कुछ लग सकता है, इसलिए यहां एक डायग्राम है जो तुम्हें यह विज़ुअलाइज़ करने में मदद करेगा कि सब कुछ कैसे फिट होता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/module_comparison.svg"
</figure>

तुम देख सकते हो कि दोनों मॉड्यूल्स की सामग्री के संदर्भ में समान इनपुट आवश्यकताएं हैं (इनपुट फ़ाइलों का एक सेट प्लस कुछ मेटाडेटा) लेकिन उस सामग्री को कैसे पैकेज किया जाता है इसके लिए बहुत अलग अपेक्षाएं हैं।
अभी के लिए versions फ़ाइल को अनदेखा करते हुए, उनका मुख्य आउटपुट भी समतुल्य है (एक कॉन्कैटेनेटेड फ़ाइल), सिवाय इसके कि CAT_CAT आउटपुट फ़ाइल के साथ मेटामैप को भी emit करता है।

पैकेजिंग अंतर से निपटना काफी आसान होगा, जैसा कि तुम थोड़ी देर में देखोगे।
हालांकि, मेटामैप भाग को समझने के लिए, हमें तुम्हें कुछ अतिरिक्त संदर्भ से परिचित कराना होगा।

### 2.2. मेटामैप्स को समझना

हमने अभी तुम्हें बताया कि CAT_CAT मॉड्यूल अपने इनपुट टपल के हिस्से के रूप में एक मेटाडेटा मैप की अपेक्षा करता है।
आओ इस पर करीब से नज़र डालने के लिए कुछ मिनट लें कि वह क्या है।

**मेटाडेटा मैप**, जिसे अक्सर संक्षेप में **मेटामैप** कहा जाता है, एक Groovy-स्टाइल मैप है जिसमें डेटा की इकाइयों के बारे में जानकारी होती है।
Nextflow पाइपलाइन्स के संदर्भ में, डेटा की इकाइयां कुछ भी हो सकती हैं जो तुम चाहते हो: व्यक्तिगत नमूने, नमूनों के बैच, या संपूर्ण डेटासेट।

परंपरा के अनुसार, एक nf-core मेटामैप को `meta` नाम दिया जाता है और इसमें आवश्यक फ़ील्ड `id` होता है, जिसका उपयोग आउटपुट्स को नाम देने और डेटा की इकाइयों को ट्रैक करने के लिए किया जाता है।

उदाहरण के लिए, एक विशिष्ट मेटाडेटा मैप इस तरह दिख सकता है:

```groovy title="Example of sample-level metamap"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

या एक मामले में जहां मेटाडेटा बैच स्तर पर संलग्न है:

```groovy title="Example of batch-level metamap"
[id: 'batch1', date: '25.10.01']
```

अब आओ इसे `CAT_CAT` प्रोसेस के संदर्भ में रखें, जो इनपुट फ़ाइलों को मेटामैप के साथ एक टपल में पैकेज होने की अपेक्षा करता है, और आउटपुट टपल के हिस्से के रूप में मेटामैप को भी आउटपुट करता है।

```groovy title="modules/nf-core/cat/cat/main.nf (excerpt)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

परिणामस्वरूप, डेटा की हर इकाई संबंधित मेटाडेटा के साथ संलग्न होकर पाइपलाइन के माध्यम से यात्रा करती है।
बाद की प्रोसेसेज़ फिर उस मेटाडेटा को भी आसानी से एक्सेस कर सकती हैं।

याद है हमने तुम्हें कैसे बताया था कि `CAT_CAT` द्वारा आउटपुट की गई फ़ाइल को मेटाडेटा का हिस्सा एक आइडेंटिफ़ायर के आधार पर नाम दिया जाएगा?
यह प्रासंगिक कोड है:

```groovy title="modules/nf-core/cat/cat/main.nf (excerpt)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

यह मोटे तौर पर इस प्रकार अनुवादित होता है: अगर एक्सटर्नल टास्क पैरामीटर सिस्टम (`task.ext`) के माध्यम से एक `prefix` प्रदान किया जाता है, तो आउटपुट फ़ाइल को नाम देने के लिए उसका उपयोग करो; अन्यथा `${meta.id}` का उपयोग करके एक बनाओ, जो मेटामैप में `id` फ़ील्ड से मेल खाता है।

तुम इस मॉड्यूल में आने वाले इनपुट चैनल की कल्पना इस तरह की सामग्री के साथ कर सकते हो:

```groovy title="Example input channel contents"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

फिर आउटपुट चैनल सामग्री इस तरह बाहर आती है:

```groovy title="Example output channel contents"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

जैसा कि पहले उल्लेख किया गया है, `tuple val(meta), path(files_in)` इनपुट सेटअप सभी nf-core मॉड्यूल्स में उपयोग किया जाने वाला एक मानक पैटर्न है।

उम्मीद है कि तुम यह देखना शुरू कर सकते हो कि यह कितना उपयोगी हो सकता है।
न केवल यह तुम्हें मेटाडेटा के आधार पर आउटपुट्स को नाम देने की अनुमति देता है, बल्कि तुम विभिन्न पैरामीटर वैल्यूज़ लागू करने के लिए इसका उपयोग जैसी चीजें भी कर सकते हो, और विशिष्ट ऑपरेटर्स के संयोजन में, तुम डेटा को पाइपलाइन के माध्यम से प्रवाहित होते समय ग्रुप, सॉर्ट या फ़िल्टर भी कर सकते हो।

!!! note "मेटाडेटा के बारे में और जानें"

    Nextflow workflows में मेटाडेटा के साथ काम करने के लिए एक व्यापक परिचय के लिए, जिसमें samplesheets से मेटाडेटा कैसे पढ़ें और प्रोसेसिंग को कस्टमाइज़ करने के लिए इसका उपयोग कैसे करें, [Metadata in workflows](../side_quests/metadata) साइड क्वेस्ट देखो।

### 2.3. किए जाने वाले परिवर्तनों का सारांश दें

हमने जो समीक्षा की है उसके आधार पर, ये प्रमुख परिवर्तन हैं जो हमें `cat/cat` मॉड्यूल का उपयोग करने के लिए अपनी पाइपलाइन में करने की आवश्यकता है:

- बैच नाम युक्त एक मेटामैप बनाएं;
- मेटामैप को कॉन्कैटेनेट करने के लिए इनपुट फ़ाइलों के सेट के साथ एक टपल में पैकेज करें (`convertToUpper` से बाहर आ रहा है);
- `collectGreetings()` से `CAT_CAT` में कॉल स्विच करें;
- `CAT_CAT` प्रोसेस द्वारा उत्पादित टपल से आउटपुट फ़ाइल को `cowpy` में पास करने से पहले निकालें।

यह काम करना चाहिए! अब जब हमारे पास एक योजना है, हम गोता लगाने के लिए तैयार हैं।

### सारांश

तुम जानते हो कि एक नए मॉड्यूल के इनपुट और आउटपुट इंटरफ़ेस का आकलन कैसे करें ताकि इसकी आवश्यकताओं की पहचान की जा सके, और तुमने सीखा है कि nf-core पाइपलाइन्स द्वारा मेटामैप्स का उपयोग कैसे किया जाता है ताकि मेटाडेटा को डेटा के साथ निकटता से जुड़ा रखा जा सके जैसे यह पाइपलाइन के माध्यम से प्रवाहित होता है।

### आगे क्या है?

नए मॉड्यूल को एक workflow में एकीकृत करो।

---

## 3. CAT_CAT को `hello.nf` workflow में एकीकृत करें

अब जब तुम मेटामैप्स के बारे में सब कुछ जानते हो (या इस कोर्स के उद्देश्यों के लिए पर्याप्त, वैसे भी), यह वास्तव में उन परिवर्तनों को लागू करने का समय है जिन्हें हमने ऊपर रेखांकित किया था।

स्पष्टता के लिए, हम इसे तोड़ेंगे और प्रत्येक चरण को अलग से कवर करेंगे।

!!! note

    नीचे दिखाए गए सभी परिवर्तन `core-hello/workflows/hello.nf` workflow फ़ाइल में `main` ब्लॉक में workflow लॉजिक में किए गए हैं।

### 3.1. एक मेटाडेटा मैप बनाएं

पहले, हमें `CAT_CAT` के लिए एक मेटाडेटा मैप बनाना होगा, यह ध्यान में रखते हुए कि nf-core मॉड्यूल्स को मेटामैप में कम से कम एक `id` फ़ील्ड की आवश्यकता होती है।

चूंकि हमें किसी अन्य मेटाडेटा की आवश्यकता नहीं है, हम इसे सरल रख सकते हैं और इस तरह कुछ उपयोग कर सकते हैं:

```groovy title="Syntax example"
def cat_meta = [id: 'test']
```

सिवाय इसके कि हम `id` वैल्यू को हार्डकोड नहीं करना चाहते; हम `params.batch` पैरामीटर की वैल्यू का उपयोग करना चाहते हैं।
तो कोड बन जाता है:

```groovy title="Syntax example"
def cat_meta = [id: params.batch]
```

हां, एक बुनियादी मेटामैप बनाना सचमुच इतना सरल है।

आओ `convertToUpper` कॉल के बाद इन लाइनों को जोड़ें, `collectGreetings` कॉल को हटाते हुए:

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

यह एक सरल मेटाडेटा मैप बनाता है जहां `id` हमारे बैच नाम पर सेट है (जो टेस्ट प्रोफ़ाइल का उपयोग करते समय `test` होगा)।

### 3.2. मेटाडेटा टपल्स के साथ एक चैनल बनाएं

अगला, फ़ाइलों के चैनल को मेटाडेटा और फ़ाइलों युक्त टपल्स के चैनल में ट्रांसफ़ॉर्म करो:

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // create a channel with metadata and files in tuple format
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

हमने जो लाइन जोड़ी है वह दो चीजें हासिल करती है:

- `.collect()` `convertToUpper` आउटपुट से सभी फ़ाइलों को एक सिंगल लिस्ट में इकट्ठा करता है
- `.map { files -> tuple(cat_meta, files) }` `[metadata, files]` के फॉर्मेट में एक टपल बनाता है जो `CAT_CAT` अपेक्षा करता है

यह सब कुछ है जो हमें `CAT_CAT` के लिए इनपुट टपल सेट अप करने के लिए करने की आवश्यकता है।

### 3.3. CAT_CAT मॉड्यूल को कॉल करें

अब नए बनाए गए चैनल पर `CAT_CAT` को कॉल करो:

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // create a channel with metadata and files in tuple format
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenate files using the nf-core cat/cat module
        CAT_CAT(ch_for_cat)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // create a channel with metadata and files in tuple format
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

यह इस प्रतिस्थापन के सबसे मुश्किल हिस्से को पूरा करता है, लेकिन हम अभी पूरी तरह से नहीं हुए हैं: हमें अभी भी यह अपडेट करने की आवश्यकता है कि हम कॉन्कैटेनेटेड आउटपुट को `cowpy` प्रोसेस में कैसे पास करते हैं।

### 3.4. `cowpy` के लिए टपल से आउटपुट फ़ाइल निकालें

पहले, `collectGreetings` प्रोसेस ने सिर्फ एक फ़ाइल उत्पन्न की जिसे हम सीधे `cowpy` में पास कर सकते थे।
हालांकि, `CAT_CAT` प्रोसेस एक टपल उत्पन्न करती है जिसमें आउटपुट फ़ाइल के अलावा मेटामैप भी शामिल है।

चूंकि `cowpy` अभी तक मेटाडेटा टपल्स को स्वीकार नहीं करता है (हम इसे कोर्स के अगले भाग में ठीक करेंगे), हमें `CAT_CAT` द्वारा उत्पादित टपल से आउटपुट फ़ाइल को `cowpy` को सौंपने से पहले निकालने की आवश्यकता है:

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // create a channel with metadata and files in tuple format
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenate the greetings
        CAT_CAT(ch_for_cat)

        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generate ASCII art of the greetings with cowpy
        cowpy(ch_for_cowpy, params.character)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // create a channel with metadata and files in tuple format
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenate the greetings
        CAT_CAT(ch_for_cat)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

`.map{ meta, file -> file }` ऑपरेशन `CAT_CAT` द्वारा उत्पादित `[metadata, file]` टपल से फ़ाइल को एक नए चैनल, `ch_for_cowpy` में निकालता है।

फिर यह सिर्फ उस अंतिम लाइन में `collectGreetings.out.outfile` के बजाय `ch_for_cowpy` को `cowpy` में पास करने की बात है।

!!! note

    कोर्स के अगले भाग में, हम `cowpy` को सीधे मेटाडेटा टपल्स के साथ काम करने के लिए अपडेट करेंगे, इसलिए यह एक्सट्रैक्शन स्टेप अब आवश्यक नहीं होगा।

### 3.5. workflow को टेस्ट करें

आओ टेस्ट करें कि workflow नए एकीकृत `cat/cat` मॉड्यूल के साथ काम करता है:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

यह काफी जल्दी चलना चाहिए।

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W ~ version 25.04.3

        Launching `./main.nf` [evil_pike] DSL2 - revision: b9e9b3b8de

        Input/output options
          input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
          outdir                    : core-hello-results

        Institutional config options
          config_profile_name       : Test profile
          config_profile_description: Minimal test dataset to check pipeline function

        Generic options
          validate_params           : false
          trace_report_suffix       : 2025-10-30_18-50-58

        Core Nextflow options
          runName                   : evil_pike
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
        [b3/f005fd] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
        [08/f923d0] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
        [34/3729a9] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
        [24/df918a] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
        -[core/hello] Pipeline completed successfully-
    ```

ध्यान दो कि `CAT_CAT` अब `collectGreetings` के बजाय प्रोसेस एक्ज़ीक्यूशन लिस्ट में दिखाई देता है।

और बस! अब हम पाइपलाइन में उस स्टेप के लिए कस्टम प्रोटोटाइप-ग्रेड कोड के बजाय एक मजबूत कम्युनिटी-क्यूरेटेड मॉड्यूल का उपयोग कर रहे हैं।

### सारांश

अब तुम जानते हो कि कैसे:

- nf-core मॉड्यूल्स खोजें और इंस्टॉल करें
- एक nf-core मॉड्यूल की आवश्यकताओं का आकलन करें
- एक nf-core मॉड्यूल के साथ उपयोग के लिए एक सरल मेटाडेटा मैप बनाएं
- एक nf-core मॉड्यूल को अपने workflow में एकीकृत करें

### आगे क्या है?

अपने लोकल मॉड्यूल्स को nf-core परंपराओं का पालन करने के लिए अनुकूलित करना सीखो।
हम तुम्हें यह भी दिखाएंगे कि nf-core टूलिंग का उपयोग करके टेम्पलेट से नए nf-core मॉड्यूल कैसे बनाएं।
