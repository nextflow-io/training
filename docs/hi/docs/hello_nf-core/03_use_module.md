# भाग 3: एक nf-core मॉड्यूल का उपयोग करें

Hello nf-core प्रशिक्षण पाठ्यक्रम के इस तीसरे भाग में, हम आपको दिखाते हैं कि अपनी pipeline में मौजूदा nf-core मॉड्यूल को कैसे खोजें, इंस्टॉल करें, और उपयोग करें।

nf-core के साथ काम करने के महान लाभों में से एक [nf-core/modules](https://github.com/nf-core/modules) रिपॉजिटरी से पूर्व-निर्मित, परीक्षित मॉड्यूल का लाभ उठाने की क्षमता है।
हर process को शुरू से लिखने के बजाय, आप community द्वारा maintain किए गए मॉड्यूल इंस्टॉल और उपयोग कर सकते हैं जो सर्वोत्तम प्रथाओं का पालन करते हैं।

यह कैसे काम करता है यह दिखाने के लिए, हम `core-hello` pipeline में nf-core/modules से `cat/cat` मॉड्यूल के साथ कस्टम `collectGreetings` मॉड्यूल को बदलेंगे।

??? info "इस सेक्शन से कैसे शुरू करें"

    पाठ्यक्रम का यह सेक्शन मानता है कि आपने [भाग 2: nf-core के लिए Hello को फिर से लिखें](./02_rewrite_hello.md) पूरा कर लिया है और आपके पास एक कार्यशील `core-hello` pipeline है।

    यदि आपने भाग 2 पूरा नहीं किया है या इस भाग के लिए नए सिरे से शुरू करना चाहते हैं, तो आप अपने शुरुआती बिंदु के रूप में `core-hello-part2` समाधान का उपयोग कर सकते हैं।
    `hello-nf-core/` डायरेक्टरी के अंदर से यह कमांड चलाएं:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    यह आपको मॉड्यूल जोड़ने के लिए तैयार एक पूरी तरह से कार्यात्मक nf-core pipeline देता है।
    आप निम्नलिखित कमांड चलाकर परीक्षण कर सकते हैं कि यह सफलतापूर्वक चलता है:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. एक उपयुक्त nf-core मॉड्यूल खोजें और इंस्टॉल करें

सबसे पहले, आइए सीखें कि मौजूदा nf-core मॉड्यूल कैसे खोजें और इसे अपनी pipeline में इंस्टॉल करें।

हम `collectGreetings` process को बदलने का लक्ष्य रखेंगे, जो एक फ़ाइल में कई greeting फ़ाइलों को जोड़ने के लिए Unix `cat` कमांड का उपयोग करता है।
फ़ाइलों को जोड़ना एक बहुत ही सामान्य ऑपरेशन है, इसलिए यह तर्कसंगत है कि nf-core में पहले से ही उस उद्देश्य के लिए डिज़ाइन किया गया एक मॉड्यूल हो सकता है।

आइए शुरू करें।

### 1.1. nf-core वेबसाइट पर उपलब्ध मॉड्यूल ब्राउज़ करें

nf-core प्रोजेक्ट [https://nf-co.re/modules](https://nf-co.re/modules) पर मॉड्यूल की एक केंद्रीकृत सूची maintain करता है।

अपने web browser में मॉड्यूल पेज पर जाएं और 'concatenate' खोजने के लिए search bar का उपयोग करें।

![module search results](./img/module-search-results.png)

जैसा कि आप देख सकते हैं, काफी कुछ परिणाम हैं, उनमें से कई बहुत विशिष्ट प्रकार की फ़ाइलों को जोड़ने के लिए डिज़ाइन किए गए मॉड्यूल हैं।
उनमें से, आपको `cat_cat` नामक एक सामान्य-उद्देश्य मॉड्यूल दिखाई देना चाहिए।

!!! note "मॉड्यूल नामकरण परंपरा"

    underscore (`_`) का उपयोग मॉड्यूल नामों में slash (`/`) वर्ण के स्थान पर किया जाता है।

    nf-core मॉड्यूल `software/command` नामकरण परंपरा का पालन करते हैं जब कोई tool कई कमांड प्रदान करता है, जैसे `samtools/view` (samtools पैकेज, view कमांड) या `gatk/haplotypecaller` (GATK पैकेज, HaplotypeCaller कमांड)।
    उन tools के लिए जो केवल एक मुख्य कमांड प्रदान करते हैं, मॉड्यूल `fastqc` या `multiqc` जैसे एकल स्तर का उपयोग करते हैं।

मॉड्यूल दस्तावेज़ देखने के लिए `cat_cat` मॉड्यूल बॉक्स पर क्लिक करें।

मॉड्यूल पेज दिखाता है:

- एक संक्षिप्त विवरण: "A module for concatenation of gzipped or uncompressed files"
- इंस्टॉलेशन कमांड: `nf-core modules install cat/cat`
- Input और output channel संरचना
- उपलब्ध पैरामीटर

### 1.2. कमांड लाइन से उपलब्ध मॉड्यूल की सूची बनाएं

वैकल्पिक रूप से, आप nf-core tools का उपयोग करके सीधे कमांड लाइन से भी मॉड्यूल खोज सकते हैं।

```bash
nf-core modules list remote
```

यह nf-core/modules रिपॉजिटरी में सभी उपलब्ध मॉड्यूल की एक सूची प्रदर्शित करेगा, हालांकि यह थोड़ा कम सुविधाजनक है यदि आप पहले से उस मॉड्यूल का नाम नहीं जानते हैं जिसे आप खोज रहे हैं।
हालाँकि, यदि आप जानते हैं, तो आप विशिष्ट मॉड्यूल खोजने के लिए सूची को `grep` में pipe कर सकते हैं:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "कमांड आउटपुट"

    ```console
    │ cat/cat
    ```

बस ध्यान रखें कि `grep` दृष्टिकोण केवल उन परिणामों को निकालेगा जिनके नाम में खोज शब्द है, जो `cat_cat` के लिए काम नहीं करेगा।

### 1.3. मॉड्यूल के बारे में विस्तृत जानकारी प्राप्त करें

कमांड लाइन से किसी विशिष्ट मॉड्यूल के बारे में विस्तृत जानकारी देखने के लिए, `info` कमांड का उपयोग करें:

```bash
nf-core modules info cat/cat
```

यह मॉड्यूल के बारे में दस्तावेज़ प्रदर्शित करता है, जिसमें इसके inputs, outputs, और बुनियादी उपयोग जानकारी शामिल है।

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

यह वही जानकारी है जो आप वेबसाइट पर पा सकते हैं।

### 1.4. cat/cat मॉड्यूल इंस्टॉल करें

अब जब हमें वह मॉड्यूल मिल गया है जो हम चाहते हैं, तो हमें इसे अपनी pipeline के source code में जोड़ना होगा।

अच्छी खबर यह है कि nf-core प्रोजेक्ट में कुछ tooling शामिल है जो इस भाग को आसान बनाती है।
विशेष रूप से, `nf-core modules install` कमांड code को retrieve करने और इसे एक ही चरण में आपके प्रोजेक्ट के लिए उपलब्ध कराने को automate करना संभव बनाता है।

अपनी pipeline डायरेक्टरी में जाएं और इंस्टॉलेशन कमांड चलाएं:

```bash
cd core-hello
nf-core modules install cat/cat
```

tool पहले आपको repository प्रकार निर्दिष्ट करने के लिए कह सकता है।
(यदि नहीं, तो "अंत में, tool मॉड्यूल इंस्टॉल करने के लिए आगे बढ़ेगा।" पर नीचे जाएं।)

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

यदि ऐसा है, तो default प्रतिक्रिया (`Pipeline`) स्वीकार करने के लिए enter दबाएं और जारी रखें।

tool तब भविष्य में इस prompt से बचने के लिए आपके प्रोजेक्ट के configuration में संशोधन करने की पेशकश करेगा।

??? success "कमांड आउटपुट"

    ```console
        INFO     To avoid this prompt in the future, add the 'repository_type' key to your .nf-core.yml file.
        ? Would you like me to add this config now? [y/n] (y):
    ```

इस सुविधाजनक tooling का लाभ उठाना भी उचित है!
default प्रतिक्रिया (हाँ) स्वीकार करने के लिए enter दबाएं।

अंत में, tool मॉड्यूल इंस्टॉल करने के लिए आगे बढ़ेगा।

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
- आपको अपने workflow में उपयोग करने के लिए सही `include` statement प्रदान करता है

!!! tip

    मॉड्यूल इंस्टॉलेशन कमांड चलाने से पहले हमेशा सुनिश्चित करें कि आपकी वर्तमान working डायरेक्टरी आपके pipeline प्रोजेक्ट की root है।

आइए जांचें कि मॉड्यूल सही तरीके से इंस्टॉल किया गया था:

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

आप स्थानीय रूप से इंस्टॉल किए गए मॉड्यूल को सूचीबद्ध करने के लिए nf-core utility से पूछकर भी इंस्टॉलेशन को सत्यापित कर सकते हैं:

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

यह पुष्टि करता है कि `cat/cat` मॉड्यूल अब आपके प्रोजेक्ट के source code का हिस्सा है।

हालाँकि, वास्तव में नए मॉड्यूल का उपयोग करने के लिए, हमें इसे अपनी pipeline में import करना होगा।

### 1.5. मॉड्यूल imports को अपडेट करें

आइए `workflows/hello.nf` workflow के imports सेक्शन में `collectGreetings` मॉड्यूल के लिए `include` statement को `CAT_CAT` के लिए वाले से बदलें।

याद दिलाने के लिए, मॉड्यूल install tool ने हमें उपयोग करने के लिए सटीक statement दिया:

```groovy title="install कमांड द्वारा निर्मित Import statement"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

ध्यान दें कि nf-core परंपरा मॉड्यूल को import करते समय मॉड्यूल नामों के लिए uppercase का उपयोग करना है।

[core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) खोलें और निम्नलिखित प्रतिस्थापन करें:

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

ध्यान दें कि nf-core मॉड्यूल के लिए path local मॉड्यूल से कैसे भिन्न है:

- **nf-core मॉड्यूल**: `'../modules/nf-core/cat/cat/main'` (`main.nf` को संदर्भित करता है)
- **Local मॉड्यूल**: `'../modules/local/collectGreetings.nf'` (एकल फ़ाइल संदर्भ)

मॉड्यूल अब workflow के लिए उपलब्ध है, इसलिए हमें बस `collectGreetings` की call को `CAT_CAT` का उपयोग करने के लिए swap करना है। सही?

इतनी जल्दी नहीं।

इस बिंदु पर, आप कूदने और code संपादित करना शुरू करने के लिए प्रेरित हो सकते हैं, लेकिन यह ध्यान से जांचने के लिए एक पल लेने लायक है कि नया मॉड्यूल क्या अपेक्षा करता है और यह क्या उत्पन्न करता है।

हम इसे एक अलग सेक्शन के रूप में निपटाने जा रहे हैं क्योंकि इसमें एक नया तंत्र शामिल है जिसे हमने अभी तक कवर नहीं किया है: metadata maps।

!!! note

    आप वैकल्पिक रूप से `collectGreetings.nf` फ़ाइल को हटा सकते हैं:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    हालाँकि, आप इसे local और nf-core मॉड्यूल के बीच अंतर को समझने के लिए एक संदर्भ के रूप में रखना चाह सकते हैं।

### निष्कर्ष

आप जानते हैं कि nf-core मॉड्यूल कैसे खोजें और इसे अपने प्रोजेक्ट के लिए उपलब्ध कैसे कराएं।

### आगे क्या है?

आकलन करें कि एक नया मॉड्यूल क्या आवश्यक है और इसे pipeline में integrate करने के लिए किसी भी महत्वपूर्ण परिवर्तन की पहचान करें।

---

## 2. नए मॉड्यूल की आवश्यकताओं का आकलन करें

विशेष रूप से, हमें मॉड्यूल के **interface** की जांच करने की आवश्यकता है, अर्थात इसकी input और output definitions, और इसकी तुलना उस मॉड्यूल के interface से करनी होगी जिसे हम बदलना चाह रहे हैं।
यह हमें यह निर्धारित करने की अनुमति देगा कि क्या हम नए मॉड्यूल को केवल एक drop-in replacement के रूप में treat कर सकते हैं या क्या हमें wiring के कुछ हिस्सों को adapt करने की आवश्यकता होगी।

आदर्श रूप से यह कुछ ऐसा है जो आपको मॉड्यूल इंस्टॉल करने से _पहले_ भी करना चाहिए, लेकिन अरे, देर से बेहतर कभी नहीं।
(क्या यह मायने रखता है, उन मॉड्यूल से छुटकारा पाने के लिए एक `uninstall` कमांड है जिन्हें आप तय करते हैं कि आप अब नहीं चाहते हैं।)

!!! note

    CAT_CAT process में विभिन्न compression प्रकारों, फ़ाइल extensions इत्यादि की कुछ चतुर handling शामिल है जो सख्ती से हम आपको यहां दिखाने की कोशिश कर रहे हैं उससे संबंधित नहीं हैं, इसलिए हम इसमें से अधिकांश को ignore करेंगे और केवल उन हिस्सों पर ध्यान केंद्रित करेंगे जो महत्वपूर्ण हैं।

### 2.1. दोनों मॉड्यूल के interfaces की तुलना करें

याद दिलाने के लिए, यह है कि हमारे `collectGreetings` मॉड्यूल का interface कैसा दिखता है:

```groovy title="modules/local/collectGreetings.nf (अंश)" linenums="1" hl_lines="6-7 10"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
```

`collectGreetings` मॉड्यूल दो inputs लेता है:

- `input_files` में process करने के लिए एक या अधिक input फ़ाइलें होती हैं;
- `batch_name` एक value है जिसका उपयोग हम output फ़ाइल को एक run-specific नाम देने के लिए करते हैं, जो metadata का एक रूप है।

पूर्ण होने पर, `collectGreetings` एकल file path आउटपुट करता है, जो `outfile` tag के साथ emit किया जाता है।

तुलना में, `cat/cat` मॉड्यूल का interface अधिक जटिल है:

```groovy title="modules/nf-core/cat/cat/main.nf (अंश)" linenums="1" hl_lines="11 14"
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

CAT_CAT मॉड्यूल एक एकल input लेता है, लेकिन वह input एक tuple है जिसमें दो चीजें हैं:

- `meta` एक संरचना है जिसमें metadata होता है, जिसे metamap कहा जाता है;
- `files_in` में process करने के लिए एक या अधिक input फ़ाइलें होती हैं, जो `collectGreetings` के `input_files` के बराबर हैं।

पूर्ण होने पर, CAT_CAT अपने outputs दो भागों में देता है:

- एक और tuple जिसमें metamap और concatenated output फ़ाइल होती है, `file_out` tag के साथ emit की गई;
- एक `versions.yml` फ़ाइल जो उपयोग किए गए software version के बारे में जानकारी capture करती है, `versions` tag के साथ emit की गई।

यह भी ध्यान दें कि default रूप से, output फ़ाइल का नाम metadata का हिस्सा है एक identifier के आधार पर होगा (code यहाँ नहीं दिखाया गया है)।

केवल code को देखते हुए यह बहुत कुछ ट्रैक करने के लिए लग सकता है, इसलिए यहां एक diagram है जो आपको सब कुछ एक साथ कैसे फिट होता है यह visualize करने में मदद करने के लिए है।

<figure class="excalidraw">
--8<-- "docs/hello_nf-core/img/module_comparison.svg"
</figure>

आप देख सकते हैं कि दोनों मॉड्यूल की सामग्री के संदर्भ में समान input आवश्यकताएं हैं (input फ़ाइलों का एक सेट प्लस कुछ metadata) लेकिन उस सामग्री को कैसे package किया जाता है, इसके लिए बहुत अलग अपेक्षाएं हैं।
अभी के लिए versions फ़ाइल को ignore करते हुए, उनका मुख्य output भी समतुल्य है (एक concatenated फ़ाइल), सिवाय CAT_CAT output फ़ाइल के संयोजन में metamap को भी emit करता है।

packaging के अंतर को संभालना काफी आसान होगा, जैसा कि आप थोड़ी देर में देखेंगे।
हालाँकि, metamap भाग को समझने के लिए, हमें आपको कुछ अतिरिक्त संदर्भ से परिचित कराने की आवश्यकता है।

### 2.2. Metamaps को समझना

हमने अभी आपको बताया कि CAT_CAT मॉड्यूल अपने input tuple के हिस्से के रूप में एक metadata map की अपेक्षा करता है।
आइए कुछ मिनट लें और गहराई से देखें कि वह क्या है।

**metadata map**, जिसे अक्सर संक्षेप में **metamap** कहा जाता है, data की units के बारे में जानकारी वाला एक Groovy-style map है।
Nextflow pipelines के संदर्भ में, data की units कुछ भी हो सकती हैं जो आप चाहते हैं: व्यक्तिगत नमूने, नमूनों के batches, या संपूर्ण datasets।

परंपरा के अनुसार, एक nf-core metamap का नाम `meta` है और इसमें आवश्यक field `id` होता है, जिसका उपयोग outputs का नामकरण और data की units को ट्रैक करने के लिए किया जाता है।

उदाहरण के लिए, एक विशिष्ट metadata map इस तरह दिख सकता है:

```groovy title="sample-level metamap का उदाहरण"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

या एक मामले में जहां metadata batch स्तर पर attached है:

```groovy title="batch-level metamap का उदाहरण"
[id: 'batch1', date: '25.10.01']
```

अब आइए इसे `CAT_CAT` process के संदर्भ में रखें, जो input फ़ाइलों को metamap के साथ tuple में package किए जाने की अपेक्षा करता है, और output tuple के हिस्से के रूप में metamap को भी आउटपुट करता है।

```groovy title="modules/nf-core/cat/cat/main.nf (अंश)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

परिणामस्वरूप, data की प्रत्येक unit संबंधित metadata के साथ जुड़ी हुई pipeline के माध्यम से यात्रा करती है।
बाद की processes तब उस metadata को भी आसानी से access कर सकती हैं।

याद है कि कैसे हमने आपको बताया कि `CAT_CAT` द्वारा आउटपुट की गई फ़ाइल का नाम metadata का हिस्सा है एक identifier के आधार पर होगा?
यह संबंधित code है:

```groovy title="modules/nf-core/cat/cat/main.nf (अंश)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

यह मोटे तौर पर इस प्रकार अनुवादित होता है: यदि external task पैरामीटर प्रणाली (`task.ext`) के माध्यम से एक `prefix` प्रदान किया जाता है, तो output फ़ाइल का नाम रखने के लिए उसका उपयोग करें; अन्यथा `${meta.id}` का उपयोग करके एक बनाएं, जो metamap में `id` field से मेल खाता है।

आप इस तरह की सामग्री के साथ इस मॉड्यूल में आने वाली input channel की कल्पना कर सकते हैं:

```groovy title="उदाहरण input channel सामग्री"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

फिर output channel सामग्री इस तरह बाहर आती है:

```groovy title="उदाहरण output channel सामग्री"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

जैसा कि पहले उल्लेख किया गया है, `tuple val(meta), path(files_in)` input setup सभी nf-core मॉड्यूल में उपयोग किया जाने वाला एक मानक pattern है।

उम्मीद है कि आप यह देखना शुरू कर सकते हैं कि यह कितना उपयोगी हो सकता है।
न केवल यह आपको metadata के आधार पर outputs का नाम रखने की अनुमति देता है, बल्कि आप विभिन्न पैरामीटर values लागू करने जैसी चीजें भी कर सकते हैं, और specific operators के संयोजन में, आप data को group, sort या filter भी कर सकते हैं क्योंकि यह pipeline के माध्यम से flows करता है।

!!! note "Metadata के बारे में अधिक जानें"

    Nextflow workflows में metadata के साथ काम करने के लिए एक व्यापक परिचय के लिए, जिसमें samplesheets से metadata को कैसे पढ़ें और processing को customize करने के लिए इसका उपयोग कैसे करें, [workflows में Metadata](../side_quests/metadata) side quest देखें।

### 2.3. किए जाने वाले परिवर्तनों का सारांश दें

जो हमने समीक्षा की है उसके आधार पर, ये वे प्रमुख परिवर्तन हैं जो हमें `cat/cat` मॉड्यूल का उपयोग करने के लिए अपनी pipeline में करने की आवश्यकता है:

- Batch name वाला एक metamap बनाएं;
- Concatenate करने के लिए input फ़ाइलों के सेट (convertToUpper से आ रही) के साथ metamap को एक tuple में package करें;
- `collectGreetings()` से `CAT_CAT` में call switch करें;
- `CAT_CAT` process द्वारा निर्मित tuple से output फ़ाइल को `cowpy` में पास करने से पहले extract करें।

बस इतना ही काफी होना चाहिए! अब जब हमारे पास एक योजना है, तो हम गोता लगाने के लिए तैयार हैं।

### निष्कर्ष

आप जानते हैं कि एक नए मॉड्यूल के input और output interface का आकलन कैसे करें ताकि इसकी आवश्यकताओं की पहचान की जा सके, और आपने सीखा है कि nf-core pipelines द्वारा metamaps का उपयोग कैसे किया जाता है ताकि metadata को data के साथ निकटता से जुड़ा रखा जा सके क्योंकि यह एक pipeline के माध्यम से flows करता है।

### आगे क्या है?

नए मॉड्यूल को workflow में integrate करें।

---

## 3. `hello.nf` workflow में CAT_CAT को integrate करें

अब जब आप metamaps के बारे में सब कुछ जानते हैं (या इस पाठ्यक्रम के उद्देश्यों के लिए पर्याप्त, वैसे भी), तो वास्तव में उन परिवर्तनों को लागू करने का समय आ गया है जिन्हें हमने ऊपर रेखांकित किया है।

स्पष्टता के लिए, हम इसे तोड़ेंगे और प्रत्येक चरण को अलग से कवर करेंगे।

!!! note

    नीचे दिखाए गए सभी परिवर्तन `core-hello/workflows/hello.nf` workflow फ़ाइल में `main` ब्लॉक में workflow logic में किए गए हैं।

### 3.1. एक metadata map बनाएं

सबसे पहले, हमें `CAT_CAT` के लिए एक metadata map बनाने की आवश्यकता है, यह ध्यान में रखते हुए कि nf-core मॉड्यूल को metamap में कम से कम एक `id` field की आवश्यकता होती है।

चूंकि हमें किसी अन्य metadata की आवश्यकता नहीं है, हम इसे सरल रख सकते हैं और इस तरह कुछ उपयोग कर सकते हैं:

```groovy title="Syntax उदाहरण"
def cat_meta = [id: 'test']
```

सिवाय हम `id` value को hardcode नहीं करना चाहते हैं; हम `params.batch` पैरामीटर की value का उपयोग करना चाहते हैं।
तो code बन जाता है:

```groovy title="Syntax उदाहरण"
def cat_meta = [id: params.batch]
```

हाँ, एक बुनियादी metamap बनाना वास्तव में इतना सरल है।

आइए इन पंक्तियों को `convertToUpper` call के बाद जोड़ें, `collectGreetings` call को हटाते हुए:

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // एक अभिवादन emit करें
        sayHello(ch_samplesheet)

        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)

        // batch name के साथ ID के रूप में metadata map बनाएं
        def cat_meta = [ id: params.batch ]

        // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // एक अभिवादन emit करें
        sayHello(ch_samplesheet)

        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)

        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
        cowpy(collectGreetings.out.outfile, params.character)
    ```

यह एक सरल metadata map बनाता है जहां `id` हमारे batch नाम पर सेट है (जो test profile का उपयोग करते समय `test` होगा)।

### 3.2. Metadata tuples के साथ एक channel बनाएं

अगला, फ़ाइलों की channel को metadata और फ़ाइलें युक्त tuples की channel में transform करें:

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // एक अभिवादन emit करें
        sayHello(ch_samplesheet)

        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)

        // batch name के साथ ID के रूप में metadata map बनाएं
        def cat_meta = [ id: params.batch ]

        // tuple format में metadata और files के साथ एक channel बनाएं
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // एक अभिवादन emit करें
        sayHello(ch_samplesheet)

        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)

        // batch name के साथ ID के रूप में metadata map बनाएं
        def cat_meta = [ id: params.batch ]

        // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
        cowpy(collectGreetings.out.outfile, params.character)
    ```

हमने जो पंक्ति जोड़ी है वह दो चीजें प्राप्त करती है:

- `.collect()` `convertToUpper` output से सभी फ़ाइलों को एक एकल list में gather करता है
- `.map { files -> tuple(cat_meta, files) }` `[metadata, files]` के format में एक tuple बनाता है जो `CAT_CAT` अपेक्षा करता है

यह सब कुछ है जो हमें `CAT_CAT` के लिए input tuple को सेट अप करने के लिए करने की आवश्यकता है।

### 3.3. CAT_CAT मॉड्यूल को call करें

अब नए बनाए गए channel पर `CAT_CAT` को call करें:

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // एक अभिवादन emit करें
        sayHello(ch_samplesheet)

        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)

        // batch name के साथ ID के रूप में metadata map बनाएं
        def cat_meta = [ id: params.batch ]

        // tuple format में metadata और files के साथ एक channel बनाएं
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // nf-core cat/cat मॉड्यूल का उपयोग करके फ़ाइलों को concatenate करें
        CAT_CAT(ch_for_cat)

        // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // एक अभिवादन emit करें
        sayHello(ch_samplesheet)

        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)

        // batch name के साथ ID के रूप में metadata map बनाएं
        def cat_meta = [ id: params.batch ]

        // tuple format में metadata और files के साथ एक channel बनाएं
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
        cowpy(collectGreetings.out.outfile, params.character)
    ```

यह इस प्रतिस्थापन का सबसे मुश्किल हिस्सा पूरा करता है, लेकिन हम अभी तक पूरी तरह से नहीं हुए हैं: हमें अभी भी अपडेट करने की आवश्यकता है कि हम concatenated output को `cowpy` process में कैसे पास करते हैं।

### 3.4. `cowpy` के लिए tuple से output file को extract करें

पहले, `collectGreetings` process ने बस एक फ़ाइल उत्पन्न की जिसे हम सीधे `cowpy` में पास कर सकते थे।
हालाँकि, `CAT_CAT` process एक tuple उत्पन्न करता है जिसमें output फ़ाइल के अलावा metamap शामिल है।

चूंकि `cowpy` अभी तक metadata tuples स्वीकार नहीं करता है (हम इसे पाठ्यक्रम के अगले भाग में ठीक करेंगे), हमें `cowpy` को सौंपने से पहले `CAT_CAT` द्वारा निर्मित tuple से output फ़ाइल को extract करने की आवश्यकता है:

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // एक अभिवादन emit करें
        sayHello(ch_samplesheet)

        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)

        // batch name के साथ ID के रूप में metadata map बनाएं
        def cat_meta = [ id: params.batch ]

        // tuple format में metadata और files के साथ एक channel बनाएं
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenate the greetings
        CAT_CAT(ch_for_cat)

        // tuple से file को extract करें क्योंकि cowpy अभी तक metadata का उपयोग नहीं करता है
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // cowpy के साथ अभिवादनों का ASCII art generate करें
        cowpy(ch_for_cowpy, params.character)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // एक अभिवादन emit करें
        sayHello(ch_samplesheet)

        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)

        // batch name के साथ ID के रूप में metadata map बनाएं
        def cat_meta = [ id: params.batch ]

        // tuple format में metadata और files के साथ एक channel बनाएं
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenate the greetings
        CAT_CAT(ch_for_cat)

        // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
        cowpy(collectGreetings.out.outfile, params.character)
    ```

`.map{ meta, file -> file }` operation `CAT_CAT` द्वारा निर्मित `[metadata, file]` tuple से file को एक नए channel, `ch_for_cowpy`, में extract करता है।

फिर यह बस उस अंतिम पंक्ति में `collectGreetings.out.outfile` के बजाय `cowpy` को `ch_for_cowpy` पास करने की बात है।

!!! note

    पाठ्यक्रम के अगले भाग में, हम `cowpy` को सीधे metadata tuples के साथ काम करने के लिए अपडेट करेंगे, इसलिए यह extraction चरण आवश्यक नहीं रहेगा।

### 3.5. Workflow का परीक्षण करें

आइए परीक्षण करें कि workflow नए integrated `cat/cat` मॉड्यूल के साथ काम करता है:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

यह यथोचित रूप से जल्दी चलना चाहिए।

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

ध्यान दें कि `collectGreetings` के बजाय अब process execution list में `CAT_CAT` दिखाई देता है।

और बस! हम अब pipeline में उस चरण के लिए custom prototype-grade code के बजाय एक मजबूत community-curated मॉड्यूल का उपयोग कर रहे हैं।

### निष्कर्ष

अब आप जानते हैं कि कैसे:

- nf-core मॉड्यूल खोजें और इंस्टॉल करें
- एक nf-core मॉड्यूल की आवश्यकताओं का आकलन करें
- nf-core मॉड्यूल के साथ उपयोग के लिए एक सरल metadata map बनाएं
- अपने workflow में एक nf-core मॉड्यूल को integrate करें

### आगे क्या है?

अपने local मॉड्यूल को nf-core परंपराओं का पालन करने के लिए adapt करना सीखें।
हम आपको यह भी दिखाएंगे कि nf-core tooling का उपयोग करके template से नए nf-core मॉड्यूल कैसे बनाएं।
