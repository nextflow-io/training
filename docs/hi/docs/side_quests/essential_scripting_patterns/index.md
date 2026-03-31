# आवश्यक Nextflow स्क्रिप्टिंग पैटर्न

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow एक प्रोग्रामिंग भाषा है जो Java Virtual Machine पर चलती है। Nextflow [Groovy](http://groovy-lang.org/) पर बनी है और उसका काफी सिंटैक्स शेयर करती है, लेकिन Nextflow सिर्फ "Groovy with extensions" नहीं है -- यह एक स्वतंत्र भाषा है जिसका पूरी तरह से निर्दिष्ट [syntax](https://nextflow.io/docs/latest/reference/syntax.html) और [standard library](https://nextflow.io/docs/latest/reference/stdlib.html) है।

तुम variables, maps, और lists के बेसिक सिंटैक्स से आगे गए बिना भी काफी Nextflow लिख सकते हो। ज़्यादातर Nextflow ट्यूटोरियल workflow orchestration (channels, processes, और data flow) पर फोकस करते हैं, और सिर्फ उसी से तुम काफी आगे जा सकते हो।

हालांकि, जब तुम्हें डेटा मैनिपुलेट करना हो, जटिल फ़ाइल नाम पार्स करने हों, conditional logic लागू करनी हो, या मजबूत प्रोडक्शन वर्कफ़्लो बनाने हों, तो अपने कोड के दो अलग पहलुओं के बारे में सोचना मददगार होता है: **dataflow** (channels, operators, processes, और workflows) और **scripting** (closures, functions, और process scripts के अंदर का कोड)। यह अंतर कुछ हद तक मनमाना है—यह सब Nextflow कोड ही है—लेकिन यह एक उपयोगी mental model देता है जो समझने में मदद करता है कि तुम कब अपनी पाइपलाइन orchestrate कर रहे हो और कब डेटा मैनिपुलेट कर रहे हो। दोनों में महारत हासिल करने से तुम्हारी स्पष्ट, maintainable वर्कफ़्लो लिखने की क्षमता में काफी सुधार होता है।

### सीखने के लक्ष्य

यह side quest तुम्हें बेसिक concepts से प्रोडक्शन-रेडी पैटर्न तक एक hands-on यात्रा पर ले जाता है।
हम एक सरल CSV-reading वर्कफ़्लो को एक sophisticated bioinformatics पाइपलाइन में बदलेंगे, इसे realistic चुनौतियों के ज़रिए step-by-step विकसित करते हुए:

- **सीमाओं को समझना:** dataflow operations और scripting के बीच अंतर करना, और यह समझना कि वे एक साथ कैसे काम करते हैं
- **डेटा मैनिपुलेशन:** शक्तिशाली operators का उपयोग करके maps और collections को extract, transform, और subset करना
- **स्ट्रिंग प्रोसेसिंग:** regex patterns से जटिल फ़ाइल naming schemes पार्स करना और variable interpolation में महारत हासिल करना
- **Reusable functions:** cleaner, अधिक maintainable वर्कफ़्लो के लिए जटिल logic को named functions में extract करना
- **Dynamic logic:** ऐसे processes बनाना जो अलग-अलग input types के अनुसार adapt हों और dynamic resource allocation के लिए closures का उपयोग करना
- **Conditional routing:** metadata characteristics के आधार पर samples को intelligently अलग-अलग processes में route करना
- **Safe operations:** null-safe operators से missing data को gracefully handle करना और clear error messages के साथ inputs validate करना
- **Configuration-based handlers:** logging, notifications, और lifecycle management के लिए workflow event handlers का उपयोग करना

### पूर्वापेक्षाएँ

इस side quest को शुरू करने से पहले, तुम्हें:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष beginner's course पूरा करना चाहिए।
- बेसिक Nextflow concepts और mechanisms (processes, channels, operators, फ़ाइलों के साथ काम करना, meta data) का उपयोग करने में सहज होना चाहिए।
- सामान्य programming constructs (variables, maps, lists) से बुनियादी परिचय होना चाहिए।

यह ट्यूटोरियल programming concepts को जैसे-जैसे हम उनसे मिलते हैं समझाएगा, इसलिए तुम्हें व्यापक programming अनुभव की ज़रूरत नहीं है।
हम fundamental concepts से शुरू करेंगे और advanced patterns तक पहुँचेंगे।

---

## 0. शुरू करना

#### Training codespace खोलो

अगर तुमने अभी तक नहीं किया है, तो [Environment Setup](../envsetup/index.md) में बताए अनुसार training environment खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Project directory में जाओ

चलो उस directory में जाते हैं जहाँ इस ट्यूटोरियल की फ़ाइलें हैं।

```bash
cd side-quests/essential_scripting_patterns
```

#### सामग्री की समीक्षा करो

तुम्हें एक main workflow फ़ाइल और एक `data` directory मिलेगी जिसमें example data फ़ाइलें हैं।

```console title="Directory contents"
.
├── collect.nf
├── data
│   ├── samples.csv
│   └── sequences
│       ├── SAMPLE_001_S1_L001_R1_001.fastq
│       ├── SAMPLE_002_S2_L001_R1_001.fastq
│       └── SAMPLE_003_S3_L001_R1_001.fastq
├── main.nf
├── modules
│   ├── fastp.nf
│   ├── generate_report.nf
│   └── trimgalore.nf
└── nextflow.config
```

हमारे sample CSV में biological samples के बारे में जानकारी है जिन्हें उनकी विशेषताओं के आधार पर अलग-अलग processing की ज़रूरत है:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

हम इस realistic dataset का उपयोग practical programming techniques explore करने के लिए करेंगे जो तुम्हें real bioinformatics वर्कफ़्लो में मिलेंगी।

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### तैयारी की जाँच

क्या तुम dive in करने के लिए तैयार हो?

- [ ] मैं इस course का लक्ष्य और इसकी पूर्वापेक्षाएँ समझता/समझती हूँ
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी working directory सही तरह से सेट की है
<!-- - [ ] I understand the assignment -->

अगर तुम सभी boxes check कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. Dataflow बनाम Scripting: सीमाओं को समझना

### 1.1. क्या है क्या, यह पहचानना

Nextflow वर्कफ़्लो लिखते समय, **dataflow** (डेटा channels और processes के ज़रिए कैसे चलता है) और **scripting** (वह कोड जो डेटा मैनिपुलेट करता है और निर्णय लेता है) के बीच अंतर करना ज़रूरी है। चलो एक वर्कफ़्लो बनाते हैं जो दिखाता है कि वे एक साथ कैसे काम करते हैं।

#### 1.1.1. बेसिक Nextflow Workflow

एक सरल वर्कफ़्लो से शुरू करो जो सिर्फ CSV फ़ाइल पढ़ता है (हमने यह तुम्हारे लिए `main.nf` में पहले से कर दिया है):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

`workflow` block हमारी पाइपलाइन की संरचना define करता है, जबकि `channel.fromPath()` एक file path से channel बनाता है। `.splitCsv()` operator CSV फ़ाइल को process करता है और प्रत्येक row को एक map data structure में convert करता है।

raw CSV डेटा देखने के लिए यह वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    Launching `main.nf` [marvelous_tuckerman] DSL2 - revision: 6113e05c17

    [sample_id:SAMPLE_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/SAMPLE_001_S1_L001_R1_001.fastq, quality_score:38.5]
    [sample_id:SAMPLE_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/SAMPLE_002_S2_L001_R1_001.fastq, quality_score:35.2]
    [sample_id:SAMPLE_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/SAMPLE_003_S3_L001_R1_001.fastq, quality_score:42.1]
    ```

#### 1.1.2. Map Operator जोड़ना

अब हम `.map()` operator का उपयोग करके डेटा transform करने के लिए scripting जोड़ेंगे, जिससे तुम शायद पहले से परिचित हो। यह operator एक 'closure' लेता है जहाँ हम प्रत्येक item को transform करने के लिए कोड लिख सकते हैं।

!!! note "नोट"

    एक **closure** कोड का एक block है जिसे pass around किया जा सकता है और बाद में execute किया जा सकता है। इसे एक ऐसे function के रूप में सोचो जिसे तुम inline define करते हो। Closures curly braces `{ }` से लिखे जाते हैं और parameters ले सकते हैं। ये Nextflow operators के काम करने के तरीके के लिए fundamental हैं और अगर तुम कुछ समय से Nextflow लिख रहे हो, तो शायद तुम पहले से इनका उपयोग कर रहे हो बिना यह जाने!

यहाँ वह map operation कैसी दिखती है:

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

यह हमारा पहला **closure** है - एक anonymous function जिसे तुम argument के रूप में pass कर सकते हो (Python में lambdas या JavaScript में arrow functions की तरह)। Closures Nextflow operators के साथ काम करने के लिए ज़रूरी हैं।

Closure `{ row -> return row }` एक parameter `row` लेता है (कोई भी नाम हो सकता है: `item`, `sample`, आदि)।

जब `.map()` operator प्रत्येक channel item को process करता है, तो वह उस item को तुम्हारे closure में pass करता है। यहाँ, `row` एक बार में एक CSV row रखता है।

यह बदलाव apply करो और वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

तुम पहले जैसा ही output देखोगे, क्योंकि हम बस input को unchanged return कर रहे हैं। यह confirm करता है कि map operator सही तरह से काम कर रहा है। अब चलो डेटा transform करना शुरू करते हैं।

#### 1.1.3. Map Data Structure बनाना

अब हम अपने closure के अंदर **scripting** logic लिखेंगे जो डेटा के प्रत्येक row को transform करे। यहाँ हम individual data items process करते हैं न कि data flow orchestrate करते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // डेटा transformation के लिए scripting
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

`sample_meta` map एक key-value data structure है (Python में dictionaries, JavaScript में objects, या Ruby में hashes की तरह) जो संबंधित जानकारी store करता है: sample ID, organism, tissue type, sequencing depth, और quality score।

हम अपने डेटा को clean करने के लिए `.toLowerCase()` और `.replaceAll()` जैसे string manipulation methods का उपयोग करते हैं, और CSV से string डेटा को उचित numeric types में convert करने के लिए `.toInteger()` और `.toDouble()` जैसे type conversion methods का उपयोग करते हैं।

यह बदलाव apply करो और वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. Conditional Logic जोड़ना

अब चलो और scripting जोड़ते हैं - इस बार data values के आधार पर निर्णय लेने के लिए ternary operator का उपयोग करते हुए।

निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="11-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2" hl_lines="11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

Ternary operator if/else statement का एक shorthand है जो `condition ? value_if_true : value_if_false` pattern follow करता है। इस line का मतलब है: "अगर quality 40 से ज़्यादा है, तो 'high' use करो, नहीं तो 'normal' use करो"। इसका cousin, **Elvis operator** (`?:`), तब default values देता है जब कुछ null या empty हो - हम उस pattern को इस ट्यूटोरियल में बाद में explore करेंगे।

Map addition operator `+` existing map को modify करने की बजाय एक **नया map** बनाता है। यह line एक नया map बनाती है जिसमें `sample_meta` के सभी key-value pairs और नई `priority` key होती है।

!!! Note "नोट"

    Closures में pass किए गए maps को कभी modify मत करो - हमेशा `+` (उदाहरण के लिए) का उपयोग करके नए बनाओ। Nextflow में, एक ही डेटा अक्सर एक साथ कई operations से होकर गुज़रता है। किसी map को in-place modify करने से unpredictable side effects हो सकते हैं जब दूसरे operations उसी object को reference करते हैं। नए maps बनाने से यह सुनिश्चित होता है कि प्रत्येक operation का अपना clean copy हो।

Modified वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

हमने quality scores के आधार पर priority level के साथ अपने metadata को enrich करने के लिए conditional logic सफलतापूर्वक जोड़ी है।

#### 1.1.5. `.subMap()` से Maps को Subset करना

जहाँ `+` operator एक map में keys जोड़ता है, वहीं कभी-कभी तुम्हें इसका उल्टा करना होता है - केवल specific keys extract करना। `.subMap()` method इसके लिए perfect है।

चलो एक line जोड़ते हैं जो हमारे metadata का एक simplified version बनाए जिसमें केवल identification fields हों:

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // डेटा transformation के लिए scripting
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "केवल ID fields: ${id_only}"

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // डेटा transformation के लिए scripting
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Modified वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    ID fields only: [id:sample_001, organism:human, tissue:liver]
    ID fields only: [id:sample_002, organism:mouse, tissue:brain]
    ID fields only: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

यह `view()` operation द्वारा display किया गया पूरा metadata और `println` से print किया गया extracted subset दोनों दिखाता है।

`.subMap()` method keys की एक list लेती है और केवल उन keys वाला एक नया map return करती है। अगर original map में कोई key नहीं है, तो वह result में simply शामिल नहीं होती।

यह विशेष रूप से उपयोगी है जब तुम्हें अलग-अलग processes के लिए अलग-अलग metadata versions बनाने की ज़रूरत हो - कुछ को पूरा metadata चाहिए हो सकता है जबकि दूसरों को केवल minimal identification fields।

अब उन println statements को हटाओ ताकि तुम्हारा वर्कफ़्लो अपनी पिछली स्थिति में वापस आ जाए, क्योंकि हमें आगे उनकी ज़रूरत नहीं है।

!!! tip "सुझाव: Map Operations का सारांश"

    - **Keys जोड़ना**: `map1 + [new_key: value]` - अतिरिक्त keys के साथ नया map बनाता है
    - **Keys extract करना**: `map1.subMap(['key1', 'key2'])` - केवल specified keys के साथ नया map बनाता है
    - **दोनों operations नए maps बनाते हैं** - Original maps unchanged रहते हैं

#### 1.1.6. Maps को Combine करना और Results Return करना

अब तक, हम केवल वही return कर रहे थे जिसे Nextflow community 'meta map' कहती है, और हम उन फ़ाइलों को ignore कर रहे थे जिनसे वह metadata संबंधित है। लेकिन अगर तुम Nextflow वर्कफ़्लो लिख रहे हो, तो शायद तुम उन फ़ाइलों के साथ कुछ करना चाहते हो।

चलो 2 elements का एक tuple वाला channel structure output करते हैं: enriched metadata map और corresponding file path। यह Nextflow में processes को डेटा pass करने का एक common pattern है।

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple( sample_meta + [priority: priority], file(row.file_path) )
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

यह बदलाव apply करो और वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

यह `[meta, file]` tuple structure Nextflow में processes को metadata और associated फ़ाइलें दोनों pass करने का एक common pattern है।

!!! note "नोट"

    **Maps और Metadata**: Maps Nextflow में metadata के साथ काम करने के लिए fundamental हैं। metadata maps के साथ काम करने की अधिक विस्तृत व्याख्या के लिए, [Working with metadata](../metadata/) side quest देखो।

हमारा वर्कफ़्लो core pattern demonstrate करता है: **dataflow operations** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orchestrate करते हैं कि डेटा पाइपलाइन से कैसे गुज़रता है, जबकि `.map()` closure के अंदर **scripting** (maps `[key: value]`, string methods, type conversions, ternary operators) individual data items के transformation को handle करती है।

### 1.2. अलग-अलग Types को समझना: Channel बनाम List

अब तक, ठीक है, हम dataflow operations और scripting के बीच अंतर कर सकते हैं। लेकिन जब एक ही method name दोनों contexts में exist करे तो क्या होगा?

एक perfect example है `collect` method, जो Nextflow standard library में channel types और List types दोनों के लिए exist करती है। List पर `collect()` method प्रत्येक element को transform करती है, जबकि channel पर `collect()` operator सभी channel emissions को एक single-item channel में gather करता है।

चलो इसे कुछ sample data के साथ demonstrate करते हैं, पहले यह refresh करते हुए कि channel `collect()` operator क्या करता है। `collect.nf` देखो:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - कई channel emissions को एक में group करता है
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

Steps:

- Sample IDs की एक List define करो
- `fromList()` से एक channel बनाओ जो प्रत्येक sample ID को अलग-अलग emit करे
- `view()` से प्रत्येक item को print करो जैसे वह flow करता है
- channel के `collect()` operator से सभी items को एक single list में gather करो
- दूसरे `view()` से collected result (सभी sample IDs वाला single item) print करो

हमने channel की structure बदली है, लेकिन डेटा खुद नहीं बदला।

इसे confirm करने के लिए वर्कफ़्लो चलाओ:

```bash
nextflow run collect.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

`view()` हर channel emission के लिए एक output return करता है, इसलिए हम जानते हैं कि इस single output में सभी 3 original items एक list में grouped हैं।

अब List पर `collect` method को action में देखते हैं। `collect.nf` को modify करो ताकि original list of sample IDs पर List का `collect` method apply हो:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - कई channel emissions को एक में group करता है
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - प्रत्येक element को transform करता है, structure preserve करता है
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - कई channel emissions को एक में group करता है
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

इस नए snippet में हम:

- एक नया variable `formatted_ids` define करते हैं जो original list में प्रत्येक sample ID को transform करने के लिए List के `collect` method का उपयोग करता है
- `println` से result print करते हैं

Modified वर्कफ़्लो चलाओ:

```bash
nextflow run collect.nf
```

??? success "कमांड आउटपुट"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

इस बार, हमने डेटा की structure नहीं बदली, list में अभी भी 3 items हैं, लेकिन हमने List के `collect` method का उपयोग करके प्रत्येक item को transform किया है जिससे modified values वाली एक नई list बनी। यह channel पर `map` operator का उपयोग करने जैसा है, लेकिन यह channel की बजाय List data structure पर operate कर रहा है।

`collect` एक extreme case है जिसे हम यहाँ एक बात समझाने के लिए use कर रहे हैं। मुख्य सबक यह है कि जब तुम वर्कफ़्लो लिख रहे हो, हमेशा **data structures** (Lists, Maps, आदि) और **channels** (dataflow constructs) के बीच अंतर करो। Operations के नाम एक जैसे हो सकते हैं लेकिन वे जिस type पर call किए जाते हैं उसके आधार पर बिल्कुल अलग तरह से behave करते हैं।

### 1.3. Spread Operator (`*.`) - Property Extraction का Shorthand

List के `collect` method से संबंधित spread operator (`*.`) है, जो collections से properties extract करने का एक concise तरीका देता है। यह essentially एक common `collect` pattern के लिए syntactic sugar है।

चलो अपनी `collect.nf` फ़ाइल में एक demonstration जोड़ते हैं:

=== "बाद में"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - कई channel emissions को एक में group करता है
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - प्रत्येक element को transform करता है, structure preserve करता है
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // Spread operator - concise property access
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operator result: ${all_ids}"
    ```

=== "पहले"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - कई channel emissions को एक में group करता है
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - प्रत्येक element को transform करता है, structure preserve करता है
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

Updated वर्कफ़्लो चलाओ:

```bash title="Spread operator test करो"
nextflow run collect.nf
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Spread operator result: [s1, s2, s3]
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

Spread operator `*.` एक common collect pattern का shorthand है:

```groovy
// ये equivalent हैं:
def ids = samples*.id
def ids = samples.collect { it.id }

// Method calls के साथ भी काम करता है:
def names = files*.getName()
def names = files.collect { it.getName() }
```

Spread operator विशेष रूप से उपयोगी है जब तुम्हें objects की list से एक single property extract करनी हो - यह पूरा `collect` closure लिखने से ज़्यादा readable है।

!!! tip "सुझाव: Spread बनाम Collect कब Use करें"

    - **Spread (`*.`) use करो** simple property access के लिए: `samples*.id`, `files*.name`
    - **Collect use करो** transformations या complex logic के लिए: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### सारांश

इस section में, तुमने सीखा:

- **Dataflow बनाम scripting**: Channel operators orchestrate करते हैं कि डेटा तुम्हारी पाइपलाइन से कैसे flow करता है, जबकि scripting individual data items को transform करती है
- **Types को समझना**: एक ही method name (जैसे `collect`) उस type के आधार पर अलग तरह से behave कर सकता है जिस पर इसे call किया जाता है (Channel बनाम List)
- **Context मायने रखता है**: हमेशा aware रहो कि तुम channels (dataflow) या data structures (scripting) के साथ काम कर रहे हो

इन सीमाओं को समझना debugging, documentation, और maintainable वर्कफ़्लो लिखने के लिए ज़रूरी है।

अगला हम string processing capabilities में deeper dive करेंगे, जो real-world डेटा handle करने के लिए ज़रूरी हैं।

---

## 2. String Processing और Dynamic Script Generation

String processing में महारत हासिल करना brittle वर्कफ़्लो को robust पाइपलाइन से अलग करता है। यह section जटिल फ़ाइल नाम parse करना, dynamic script generation, और variable interpolation cover करता है।

### 2.1. Pattern Matching और Regular Expressions

Bioinformatics फ़ाइलों में अक्सर metadata encode करने वाले complex naming conventions होते हैं। चलो regular expressions के साथ pattern matching का उपयोग करके इसे automatically extract करते हैं।

हम अपने `main.nf` वर्कफ़्लो पर वापस जाएंगे और फ़ाइल नामों से additional sample information extract करने के लिए कुछ pattern matching logic जोड़ेंगे। हमारे dataset में FASTQ फ़ाइलें Illumina-style naming conventions follow करती हैं जैसे `SAMPLE_001_S1_L001_R1_001.fastq.gz`। ये cryptic लग सकती हैं, लेकिन ये actually sample ID, lane number, और read direction जैसे useful metadata encode करती हैं। हम इन नामों को parse करने के लिए regex capabilities का उपयोग करेंगे।

अपने existing `main.nf` वर्कफ़्लो में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // डेटा transformation के लिए scripting
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="4" hl_lines="10-11"
            .map { row ->
                // डेटा transformation के लिए scripting
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + [priority: priority], file(row.file_path))
            }
    ```

यह key **string processing concepts** demonstrate करता है:

1. **Regular expression literals** `~/pattern/` syntax का उपयोग करके - यह backslashes escape करने की ज़रूरत के बिना एक regex pattern बनाता है
2. **Pattern matching** `=~` operator के साथ - यह एक string को regex pattern से match करने की कोशिश करता है
3. **Matcher objects** जो `[0][1]`, `[0][2]`, आदि से groups capture करते हैं - `[0]` पूरे match को refer करता है, `[1]`, `[2]`, आदि parentheses में captured groups को refer करते हैं

चलो regex pattern `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$` को break down करते हैं:

| Pattern             | Matches करता है                        | Captures करता है                        |
| ------------------- | -------------------------------------- | --------------------------------------- |
| `^(.+)`             | शुरू से sample name                    | Group 1: sample name                    |
| `_S(\d+)`           | Sample number `_S1`, `_S2`, आदि        | Group 2: sample number                  |
| `_L(\d{3})`         | Lane number `_L001`                    | Group 3: lane (3 digits)                |
| `_(R[12])`          | Read direction `_R1` या `_R2`          | Group 4: read direction                 |
| `_(\d{3})`          | Chunk number `_001`                    | Group 5: chunk (3 digits)               |
| `\.fastq(?:\.gz)?$` | File extension `.fastq` या `.fastq.gz` | Capture नहीं होता (?: non-capturing है) |

यह metadata automatically extract करने के लिए Illumina-style naming conventions parse करता है।

Modified वर्कफ़्लो चलाओ:

```bash title="Pattern matching test करो"
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

यह फ़ाइल नामों से enriched metadata दिखाता है।

### 2.2. Processes में Dynamic Script Generation

Process script blocks essentially multi-line strings हैं जो shell को pass किए जाते हैं। तुम input characteristics के आधार पर अलग-अलग script strings dynamically generate करने के लिए **conditional logic** (if/else, ternary operators) का उपयोग कर सकते हो। यह process definitions duplicate किए बिना diverse input types—जैसे single-end बनाम paired-end sequencing reads—handle करने के लिए ज़रूरी है।

चलो अपने वर्कफ़्लो में एक process जोड़ते हैं जो यह pattern demonstrate करे। `modules/fastp.nf` खोलो और देखो:

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz \\
        --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread $task.cpus
    """
}
```

यह process FASTQ फ़ाइलें input के रूप में लेता है और adapters trim करने और low-quality reads filter करने के लिए `fastp` tool चलाता है। दुर्भाग्य से, जिस व्यक्ति ने यह process लिखा उसने हमारे example dataset में single-end reads के लिए allow नहीं किया। चलो इसे अपने वर्कफ़्लो में जोड़ते हैं और देखते हैं क्या होता है:

पहले, अपने `main.nf` वर्कफ़्लो की पहली line पर module include करो:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

फिर `workflow` block को modify करो ताकि `ch_samples` channel को `FASTP` process से connect किया जा सके:

=== "बाद में"

    ```groovy title="main.nf" linenums="25" hl_lines="27"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="25" hl_lines="26"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + file_meta + [priority: priority], file(row.file_path)]
            }
            .view()
    }
    ```

यह modified वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? failure "कमांड आउटपुट"

    ```console
    ERROR ~ Error executing process > 'FASTP (3)'

    Caused by:
      Process `FASTP (3)` terminated with an error exit status (255)


    Command executed:

      fastp \
          --in1 SAMPLE_003_S3_L001_R1_001.fastq \
          --in2 null \
          --out1 sample_003_trimmed_R1.fastq.gz \
          --out2 sample_003_trimmed_R2.fastq.gz \
          --json sample_003.fastp.json \
          --html sample_003.fastp.html \
          --thread 2

    Command exit status:
      255

    Command output:
      (empty)
    ```

तुम देख सकते हो कि process दूसरी input फ़ाइल के लिए `null` value के साथ `fastp` चलाने की कोशिश कर रहा है, जिससे यह fail हो रहा है। यह इसलिए है क्योंकि हमारे dataset में single-end reads हैं, लेकिन process hardcoded है कि paired-end reads (एक बार में दो input फ़ाइलें) expect करे।

इसे fix करने के लिए `FASTP` process के `script:` block में conditional logic जोड़ो। एक if/else statement read file count check करता है और command accordingly adjust करता है।

=== "बाद में"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Simple single-end बनाम paired-end detection
        def is_single = reads instanceof List ? reads.size() == 1 : true

        if (is_single) {
            def input_file = reads instanceof List ? reads[0] : reads
            """
            fastp \\
                --in1 ${input_file} \\
                --out1 ${meta.id}_trimmed.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        } else {
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="10" hl_lines="2-11"
            script:
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

अब वर्कफ़्लो single-end और paired-end reads दोनों को gracefully handle कर सकता है। Conditional logic input फ़ाइलों की संख्या check करता है और `fastp` के लिए appropriate command construct करता है। चलो देखते हैं कि यह काम करता है या नहीं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_rosalind] DSL2 - revision: 04b1cd93e9

    executor >  local (3)
    [31/a8ad4d] process > FASTP (3) [100%] 3 of 3 ✔
    ```

बढ़िया! अगर हम actual commands check करें जो run हुए (अपने task hash के लिए customize करो):

```console title="Execute किए गए commands check करो"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

हम देख सकते हैं कि Nextflow ने single-end reads के लिए सही command चुना:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Dynamic script logic का एक और common usage [Nextflow for Science Genomics module](../../nf4science/genomics/02_joint_calling) में देखा जा सकता है। उस module में, call किया जाने वाला GATK process कई input फ़ाइलें ले सकता है, लेकिन प्रत्येक को एक correct command line बनाने के लिए `-V` से prefix करना होगा। Process input फ़ाइलों के collection (`all_gvcfs`) को correct command arguments में transform करने के लिए scripting का उपयोग करता है:

```groovy title="GATK के लिए command line manipulation" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Process script blocks में scripting use करने के ये patterns बेहद powerful हैं और कई scenarios में apply किए जा सकते हैं - variable input types handle करने से लेकर file collections से complex command-line arguments बनाने तक, जो तुम्हारे processes को real-world डेटा की diverse requirements के लिए truly adaptable बनाते हैं।

### 2.3. Variable Interpolation: Nextflow और Shell Variables

Process scripts Nextflow variables, shell variables, और command substitutions को mix करते हैं, प्रत्येक के अलग-अलग interpolation syntax के साथ। गलत syntax use करने से errors होती हैं। चलो इन्हें एक process के साथ explore करते हैं जो एक processing report बनाता है।

Module फ़ाइल `modules/generate_report.nf` देखो:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

यह process sample ID और filename के साथ एक simple report लिखता है। अब चलो इसे run करते हैं और देखते हैं कि जब हमें अलग-अलग types के variables mix करने की ज़रूरत हो तो क्या होता है।

Process को अपने `main.nf` में include करो और वर्कफ़्लो में जोड़ो:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="2 30"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1" hl_lines="1 10-29"
    include { FASTP } from './modules/fastp.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

अब वर्कफ़्लो चलाओ और `results/reports/` में generated reports check करो। उनमें प्रत्येक sample के बारे में basic information होनी चाहिए।

<!-- TODO: add the run command -->

??? success "कमांड आउटपुट"

    ```console
    <!-- TODO: output -->
    ```

लेकिन क्या होगा अगर हम processing कब और कहाँ हुई इसके बारे में जानकारी जोड़ना चाहें? चलो process को modify करते हैं ताकि report में current user, hostname, और date शामिल करने के लिए **shell** variables और थोड़ी command substitution का उपयोग किया जा सके:

=== "बाद में"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "पहले"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

अगर तुम यह run करते हो, तो तुम्हें एक error दिखेगी - Nextflow `${USER}` को एक Nextflow variable के रूप में interpret करने की कोशिश करता है जो exist नहीं करता।

??? failure "कमांड आउटपुट"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

हमें इसे escape करना होगा ताकि Bash इसे handle कर सके।

Shell variables और command substitutions को backslash (`\`) से escape करके इसे fix करो:

=== "बाद में"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: \${USER}" >> ${meta.id}_report.txt
        echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
        echo "Date: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "पहले"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

अब यह काम करता है! Backslash (`\`) Nextflow को बताता है "इसे interpret मत करो, इसे Bash को pass करो।"

### सारांश

इस section में, तुमने **string processing** techniques सीखी:

- **फ़ाइल parsing के लिए regular expressions**: जटिल फ़ाइल naming conventions से metadata extract करने के लिए `=~` operator और regex patterns (`~/pattern/`) का उपयोग करना
- **Dynamic script generation**: Input characteristics के आधार पर अलग-अलग script strings generate करने के लिए conditional logic (if/else, ternary operators) का उपयोग करना
- **Variable interpolation**: यह समझना कि Nextflow strings कब interpret करता है बनाम shell कब करता है
  - `${var}` - Nextflow variables (workflow compile time पर Nextflow द्वारा interpolated)
  - `\${var}` - Shell environment variables (escaped, runtime पर bash को pass)
  - `\$(cmd)` - Shell command substitution (escaped, runtime पर bash द्वारा execute)

ये string processing और generation patterns real-world bioinformatics वर्कफ़्लो में मिलने वाले diverse file formats और naming conventions handle करने के लिए ज़रूरी हैं।

---

## 3. Reusable Functions बनाना

Channel operators या process definitions में inline complex workflow logic readability और maintainability कम करती है। **Functions** तुम्हें इस logic को named, reusable components में extract करने देती हैं।

हमारा map operation लंबा और complex हो गया है। चलो इसे `def` keyword का उपयोग करके एक reusable function में extract करते हैं।

यह हमारे existing वर्कफ़्लो के साथ कैसा दिखता है यह illustrate करने के लिए, नीचे दिया गया modification करो, `separateMetadata` नाम का एक reusable function define करने के लिए `def` का उपयोग करते हुए:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def fastq_path = file(row.file_path)

        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [
            sample_num: m[0][2].toInteger(),
            lane: m[0][3],
            read: m[0][4],
            chunk: m[0][5]
        ] : [:]

        def priority = sample_meta.quality > 40 ? 'high' : 'normal'
        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

इस logic को एक function में extract करके, हमने actual workflow logic को कुछ बहुत cleaner में reduce कर दिया है:

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

यह workflow logic को एक नज़र में पढ़ना और समझना बहुत आसान बनाता है। Function `separateMetadata` metadata parse करने और enrich करने के लिए सभी complex logic encapsulate करता है, जिससे यह reusable और testable बनता है।

वर्कफ़्लो चलाओ यह सुनिश्चित करने के लिए कि यह अभी भी काम करता है:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

Output दोनों processes को successfully complete होते दिखाना चाहिए। वर्कफ़्लो अब बहुत cleaner और maintain करने में आसान है, सभी complex metadata processing logic `separateMetadata` function में encapsulated है।

### सारांश

इस section में, तुमने **function creation** सीखी:

- **`def` से functions define करना**: Named functions बनाने के लिए keyword (Python में `def` या JavaScript में `function` की तरह)
- **Function scope**: Script level पर define किए गए functions तुम्हारे पूरे Nextflow वर्कफ़्लो में accessible हैं
- **Return values**: Functions automatically last expression return करते हैं, या explicit `return` use करते हैं
- **Cleaner code**: Complex logic को functions में extract करना किसी भी language में fundamental software engineering practice है

अगला, हम process directives में dynamic resource allocation के लिए closures का उपयोग करना explore करेंगे।

---

## 4. Closures के साथ Dynamic Resource Directives

अब तक हमने processes के `script` block में scripting का उपयोग किया है। लेकिन **closures** (Section 1.1 में introduce किए गए) process directives में भी बेहद उपयोगी हैं, विशेष रूप से dynamic resource allocation के लिए। चलो अपने FASTP process में resource directives जोड़ते हैं जो sample characteristics के आधार पर adapt करें।

### 4.1. Sample-specific resource allocation

वर्तमान में, हमारा FASTP process default resources use करता है। चलो इसे smarter बनाते हैं high-depth samples के लिए अधिक CPUs allocate करके। `modules/fastp.nf` edit करो ताकि एक dynamic `cpus` directive और एक static `memory` directive शामिल हो:

=== "बाद में"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "पहले"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

Closure `{ meta.depth > 40000000 ? 2 : 1 }` **ternary operator** (Section 1.1 में cover किया गया) का उपयोग करता है और प्रत्येक task के लिए evaluate किया जाता है, जिससे per-sample resource allocation possible होती है। High-depth samples (>40M reads) को 2 CPUs मिलते हैं, जबकि बाकी को 1 CPU मिलता है।

!!! note "नोट: Directives में Input Variables Access करना"

    Closure किसी भी input variables (जैसे यहाँ `meta`) access कर सकता है क्योंकि Nextflow इन closures को प्रत्येक task execution के context में evaluate करता है।

वर्कफ़्लो को फिर से `-ansi-log false` option के साथ चलाओ ताकि task hashes देखना आसान हो।

```bash
nextflow run main.nf -ansi-log false
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
    [bd/ff3d41] Submitted process > FASTP (2)
    [a4/a3aab2] Submitted process > FASTP (1)
    [48/6db0c9] Submitted process > FASTP (3)
    [ec/83439d] Submitted process > GENERATE_REPORT (3)
    [bd/15d7cc] Submitted process > GENERATE_REPORT (2)
    [42/699357] Submitted process > GENERATE_REPORT (1)
    ```

किसी भी task के लिए CPU allocation देखने के लिए exact `docker` command check कर सकते हो:

```console title="Docker command check करो"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

तुम्हें कुछ ऐसा दिखना चाहिए:

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

इस example में हमने एक ऐसा example चुना है जिसने 2 CPUs request किए (`--cpu-shares 2048`), क्योंकि यह एक high-depth sample था, लेकिन तुम्हें sample depth के आधार पर अलग-अलग CPU allocations दिखनी चाहिए। दूसरे tasks के लिए भी यह try करो।

### 4.2. Retry strategies

एक और powerful pattern है `task.attempt` का retry strategies के लिए उपयोग करना। यह क्यों उपयोगी है यह दिखाने के लिए, हम पहले FASTP को memory allocation उससे कम करके शुरू करेंगे जितनी उसे ज़रूरत है। `modules/fastp.nf` में `memory` directive को `1.GB` में बदलो:

=== "बाद में"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "पहले"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... और वर्कफ़्लो फिर से चलाओ:

```bash
nextflow run main.nf
```

??? failure "कमांड आउटपुट"

    ```console hl_lines="2 11"
    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      Detecting adapter sequence for read1...
      No adapter detected for read1

      .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
    ```

यह indicate करता है कि process memory limits exceed करने के कारण kill हो गया।

यह real-world वर्कफ़्लो में एक बहुत common scenario है - कभी-कभी तुम्हें पता नहीं होता कि एक task को कितनी memory चाहिए जब तक तुम इसे run नहीं करते।

अपने वर्कफ़्लो को अधिक robust बनाने के लिए, हम एक retry strategy implement कर सकते हैं जो प्रत्येक attempt पर memory allocation बढ़ाती है, एक बार फिर Groovy closure का उपयोग करते हुए। `memory` directive को modify करो ताकि base memory को `task.attempt` से multiply किया जाए, और `errorStrategy 'retry'` और `maxRetries 2` directives जोड़ो:

=== "बाद में"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { 1.GB * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "पहले"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

अब अगर process insufficient memory के कारण fail होता है, तो Nextflow अधिक memory के साथ retry करेगा:

- पहला attempt: 1 GB (task.attempt = 1)
- दूसरा attempt: 2.GB (task.attempt = 2)

... और इसी तरह, `maxRetries` limit तक।

### सारांश

Closures के साथ dynamic directives तुम्हें allow करते हैं:

- Input characteristics के आधार पर resources allocate करना
- बढ़ते resources के साथ automatic retry strategies implement करना
- कई factors (metadata, attempt number, priorities) combine करना
- Complex resource calculations के लिए conditional logic use करना

यह तुम्हारे वर्कफ़्लो को अधिक efficient (over-allocating नहीं) और अधिक robust (अधिक resources के साथ automatic retry) दोनों बनाता है।

---

## 5. Conditional Logic और Process Control

पहले, हमने channel data transform करने के लिए `.map()` के साथ scripting का उपयोग किया। अब हम डेटा के आधार पर कौन से processes execute होते हैं यह control करने के लिए conditional logic का उपयोग करेंगे—यह अलग-अलग sample types के अनुसार adapt होने वाले flexible वर्कफ़्लो के लिए ज़रूरी है।

Nextflow के [dataflow operators](https://www.nextflow.io/docs/latest/reference/operator.html) runtime पर evaluate किए गए closures लेते हैं, जो channel content के आधार पर workflow decisions drive करने के लिए conditional logic enable करते हैं।

### 5.1. `.branch()` से Routing

उदाहरण के लिए, मान लो कि हमारे sequencing samples को FASTP से तभी trim करना है जब वे human samples हों और एक निश्चित threshold से ऊपर coverage हो। Mouse samples या low-coverage samples को Trimgalore से run करना चाहिए (यह एक contrived example है, लेकिन यह बात illustrate करता है)।

हमने `modules/trimgalore.nf` में एक simple Trimgalore process provide किया है, अगर चाहो तो देख सकते हो, लेकिन इस exercise के लिए details ज़रूरी नहीं हैं। मुख्य बात यह है कि हम samples को उनके metadata के आधार पर route करना चाहते हैं।

`modules/trimgalore.nf` से नया module include करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... और फिर अपने `main.nf` वर्कफ़्लो को modify करो ताकि samples को उनके metadata के आधार पर branch किया जाए और appropriate trimming process के ज़रिए route किया जाए, इस तरह:

=== "बाद में"

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

यह modified वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

यहाँ, हमने `.branch{}` operator के अंदर छोटे लेकिन powerful conditional expressions का उपयोग करके samples को उनके metadata के आधार पर route किया है। High coverage वाले human samples `FASTP` से गुज़रते हैं, जबकि बाकी सभी samples `TRIMGALORE` से गुज़रते हैं।

### 5.2. Truthiness के साथ `.filter()` का उपयोग

Workflow execution control करने का एक और powerful pattern `.filter()` operator है, जो यह determine करने के लिए closure का उपयोग करता है कि कौन से items pipeline में आगे जाने चाहिए। Filter closure के अंदर, तुम **boolean expressions** लिखोगे जो decide करती हैं कि कौन से items pass होते हैं।

Nextflow (कई dynamic languages की तरह) में **"truthiness"** का concept है जो determine करता है कि boolean contexts में कौन से values `true` या `false` evaluate होते हैं:

- **Truthy**: Non-null values, non-empty strings, non-zero numbers, non-empty collections
- **Falsy**: `null`, empty strings `""`, zero `0`, empty collections `[]` या `[:]`, `false`

इसका मतलब है `meta.id` अकेले (explicit `!= null` के बिना) check करता है कि ID exist करती है और empty नहीं है। चलो इसका उपयोग उन samples को filter out करने के लिए करते हैं जो हमारी quality requirements पूरी नहीं करते।

Branch operation से पहले निम्नलिखित जोड़ो:

=== "बाद में"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Invalid या low-quality samples filter out करो
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

वर्कफ़्लो फिर से चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
    [94/b48eac] Submitted process > FASTP (2)
    [2c/d2b28f] Submitted process > GENERATE_REPORT (2)
    [65/2e3be4] Submitted process > GENERATE_REPORT (1)
    [94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
    [3e/0d8664] Submitted process > TRIMGALORE (1)
    [6a/9137b0] Submitted process > FASTP (1)
    [6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
    [83/577ac0] Submitted process > GENERATE_REPORT (3)
    [a2/5117de] Re-submitted process > FASTP (1)
    [1f/a1a4ca] Re-submitted process > FASTP (2)
    ```

क्योंकि हमने एक filter चुना है जो कुछ samples को exclude करता है, कम tasks execute हुए।

Filter expression `meta.id && meta.organism && meta.depth >= 25000000` truthiness को explicit comparisons के साथ combine करता है:

- `meta.id && meta.organism` check करता है कि दोनों fields exist करती हैं और non-empty हैं (truthiness का उपयोग करके)
- `meta.depth >= 25000000` explicit comparison के साथ sufficient sequencing depth ensure करता है

!!! note "नोट: Truthiness in Practice"

    Expression `meta.id && meta.organism` लिखना इससे ज़्यादा concise है:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    यह filtering logic को बहुत cleaner और पढ़ने में आसान बनाता है।

### सारांश

इस section में, तुमने Nextflow operators जैसे `.branch{}` और `.filter{}` के closure interfaces का उपयोग करके workflow execution control करने के लिए conditional logic use करना सीखा, concise conditional expressions लिखने के लिए truthiness का लाभ उठाते हुए।

हमारी पाइपलाइन अब intelligently samples को appropriate processes में route करती है, लेकिन production वर्कफ़्लो को invalid data gracefully handle करने की ज़रूरत है। चलो अपने वर्कफ़्लो को missing या null values के खिलाफ robust बनाते हैं।

---

## 6. Safe Navigation और Elvis Operators

हमारा `separateMetadata` function वर्तमान में assume करता है कि सभी CSV fields present और valid हैं। लेकिन incomplete data के साथ क्या होगा? चलो पता करते हैं।

### 6.1. समस्या: ऐसी Properties Access करना जो Exist नहीं करतीं

मान लो हम optional sequencing run information के लिए support जोड़ना चाहते हैं। कुछ labs में, samples में sequencing run ID या batch number के लिए एक additional field हो सकता है, लेकिन हमारे current CSV में यह column नहीं है। चलो फिर भी इसे access करने की कोशिश करते हैं।

`separateMetadata` function को modify करो ताकि एक run_id field शामिल हो:

=== "बाद में"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

अब वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

यह NullPointerException के साथ crash होता है।

समस्या यह है कि `row.run_id` `null` return करता है क्योंकि हमारे CSV में `run_id` column exist नहीं करता। जब हम `null` पर `.toUpperCase()` call करने की कोशिश करते हैं, तो यह crash होता है। यहीं safe navigation operator काम आता है।

### 6.2. Safe Navigation Operator (`?.`)

Safe navigation operator (`?.`) null value पर call किए जाने पर exception throw करने की बजाय `null` return करता है। अगर `?.` से पहले का object `null` है, तो पूरा expression method execute किए बिना `null` evaluate होता है।

Safe navigation use करने के लिए function update करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

फिर से चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    <!-- TODO: output -->
    ```

Crash नहीं हुआ! वर्कफ़्लो अब missing field को gracefully handle करता है। जब `row.run_id` `null` होता है, तो `?.` operator `.toUpperCase()` call को prevent करता है, और `run_id` exception cause करने की बजाय `null` बन जाता है।

### 6.3. Default Values के लिए Elvis Operator (`?:`)

Elvis operator (`?:`) तब default values देता है जब left side "falsy" हो (जैसा पहले explain किया गया)। इसका नाम Elvis Presley के नाम पर रखा गया है क्योंकि `?:` sideways देखने पर उनके famous hair और eyes जैसा दिखता है!

अब जब हम safe navigation use कर रहे हैं, तो उस field के बिना samples के लिए `run_id` `null` होगा। चलो Elvis operator का उपयोग करके एक default value provide करते हैं और इसे अपने `sample_meta` map में जोड़ते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
        sample_meta.run = run_id
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

Results देखने के लिए वर्कफ़्लो में एक `view()` operator भी जोड़ो:

=== "बाद में"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

और वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

बढ़िया! अब सभी samples में एक `run` field है जिसमें या तो उनका actual run ID (uppercase में) है या default value 'UNSPECIFIED' है। `?.` और `?:` का combination दोनों safety (कोई crash नहीं) और sensible defaults देता है।

`.view()` operator को अब हटा दो क्योंकि हमने confirm कर लिया है कि यह काम करता है।

!!! tip "सुझाव: Safe Navigation और Elvis को Combine करना"

    Pattern `value?.method() ?: 'default'` production वर्कफ़्लो में common है:

    - `value?.method()` - Safely method call करता है, अगर `value` `null` है तो `null` return करता है
    - `?: 'default'` - अगर result `null` है तो fallback देता है

    यह pattern missing/incomplete data को gracefully handle करता है।

इन operators को functions, operator closures (`.map{}`, `.filter{}`), process scripts, और config files में consistently use करो। ये real-world data handle करते समय crashes prevent करते हैं।

### सारांश

- **Safe navigation (`?.`)**: Null values पर crashes prevent करता है - exception throw करने की बजाय null return करता है
- **Elvis operator (`?:`)**: Defaults देता है - `value ?: 'default'`
- **Combining**: `value?.method() ?: 'default'` common pattern है

ये operators वर्कफ़्लो को incomplete data के प्रति resilient बनाते हैं - real-world काम के लिए ज़रूरी।

---

## 7. `error()` और `log.warn` से Validation

कभी-कभी तुम्हें workflow को immediately रोकने की ज़रूरत होती है अगर input parameters invalid हों। Nextflow में, तुम validation logic implement करने के लिए `error()` और `log.warn` जैसे built-in functions, साथ ही `if` statements और boolean logic जैसे standard programming constructs का उपयोग कर सकते हो। चलो अपने वर्कफ़्लो में validation जोड़ते हैं।

अपने workflow block से पहले एक validation function बनाओ, इसे workflow से call करो, और CSV file path के लिए एक parameter use करने के लिए channel creation बदलो। अगर parameter missing है या file exist नहीं करती, तो clear message के साथ execution रोकने के लिए `error()` call करो।

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Check करो कि input parameter provide किया गया है
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // Check करो कि CSV file exist करती है
        if (!file(params.input).exists()) {
            error("Input CSV file not found: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

अब CSV file के बिना run करने की कोशिश करो:

```bash
nextflow run main.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Input CSV file path not provided. Please specify --input <file.csv>
    ```

वर्कफ़्लो बाद में mysteriously fail होने की बजाय immediately clear error message के साथ रुक जाता है।

अब एक non-existent file के साथ run करो:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Input CSV file not found: ./data/nonexistent.csv
    ```

अंत में, correct file के साथ run करो:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "कमांड आउटपुट"

    ```console
    <!-- TODO: output -->
    ```

इस बार यह successfully run होता है।

तुम `separateMetadata` function के अंदर भी validation जोड़ सकते हो। चलो non-fatal `log.warn` का उपयोग करके low sequencing depth वाले samples के लिए warnings issue करते हैं, लेकिन फिर भी वर्कफ़्लो को continue करने देते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validate करो कि data sense बनाता है
        if (sample_meta.depth < 30000000) {
            log.warn "Low sequencing depth for ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Original CSV के साथ वर्कफ़्लो फिर से चलाओ:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Low sequencing depth for sample_002: 25000000
    ```

हम samples में से एक के low sequencing depth के बारे में एक warning देखते हैं।

### सारांश

- **`error()`**: Clear message के साथ workflow immediately रोकता है
- **`log.warn`**: Workflow रोके बिना warnings issue करता है
- **Early validation**: Processing से पहले inputs check करो ताकि helpful errors के साथ fast fail हो
- **Validation functions**: Reusable validation logic बनाओ जिसे workflow start पर call किया जा सके

Proper validation workflows को अधिक robust और user-friendly बनाती है problems को early clear error messages के साथ catch करके।

---

## 8. Workflow Event Handlers

अब तक, हम अपने workflow scripts और process definitions में code लिख रहे थे। लेकिन एक और important feature है जो तुम्हें जाननी चाहिए: workflow event handlers।

Event handlers closures हैं जो तुम्हारे workflow के lifecycle में specific points पर run होते हैं। ये logging, notifications, या cleanup operations जोड़ने के लिए perfect हैं। इन handlers को तुम्हारे workflow script में तुम्हारी workflow definition के साथ define किया जाना चाहिए।

### 8.1. `onComplete` Handler

सबसे commonly used event handler `onComplete` है, जो तुम्हारा workflow finish होने पर run होता है (चाहे वह succeed हुआ हो या fail)। चलो अपनी pipeline results summarize करने के लिए एक जोड़ते हैं।

Event handler को अपनी `main.nf` file में, अपनी workflow definition के अंदर जोड़ो:

=== "बाद में"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

यह closure workflow complete होने पर run होता है। अंदर, तुम्हारे पास `workflow` object तक access है जो execution के बारे में useful properties देता है।

अपना वर्कफ़्लो चलाओ और तुम्हें अंत में यह summary दिखेगी!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:14:24.885384+01:00
    Duration    : 2.9s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0
    ```

चलो इसे conditional logic जोड़कर और उपयोगी बनाते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline completed successfully!"
            } else {
                println "❌ Pipeline failed!"
                println "Error: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

अब हमें एक और informative summary मिलती है, जिसमें success/failure message और specified होने पर output directory शामिल है:

<!-- TODO: add run command -->

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:16:00.522569+01:00
    Duration    : 3.6s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0

    ✅ Pipeline completed successfully!
    ```

तुम file operations का उपयोग करके summary को एक file में भी लिख सकते हो:

````groovy title="main.nf - Summary को file में लिखना"
workflow {
    // ... तुम्हारा workflow code ...

    workflow.onComplete = {
        def summary = """
        Pipeline Execution Summary
        ===========================
        Completed: ${workflow.complete}
        Duration : ${workflow.duration}```groovy title="main.nf - Summary को file में लिखना"
workflow {
    // ... तुम्हारा workflow code ...

    workflow.onComplete = {
        def summary = """
        Pipeline Execution Summary
        ===========================
        Completed: ${workflow.complete}
        Duration : ${workflow.duration}
        Success  : ${workflow.success}
        Command  : ${workflow.commandLine}
        """

        println summary

        // एक log file में लिखो
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
````

### 8.2. `onError` Handler

`onComplete` के अलावा, एक और event handler है जिसे तुम use कर सकते हो: `onError`, जो केवल तभी run होता है जब workflow fail हो:

```groovy title="main.nf - onError handler"
workflow {
    // ... तुम्हारा workflow code ...

    workflow.onError = {
        println "="* 50
        println "Pipeline execution failed!"
        println "Error message: ${workflow.errorMessage}"
        println "="* 50

        // Detailed error log लिखो
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Workflow Error Report
        =====================
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
}
```

तुम अपने workflow script में कई handlers एक साथ use कर सकते हो:

```groovy title="main.nf - Combined handlers"
workflow {
    // ... तुम्हारा workflow code ...

    workflow.onError = {
        println "Workflow failed: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUCCESS ✅" : "FAILED ❌"

        println """
        Pipeline finished: ${status}
        Duration: ${duration_mins} minutes
        """
    }
}
```

### सारांश

इस section में, तुमने सीखा:

- **Event handler closures**: तुम्हारे workflow script में closures जो अलग-अलग lifecycle points पर run होते हैं
- **`onComplete` handler**: Execution summaries और result reporting के लिए
- **`onError` handler**: Error handling और failures logging के लिए
- **Workflow object properties**: `workflow.success`, `workflow.duration`, `workflow.errorMessage`, आदि access करना

Event handlers दिखाते हैं कि तुम अपने workflow scripts में Nextflow language की पूरी power का उपयोग कैसे कर सकते हो sophisticated logging और notification capabilities जोड़ने के लिए।

---

## सारांश

बधाई हो, तुमने यह पूरा कर लिया!

इस side quest में, तुमने एक comprehensive sample processing pipeline बनाई जो basic metadata handling से एक sophisticated, production-ready वर्कफ़्लो तक विकसित हुई।
प्रत्येक section ने पिछले पर build किया, यह demonstrate करते हुए कि programming constructs simple वर्कफ़्लो को powerful data processing systems में कैसे transform करते हैं, निम्नलिखित benefits के साथ:

- **Clearer code**: Dataflow बनाम scripting को समझना तुम्हें अधिक organized वर्कफ़्लो लिखने में मदद करता है
- **Robust handling**: Safe navigation और Elvis operators वर्कफ़्लो को missing data के प्रति resilient बनाते हैं
- **Flexible processing**: Conditional logic तुम्हारे वर्कफ़्लो को अलग-अलग sample types को appropriately process करने देती है
- **Adaptive resources**: Dynamic directives input characteristics के आधार पर resource usage optimize करते हैं

यह progression bioinformatics pipelines के real-world evolution को mirror करता है, कुछ samples handle करने वाले research prototypes से लेकर laboratories और institutions में हज़ारों samples process करने वाले production systems तक।
तुमने जो हर challenge solve किया और pattern सीखा वह actual problems reflect करता है जो developers Nextflow वर्कफ़्लो scale करते समय face करते हैं।

अपने काम में इन patterns को apply करने से तुम robust, production-ready वर्कफ़्लो बनाने में सक्षम होगे।

### मुख्य patterns

1.  **Dataflow बनाम Scripting:** तुमने dataflow operations (channel orchestration) और scripting (डेटा manipulate करने वाला code) के बीच अंतर करना सीखा, जिसमें Channel बनाम List पर `collect` जैसे अलग-अलग types पर operations के बीच crucial differences शामिल हैं।

    - Dataflow: channel orchestration

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: collections पर data processing

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Advanced String Processing**: तुमने file names parse करने के लिए regular expressions, processes में dynamic script generation, और variable interpolation (Nextflow बनाम Bash बनाम Shell) में महारत हासिल की।

    - Pattern matching

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Conditional return के साथ function

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - File collection से command arguments (process script block में)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Reusable Functions बनाना**: तुमने complex logic को named functions में extract करना सीखा जिन्हें channel operators से call किया जा सकता है, जिससे वर्कफ़्लो अधिक readable और maintainable बनते हैं।

    - एक named function define करो

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* संक्षिप्तता के लिए code छुपाया गया */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* संक्षिप्तता के लिए code छुपाया गया */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Workflow में named function call करो

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Closures के साथ Dynamic Resource Directives**: तुमने input characteristics के आधार पर adaptive resource allocation के लिए process directives में closures का उपयोग explore किया।

    - Named closures और composition

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Scope access के साथ closures

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Conditional Logic और Process Control**: तुमने `.branch()` और `.filter()` operators का उपयोग करके intelligent routing जोड़ी, concise conditional expressions के लिए truthiness का लाभ उठाते हुए।

    - Data को अलग-अलग workflow branches में route करने के लिए `.branch()` use करो

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Groovy Truth के साथ boolean evaluation

    ```groovy
    if (sample.files) println "Has files"
    ```

    - 'Truthiness' के साथ data subset करने के लिए `filter()` use करो

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Safe Navigation और Elvis Operators**: तुमने null-safe property access के लिए `?.` और default values देने के लिए `?:` का उपयोग करके pipeline को missing data के खिलाफ robust बनाया।

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **error() और log.warn से Validation**: तुमने inputs को early validate करना और clear error messages के साथ fast fail करना सीखा।

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Configuration Event Handlers**: तुमने logging, notifications, और lifecycle management के लिए workflow event handlers (`onComplete` और `onError`) का उपयोग करना सीखा।

    - Log और notify करने के लिए `onComplete` का उपयोग

    ```groovy
    workflow.onComplete = {
        println "Success     : ${workflow.success}"
        println "exit status : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline completed successfully!"
        } else {
            println "❌ Pipeline failed!"
            println "Error: ${workflow.errorMessage}"
        }
    }
    ```

    - Failure के case में specifically action लेने के लिए `onError` का उपयोग

    ```groovy
    workflow.onError = {
        // Detailed error log लिखो
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
    ```

### अतिरिक्त संसाधन

- [Nextflow Language Reference](https://nextflow.io/docs/latest/reference/syntax.html)
- [Nextflow Operators](https://www.nextflow.io/docs/latest/operator.html)
- [Nextflow Script Syntax](https://www.nextflow.io/docs/latest/script.html)
- [Nextflow Standard Library](https://nextflow.io/docs/latest/reference/stdlib.html)

जब तुम्हें अधिक advanced features explore करने की ज़रूरत हो तो इन resources को ज़रूर देखो।

अपनी skills practice करने और expand करने से तुम्हें फायदा होगा ताकि तुम:

- Dataflow और scripting के बीच proper separation के साथ cleaner वर्कफ़्लो लिख सको
- Nextflow, Bash, और shell variables के साथ common pitfalls से बचने के लिए variable interpolation में महारत हासिल करो
- Efficient, adaptive वर्कफ़्लो के लिए dynamic resource directives use करो
- File collections को properly formatted command-line arguments में transform करो
- Regex और string processing का उपयोग करके अलग-अलग file naming conventions और input formats gracefully handle करो
- Advanced closure patterns और functional programming का उपयोग करके reusable, maintainable code बनाओ
- Collection operations का उपयोग करके complex datasets process और organize करो
- अपने वर्कफ़्लो को production-ready बनाने के लिए validation, error handling, और logging जोड़ो
- Event handlers के साथ workflow lifecycle management implement करो

---

## आगे क्या है?

[Side Quests के menu](../) पर वापस जाओ या list में अगले topic पर जाने के लिए page के bottom right में button click करो।
