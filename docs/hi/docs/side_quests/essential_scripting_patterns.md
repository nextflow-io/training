# आवश्यक Nextflow स्क्रिप्टिंग पैटर्न

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow एक प्रोग्रामिंग भाषा है जो Java Virtual Machine पर चलती है। जबकि Nextflow [Groovy](http://groovy-lang.org/) पर बनाई गई है और इसके अधिकांश सिंटैक्स को साझा करती है, Nextflow सिर्फ "एक्सटेंशन के साथ Groovy" नहीं है -- यह पूर्णतः निर्दिष्ट [सिंटैक्स](https://nextflow.io/docs/latest/reference/syntax.html) और [standard library](https://nextflow.io/docs/latest/reference/stdlib.html) के साथ एक स्वतंत्र भाषा है।

आप वेरिएबल, maps, और lists के लिए बुनियादी सिंटैक्स से आगे बढ़े बिना बहुत सारी Nextflow लिख सकते हैं। अधिकांश Nextflow ट्यूटोरियल workflow orchestration (channels, processes, और data flow) पर केंद्रित होते हैं, और आप सिर्फ इतने से ही आश्चर्यजनक रूप से बहुत आगे जा सकते हैं।

हालांकि, जब आपको डेटा में हेरफेर करने, जटिल फ़ाइल नामों को parse करने, conditional logic को लागू करने, या मजबूत production workflows बनाने की आवश्यकता होती है, तो यह आपके कोड के दो अलग पहलुओं के बारे में सोचने में मदद करता है: **dataflow** (channels, operators, processes, और workflows) और **scripting** (closures, functions, और process scripts के अंदर का कोड)। जबकि यह भेद कुछ हद तक मनमाना है—यह सब Nextflow कोड है—यह यह समझने के लिए एक उपयोगी मानसिक मॉडल प्रदान करता है कि आप कब अपनी pipeline को orchestrate कर रहे हैं बनाम कब आप डेटा में हेरफेर कर रहे हैं। दोनों में महारत हासिल करने से स्पष्ट, maintainable workflows लिखने की आपकी क्षमता में नाटकीय रूप से सुधार होता है।

### सीखने के लक्ष्य

यह side quest आपको बुनियादी अवधारणाओं से लेकर production-ready पैटर्न तक एक व्यावहारिक यात्रा पर ले जाता है।
हम एक साधारण CSV-reading workflow को एक sophisticated bioinformatics pipeline में बदल देंगे, इसे realistic चुनौतियों के माध्यम से चरण-दर-चरण विकसित करते हुए:

- **सीमाओं को समझना:** Dataflow operations और scripting के बीच अंतर करें, और समझें कि वे एक साथ कैसे काम करते हैं
- **डेटा में हेरफेर:** शक्तिशाली ऑपरेटर्स का उपयोग करके maps और collections को extract, transform, और subset करें
- **स्ट्रिंग प्रोसेसिंग:** Regex पैटर्न के साथ जटिल फ़ाइल नामकरण योजनाओं को parse करें और वेरिएबल interpolation में महारत हासिल करें
- **पुन: प्रयोज्य फ़ंक्शन:** साफ़, अधिक maintainable workflows के लिए जटिल logic को named functions में extract करें
- **Dynamic logic:** ऐसे processes बनाएं जो विभिन्न इनपुट प्रकारों के अनुकूल हों और dynamic resource allocation के लिए closures का उपयोग करें
- **Conditional routing:** उनकी metadata विशेषताओं के आधार पर विभिन्न processes के माध्यम से बुद्धिमानी से नमूनों को route करें
- **सुरक्षित operations:** Null-safe operators के साथ लापता डेटा को gracefully handle करें और स्पष्ट त्रुटि संदेशों के साथ inputs को validate करें
- **Configuration-based handlers:** Logging, notifications, और lifecycle management के लिए workflow event handlers का उपयोग करें

### पूर्वापेक्षाएँ

इस side quest को लेने से पहले, आपको चाहिए:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष beginner's कोर्स पूरा किया हो।
- बुनियादी Nextflow अवधारणाओं और mechanisms (processes, channels, operators, फ़ाइलों के साथ काम करना, meta data) के साथ सहज हों
- सामान्य प्रोग्रामिंग constructs (वेरिएबल, maps, lists) के साथ बुनियादी परिचय हो

यह ट्यूटोरियल प्रोग्रामिंग अवधारणाओं को समझाएगा जैसे हम उनका सामना करते हैं, इसलिए आपको व्यापक प्रोग्रामिंग अनुभव की आवश्यकता नहीं है।
हम मौलिक अवधारणाओं से शुरू करेंगे और उन्नत पैटर्न तक बनाएंगे।

---

## 0. शुरू करें

#### ट्रेनिंग codespace खोलें

यदि आपने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित अनुसार ट्रेनिंग वातावरण को खोलना सुनिश्चित करें।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाएँ

आइए उस डायरेक्टरी में जाएँ जहाँ इस ट्यूटोरियल के लिए फ़ाइलें स्थित हैं।

```bash
cd side-quests/essential_scripting_patterns
```

#### सामग्री की समीक्षा करें

आपको एक मुख्य workflow फ़ाइल और एक `data` डायरेक्टरी मिलेगी जिसमें उदाहरण डेटा फ़ाइलें हैं।

```console title="डायरेक्टरी सामग्री"
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

हमारे नमूना CSV में जैविक नमूनों के बारे में जानकारी है जिन्हें उनकी विशेषताओं के आधार पर विभिन्न प्रसंस्करण की आवश्यकता है:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

हम इस realistic डेटासेट का उपयोग व्यावहारिक प्रोग्रामिंग तकनीकों का पता लगाने के लिए करेंगे जो आप वास्तविक bioinformatics workflows में सामना करेंगे।

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### तैयारी चेकलिस्ट

लगता है कि आप शुरू करने के लिए तैयार हैं?

- [ ] मैं इस कोर्स के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता/समझती हूँ
- [ ] मेरा codespace चालू है और चल रहा है
- [ ] मैंने अपनी कार्य डायरेक्टरी को उचित रूप से सेट कर लिया है
<!-- - [ ] मैं असाइनमेंट को समझता/समझती हूँ -->

यदि आप सभी बॉक्स को चेक कर सकते हैं, तो आप जाने के लिए तैयार हैं।

---

## 1. Dataflow बनाम Scripting: सीमाओं को समझना

### 1.1. क्या क्या है इसकी पहचान करना

Nextflow workflows लिखते समय, **dataflow** (channels और processes के माध्यम से डेटा कैसे चलता है) और **scripting** (वह कोड जो डेटा में हेरफेर करता है और निर्णय लेता है) के बीच अंतर करना महत्वपूर्ण है। आइए एक workflow बनाएं जो दर्शाता है कि वे एक साथ कैसे काम करते हैं।

#### 1.1.1. बुनियादी Nextflow Workflow

एक साधारण workflow से शुरू करें जो सिर्फ CSV फ़ाइल को पढ़ता है (हमने पहले से ही यह आपके लिए `main.nf` में किया है):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

`workflow` ब्लॉक हमारी pipeline संरचना को परिभाषित करता है, जबकि `channel.fromPath()` एक फ़ाइल पथ से एक channel बनाता है। `.splitCsv()` ऑपरेटर CSV फ़ाइल को प्रोसेस करता है और प्रत्येक पंक्ति को एक map डेटा संरचना में परिवर्तित करता है।

कच्चे CSV डेटा को देखने के लिए इस workflow को चलाएं:

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

अब हम डेटा को transform करने के लिए scripting जोड़ने जा रहे हैं, `.map()` ऑपरेटर का उपयोग करते हुए जिससे आप शायद पहले से परिचित हैं। यह ऑपरेटर एक 'closure' लेता है जहाँ हम प्रत्येक आइटम को transform करने के लिए कोड लिख सकते हैं।

!!! note

    एक **closure** कोड का एक ब्लॉक है जिसे पास किया जा सकता है और बाद में निष्पादित किया जा सकता है। इसे एक ऐसे फ़ंक्शन के रूप में सोचें जिसे आप inline परिभाषित करते हैं। Closures को curly braces `{ }` के साथ लिखा जाता है और पैरामीटर ले सकते हैं। वे Nextflow operators के काम करने के तरीके के लिए मौलिक हैं और यदि आप कुछ समय से Nextflow लिख रहे हैं, तो आप शायद पहले से ही उनका उपयोग कर रहे हैं बिना यह महसूस किए!

यहाँ वह map operation कैसा दिखता है:

=== "बाद"

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

यह हमारा पहला **closure** है - एक anonymous फ़ंक्शन जिसे आप argument के रूप में पास कर सकते हैं (Python में lambdas या JavaScript में arrow functions के समान)। Closures Nextflow operators के साथ काम करने के लिए आवश्यक हैं।

Closure `{ row -> return row }` एक पैरामीटर `row` लेता है (कोई भी नाम हो सकता है: `item`, `sample`, आदि)।

जब `.map()` ऑपरेटर प्रत्येक channel आइटम को प्रोसेस करता है, तो यह उस आइटम को आपके closure को पास करता है। यहाँ, `row` एक बार में एक CSV पंक्ति रखता है।

यह परिवर्तन लागू करें और workflow चलाएं:

```bash
nextflow run main.nf
```

आप पहले जैसा ही आउटपुट देखेंगे, क्योंकि हम बस इनपुट को अपरिवर्तित रूप से वापस कर रहे हैं। यह पुष्टि करता है कि map ऑपरेटर सही ढंग से काम कर रहा है। अब आइए डेटा को transform करना शुरू करें।

#### 1.1.3. एक Map डेटा संरचना बनाना

अब हम डेटा की प्रत्येक पंक्ति को transform करने के लिए अपने closure के अंदर **scripting** logic लिखने जा रहे हैं। यह वह जगह है जहाँ हम data flow को orchestrate करने के बजाय व्यक्तिगत डेटा आइटम को प्रोसेस करते हैं।

=== "बाद"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // डेटा transformation के लिए Scripting
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

`sample_meta` map एक key-value डेटा संरचना है (Python में dictionaries, JavaScript में objects, या Ruby में hashes की तरह) जो संबंधित जानकारी संग्रहीत करती है: sample ID, organism, tissue type, sequencing depth, और quality score।

हम अपने डेटा को साफ़ करने के लिए string manipulation methods जैसे `.toLowerCase()` और `.replaceAll()` का उपयोग करते हैं, और CSV से string डेटा को उपयुक्त numeric types में परिवर्तित करने के लिए type conversion methods जैसे `.toInteger()` और `.toDouble()` का उपयोग करते हैं।

यह परिवर्तन लागू करें और workflow चलाएं:

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

अब आइए और scripting जोड़ें - इस बार डेटा मानों के आधार पर निर्णय लेने के लिए एक ternary ऑपरेटर का उपयोग करते हुए।

निम्नलिखित परिवर्तन करें:

=== "बाद"

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

Ternary ऑपरेटर एक if/else statement के लिए एक shorthand है जो पैटर्न `condition ? value_if_true : value_if_false` का अनुसरण करता है। इस पंक्ति का अर्थ है: "यदि quality 40 से अधिक है, तो 'high' का उपयोग करें, अन्यथा 'normal' का उपयोग करें"। इसका cousin, **Elvis operator** (`?:`), null या empty होने पर default values प्रदान करता है - हम इस ट्यूटोरियल में बाद में उस पैटर्न का पता लगाएंगे।

Map addition ऑपरेटर `+` मौजूदा को संशोधित करने के बजाय एक **नया map** बनाता है। यह पंक्ति एक नया map बनाती है जिसमें `sample_meta` के सभी key-value pairs प्लस नई `priority` key होती है।

!!! Note

    Closures में पास किए गए maps को कभी भी संशोधित न करें - हमेशा `+` (उदाहरण के लिए) का उपयोग करके नए बनाएं। Nextflow में, एक ही डेटा अक्सर एक साथ कई operations के माध्यम से बहता है। In-place में एक map को संशोधित करने से अप्रत्याशित side effects हो सकते हैं जब अन्य operations उसी object को reference करते हैं। नए maps बनाना सुनिश्चित करता है कि प्रत्येक operation की अपनी साफ़ copy है।

संशोधित workflow चलाएं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

हमने quality scores के आधार पर एक priority level के साथ हमारे metadata को enrich करने के लिए conditional logic को सफलतापूर्वक जोड़ा है।

#### 1.1.5. `.subMap()` के साथ Maps को Subset करना

जबकि `+` ऑपरेटर एक map में keys जोड़ता है, कभी-कभी आपको विपरीत करने की आवश्यकता होती है - केवल विशिष्ट keys को extract करें। `.subMap()` method इसके लिए बिल्कुल सही है।

आइए हमारे metadata का एक सरलीकृत संस्करण बनाने के लिए एक पंक्ति जोड़ें जिसमें केवल पहचान फ़ील्ड हों:

=== "बाद"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // डेटा transformation के लिए Scripting
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "केवल ID फ़ील्ड: ${id_only}"

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
                // डेटा transformation के लिए Scripting
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

संशोधित workflow चलाएं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    केवल ID फ़ील्ड: [id:sample_001, organism:human, tissue:liver]
    केवल ID फ़ील्ड: [id:sample_002, organism:mouse, tissue:brain]
    केवल ID फ़ील्ड: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

यह `view()` operation द्वारा प्रदर्शित पूर्ण metadata और extracted subset दोनों को दिखाता है जिसे हमने `println` के साथ प्रिंट किया।

`.subMap()` method keys की एक सूची लेता है और केवल उन keys वाला एक नया map लौटाता है। यदि मूल map में कोई key मौजूद नहीं है, तो यह बस परिणाम में शामिल नहीं है।

यह विशेष रूप से उपयोगी है जब आपको विभिन्न processes के लिए विभिन्न metadata versions बनाने की आवश्यकता है - कुछ को पूर्ण metadata की आवश्यकता हो सकती है जबकि अन्य को केवल न्यूनतम पहचान फ़ील्ड की आवश्यकता है।

अब उन println statements को हटा दें ताकि आपका workflow अपनी पिछली स्थिति में वापस आ जाए, क्योंकि हमें उन्हें आगे बढ़ने की आवश्यकता नहीं है।

!!! tip "Map Operations सारांश"

    - **Keys जोड़ें**: `map1 + [new_key: value]` - अतिरिक्त keys के साथ नया map बनाता है
    - **Keys extract करें**: `map1.subMap(['key1', 'key2'])` - केवल निर्दिष्ट keys के साथ नया map बनाता है
    - **दोनों operations नए maps बनाते हैं** - मूल maps अपरिवर्तित रहते हैं

#### 1.1.6. Maps को संयोजित करना और परिणाम लौटाना

अब तक, हम केवल वही लौटा रहे थे जिसे Nextflow community 'meta map' कहती है, और हम उन फ़ाइलों को नज़रअंदाज कर रहे थे जिनसे वह metadata संबंधित है। लेकिन यदि आप Nextflow workflows लिख रहे हैं, तो आप शायद उन फ़ाइलों के साथ कुछ करना चाहते हैं।

आइए 2 तत्वों के एक tuple वाली एक channel संरचना आउटपुट करें: enriched metadata map और संबंधित फ़ाइल पथ। यह processes को डेटा पास करने के लिए Nextflow में एक सामान्य पैटर्न है।

=== "बाद"

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

यह परिवर्तन लागू करें और workflow चलाएं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

यह `[meta, file]` tuple संरचना processes को metadata और संबंधित फ़ाइलों दोनों को पास करने के लिए Nextflow में एक सामान्य पैटर्न है।

!!! note

    **Maps और Metadata**: Maps Nextflow में metadata के साथ काम करने के लिए मौलिक हैं। Metadata maps के साथ काम करने के अधिक विस्तृत स्पष्टीकरण के लिए, [Working with metadata](./metadata.md) side quest देखें।

हमारा workflow मुख्य पैटर्न को प्रदर्शित करता है: **dataflow operations** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) यह orchestrate करते हैं कि पाइपलाइन के माध्यम से डेटा कैसे चलता है, जबकि **scripting** (maps `[key: value]`, string methods, type conversions, ternary operators) `.map()` closure के अंदर व्यक्तिगत डेटा आइटम के transformation को संभालता है।

### 1.2. विभिन्न प्रकारों को समझना: Channel बनाम List

अब तक, बहुत अच्छा, हम dataflow operations और scripting के बीच अंतर कर सकते हैं। लेकिन क्या होगा जब Nextflow standard library में दोनों संदर्भों में एक ही method का नाम मौजूद हो?

एक बिल्कुल सही उदाहरण `collect` method है, जो channel types और List types दोनों के लिए मौजूद है। List पर `collect()` method प्रत्येक तत्व को transform करती है, जबकि channel पर `collect()` ऑपरेटर सभी channel emissions को एक single-item channel में इकट्ठा करता है।

आइए कुछ नमूना डेटा के साथ इसे प्रदर्शित करें, यह देखने के लिए कि channel `collect()` ऑपरेटर क्या करता है। `collect.nf` देखें:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - कई channel emissions को एक में group करता है
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

चरण:

- Sample IDs की एक List परिभाषित करें
- `fromList()` के साथ एक channel बनाएं जो प्रत्येक sample ID को अलग से emit करता है
- जैसे यह बहता है, `view()` के साथ प्रत्येक आइटम को प्रिंट करें
- Channel के `collect()` ऑपरेटर के साथ सभी आइटम को एक single list में इकट्ठा करें
- Collected result (सभी sample IDs वाला single आइटम) को एक दूसरे `view()` के साथ प्रिंट करें

हमने channel की संरचना को बदल दिया है, लेकिन हमने डेटा को स्वयं नहीं बदला है।

इसकी पुष्टि करने के लिए workflow चलाएं:

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

`view()` हर channel emission के लिए एक आउटपुट लौटाता है, इसलिए हम जानते हैं कि यह single आउटपुट एक list में grouped सभी 3 मूल आइटम शामिल हैं।

अब आइए List पर `collect` method को action में देखें। मूल sample IDs की सूची में List के `collect` method को लागू करने के लिए `collect.nf` को संशोधित करें:

=== "बाद"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - कई channel emissions को एक में group करता है
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - प्रत्येक तत्व को transform करता है, संरचना को preserve करता है
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

- एक नया वेरिएबल `formatted_ids` परिभाषित करते हैं जो मूल list में प्रत्येक sample ID को transform करने के लिए List के `collect` method का उपयोग करता है
- `println` का उपयोग करके परिणाम को प्रिंट करते हैं

संशोधित workflow चलाएं:

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

इस बार, हमने डेटा की संरचना को नहीं बदला है, list में अभी भी 3 आइटम हैं, लेकिन हमने संशोधित मानों के साथ एक नई list बनाने के लिए List के `collect` method का उपयोग करके प्रत्येक आइटम को transform किया है। यह एक channel पर `map` ऑपरेटर का उपयोग करने के समान है, लेकिन यह channel के बजाय एक List डेटा संरचना पर काम कर रहा है।

`collect` एक चरम मामला है जिसका हम यहाँ एक बिंदु बनाने के लिए उपयोग कर रहे हैं। मुख्य सबक यह है कि जब आप workflows लिख रहे हों, तो हमेशा **डेटा संरचनाओं** (Lists, Maps, आदि) और **channels** (dataflow constructs) के बीच अंतर करें। Operations नाम साझा कर सकते हैं लेकिन जिस type पर उन्हें बुलाया जाता है, उसके आधार पर पूरी तरह से अलग तरीके से व्यवहार कर सकते हैं।

### 1.3. Spread Operator (`*.`) - प्रॉपर्टी निष्कर्षण के लिए Shorthand

List के `collect` method से संबंधित spread ऑपरेटर (`*.`) है, जो collections से properties निकालने का एक संक्षिप्त तरीका प्रदान करता है। यह अनिवार्य रूप से एक सामान्य `collect` पैटर्न के लिए syntactic sugar है।

आइए अपनी `collect.nf` फ़ाइल में एक प्रदर्शन जोड़ें:

=== "बाद"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - कई channel emissions को एक में group करता है
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - प्रत्येक तत्व को transform करता है, संरचना को preserve करता है
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // Spread operator - संक्षिप्त प्रॉपर्टी access
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

    // List.collect() - प्रत्येक तत्व को transform करता है, संरचना को preserve करता है
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

अपडेट किया गया workflow चलाएं:

```bash title="spread ऑपरेटर को test करें"
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

Spread ऑपरेटर `*.` एक सामान्य collect पैटर्न के लिए एक shorthand है:

```groovy
// ये समकक्ष हैं:
def ids = samples*.id
def ids = samples.collect { it.id }

// Method calls के साथ भी काम करता है:
def names = files*.getName()
def names = files.collect { it.getName() }
```

Spread ऑपरेटर विशेष रूप से उपयोगी है जब आपको objects की एक list से एक single प्रॉपर्टी extract करनी हो - यह पूर्ण `collect` closure लिखने की तुलना में अधिक पठनीय है।

!!! tip "Spread बनाम Collect का उपयोग कब करें"

    - **Spread (`*.`) का उपयोग करें** साधारण प्रॉपर्टी access के लिए: `samples*.id`, `files*.name`
    - **Collect का उपयोग करें** transformations या complex logic के लिए: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### निष्कर्ष

इस खंड में, आपने सीखा:

- **Dataflow बनाम scripting**: Channel operators यह orchestrate करते हैं कि आपकी pipeline के माध्यम से डेटा कैसे बहता है, जबकि scripting व्यक्तिगत डेटा आइटम को transform करता है
- **प्रकारों को समझना**: एक ही method का नाम (जैसे `collect`) जिस type पर इसे बुलाया जाता है (Channel बनाम List) के आधार पर अलग तरह से व्यवहार कर सकता है
- **संदर्भ महत्वपूर्ण है**: हमेशा जागरूक रहें कि आप channels (dataflow) या डेटा संरचनाओं (scripting) के साथ काम कर रहे हैं या नहीं

इन सीमाओं को समझना debugging, documentation, और maintainable workflows लिखने के लिए आवश्यक है।

अगला हम string processing क्षमताओं में गहराई से जाएंगे, जो वास्तविक दुनिया के डेटा को संभालने के लिए आवश्यक हैं।

---

## 2. स्ट्रिंग प्रोसेसिंग और Dynamic Script Generation

String processing में महारत हासिल करना नाजुक workflows को मजबूत pipelines से अलग करता है। यह खंड जटिल फ़ाइल नामों को parsing करने, dynamic script generation, और वेरिएबल interpolation को कवर करता है।

### 2.1. Pattern Matching और Regular Expressions

Bioinformatics फ़ाइलों में अक्सर जटिल नामकरण सम्मेलन होते हैं जो metadata को encode करते हैं। आइए pattern matching का उपयोग करके regular expressions के साथ इसे automatically extract करें।

हम अपने `main.nf` workflow में वापस जाने वाले हैं और फ़ाइल नामों से अतिरिक्त नमूना जानकारी extract करने के लिए कुछ pattern matching logic जोड़ने जा रहे हैं। हमारे डेटासेट में FASTQ फ़ाइलें Illumina-style नामकरण सम्मेलनों का पालन करती हैं जिनके नाम `SAMPLE_001_S1_L001_R1_001.fastq.gz` जैसे हैं। ये रहस्यमय लग सकते हैं, लेकिन वे वास्तव में sample ID, lane number, और read direction जैसी उपयोगी metadata को encode करते हैं। हम इन नामों को parse करने के लिए regex क्षमताओं का उपयोग करने जा रहे हैं।

अपने मौजूदा `main.nf` workflow में निम्नलिखित परिवर्तन करें:

=== "बाद"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // डेटा transformation के लिए Scripting
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
                // डेटा transformation के लिए Scripting
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

यह महत्वपूर्ण **string processing अवधारणाओं** को प्रदर्शित करता है:

1. `~/pattern/` सिंटैक्स का उपयोग करके **Regular expression literals** - यह backslashes को escape करने की आवश्यकता के बिना एक regex पैटर्न बनाता है
2. `=~` ऑपरेटर के साथ **Pattern matching** - यह एक regex पैटर्न के विरुद्ध एक string को match करने का प्रयास करता है
3. **Matcher objects** जो `[0][1]`, `[0][2]`, आदि के साथ groups को capture करते हैं - `[0]` पूरे match को संदर्भित करता है, `[1]`, `[2]`, आदि parentheses में captured groups को संदर्भित करते हैं

आइए regex पैटर्न `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$` को break down करें:

| पैटर्न              | Matches                                 | Captures                            |
| ------------------- | --------------------------------------- | ----------------------------------- |
| `^(.+)`             | शुरुआत से sample नाम                    | Group 1: sample name                |
| `_S(\d+)`           | Sample number `_S1`, `_S2`, आदि         | Group 2: sample number              |
| `_L(\d{3})`         | Lane number `_L001`                     | Group 3: lane (3 digits)            |
| `_(R[12])`          | Read direction `_R1` या `_R2`           | Group 4: read direction             |
| `_(\d{3})`          | Chunk number `_001`                     | Group 5: chunk (3 digits)           |
| `\.fastq(?:\.gz)?$` | फ़ाइल extension `.fastq` या `.fastq.gz` | Captured नहीं (?: non-capturing है) |

यह metadata को automatically extract करने के लिए Illumina-style नामकरण सम्मेलनों को parse करता है।

संशोधित workflow चलाएं:

```bash title="pattern matching को test करें"
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

यह फ़ाइल नामों से enriched metadata को दिखाता है।

### 2.2. Processes में Dynamic Script Generation

Process script blocks अनिवार्य रूप से multi-line strings हैं जो shell को pass हो जाते हैं। आप इनपुट विशेषताओं के आधार पर dynamically विभिन्न script strings generate करने के लिए **conditional logic** (if/else, ternary operators) का उपयोग कर सकते हैं। यह विभिन्न इनपुट प्रकारों—जैसे single-end बनाम paired-end sequencing reads—को process definitions को duplicate किए बिना संभालने के लिए आवश्यक है।

आइए अपने workflow में एक process जोड़ें जो इस पैटर्न को प्रदर्शित करता है। `modules/fastp.nf` खोलें और देखें:

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

Process इनपुट के रूप में FASTQ फ़ाइलें लेता है और adapters को trim करने और low-quality reads को filter करने के लिए `fastp` tool चलाता है। दुर्भाग्य से, जिस व्यक्ति ने इस process को लिखा, उसने हमारे उदाहरण डेटासेट में single-end reads की अनुमति नहीं दी। आइए इसे अपने workflow में जोड़ें और देखें कि क्या होता है:

सबसे पहले, अपने `main.nf` workflow की पहली पंक्ति में module को include करें:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

फिर `ch_samples` channel को `FASTP` process से connect करने के लिए `workflow` ब्लॉक को संशोधित करें:

=== "बाद"

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

यह संशोधित workflow चलाएं:

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

आप देख सकते हैं कि process दूसरी इनपुट फ़ाइल के लिए `null` मान के साथ `fastp` चलाने की कोशिश कर रहा है, जिससे यह fail हो रहा है। ऐसा इसलिए है क्योंकि हमारे डेटासेट में single-end reads हैं, लेकिन process paired-end reads (एक बार में दो इनपुट फ़ाइलें) की अपेक्षा करने के लिए hardcoded है।

`FASTP` process `script:` ब्लॉक में conditional logic जोड़कर इसे ठीक करें। एक if/else statement read फ़ाइलों की संख्या की जांच करता है और उसके अनुसार command को समायोजित करता है।

=== "बाद"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // सरल single-end बनाम paired-end detection
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

अब workflow single-end और paired-end reads दोनों को gracefully handle कर सकता है। Conditional logic इनपुट फ़ाइलों की संख्या की जांच करती है और `fastp` के लिए उचित command बनाती है। आइए देखें कि यह काम करता है या नहीं:

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

अच्छा लग रहा है! यदि हम उन वास्तविक commands की जांच करें जो चलाई गई थीं (अपने task hash के लिए customize करें):

```console title="निष्पादित commands की जांच करें"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

हम देख सकते हैं कि Nextflow ने single-end reads के लिए सही command को सही ढंग से चुना:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Dynamic script logic का एक अन्य सामान्य उपयोग [Nextflow for Science Genomics module](../../nf4science/genomics/02_joint_calling) में देखा जा सकता है। उस module में, जिस GATK process को call किया जा रहा है वह कई इनपुट फ़ाइलें ले सकता है, लेकिन प्रत्येक को एक सही command line बनाने के लिए `-V` के साथ prefix किया जाना चाहिए। Process इनपुट फ़ाइलों के collection (`all_gvcfs`) को सही command arguments में transform करने के लिए scripting का उपयोग करता है:

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

Process script blocks में scripting का उपयोग करने के ये पैटर्न अत्यंत शक्तिशाली हैं और कई परिदृश्यों में लागू किए जा सकते हैं - variable input types को handle करने से लेकर file collections से complex command-line arguments बनाने तक, जो आपके processes को वास्तविक दुनिया के डेटा की विविध आवश्यकताओं के लिए सचमुच अनुकूलनीय बनाते हैं।

### 2.3. Variable Interpolation: Nextflow और Shell Variables

Process scripts Nextflow variables, shell variables, और command substitutions को मिलाते हैं, प्रत्येक अलग interpolation syntax के साथ। गलत syntax का उपयोग करने से errors होती हैं। आइए इन्हें एक processing report बनाने वाले process के साथ explore करें।

`modules/generate_report.nf` module फ़ाइल पर एक नज़र डालें:

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

यह process sample ID और filename के साथ एक सरल report लिखता है। अब इसे चलाएं यह देखने के लिए कि जब हमें विभिन्न प्रकार के variables को mix करने की आवश्यकता होती है तो क्या होता है।

अपने `main.nf` में process को include करें और इसे workflow में जोड़ें:

=== "बाद"

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

अब workflow चलाएं और `results/reports/` में generated reports की जांच करें। उनमें प्रत्येक sample के बारे में बुनियादी जानकारी होनी चाहिए।

<!-- TODO: add the run command -->

??? success "कमांड आउटपुट"

    ```console
    <!-- TODO: output -->
    ```

लेकिन क्या होगा यदि हम processing कब और कहाँ हुई इसके बारे में जानकारी जोड़ना चाहते हैं? आइए report में current user, hostname, और date शामिल करने के लिए **shell** variables और कुछ command substitution का उपयोग करके process को संशोधित करें:

=== "बाद"

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

यदि आप इसे चलाते हैं, तो आपको एक error दिखेगा - Nextflow `${USER}` को एक Nextflow variable के रूप में interpret करने की कोशिश करता है जो मौजूद नहीं है।

??? failure "कमांड आउटपुट"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

हमें इसे escape करने की आवश्यकता है ताकि Bash इसे handle कर सके।

Shell variables और command substitutions को backslash (`\`) के साथ escape करके इसे ठीक करें:

=== "बाद"

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

अब यह काम करता है! Backslash (`\`) Nextflow को बताता है "इसे interpret मत करो, इसे Bash को pass कर दो।"

### निष्कर्ष

इस खंड में, आपने **string processing** तकनीकें सीखीं:

- **फ़ाइल parsing के लिए Regular expressions**: फ़ाइल नामकरण सम्मेलनों से metadata extract करने के लिए `=~` ऑपरेटर और regex पैटर्न (`~/pattern/`) का उपयोग
- **Dynamic script generation**: इनपुट विशेषताओं के आधार पर विभिन्न script strings generate करने के लिए conditional logic (if/else, ternary operators) का उपयोग
- **Variable interpolation**: यह समझना कि Nextflow कब strings को interpret करता है बनाम shell कब करता है
  - `${var}` - Nextflow variables (Nextflow द्वारा workflow compile time पर interpolated)
  - `\${var}` - Shell environment variables (escaped, runtime पर bash को passed)
  - `\$(cmd)` - Shell command substitution (escaped, runtime पर bash द्वारा executed)

ये string processing और generation पैटर्न वास्तविक bioinformatics workflows में आपको मिलने वाले विविध फ़ाइल formats और नामकरण सम्मेलनों को handle करने के लिए आवश्यक हैं।

---

## 3. पुन: प्रयोज्य फ़ंक्शन बनाना

Channel operators या process definitions में inline जटिल workflow logic पठनीयता और maintainability को कम करता है। **Functions** आपको इस logic को named, reusable components में extract करने देते हैं।

हमारा map operation लंबा और जटिल हो गया है। आइए इसे `def` keyword का उपयोग करके एक reusable function में extract करें।

यह दिखाने के लिए कि यह हमारे मौजूदा workflow के साथ कैसा दिखता है, नीचे दिए गए modification को करें, `separateMetadata` नामक एक reusable function परिभाषित करने के लिए `def` का उपयोग करते हुए:

=== "बाद"

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

इस logic को एक function में extract करके, हमने वास्तविक workflow logic को कुछ बहुत साफ़ में कम कर दिया है:

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

यह workflow logic को एक नज़र में पढ़ने और समझने में बहुत आसान बनाता है। `separateMetadata` function metadata को parse और enrich करने के लिए सभी जटिल logic को encapsulate करता है, जिससे यह reusable और testable बनता है।

यह सुनिश्चित करने के लिए workflow चलाएं कि यह अभी भी काम करता है:

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

आउटपुट को दोनों processes के सफलतापूर्वक पूरा होना दिखाना चाहिए। Workflow अब बहुत साफ़ और maintain करने में आसान है, सभी जटिल metadata processing logic `separateMetadata` function में encapsulated है।

### निष्कर्ष

इस खंड में, आपने **function creation** सीखा:

- **`def` के साथ functions परिभाषित करना**: Named functions बनाने के लिए keyword (Python में `def` या JavaScript में `function` की तरह)
- **Function scope**: Script level पर परिभाषित functions आपके पूरे Nextflow workflow में accessible हैं
- **Return values**: Functions automatically अंतिम expression return करते हैं, या explicit `return` का उपयोग करें
- **साफ़ कोड**: जटिल logic को functions में extract करना किसी भी भाषा में एक मौलिक software engineering practice है

अगला, हम dynamic resource allocation के लिए process directives में closures का उपयोग कैसे करें यह explore करेंगे।

---

## 4. Closures के साथ Dynamic Resource Directives

अब तक हमने processes के `script` block में scripting का उपयोग किया है। लेकिन **closures** (Section 1.1 में पेश किए गए) process directives में भी अविश्वसनीय रूप से उपयोगी हैं, विशेष रूप से dynamic resource allocation के लिए। आइए अपने FASTP process में resource directives जोड़ें जो sample विशेषताओं के आधार पर adapt होते हैं।

### 4.1. Sample-specific resource allocation

वर्तमान में, हमारा FASTP process default resources का उपयोग करता है। आइए इसे smarter बनाएं high-depth samples के लिए अधिक CPUs allocate करके। Dynamic `cpus` directive और एक static `memory` directive शामिल करने के लिए `modules/fastp.nf` को edit करें:

=== "बाद"

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

Closure `{ meta.depth > 40000000 ? 2 : 1 }` **ternary operator** (Section 1.1 में कवर किया गया) का उपयोग करता है और प्रत्येक task के लिए evaluate किया जाता है, जो per-sample resource allocation की अनुमति देता है। High-depth samples (>40M reads) को 2 CPUs मिलती हैं, जबकि अन्य को 1 CPU मिलती है।

!!! note "Directives में Input Variables को Access करना"

    Closure किसी भी input variable (जैसे यहाँ `meta`) को access कर सकता है क्योंकि Nextflow इन closures को प्रत्येक task execution के context में evaluate करता है।

`-ansi-log false` option के साथ workflow को फिर से चलाएं ताकि task hashes देखना आसान हो।

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

आप किसी भी दिए गए task के लिए CPU allocation देखने के लिए चलाए गए वास्तविक `docker` command की जांच कर सकते हैं:

```console title="docker command की जांच करें"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

आपको कुछ इस तरह दिखना चाहिए:

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

इस उदाहरण में हमने एक ऐसा उदाहरण चुना है जिसने 2 CPUs (`--cpu-shares 2048`) का अनुरोध किया, क्योंकि यह एक high-depth sample था, लेकिन आपको sample depth के आधार पर विभिन्न CPU allocations दिखनी चाहिए। अन्य tasks के लिए भी यह आज़माएं।

### 4.2. Retry strategies

एक और शक्तिशाली पैटर्न retry strategies के लिए `task.attempt` का उपयोग करना है। यह दिखाने के लिए कि यह क्यों उपयोगी है, हम FASTP को memory allocation को उसकी ज़रूरत से कम करके शुरू करेंगे। `modules/fastp.nf` में `memory` directive को `1.GB` में बदलें:

=== "बाद"

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

... और workflow को फिर से चलाएं:

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

यह indicate करता है कि process को memory limits से अधिक होने के कारण kill कर दिया गया था।

यह वास्तविक workflows में एक बहुत सामान्य परिदृश्य है - कभी-कभी आप नहीं जानते कि किसी task को कितनी memory की आवश्यकता होगी जब तक आप इसे चला नहीं लेते।

हमारे workflow को अधिक robust बनाने के लिए, हम एक retry strategy implement कर सकते हैं जो प्रत्येक attempt पर memory allocation बढ़ाती है, फिर से एक Groovy closure का उपयोग करते हुए। Base memory को `task.attempt` से multiply करने के लिए `memory` directive को संशोधित करें, और `errorStrategy 'retry'` और `maxRetries 2` directives जोड़ें:

=== "बाद"

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

अब यदि process अपर्याप्त memory के कारण fail होता है, तो Nextflow अधिक memory के साथ retry करेगा:

- पहला attempt: 1 GB (task.attempt = 1)
- दूसरा attempt: 2.GB (task.attempt = 2)

... और इसी तरह, `maxRetries` limit तक।

### निष्कर्ष

Closures के साथ dynamic directives आपको अनुमति देते हैं:

- इनपुट विशेषताओं के आधार पर resources allocate करें
- बढ़ते resources के साथ automatic retry strategies implement करें
- कई factors (metadata, attempt number, priorities) को combine करें
- जटिल resource calculations के लिए conditional logic का उपयोग करें

यह आपके workflows को अधिक efficient (over-allocating नहीं) और अधिक robust (अधिक resources के साथ automatic retry) दोनों बनाता है।

---

## 5. Conditional Logic और Process Control

पहले, हमने channel data को transform करने के लिए scripting के साथ `.map()` का उपयोग किया। अब हम data के आधार पर कौन से processes execute होंगे यह control करने के लिए conditional logic का उपयोग करेंगे—विभिन्न sample types के अनुकूल flexible workflows के लिए आवश्यक।

Nextflow के [dataflow operators](https://www.nextflow.io/docs/latest/reference/operator.html) runtime पर evaluate होने वाले closures लेते हैं, जो channel content के आधार पर workflow decisions को drive करने के लिए conditional logic को सक्षम करते हैं।

### 5.1. `.branch()` के साथ Routing

उदाहरण के लिए, मान लें कि हमारे sequencing samples को केवल तभी FASTP से trim करने की आवश्यकता है जब वे एक निश्चित threshold से ऊपर coverage वाले human samples हों। Mouse samples या low-coverage samples को इसके बजाय Trimgalore के साथ चलाना चाहिए (यह एक काल्पनिक उदाहरण है, लेकिन यह बात को illustrate करता है)।

हमने `modules/trimgalore.nf` में एक सरल Trimgalore process प्रदान किया है, यदि आप चाहें तो देखें, लेकिन इस exercise के लिए details महत्वपूर्ण नहीं हैं। मुख्य बात यह है कि हम samples को उनके metadata के आधार पर route करना चाहते हैं।

`modules/trimgalore.nf` से नए module को include करें:

=== "बाद"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... और फिर अपने `main.nf` workflow को samples को उनके metadata के आधार पर branch करने और उन्हें उचित trimming process के माध्यम से route करने के लिए संशोधित करें, इस तरह:

=== "बाद"

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

यह संशोधित workflow चलाएं:

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

यहाँ, हमने `.branch{}` operator के अंदर छोटे लेकिन शक्तिशाली conditional expressions का उपयोग किया है ताकि samples को उनके metadata के आधार पर route किया जा सके। High coverage वाले human samples `FASTP` से गुज़रते हैं, जबकि अन्य सभी samples `TRIMGALORE` से गुज़रते हैं।

### 5.2. Truthiness के साथ `.filter()` का उपयोग

Workflow execution को control करने का एक और शक्तिशाली पैटर्न `.filter()` operator है, जो यह निर्धारित करने के लिए एक closure का उपयोग करता है कि कौन से items pipeline में आगे बढ़ने चाहिए। Filter closure के अंदर, आप **boolean expressions** लिखेंगे जो तय करते हैं कि कौन से items pass होते हैं।

Nextflow (कई dynamic languages की तरह) में **"truthiness"** की अवधारणा है जो यह निर्धारित करती है कि boolean contexts में कौन से values `true` या `false` evaluate होते हैं:

- **Truthy**: Non-null values, non-empty strings, non-zero numbers, non-empty collections
- **Falsy**: `null`, empty strings `""`, zero `0`, empty collections `[]` या `[:]`, `false`

इसका मतलब है कि `meta.id` अकेले (explicit `!= null` के बिना) जांच करता है कि ID मौजूद है और empty नहीं है। आइए इसका उपयोग उन samples को filter out करने के लिए करें जो हमारी quality requirements को पूरा नहीं करते।

Branch operation से पहले निम्नलिखित जोड़ें:

=== "बाद"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // अमान्य या low-quality samples को filter out करें
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

Workflow को फिर से चलाएं:

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

चूंकि हमने एक ऐसा filter चुना है जो कुछ samples को exclude करता है, कम tasks execute हुए।

Filter expression `meta.id && meta.organism && meta.depth >= 25000000` truthiness को explicit comparisons के साथ combine करता है:

- `meta.id && meta.organism` जांचता है कि दोनों fields मौजूद हैं और non-empty हैं (truthiness का उपयोग करते हुए)
- `meta.depth >= 25000000` explicit comparison के साथ पर्याप्त sequencing depth सुनिश्चित करता है

!!! note "व्यवहार में Truthiness"

    Expression `meta.id && meta.organism` यह लिखने से अधिक संक्षिप्त है:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    यह filtering logic को बहुत साफ़ और पढ़ने में आसान बनाता है।

### निष्कर्ष

इस खंड में, आपने Nextflow operators जैसे `.branch{}` और `.filter{}` के closure interfaces का उपयोग करके workflow execution को control करने के लिए conditional logic का उपयोग करना सीखा, संक्षिप्त conditional expressions लिखने के लिए truthiness का लाभ उठाते हुए।

हमारी pipeline अब बुद्धिमानी से samples को उचित processes के माध्यम से route करती है, लेकिन production workflows को अमान्य data को gracefully handle करने की आवश्यकता है। आइए हमारे workflow को missing या null values के विरुद्ध robust बनाएं।

---

## 6. Safe Navigation और Elvis Operators

हमारा `separateMetadata` function वर्तमान में मानता है कि सभी CSV fields मौजूद और valid हैं। लेकिन अपूर्ण data के साथ क्या होता है? आइए पता करें।

### 6.1. समस्या: ऐसी Properties को Access करना जो मौजूद नहीं हैं

मान लें कि हम optional sequencing run information के लिए support जोड़ना चाहते हैं। कुछ labs में, samples में sequencing run ID या batch number के लिए एक अतिरिक्त field हो सकता है, लेकिन हमारे वर्तमान CSV में यह column नहीं है। आइए इसे किसी भी तरह access करने का प्रयास करें।

`separateMetadata` function को एक run_id field शामिल करने के लिए संशोधित करें:

=== "बाद"

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

अब workflow चलाएं:

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

यह NullPointerException के साथ crash हो जाता है।

समस्या यह है कि `row.run_id` `null` return करता है क्योंकि हमारे CSV में `run_id` column मौजूद नहीं है। जब हम `null` पर `.toUpperCase()` call करने की कोशिश करते हैं, तो यह crash हो जाता है। यहीं safe navigation operator बचाव में आता है।

### 6.2. Safe Navigation Operator (`?.`)

Safe navigation operator (`?.`) `null` value पर call किए जाने पर exception throw करने के बजाय `null` return करता है। यदि `?.` से पहले object `null` है, तो method को execute किए बिना पूरी expression `null` evaluate होती है।

Safe navigation का उपयोग करने के लिए function को update करें:

=== "बाद"

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

फिर से चलाएं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    <!-- TODO: output -->
    ```

कोई crash नहीं! Workflow अब missing field को gracefully handle करता है। जब `row.run_id` `null` है, तो `?.` operator `.toUpperCase()` call को रोकता है, और `run_id` exception cause करने के बजाय `null` बन जाता है।

### 6.3. Defaults के लिए Elvis Operator (`?:`)

Elvis operator (`?:`) जब left side "falsy" हो (जैसा कि पहले बताया गया) तो default values प्रदान करता है। इसका नाम Elvis Presley के नाम पर रखा गया है क्योंकि `?:` बगल से देखने पर उनके प्रसिद्ध बालों और आँखों जैसा दिखता है!

अब जब हम safe navigation का उपयोग कर रहे हैं, तो उस field के बिना samples के लिए `run_id` `null` होगा। आइए एक default value प्रदान करने और इसे हमारे `sample_meta` map में जोड़ने के लिए Elvis operator का उपयोग करें:

=== "बाद"

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

Results देखने के लिए workflow में एक `view()` operator भी जोड़ें:

=== "बाद"

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

और workflow चलाएं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

बिल्कुल सही! अब सभी samples में एक `run` field है जिसमें या तो उनका वास्तविक run ID (uppercase में) या default value 'UNSPECIFIED' है। `?.` और `?:` का combination safety (कोई crashes नहीं) और समझदार defaults दोनों प्रदान करता है।

अब `.view()` operator को हटा दें क्योंकि हमने confirm कर लिया है कि यह काम करता है।

!!! tip "Safe Navigation और Elvis को Combine करना"

    Pattern `value?.method() ?: 'default'` production workflows में सामान्य है:

    - `value?.method()` - Safely method को call करता है, `value` `null` होने पर `null` return करता है
    - `?: 'default'` - Result `null` होने पर fallback प्रदान करता है

    यह pattern missing/incomplete data को gracefully handle करता है।

इन operators का उपयोग functions, operator closures (`.map{}`, `.filter{}`), process scripts, और config files में consistently करें। वे वास्तविक डेटा handle करते समय crashes को रोकते हैं।

### निष्कर्ष

- **Safe navigation (`?.`)**: Null values पर crashes को रोकता है - exception throw करने के बजाय null return करता है
- **Elvis operator (`?:`)**: Defaults प्रदान करता है - `value ?: 'default'`
- **Combining**: `value?.method() ?: 'default'` सामान्य pattern है

ये operators workflows को अपूर्ण data के लिए resilient बनाते हैं - वास्तविक कार्य के लिए आवश्यक।

---

## 7. `error()` और `log.warn` के साथ Validation

कभी-कभी आपको input parameters अमान्य होने पर workflow को तुरंत रोकने की आवश्यकता होती है। Nextflow में, आप validation logic implement करने के लिए `error()` और `log.warn` जैसे built-in functions, साथ ही `if` statements और boolean logic जैसे standard programming constructs का उपयोग कर सकते हैं। आइए अपने workflow में validation जोड़ें।

अपने workflow block से पहले एक validation function बनाएं, इसे workflow से call करें, और CSV file path के लिए एक parameter का उपयोग करने के लिए channel creation को बदलें। यदि parameter missing है या file मौजूद नहीं है, तो एक स्पष्ट message के साथ execution रोकने के लिए `error()` call करें।

=== "बाद"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // जांचें कि input parameter प्रदान किया गया है
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // जांचें कि CSV file मौजूद है
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

अब CSV file के बिना चलाने का प्रयास करें:

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

Workflow बाद में रहस्यमय तरीके से fail होने के बजाय एक स्पष्ट error message के साथ तुरंत रुक जाता है।

अब एक non-existent file के साथ चलाएं:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Input CSV file not found: ./data/nonexistent.csv
    ```

अंत में, सही file के साथ चलाएं:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "कमांड आउटपुट"

    ```console
    <!-- TODO: output -->
    ```

इस बार यह सफलतापूर्वक चलता है।

आप `separateMetadata` function के अंदर भी validation जोड़ सकते हैं। आइए non-fatal `log.warn` का उपयोग करके low sequencing depth वाले samples के लिए warnings issue करें, लेकिन workflow को जारी रहने दें:

=== "बाद"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validate करें कि data सही है
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

मूल CSV के साथ workflow को फिर से चलाएं:

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

हमें एक sample के low sequencing depth के बारे में एक warning दिखाई देती है।

### निष्कर्ष

- **`error()`**: स्पष्ट message के साथ workflow को तुरंत रोकता है
- **`log.warn`**: Workflow को रोके बिना warnings issue करता है
- **Early validation**: Helpful errors के साथ तेज़ी से fail करने के लिए processing से पहले inputs की जांच करें
- **Validation functions**: Reusable validation logic बनाएं जिसे workflow start पर call किया जा सके

उचित validation स्पष्ट error messages के साथ जल्दी problems पकड़कर workflows को अधिक robust और user-friendly बनाती है।

---

## 8. Workflow Event Handlers

अब तक, हम अपने workflow scripts और process definitions में code लिख रहे थे। लेकिन एक और महत्वपूर्ण feature है जो आपको पता होनी चाहिए: workflow event handlers।

Event handlers closures हैं जो आपके workflow के lifecycle में विशिष्ट बिंदुओं पर चलते हैं। वे logging, notifications, या cleanup operations जोड़ने के लिए बिल्कुल सही हैं। ये handlers आपके workflow script में आपकी workflow definition के साथ परिभाषित किए जाने चाहिए।

### 8.1. `onComplete` Handler

सबसे अधिक उपयोग किया जाने वाला event handler `onComplete` है, जो आपके workflow के समाप्त होने पर (चाहे यह सफल हो या fail) चलता है। आइए हमारे pipeline results को summarize करने के लिए एक जोड़ें।

अपने `main.nf` file में, अपनी workflow definition के अंदर event handler जोड़ें:

=== "बाद"

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

यह closure workflow पूरा होने पर चलता है। अंदर, आपके पास `workflow` object तक access है जो execution के बारे में उपयोगी properties प्रदान करता है।

अपना workflow चलाएं और आपको अंत में यह summary दिखाई देगी!

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

Conditional logic जोड़कर इसे और उपयोगी बनाएं:

=== "बाद"

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

अब हमें एक और भी अधिक जानकारीपूर्ण summary मिलती है, जिसमें एक success/failure message और यदि specified हो तो output directory शामिल है:

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

आप file operations का उपयोग करके summary को एक file में भी लिख सकते हैं:

```groovy title="main.nf - Summary को file में लिखना"
workflow {
    // ... आपका workflow code ...

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

        // एक log file में लिखें
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. `onError` Handler

`onComplete` के अलावा, एक और event handler है जिसका आप उपयोग कर सकते हैं: `onError`, जो केवल workflow fail होने पर चलता है:

```groovy title="main.nf - onError handler"
workflow {
    // ... आपका workflow code ...

    workflow.onError = {
        println "="* 50
        println "Pipeline execution failed!"
        println "Error message: ${workflow.errorMessage}"
        println "="* 50

        // विस्तृत error log लिखें
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

आप अपने workflow script में कई handlers को एक साथ उपयोग कर सकते हैं:

```groovy title="main.nf - Combined handlers"
workflow {
    // ... आपका workflow code ...

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

### निष्कर्ष

इस खंड में, आपने सीखा:

- **Event handler closures**: आपके workflow script में closures जो lifecycle के विभिन्न बिंदुओं पर चलते हैं
- **`onComplete` handler**: Execution summaries और result reporting के लिए
- **`onError` handler**: Error handling और failures logging के लिए
- **Workflow object properties**: `workflow.success`, `workflow.duration`, `workflow.errorMessage`, आदि को access करना

Event handlers दिखाते हैं कि आप अपने workflow scripts में sophisticated logging और notification capabilities जोड़ने के लिए Nextflow भाषा की पूरी शक्ति का उपयोग कैसे कर सकते हैं।

---

## सारांश

बधाई हो, आपने इसे पूरा कर लिया!

इस side quest में, आपने एक व्यापक sample processing pipeline बनाई जो बुनियादी metadata handling से एक sophisticated, production-ready workflow तक विकसित हुई।
प्रत्येक खंड पिछले पर बना, यह प्रदर्शित करते हुए कि programming constructs कैसे सरल workflows को शक्तिशाली data processing systems में transform करते हैं, निम्नलिखित लाभों के साथ:

- **स्पष्ट कोड**: Dataflow बनाम scripting को समझना आपको अधिक organized workflows लिखने में मदद करता है
- **Robust handling**: Safe navigation और Elvis operators workflows को missing data के लिए resilient बनाते हैं
- **Flexible processing**: Conditional logic आपके workflows को विभिन्न sample types को उचित रूप से process करने देती है
- **Adaptive resources**: Dynamic directives input विशेषताओं के आधार पर resource usage को optimize करते हैं

यह progression bioinformatics pipelines के वास्तविक विश्व evolution को mirror करती है, कुछ samples handle करने वाले research prototypes से लेकर laboratories और institutions में हज़ारों samples process करने वाले production systems तक।
आपने जो भी चुनौती हल की और पैटर्न सीखा, वह Nextflow workflows को scale करते समय developers को आने वाली वास्तविक समस्याओं को reflect करता है।

इन patterns को अपने काम में लागू करना आपको robust, production-ready workflows बनाने में सक्षम करेगा।

### मुख्य पैटर्न

1.  **Dataflow बनाम Scripting:** आपने dataflow operations (channel orchestration) और scripting (data manipulate करने वाला code) के बीच अंतर करना सीखा, जिसमें Channel बनाम List पर `collect` जैसे विभिन्न types पर operations के बीच महत्वपूर्ण अंतर शामिल हैं।

    - Dataflow: channel orchestration

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: collections पर data processing

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Advanced String Processing**: आपने file names parsing के लिए regular expressions, processes में dynamic script generation, और variable interpolation (Nextflow बनाम Bash बनाम Shell) में महारत हासिल की।

    - Pattern matching

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Conditional return के साथ Function

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - File collection को command arguments में (process script block में)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Reusable Functions बनाना**: आपने जटिल logic को named functions में extract करना सीखा जिन्हें channel operators से call किया जा सकता है, workflows को अधिक readable और maintainable बनाते हुए।

    - एक named function परिभाषित करें

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* संक्षिप्तता के लिए code hidden */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* संक्षिप्तता के लिए code hidden */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Workflow में named function को call करें

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Closures के साथ Dynamic Resource Directives**: आपने input विशेषताओं के आधार पर adaptive resource allocation के लिए process directives में closures का उपयोग करना explore किया।

    - Named closures और composition

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Scope access के साथ Closures

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Conditional Logic और Process Control**: आपने `.branch()` और `.filter()` operators का उपयोग करके intelligent routing जोड़ा, संक्षिप्त conditional expressions के लिए truthiness का लाभ उठाते हुए।

    - विभिन्न workflow branches के माध्यम से data route करने के लिए `.branch()` का उपयोग करें

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Groovy Truth के साथ Boolean evaluation

    ```groovy
    if (sample.files) println "Has files"
    ```

    - 'Truthiness' के साथ data subset करने के लिए `filter()` का उपयोग करें

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Safe Navigation और Elvis Operators**: आपने null-safe property access के लिए `?.` और default values प्रदान करने के लिए `?:` का उपयोग करके pipeline को missing data के विरुद्ध robust बनाया।

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **error() और log.warn के साथ Validation**: आपने inputs को जल्दी validate करना और स्पष्ट error messages के साथ तेज़ी से fail करना सीखा।

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Configuration Event Handlers**: आपने logging, notifications, और lifecycle management के लिए workflow event handlers (`onComplete` और `onError`) का उपयोग करना सीखा।

    - Logging और notify करने के लिए `onComplete` का उपयोग

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

    - Failure के मामले में विशेष रूप से action लेने के लिए `onError` का उपयोग

    ```groovy
    workflow.onError = {
        // विस्तृत error log लिखें
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

जब आपको अधिक advanced features explore करने की आवश्यकता हो तो इन resources की जांच अवश्य करें।

आप अपने skills का अभ्यास और विस्तार करके लाभ उठाएंगे:

- Dataflow और scripting के बीच उचित separation के साथ साफ़ workflows लिखें
- Nextflow, Bash, और shell variables के साथ सामान्य pitfalls से बचने के लिए variable interpolation में महारत हासिल करें
- Efficient, adaptive workflows के लिए dynamic resource directives का उपयोग करें
- File collections को properly formatted command-line arguments में transform करें
- Regex और string processing का उपयोग करके विभिन्न file naming conventions और input formats को gracefully handle करें
- Advanced closure patterns और functional programming का उपयोग करके reusable, maintainable code बनाएं
- Collection operations का उपयोग करके complex datasets को process और organize करें
- अपने workflows को production-ready बनाने के लिए validation, error handling, और logging जोड़ें
- Event handlers के साथ workflow lifecycle management implement करें

---

## आगे क्या?

[Side Quests के मेनू](./index.md) पर वापस जाएं या सूची में अगले विषय पर जाने के लिए पेज के नीचे दाईं ओर बटन पर क्लिक करें।
