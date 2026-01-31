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

                def m = (fastq_path.name =~ /^(.
