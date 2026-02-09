# आवश्यक Nextflow स्क्रिप्टिंग पैटर्न

Nextflow एक प्रोग्रामिंग भाषा है जो Java Virtual Machine पर चलती है। जबकि Nextflow [Groovy](http://groovy-lang.org/) पर बनी है और इसका अधिकांश syntax साझा करती है, Nextflow सिर्फ "Groovy with extensions" से कहीं अधिक है -- यह एक स्वतंत्र भाषा है जिसमें पूरी तरह से निर्दिष्ट [syntax](https://nextflow.io/docs/latest/reference/syntax.html) और [standard library](https://nextflow.io/docs/latest/reference/stdlib.html) है।

तुम variables, maps, और lists के लिए बुनियादी syntax से आगे बढ़े बिना भी बहुत सारा Nextflow लिख सकते हो। अधिकांश Nextflow tutorials workflow orchestration (channels, processes, और data flow) पर ध्यान केंद्रित करते हैं, और तुम सिर्फ उसी के साथ आश्चर्यजनक रूप से दूर तक जा सकते हो।

हालांकि, जब तुम्हें data को manipulate करने, जटिल filenames को parse करने, conditional logic को implement करने, या robust production workflows बनाने की आवश्यकता होती है, तो अपने code के दो अलग-अलग पहलुओं के बारे में सोचना मददगार होता है: **dataflow** (channels, operators, processes, और workflows) और **scripting** (closures, functions, और process scripts के अंदर का code)। जबकि यह अंतर कुछ हद तक मनमाना है—यह सब Nextflow code है—यह यह समझने के लिए एक उपयोगी mental model प्रदान करता है कि तुम कब अपनी pipeline को orchestrate कर रहे हो बनाम कब तुम data को manipulate कर रहे हो। दोनों में महारत हासिल करना स्पष्ट, maintainable workflows लिखने की तुम्हारी क्षमता को नाटकीय रूप से बेहतर बनाता है।

### सीखने के लक्ष्य

यह side quest तुम्हें बुनियादी अवधारणाओं से production-ready पैटर्न तक एक hands-on यात्रा पर ले जाता है।
हम एक सरल CSV-reading workflow को एक sophisticated bioinformatics pipeline में बदलेंगे, इसे realistic challenges के माध्यम से step-by-step विकसित करते हुए:

- **सीमाओं को समझना:** Dataflow operations और scripting के बीच अंतर करना, और समझना कि वे एक साथ कैसे काम करते हैं
- **Data manipulation:** Powerful operators का उपयोग करके maps और collections को extract, transform, और subset करना
- **String processing:** Regex patterns के साथ जटिल file naming schemes को parse करना और variable interpolation में महारत हासिल करना
- **Reusable functions:** Cleaner, अधिक maintainable workflows के लिए जटिल logic को named functions में extract करना
- **Dynamic logic:** ऐसे processes बनाना जो विभिन्न input types के अनुकूल हों और dynamic resource allocation के लिए closures का उपयोग करें
- **Conditional routing:** Samples को उनकी metadata विशेषताओं के आधार पर विभिन्न processes के माध्यम से intelligently route करना
- **Safe operations:** Null-safe operators के साथ missing data को gracefully handle करना और स्पष्ट error messages के साथ inputs को validate करना
- **Configuration-based handlers:** Logging, notifications, और lifecycle management के लिए workflow event handlers का उपयोग करना

### पूर्वापेक्षाएँ

इस side quest को शुरू करने से पहले, तुम्हें चाहिए:

- [Hello Nextflow](../hello_nextflow/README.md) tutorial या समकक्ष beginner's course पूरा किया हो।
- बुनियादी Nextflow अवधारणाओं और mechanisms (processes, channels, operators, files के साथ काम करना, meta data) का उपयोग करने में सहज हो
- सामान्य programming constructs (variables, maps, lists) से बुनियादी परिचित हो

यह tutorial programming अवधारणाओं को समझाएगा जैसे-जैसे हम उनका सामना करते हैं, इसलिए तुम्हें व्यापक programming अनुभव की आवश्यकता नहीं है।
हम मौलिक अवधारणाओं से शुरू करेंगे और advanced पैटर्न तक बनाएंगे।

---

## 0. शुरू करना

#### Training codespace खोलें

यदि तुमने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित अनुसार training environment खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Project directory में जाओ

चलो उस directory में चलते हैं जहाँ इस tutorial के लिए files स्थित हैं।

```bash
cd side-quests/essential_scripting_patterns
```

#### Materials की समीक्षा करो

तुम्हें एक main workflow file और example data files वाली एक `data` directory मिलेगी।

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

हमारी sample CSV में biological samples के बारे में जानकारी है जिन्हें उनकी विशेषताओं के आधार पर अलग-अलग processing की आवश्यकता है:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

हम इस realistic dataset का उपयोग practical programming techniques का पता लगाने के लिए करेंगे जिनका तुम real bioinformatics workflows में सामना करोगे।

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### तैयारी checklist

लगता है कि तुम dive करने के लिए तैयार हो?

- [ ] मैं इस course के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता हूँ
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी working directory उचित रूप से सेट कर ली है
<!-- - [ ] I understand the assignment -->

यदि तुम सभी boxes को check कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. Dataflow vs Scripting: सीमाओं को समझना

### 1.1. क्या है क्या की पहचान करना

Nextflow workflows लिखते समय, **dataflow** (data channels और processes के माध्यम से कैसे चलता है) और **scripting** (वह code जो data को manipulate करता है और निर्णय लेता है) के बीच अंतर करना महत्वपूर्ण है। चलो एक workflow बनाते हैं जो दर्शाता है कि वे एक साथ कैसे काम करते हैं।

#### 1.1.1. बुनियादी Nextflow Workflow

एक सरल workflow से शुरू करो जो सिर्फ CSV file को पढ़ता है (हमने पहले से ही `main.nf` में तुम्हारे लिए यह किया है):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

`workflow` block हमारी pipeline structure को define करता है, जबकि `channel.fromPath()` एक file path से एक channel बनाता है। `.splitCsv()` operator CSV file को process करता है और प्रत्येक row को एक map data structure में convert करता है।

Raw CSV data देखने के लिए इस workflow को चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    Launching `main.nf` [marvelous_tuckerman] DSL2 - revision: 6113e05c17

    [sample_id:SAMPLE_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/SAMPLE_001_S1_L001_R1_001.fastq, quality_score:38.5]
    [sample_id:SAMPLE_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/SAMPLE_002_S2_L001_R1_001.fastq, quality_score:35.2]
    [sample_id:SAMPLE_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/SAMPLE_003_S3_L001_R1_001.fastq, quality_score:42.1]
    ```

#### 1.1.2. Map Operator जोड़ना

अब हम `.map()` operator का उपयोग करके data को transform करने के लिए scripting जोड़ने जा रहे हैं, जिससे तुम शायद पहले से ही परिचित होगे। यह operator एक 'closure' लेता है जहाँ हम प्रत्येक item को transform करने के लिए code लिख सकते हैं।

!!! note

    एक **closure** code का एक block है जिसे pass किया जा सकता है और बाद में execute किया जा सकता है। इसे एक ऐसे function के रूप में सोचो जिसे तुम inline define करते हो। Closures को curly braces `{ }` के साथ लिखा जाता है और parameters ले सकते हैं। वे Nextflow operators के काम करने के तरीके के लिए मौलिक हैं और यदि तुम कुछ समय से Nextflow लिख रहे हो, तो तुम पहले से ही उनका उपयोग कर रहे हो सकते हो बिना यह महसूस किए!

यहाँ वह map operation कैसा दिखता है:

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

यह हमारा पहला **closure** है - एक anonymous function जिसे तुम एक argument के रूप में pass कर सकते हो (Python में lambdas या JavaScript में arrow functions के समान)। Closures Nextflow operators के साथ काम करने के लिए आवश्यक हैं।

Closure `{ row -> return row }` एक parameter `row` लेता है (कोई भी नाम हो सकता है: `item`, `sample`, आदि)।

जब `.map()` operator प्रत्येक channel item को process करता है, तो यह उस item को तुम्हारे closure में pass करता है। यहाँ, `row` एक समय में एक CSV row रखता है।

इस परिवर्तन को लागू करो और workflow चलाओ:

```bash
nextflow run main.nf
```

तुम पहले जैसा ही output देखोगे, क्योंकि हम बस input को unchanged return कर रहे हैं। यह confirm करता है कि map operator सही तरीके से काम कर रहा है। अब चलो data को transform करना शुरू करते हैं।

#### 1.1.3. एक Map Data Structure बनाना

अब हम प्रत्येक row of data को transform करने के लिए अपने closure के अंदर **scripting** logic लिखने जा रहे हैं। यहीं पर हम data flow को orchestrate करने के बजाय individual data items को process करते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting for data transformation
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

हम अपने data को साफ करने के लिए `.toLowerCase()` और `.replaceAll()` जैसे string manipulation methods का उपयोग करते हैं, और CSV से string data को उपयुक्त numeric types में convert करने के लिए `.toInteger()` और `.toDouble()` जैसे type conversion methods का उपयोग करते हैं।

इस परिवर्तन को लागू करो और workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. Conditional Logic जोड़ना

अब चलो और scripting जोड़ते हैं - इस बार data values के आधार पर निर्णय लेने के लिए एक ternary operator का उपयोग करते हुए।

निम्नलिखित परिवर्तन करो:

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

Ternary operator एक if/else statement के लिए एक shorthand है जो `condition ? value_if_true : value_if_false` pattern का पालन करता है। इस line का मतलब है: "यदि quality 40 से अधिक है, तो 'high' का उपयोग करो, अन्यथा 'normal' का उपयोग करो"। इसका cousin, **Elvis operator** (`?:`), default values प्रदान करता है जब कुछ null या empty होता है - हम इस tutorial में बाद में उस pattern का पता लगाएंगे।

Map addition operator `+` existing map को modify करने के बजाय एक **नया map** बनाता है। यह line एक नया map बनाती है जिसमें `sample_meta` से सभी key-value pairs और नई `priority` key शामिल है।

!!! Note

    Closures में pass किए गए maps को कभी भी modify न करो - हमेशा `+` का उपयोग करके नए बनाओ (उदाहरण के लिए)। Nextflow में, एक ही data अक्सर एक साथ कई operations के माध्यम से flow करता है। In-place map को modify करना unpredictable side effects का कारण बन सकता है जब अन्य operations उसी object को reference करते हैं। नए maps बनाना सुनिश्चित करता है कि प्रत्येक operation की अपनी clean copy है।

Modified workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

हमने quality scores के आधार पर priority level के साथ अपने metadata को enrich करने के लिए सफलतापूर्वक conditional logic जोड़ा है।

#### 1.1.5. `.subMap()` के साथ Maps को Subset करना

जबकि `+` operator एक map में keys जोड़ता है, कभी-कभी तुम्हें विपरीत करने की आवश्यकता होती है - केवल specific keys को extract करना। `.subMap()` method इसके लिए perfect है।

चलो हमारे metadata का एक simplified version बनाने के लिए एक line जोड़ते हैं जिसमें केवल identification fields हों:

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting for data transformation
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "ID fields only: ${id_only}"

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
                // Scripting for data transformation
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

Modified workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

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

यह `view()` operation द्वारा displayed पूर्ण metadata और हमने `println` के साथ print किए गए extracted subset दोनों को दिखाता है।

`.subMap()` method keys की एक list लेता है और केवल उन keys वाला एक नया map return करता है। यदि कोई key original map में मौजूद नहीं है, तो इसे बस result में शामिल नहीं किया जाता है।

यह विशेष रूप से उपयोगी है जब तुम्हें विभिन्न processes के लिए अलग-अलग metadata versions बनाने की आवश्यकता होती है - कुछ को पूर्ण metadata की आवश्यकता हो सकती है जबकि अन्य को केवल minimal identification fields की आवश्यकता होती है।

अब उन println statements को हटा दो ताकि तुम्हारा workflow अपनी पिछली स्थिति में restore हो जाए, क्योंकि हमें आगे बढ़ने के लिए उनकी आवश्यकता नहीं है।

!!! tip "Map Operations सारांश"

    - **Keys जोड़ें**: `map1 + [new_key: value]` - अतिरिक्त keys के साथ नया map बनाता है
    - **Keys extract करें**: `map1.subMap(['key1', 'key2'])` - केवल निर्दिष्ट keys के साथ नया map बनाता है
    - **दोनों operations नए maps बनाते हैं** - Original maps unchanged रहते हैं

#### 1.1.6. Maps को Combine करना और Results Return करना

अब तक, हम केवल वही return कर रहे हैं जिसे Nextflow community 'meta map' कहती है, और हम उन files को ignore कर रहे हैं जिनसे वे metadata संबंधित हैं। लेकिन यदि तुम Nextflow workflows लिख रहे हो, तो तुम शायद उन files के साथ कुछ करना चाहते हो।

चलो एक channel structure output करते हैं जिसमें 2 elements का एक tuple हो: enriched metadata map और corresponding file path। यह Nextflow में processes को data pass करने के लिए एक सामान्य pattern है।

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

इस परिवर्तन को लागू करो और workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

यह `[meta, file]` tuple structure Nextflow में processes को metadata और associated files दोनों pass करने के लिए एक सामान्य pattern है।

!!! note

    **Maps और Metadata**: Maps Nextflow में metadata के साथ काम करने के लिए मौलिक हैं। Metadata maps के साथ काम करने की अधिक विस्तृत व्याख्या के लिए, [Working with metadata](./metadata.md) side quest देखो।

हमारा workflow core pattern को demonstrate करता है: **dataflow operations** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orchestrate करते हैं कि data pipeline के माध्यम से कैसे चलता है, जबकि **scripting** (maps `[key: value]`, string methods, type conversions, ternary operators) `.map()` closure के अंदर individual data items के transformation को handle करता है।

### 1.2. विभिन्न Types को समझना: Channel vs List

अब तक, बहुत अच्छा, हम dataflow operations और scripting के बीच अंतर कर सकते हैं। लेकिन क्या होगा जब दोनों contexts में एक ही method name मौजूद हो?

एक perfect उदाहरण `collect` method है, जो Nextflow standard library में channel types और List types दोनों के लिए मौजूद है। List पर `collect()` method प्रत्येक element को transform करता है, जबकि channel पर `collect()` operator सभी channel emissions को एक single-item channel में gather करता है।

चलो कुछ sample data के साथ इसे demonstrate करते हैं, यह refresh करते हुए कि channel `collect()` operator क्या करता है। `collect.nf` देखो:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - groups multiple channel emissions into one
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

Steps:

- Sample IDs की एक List define करो
- `fromList()` के साथ एक channel बनाओ जो प्रत्येक sample ID को अलग से emit करता है
- प्रत्येक item को `view()` के साथ print करो जैसे यह flow करता है
- सभी items को channel के `collect()` operator के साथ एक single list में gather करो
- Collected result (सभी sample IDs वाला single item) को दूसरे `view()` के साथ print करो

हमने channel की structure को बदल दिया है, लेकिन हमने data को नहीं बदला है।

इसे confirm करने के लिए workflow चलाओ:

```bash
nextflow run collect.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

`view()` हर channel emission के लिए एक output return करता है, इसलिए हम जानते हैं कि यह single output में सभी 3 original items एक list में grouped हैं।

अब चलो List के `collect` method को action में देखते हैं। `collect.nf` को modify करो ताकि List के `collect` method को original list of sample IDs पर apply किया जा सके:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - groups multiple channel emissions into one
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforms each element, preserves structure
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - groups multiple channel emissions into one
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

इस नए snippet में हम:

- एक नया variable `formatted_ids` define करते हैं जो original list में प्रत्येक sample ID को transform करने के लिए List के `collect` method का उपयोग करता है
- `println` का उपयोग करके result को print करते हैं

Modified workflow चलाओ:

```bash
nextflow run collect.nf
```

??? success "Command output"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

इस बार, हमने data की structure को नहीं बदला है, हमारे पास अभी भी list में 3 items हैं, लेकिन हमने modified values के साथ एक नई list produce करने के लिए List के `collect` method का उपयोग करके प्रत्येक item को transform किया है। यह channel पर `map` operator का उपयोग करने के समान है, लेकिन यह channel के बजाय एक List data structure पर operate कर रहा है।

`collect` एक extreme case है जिसका हम यहाँ एक point बनाने के लिए उपयोग कर रहे हैं। मुख्य सबक यह है कि जब तुम workflows लिख रहे हो, तो हमेशा **data structures** (Lists, Maps, आदि) और **channels** (dataflow constructs) के बीच अंतर करो। Operations नाम साझा कर सकते हैं लेकिन उस type के आधार पर पूरी तरह से अलग तरीके से व्यवहार करते हैं जिस पर उन्हें call किया जाता है।

### 1.3. Spread Operator (`*.`) - Property Extraction के लिए Shorthand

List के `collect` method से संबंधित spread operator (`*.`) है, जो collections से properties को extract करने का एक concise तरीका प्रदान करता है। यह अनिवार्य रूप से एक सामान्य `collect` pattern के लिए syntactic sugar है।

चलो अपनी `collect.nf` file में एक demonstration जोड़ते हैं:

=== "बाद में"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - groups multiple channel emissions into one
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforms each element, preserves structure
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

    // channel.collect() - groups multiple channel emissions into one
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforms each element, preserves structure
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

Updated workflow चलाओ:

```bash title="Test spread operator"
nextflow run collect.nf
```

??? success "Command output"

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

Spread operator `*.` एक सामान्य collect pattern के लिए एक shorthand है:

```groovy
// These are equivalent:
def ids = samples*.id
def ids = samples.collect { it.id }

// Also works with method calls:
def names = files*.getName()
def names = files.collect { it.getName() }
```

Spread operator विशेष रूप से उपयोगी है जब तुम्हें objects की list से एक single property को extract करने की आवश्यकता होती है - यह पूर्ण `collect` closure लिखने से अधिक readable है।

!!! tip "Spread vs Collect का उपयोग कब करें"

    - **Spread (`*.`) का उपयोग करो** simple property access के लिए: `samples*.id`, `files*.name`
    - **Collect का उपयोग करो** transformations या complex logic के लिए: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### सारांश

इस section में, तुमने सीखा:

- **Dataflow vs scripting**: Channel operators orchestrate करते हैं कि data तुम्हारी pipeline के माध्यम से कैसे flow करता है, जबकि scripting individual data items को transform करता है
- **Types को समझना**: एक ही method name (जैसे `collect`) उस type के आधार पर अलग तरीके से व्यवहार कर सकता है जिस पर इसे call किया जाता है (Channel vs List)
- **Context matters**: हमेशा aware रहो कि तुम channels (dataflow) या data structures (scripting) के साथ काम कर रहे हो

इन सीमाओं को समझना debugging, documentation, और maintainable workflows लिखने के लिए आवश्यक है।

अगले हम string processing capabilities में गहराई से dive करेंगे, जो real-world data को handle करने के लिए आवश्यक हैं।

---

## 2. String Processing और Dynamic Script Generation

String processing में महारत हासिल करना brittle workflows को robust pipelines से अलग करता है। यह section complex file names को parse करने, dynamic script generation, और variable interpolation को cover करता है।

### 2.1. Pattern Matching और Regular Expressions

Bioinformatics files में अक्सर जटिल naming conventions होते हैं जो metadata को encode करते हैं। चलो regular expressions के साथ pattern matching का उपयोग करके इसे automatically extract करते हैं।

हम अपने `main.nf` workflow पर वापस जाने वाले हैं और file names से अतिरिक्त sample information को extract करने के लिए कुछ pattern matching logic जोड़ने वाले हैं। हमारे dataset में FASTQ files Illumina-style naming conventions का पालन करती हैं जिनके नाम `SAMPLE_001_S1_L001_R1_001.fastq.gz` जैसे हैं। ये cryptic लग सकते हैं, लेकिन वे वास्तव में sample ID, lane number, और read direction जैसे उपयोगी metadata को encode करते हैं। हम इन names को parse करने के लिए regex capabilities का उपयोग करने जा रहे हैं।

अपने existing `main.nf` workflow में निम्नलिखित परिवर्तन करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Scripting for data transformation
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
                // Scripting for data transformation
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

यह मुख्य **string processing concepts** को demonstrate करता है:

1. **Regular expression literals** `~/pattern/` syntax का उपयोग करते हुए - यह backslashes को escape करने की आवश्यकता के बिना एक regex pattern बनाता है
2. **Pattern matching** `=~` operator के साथ - यह एक string को एक regex pattern के साथ match करने का प्रयास करता है
3. **Matcher objects** जो `[0][1]`, `[0][2]`, आदि के साथ groups को capture करते हैं - `[0]` पूरे match को refer करता है, `[1]`, `[2]`, आदि parentheses में captured groups को refer करते हैं

चलो regex pattern `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$` को break down करते हैं:

| Pattern             | Matches                                | Captures                           |
| ------------------- | -------------------------------------- | ---------------------------------- |
| `^(.+)`             | Start से sample name                   | Group 1: sample name               |
| `_S(\d+)`           | Sample number `_S1`, `_S2`, आदि        | Group 2: sample number             |
| `_L(\d{3})`         | Lane number `_L001`                    | Group 3: lane (3 digits)           |
| `_(R[12])`          | Read direction `_R1` या `_R2`          | Group 4: read direction            |
| `_(\d{3})`          | Chunk number `_001`                    | Group 5: chunk (3 digits)          |
| `\.fastq(?:\.gz)?$` | File extension `.fastq` या `.fastq.gz` | Captured नहीं (?: non-capturing है) |

यह metadata को automatically extract करने के लिए Illumina-style naming conventions को parse करता है।

Modified workflow चलाओ:

```bash title="Test pattern matching"
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

यह file names से enriched metadata को दिखाता है।

### 2.2. Processes में Dynamic Script Generation

Process script blocks अनिवार्य रूप से multi-line strings हैं जो shell को pass की जाती हैं। तुम input विशेषताओं के आधार पर अलग-अलग script strings को dynamically generate करने के लिए **conditional logic** (if/else, ternary operators) का उपयोग कर सकते हो। यह diverse input types—जैसे single-end vs paired-end sequencing reads—को handle करने के लिए आवश्यक है बिना process definitions को duplicate किए।

चलो अपने workflow में एक process जोड़ते हैं जो इस pattern को demonstrate करता है। `modules/fastp.nf` खोलो और देखो:

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

Process FASTQ files को input के रूप में लेता है और adapters को trim करने और low-quality reads को filter करने के लिए `fastp` tool चलाता है। दुर्भाग्य से, जिस व्यक्ति ने यह process लिखा था उसने हमारे example dataset में single-end reads के लिए allow नहीं किया। चलो इसे अपने workflow में जोड़ते हैं और देखते हैं कि क्या होता है:

पहले, अपनी `main.nf` workflow की पहली line पर module को include करो:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

फिर `ch_samples` channel को `FASTP` process से connect करने के लिए `workflow` block को modify करो:

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

इस modified workflow को चलाओ:

```bash
nextflow run main.nf
```

??? failure "Command output"

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

तुम देख सकते हो कि process दूसरी input file के लिए `null` value के साथ `fastp` चलाने की कोशिश कर रहा है, जो इसे fail करने का कारण बन रहा है। यह इसलिए है क्योंकि हमारे dataset में single-end reads हैं, लेकिन process paired-end reads (एक समय में दो input files) की उम्मीद करने के लिए hardcoded है।

`FASTP` process `script:` block में conditional logic जोड़कर इसे ठीक करो। एक if/else statement read file count को check करता है और तदनुसार command को adjust करता है।

=== "बाद में"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Simple single-end vs paired-end detection
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

अब workflow single-end और paired-end reads दोनों को gracefully handle कर सकता है। Conditional logic input files की संख्या को check करता है और `fastp` के लिए उपयुक्त command construct करता है। चलो देखते हैं कि यह काम करता है या नहीं:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_rosalind] DSL2 - revision: 04b1cd93e9

    executor >  local (3)
    [31/a8ad4d] process > FASTP (3) [100%] 3 of 3 ✔
    ```

अच्छा लग रहा है! यदि हम actual commands को check करें जो run किए गए थे (अपने task hash के लिए customize करो):

```console title="Check commands executed"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

हम देख सकते हैं कि Nextflow ने single-end reads के लिए सही command को correctly pick किया:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Dynamic script logic का एक और सामान्य उपयोग [Nextflow for Science Genomics module](../../nf4science/genomics/02_joint_calling) में देखा जा सकता है। उस module में, GATK process जिसे call किया जा रहा है वह कई input files ले सकता है, लेकिन प्रत्येक को एक सही command line बनाने के लिए `-V` के साथ prefix किया जाना चाहिए। Process input files (`all_gvcfs`) के collection को सही command arguments में transform करने के लिए scripting का उपयोग करता है:

```groovy title="command line manipulation for GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Process script blocks में scripting का उपयोग करने के ये patterns अत्यंत शक्तिशाली हैं और कई scenarios में लागू किए जा सकते हैं - variable input types को handle करने से लेकर file collections से complex command-line arguments बनाने तक, तुम्हारे processes को real-world data की विविध आवश्यकताओं के लिए वास्तव में adaptable बनाते हैं।

### 2.3. Variable Interpolation: Nextflow और Shell Variables

Process scripts Nextflow variables, shell variables, और command substitutions को mix करते हैं, प्रत्येक अलग interpolation syntax के साथ। गलत syntax का उपयोग करने से errors होती हैं। चलो एक process के साथ इन्हें explore करते हैं जो एक processing report बनाता है।

Module file `modules/generate_report.nf` देखो:

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

यह process sample ID और filename के साथ एक simple report लिखता है। अब चलो इसे चलाते हैं यह देखने के लिए कि जब हमें विभिन्न प्रकार के variables को mix करने की आवश्यकता होती है तो क्या होता है।

अपने `main.nf` में process को include करो और इसे workflow में जोड़ो:

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

अब workflow चलाओ और `results/reports/` में generated reports को check करो। उनमें प्रत्येक sample के बारे में बुनियादी जानकारी होनी चाहिए।

<!-- TODO: add the run command -->

??? success "Command output"

    ```console
    <!-- TODO: output -->
    ```

लेकिन क्या होगा यदि हम processing कब और कहाँ हुई इसके बारे में जानकारी जोड़ना चाहते हैं? चलो report में current user, hostname, और date को शामिल करने के लिए **shell** variables और थोड़ा command substitution का उपयोग करने के लिए process को modify करते हैं:

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

यदि तुम इसे चलाते हो, तो तुम्हें एक error दिखाई देगी - Nextflow `${USER}` को एक Nextflow variable के रूप में interpret करने की कोशिश करता है जो मौजूद नहीं है।

??? failure "Command output"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

हमें इसे escape करने की आवश्यकता है ताकि Bash इसे handle कर सके।

Shell variables और command substitutions को एक backslash (`\`) के साथ escape करके इसे ठीक करो:

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

इस section में, तुमने **string processing** techniques सीखीं:

- **File parsing के लिए regular expressions**: Complex file naming conventions से metadata को extract करने के लिए `=~` operator और regex patterns (`~/pattern/`) का उपयोग करना
- **Dynamic script generation**: Input विशेषताओं के आधार पर अलग-अलग script strings generate करने के लिए conditional logic (if/else, ternary operators) का उपयोग करना
- **Variable interpolation**: समझना कि Nextflow strings को कब interpret करता है बनाम shell कब करता है
  - `${var}` - Nextflow variables (workflow compile time पर Nextflow द्वारा interpolated)
  - `\${var}` - Shell environment variables (escaped, runtime पर bash को passed)
  - `\$(cmd)` - Shell command substitution (escaped, runtime पर bash द्वारा executed)

ये string processing और generation patterns उन विविध file formats और naming conventions को handle करने के लिए आवश्यक हैं जिनका तुम real-world bioinformatics workflows में सामना करोगे।

---

## 3. Reusable Functions बनाना

Channel operators या process definitions में inline complex workflow logic readability और maintainability को कम करता है। **Functions** तुम्हें इस logic को named, reusable components में extract करने देते हैं।

हमारा map operation लंबा और जटिल हो गया है। चलो `def` keyword का उपयोग करके इसे एक reusable function में extract करते हैं।

यह illustrate करने के लिए कि हमारे existing workflow के साथ यह कैसा दिखता है, नीचे दिए गए modification को करो, `separateMetadata` नामक एक reusable function को define करने के लिए `def` का उपयोग करते हुए:

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

यह workflow logic को एक नज़र में पढ़ना और समझना बहुत आसान बनाता है। Function `separateMetadata` metadata को parse करने और enrich करने के लिए सभी complex logic को encapsulate करता है, इसे reusable और testable बनाता है।

यह सुनिश्चित करने के लिए workflow चलाओ कि यह अभी भी काम करता है:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

Output को दोनों processes को successfully complete होते हुए दिखाना चाहिए। Workflow अब बहुत cleaner और maintain करने में आसान है, सभी complex metadata processing logic `separateMetadata` function में encapsulated है।

### सारांश

इस section में, तुमने **function creation** सीखा:

- **`def` के साथ functions को define करना**: Named functions बनाने के लिए keyword (Python में `def` या JavaScript में `function` की तरह)
- **Function scope**: Script level पर defined functions तुम्हारे पूरे Nextflow workflow में accessible हैं
- **Return values**: Functions automatically अंतिम expression को return करते हैं, या explicit `return` का उपयोग करते हैं
- **Cleaner code**: Complex logic को functions में extract करना किसी भी भाषा में एक मौलिक software engineering practice है

अगले, हम dynamic resource allocation के लिए process directives में closures का उपयोग करने का पता लगाएंगे।

---

## 4. Closures के साथ Dynamic Resource Directives

अब तक हमने processes के `script` block में scripting का उपयोग किया है। लेकिन **closures** (Section 1.1 में introduced) process directives में भी अविश्वसनीय रूप से उपयोगी हैं, विशेष रूप से dynamic resource allocation के लिए। चलो अपने FASTP process में resource directives जोड़ते हैं जो sample विशेषताओं के आधार पर adapt करते हैं।

### 4.1. Sample-specific resource allocation

वर्तमान में, हमारा FASTP process default resources का उपयोग करता है। चलो इसे high-depth samples के लिए अधिक CPUs allocate करके smarter बनाते हैं। `modules/fastp.nf` को edit करो ताकि एक dynamic `cpus` directive और एक static `memory` directive शामिल हो:

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

Closure `{ meta.depth > 40000000 ? 2 : 1 }` **ternary operator** (Section 1.1 में covered) का उपयोग करता है और प्रत्येक task के लिए evaluated किया जाता है, per-sample resource allocation की अनुमति देता है। High-depth samples (>40M reads) को 2 CPUs मिलते हैं, जबकि अन्य को 1 CPU मिलता है।

!!! note "Directives में Input Variables को Access करना"

    Closure किसी भी input variables (जैसे यहाँ `meta`) को access कर सकता है क्योंकि Nextflow इन closures को प्रत्येक task execution के context में evaluate करता है।

`-ansi-log false` option के साथ workflow को फिर से चलाओ ताकि task hashes को देखना आसान हो।

```bash
nextflow run main.nf -ansi-log false
```

??? success "Command output"

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

तुम किसी भी दिए गए task के लिए CPU allocation देखने के लिए run किए गए exact `docker` command को check कर सकते हो:

```console title="Check docker command"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

तुम्हें कुछ ऐसा दिखना चाहिए:

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

इस उदाहरण में हमने एक ऐसा उदाहरण चुना है जिसने 2 CPUs (`--cpu-shares 2048`) का request किया, क्योंकि यह एक high-depth sample था, लेकिन तुम्हें sample depth के आधार पर अलग-अलग CPU allocations दिखने चाहिए। अन्य tasks के लिए भी यह try करो।

### 4.2. Retry strategies

एक और powerful pattern retry strategies के लिए `task.attempt` का उपयोग करना है। यह दिखाने के लिए कि यह क्यों उपयोगी है, हम FASTP को memory allocation को उससे कम करके शुरू करने जा रहे हैं जितनी इसे आवश्यकता है। `modules/fastp.nf` में `memory` directive को `1.GB` में बदलो:

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

... और workflow को फिर से चलाओ:

```bash
nextflow run main.nf
```

??? failure "Command output"

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

यह indicate करता है कि process को memory limits exceed करने के लिए kill कर दिया गया था।

यह real-world workflows में एक बहुत सामान्य scenario है - कभी-कभी तुम बस नहीं जानते कि एक task को कितनी memory की आवश्यकता होगी जब तक तुम इसे run नहीं करते।

अपने workflow को अधिक robust बनाने के लिए, हम एक retry strategy implement कर सकते हैं जो प्रत्येक attempt पर memory allocation को बढ़ाता है, एक बार फिर Groovy closure का उपयोग करते हुए। `memory` directive को `task.attempt` से base memory को multiply करने के लिए modify करो, और `errorStrategy 'retry'` और `maxRetries 2` directives जोड़ो:

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

अब यदि process insufficient memory के कारण fail होता है, तो Nextflow अधिक memory के साथ retry करेगा:

- पहला attempt: 1 GB (task.attempt = 1)
- दूसरा attempt: 2.GB (task.attempt = 2)

... और इसी तरह, `maxRetries` limit तक।

### सारांश

Closures के साथ dynamic directives तुम्हें allow करते हैं:

- Input विशेषताओं के आधार पर resources को allocate करना
- बढ़ते resources के साथ automatic retry strategies को implement करना
- कई factors (metadata, attempt number, priorities) को combine करना
- Complex resource calculations के लिए conditional logic का उपयोग करना

यह तुम्हारे workflows को अधिक efficient (over-allocating नहीं) और अधिक robust (अधिक resources के साथ automatic retry) दोनों बनाता है।

---

## 5. Conditional Logic और Process Control

पहले, हमने channel data को transform करने के लिए scripting के साथ `.map()` का उपयोग किया। अब हम data के आधार पर यह control करने के लिए conditional logic का उपयोग करेंगे कि कौन से processes execute होते हैं—विभिन्न sample types के अनुकूल flexible workflows के लिए आवश्यक।

Nextflow के [dataflow operators](https://www.nextflow.io/docs/latest/reference/operator.html) closures लेते हैं जो runtime पर evaluated किए जाते हैं, channel content के आधार पर workflow decisions को drive करने के लिए conditional logic को enable करते हैं।

### 5.1. `.branch()` के साथ Routing

उदाहरण के लिए, चलो pretend करते हैं कि हमारे sequencing samples को FASTP के साथ trim करने की आवश्यकता है केवल तभी जब वे एक निश्चित threshold से ऊपर coverage वाले human samples हों। Mouse samples या low-coverage samples को Trimgalore के साथ run किया जाना चाहिए (यह एक contrived उदाहरण है, लेकिन यह point को illustrate करता है)।

हमने `modules/trimgalore.nf` में एक simple Trimgalore process प्रदान किया है, यदि तुम चाहो तो देखो, लेकिन details इस exercise के लिए महत्वपूर्ण नहीं हैं। मुख्य point यह है कि हम samples को उनके metadata के आधार पर route करना चाहते हैं।

`modules/trimgalore.nf` से नए को include करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... और फिर अपने `main.nf` workflow को modify करो ताकि samples को उनके metadata के आधार पर branch किया जा सके और उन्हें उपयुक्त trimming process के माध्यम से route किया जा सके, इस तरह:

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

इस modified workflow को चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

यहाँ, हमने samples को उनके metadata के आधार पर route करने के लिए `.branch{}` operator के अंदर छोटे लेकिन mighty conditional expressions का उपयोग किया है। High coverage वाले human samples `FASTP` के माध्यम से जाते हैं, जबकि अन्य सभी samples `TRIMGALORE` के माध्यम से जाते हैं।

### 5.2. Truthiness के साथ `.filter()` का उपयोग करना

Workflow execution को control करने के लिए एक और powerful pattern `.filter()` operator है, जो यह निर्धारित करने के लिए एक closure का उपयोग करता है कि कौन से items pipeline के नीचे जारी रहने चाहिए। Filter closure के अंदर, तुम **boolean expressions** लिखोगे जो decide करते हैं कि कौन से items pass होते हैं।

Nextflow (कई dynamic languages की तरह) में **"truthiness"** की एक अवधारणा है जो निर्धारित करती है कि boolean contexts में कौन से values `true` या `false` के रूप में evaluate होते हैं:

- **Truthy**: Non-null values, non-empty strings, non-zero numbers, non-empty collections
- **Falsy**: `null`, empty strings `""`, zero `0`, empty collections `[]` या `[:]`, `false`

इसका मतलब है कि `meta.id` अकेले (explicit `!= null` के बिना) check करता है कि ID मौजूद है और empty नहीं है। चलो हमारी quality requirements को पूरा नहीं करने वाले samples को filter करने के लिए इसका उपयोग करते हैं।

Branch operation से पहले निम्नलिखित जोड़ो:

=== "बाद में"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Filter out invalid or low-quality samples
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

Workflow को फिर से चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

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

क्योंकि हमने एक ऐसा filter चुना है जो कुछ samples को exclude करता है, कम tasks execute किए गए।

Filter expression `meta.id && meta.organism && meta.depth >= 25000000` truthiness को explicit comparisons के साथ combine करता है:

- `meta.id && meta.organism` check करता है कि दोनों fields मौजूद हैं और non-empty हैं (truthiness का उपयोग करते हुए)
- `meta.depth >= 25000000` एक explicit comparison के साथ पर्याप्त sequencing depth सुनिश्चित करता है

!!! note "Practice में Truthiness"

    Expression `meta.id && meta.organism` लिखने से अधिक concise है:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    यह filtering logic को बहुत cleaner और पढ़ने में आसान बनाता है।

### सारांश

इस section में, तुमने Nextflow operators जैसे `.branch{}` और `.filter{}` के closure interfaces का उपयोग करके workflow execution को control करने के लिए conditional logic का उपयोग करना सीखा, concise conditional expressions लिखने के लिए truthiness का लाभ उठाते हुए।

हमारी pipeline अब intelligently samples को उपयुक्त processes के माध्यम से route करती है, लेकिन production workflows को invalid data को gracefully handle करने की आवश्यकता होती है। चलो अपने workflow को missing या null values के खिलाफ robust बनाते हैं।

---

## 6. Safe Navigation और Elvis Operators

हमारा `separateMetadata` function वर्तमान में मानता है कि सभी CSV fields मौजूद और valid हैं। लेकिन incomplete data के साथ क्या होता है? चलो पता लगाते हैं।

### 6.1. समस्या: Properties को Access करना जो मौजूद नहीं हैं

चलो कहते हैं कि हम optional sequencing run information के लिए support जोड़ना चाहते हैं। कुछ labs में, samples में sequencing run ID या batch number के लिए एक अतिरिक्त field हो सकता है, लेकिन हमारे current CSV में यह column नहीं है। चलो वैसे भी इसे access करने की कोशिश करते हैं।

`separateMetadata` function को एक run_id field शामिल करने के लिए modify करो:

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

अब workflow चलाओ:

```bash
nextflow run main.nf
```

??? failure "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

यह NullPointerException के साथ crash होता है।

समस्या यह है कि `row.run_id` `null` return करता है क्योंकि `run_id` column हमारे CSV में मौजूद नहीं है। जब हम `null` पर `.toUpperCase()` को call करने की कोशिश करते हैं, तो यह crash हो जाता है। यहीं पर safe navigation operator दिन बचाता है।

### 6.2. Safe Navigation Operator (`?.`)

Safe navigation operator (`?.`) एक `null` value पर call किए जाने पर exception throw करने के बजाय `null` return करता है। यदि `?.` से पहले object `null` है, तो पूरा expression method को execute किए बिना `null` के रूप में evaluate होता है।

Safe navigation का उपयोग करने के लिए function को update करो:

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

??? success "Command output"

    ```console
    <!-- TODO: output -->
    ```

कोई crash नहीं! Workflow अब missing field को gracefully handle करता है। जब `row.run_id` `null` है, तो `?.` operator `.toUpperCase()` call को prevent करता है, और `run_id` exception cause करने के बजाय `null` बन जाता है।

### 6.3. Defaults के लिए Elvis Operator (`?:`)

Elvis operator (`?:`) default values प्रदान करता है जब left side "falsy" (जैसा कि पहले समझाया गया) होता है। इसका नाम Elvis Presley के नाम पर रखा गया है क्योंकि `?:` sideways से देखने पर उनके प्रसिद्ध बालों और आँखों की तरह दिखता है!

अब जब हम safe navigation का उपयोग कर रहे हैं, तो `run_id` उन samples के लिए `null` होगा जिनमें वह field नहीं है। चलो एक default value प्रदान करने के लिए Elvis operator का उपयोग करते हैं और इसे अपने `sample_meta` map में जोड़ते हैं:

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

Results देखने के लिए workflow में एक `view()` operator भी जोड़ो:

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

और workflow चलाओ:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Perfect! अब सभी samples में एक `run` field है जिसमें या तो उनकी actual run ID (uppercase में) या default value 'UNSPECIFIED' है। `?.` और `?:` का combination दोनों safety (कोई crashes नहीं) और sensible defaults प्रदान करता है।

अब `.view()` operator को हटा दो क्योंकि हमने confirm कर लिया है कि यह काम करता है।

!!! tip "Safe Navigation और Elvis को Combine करना"

    Pattern `value?.method() ?: 'default'` production workflows में सामान्य है:

    - `value?.method()` - Safely method को call करता है, यदि `value` `null` है तो `null` return करता है
    - `?: 'default'` - यदि result `null` है तो fallback प्रदान करता है

    यह pattern missing/incomplete data को gracefully handle करता है।

इन operators का उपयोग functions, operator closures (`.map{}`, `.filter{}`), process scripts, और config files में consistently करो। वे real-world data को handle करते समय crashes को prevent करते हैं।

### सारांश

- **Safe navigation (`?.`)**: Null values पर crashes को prevent करता है - exception throw करने के बजाय null return करता है
- **Elvis operator (`?:`)**: Defaults प्रदान करता है - `value ?: 'default'`
- **Combining**: `value?.method() ?: 'default'` सामान्य pattern है

ये operators workflows को incomplete data के लिए resilient बनाते हैं - real-world काम के लिए आवश्यक।

---

## 7. `error()` और `log.warn` के साथ Validation

कभी-कभी तुम्हें workflow को तुरंत रोकने की आवश्यकता होती है यदि input parameters invalid हैं। Nextflow में, तुम validation logic को implement करने के लिए `error()` और `log.warn` जैसे built-in functions के साथ-साथ `if` statements और boolean logic जैसे standard programming constructs का उपयोग कर सकते हो। चलो अपने workflow में validation जोड़ते हैं।

अपने workflow block से पहले एक validation function बनाओ, इसे workflow से call करो, और CSV file path के लिए एक parameter का उपयोग करने के लिए channel creation को बदलो। यदि parameter missing है या file मौजूद नहीं है, तो एक स्पष्ट message के साथ execution को रोकने के लिए `error()` को call करो।

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Check input parameter is provided
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // Check CSV file exists
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

अब CSV file के बिना चलाने की कोशिश करो:

```bash
nextflow run main.nf
```

??? failure "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Input CSV file path not provided. Please specify --input <file.csv>
    ```

Workflow बाद में mysteriously fail होने के बजाय एक स्पष्ट error message के साथ तुरंत रुक जाता है

अब इसे एक non-existent file के साथ चलाओ:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Input CSV file not found: ./data/nonexistent.csv
    ```

अंत में, इसे सही file के साथ चलाओ:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Command output"

    ```console
    <!-- TODO: output -->
    ```

इस बार यह successfully चलता है।

तुम `separateMetadata` function के भीतर भी validation जोड़ सकते हो। चलो low sequencing depth वाले samples के लिए warnings issue करने के लिए non-fatal `log.warn` का उपयोग करते हैं, लेकिन फिर भी workflow को जारी रखने की अनुमति देते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validate data makes sense
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

Original CSV के साथ workflow को फिर से चलाओ:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Low sequencing depth for sample_002: 25000000
    ```

हम एक sample के लिए low sequencing depth के बारे में एक warning देखते हैं।

### सारांश

- **`error()`**: स्पष्ट message के साथ workflow को तुरंत रोकता है
- **`log.warn`**: Workflow को रोके बिना warnings issue करता है
- **Early validation**: Processing से पहले inputs को check करो ताकि helpful errors के साथ fast fail हो
- **Validation functions**: Reusable validation logic बनाओ जिसे workflow start पर call किया जा सके

Proper validation workflows को अधिक robust और user-friendly बनाता है स्पष्ट error messages के साथ problems को जल्दी catch करके।

---

## 8. Workflow Event Handlers

अब तक, हम अपने workflow scripts और process definitions में code लिख रहे हैं। लेकिन एक और महत्वपूर्ण feature है जिसके बारे में तुम्हें पता होना चाहिए: workflow event handlers।

Event handlers closures हैं जो तुम्हारे workflow के lifecycle में specific points पर run होते हैं। वे logging, notifications, या cleanup operations जोड़ने के लिए perfect हैं। इन handlers को तुम्हारे workflow script में तुम्हारे workflow definition के साथ define किया जाना चाहिए।

### 8.1. `onComplete` Handler

सबसे अधिक उपयोग किया जाने वाला event handler `onComplete` है, जो तब run होता है जब तुम्हारा workflow finish होता है (चाहे यह succeed हो या fail)। चलो हमारे pipeline results को summarize करने के लिए एक जोड़ते हैं।

अपनी `main.nf` file में event handler जोड़ो, अपनी workflow definition के अंदर:

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

यह closure तब run होता है जब workflow complete होता है। अंदर, तुम्हारे पास `workflow` object तक access है जो execution के बारे में उपयोगी properties प्रदान करता है।

अपना workflow चलाओ और तुम यह summary अंत में दिखाई देगी!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Command output"

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

चलो conditional logic जोड़कर इसे और उपयोगी बनाते हैं:

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

अब हमें एक और भी informative summary मिलता है, जिसमें एक success/failure message और यदि निर्दिष्ट किया गया हो तो output directory शामिल है:

<!-- TODO: add run command -->

??? success "Command output"

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

```groovy title="main.nf - Writing summary to file"
workflow {
    // ... your workflow code ...

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

        // Write to a log file
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. `onError` Handler

`onComplete` के अलावा, एक और event handler है जिसका तुम उपयोग कर सकते हो: `onError`, जो केवल तभी run होता है जब workflow fail होता है:

```groovy title="main.nf - onError handler"
workflow {
    // ... your workflow code ...

    workflow.onError = {
        println "="* 50
        println "Pipeline execution failed!"
        println "Error message: ${workflow.errorMessage}"
        println "="* 50

        // Write detailed error log
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

तुम अपने workflow script में कई handlers को एक साथ उपयोग कर सकते हो:

```groovy title="main.nf - Combined handlers"
workflow {
    // ... your workflow code ...

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

- **Event handler closures**: तुम्हारे workflow script में closures जो विभिन्न lifecycle points पर run होते हैं
- **`onComplete` handler**: Execution summaries और result reporting के लिए
- **`onError` handler**: Error handling और failures को log करने के लिए
- **Workflow object properties**: `workflow.success`, `workflow.duration`, `workflow.errorMessage`, आदि को access करना

Event handlers दिखाते हैं कि तुम sophisticated logging और notification capabilities जोड़ने के लिए अपने workflow scripts के भीतर Nextflow भाषा की पूरी शक्ति का उपयोग कैसे कर सकते हो।

---

## सारांश

बधाई हो, तुमने इसे पूरा कर लिया!

इस side quest के दौरान, तुमने एक comprehensive sample processing pipeline बनाया है जो बुनियादी metadata handling से एक sophisticated, production-ready workflow में विकसित हुआ।
प्रत्येक section पिछले एक पर बना, यह demonstrate करते हुए कि programming constructs कैसे simple workflows को powerful data processing systems में transform करते हैं, निम्नलिखित लाभों के साथ:

- **Clearer code**: Dataflow vs scripting को समझना तुम्हें अधिक organized workflows लिखने में मदद करता है
- **Robust handling**: Safe navigation और Elvis operators workflows को missing data के लिए resilient बनाते हैं
- **Flexible processing**: Conditional logic तुम्हारे workflows को विभिन्न sample types को appropriately process करने देता है
- **Adaptive resources**: Dynamic directives input विशेषताओं के आधार पर resource usage को optimize करते हैं

यह progression real-world bioinformatics pipelines के evolution को mirror करता है, कुछ samples को handle करने वाले research prototypes से लेकर laboratories और institutions में हजारों samples को process करने वाली production systems तक।
हर challenge जिसे तुमने solve किया और हर pattern जो तुमने सीखा वह actual problems को reflect करता है जिनका developers Nextflow workflows को scale करते समय सामना करते हैं।

अपने खुद के काम में इन patterns को apply करना तुम्हें robust, production-ready workflows बनाने में सक्षम बनाएगा।

### मुख्य patterns

1.  **Dataflow vs Scripting:** तुमने dataflow operations (channel orchestration) और scripting (code जो data को manipulate करता है) के बीच अंतर करना सीखा, जिसमें विभिन्न types पर operations के बीच महत्वपूर्ण अंतर शामिल हैं जैसे Channel vs List पर `collect`।

    - Dataflow: channel orchestration

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: collections पर data processing

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Advanced String Processing**: तुमने file names को parse करने के लिए regular expressions, processes में dynamic script generation, और variable interpolation (Nextflow vs Bash vs Shell) में महारत हासिल की।

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

    - File collection को command arguments में (process script block में)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Reusable Functions बनाना**: तुमने complex logic को named functions में extract करना सीखा जिन्हें channel operators से call किया जा सकता है, workflows को अधिक readable और maintainable बनाते हुए।

    - एक named function define करो

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* code hidden for brevity */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* code hidden for brevity */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - एक workflow में named function को call करो

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Closures के साथ Dynamic Resource Directives**: तुमने input विशेषताओं के आधार पर adaptive resource allocation के लिए process directives में closures का उपयोग करना explored किया।

    - Named closures और composition

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Scope access के साथ closures

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Conditional Logic और Process Control**: तुमने `.branch()` और `.filter()` operators का उपयोग करके intelligent routing जोड़ा, concise conditional expressions के लिए truthiness का लाभ उठाते हुए।

    - विभिन्न workflow branches के माध्यम से data को route करने के लिए `.branch()` का उपयोग करो

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

    - 'Truthiness' के साथ data को subset करने के लिए `filter()` का उपयोग करो

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Safe Navigation और Elvis Operators**: तुमने null-safe property access के लिए `?.` और default values प्रदान करने के लिए `?:` का उपयोग करके pipeline को missing data के खिलाफ robust बनाया।

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **error() और log.warn के साथ Validation**: तुमने inputs को जल्दी validate करना और स्पष्ट error messages के साथ fast fail करना सीखा।

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Configuration Event Handlers**: तुमने logging, notifications, और lifecycle management के लिए workflow event handlers (`onComplete` और `onError`) का उपयोग करना सीखा।

    - Log और notify करने के लिए `onComplete` का उपयोग करना

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

    - Failure के मामले में विशेष रूप से action लेने के लिए `onError` का उपयोग करना

    ```groovy
    workflow.onError = {
        // Write detailed error log
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
    ```

### अतिरिक्त resources

- [Nextflow Language Reference](https://nextflow.io/docs/latest/reference/syntax.html)
- [Nextflow Operators](https://www.nextflow.io/docs/latest/operator.html)
- [Nextflow Script Syntax](https://www.nextflow.io/docs/latest/script.html)
- [Nextflow Standard Library](https://nextflow.io/docs/latest/reference/stdlib.html)

जब तुम्हें अधिक advanced features को explore करने की आवश्यकता हो तो इन resources को check करना सुनिश्चित करो।

तुम्हें अपने skills को practice करने और expand करने से लाभ होगा ताकि:

- Dataflow और scripting के बीच उचित separation के साथ cleaner workflows लिखो
- Nextflow, Bash, और shell variables के साथ सामान्य pitfalls से बचने के लिए variable interpolation में महारत हासिल करो
- Efficient, adaptive workflows के लिए dynamic resource directives का उपयोग करो
- File collections को properly formatted command-line arguments में transform करो
- Regex और string processing का उपयोग करके विभिन्न file naming conventions और input formats को gracefully handle करो
- Advanced closure patterns और functional programming का उपयोग करके reusable, maintainable code बनाओ
- Collection operations का उपयोग करके complex datasets को process और organize करो
- अपने workflows को production-ready बनाने के लिए validation, error handling, और logging जोड़ो
- Event handlers के साथ workflow lifecycle management को implement करो

---

## आगे क्या है?

[Side Quests के menu](./index.md) पर वापस जाओ या list में अगले topic पर जाने के लिए page के नीचे दाईं ओर button पर click करो।
