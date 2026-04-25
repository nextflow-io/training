# मेटाडेटा और मेटा मैप्स

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

किसी भी वैज्ञानिक विश्लेषण में, हम केवल raw डेटा फ़ाइलों के साथ काम नहीं करते।
हर फ़ाइल के साथ अपनी अतिरिक्त जानकारी होती है: वह क्या है, कहाँ से आई है, और उसे क्या खास बनाता है।
इस अतिरिक्त जानकारी को हम मेटाडेटा कहते हैं।

मेटाडेटा वह डेटा है जो दूसरे डेटा का वर्णन करता है।
मेटाडेटा फ़ाइलों और प्रयोगात्मक परिस्थितियों के बारे में महत्वपूर्ण विवरण ट्रैक करता है, और प्रत्येक डेटासेट की अनूठी विशेषताओं के अनुसार विश्लेषण को अनुकूलित करने में मदद करता है।

इसे एक लाइब्रेरी कैटलॉग की तरह सोचो: जबकि किताबों में वास्तविक सामग्री (raw डेटा) होती है, कैटलॉग कार्ड प्रत्येक किताब के बारे में आवश्यक जानकारी प्रदान करते हैं—कब प्रकाशित हुई, किसने लिखी, कहाँ मिलेगी (मेटाडेटा)।
Nextflow पाइपलाइनों में, मेटाडेटा का उपयोग इन कामों के लिए किया जा सकता है:

- पूरे वर्कफ़्लो में फ़ाइल-विशिष्ट जानकारी ट्रैक करना
- फ़ाइल की विशेषताओं के आधार पर प्रोसेस को कॉन्फ़िगर करना
- संयुक्त विश्लेषण के लिए संबंधित फ़ाइलों को ग्रुप करना

### सीखने के लक्ष्य

इस side quest में, हम वर्कफ़्लो में मेटाडेटा को संभालने का तरीका सीखेंगे।
एक सरल datasheet (जिसे bioinformatics में अक्सर samplesheet कहा जाता है) से शुरू करते हुए जिसमें बुनियादी फ़ाइल जानकारी होती है, तुम सीखोगे कि कैसे:

- CSV फ़ाइलों से फ़ाइल मेटाडेटा पढ़ें और parse करें
- मेटाडेटा मैप्स बनाएं और उनमें बदलाव करें
- वर्कफ़्लो execution के दौरान नए मेटाडेटा फ़ील्ड जोड़ें
- प्रोसेस के व्यवहार को कस्टमाइज़ करने के लिए मेटाडेटा का उपयोग करें

ये कौशल तुम्हें अधिक मज़बूत और लचीली पाइपलाइनें बनाने में मदद करेंगे जो जटिल फ़ाइल संबंधों और प्रोसेसिंग आवश्यकताओं को संभाल सकती हैं।

### पूर्वापेक्षाएँ

इस side quest को शुरू करने से पहले, तुम्हें:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष शुरुआती कोर्स पूरा करना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (प्रोसेस, चैनल, ऑपरेटर) का उपयोग करने में सहज होना चाहिए।

---

## 0. शुरू करना

#### ट्रेनिंग codespace खोलो

अगर तुमने अभी तक नहीं किया है, तो [Environment Setup](../envsetup/index.md) में बताए अनुसार ट्रेनिंग वातावरण खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाओ

चलो उस डायरेक्टरी में जाते हैं जहाँ इस ट्यूटोरियल की फ़ाइलें हैं।

```bash
cd side-quests/metadata
```

तुम VSCode को इस डायरेक्टरी पर फ़ोकस करने के लिए सेट कर सकते हो:

```bash
code .
```

#### सामग्री की समीक्षा करो

तुम्हें एक main workflow फ़ाइल और एक `data` डायरेक्टरी मिलेगी जिसमें एक datasheet और कुछ डेटा फ़ाइलें हैं।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    .
    ├── data
    │   ├── bonjour.txt
    │   ├── ciao.txt
    │   ├── guten_tag.txt
    │   ├── hallo.txt
    │   ├── hello.txt
    │   ├── hola.txt
    │   ├── salut.txt
    │   └── datasheet.csv
    ├── main.nf
    └── nextflow.config
    ```

`main.nf` फ़ाइल में वर्कफ़्लो एक stub है जिसे तुम धीरे-धीरे एक पूरी तरह से काम करने वाले वर्कफ़्लो में विस्तारित करोगे।

datasheet डेटा फ़ाइलों के पाथ और कुछ संबंधित मेटाडेटा सूचीबद्ध करती है, जो 3 कॉलम में व्यवस्थित है:

- `id`: स्व-व्याख्यात्मक, फ़ाइल को दिया गया एक ID
- `character`: एक character का नाम, जिसका उपयोग हम बाद में अलग-अलग प्राणी बनाने के लिए करेंगे
- `data`: `.txt` फ़ाइलों के पाथ जिनमें विभिन्न भाषाओं में अभिवादन हैं

```console title="datasheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

प्रत्येक डेटा फ़ाइल में पाँच भाषाओं में से एक में कुछ अभिवादन टेक्स्ट है (fr: फ्रेंच, de: जर्मन, es: स्पेनिश, it: इतालवी, en: अंग्रेज़ी)।

हम तुम्हें `langid` नामक एक containerized भाषा विश्लेषण टूल भी प्रदान करेंगे।

#### असाइनमेंट की समीक्षा करो

तुम्हारी चुनौती एक Nextflow वर्कफ़्लो लिखना है जो:

1. प्रत्येक फ़ाइल में भाषा को स्वचालित रूप से **पहचाने**
2. फ़ाइलों को भाषा परिवार के अनुसार **ग्रुप करे** (Germanic बनाम Romance भाषाएँ)
3. प्रत्येक फ़ाइल की भाषा और मेटाडेटा के आधार पर प्रोसेसिंग को **कस्टमाइज़ करे**
4. आउटपुट को भाषा समूह के अनुसार **व्यवस्थित करे**

यह एक विशिष्ट वर्कफ़्लो पैटर्न है जहाँ फ़ाइल-विशिष्ट मेटाडेटा प्रोसेसिंग निर्णयों को नियंत्रित करता है; बिल्कुल वही समस्या जिसे मेटा मैप्स सुंदर तरीके से हल करते हैं।

#### तैयारी की जाँच सूची

क्या तुम शुरू करने के लिए तैयार हो?

- [ ] मैं इस कोर्स के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता/समझती हूँ
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी working directory उचित रूप से सेट की है
- [ ] मैं असाइनमेंट समझता/समझती हूँ

अगर तुम सभी बॉक्स चेक कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. Datasheet से मेटाडेटा लोड करना

`main.nf` वर्कफ़्लो फ़ाइल खोलो और उस वर्कफ़्लो stub की जाँच करो जो हम तुम्हें शुरुआती बिंदु के रूप में दे रहे हैं।

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

तुम देख सकते हो कि हमने example datasheet को एक फ़ाइल के रूप में लोड करने के लिए एक बुनियादी channel factory सेट किया है, लेकिन यह अभी तक फ़ाइल की सामग्री नहीं पढ़ेगा।
चलो इसे जोड़ने से शुरू करते हैं।

### 1.1. `splitCsv` से सामग्री पढ़ना

हमें एक ऐसा ऑपरेटर चुनना होगा जो हमारी ओर से न्यूनतम प्रयास के साथ फ़ाइल की सामग्री को उचित रूप से parse करे।
चूँकि हमारी datasheet CSV फॉर्मेट में है, यह [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) ऑपरेटर का काम है, जो फ़ाइल की प्रत्येक पंक्ति को चैनल में एक element के रूप में लोड करता है।

channel construction कोड में `splitCsv()` ऑपरेशन जोड़ने के लिए निम्नलिखित बदलाव करो, साथ ही एक `view()` ऑपरेशन भी जोड़ो ताकि यह जाँच सको कि फ़ाइल की सामग्री चैनल में सही तरीके से लोड हो रही है।

=== "बाद में"

    ```groovy title="main.nf" linenums="3" hl_lines="4-5"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")

    }
    ```

ध्यान दो कि हम `header: true` विकल्प का उपयोग कर रहे हैं ताकि Nextflow को CSV फ़ाइल की पहली पंक्ति को header पंक्ति के रूप में पढ़ने के लिए कहा जा सके।

चलो देखते हैं इससे क्या निकलता है, है ना?
वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

हम देख सकते हैं कि ऑपरेटर ने CSV फ़ाइल की प्रत्येक पंक्ति के लिए key-value pairs का एक map बनाया है, जिसमें column headers संबंधित values की keys हैं।

प्रत्येक map entry हमारी datasheet के एक कॉलम से मेल खाती है:

- `id`
- `character`
- `recording`

यह बहुत अच्छा है! इससे प्रत्येक फ़ाइल से विशिष्ट फ़ील्ड तक पहुँचना आसान हो जाता है।
उदाहरण के लिए, हम `id` से फ़ाइल ID और `recording` से txt फ़ाइल पाथ तक पहुँच सकते हैं।

??? info "(वैकल्पिक) मैप्स के बारे में अधिक जानकारी"

    Groovy में, जो प्रोग्रामिंग भाषा Nextflow पर आधारित है, एक map एक key-value डेटा संरचना है जो Python में dictionaries, JavaScript में objects, या Ruby में hashes के समान है।

    यहाँ एक चलाने योग्य script है जो दिखाती है कि तुम व्यवहार में एक map कैसे define कर सकते हो और उसकी सामग्री तक कैसे पहुँच सकते हो:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // एक सरल map बनाएं
    def my_map = [id:'sampleA', character:'squirrel']

    // पूरा map प्रिंट करें
    println "map: ${my_map}"

    // dot notation का उपयोग करके अलग-अलग values तक पहुँचें
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    भले ही इसमें एक उचित `workflow` block नहीं है, Nextflow इसे एक वर्कफ़्लो की तरह चला सकता है:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    और यहाँ है जो तुम आउटपुट में देखने की उम्मीद कर सकते हो:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. `map` से विशिष्ट फ़ील्ड चुनना

मान लो हम datasheet से `character` कॉलम तक पहुँचना और उसे प्रिंट करना चाहते हैं।
हम Nextflow `map` ऑपरेटर का उपयोग करके अपने चैनल में प्रत्येक item पर iterate कर सकते हैं और map object से विशेष रूप से `character` entry चुन सकते हैं।

वर्कफ़्लो में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="3" hl_lines="5-7"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()

    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

अब वर्कफ़्लो फिर से चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

सफलता! हमने अपनी datasheet से प्राप्त map संरचना का लाभ उठाकर प्रत्येक पंक्ति के लिए अलग-अलग कॉलम से values तक पहुँचा।

अब जब हमने datasheet को सफलतापूर्वक पढ़ लिया है और प्रत्येक पंक्ति के डेटा तक पहुँच है, तो हम अपनी पाइपलाइन logic को लागू करना शुरू कर सकते हैं।

### 1.3. मेटाडेटा को 'meta map' में व्यवस्थित करना

वर्कफ़्लो की वर्तमान स्थिति में, इनपुट फ़ाइलें (`recording` key के अंतर्गत) और संबंधित मेटाडेटा (`id`, `character`) सभी एक ही स्तर पर हैं, जैसे वे सब एक बड़े थैले में हों।
इसका व्यावहारिक परिणाम यह है कि इस चैनल का उपयोग करने वाले प्रत्येक प्रोसेस को इस संरचना को ध्यान में रखकर कॉन्फ़िगर करना होगा:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

यह तब तक ठीक है जब तक datasheet में कॉलम की संख्या नहीं बदलती।
हालाँकि, अगर तुम datasheet में एक भी कॉलम जोड़ते हो, तो चैनल का आकार प्रोसेस की अपेक्षाओं से मेल नहीं खाएगा, और वर्कफ़्लो में त्रुटियाँ आएंगी।
इससे प्रोसेस को दूसरों के साथ share करना भी मुश्किल हो जाता है जिनके पास थोड़ा अलग इनपुट डेटा हो सकता है, और तुम्हें script block में ऐसे variables hard-code करने पड़ सकते हैं जिनकी ज़रूरत नहीं है।

इस समस्या से बचने के लिए, हमें tuple के भीतर एक item में सभी मेटाडेटा को इकट्ठा करने का तरीका खोजना होगा, जिसे हम metadata map, या अधिक सरल रूप से 'meta map' कहेंगे।

`map` ऑपरेशन में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

हमने अपने चैनल elements को दो elements वाले tuple में पुनर्गठित किया है: meta map और संबंधित फ़ाइल object।

वर्कफ़्लो चलाते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

अब, चैनल में प्रत्येक element में पहले metadata map और दूसरे स्थान पर संबंधित फ़ाइल object है:

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

परिणामस्वरूप, datasheet में अधिक कॉलम जोड़ने से `meta` map में अधिक मेटाडेटा उपलब्ध होगा, लेकिन चैनल का आकार नहीं बदलेगा।
इससे हम ऐसे प्रोसेस लिख सकते हैं जो चैनल का उपयोग करते हैं बिना input specification में मेटाडेटा items को hard-code किए:

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

यह Nextflow वर्कफ़्लो में मेटाडेटा व्यवस्थित करने के लिए एक व्यापक रूप से उपयोग किया जाने वाला convention है।

### सारांश

इस section में, तुमने सीखा:

- **मेटाडेटा क्यों महत्वपूर्ण है:** अपने डेटा के साथ मेटाडेटा रखने से पूरे वर्कफ़्लो में महत्वपूर्ण फ़ाइल जानकारी सुरक्षित रहती है।
- **Datasheet कैसे पढ़ें:** header जानकारी के साथ CSV फ़ाइलें पढ़ने और पंक्तियों को संरचित डेटा में बदलने के लिए `splitCsv` का उपयोग करना।
- **Meta map कैसे बनाएं:** tuple संरचना `[ [id:value, ...], file ]` का उपयोग करके मेटाडेटा को फ़ाइल डेटा से अलग करना।

---

## 2. मेटाडेटा में बदलाव करना

अब जब हमारा मेटाडेटा लोड हो गया है, तो चलो इसके साथ कुछ करते हैं!

हम [`langid`](https://github.com/saffsd/langid.py) नामक एक टूल का उपयोग करेंगे ताकि प्रत्येक creature की recording फ़ाइल में मौजूद भाषा की पहचान की जा सके।
यह टूल भाषाओं के एक सेट पर pre-trained है, और टेक्स्ट का एक अंश दिए जाने पर, यह एक भाषा prediction और एक संबंधित probability score दोनों को `stdout` पर आउटपुट करेगा।

### 2.1. प्रोसेस import करना और कोड की जाँच करना

हम तुम्हें `IDENTIFY_LANGUAGE` नामक एक pre-written process module प्रदान करते हैं जो `langid` टूल को wrap करता है, इसलिए तुम्हें बस workflow block से पहले एक include statement जोड़ना है।

वर्कफ़्लो में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

तुम module फ़ाइल खोलकर उसके कोड की जाँच कर सकते हो:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// प्रत्येक इनपुट फ़ाइल की भाषा predict करने के लिए langid का उपयोग करें
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}
```

जैसा कि तुम देख सकते हो, input definition वही `tuple val(meta), path(file)` संरचना का उपयोग करती है जो हमने अभी अपने इनपुट चैनल पर लागू की है।

output definition एक tuple के रूप में संरचित है जो input की तरह है, सिवाय इसके कि इसमें तीसरे element के रूप में `stdout` भी है।
यह `tuple val(meta), path(file), <output>` पैटर्न मेटाडेटा को इनपुट डेटा और आउटपुट दोनों के साथ जोड़े रखता है जैसे-जैसे यह पाइपलाइन से गुज़रता है।

ध्यान दो कि हम यहाँ Nextflow के [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) output qualifier का उपयोग कर रहे हैं क्योंकि टूल अपना आउटपुट सीधे console पर प्रिंट करता है न कि किसी फ़ाइल में; और हम command line में `sed` का उपयोग probability score हटाने, newline characters हटाकर string को साफ करने, और केवल भाषा prediction वापस करने के लिए करते हैं।

### 2.2. `IDENTIFY_LANGUAGE` को call करना

अब जब प्रोसेस वर्कफ़्लो के लिए उपलब्ध है, तो हम data channel पर इसे चलाने के लिए `IDENTIFY_LANGUAGE` प्रोसेस को call कर सकते हैं।

वर्कफ़्लो में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएं
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

ध्यान दो कि हमने channel construction में मूल `.view()` ऑपरेशन हटा दिया है।

अब हम वर्कफ़्लो चला सकते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (7)
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

बढ़िया! अब हमारे पास prediction है कि प्रत्येक character कौन सी भाषा बोलता है।

और जैसा कि पहले उल्लेख किया गया था, हमने आउटपुट में इनपुट फ़ाइल और meta map भी शामिल किया है, जिसका अर्थ है कि दोनों हमारे द्वारा अभी उत्पन्न नई जानकारी के साथ जुड़े रहते हैं।
यह अगले चरण में उपयोगी साबित होगा।

!!! note "नोट"

    अधिक सामान्य रूप से, meta map को results के साथ जोड़े रखने का यह पैटर्न उन संबंधित results को जोड़ना आसान बनाता है जो समान identifiers share करते हैं।

    जैसा कि तुम पहले ही सीख चुके होगे, तुम results को उनके पार match करने के लिए चैनलों में items के क्रम पर भरोसा नहीं कर सकते।
    इसके बजाय, तुम्हें डेटा को सही तरीके से जोड़ने के लिए keys का उपयोग करना होगा, और meta maps इस उद्देश्य के लिए एक आदर्श संरचना प्रदान करते हैं।

    हम इस use case को [Splitting & Grouping](../splitting_and_grouping/) side quest में विस्तार से explore करते हैं।

### 2.3. प्रोसेस आउटपुट से मेटाडेटा बढ़ाना

चूँकि हमने जो results अभी उत्पन्न किए हैं वे स्वयं फ़ाइलों की सामग्री के बारे में मेटाडेटा का एक रूप हैं, इसलिए उन्हें हमारे meta map में जोड़ना उपयोगी होगा।

हालाँकि, हम मौजूदा meta map को in-place modify नहीं करना चाहते।
तकनीकी दृष्टिकोण से, ऐसा करना _संभव_ है, लेकिन यह असुरक्षित है।

इसलिए, हम `+` ऑपरेटर (एक Groovy feature) का उपयोग करके एक नया meta map बनाएंगे जिसमें मौजूदा meta map की सामग्री के साथ-साथ नई जानकारी रखने वाला एक नया `lang: lang_id` key-value pair होगा।
और हम इसे एक [`map`](https://www.nextflow.io/docs/latest/operator.html#map) ऑपरेशन के साथ जोड़ेंगे ताकि पुराने map को नए से बदला जा सके।

वर्कफ़्लो में ये बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएं
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएं
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

अगर तुम अभी तक `+` ऑपरेटर से परिचित नहीं हो, या अगर यह भ्रामक लगता है, तो नीचे दिए गए विस्तृत स्पष्टीकरण को पढ़ने में कुछ मिनट लगाओ।

??? info "`+` ऑपरेटर का उपयोग करके नया meta map बनाना"

    **पहले, तुम्हें यह जानना होगा कि हम Groovy ऑपरेटर `+` का उपयोग करके दो maps की सामग्री को merge कर सकते हैं।**

    मान लो हमारे पास निम्नलिखित maps हैं:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    हम उन्हें इस तरह merge कर सकते हैं:

    ```groovy
    new_map = map1 + map2
    ```

    `new_map` की सामग्री होगी:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    बढ़िया!

    **लेकिन क्या होगा अगर तुम्हें एक ऐसा field जोड़ना हो जो पहले से किसी map का हिस्सा नहीं है?**

    मान लो तुम फिर से `map1` से शुरू करते हो, लेकिन भाषा prediction अपने map में नहीं है (कोई `map2` नहीं है)।
    इसके बजाय, यह `lang_id` नामक एक variable में है, और तुम जानते हो कि तुम इसकी value (`'fr'`) को `lang` key के साथ store करना चाहते हो।

    तुम वास्तव में निम्नलिखित कर सकते हो:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    यहाँ, `[lang: new_info]` on the fly एक नया unnamed map बनाता है, और `map1 + ` `map1` को नए unnamed map के साथ merge करता है, जो पहले जैसी ही `new_map` सामग्री उत्पन्न करता है।

    अच्छा है, है ना?

    **अब इसे Nextflow `channel.map()` ऑपरेशन के संदर्भ में transpose करते हैं।**

    कोड बन जाता है:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    यह निम्नलिखित करता है:

    - `map1, lang_id ->` tuple में दो items लेता है
    - `[map1 + [lang: lang_id]]` ऊपर बताए अनुसार नया map बनाता है

    आउटपुट हमारे उदाहरण में `new_map` जैसी ही सामग्री वाला एक single unnamed map है।
    इसलिए हमने effectively transform किया है:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    को:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    उम्मीद है कि तुम देख सकते हो कि अगर हम `map1` को `meta` में बदलते हैं, तो यह मूल रूप से वह सब है जो हमें अपने वर्कफ़्लो में meta map में भाषा prediction जोड़ने के लिए चाहिए।

    सिवाय एक चीज़ के!

    हमारे वर्कफ़्लो के मामले में, **हमें tuple में `file` object की उपस्थिति का भी ध्यान रखना होगा**, जो `meta, file, lang_id` से बना है।

    तो यहाँ कोड बन जाएगा:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    अगर तुम्हें यह समझने में कठिनाई हो रही है कि `map` ऑपरेशन में `file` क्यों इधर-उधर होता दिखता है, तो कल्पना करो कि `[meta + [lang: lang_id], file]` के बजाय वह लाइन `[new_map, file]` पढ़ती है।
    इससे यह स्पष्ट होना चाहिए कि हम बस `file` को tuple में दूसरे स्थान पर उसकी मूल जगह पर छोड़ रहे हैं। हमने बस `new_info` value ली और उसे पहले स्थान पर मौजूद map में fold कर दिया।

    **और यह हमें वापस `tuple val(meta), path(file)` चैनल संरचना पर लाता है!**

एक बार जब तुम्हें विश्वास हो जाए कि तुम समझ गए हो यह कोड क्या कर रहा है, तो वर्कफ़्लो चलाओ और देखो कि यह काम किया या नहीं:

```bash
nextflow run main.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
    [[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
    [[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
    [[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
    [[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
    [[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
    ```

हाँ, यह सही है!
हमने प्रोसेस के आउटपुट को `meta, file, lang_id` से सुंदर तरीके से पुनर्गठित किया है ताकि `lang_id` अब meta map में keys में से एक हो, और चैनल के tuples फिर से `meta, file` मॉडल में फिट हों।

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. Conditionals का उपयोग करके भाषा समूह असाइन करना

अब जब हमारे पास भाषा predictions हैं, तो चलो इस जानकारी का उपयोग कुछ नए groupings असाइन करने के लिए करते हैं।

हमारे उदाहरण डेटा में, हमारे characters द्वारा उपयोग की जाने वाली भाषाओं को Germanic भाषाओं (अंग्रेज़ी, जर्मन) और Romance भाषाओं (फ्रेंच, स्पेनिश, इतालवी) में ग्रुप किया जा सकता है।
पाइपलाइन में बाद में उस classification को आसानी से उपलब्ध रखना उपयोगी हो सकता है, तो चलो उस जानकारी को meta map में जोड़ते हैं।

और, अच्छी खबर, यह एक और मामला है जो `map` ऑपरेटर का उपयोग करने के लिए बिल्कुल उपयुक्त है!

विशेष रूप से, हम `lang_group` नामक एक variable define करेंगे, और प्रत्येक डेटा के लिए `lang_group` को क्या value असाइन करनी है यह निर्धारित करने के लिए कुछ सरल conditional logic का उपयोग करेंगे।

सामान्य syntax इस तरह दिखेगी:

```groovy
.map { meta, file ->

    // lang_group define करने वाला conditional logic यहाँ जाता है

    [meta + [lang_group: lang_group], file]
}
```

तुम देख सकते हो यह पिछले चरण में उपयोग किए गए on-the-fly map merging ऑपरेशन के बहुत समान है।
हमें बस conditional statements लिखने हैं।

यहाँ वह conditional logic है जिसे हम लागू करना चाहते हैं:

- `'unknown'` default value के साथ `lang_group` नामक एक variable define करो।
- अगर `lang` German (`'de'`) या English (`'en'`) है, तो `lang_group` को `germanic` में बदलो।
- Else if `lang` French (`'fr'`), Spanish (`'es'`) और Italian (`'it'`) वाली list में शामिल है, तो `lang_group` को `romance` में बदलो।

अगर तुम पहले से Nextflow में conditional statements लिखना जानते हो तो खुद इसे लिखने की कोशिश करो।

!!! tip "सुझाव"

    तुम map ऑपरेशन के भीतर `meta.lang` से `lang` की value तक पहुँच सकते हो।

तुम्हें वर्कफ़्लो में निम्नलिखित बदलाव करने चाहिए:

=== "बाद में"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएं
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="13" hl_lines="7"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएं
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

यहाँ मुख्य बिंदु हैं:

- हम `def lang_group = "unknown"` का उपयोग `lang_group` variable को `unknown` default value के साथ बनाने के लिए करते हैं।
- हम conditional logic के लिए `if {} else if {}` संरचना का उपयोग करते हैं, दो Germanic भाषाओं के लिए वैकल्पिक `.equals()` tests के साथ, और तीन Romance भाषाओं के लिए list में existence की जाँच के साथ।
- हम पहले की तरह updated meta map generate करने के लिए `meta + [lang_group:lang_group]` merge ऑपरेशन का उपयोग करते हैं।

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

एक बार जब यह सब समझ में आ जाए, तो परिणाम देखने के लिए वर्कफ़्लो फिर से चलाओ:

```bash
nextflow run main.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

जैसा कि तुम देख सकते हो, चैनल elements अपनी `[meta, file]` संरचना बनाए रखते हैं, लेकिन meta map में अब यह नई classification शामिल है।

### सारांश

इस section में, तुमने सीखा कि कैसे:

- **इनपुट मेटाडेटा को आउटपुट चैनलों पर लागू करें**: इस तरह मेटाडेटा copy करने से हम बाद में मेटाडेटा सामग्री के आधार पर results को जोड़ सकते हैं।
- **Custom keys बनाएं**: तुमने अपने meta map में दो नई keys बनाईं, उन्हें `meta + [new_key:value]` से मौजूदा meta map में merge किया। एक प्रोसेस से computed value पर आधारित, और एक `map` ऑपरेटर में तुम्हारे द्वारा सेट की गई condition पर आधारित।

ये तुम्हें अपनी पाइपलाइन में आगे बढ़ते हुए फ़ाइलों के साथ नए और मौजूदा मेटाडेटा को जोड़ने की अनुमति देते हैं।
भले ही तुम किसी प्रोसेस के हिस्से के रूप में मेटाडेटा का उपयोग नहीं कर रहे हो, इस तरह meta map को डेटा के साथ जोड़े रखने से सभी प्रासंगिक जानकारी को एक साथ रखना आसान हो जाता है।

---

## 3. प्रोसेस में meta map जानकारी का उपयोग करना

अब जब तुम जानते हो कि meta map कैसे बनाएं और update करें, तो हम वास्तव में मज़ेदार हिस्से पर आ सकते हैं: प्रोसेस में मेटाडेटा का वास्तव में उपयोग करना।

विशेष रूप से, हम अपने वर्कफ़्लो में एक दूसरा चरण जोड़ने जा रहे हैं ताकि प्रत्येक जानवर को ASCII art के रूप में बनाया जा सके और उसे speech bubble में recorded text बोलते हुए दिखाया जा सके।
हम इसे [`cowpy`](https://github.com/jeffbuttars/cowpy) नामक एक टूल का उपयोग करके करेंगे।

??? info "`cowpy` क्या करता है?"

    `cowpy` एक command-line टूल है जो arbitrary text inputs को मज़ेदार तरीके से प्रदर्शित करने के लिए ASCII art generate करता है।
    यह Tony Monroe के classic [cowsay](https://en.wikipedia.org/wiki/Cowsay) टूल का python implementation है।

    ```console
    cowpy "Hello Nextflow"
    ```

    ```console
    ______________________________________________________
    < Hello Nextflow >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

    वैकल्पिक रूप से, तुम default cow के बजाय उपयोग करने के लिए एक character (या 'cowacter') चुन सकते हो।

    ```console
    cowpy "Hello Nextflow" -c tux
    ```

    ```console
    __________________
    < Hello Nextflow >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

अगर तुमने Hello Nextflow कोर्स पूरा किया है, तो तुमने इस टूल को पहले से काम करते देखा है।
अगर नहीं, तो चिंता मत करो; हम जैसे-जैसे आगे बढ़ेंगे, तुम्हें जो कुछ जानना है वह सब cover करेंगे।

### 3.1. प्रोसेस import करना और कोड की जाँच करना

हम तुम्हें `COWPY` नामक एक pre-written process module प्रदान करते हैं जो `cowpy` टूल को wrap करता है, इसलिए तुम्हें बस workflow block से पहले एक include statement जोड़ना है।

वर्कफ़्लो में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="4"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
    include { COWPY } from './modules/cowpy.nf'

    workflow {
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

तुम module फ़ाइल खोलकर उसके कोड की जाँच कर सकते हो:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// cowpy से ASCII art generate करें
process COWPY {

    publishDir "results/", mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
}
```

जैसा कि तुम देख सकते हो, यह प्रोसेस वर्तमान में एक इनपुट फ़ाइल (प्रदर्शित किए जाने वाले टेक्स्ट वाली) और ASCII art में बनाए जाने वाले character को specify करने वाली एक value लेने के लिए designed है, जो आमतौर पर workflow level पर command-line parameter द्वारा प्रदान की जाती है।

### 3.2. Meta map field को input के रूप में pass करना

जब हमने Hello Nextflow कोर्स में `cowpy` टूल का उपयोग किया था, तो हमने final image बनाने के लिए किस character का उपयोग करना है यह निर्धारित करने के लिए एक command-line parameter का उपयोग किया था।
यह समझ में आता था, क्योंकि हम पाइपलाइन के प्रत्येक run में केवल एक image generate कर रहे थे।

हालाँकि, इस ट्यूटोरियल में, हम प्रत्येक subject के लिए एक उपयुक्त image generate करना चाहते हैं जिसे हम process कर रहे हैं, इसलिए command-line parameter का उपयोग करना बहुत सीमित होगा।

अच्छी खबर: हमारे datasheet में और इसलिए हमारे meta map में एक `character` कॉलम है।
चलो उसका उपयोग करते हैं ताकि प्रत्येक entry के लिए प्रोसेस को किस character का उपयोग करना चाहिए यह सेट किया जा सके।

इसके लिए, हमें तीन काम करने होंगे:

1. पिछले प्रोसेस से आने वाले output channel को एक नाम दो ताकि हम उस पर अधिक सुविधाजनक तरीके से काम कर सकें।
2. रुचि की जानकारी तक कैसे पहुँचें यह निर्धारित करो।
3. दूसरे प्रोसेस को call करो और जानकारी उचित तरीके से feed करो।

चलो शुरू करते हैं।

#### 3.2.1. पिछले output channel को नाम देना

हमने पिछले manipulations को सीधे पहले प्रोसेस `IDENTIFY_LANGUAGE.out` के output channel पर लागू किया।
उस channel की सामग्री को अगले प्रोसेस में feed करने के लिए (और इसे स्पष्ट और पढ़ने में आसान तरीके से करने के लिए) हम इसे अपना नाम `ch_languages` देना चाहते हैं।

हम [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) ऑपरेटर का उपयोग करके ऐसा कर सकते हैं।

main workflow में, `.view()` ऑपरेटर को `.set { ch_languages }` से बदलो, और एक लाइन जोड़ो जो test करे कि हम channel को नाम से refer कर सकते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएं
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .set { ch_languages }

        // अस्थायी: ch_languages में झाँकें
        ch_languages.view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएं
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

चलो यह चलाते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [friendly_austin] DSL2 - revision: 3dbe460fd6

    [36/cca6a7] IDENTIFY_LANGUAGE (7) | 7 of 7 ✔
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/e2/6db2402d83cf72081bcd2d11784714/guten_tag.txt]
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/6c/114c818317d169457d6e7336d5d55b/bonjour.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/55/68c69c5efb527f3604ddb3daab8057/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/work/2a/4752055ccb5d1370b0ef9da41d3993/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/work/f4/fcd3186dc666d5d239ffa6c37d125d/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/c3/3b2627f733f278a7088332a5806108/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/work/36/cca6a7dbfa26ac24f9329787a32e9d/ciao.txt]
    ```

यह confirm करता है कि हम अब channel को नाम से refer कर सकते हैं।

#### 3.2.2. फ़ाइल और character मेटाडेटा तक पहुँचना

हम module कोड देखकर जानते हैं कि `COWPY` प्रोसेस को एक text file और एक `character` value दिए जाने की उम्मीद है।
`COWPY` प्रोसेस call लिखने के लिए, हमें बस यह जानना है कि channel में प्रत्येक element से संबंधित फ़ाइल object और मेटाडेटा कैसे निकालें।

जैसा कि अक्सर होता है, ऐसा करने का सबसे सरल तरीका `map` ऑपरेशन का उपयोग करना है।

हमारे channel में `[meta, file]` के रूप में structured tuples हैं, इसलिए हम `file` object तक सीधे पहुँच सकते हैं, और meta map के अंदर stored `character` value तक `meta.character` के रूप में refer करके पहुँच सकते हैं।

main workflow में, निम्नलिखित code बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="34"
        // अस्थायी: फ़ाइल और character तक पहुँचें
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="34"
        // अस्थायी: ch_languages में झाँकें
        ch_languages.view()
    ```

ध्यान दो कि हम closures (जैसे `{ file -> "File: " + file }`) का उपयोग `.view` ऑपरेशन के आउटपुट को अधिक पठनीय बनाने के लिए कर रहे हैं।

चलो यह चलाते हैं:

```bash
nextflow run main.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Character: squirrel
    File: /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    File: /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Character: tux
    Character: turkey
    File: /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    File: /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Character: sheep
    Character: moose
    Character: stegosaurus
    File: /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    File: /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    File: /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Character: turtle
    ```

_फ़ाइल पाथ और character values तुम्हारे आउटपुट में अलग क्रम में आ सकते हैं।_

यह confirm करता है कि हम channel में प्रत्येक element के लिए फ़ाइल और character तक पहुँचने में सक्षम हैं।

#### 3.2.3. `COWPY` प्रोसेस को call करना

अब चलो सब कुछ एक साथ रखते हैं और वास्तव में `ch_languages` channel पर `COWPY` प्रोसेस को call करते हैं।

main workflow में, निम्नलिखित code बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="34"
        // ASCII art generate करने के लिए cowpy चलाएं
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="34"
        // अस्थायी: फ़ाइल और character तक पहुँचें
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

तुम देख सकते हो कि हम बस दो map ऑपरेशन (`.view()` statements के बिना) को प्रोसेस call के inputs के रूप में copy करते हैं।
बस यह सुनिश्चित करो कि उनके बीच comma न भूलो!

यह थोड़ा clunky है, लेकिन हम अगले section में देखेंगे कि इसे बेहतर कैसे बनाया जाए।

चलो यह चलाते हैं:

```bash
nextflow run main.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (7)
    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [e7/317c18] COWPY (6)             [100%] 7 of 7 ✔
    ```

अगर तुम results डायरेक्टरी में देखो, तो तुम्हें संबंधित character द्वारा बोले गए प्रत्येक अभिवादन की ASCII art वाली individual files दिखनी चाहिए।

??? abstract "डायरेक्टरी और उदाहरण फ़ाइल सामग्री"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

    ```text title="results/cowpy-bonjour.txt"
     _________________
    / Bonjour         \
    \ Salut, à demain /
    -----------------
      \
        \
                      _ _
          | \__/|  .~    ~.
          /oo `./      .'
          {o__,   \    {
            / .  . )    \
            `-` '-' \    }
          .(   _(   )_.'
          '---.~_ _ _|
    ```

इससे पता चलता है कि हम पाइपलाइन के दूसरे चरण में command को parameterize करने के लिए meta map में जानकारी का उपयोग करने में सक्षम थे।

हालाँकि, जैसा कि ऊपर उल्लेख किया गया है, इसमें शामिल कुछ कोड थोड़ा clunky था, क्योंकि हमें workflow body के context में ही meta data unpack करना पड़ा।
यह approach meta map से कम संख्या में fields का उपयोग करने के लिए ठीक काम करती है, लेकिन अगर हम बहुत अधिक fields का उपयोग करना चाहते तो यह poorly scale होती।

`multiMap()` नामक एक और ऑपरेटर है जो इसे थोड़ा streamline करने की अनुमति देता है, लेकिन तब भी यह ideal नहीं है।

??? info "(वैकल्पिक) `multiMap()` के साथ वैकल्पिक version"

    अगर तुम सोच रहे हो, तो हम एक single `map()` ऑपरेशन नहीं लिख सकते थे जो `file` और `character` दोनों output करे, क्योंकि वह उन्हें tuple के रूप में return करता।
    हमें `file` और `character` elements को प्रोसेस में अलग-अलग feed करने के लिए दो अलग `map()` ऑपरेशन लिखने पड़े।

    तकनीकी रूप से `multiMap()` ऑपरेटर का उपयोग करके एक single mapping ऑपरेशन के माध्यम से ऐसा करने का एक और तरीका है, जो multiple channels emit करने में सक्षम है।
    उदाहरण के लिए, तुम ऊपर `COWPY` के call को निम्नलिखित कोड से बदल सकते हो:

    === "बाद में"

        ```groovy title="main.nf" linenums="34"
            // ASCII art generate करने के लिए cowpy चलाएं
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "पहले"

        ```groovy title="main.nf" linenums="34"
            // ASCII art generate करने के लिए cowpy चलाएं
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    यह बिल्कुल वही परिणाम देता है।

किसी भी मामले में, यह awkward है कि हमें workflow level पर कुछ unpacking करनी पड़ती है।

बेहतर होगा अगर हम पूरा meta map प्रोसेस में feed कर सकें और वहाँ जो चाहिए वह pick कर सकें।

### 3.3. पूरा meta map pass करना और उपयोग करना

meta map का उद्देश्य आखिरकार सभी मेटाडेटा को एक bundle के रूप में pass करना है।
एकमात्र कारण जिससे हम ऊपर ऐसा नहीं कर सके वह यह है कि प्रोसेस meta map accept करने के लिए set up नहीं है।
लेकिन चूँकि हम प्रोसेस कोड को control करते हैं, हम इसे बदल सकते हैं।

चलो `COWPY` प्रोसेस को `[meta, file]` tuple संरचना accept करने के लिए modify करते हैं जो हमने पहले प्रोसेस में उपयोग की थी ताकि हम वर्कफ़्लो को streamline कर सकें।

इसके लिए, हमें तीन काम करने होंगे:

1. `COWPY` process module की input definitions modify करो
2. meta map का उपयोग करने के लिए process command update करो
3. workflow body में process call update करो

तैयार हो? चलो शुरू करते हैं!

#### 3.3.1. `COWPY` module input modify करना

`cowpy.nf` module फ़ाइल में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2"
    input:
    tuple val(meta), path(input_file)
    ```

=== "पहले"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2-3"
    input:
    path(input_file)
    val character
    ```

यह हमें ट्यूटोरियल में पहले cover की गई `[meta, file]` tuple संरचना का उपयोग करने में सक्षम बनाता है।

ध्यान दो कि हमने ट्यूटोरियल को streamlined रखने के लिए meta map output करने के लिए process output definition update नहीं की, लेकिन `IDENTIFY_LANGUAGE` प्रोसेस के model का अनुसरण करते हुए इसे खुद exercise के रूप में करने के लिए स्वतंत्र महसूस करो।

#### 3.3.2. Meta map field का उपयोग करने के लिए command update करना

पूरा meta map अब प्रोसेस के अंदर उपलब्ध है, इसलिए हम command block के अंदर से सीधे इसमें मौजूद जानकारी को refer कर सकते हैं।

`cowpy.nf` module फ़ाइल में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
    """
    ```

=== "पहले"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
    ```

हमने standalone input के रूप में पहले pass किए गए `character` value के reference को meta map में held value से बदल दिया है, जिसे हम `meta.character` का उपयोग करके refer करते हैं।

अब process call को उसी के अनुसार update करते हैं।

#### 3.3.3. Process call update करना और चलाना

प्रोसेस अब अपने input के लिए `[meta, file]` tuple संरचना की उम्मीद करता है, जो पिछला प्रोसेस output करता है, इसलिए हम बस पूरे `ch_languages` channel को `COWPY` प्रोसेस में pass कर सकते हैं।

main workflow में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // ASCII art generate करने के लिए cowpy चलाएं
    COWPY(ch_languages)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // ASCII art generate करने के लिए cowpy चलाएं
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

इससे call काफी सरल हो जाती है!

चलो पिछले execution के results delete करते हैं और इसे चलाते हैं:

```bash
rm -r results
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (14)
    [5d/dffd4e] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [25/9243df] process > COWPY (7)             [100%] 7 of 7 ✔
    ```

अगर तुम results डायरेक्टरी में देखो, तो तुम्हें पहले जैसे ही आउटपुट दिखने चाहिए, _यानी_ संबंधित character द्वारा बोले गए प्रत्येक अभिवादन की ASCII art वाली individual files।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    ./results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

तो यह सरल कोड के साथ पहले जैसे ही results देता है।

बेशक, यह मानता है कि तुम process कोड modify करने में सक्षम हो।
कुछ मामलों में, तुम्हें मौजूदा processes पर निर्भर रहना पड़ सकता है जिन्हें तुम modify करने की स्थिति में नहीं हो, जो तुम्हारे विकल्पों को सीमित करता है।
अच्छी खबर, अगर तुम [nf-core](https://nf-co.re/) project के modules का उपयोग करने की योजना बना रहे हो, तो nf-core modules सभी standard के रूप में `[meta, file]` tuple संरचना का उपयोग करने के लिए set up हैं।

### 3.4. Missing required inputs की troubleshooting

`COWPY` प्रोसेस को सफलतापूर्वक चलाने के लिए `character` value आवश्यक है।
अगर हम configuration file में इसके लिए default value सेट नहीं करते, तो हमें datasheet में इसके लिए एक value ज़रूर प्रदान करनी होगी।

**क्या होगा अगर हम नहीं करते?**
यह इस बात पर निर्भर करता है कि input datasheet में क्या है और हम वर्कफ़्लो का कौन सा version चला रहे हैं।

#### 3.4.1. Character column मौजूद है लेकिन खाली है

मान लो हम data collection error simulate करने के लिए अपनी datasheet में एक entry के लिए character value delete करते हैं:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

ऊपर उपयोग किए गए वर्कफ़्लो के किसी भी version के लिए, datasheet पढ़े जाने पर सभी entries के लिए `character` key बनाई जाएगी, लेकिन `sampleA` के लिए value एक empty string होगी।

इससे एक error आएगी।

??? failure "कमांड आउटपुट"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > cowpy-bonjour.txt

    Command exit status:
      2

    Command output:
      (empty)

    Command error:
      usage: cowpy [-h] [-l] [-L] [-t] [-u] [-e EYES] [-c COWACTER] [-E] [-r] [-x]
                  [-C]
                  [msg ...]
      cowpy: error: argument -c/--cowacter: expected one argument

    Work dir:
      /workspaces/training/side-quests/metadata/work/ca/9d49796612a54dec5ed466063c809b

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

जब Nextflow उस नमूने के लिए `cowpy` command line चलाता है, तो `${meta.character}` `cowpy` command line में एक empty string से भर जाता है, इसलिए `cowpy` टूल एक error throw करता है जो कहता है कि `-c` argument के लिए कोई value प्रदान नहीं की गई।

#### 3.4.2. Datasheet में character column मौजूद नहीं है

अब मान लो हम अपनी datasheet से `character` column पूरी तरह delete करते हैं:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,/workspaces/training/side-quests/metadata/data/ciao.txt
```

इस मामले में datasheet पढ़े जाने पर `character` key बिल्कुल नहीं बनाई जाएगी।

##### 3.4.2.1. Workflow level पर access की गई value

अगर हम section 3.2 में लिखे गए कोड के version का उपयोग कर रहे हैं, तो Nextflow `COWPY` प्रोसेस को call करने से पहले meta map में `character` key access करने का प्रयास करेगा।

उसे कोई ऐसे elements नहीं मिलेंगे जो instruction से match करें, इसलिए वह `COWPY` बिल्कुल नहीं चलाएगा।

??? success "कमांड आउटपुट"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Nextflow की दृष्टि से, यह वर्कफ़्लो सफलतापूर्वक चला!
हालाँकि, हम जो आउटपुट चाहते हैं उनमें से कोई भी produce नहीं होगा।

##### 3.4.2.2. Process level पर access की गई value

अगर हम section 3.3 के version का उपयोग कर रहे हैं, तो Nextflow पूरा meta map `COWPY` प्रोसेस में pass करेगा और command चलाने का प्रयास करेगा।

इससे एक error आएगी, लेकिन पहले मामले की तुलना में अलग।

??? failure "कमांड आउटपुट"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > cowpy-guten_tag.txt

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Traceback (most recent call last):
        File "/opt/conda/bin/cowpy", line 10, in <module>
          sys.exit(main())
                  ~~~~^^
        File "/opt/conda/lib/python3.13/site-packages/cowpy/cow.py", line 1215, in main
          print(cow(eyes=args.eyes,
                ~~~^^^^^^^^^^^^^^^^
                tongue=args.tongue,
                ^^^^^^^^^^^^^^^^^^^
                thoughts=args.thoughts
                ^^^^^^^^^^^^^^^^^^^^^^
                    ).milk(msg)
                    ^
      TypeError: 'str' object is not callable

    Work dir:
      /workspaces/training/side-quests/metadata/work/06/28065f7d9fd7d22bba084aa941b6d6

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

यह इसलिए होता है क्योंकि `meta.character` मौजूद नहीं है, इसलिए इसे access करने का हमारा प्रयास `null` return करता है। परिणामस्वरूप, Nextflow literally command-line में `null` plug in करता है, जिसे `cowpy` टूल naturally नहीं पहचानता।

#### 3.4.3. समाधान

Workflow configuration के हिस्से के रूप में default value supply करने के अलावा, इसे अधिक robustly handle करने के लिए हम दो काम कर सकते हैं:

1. अपने वर्कफ़्लो में input validation implement करो ताकि यह सुनिश्चित हो सके कि datasheet में सभी आवश्यक जानकारी है। तुम Hello nf-core training course में [input validation का परिचय](../hello_nf-core/05_input_validation.md) पा सकते हो। <!-- TODO (future) pending a proper Validation side quest -->

2. अगर तुम यह सुनिश्चित करना चाहते हो कि तुम्हारे process module का उपयोग करने वाला कोई भी व्यक्ति तुरंत required inputs की पहचान कर सके, तो तुम required metadata property को एक explicit input भी बना सकते हो।

यहाँ एक उदाहरण है कि यह कैसे काम करेगा।

पहले, process level पर, input definition को इस प्रकार update करो:

=== "बाद में"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "पहले"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

फिर, workflow level पर, `character` property को metadata से extract करने और इसे input tuple का explicit component बनाने के लिए एक mapping ऑपरेशन का उपयोग करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

इस approach का यह फायदा है कि यह स्पष्ट रूप से दिखाता है कि `character` आवश्यक है, और प्रोसेस को अन्य contexts में redeploy करना आसान बनाता है।

यह एक महत्वपूर्ण design principle को उजागर करता है:

**Meta map का उपयोग optional, descriptive जानकारी के लिए करो, लेकिन required values को explicit inputs के रूप में extract करो।**

Meta map channel structures को clean रखने और arbitrary channel structures को रोकने के लिए उत्कृष्ट है, लेकिन mandatory elements के लिए जो किसी प्रोसेस में directly referenced हैं, उन्हें explicit inputs के रूप में extract करने से अधिक robust और maintainable कोड बनता है।

### सारांश

इस section में, तुमने सीखा कि किसी प्रोसेस के execution को customize करने के लिए मेटाडेटा का उपयोग कैसे करें, इसे workflow level पर या process level पर access करते हुए।

---

## अनुपूरक अभ्यास

अगर तुम किसी प्रोसेस के अंदर से meta map जानकारी का उपयोग करने का अभ्यास करना चाहते हो, तो meta map से अन्य जानकारी जैसे `lang` और `lang_group` का उपयोग करके आउटपुट के नाम और/या व्यवस्था को customize करने की कोशिश करो।

उदाहरण के लिए, इस परिणाम को produce करने के लिए कोड modify करने की कोशिश करो:

```console title="Results directory contents"
results/
├── germanic
│   ├── de-guten_tag.txt
│   ├── de-hallo.txt
│   └── en-hello.txt
└── romance
    ├── es-hola.txt
    ├── fr-bonjour.txt
    ├── fr-salut.txt
    └── it-ciao.txt
```

<!-- TODO (future) Provide worked out solution -->
<!-- the renaming should use the meta inside the process -->
<!-- the output org should use the meta in the workflow outputs -->

---

## सारांश

इस side quest में, तुमने Nextflow वर्कफ़्लो में मेटाडेटा के साथ प्रभावी ढंग से काम करने का तरीका explore किया।

मेटाडेटा को explicit रखने और डेटा के साथ जोड़े रखने का यह पैटर्न Nextflow में एक core best practice है, जो file information को hard-code करने की तुलना में कई फायदे प्रदान करता है:

- File metadata पूरे वर्कफ़्लो में files के साथ जुड़ी रहती है
- प्रत्येक फ़ाइल के लिए process behavior customize किया जा सकता है
- Output organization file metadata को reflect कर सकती है
- Pipeline execution के दौरान file information expand की जा सकती है

अपने काम में इस पैटर्न को लागू करने से तुम robust, maintainable bioinformatics वर्कफ़्लो बनाने में सक्षम होगे।

### मुख्य पैटर्न

1.  **मेटाडेटा पढ़ना और संरचित करना:** CSV फ़ाइलें पढ़ना और organized metadata maps बनाना जो तुम्हारी data files के साथ जुड़े रहते हैं।

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **वर्कफ़्लो के दौरान मेटाडेटा expand करना:** Process outputs जोड़कर और conditional logic के माध्यम से values derive करके अपनी पाइपलाइन के आगे बढ़ने पर अपने मेटाडेटा में नई जानकारी जोड़ना।

    - Process output के आधार पर नई keys जोड़ना

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - Conditional clause का उपयोग करके नई keys जोड़ना

    ```groovy
    .map{ meta, file ->
        if ( meta.lang.equals("de") || meta.lang.equals('en') ){
            lang_group = "germanic"
        } else if ( meta.lang in ["fr", "es", "it"] ) {
            lang_group = "romance"
        } else {
            lang_group = "unknown"
        }
    }
    ```

3.  **Process Behavior Customize करना:** प्रोसेस के अंदर मेटाडेटा का उपयोग करना।

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### अतिरिक्त संसाधन

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## आगे क्या है?

[Side Quests के मेनू](../) पर वापस जाओ या सूची में अगले विषय पर जाने के लिए पृष्ठ के नीचे दाईं ओर बटन पर क्लिक करो।
