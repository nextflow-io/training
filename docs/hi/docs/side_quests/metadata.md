# मेटाडेटा और मेटा मैप्स

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

किसी भी वैज्ञानिक विश्लेषण में, हम शायद ही कभी केवल कच्ची डेटा फ़ाइलों के साथ काम करते हैं।
प्रत्येक फ़ाइल अपनी अतिरिक्त जानकारी के साथ आती है: यह क्या है, यह कहाँ से आई, और इसे क्या खास बनाता है।
इस अतिरिक्त जानकारी को हम मेटाडेटा कहते हैं।

मेटाडेटा अन्य डेटा का वर्णन करने वाला डेटा है।
मेटाडेटा फ़ाइलों और प्रयोगात्मक स्थितियों के बारे में महत्वपूर्ण विवरण ट्रैक करता है, और प्रत्येक डेटासेट की अनूठी विशेषताओं के अनुसार विश्लेषण को तैयार करने में मदद करता है।

इसे एक लाइब्रेरी कैटलॉग की तरह समझें: जबकि पुस्तकों में वास्तविक सामग्री (कच्चा डेटा) होती है, कैटलॉग कार्ड प्रत्येक पुस्तक के बारे में आवश्यक जानकारी प्रदान करते हैं—यह कब प्रकाशित हुई थी, इसे किसने लिखा, इसे कहाँ खोजें (मेटाडेटा)।
Nextflow pipelines में, मेटाडेटा का उपयोग निम्न के लिए किया जा सकता है:

- पूरे workflow में फ़ाइल-विशिष्ट जानकारी को ट्रैक करना
- फ़ाइल विशेषताओं के आधार पर processes को कॉन्फ़िगर करना
- संयुक्त विश्लेषण के लिए संबंधित फ़ाइलों को समूहित करना

### सीखने के लक्ष्य

इस side quest में, हम workflows में मेटाडेटा को संभालने का तरीका जानेंगे।
बुनियादी फ़ाइल जानकारी वाली एक साधारण datasheet (जिसे अक्सर bioinformatics में samplesheet कहा जाता है) से शुरू करते हुए, आप सीखेंगे कि कैसे:

- CSV फ़ाइलों से फ़ाइल मेटाडेटा पढ़ें और पार्स करें
- मेटाडेटा मैप्स बनाएं और उनमें परिवर्तन करें
- workflow निष्पादन के दौरान नए मेटाडेटा फ़ील्ड जोड़ें
- process व्यवहार को अनुकूलित करने के लिए मेटाडेटा का उपयोग करें

ये कौशल आपको अधिक मजबूत और लचीली pipelines बनाने में मदद करेंगे जो जटिल फ़ाइल संबंधों और प्रोसेसिंग आवश्यकताओं को संभाल सकते हैं।

### पूर्वापेक्षाएँ

इस side quest को शुरू करने से पहले, आपको:

- [Hello Nextflow](../hello_nextflow/README.md) tutorial या समकक्ष beginner's course पूरा किया होना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (processes, channels, operators) का उपयोग करने में सहज होना चाहिए

---

## 0. शुरुआत करें

#### प्रशिक्षण codespace खोलें

यदि आपने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित प्रशिक्षण वातावरण को खोलना सुनिश्चित करें।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाएं

आइए उस डायरेक्टरी में चलें जहाँ इस tutorial के लिए फ़ाइलें स्थित हैं।

```bash
cd side-quests/metadata
```

आप VSCode को इस डायरेक्टरी पर फ़ोकस करने के लिए सेट कर सकते हैं:

```bash
code .
```

#### सामग्री की समीक्षा करें

आपको एक मुख्य workflow फ़ाइल और एक `data` डायरेक्टरी मिलेगी जिसमें एक datasheet और कुछ डेटा फ़ाइलें हैं।

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

`main.nf` फ़ाइल में workflow एक stub है जिसे आप धीरे-धीरे एक पूर्ण रूप से कार्यशील workflow में विस्तारित करेंगे।

datasheet डेटा फ़ाइलों के paths और कुछ संबंधित मेटाडेटा को सूचीबद्ध करती है, जो 3 columns में व्यवस्थित हैं:

- `id`: स्व-व्याख्यात्मक, फ़ाइल को दी गई एक ID
- `character`: एक character नाम, जिसका उपयोग हम बाद में विभिन्न creatures बनाने के लिए करेंगे
- `data`: `.txt` फ़ाइलों के paths जिनमें विभिन्न भाषाओं में अभिवादन शामिल हैं

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

प्रत्येक डेटा फ़ाइल में पाँच भाषाओं (fr: French, de: German, es: Spanish, it: Italian, en: English) में से एक में कुछ अभिवादन text शामिल है।

हम आपको `langid` नामक एक containerized भाषा विश्लेषण tool भी प्रदान करेंगे।

#### असाइनमेंट की समीक्षा करें

आपकी चुनौती एक Nextflow workflow लिखना है जो:

1. प्रत्येक फ़ाइल में भाषा को स्वचालित रूप से **पहचानेगा**
2. फ़ाइलों को भाषा परिवार (Germanic बनाम Romance languages) द्वारा **समूहित** करेगा
3. प्रत्येक फ़ाइल के लिए उसकी भाषा और मेटाडेटा के आधार पर processing को **अनुकूलित** करेगा
4. outputs को भाषा समूह द्वारा **व्यवस्थित** करेगा

यह एक विशिष्ट workflow pattern का प्रतिनिधित्व करता है जहाँ फ़ाइल-विशिष्ट मेटाडेटा processing निर्णयों को संचालित करता है; बिल्कुल वैसी ही समस्या जिसे मेटाडेटा मैप्स सुंदरता से हल करते हैं।

#### तैयारी checklist

क्या आपको लगता है कि आप शुरू करने के लिए तैयार हैं?

- [ ] मैं इस course के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता हूँ
- [ ] मेरा codespace चालू है और चल रहा है
- [ ] मैंने अपनी working डायरेक्टरी को उचित रूप से सेट कर लिया है
- [ ] मैं असाइनमेंट को समझता हूँ

यदि आप सभी boxes को चेक कर सकते हैं, तो आप जाने के लिए तैयार हैं।

---

## 1. datasheet से मेटाडेटा लोड करें

`main.nf` workflow फ़ाइल खोलें और workflow stub की जाँच करें जो हम आपको शुरुआती बिंदु के रूप में दे रहे हैं।

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

आप देख सकते हैं कि हमने उदाहरण datasheet को फ़ाइल के रूप में लोड करने के लिए एक बुनियादी channel factory सेट अप किया है, लेकिन यह अभी तक फ़ाइल की सामग्री नहीं पढ़ेगा।
आइए इसे जोड़कर शुरू करें।

### 1.1. `splitCsv` के साथ सामग्री पढ़ें

हमें एक ऑपरेटर चुनने की आवश्यकता है जो हमारी ओर से न्यूनतम प्रयास के साथ फ़ाइल सामग्री को उचित रूप से पार्स करे।
चूंकि हमारी datasheet CSV format में है, यह [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) ऑपरेटर के लिए एक काम है, जो फ़ाइल में प्रत्येक पंक्ति को channel में एक तत्व के रूप में लोड करता है।

channel निर्माण कोड में `splitCsv()` operation जोड़ने के लिए निम्नलिखित परिवर्तन करें, साथ ही यह जांचने के लिए `view()` operation भी कि फ़ाइल की सामग्री channel में सही ढंग से लोड हो रही है।

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

ध्यान दें कि हम Nextflow को CSV फ़ाइल की पहली पंक्ति को header पंक्ति के रूप में पढ़ने के लिए `header: true` option का उपयोग कर रहे हैं।

आइए देखें कि इससे क्या निकलता है, ठीक है?
workflow चलाएँ:

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

हम देख सकते हैं कि ऑपरेटर ने CSV फ़ाइल में प्रत्येक पंक्ति के लिए key-value pairs का एक map बनाया है, जिसमें संबंधित values के लिए keys के रूप में column headers हैं।

प्रत्येक map entry हमारी datasheet में एक column से मेल खाती है:

- `id`
- `character`
- `recording`

यह बहुत बढ़िया है! यह प्रत्येक फ़ाइल से विशिष्ट फ़ील्ड्स तक पहुँचना आसान बनाता है।
उदाहरण के लिए, हम `id` के साथ फ़ाइल ID या `recording` के साथ txt फ़ाइल path तक पहुँच सकते हैं।

??? info "(वैकल्पिक) मैप्स के बारे में अधिक"

    Groovy में, प्रोग्रामिंग भाषा जिस पर Nextflow बना है, एक map एक key-value डेटा structure है जो Python में dictionaries, JavaScript में objects, या Ruby में hashes के समान है।

    यहाँ एक runnable script है जो दिखाता है कि आप व्यवहार में एक map कैसे परिभाषित कर सकते हैं और इसकी सामग्री तक कैसे पहुँच सकते हैं:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // एक साधारण map बनाएं
    def my_map = [id:'sampleA', character:'squirrel']

    // पूरा map प्रिंट करें
    println "map: ${my_map}"

    // dot notation का उपयोग करके व्यक्तिगत values तक पहुँचें
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    भले ही इसमें एक उचित `workflow` block नहीं है, Nextflow इसे workflow के रूप में चला सकता है:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    और यहाँ आउटपुट में आप क्या देखने की उम्मीद कर सकते हैं:

    ```console title="आउटपुट"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. `map` के साथ विशिष्ट फ़ील्ड्स चुनें

मान लीजिए हम datasheet से `character` column तक पहुँचना और इसे प्रिंट करना चाहते हैं।
हम अपने channel में प्रत्येक item पर iterate करने और विशेष रूप से map object से `character` entry चुनने के लिए Nextflow `map` ऑपरेटर का उपयोग कर सकते हैं।

workflow में निम्नलिखित संपादन करें:

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

अब फिर से workflow चलाएँ:

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

सफलता! हमने प्रत्येक पंक्ति के लिए व्यक्तिगत columns से values तक पहुँचने के लिए अपनी datasheet से प्राप्त map structure का लाभ उठाया है।

अब जब हमने सफलतापूर्वक datasheet पढ़ ली है और प्रत्येक पंक्ति में डेटा तक पहुँच है, तो हम अपनी pipeline logic को लागू करना शुरू कर सकते हैं।

### 1.3. मेटाडेटा को 'meta map' में व्यवस्थित करें

workflow की वर्तमान स्थिति में, इनपुट फ़ाइलें (`recording` key के तहत) और संबंधित मेटाडेटा (`id`, `character`) सभी एक ही स्तर पर हैं, जैसे वे सभी एक बड़े bag में हैं।
व्यावहारिक परिणाम यह है कि इस channel का उपभोग करने वाले प्रत्येक process को इस structure को ध्यान में रखते हुए configure करने की आवश्यकता होगी:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

यह ठीक है जब तक datasheet में columns की संख्या नहीं बदलती।
हालांकि, यदि आप datasheet में केवल एक column जोड़ते हैं, तो channel का आकार उस से मेल नहीं खाएगा जो process उम्मीद करता है, और workflow त्रुटियाँ उत्पन्न करेगा।
यह process को दूसरों के साथ साझा करना भी मुश्किल बनाता है जिनके पास थोड़ा अलग इनपुट डेटा हो सकता है, और आपको variables को process में hard-code करना पड़ सकता है जो script block द्वारा आवश्यक नहीं हैं।

इस समस्या से बचने के लिए, हमें channel structure को consistent रखने का एक तरीका खोजने की आवश्यकता है, भले ही उस datasheet में कितने columns हों।

हम सभी मेटाडेटा को tuple के भीतर एक item में एकत्र करके ऐसा कर सकते हैं, जिसे हम मेटाडेटा map, या अधिक सरलता से 'meta map' कहेंगे।

`map` operation में निम्नलिखित संपादन करें:

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

हमने अपने channel elements को दो तत्वों, meta map और संबंधित फ़ाइल object से मिलकर बने tuple में पुनर्गठित किया है।

आइए workflow चलाएँ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console title="meta map देखें"
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

अब, channel में प्रत्येक तत्व में मेटाडेटा map पहले और संबंधित फ़ाइल object दूसरे स्थान पर है:

```console title="उदाहरण आउटपुट structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

परिणामस्वरूप, datasheet में अधिक columns जोड़ने से `meta` map में अधिक मेटाडेटा उपलब्ध होगा, लेकिन channel का आकार नहीं बदलेगा।
यह हमें ऐसी processes लिखने में सक्षम बनाता है जो channel का उपभोग करते हैं बिना मेटाडेटा items को इनपुट specification में hard-code किए:

```groovy title="Syntax उदाहरण"
    input:
    tuple val(meta), file(recording)
```

यह Nextflow workflows में मेटाडेटा को व्यवस्थित करने के लिए व्यापक रूप से उपयोग किया जाने वाला convention है।

### मुख्य बातें

इस खंड में, आपने सीखा:

- **मेटाडेटा क्यों महत्वपूर्ण है:** अपने डेटा के साथ मेटाडेटा रखने से पूरे workflow में महत्वपूर्ण फ़ाइल जानकारी संरक्षित होती है।
- **datasheets को कैसे पढ़ें:** header जानकारी के साथ CSV फ़ाइलों को पढ़ने और पंक्तियों को structured डेटा में transform करने के लिए `splitCsv` का उपयोग करना
- **meta map कैसे बनाएं:** tuple structure `[ [id:value, ...], file ]` का उपयोग करके फ़ाइल डेटा से मेटाडेटा को अलग करना

---

## 2. मेटाडेटा में परिवर्तन करना

अब जब हमारा मेटाडेटा लोड हो गया है, तो आइए इसके साथ कुछ करें!

हम प्रत्येक creature की recording फ़ाइल में निहित भाषा की पहचान करने के लिए [`langid`](https://github.com/saffsd/langid.py) नामक एक tool का उपयोग करने जा रहे हैं।
यह tool भाषाओं के एक सेट पर पूर्व-प्रशिक्षित आता है, और text का एक snippet दिए जाने पर, यह एक भाषा prediction और एक संबंधित probability score दोनों को `stdout` में आउटपुट करेगा।

### 2.1. process को import करें और कोड की जाँच करें

हम आपको `IDENTIFY_LANGUAGE` नामक एक पूर्व-लिखित process मॉड्यूल प्रदान करते हैं जो `langid` tool को wrap करता है, इसलिए आपको workflow block से पहले केवल एक include statement जोड़ने की आवश्यकता है।

workflow में निम्नलिखित संपादन करें:

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

आप इसके कोड की जाँच करने के लिए मॉड्यूल फ़ाइल खोल सकते हैं:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// प्रत्येक इनपुट फ़ाइल की भाषा की भविष्यवाणी करने के लिए langid का उपयोग करें
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

जैसा कि आप देख सकते हैं, इनपुट definition उसी `tuple val(meta), path(file)` structure का उपयोग करता है जिसे हमने अभी अपने इनपुट channel पर लागू किया था।

आउटपुट definition एक tuple के रूप में structured है जिसमें इनपुट के समान structure है, सिवाय इसके कि इसमें तीसरे तत्व के रूप में `stdout` भी है।
यह `tuple val(meta), path(file), <output>` pattern मेटाडेटा को इनपुट डेटा और outputs दोनों से जुड़ा रखता है क्योंकि यह pipeline के माध्यम से प्रवाहित होता है।

ध्यान दें कि हम यहाँ Nextflow के [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) आउटपुट qualifier का उपयोग कर रहे हैं क्योंकि tool अपने आउटपुट को फ़ाइल लिखने के बजाय सीधे console पर print करता है; और हम probability score को हटाने, newline characters को हटाकर string को साफ करने, और केवल भाषा prediction वापस करने के लिए कमांड लाइन में `sed` का उपयोग करते हैं।

### 2.2. `IDENTIFY_LANGUAGE` की call जोड़ें

अब जब process workflow के लिए उपलब्ध है, तो हम डेटा channel पर इसे चलाने के लिए `IDENTIFY_LANGUAGE` process में एक call जोड़ सकते हैं।

workflow में निम्नलिखित संपादन करें:

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएँ
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

ध्यान दें कि हमने channel निर्माण में मूल `.view()` operation हटा दिया है।

अब हम workflow चला सकते हैं:

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

उत्कृष्ट! अब हमारे पास इस बात की prediction है कि प्रत्येक character कौन सी भाषा बोलता है।

और जैसा कि पहले उल्लेख किया गया है, हमने आउटपुट में इनपुट फ़ाइल और meta map भी शामिल किया है, जिसका अर्थ है कि दोनों नई जानकारी के साथ जुड़े रहते हैं जो हमने अभी उत्पन्न की है।
यह अगले चरण में उपयोगी साबित होगा।

!!! note

    अधिक सामान्यतः, meta map को results से जुड़ा रखने का यह pattern संबंधित results को जोड़ना आसान बनाता है जो समान identifiers साझा करते हैं।

    जैसा कि आप पहले ही सीख चुके होंगे, आप results को उनके बीच match करने के लिए channels में items के क्रम पर भरोसा नहीं कर सकते।
    इसके बजाय, आपको डेटा को सही ढंग से जोड़ने के लिए keys का उपयोग करना चाहिए, और meta maps इस उद्देश्य के लिए एक आदर्श structure प्रदान करते हैं।

    हम [Splitting & Grouping](./splitting_and_grouping.md) side quest में इस use case का विस्तार से पता लगाते हैं।

### 2.3. process outputs के साथ मेटाडेटा बढ़ाएं

यह देखते हुए कि हमने अभी जो परिणाम उत्पन्न किए हैं वे स्वयं फ़ाइलों की सामग्री के बारे में मेटाडेटा का एक रूप हैं, उन्हें हमारे meta map में जोड़ना उपयोगी होगा।

हालाँकि, हम मौजूदा meta map को स्थान पर संशोधित नहीं करना चाहते।
तकनीकी दृष्टिकोण से, ऐसा करना _संभव_ है, लेकिन यह असुरक्षित है।

इसलिए इसके बजाय, हम एक नया meta map बनाएंगे जिसमें मौजूदा meta map की सामग्री और एक नई `lang: lang_id` key-value pair होगी जो नई जानकारी रखेगी, `+` ऑपरेटर (एक Groovy feature) का उपयोग करते हुए।
और हम इसे एक [`map`](https://www.nextflow.io/docs/latest/operator.html#map) operation के साथ जोड़ेंगे ताकि पुराने map को नए से बदल सकें।

यहाँ वे संपादन हैं जो आपको workflow में करने की आवश्यकता है:

=== "बाद में"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएँ
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएँ
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

यदि आप अभी तक `+` ऑपरेटर से परिचित नहीं हैं, या यदि यह भ्रमित लगता है, तो नीचे दिए गए विस्तृत स्पष्टीकरण को पढ़ने में कुछ मिनट लगाएं।

??? info "`+` ऑपरेटर का उपयोग करके नए meta map का निर्माण"

    **सबसे पहले, आपको जानना होगा कि हम Groovy ऑपरेटर `+` का उपयोग करके दो maps की सामग्री को merge कर सकते हैं।**

    मान लीजिए हमारे पास निम्नलिखित maps हैं:

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

    **लेकिन अगर आपको एक ऐसा फ़ील्ड जोड़ने की आवश्यकता है जो पहले से ही map का हिस्सा नहीं है तो क्या होगा?**

    मान लीजिए आप फिर से `map1` से शुरू करते हैं, लेकिन भाषा prediction अपने map में नहीं है (कोई `map2` नहीं है)।
    इसके बजाय, यह `lang_id` नामक एक variable में है, और आप जानते हैं कि आप इसके value (`'fr'`) को key `lang` के साथ store करना चाहते हैं।

    आप वास्तव में निम्नलिखित कर सकते हैं:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    यहाँ, `[lang: new_info]` तुरंत एक नया unnamed map बनाता है, और `map1 + ` `map1` को नए unnamed map के साथ merge करता है, जो पहले की तरह समान `new_map` सामग्री उत्पन्न करता है।

    साफ़, है ना?

    **अब आइए इसे Nextflow `channel.map()` operation के संदर्भ में transpose करें।**

    कोड बन जाता है:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    यह निम्नलिखित करता है:

    - `map1, lang_id ->` tuple में दो items लेता है
    - `[map1 + [lang: lang_id]]` ऊपर विस्तृत रूप से नया map बनाता है

    आउटपुट एक single unnamed map है जिसकी सामग्री हमारे ऊपर के उदाहरण में `new_map` के समान है।
    इसलिए हमने प्रभावी रूप से परिवर्तित किया है:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'it'
    ```

    में:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    उम्मीद है कि आप देख सकते हैं कि यदि हम `map1` को `meta` में बदलते हैं, तो मूल रूप से हमारे workflow में हमारे meta map में भाषा predication जोड़ने के लिए यही सब हमें चाहिए।

    एक चीज़ के अलावा!

    हमारे workflow के मामले में, **हमें tuple में `file` object की उपस्थिति के लिए भी account करने की आवश्यकता है**, जो `meta, file, lang_id` से बना है।

    तो यहाँ कोड बन जाएगा:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    यदि आपको यह समझने में कठिनाई हो रही है कि `file` `map` operation में क्यों घूम रहा है, तो कल्पना करें कि `[meta + [lang: lang_id], file]` के बजाय, वह लाइन `[new_map, file]` पढ़ती है।
    इससे यह अधिक स्पष्ट होना चाहिए कि हम केवल `file` को tuple में दूसरी position में अपनी मूल जगह पर छोड़ रहे हैं। हमने बस `new_info` value लिया है और इसे map में fold कर दिया है जो पहली position में है।

    **और यह हमें वापस `tuple val(meta), path(file)` channel structure पर लाता है!**

एक बार जब आपको विश्वास हो जाए कि आप समझ गए हैं कि यह कोड क्या कर रहा है, तो यह देखने के लिए workflow चलाएं कि यह काम किया या नहीं:

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
हमने process के आउटपुट को `meta, file, lang_id` से साफ़ तरीके से पुनर्गठित कर दिया है ताकि `lang_id` अब meta map में keys में से एक हो, और channel के tuples फिर से `meta, file` model के अनुरूप हों।

<!-- TODO (future) क्या हमें यह भी दिखाना चाहिए कि subMap का उपयोग करके key कैसे हटाएं?! या नोट करें कि यह कहाँ मिलेगा। -->

### 2.4. conditionals का उपयोग करके भाषा समूह असाइन करें

अब जब हमारे पास हमारी भाषा predictions हैं, तो आइए कुछ नए groupings असाइन करने के लिए जानकारी का उपयोग करें।

हमारे उदाहरण डेटा में, हमारे characters द्वारा उपयोग की जाने वाली भाषाओं को germanic languages (English, German) और romance language (French, Spanish, Italian) में समूहीकृत किया जा सकता है।
यह classification बाद में pipeline में कहीं आसानी से उपलब्ध होना उपयोगी हो सकता है, तो आइए उस जानकारी को meta map में जोड़ें।

और, अच्छी खबर है, यह एक और मामला है जो `map` ऑपरेटर का उपयोग करने के लिए पूरी तरह से उपयुक्त है!

विशेष रूप से, हम `lang_group` नामक एक variable परिभाषित करने जा रहे हैं, डेटा के प्रत्येक टुकड़े के लिए `lang_group` को क्या value असाइन करनी है यह निर्धारित करने के लिए कुछ साधारण conditional logic का उपयोग करेंगे।

सामान्य syntax इस तरह दिखने वाली है:

```groovy
.map { meta, file ->

    // यहाँ lang_group को परिभाषित करने वाली conditional logic जाती है

    [meta + [lang_group: lang_group], file]
}
```

आप देख सकते हैं कि यह पिछले चरण में हमने उपयोग किए गए on-the-fly map merging operation के समान है।
हमें बस conditional statements लिखने की आवश्यकता है।

यहाँ conditional logic है जो हम लागू करना चाहते हैं:

- default value `'unknown'` के साथ `lang_group` नामक एक variable परिभाषित करें।
- यदि `lang` German (`'de'`) या English (`'en'`) है, तो `lang_group` को `germanic` में बदलें।
- अन्यथा यदि `lang` French (`'fr'`), Spanish (`'es'`) और Italian (`'it'`) युक्त सूची में शामिल है, तो `lang_group` को `romance` में बदलें।

यदि आप पहले से ही जानते हैं कि Nextflow में conditional statements कैसे लिखें तो इसे स्वयं लिखने का प्रयास करें।

!!! tip

    आप map operation के भीतर `meta.lang` के साथ `lang` के value तक पहुँच सकते हैं।

आपको workflow में निम्नलिखित परिवर्तन करने चाहिए:

=== "बाद में"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएँ
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
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएँ
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

यहाँ मुख्य बिंदु हैं:

- हम `lang_group` variable बनाने के लिए `def lang_group = "unknown"` का उपयोग करते हैं जिसकी default value `unknown` पर सेट है।
- हम conditional logic के लिए `if {} else if {}` structure का उपयोग करते हैं, दो germanic भाषाओं के लिए वैकल्पिक `.equals()` tests के साथ, और तीन romance भाषाओं के लिए सूची में अस्तित्व के लिए एक test के साथ।
- हम updated meta map generate करने के लिए पहले की तरह `meta + [lang_group:lang_group]` merge operation का उपयोग करते हैं।

<!-- TODO (future) अतिरिक्त संसाधन अनुभाग में संबंधित docs के लिए नोट/लिंक जोड़ें -->

एक बार जब यह सब समझ में आ जाए, तो परिणाम देखने के लिए फिर से workflow चलाएँ:

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

जैसा कि आप देख सकते हैं, channel elements अपनी `[meta, file]` structure बनाए रखते हैं, लेकिन meta map में अब यह नया classification शामिल है।

### मुख्य बातें

इस खंड में, आपने सीखा कि कैसे:

- **आउटपुट channels में इनपुट मेटाडेटा लागू करें**: इस तरह से मेटाडेटा कॉपी करने से हम बाद में मेटाडेटा सामग्री के आधार पर results को जोड़ सकते हैं।
- **कस्टम keys बनाएँ**: आपने अपने meta map में दो नई keys बनाईं, उन्हें `meta + [new_key:value]` के साथ मौजूदा meta map में merge करते हुए। एक process से computed value के आधार पर, और एक आपके द्वारा `map` ऑपरेटर में सेट की गई condition के आधार पर।

ये आपको अपनी pipeline में आगे बढ़ते समय फ़ाइलों के साथ नए और मौजूदा मेटाडेटा को जोड़ने की अनुमति देते हैं।
भले ही आप process के हिस्से के रूप में मेटाडेटा का उपयोग नहीं कर रहे हों, इस तरह से meta map को डेटा से जुड़ा रखने से सभी प्रासंगिक जानकारी को एक साथ रखना आसान हो जाता है।

---

## 3. process में meta map जानकारी का उपयोग करना

अब जब आप जानते हैं कि meta map कैसे बनाएं और update करें, तो हम वास्तव में मज़ेदार बिट पर आ सकते हैं: वास्तव में process में मेटाडेटा का उपयोग करना।

अधिक विशेष रूप से, हम अपने workflow में एक दूसरा चरण जोड़ने जा रहे हैं ताकि प्रत्येक animal को ASCII art के रूप में draw किया जा सके और उसे speech bubble में recorded text बोलने के लिए कहा जा सके।
हम [`cowpy`](https://github.com/jeffbuttars/cowpy) नामक एक tool का उपयोग करके ऐसा करने जा रहे हैं।

??? info "`cowpy` क्या करता है?"

    `cowpy` एक command-line tool है जो मनमाना text inputs को एक मज़ेदार तरीके से प्रदर्शित करने के लिए ASCII art generate करता है।
    यह Tony Monroe के classic [cowsay](https://en.wikipedia.org/wiki/Cowsay) tool का python implementation है।

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

    वैकल्पिक रूप से, आप default cow के बजाय उपयोग करने के लिए एक character (या 'cowacter') चुन सकते हैं।

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

यदि आपने Hello Nextflow course किया है, तो आप पहले ही इस tool को action में देख चुके हैं।
यदि नहीं, तो चिंता न करें; हम सब कुछ कवर करेंगे जो आपको जानने की आवश्यकता है।

### 3.1. process को import करें और कोड की जाँच करें

हम आपको `COWPY` नामक एक पूर्व-लिखित process मॉड्यूल प्रदान करते हैं जो `cowpy` tool को wrap करता है, इसलिए आपको workflow block से पहले केवल एक include statement जोड़ने की आवश्यकता है।

workflow में निम्नलिखित संपादन करें:

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

आप इसके कोड की जाँच करने के लिए मॉड्यूल फ़ाइल खोल सकते हैं:

```groovy title="modules/cowpy.nf" linen

```
