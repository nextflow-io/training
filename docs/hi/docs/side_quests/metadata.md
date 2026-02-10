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

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// cowpy के साथ ASCII art generate करें
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

जैसा कि आप देख सकते हैं, यह process वर्तमान में एक इनपुट फ़ाइल (प्रदर्शित करने के लिए text युक्त) और एक value लेने के लिए डिज़ाइन किया गया है जो ASCII art में draw किए जाने वाले character को निर्दिष्ट करता है, आमतौर पर workflow स्तर पर एक command-line parameter द्वारा प्रदान किया जाता है।

### 3.2. meta map फ़ील्ड को इनपुट के रूप में पास करें

जब हमने Hello Nextflow course में `cowpy` tool का उपयोग किया था, तो हमने अंतिम image बनाने के लिए किस character का उपयोग करना है यह निर्धारित करने के लिए एक command-line parameter का उपयोग किया था।
यह समझ में आया, क्योंकि हम pipeline के प्रति run केवल एक image generate कर रहे थे।

हालांकि, इस tutorial में, हम प्रत्येक subject के लिए एक उपयुक्त image generate करना चाहते हैं जिसे हम process कर रहे हैं, इसलिए command-line parameter का उपयोग बहुत सीमित होगा।

अच्छी खबर: हमारे datasheet में और इसलिए हमारे meta map में एक `character` column है।
आइए इसका उपयोग करें प्रत्येक entry के लिए process द्वारा उपयोग किए जाने वाले character को set करने के लिए।

इसके लिए, हमें तीन काम करने होंगे:

1. पिछले process से आने वाले आउटपुट channel को एक नाम दें ताकि हम इसे अधिक सुविधाजनक रूप से operate कर सकें।
2. रुचि की जानकारी तक पहुँचने का तरीका निर्धारित करें
3. दूसरे process को call जोड़ें और जानकारी को उचित रूप से feed करें।

आइए शुरू करें।

#### 3.2.1. पिछले आउटपुट channel को नाम दें

हमने पिछली manipulations सीधे पहले process के आउटपुट channel, `IDENTIFY_LANGUAGE.out` पर लागू कीं।
उस channel की सामग्री को अगले process को feed करने के लिए (और ऐसा इस तरह से करें जो स्पष्ट और पढ़ने में आसान हो) हम इसे अपना नाम देना चाहते हैं, `ch_languages`।

हम [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) ऑपरेटर का उपयोग करके ऐसा कर सकते हैं।

मुख्य workflow में, `.view()` ऑपरेटर को `.set { ch_languages }` से बदलें, और एक लाइन जोड़ें जो test करती है कि हम channel को नाम से refer कर सकते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
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
            .set { ch_languages }

        // अस्थायी: ch_languages में झाँकें
        ch_languages.view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
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

आइए इसे चलाएँ:

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

यह पुष्टि करता है कि अब हम channel को नाम से refer कर सकते हैं।

#### 3.2.2. फ़ाइल और character मेटाडेटा तक पहुँचें

हम module code देखकर जानते हैं कि `COWPY` process एक text फ़ाइल और एक `character` value प्राप्त करने की उम्मीद करता है।
`COWPY` process call लिखने के लिए, हमें बस यह जानने की आवश्यकता है कि channel में प्रत्येक element से संबंधित फ़ाइल object और मेटाडेटा कैसे extract करें।

जैसा कि अक्सर होता है, ऐसा करने का सबसे सरल तरीका `map` operation का उपयोग करना है।

हमारे channel में `[meta, file]` के रूप में structured tuples हैं, इसलिए हम `file` object तक सीधे पहुँच सकते हैं, और हम meta map के अंदर stored `character` value तक `meta.character` के रूप में refer करके पहुँच सकते हैं।

मुख्य workflow में, निम्नलिखित code परिवर्तन करें:

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

ध्यान दें कि हम `.view` operations के आउटपुट को अधिक readable बनाने के लिए closures (जैसे `{ file -> "File: " + file }`) का उपयोग कर रहे हैं।

आइए इसे चलाएँ:

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

_फ़ाइल paths और character values आपके आउटपुट में अलग क्रम में आ सकते हैं।_

यह पुष्टि करता है कि हम channel में प्रत्येक element के लिए फ़ाइल और character तक पहुँच सकते हैं।

#### 3.2.3. `COWPY` process को call करें

अब सब कुछ एक साथ रखें और वास्तव में `ch_languages` channel पर `COWPY` process को call करें।

मुख्य workflow में, निम्नलिखित code परिवर्तन करें:

=== "बाद में"

    ```groovy title="main.nf" linenums="34"
        // ASCII art generate करने के लिए cowpy चलाएँ
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

आप देखते हैं कि हम बस दो map operations (`.view()` statements को हटाकर) को process call के inputs के रूप में copy करते हैं।
बस उनके बीच comma न भूलना सुनिश्चित करें!

यह थोड़ा बोझिल है, लेकिन हम अगले section में देखेंगे कि इसे कैसे बेहतर बनाया जाए।

आइए इसे चलाएँ:

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

यदि आप results डायरेक्टरी में देखें, तो आपको प्रत्येक अभिवादन का ASCII art वाली अलग-अलग फ़ाइलें दिखनी चाहिए, जो संबंधित character द्वारा बोली गई हैं।

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

यह दिखाता है कि हम pipeline के दूसरे चरण में command को parameterize करने के लिए meta map में जानकारी का उपयोग कर सके।

हालांकि, जैसा कि ऊपर उल्लेख किया गया है, इसमें शामिल कुछ code थोड़ा बोझिल था, क्योंकि हमें workflow body के context में रहते हुए meta data को unpack करना पड़ा।
यह approach meta map से कम संख्या में fields उपयोग करने के लिए ठीक काम करता है, लेकिन यदि हम बहुत अधिक उपयोग करना चाहते तो यह अच्छी तरह scale नहीं करता।

`multiMap()` नामक एक अन्य operator है जो हमें इसे थोड़ा streamline करने की अनुमति देता है, लेकिन तब भी यह आदर्श नहीं है।

??? info "(वैकल्पिक) `multiMap()` के साथ वैकल्पिक version"

    यदि आप सोच रहे हैं, तो हम केवल एक single `map()` operation नहीं लिख सकते थे जो `file` और `character` दोनों output करे, क्योंकि वह उन्हें tuple के रूप में return करता।
    हमें `file` और `character` elements को process को अलग-अलग feed करने के लिए दो separate `map()` operations लिखनी पड़ीं।

    तकनीकी रूप से एक single mapping operation के माध्यम से ऐसा करने का एक और तरीका है, `multiMap()` operator का उपयोग करके, जो multiple channels emit करने में सक्षम है।
    उदाहरण के लिए, आप ऊपर `COWPY` की call को निम्नलिखित code से बदल सकते हैं:

    === "बाद में"

        ```groovy title="main.nf" linenums="34"
            // ASCII art generate करने के लिए cowpy चलाएँ
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "पहले"

        ```groovy title="main.nf" linenums="34"
            // ASCII art generate करने के लिए cowpy चलाएँ
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    यह बिल्कुल वही परिणाम उत्पन्न करता है।

किसी भी स्थिति में, यह अजीब है कि हमें workflow स्तर पर कुछ unpacking करनी पड़ती है।

बेहतर होगा यदि हम पूरी meta map को process में feed कर सकें और वहाँ जाकर जो चाहिए वह चुन सकें।

### 3.3. पूरी meta map पास और उपयोग करें

meta map का मतलब आखिरकार सभी मेटाडेटा को एक bundle के रूप में एक साथ पास करना है।
ऊपर हम ऐसा क्यों नहीं कर सके इसका एकमात्र कारण यह है कि process meta map स्वीकार करने के लिए set up नहीं है।
लेकिन चूंकि हम process code को control करते हैं, हम इसे बदल सकते हैं।

आइए `COWPY` process को `[meta, file]` tuple structure स्वीकार करने के लिए modify करें जो हमने पहले process में उपयोग किया था ताकि हम workflow को streamline कर सकें।

इसके लिए, हमें तीन काम करने होंगे:

1. `COWPY` process module की input definitions को modify करें
2. meta map का उपयोग करने के लिए process command को update करें
3. workflow body में process call को update करें

तैयार? चलिए शुरू करते हैं!

#### 3.3.1. `COWPY` module input को modify करें

`cowpy.nf` module फ़ाइल में निम्नलिखित संपादन करें:

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

यह हमें `[meta, file]` tuple structure का उपयोग करने में सक्षम बनाता है जो हमने tutorial में पहले cover किया था।

ध्यान दें कि हमने meta map output करने के लिए process output definition को update नहीं किया, tutorial को streamlined रखने के लिए, लेकिन `IDENTIFY_LANGUAGE` process के model का अनुसरण करते हुए इसे स्वयं एक exercise के रूप में करने के लिए स्वतंत्र महसूस करें।

#### 3.3.2. meta map field उपयोग करने के लिए command update करें

पूरा meta map अब process के अंदर उपलब्ध है, इसलिए हम command block के अंदर से सीधे इसमें निहित जानकारी को refer कर सकते हैं।

`cowpy.nf` module फ़ाइल में निम्नलिखित संपादन करें:

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

हमने पहले standalone input के रूप में पास किए गए `character` value के reference को meta map में रखे गए value से बदल दिया है, जिसे हम `meta.character` का उपयोग करके refer करते हैं।

अब process call को तदनुसार update करते हैं।

#### 3.3.3. process call update करें और चलाएँ

process अब अपने input के लिए `[meta, file]` tuple structure उपयोग करने की उम्मीद करता है, जो पिछला process output करता है, इसलिए हम बस पूरा `ch_languages` channel `COWPY` process को पास कर सकते हैं।

मुख्य workflow में निम्नलिखित संपादन करें:

=== "बाद में"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // ASCII art generate करने के लिए cowpy चलाएँ
    COWPY(ch_languages)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // ASCII art generate करने के लिए cowpy चलाएँ
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

यह call को काफी सरल बनाता है!

आइए पिछले execution के results हटाएँ और चलाएँ:

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

यदि आप results डायरेक्टरी में देखें, तो आपको पहले जैसे ही outputs दिखने चाहिए, अर्थात् प्रत्येक अभिवादन का ASCII art वाली अलग-अलग फ़ाइलें, संबंधित character द्वारा बोली गई।

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

तो यह सरल code के साथ पहले जैसे ही results produce करता है।

निश्चित रूप से, यह मानता है कि आप process code को modify करने में सक्षम हैं।
कुछ मामलों में, आपको मौजूदा processes पर निर्भर रहना पड़ सकता है जिन्हें आप modify करने के लिए स्वतंत्र नहीं हैं, जो आपके options को सीमित करता है।
अच्छी खबर, यदि आप [nf-core](https://nf-co.re/) project से modules उपयोग करने की योजना बना रहे हैं, तो nf-core modules सभी standard के रूप में `[meta, file]` tuple structure उपयोग करने के लिए set up हैं।

### 3.4. गायब आवश्यक inputs की समस्या निवारण

`COWPY` process के सफलतापूर्वक चलने के लिए `character` value आवश्यक है।
यदि हम configuration फ़ाइल में इसके लिए default value set नहीं करते हैं, तो हमें datasheet में एक value प्रदान करनी होगी।

**क्या होता है यदि हम नहीं करते?**
यह इस पर निर्भर करता है कि input datasheet में क्या है और हम workflow का कौन सा version चला रहे हैं।

#### 3.4.1. character column मौजूद है लेकिन खाली है

मान लीजिए हम data collection error simulate करने के लिए अपनी datasheet में एक entry का character value हटा देते हैं:

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

ऊपर उपयोग किए गए workflow के किसी भी version के लिए, datasheet पढ़ने पर सभी entries के लिए `character` key बनाई जाएगी, लेकिन `sampleA` के लिए value एक empty string होगी।

यह एक error पैदा करेगा।

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

जब Nextflow उस sample के लिए `cowpy` command line चलाता है, तो `${meta.character}` `cowpy` command line में एक empty string से भर जाता है, इसलिए `cowpy` tool एक error फेंकता है जो कहता है कि `-c` argument के लिए कोई value प्रदान नहीं की गई।

#### 3.4.2. character column datasheet में मौजूद नहीं है

अब मान लीजिए हम अपनी datasheet से `character` column को पूरी तरह हटा देते हैं:

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

इस मामले में datasheet पढ़ने पर `character` key बिल्कुल भी नहीं बनाई जाएगी।

##### 3.4.2.1. Workflow स्तर पर value access किया गया

यदि हम section 3.2 में लिखे गए code का version उपयोग कर रहे हैं, तो Nextflow `COWPY` process को call करने से पहले meta map में `character` key तक पहुँचने का प्रयास करेगा।

उसे कोई element नहीं मिलेगा जो instruction से मेल खाता हो, इसलिए वह `COWPY` को बिल्कुल भी नहीं चलाएगा।

??? success "कमांड आउटपुट"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Nextflow के अनुसार, यह workflow सफलतापूर्वक चला!
हालांकि, हमारे वांछित outputs में से कोई भी produce नहीं होगा।

##### 3.4.2.2. Process स्तर पर value access किया गया

यदि हम section 3.3 का version उपयोग कर रहे हैं, तो Nextflow पूरी meta map को `COWPY` process को पास करेगा और command चलाने का प्रयास करेगा।

यह एक error पैदा करेगा, लेकिन पहले मामले की तुलना में एक अलग error।

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

यह इसलिए होता है क्योंकि `meta.character` मौजूद नहीं है, इसलिए इसे access करने का हमारा प्रयास `null` लौटाता है। परिणामस्वरूप, Nextflow शाब्दिक रूप से command-line में `null` डाल देता है, जो निश्चित रूप से `cowpy` tool द्वारा पहचाना नहीं जाता।

#### 3.4.3. समाधान

Workflow configuration के हिस्से के रूप में default value प्रदान करने के अलावा, इसे अधिक robustly handle करने के लिए हम दो काम कर सकते हैं:

1. अपने workflow में input validation implement करें ताकि यह सुनिश्चित हो सके कि datasheet में सभी आवश्यक जानकारी है। आप Hello nf-core training course में [input validation का परिचय](../hello_nf-core/05_input_validation.md) पा सकते हैं। <!-- TODO (future) pending a proper Validation side quest -->

2. यदि आप यह सुनिश्चित करना चाहते हैं कि आपके process module का उपयोग करने वाला कोई भी व्यक्ति तुरंत आवश्यक inputs की पहचान कर सके, तो आप आवश्यक metadata property को एक explicit input भी बना सकते हैं।

यहाँ एक उदाहरण है कि यह कैसे काम करेगा।

सबसे पहले, process स्तर पर, input definition को इस प्रकार update करें:

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

फिर, workflow स्तर पर, metadata से `character` property extract करने और इसे input tuple का एक explicit component बनाने के लिए mapping operation का उपयोग करें:

=== "बाद में"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

इस approach का लाभ यह है कि यह स्पष्ट रूप से दिखाता है कि `character` आवश्यक है, और process को अन्य contexts में redeploy करना आसान बनाता है।

यह एक महत्वपूर्ण design principle को उजागर करता है:

**वैकल्पिक, वर्णनात्मक जानकारी के लिए meta map का उपयोग करें, लेकिन आवश्यक values को explicit inputs के रूप में extract करें।**

meta map channel structures को साफ रखने और मनमानी channel structures को रोकने के लिए उत्कृष्ट है, लेकिन अनिवार्य elements के लिए जो process में सीधे referenced हैं, उन्हें explicit inputs के रूप में extract करने से अधिक robust और maintainable code बनता है।

### मुख्य बातें

इस section में, आपने सीखा कि कैसे मेटाडेटा का उपयोग करके process के execution को customize किया जाता है, इसे workflow स्तर पर या process स्तर पर access करके।

---

## पूरक exercise

यदि आप process के अंदर से meta map जानकारी का उपयोग करने का अभ्यास करना चाहते हैं, तो meta map से `lang` और `lang_group` जैसी अन्य जानकारी का उपयोग करके outputs का नाम और/या संगठन customize करने का प्रयास करें।

उदाहरण के लिए, इस result produce करने के लिए code को modify करने का प्रयास करें:

```console title="Results डायरेक्टरी सामग्री"
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

इस side quest में, आपने Nextflow workflows में मेटाडेटा के साथ प्रभावी ढंग से काम करने का तरीका जाना।

मेटाडेटा को स्पष्ट और डेटा से जुड़ा रखने का यह pattern Nextflow में एक मुख्य best practice है, जो फ़ाइल जानकारी को hardcoding करने पर कई फायदे प्रदान करता है:

- फ़ाइल मेटाडेटा पूरे workflow में फ़ाइलों से जुड़ा रहता है
- प्रत्येक फ़ाइल के लिए process व्यवहार customize किया जा सकता है
- output organization फ़ाइल मेटाडेटा को reflect कर सकता है
- pipeline execution के दौरान फ़ाइल जानकारी expand की जा सकती है

अपने काम में इस pattern को लागू करने से आप robust, maintainable bioinformatics workflows बना सकेंगे।

### मुख्य patterns

1.  **मेटाडेटा पढ़ना और structuring:** CSV फ़ाइलें पढ़ना और organized मेटाडेटा maps बनाना जो आपकी डेटा फ़ाइलों से जुड़ी रहती हैं।

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **workflow के दौरान मेटाडेटा expand करना:** process outputs जोड़कर और conditional logic के माध्यम से values derive करके अपनी pipeline आगे बढ़ने पर अपने मेटाडेटा में नई जानकारी जोड़ना।

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

3.  **Process व्यवहार customize करना:** Process के अंदर मेटाडेटा का उपयोग करना।

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### अतिरिक्त संसाधन

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## आगे क्या?

[Side Quests के मेनू](./index.md) पर वापस लौटें या सूची में अगले topic पर जाने के लिए पृष्ठ के नीचे दाईं ओर बटन पर क्लिक करें।
