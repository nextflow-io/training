# मेटाडेटा और मेटा मैप्स

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

किसी भी वैज्ञानिक विश्लेषण में, हम शायद ही कभी केवल रॉ डेटा फ़ाइलों के साथ काम करते हैं।
प्रत्येक फ़ाइल अपनी अतिरिक्त जानकारी के साथ आती है: यह क्या है, यह कहाँ से आई, और इसे क्या खास बनाता है।
इस अतिरिक्त जानकारी को हम मेटाडेटा कहते हैं।

मेटाडेटा अन्य डेटा का वर्णन करने वाला डेटा है।
मेटाडेटा फ़ाइलों और प्रयोगात्मक स्थितियों के बारे में महत्वपूर्ण विवरण ट्रैक करता है, और प्रत्येक डेटासेट की अनूठी विशेषताओं के अनुसार विश्लेषण को अनुकूलित करने में मदद करता है।

इसे एक लाइब्रेरी कैटलॉग की तरह सोचो: जबकि किताबों में वास्तविक सामग्री (रॉ डेटा) होती है, कैटलॉग कार्ड प्रत्येक किताब के बारे में आवश्यक जानकारी प्रदान करते हैं—यह कब प्रकाशित हुई, इसे किसने लिखा, इसे कहाँ खोजना है (मेटाडेटा)।
Nextflow पाइपलाइनों में, मेटाडेटा का उपयोग इसके लिए किया जा सकता है:

- पूरे workflow में फ़ाइल-विशिष्ट जानकारी ट्रैक करना
- फ़ाइल विशेषताओं के आधार पर processes को कॉन्फ़िगर करना
- संयुक्त विश्लेषण के लिए संबंधित फ़ाइलों को समूहित करना

### सीखने के लक्ष्य

इस side quest में, हम workflows में मेटाडेटा को संभालने का तरीका जानेंगे।
बुनियादी फ़ाइल जानकारी वाली एक सरल datasheet (जिसे अक्सर bioinformatics में samplesheet कहा जाता है) से शुरू करते हुए, तुम सीखोगे कि कैसे:

- CSV फ़ाइलों से फ़ाइल मेटाडेटा पढ़ना और parse करना
- मेटाडेटा maps बनाना और उनमें हेरफेर करना
- workflow execution के दौरान नए मेटाडेटा फ़ील्ड जोड़ना
- process व्यवहार को अनुकूलित करने के लिए मेटाडेटा का उपयोग करना

ये कौशल तुम्हें अधिक मजबूत और लचीली पाइपलाइनें बनाने में मदद करेंगे जो जटिल फ़ाइल संबंधों और प्रोसेसिंग आवश्यकताओं को संभाल सकती हैं।

### पूर्वापेक्षाएँ

इस side quest को शुरू करने से पहले, तुम्हें:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष beginner's कोर्स पूरा कर लेना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (processes, channels, operators) का उपयोग करने में सहज होना चाहिए

---

## 0. शुरू करना

#### ट्रेनिंग codespace खोलो

यदि तुमने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित अनुसार ट्रेनिंग environment खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाओ

चलो उस डायरेक्टरी में चलते हैं जहाँ इस ट्यूटोरियल के लिए फ़ाइलें स्थित हैं।

```bash
cd side-quests/metadata
```

तुम VSCode को इस डायरेक्टरी पर फ़ोकस करने के लिए सेट कर सकते हो:

```bash
code .
```

#### सामग्री की समीक्षा करो

तुम्हें एक मुख्य workflow फ़ाइल और एक `data` डायरेक्टरी मिलेगी जिसमें एक datasheet और कुछ डेटा फ़ाइलें हैं।

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

`main.nf` फ़ाइल में workflow एक stub है जिसे तुम धीरे-धीरे एक पूरी तरह से काम करने वाली workflow में विस्तारित करोगे।

datasheet डेटा फ़ाइलों के paths और कुछ संबंधित मेटाडेटा को सूचीबद्ध करती है, जो 3 columns में व्यवस्थित है:

- `id`: स्व-व्याख्यात्मक, फ़ाइल को दी गई एक ID
- `character`: एक character नाम, जिसका उपयोग हम बाद में विभिन्न creatures बनाने के लिए करेंगे
- `data`: `.txt` फ़ाइलों के paths जिनमें विभिन्न भाषाओं में अभिवादन हैं

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

प्रत्येक डेटा फ़ाइल में पाँच भाषाओं में से एक में कुछ अभिवादन text है (fr: French, de: German, es: Spanish, it: Italian, en: English)।

हम तुम्हें `langid` नामक एक containerized भाषा विश्लेषण tool भी प्रदान करेंगे।

#### असाइनमेंट की समीक्षा करो

तुम्हारी चुनौती एक Nextflow workflow लिखना है जो:

1. प्रत्येक फ़ाइल में भाषा को स्वचालित रूप से **पहचानेगी**
2. फ़ाइलों को भाषा परिवार (Germanic बनाम Romance भाषाएँ) द्वारा **समूहित** करेगी
3. प्रत्येक फ़ाइल के लिए प्रोसेसिंग को उसकी भाषा और मेटाडेटा के आधार पर **अनुकूलित** करेगी
4. outputs को भाषा समूह द्वारा **व्यवस्थित** करेगी

यह एक विशिष्ट workflow पैटर्न का प्रतिनिधित्व करता है जहाँ फ़ाइल-विशिष्ट मेटाडेटा प्रोसेसिंग निर्णयों को संचालित करता है; बिल्कुल उस तरह की समस्या जिसे मेटाडेटा maps सुंदरता से हल करते हैं।

#### तैयारी चेकलिस्ट

लगता है कि तुम शुरू करने के लिए तैयार हो?

- [ ] मैं इस कोर्स के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता हूँ
- [ ] मेरा codespace चालू है और चल रहा है
- [ ] मैंने अपनी working डायरेक्टरी उचित रूप से सेट कर ली है
- [ ] मैं असाइनमेंट को समझता हूँ

यदि तुम सभी बॉक्स चेक कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. datasheet से मेटाडेटा लोड करना

`main.nf` workflow फ़ाइल खोलो और workflow stub की जाँच करो जो हम तुम्हें शुरुआती बिंदु के रूप में दे रहे हैं।

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

तुम देख सकते हो कि हमने उदाहरण datasheet को एक फ़ाइल के रूप में लोड करने के लिए एक बुनियादी channel factory सेट अप किया है, लेकिन यह अभी तक फ़ाइल की सामग्री को नहीं पढ़ेगा।
चलो इसे जोड़कर शुरू करते हैं।

### 1.1. `splitCsv` के साथ सामग्री पढ़ना

हमें एक operator चुनना होगा जो हमारी ओर से न्यूनतम प्रयास के साथ फ़ाइल सामग्री को उचित रूप से parse करेगा।
चूंकि हमारी datasheet CSV format में है, यह [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) operator के लिए एक काम है, जो फ़ाइल में प्रत्येक row को channel में एक element के रूप में लोड करता है।

channel construction code में `splitCsv()` operation जोड़ने के लिए निम्नलिखित परिवर्तन करो, साथ ही यह जाँचने के लिए एक `view()` operation जोड़ो कि फ़ाइल की सामग्री channel में सही तरीके से लोड हो रही है।

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

ध्यान दो कि हम `header: true` option का उपयोग कर रहे हैं ताकि Nextflow को CSV फ़ाइल की पहली row को header row के रूप में पढ़ने के लिए कहा जा सके।

चलो देखते हैं कि इससे क्या निकलता है, ठीक है?
workflow चलाओ:

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

हम देख सकते हैं कि operator ने CSV फ़ाइल में प्रत्येक row के लिए key-value pairs का एक map बनाया है, जिसमें column headers संबंधित values के लिए keys के रूप में हैं।

प्रत्येक map entry हमारी datasheet में एक column से मेल खाती है:

- `id`
- `character`
- `recording`

यह बहुत बढ़िया है! यह प्रत्येक फ़ाइल से विशिष्ट फ़ील्ड तक पहुँचना आसान बनाता है।
उदाहरण के लिए, हम `id` के साथ फ़ाइल ID या `recording` के साथ txt फ़ाइल path तक पहुँच सकते हैं।

??? info "(वैकल्पिक) maps के बारे में अधिक"

    Groovy में, वह programming भाषा जिस पर Nextflow बनाया गया है, एक map एक key-value डेटा संरचना है जो Python में dictionaries, JavaScript में objects, या Ruby में hashes के समान है।

    यहाँ एक runnable script है जो दिखाता है कि तुम व्यवहार में एक map को कैसे परिभाषित कर सकते हो और इसकी सामग्री तक कैसे पहुँच सकते हो:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // एक सरल map बनाओ
    def my_map = [id:'sampleA', character:'squirrel']

    // पूरे map को प्रिंट करो
    println "map: ${my_map}"

    // dot notation का उपयोग करके व्यक्तिगत values तक पहुँचो
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    भले ही इसमें एक उचित `workflow` block नहीं है, Nextflow इसे एक workflow की तरह चला सकता है:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    और यहाँ वह है जो तुम output में देखने की उम्मीद कर सकते हो:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. `map` के साथ विशिष्ट फ़ील्ड चुनना

मान लो हम datasheet से `character` column तक पहुँचना और इसे प्रिंट करना चाहते हैं।
हम Nextflow `map` operator का उपयोग करके अपने channel में प्रत्येक item पर iterate कर सकते हैं और विशेष रूप से map object से `character` entry को चुन सकते हैं।

workflow में निम्नलिखित संपादन करो:

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

अब workflow फिर से चलाओ:

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

सफलता! हमने अपनी datasheet से प्राप्त map संरचना का लाभ उठाया है ताकि प्रत्येक row के लिए व्यक्तिगत columns से values तक पहुँच सकें।

अब जब हमने सफलतापूर्वक datasheet पढ़ ली है और प्रत्येक row में डेटा तक पहुँच है, तो हम अपनी pipeline logic को लागू करना शुरू कर सकते हैं।

### 1.3. मेटाडेटा को 'meta map' में व्यवस्थित करना

workflow की वर्तमान स्थिति में, input फ़ाइलें (`recording` key के तहत) और संबंधित मेटाडेटा (`id`, `character`) सभी एक ही स्तर पर हैं, जैसे वे सभी एक बड़े bag में हैं।
व्यावहारिक परिणाम यह है कि इस channel को consume करने वाली प्रत्येक process को इस संरचना को ध्यान में रखते हुए कॉन्फ़िगर करना होगा:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

यह ठीक है जब तक datasheet में columns की संख्या नहीं बदलती।
हालाँकि, यदि तुम datasheet में केवल एक column भी जोड़ते हो, तो channel का shape अब process की अपेक्षा से मेल नहीं खाएगा, और workflow errors उत्पन्न करेगी।
यह process को दूसरों के साथ साझा करना भी कठिन बनाता है जिनके पास थोड़ा अलग input डेटा हो सकता है, और तुम्हें process में variables को hard-code करना पड़ सकता है जो script block द्वारा आवश्यक नहीं हैं।

इस समस्या से बचने के लिए, हमें channel संरचना को consistent रखने का एक तरीका खोजना होगा, चाहे datasheet में कितने भी columns हों।

हम सभी मेटाडेटा को tuple के भीतर एक item में एकत्र करके ऐसा कर सकते हैं, जिसे हम metadata map, या अधिक सरलता से 'meta map' कहेंगे।

`map` operation में निम्नलिखित संपादन करो:

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

हमने अपने channel elements को एक tuple में पुनर्गठित किया है जिसमें दो elements हैं, meta map और संबंधित फ़ाइल object।

चलो workflow चलाते हैं:

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

अब, channel में प्रत्येक element में पहले metadata map और दूसरे संबंधित फ़ाइल object है:

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

परिणामस्वरूप, datasheet में अधिक columns जोड़ने से `meta` map में अधिक मेटाडेटा उपलब्ध होगा, लेकिन channel shape नहीं बदलेगा।
यह हमें ऐसी processes लिखने में सक्षम बनाता है जो channel को consume करती हैं बिना input specification में मेटाडेटा items को hard-code किए:

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

यह Nextflow workflows में मेटाडेटा को व्यवस्थित करने के लिए व्यापक रूप से उपयोग किया जाने वाला convention है।

### सारांश

इस section में, तुमने सीखा है:

- **मेटाडेटा क्यों महत्वपूर्ण है:** अपने डेटा के साथ मेटाडेटा रखने से पूरे workflow में महत्वपूर्ण फ़ाइल जानकारी संरक्षित रहती है।
- **datasheets को कैसे पढ़ें:** header जानकारी के साथ CSV फ़ाइलों को पढ़ने और rows को structured डेटा में बदलने के लिए `splitCsv` का उपयोग करना
- **meta map कैसे बनाएँ:** tuple संरचना `[ [id:value, ...], file ]` का उपयोग करके मेटाडेटा को फ़ाइल डेटा से अलग करना

---

## 2. मेटाडेटा में हेरफेर करना

अब जब हमारा मेटाडेटा लोड हो गया है, चलो इसके साथ कुछ करते हैं!

हम [`langid`](https://github.com/saffsd/langid.py) नामक एक tool का उपयोग करने जा रहे हैं ताकि प्रत्येक creature की recording फ़ाइल में निहित भाषा की पहचान की जा सके।
यह tool भाषाओं के एक set पर pre-trained आता है, और text के एक snippet को देखते हुए, यह एक भाषा prediction और एक संबंधित probability score output करेगा, दोनों `stdout` पर।

### 2.1. process को import करो और code की जाँच करो

हम तुम्हें `IDENTIFY_LANGUAGE` नामक एक pre-written process module प्रदान करते हैं जो `langid` tool को wrap करता है, इसलिए तुम्हें बस workflow block से पहले एक include statement जोड़ना होगा।

workflow में निम्नलिखित संपादन करो:

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

तुम इसके code की जाँच करने के लिए module फ़ाइल खोल सकते हो:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// प्रत्येक input फ़ाइल की भाषा predict करने के लिए langid का उपयोग करो
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

जैसा कि तुम देख सकते हो, input definition उसी `tuple val(meta), path(file)` संरचना का उपयोग करती है जिसे हमने अभी अपने input channel पर लागू किया है।

output definition एक tuple के रूप में structured है जिसकी संरचना input के समान है, सिवाय इसके कि इसमें तीसरे element के रूप में `stdout` भी है।
यह `tuple val(meta), path(file), <output>` पैटर्न मेटाडेटा को input डेटा और outputs दोनों के साथ जुड़ा रखता है क्योंकि यह pipeline के माध्यम से बहता है।

ध्यान दो कि हम यहाँ Nextflow के [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) output qualifier का उपयोग कर रहे हैं क्योंकि tool अपने output को सीधे console पर प्रिंट करता है बजाय एक फ़ाइल लिखने के; और हम command line में `sed` का उपयोग probability score को हटाने, newline characters को हटाकर string को साफ करने, और केवल भाषा prediction return करने के लिए करते हैं।

### 2.2. `IDENTIFY_LANGUAGE` के लिए एक call जोड़ो

अब जब process workflow के लिए उपलब्ध है, तो हम डेटा channel पर इसे चलाने के लिए `IDENTIFY_LANGUAGE` process के लिए एक call जोड़ सकते हैं।

workflow में निम्नलिखित संपादन करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाओ
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

ध्यान दो कि हमने channel construction में मूल `.view()` operation को हटा दिया है।

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

उत्कृष्ट! अब हमारे पास एक prediction है कि प्रत्येक character कौन सी भाषा बोलता है।

और जैसा कि पहले उल्लेख किया गया है, हमने output में input फ़ाइल और meta map भी शामिल किया है, जिसका अर्थ है कि दोनों नई जानकारी के साथ जुड़े रहते हैं जो हमने अभी उत्पन्न की है।
यह अगले step में उपयोगी साबित होगा।

!!! note "नोट"

    अधिक सामान्य रूप से, meta map को results के साथ जुड़ा रखने का यह पैटर्न संबंधित results को associate करना आसान बनाता है जो समान identifiers साझा करते हैं।

    जैसा कि तुम पहले ही सीख चुके होगे, तुम channels में items के order पर भरोसा नहीं कर सकते कि वे उनके बीच results से मेल खाएँ।
    इसके बजाय, तुम्हें डेटा को सही तरीके से associate करने के लिए keys का उपयोग करना होगा, और meta maps इस उद्देश्य के लिए एक आदर्श संरचना प्रदान करते हैं।

    हम [Splitting & Grouping](./splitting_and_grouping.md) side quest में इस use case का विस्तार से पता लगाते हैं।

### 2.3. process outputs के साथ मेटाडेटा बढ़ाना

यह देखते हुए कि हमने अभी जो results उत्पन्न किए हैं वे स्वयं फ़ाइलों की सामग्री के बारे में मेटाडेटा का एक रूप हैं, उन्हें हमारे meta map में जोड़ना उपयोगी होगा।

हालाँकि, हम मौजूदा meta map को in place modify नहीं करना चाहते।
तकनीकी दृष्टिकोण से, ऐसा करना _संभव_ है, लेकिन यह असुरक्षित है।

इसलिए इसके बजाय, हम मौजूदा meta map की सामग्री और नई जानकारी रखने वाली एक नई `lang: lang_id` key-value pair के साथ एक नया meta map बनाएँगे, `+` operator (एक Groovy feature) का उपयोग करके।
और हम इसे पुराने map को नए से बदलने के लिए एक [`map`](https://www.nextflow.io/docs/latest/operator.html#map) operation के साथ combine करेंगे।

यहाँ वे संपादन हैं जो तुम्हें workflow में करने होंगे:

=== "बाद में"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाओ
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाओ
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

यदि तुम अभी तक `+` operator से परिचित नहीं हो, या यदि यह भ्रमित करने वाला लगता है, तो नीचे विस्तृत व्याख्या के माध्यम से कुछ मिनट लो।

??? info "`+` operator का उपयोग करके नए meta map का निर्माण"

    **सबसे पहले, तुम्हें यह जानना होगा कि हम Groovy operator `+` का उपयोग करके दो maps की सामग्री को merge कर सकते हैं।**

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

    **लेकिन क्या होगा यदि तुम्हें एक फ़ील्ड जोड़ना है जो पहले से map का हिस्सा नहीं है?**

    मान लो तुम फिर से `map1` से शुरू करते हो, लेकिन भाषा prediction अपने map में नहीं है (कोई `map2` नहीं है)।
    इसके बजाय, यह `lang_id` नामक एक variable में रखा गया है, और तुम जानते हो कि तुम इसकी value (`'fr'`) को key `lang` के साथ store करना चाहते हो।

    तुम वास्तव में निम्नलिखित कर सकते हो:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    यहाँ, `[lang: new_info]` on the fly एक नया unnamed map बनाता है, और `map1 + ` `map1` को नए unnamed map के साथ merge करता है, पहले जैसी ही `new_map` सामग्री उत्पन्न करता है।

    साफ़, है ना?

    **अब चलो इसे Nextflow `channel.map()` operation के संदर्भ में transpose करते हैं।**

    code बन जाता है:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    यह निम्नलिखित करता है:

    - `map1, lang_id ->` tuple में दो items लेता है
    - `[map1 + [lang: lang_id]]` ऊपर विस्तृत रूप से नया map बनाता है

    output एक single unnamed map है जिसमें हमारे ऊपर के उदाहरण में `new_map` के समान सामग्री है।
    तो हमने प्रभावी रूप से transform किया है:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    में:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    उम्मीद है कि तुम देख सकते हो कि यदि हम `map1` को `meta` में बदलते हैं, तो यह मूल रूप से वह सब कुछ है जो हमें अपने workflow में अपने meta map में भाषा prediction जोड़ने के लिए चाहिए।

    एक चीज़ को छोड़कर!

    हमारे workflow के मामले में, **हमें tuple में `file` object की उपस्थिति के लिए भी account करना होगा**, जो `meta, file, lang_id` से composed है।

    तो यहाँ code बन जाएगा:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    यदि तुम्हें यह समझने में कठिनाई हो रही है कि `file` `map` operation में क्यों घूम रही है, तो कल्पना करो कि `[meta + [lang: lang_id], file]` के बजाय, वह line `[new_map, file]` पढ़ती है।
    इससे यह अधिक स्पष्ट होना चाहिए कि हम बस `file` को tuple में दूसरी position में अपनी मूल जगह पर छोड़ रहे हैं। हमने बस `new_info` value ली है और इसे पहली position में map में fold कर दिया है।

    **और यह हमें `tuple val(meta), path(file)` channel संरचना पर वापस लाता है!**

एक बार जब तुम्हें विश्वास हो जाए कि तुम समझते हो कि यह code क्या कर रहा है, तो देखने के लिए workflow चलाओ कि क्या यह काम किया:

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

हाँ, यह check out करता है!
हमने process के output को `meta, file, lang_id` से साफ़-सुथरे तरीके से पुनर्गठित किया है ताकि `lang_id` अब meta map में keys में से एक हो, और channel के tuples फिर से `meta, file` model में fit हों।

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. conditionals का उपयोग करके एक भाषा समूह assign करना

अब जब हमारे पास हमारी भाषा predictions हैं, चलो कुछ नए groupings assign करने के लिए जानकारी का उपयोग करते हैं।

हमारे उदाहरण डेटा में, हमारे characters द्वारा उपयोग की जाने वाली भाषाओं को germanic भाषाओं (English, German) और romance भाषाओं (French, Spanish, Italian) में grouped किया जा सकता है।
यह उपयोगी हो सकता है कि pipeline में बाद में कहीं वह classification आसानी से उपलब्ध हो, तो चलो उस जानकारी को meta map में जोड़ते हैं।

और, अच्छी खबर, यह फिर से एक case है जो `map` operator का उपयोग करने के लिए पूरी तरह से उपयुक्त है!

विशेष रूप से, हम `lang_group` नामक एक variable define करने जा रहे हैं, डेटा के प्रत्येक piece के लिए `lang_group` को कौन सी value assign करनी है यह निर्धारित करने के लिए कुछ सरल conditional logic का उपयोग करेंगे।

सामान्य syntax इस तरह दिखने वाला है:

```groovy
.map { meta, file ->

    // lang_group को define करने वाली conditional logic यहाँ जाती है

    [meta + [lang_group: lang_group], file]
}
```

तुम देख सकते हो कि यह पिछले step में हमने उपयोग किए गए on-the-fly map merging operation के समान है।
हमें बस conditional statements लिखने की ज़रूरत है।

यहाँ वह conditional logic है जिसे हम apply करना चाहते हैं:

- default value `'unknown'` के साथ `lang_group` नामक एक variable define करो।
- यदि `lang` या तो German (`'de'`) या English (`'en'`) है, तो `lang_group` को `germanic` में बदलो।
- Else यदि `lang` French (`'fr'`), Spanish (`'es'`) और Italian (`'it'`) वाली list में शामिल है, तो `lang_group` को `romance` में बदलो।

यदि तुम पहले से जानते हो कि Nextflow में conditional statements कैसे लिखें तो इसे स्वयं लिखने का प्रयास करो।

!!! tip "सुझाव"

    तुम map operation के भीतर `meta.lang` के साथ `lang` की value तक पहुँच सकते हो।

तुम्हें workflow में निम्नलिखित परिवर्तन करने चाहिए:

=== "बाद में"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाओ
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
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाओ
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

यहाँ मुख्य बिंदु हैं:

- हम default value को `unknown` पर set करके `lang_group` variable बनाने के लिए `def lang_group = "unknown"` का उपयोग करते हैं।
- हम conditional logic के लिए `if {} else if {}` संरचना का उपयोग करते हैं, दो germanic भाषाओं के लिए वैकल्पिक `.equals()` tests के साथ, और तीन romance भाषाओं के लिए एक list में existence के लिए एक test के साथ।
- हम updated meta map generate करने के लिए पहले की तरह `meta + [lang_group:lang_group]` merge operation का उपयोग करते हैं।

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

एक बार जब यह सब समझ में आ जाए, तो परिणाम देखने के लिए workflow फिर से चलाओ:

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

जैसा कि तुम देख सकते हो, channel elements अपनी `[meta, file]` संरचना बनाए रखते हैं, लेकिन meta map में अब यह नया classification शामिल है।

### सारांश

इस section में, तुमने सीखा है कि कैसे:

- **output channels पर input मेटाडेटा apply करें**: इस तरह से मेटाडेटा copy करने से हम बाद में मेटाडेटा सामग्री के आधार पर results को associate कर सकते हैं।
- **custom keys बनाएँ**: तुमने अपने meta map में दो नई keys बनाईं, उन्हें `meta + [new_key:value]` के साथ मौजूदा meta map में merge किया। एक process से computed value के आधार पर, और एक `map` operator में तुमने set की गई condition के आधार पर।

ये तुम्हें अपनी pipeline के माध्यम से progress करते समय फ़ाइलों के साथ नए और मौजूदा मेटाडेटा को associate करने की अनुमति देते हैं।
भले ही तुम किसी process के हिस्से के रूप में मेटाडेटा का उपयोग नहीं कर रहे हो, meta map को इस तरह डेटा के साथ जुड़ा रखने से सभी प्रासंगिक जानकारी को एक साथ रखना आसान हो जाता है।

---

## 3. process में meta map जानकारी का उपयोग करना

अब जब तुम जानते हो कि meta map कैसे बनाएँ और update करें, तो हम वास्तव में मज़ेदार bit पर पहुँच सकते हैं: वास्तव में एक process में मेटाडेटा का उपयोग करना।

अधिक विशेष रूप से, हम अपने workflow में एक दूसरा step जोड़ने जा रहे हैं ताकि प्रत्येक animal को ASCII art के रूप में draw किया जा सके और इसे एक speech bubble में recorded text कहने के लिए बनाया जा सके।
हम [`cowpy`](https://github.com/jeffbuttars/cowpy) नामक एक tool का उपयोग करके ऐसा करने जा रहे हैं।

??? info "`cowpy` क्या करता है?"

    `cowpy` एक command-line tool है जो मज़ेदार तरीके से arbitrary text inputs display करने के लिए ASCII art generate करता है।
    यह Tony Monroe के classic [cowsay](https://en.wikipedia.org/wiki/Cowsay) tool का एक python implementation है।

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

    वैकल्पिक रूप से, तुम default cow के बजाय उपयोग करने के लिए एक character (या 'cowacter') select कर सकते हो।

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

यदि तुमने Hello Nextflow कोर्स के माध्यम से काम किया है, तो तुमने पहले ही इस tool को action में देखा है।
यदि नहीं, तो चिंता मत करो; हम सब कुछ cover करेंगे जो तुम्हें जानना चाहिए जैसे हम आगे बढ़ते हैं।

### 3.1. process को import करो और code की जाँच करो

हम तुम्हें `COWPY` नामक एक pre-written process module प्रदान करते हैं जो `cowpy` tool को wrap करता है, इसलिए तुम्हें बस workflow block से पहले एक include statement जोड़ना होगा।

workflow में निम्नलिखित संपादन करो:

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

तुम इसके code की जाँच करने के लिए module फ़ाइल खोल सकते हो:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// cowpy के साथ ASCII art generate करो
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

जैसा कि तुम देख सकते हो, यह process वर्तमान में एक input फ़ाइल (display किए जाने वाले text वाली) और एक value लेने के लिए designed है जो उस character को specify करती है जिसे ASCII art में draw किया जाना चाहिए, आमतौर पर workflow level पर एक command-line parameter द्वारा प्रदान किया जाता है।

### 3.2. एक meta map फ़ील्ड को input के रूप में pass करना

जब हमने Hello Nextflow कोर्स में `cowpy` tool का उपयोग किया, तो हमने यह निर्धारित करने के लिए एक command-line parameter का उपयोग किया कि final image draw करने के लिए किस character का उपयोग करना है।
यह समझ में आया, क्योंकि हम pipeline के प्रति run केवल एक image generate कर रहे थे।

हालाँकि, इस ट्यूटोरियल में, हम प्रत्येक subject के लिए एक उपयुक्त image generate करना चाहते हैं जिसे हम process कर रहे हैं, इसलिए एक command-line parameter का उपयोग करना बहुत सीमित होगा।

अच्छी खबर: हमारे datasheet में एक `character` column है और इसलिए, हमारे meta map में।
चलो इसका उपयोग करके उस character को set करते हैं जिसे process को प्रत्येक entry के लिए उपयोग करना चाहिए।

इसके लिए, हमें तीन चीजें करनी होंगी:

1. पिछली process से निकलने वाले output channel को एक नाम दो ताकि हम इस पर अधिक सुविधाजनक तरीके से operate कर सकें।
2. निर्धारित करो कि interest की जानकारी तक कैसे पहुँचें
3. दूसरी process के लिए एक call जोड़ो और जानकारी को उचित रूप से feed करो।

चलो शुरू करते हैं।

#### 3.2.1. पिछले output channel को नाम दो

हमने पहली process, `IDENTIFY_LANGUAGE.out` के output channel पर सीधे पिछले manipulations apply किए।
अगली process को channel की सामग्री feed करने के लिए (और ऐसा करने के लिए जो स्पष्ट और पढ़ने में आसान हो) हम इसे अपना नाम देना चाहते हैं, `ch_languages`।

हम [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) operator का उपयोग करके ऐसा कर सकते हैं।

मुख्य workflow में, `.view()` operator को `.set { ch_languages }` से बदलो, और एक line जोड़ो जो test करे कि हम channel को नाम से refer कर सकते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाओ
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

        // अस्थायी: ch_languages में झाँको
        ch_languages.view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाओ
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

चलो इसे चलाते हैं:

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

यह confirm करता है कि अब हम channel को नाम से refer कर सकते हैं।

#### 3.2.2. फ़ाइल और character मेटाडेटा तक पहुँचो

हम module code को देखने से जानते हैं कि `COWPY` process को एक text फ़ाइल और एक `character` value दी जाने की उम्मीद है।
`COWPY` process call लिखने के लिए, हमें बस यह जानना होगा कि channel में प्रत्येक element से संबंधित फ़ाइल object और मेटाडेटा को कैसे extract करें।

जैसा कि अक्सर होता है, ऐसा करने का सबसे सरल तरीका `map` operation का उपयोग करना है।

हमारे channel में tuples हैं जो `[meta, file]` के रूप में structured हैं, इसलिए हम `file` object को सीधे access कर सकते हैं, और हम meta map के अंदर stored `character` value को `meta.character` के रूप में refer करके access कर सकते हैं।

मुख्य workflow में, निम्नलिखित code परिवर्तन करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="34"
        // अस्थायी: फ़ाइल और character तक पहुँचो
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="34"
        // अस्थायी: ch_languages में झाँको
        ch_languages.view()
    ```

ध्यान दो कि हम `.view` operations के output को अधिक readable बनाने के लिए closures (जैसे `{ file -> "File: " + file }`) का उपयोग कर रहे हैं।

चलो इसे चलाते हैं:

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

_तुम्हारे output में फ़ाइल paths और character values एक अलग order में आ सकते हैं।_

यह confirm करता है कि हम channel में प्रत्येक element के लिए फ़ाइल और character तक पहुँचने में सक्षम हैं।

#### 3.2.3. `COWPY` process को call करो

अब चलो इसे सब एक साथ रखते हैं और वास्तव में `ch_languages` channel पर `COWPY` process को call करते हैं।

मुख्य workflow में, निम्नलिखित code परिवर्तन करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="34"
        // ASCII art generate करने के लिए cowpy चलाओ
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="34"
        // अस्थायी: फ़ाइल और character तक पहुँचो
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

तुम देखते हो कि हम बस दो map operations (minus `.view()` statements) को process call के inputs के रूप में copy करते हैं।
बस सुनिश्चित करो कि तुम उनके बीच comma भूल नहीं जाते!

यह थोड़ा clunky है, लेकिन हम अगले section में देखेंगे कि इसे बेहतर कैसे बनाया जाए।

चलो इसे चलाते हैं:

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

यदि तुम results डायरेक्टरी में देखते हो, तो तुम्हें प्रत्येक अभिवादन के ASCII art वाली व्यक्तिगत फ़ाइलें दिखनी चाहिए जो संबंधित character द्वारा बोली गई हैं।

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

यह दिखाता है कि हम pipeline के दूसरे step में command को parameterize करने के लिए meta map में जानकारी का उपयोग करने में सक्षम थे।

हालाँकि, जैसा कि ऊपर उल्लेख किया गया है, शामिल कुछ code थोड़ा clunky था, क्योंकि हमें workflow body के संदर्भ में रहते हुए meta डेटा को unpack करना पड़ा।
यह approach meta map से कम संख्या में फ़ील्ड का उपयोग करने के लिए ठीक काम करता है, लेकिन यदि हम बहुत अधिक का उपयोग करना चाहते हैं तो खराब तरीके से scale होगा।

एक और operator है जिसे `multiMap()` कहा जाता है जो हमें इसे थोड़ा streamline करने की अनुमति देता है, लेकिन तब भी यह आदर्श नहीं है।

??? info "(वैकल्पिक) `multiMap()` के साथ वैकल्पिक संस्करण"

    यदि तुम सोच रहे हो, तो हम एक single `map()` operation नहीं लिख सकते थे जो `file` और `character` दोनों को output करता है, क्योंकि वह उन्हें एक tuple के रूप में return करता।
    हमें `file` और `character` elements को process को अलग से feed करने के लिए दो अलग `map()` operations लिखने पड़े।

    तकनीकी रूप से एक single mapping operation के माध्यम से ऐसा करने का एक और तरीका है, `multiMap()` operator का उपयोग करके, जो कई channels emit करने में सक्षम है।
    उदाहरण के लिए, तुम ऊपर `COWPY` के call को निम्नलिखित code से बदल सकते हो:

    === "बाद में"

        ```groovy title="main.nf" linenums="34"
            // ASCII art generate करने के लिए cowpy चलाओ
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "पहले"

        ```groovy title="main.nf" linenums="34"
            // ASCII art generate करने के लिए cowpy चलाओ
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    यह बिल्कुल वही परिणाम उत्पन्न करता है।

किसी भी मामले में, यह अजीब है कि हमें workflow level पर कुछ unpacking करना पड़ता है।

यह बेहतर होगा यदि हम पूरे meta map को process में feed कर सकें और वहाँ पहुँचने के बाद जो हमें चाहिए उसे चुन सकें।

### 3.3. पूरे meta map को pass और उपयोग करना

meta map का point आखिरकार सभी मेटाडेटा को एक bundle के रूप में एक साथ pass करना है।
एकमात्र कारण जिससे हम ऊपर ऐसा नहीं कर सके वह यह है कि process एक meta map स्वीकार करने के लिए set up नहीं है।
लेकिन चूंकि हम process code को control करते हैं, हम इसे बदल सकते हैं।

चलो `COWPY` process को `[meta, file]` tuple संरचना स्वीकार करने के लिए modify करते हैं जिसका उपयोग हमने पहली process में किया था ताकि हम workflow को streamline कर सकें।

इसके लिए, हमें तीन चीजें करनी होंगी:

1. `COWPY` process module के input definitions को modify करो
2. meta map का उपयोग करने के लिए process command को update करो
3. workflow body में process call को update करो

तैयार? चलो चलते हैं!

#### 3.3.1. `COWPY` module input को modify करो

`cowpy.nf` module फ़ाइल में निम्नलिखित संपादन करो:

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

यह हमें `[meta, file]` tuple संरचना का उपयोग करने में सक्षम बनाता है जिसे हमने ट्यूटोरियल में पहले cover किया था।

ध्यान दो कि हमने meta map को output करने के लिए process output definition को update नहीं किया, ताकि ट्यूटोरियल को streamlined रखा जा सके, लेकिन `IDENTIFY_LANGUAGE` process के model का पालन करते हुए स्वयं ऐसा करने के लिए स्वतंत्र महसूस करो।

#### 3.3.2. meta map फ़ील्ड का उपयोग करने के लिए command को update करो

पूरा meta map अब process के अंदर उपलब्ध है, इसलिए हम command block के अंदर से सीधे इसमें निहित जानकारी को refer कर सकते हैं।

`cowpy.nf` module फ़ाइल में निम्नलिखित संपादन करो:

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

हमने पहले standalone input के रूप में pass की गई `character` value के reference को meta map में held value से बदल दिया है, जिसे हम `meta.character` का उपयोग करके refer करते हैं।

अब चलो तदनुसार process call को update करते हैं।

#### 3.3.3. process call को update करो और इसे चलाओ

process अब अपने input को `[meta, file]` tuple संरचना का उपयोग करने की उम्मीद करती है, जो कि पिछली process outputs है, इसलिए हम बस पूरे `ch_languages` channel को `COWPY` process को pass कर सकते हैं।

मुख्य workflow में निम्नलिखित संपादन करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // ASCII art generate करने के लिए cowpy चलाओ
    COWPY(ch_languages)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // ASCII art generate करने के लिए cowpy चलाओ
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

यह call को काफी simplify करता है!

चलो पिछले execution के results को delete करते हैं और इसे चलाते हैं:

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

यदि तुम results डायरेक्टरी में देखते हो, तो तुम्हें पहले जैसे ही outputs दिखने चाहिए, _यानी_ प्रत्येक अभिवादन के ASCII art वाली व्यक्तिगत फ़ाइलें जो संबंधित character द्वारा बोली गई हैं।

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

तो यह सरल code के साथ पहले जैसे ही results उत्पन्न करता है।

बेशक, यह मानता है कि तुम process code को modify करने में सक्षम हो।
कुछ मामलों में, तुम्हें मौजूदा processes पर भरोसा करना पड़ सकता है जिन्हें तुम modify करने के लिए स्वतंत्र नहीं हो, जो तुम्हारे options को सीमित करता है।
अच्छी खबर, यदि तुम [nf-core](https://nf-co.re/) प्रोजेक्ट से modules का उपयोग करने की योजना बना रहे हो, तो यह है कि nf-core modules सभी `[meta, file]` tuple संरचना को एक standard के रूप में उपयोग करने के लिए set up हैं।

### 3.4. missing required inputs की troubleshooting

`COWPY` process को सफलतापूर्वक चलाने के लिए `character` value आवश्यक है।
यदि हम configuration फ़ाइल में इसके लिए एक default value set नहीं करते हैं, तो हमें datasheet में इसके लिए एक value प्रदान करनी होगी।

**क्या होता है यदि हम नहीं करते?**
यह इस बात पर निर्भर करता है कि input datasheet में क्या है और हम workflow के किस version को चला रहे हैं।

#### 3.4.1. character column मौजूद है लेकिन खाली है

मान लो हम डेटा collection error को simulate करने के लिए अपनी datasheet में entries में से एक के लिए character value को delete करते हैं:

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

हमने ऊपर उपयोग किए गए workflow के किसी भी version के लिए, जब datasheet पढ़ी जाती है तो सभी entries के लिए `character` key बनाई जाएगी, लेकिन `sampleA` के लिए value एक empty string होगी।

यह एक error का कारण बनेगा।

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

जब Nextflow उस sample के लिए `cowpy` command line चलाता है, तो `${meta.character}` को `cowpy` command line में एक empty string से भरा जाता है, इसलिए `cowpy` tool एक error throw करता है जो कहता है कि `-c` argument के लिए कोई value प्रदान नहीं की गई थी।

#### 3.4.2. character column datasheet में मौजूद नहीं है

अब मान लो हम अपनी datasheet से `character` column को पूरी तरह से delete करते हैं:

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

इस मामले में जब datasheet पढ़ी जाती है तो `character` key बिल्कुल नहीं बनाई जाएगी।

##### 3.4.2.1. workflow level पर value accessed

यदि हम section 3.2 में लिखे गए code के version का उपयोग कर रहे हैं, तो Nextflow `COWPY` process को call करने से पहले meta map में `character` key तक पहुँचने का प्रयास करेगा।

यह कोई elements नहीं पाएगा जो instruction से मेल खाते हैं, इसलिए यह `COWPY` को बिल्कुल नहीं चलाएगा।

??? success "कमांड आउटपुट"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

जहाँ तक Nextflow का सवाल है, यह workflow सफलतापूर्वक चली!
हालाँकि, हम जो outputs चाहते हैं उनमें से कोई भी उत्पन्न नहीं होगा।

##### 3.4.2.2. process level पर value accessed

यदि हम section 3.3 में version का उपयोग कर रहे हैं, तो Nextflow पूरे meta map को `COWPY` process को pass करेगा और command चलाने का प्रयास करेगा।

यह एक error का कारण बनेगा, लेकिन पहले case की तुलना में एक अलग।

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

ऐसा इसलिए होता है क्योंकि `meta.character` मौजूद नहीं है, इसलिए इसे access करने का हमारा प्रयास `null` return करता है। परिणामस्वरूप, Nextflow literally command-line में `null` plug in करता है, जो बेशक `cowpy` tool द्वारा recognized नहीं है।

#### 3.4.3. समाधान

workflow configuration के हिस्से के रूप में एक default value supply करने के अलावा, इसे अधिक robustly handle करने के लिए हम दो चीजें कर सकते हैं:

1. यह सुनिश्चित करने के लिए अपने workflow में input validation implement करो कि datasheet में सभी आवश्यक जानकारी है। तुम Hello nf-core ट्रेनिंग कोर्स में [input validation का परिचय](../hello_nf-core/05_input_validation.md) पा सकते हो। <!-- TODO (future) pending a proper Validation side quest -->

2. यदि तुम यह सुनिश्चित करना चाहते हो कि कोई भी जो तुम्हारे process module का उपयोग करता है वह तुरंत required inputs की पहचान कर सके, तो तुम required मेटाडेटा property को एक explicit input भी बना सकते हो।

यहाँ एक उदाहरण है कि यह कैसे काम करेगा।

सबसे पहले, process level पर, input definition को निम्नानुसार update करो:

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

फिर, workflow level पर, मेटाडेटा से `character` property को extract करने और इसे input tuple का एक explicit component बनाने के लिए एक mapping operation का उपयोग करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

इस approach का फायदा यह है कि यह explicitly दिखाता है कि `character` required है, और process को अन्य contexts में redeploy करना आसान बनाता है।

यह एक महत्वपूर्ण design principle को highlight करता है:

**optional, descriptive जानकारी के लिए meta map का उपयोग करो, लेकिन required values को explicit inputs के रूप में extract करो।**

meta map channel संरचनाओं को clean रखने और arbitrary channel संरचनाओं को रोकने के लिए उत्कृष्ट है, लेकिन mandatory elements के लिए जो सीधे एक process में referenced हैं, उन्हें explicit inputs के रूप में extract करना अधिक robust और maintainable code बनाता है।

### सारांश

इस section में, तुमने सीखा है कि एक process के execution को customize करने के लिए मेटाडेटा का उपयोग कैसे करें, इसे workflow level पर या process level पर access करते हुए।

---

## पूरक अभ्यास

यदि तुम एक process के अंदर से meta map जानकारी का उपयोग करने का अभ्यास करना चाहते हो, तो meta map से अन्य जानकारी के pieces जैसे `lang` और `lang_group` का उपयोग करके outputs को नाम दिए जाने और/या व्यवस्थित किए जाने के तरीके को customize करने का प्रयास करो।

उदाहरण के लिए, यह परिणाम उत्पन्न करने के लिए code को modify करने का प्रयास करो:

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

इस side quest में, तुमने Nextflow workflows में मेटाडेटा के साथ प्रभावी ढंग से काम करने का तरीका explored किया है।

मेटाडेटा को explicit और डेटा के साथ attached रखने का यह पैटर्न Nextflow में एक core best practice है, जो फ़ाइल जानकारी को hardcoding करने की तुलना में कई फायदे प्रदान करता है:

- फ़ाइल मेटाडेटा पूरे workflow में फ़ाइलों के साथ जुड़ा रहता है
- process व्यवहार को प्रति फ़ाइल customized किया जा सकता है
- output organization फ़ाइल मेटाडेटा को reflect कर सकता है
- pipeline execution के दौरान फ़ाइल जानकारी को expanded किया जा सकता है

अपने काम में इस पैटर्न को apply करने से तुम्हें robust, maintainable bioinformatics workflows बनाने में सक्षम होगा।

### मुख्य पैटर्न

1.  **मेटाडेटा को पढ़ना और संरचित करना:** CSV फ़ाइलों को पढ़ना और organized मेटाडेटा maps बनाना जो तुम्हारी डेटा फ़ाइलों के साथ जुड़े रहते हैं।

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **Workflow के दौरान मेटाडेटा का विस्तार करना** जैसे-जैसे तुम्हारी pipeline progress करती है, process outputs जोड़कर और conditional logic के माध्यम से values derive करके अपने मेटाडेटा में नई जानकारी जोड़ना।

    - process output के आधार पर नई keys जोड़ना

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - conditional clause का उपयोग करके नई keys जोड़ना

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

3.  **Process व्यवहार को अनुकूलित करना:** process के अंदर मेटाडेटा का उपयोग करना।

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### अतिरिक्त संसाधन

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## आगे क्या है?

[Side Quests के menu](./index.md) पर वापस जाओ या list में अगले topic पर जाने के लिए पृष्ठ के निचले दाएँ कोने में button पर क्लिक करो।
