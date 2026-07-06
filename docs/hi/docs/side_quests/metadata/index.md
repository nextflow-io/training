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
- समझें कि "meta map + data file" interface एक व्यापक रूप से उपयोग किया जाने वाला convention क्यों है
- वर्कफ़्लो execution के दौरान नए मेटाडेटा फ़ील्ड जोड़ें
- आउटपुट को customize करने और व्यवस्थित करने के लिए मेटाडेटा का उपयोग करें

ये कौशल तुम्हें अधिक मज़बूत और लचीली पाइपलाइनें बनाने में मदद करेंगे जो जटिल फ़ाइल संबंधों और प्रोसेसिंग आवश्यकताओं को संभाल सकती हैं।

### पूर्वापेक्षाएँ

इस side quest को शुरू करने से पहले, तुम्हें:

- [Hello Nextflow](../../hello_nextflow/index.md) ट्यूटोरियल या समकक्ष शुरुआती कोर्स पूरा करना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (प्रोसेस, चैनल, ऑपरेटर) का उपयोग करने में सहज होना चाहिए।

---

## 0. शुरू करना

#### ट्रेनिंग codespace खोलो

अगर तुमने अभी तक नहीं किया है, तो [Environment Setup](../../envsetup/index.md) में बताए अनुसार ट्रेनिंग वातावरण खोलना सुनिश्चित करो।

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

Editor प्रोजेक्ट डायरेक्टरी पर फ़ोकस करके खुलता है।

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

हम [`COWPY`](https://github.com/jeffbuttars/cowpy) नामक एक टूल का उपयोग करेंगे ताकि प्रत्येक character की ASCII art बनाई जा सके और उसे अपना recorded अभिवादन बोलते हुए दिखाया जा सके।

??? info "`COWPY` क्या करता है?"

    `COWPY` एक command-line टूल है जो arbitrary text inputs को मज़ेदार तरीके से प्रदर्शित करने के लिए ASCII art generate करता है।
    यह Tony Monroe के classic [cowsay](https://en.wikipedia.org/wiki/Cowsay) टूल का Python implementation है।

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

इसके अलावा, हम `langid` नामक एक भाषा विश्लेषण टूल का उपयोग करेंगे ताकि पहचाना जा सके कि प्रत्येक character कौन सी भाषा बोलता है और पाइपलाइन के आउटपुट को उसी के अनुसार व्यवस्थित किया जा सके।

#### असाइनमेंट की समीक्षा करो

तुम्हारी चुनौती एक Nextflow वर्कफ़्लो लिखना है जो:

1. प्रत्येक character की **ASCII art बनाए**
2. आउटपुट को भाषा परिवार के अनुसार **व्यवस्थित करे** (Germanic बनाम Romance भाषाएँ)

यह एक विशिष्ट वर्कफ़्लो पैटर्न है जहाँ फ़ाइल-विशिष्ट मेटाडेटा प्रोसेसिंग निर्णयों को नियंत्रित करता है; बिल्कुल वही समस्या जिसे metadata maps सुंदर तरीके से हल करते हैं।

#### तैयारी की जाँच सूची

क्या तुम शुरू करने के लिए तैयार हो?

- [ ] मैं इस कोर्स के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता/समझती हूँ
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी working directory उचित रूप से सेट की है
- [ ] मैं असाइनमेंट समझता/समझती हूँ

अगर तुम सभी बॉक्स चेक कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. मेटाडेटा लोड करने और उपयोग करने के बुनियादी विकल्प

`main.nf` वर्कफ़्लो फ़ाइल खोलो और उस वर्कफ़्लो stub की जाँच करो जो हम तुम्हें शुरुआती बिंदु के रूप में दे रहे हैं।

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .view()

    publish:
    cowpy_art = channel.empty()
}

output {
    cowpy_art {
    }
}
```

[`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) ऑपरेटर फ़ाइल की प्रत्येक पंक्ति को चैनल में एक element के रूप में पढ़ता है।
यह वही approach है जिसका उपयोग हम Hello Nextflow (हमारे beginner course) में CSV डेटा लोड करने के लिए करते हैं।
अगर तुम्हें याद दिलाने की ज़रूरत हो कि यह कैसे काम करता है, तो [इस section](../../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file) को देखो।

`header: true` के साथ, पहली पंक्ति को column headers के रूप में माना जाता है, इसलिए प्रत्येक element column name से keyed key-value pairs का एक map बन जाता है।

ध्यान दो कि चूँकि हम अभी तक डेटा पर कोई प्रोसेस नहीं चला रहे, `publish` और `output` blocks केवल stubs हैं।

### 1.1. वर्कफ़्लो चलाओ

वर्कफ़्लो चलाओ ताकि देख सको कि सब कुछ लोड होने के बाद चैनल की सामग्री कैसे structured है:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

जैसा कि तुम देख सकते हो, ऑपरेटर ने CSV फ़ाइल की प्रत्येक पंक्ति के लिए key-value pairs का एक map बनाया है, जिसमें column headers संबंधित values की keys हैं।

प्रत्येक map entry हमारी datasheet के एक कॉलम से मेल खाती है:

- `id`
- `character`
- `recording`

यह बहुत अच्छा है! इससे प्रत्येक पंक्ति से विशिष्ट फ़ील्ड तक पहुँचना आसान हो जाता है।
उदाहरण के लिए, हम `id` से फ़ाइल ID और `recording` से txt फ़ाइल पाथ तक पहुँच सकते हैं।

??? info "(वैकल्पिक) Groovy maps के बारे में अधिक जानकारी"

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
     N E X T F L O W   ~  version 25.10.4

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. `map` से एक विशिष्ट फ़ील्ड चुनना

हम `map` ऑपरेटर का उपयोग करके चैनल में प्रत्येक element पर iterate करेंगे और केवल character फ़ील्ड चुनेंगे, जिसे हम dot notation का उपयोग करके नाम से access कर सकते हैं।

#### 1.2.1. map operation जोड़ो

`character` कॉलम तक पहुँचने के लिए, `.view()` operation से पहले `map` operation इस प्रकार जोड़ो:

=== "बाद में"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()
    ```

किसी विशिष्ट फ़ील्ड तक पहुँचने का यह तरीका Hello Nextflow के [इस section](../../hello_nextflow/02_hello_channels.md#43-use-the-map-operator-to-extract-the-greetings) में अधिक विस्तार से समझाया गया है, अगर तुम्हें याद दिलाने की ज़रूरत हो।

#### 1.2.2. वर्कफ़्लो चलाओ

वर्कफ़्लो चलाओ ताकि verify कर सको कि तुम extracted character names देख सकते हो।

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

इससे पता चलता है कि हम प्रत्येक पंक्ति के लिए `character` कॉलम की values तक पहुँचने में सक्षम हैं।

अब चलो इस डेटा के साथ कुछ करते हैं: `character` और `recording` फ़ील्ड का एक साथ उपयोग करके `COWPY` से ASCII art generate करते हैं।

### 1.3. `multiMap` से sub-channels emit करना

हम तुम्हें एक pre-written `COWPY` process module प्रदान करते हैं, इसलिए पहले तुम्हें प्रोसेस की input requirements की जाँच करनी होगी।

तुम फ़ाइल खोलकर देख सकते हो कि प्रोसेस कैसी दिखती है:

```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 8"
// cowpy से ASCII art generate करें
process COWPY {

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

जैसा कि तुम देख सकते हो, प्रोसेस दो अलग inputs लेती है: एक recording फ़ाइल और एक character नाम।
महत्वपूर्ण बात यह है कि हमारे पास दोनों के लिए values हैं, लेकिन वे वर्तमान में चैनल के प्रत्येक element के अंदर bundled हैं।

Multiple fields को अलग-अलग channels में extract करने का एक तरीका [`multiMap`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap) ऑपरेटर है, जो एक single operation में एक channel को multiple named sub-channels में split करता है।

#### 1.3.1. multiMap operation जोड़ो

`map` operation को `multiMap` से बदलो:

=== "बाद में"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

`multiMap` block प्रत्येक पंक्ति से दो named sub-channels (`file` और `character`) define करता है, जिन्हें हम `ch_datasheet.file` और `ch_datasheet.character` के रूप में access कर सकते हैं।

#### 1.3.2. Sub-channels पर COWPY call करो

अब `COWPY` प्रोसेस को include करो और प्रत्येक sub-channel को एक अलग argument के रूप में दो:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="3 14"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

इससे हम दोनों फ़ील्ड को अलग-अलग pass कर सकते हैं जैसा `COWPY` को चाहिए।

#### 1.3.3. Output publishing सेट करो

अंत में, `COWPY` के output को `publish:` block में जोड़ो:

=== "बाद में"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = COWPY.out
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    ```

इससे हम वर्कफ़्लो द्वारा produce किए गए outputs को आसानी से देख सकेंगे।

#### 1.3.4. वर्कफ़्लो चलाओ

वर्कफ़्लो चलाओ ताकि जाँच सको कि `COWPY` हमारे दिए गए inputs पर चलता है:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_dijkstra] DSL2 - revision: a1b2c3d4e5

    executor >  local (7)
    [3a/f1c290] COWPY (7) [100%] 7 of 7 ✔
    ```

जैसा कि तुम देख सकते हो, `COWPY` प्रत्येक फ़ाइल पर सही character का उपयोग करके चला।

??? abstract "Results डायरेक्टरी सामग्री"

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

??? example "results/cowpy-guten_tag.txt की सामग्री"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

यह approach काम करती है, लेकिन इसकी एक सीमा है: हमें channel को दो अलग sub-channels में split करना पड़ा।
अगर हम प्रोसेस को और अधिक fields pass करना चाहते, तो हमें उन्हें और अधिक sub-channels में split करना पड़ता।
यह परेशान करने वाला और गड़बड़ हो सकता है।

अच्छी खबर: इसे करने का एक सरल तरीका है।

### 1.4. सब कुछ प्रोसेस में एक single input के रूप में group करना

Fields को अलग-अलग channels में split करने के बजाय, हम प्रोसेस को update कर सकते हैं ताकि वह सभी inputs को एक single tuple के रूप में receive करे, जिससे प्रोसेस को call करना सरल हो जाता है।

#### 1.4.1. COWPY प्रोसेस update करो

`COWPY` को प्रत्येक पंक्ति के तीन elements के अनुरूप एक tuple accept करने के लिए update करो:

=== "बाद में"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // cowpy से ASCII art generate करें
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

=== "पहले"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7-8 11 15"
    // cowpy से ASCII art generate करें
    process COWPY {

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

अब प्रोसेस केवल एक input लेती है जिसमें वह सब कुछ है जो हम उसे देना चाहते हैं।

#### 1.4.2. Input tuple बनाने के लिए `map()` का उपयोग करो

हमें अभी भी एक mapping operation का उपयोग करना होगा ताकि उन elements को enumerate किया जा सके जिन्हें हम tuple में प्रोसेस को pass करना चाहते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

तुम सोच सकते हो कि हम `splitCsv` से आने वाले पूरे Groovy map को ऐसे ही क्यों नहीं pass कर सकते।
इसलिए क्योंकि हमें Nextflow को explicitly बताना होगा कि recording फ़ाइल को path के रूप में handle किया जाना चाहिए (यानी उसे ठीक से stage किया जाना चाहिए)।
यह `COWPY` के input interface के स्तर पर होता है, जहाँ `recording` element को explicitly `path` के रूप में designate किया जाता है।

#### 1.4.3. प्रोसेस को call update करो

अंत में, process call में दो अलग inputs को उस single tuple से बदलो जो हमने अभी बनाई:

=== "बाद में"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)
    ```

इससे प्रोसेस को call करना थोड़ा सरल हो जाता है।

#### 1.4.4. वर्कफ़्लो चलाओ

वर्कफ़्लो चलाओ ताकि verify कर सको कि `COWPY` अभी भी डेटा को सही तरीके से process कर सकता है:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [pedantic_lovelace] DSL2 - revision: b2c3d4e5f6

    executor >  local (7)
    [5e/2a1b34] COWPY (7) [100%] 7 of 7 ✔
    ```

आउटपुट पहले जैसी ही सात `cowpy-*.txt` फ़ाइलें हैं, अब `COWPY` को एक सरल call के साथ produce की गई।

??? abstract "Results डायरेक्टरी सामग्री"

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

??? example "results/cowpy-guten_tag.txt की सामग्री"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

यह `multiMap` approach की तुलना में थोड़ा बेहतर है।
लेकिन हमें अभी भी input tuple बनाने के लिए original Groovy map को unpack करना पड़ा, और प्रोसेस और datasheet के बीच एक tight coupling है: `COWPY` की input definition अब column names `id`, `character`, और `recording` को directly reference करती है।

```groovy
input:
tuple val(id), val(character), path(recording)
```

अगर कोई collaborator अलग तरह से structured datasheet का उपयोग करे—अतिरिक्त columns के साथ, या अलग क्रम में columns के साथ—तो यह प्रोसेस बिना modification के काम नहीं करेगी।
इससे प्रोसेस fragile हो जाती है, क्योंकि इसकी input structure datasheet की exact composition से tied है।

इसे solve करने के लिए, हमें एक ऐसा तरीका चाहिए जिससे सभी मेटाडेटा को एक bundle के रूप में pass किया जा सके बिना process interface में उसकी exact structure को hard-code किए।

### 1.5. Meta map + file interface का उपयोग करना

इसका solution यह है कि channel में दो अलग-अलग concerns को अलग किया जाए: **किसी sample के बारे में मेटाडेटा**, और **data file** खुद।
सभी मेटाडेटा को एक single map में bundle करके—"meta map"—हमें एक consistent two-element tuple मिलता है चाहे datasheet में कितने भी metadata columns हों:

```groovy title="Syntax example"
input:
tuple val(meta), path(file)
```

Datasheet में columns जोड़ने या हटाने से `meta` के अंदर की सामग्री बदलती है, लेकिन tuple का shape `[meta, file]` स्थिर रहता है।
इस structure को accept करने वाले processes को यह जानने या परवाह करने की ज़रूरत नहीं कि कितने metadata fields मौजूद हैं।

#### 1.5.1. Tuple की सामग्री को meta map में पुनर्गठित करो

`map` operation को `[meta, file]` tuple produce करने के लिए restructure करते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="5" hl_lines="4 6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // अगले चरण में update करेंगे

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5" hl_lines="4 7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

तुम देखोगे कि हमने एक `view()` statement भी जोड़ा है, `COWPY` call को comment out किया है और `COWPY.out` को `channel.empty()` से बदल दिया है क्योंकि process input definition अभी नई structure से match नहीं करती।

#### 1.5.2. पुनर्गठित सामग्री देखने के लिए वर्कफ़्लो चलाओ

नई channel shape देखने के लिए वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

चैनल में प्रत्येक element अब एक two-element tuple है: पहले meta map, दूसरे स्थान पर फ़ाइल।

```console title="Example element structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

अगर हम बाद में datasheet में एक `language` column जोड़ते हैं, तो वह `meta.language` के रूप में उपलब्ध हो जाएगा बिना process input definition में कोई बदलाव किए।

#### 1.5.3. Meta map का उपयोग करने के लिए `COWPY` प्रोसेस update करो

`COWPY` को `[meta, file]` tuple structure accept करने के लिए update करो:

=== "बाद में"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // cowpy से ASCII art generate करें
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(meta), path(input_file)

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    }
    ```

=== "पहले"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // cowpy से ASCII art generate करें
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

Script block के अंदर, `meta.character` meta map से `character` फ़ील्ड access करता है।
Meta map में कोई भी फ़ील्ड इसी तरह accessible है।

#### 1.5.4. Process call update करो

`COWPY` call restore करो और publishing के लिए इसका output connect करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="5" hl_lines="7 10"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="5" hl_lines="6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // अगले चरण में update करेंगे

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

हमने output publishing भी restore कर दी है।

#### 1.5.5. वर्कफ़्लो चलाओ

वर्कफ़्लो चलाओ ताकि जाँच सको कि सब कुछ काम करता है:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (7)
    [5d/dffd4e] COWPY (7) [100%] 7 of 7 ✔
    ```

Results डायरेक्टरी में अब ASCII art फ़ाइलें हैं।

??? abstract "डायरेक्टरी सामग्री"

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

??? example "results/cowpy-guten_tag.txt की सामग्री"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

प्रोसेस अब `meta` के ज़रिए सभी मेटाडेटा को एक bundle के रूप में receive करती है, जो चाहिए उसका उपयोग करती है (`meta.character`), और बाकी को ignore करती है।

यह [nf-core](https://nf-co.re/) के सभी modules द्वारा उपयोग किया जाने वाला standard interface है।
`tuple val(meta), path(file)` pattern nf-core module library में consistently दिखता है, यही कारण है कि इस convention को adopt करने वाले workflows nf-core modules को minimal friction के साथ swap in कर सकते हैं।

### सारांश

इस section में, तुमने सीखा:

- **Datasheets कैसे पढ़ें:** Header जानकारी के साथ CSV फ़ाइलें parse करने के लिए `splitCsv` का उपयोग करना
- **Meta map convention क्यों मौजूद है:** मेटाडेटा को data files से `[meta, file]` tuples में अलग करने से channel structure stable रहती है जैसे-जैसे datasheet evolve होती है
- **Process के अंदर meta map fields का उपयोग कैसे करें:** Meta map में कोई भी फ़ील्ड script block में dot notation के ज़रिए accessible है

---

## 2. अतिरिक्त मेटाडेटा manipulations

अब जब meta map interface तैयार है, तो हम इसे enrich कर सकते हैं जैसे-जैसे डेटा पाइपलाइन से गुज़रता है।

हम [`langid`](https://github.com/saffsd/langid.py) नामक एक टूल का उपयोग करेंगे ताकि प्रत्येक recording फ़ाइल में मौजूद भाषा की पहचान की जा सके।
टेक्स्ट का एक अंश दिए जाने पर, यह `stdout` पर एक भाषा prediction और एक probability score output करता है।

### 2.1. भाषा पहचान का चरण जोड़ना

हम तुम्हें `IDENTIFY_LANGUAGE` नामक एक pre-written process module प्रदान करते हैं जो `langid` टूल को wrap करता है।

Module फ़ाइल खोलकर उसके कोड की जाँच करो:

```groovy title="modules/langid.nf" linenums="1" hl_lines="7 10"
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

Input definition वही `tuple val(meta), path(file)` structure का उपयोग करती है जो हमने Section 1 में बनाई, इसलिए `ch_datasheet` बिना किसी adaptation के directly इस प्रोसेस में feed हो सकता है।

Output में `stdout` को तीसरे element के रूप में जोड़ा गया है: यह उस language prediction को capture करता है जो `langid` console पर print करता है।
`sed` command probability score और trailing newline को strip करती है, केवल two-letter language code छोड़ती है।

#### 2.1.1. `IDENTIFY_LANGUAGE` को call जोड़ो

`IDENTIFY_LANGUAGE` process module को include करो और datasheet channel पर इसे call करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="4 14-16"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'
    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएं
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()

        COWPY(ch_datasheet)
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)
    ```

इस प्रोसेस का मुख्य output केवल एक string है, इसलिए publish करने के लिए कोई output files नहीं हैं।
इसके बजाय, हम operation के results देखने के लिए `IDENTIFY_LANGUAGE.out.view()` का उपयोग करते हैं।

#### 2.1.2. वर्कफ़्लो चलाओ

Language identification produce करने के लिए वर्कफ़्लो चलाओ, `COWPY` tasks को फिर से चलाने से बचने के लिए `-resume` का उपयोग करो:

```bash
nextflow run main.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (14)
    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

अब हमारे पास dataset में प्रत्येक फ़ाइल के लिए language prediction है।

ध्यान दो कि output tuple `[meta, file, lang_id]` से बना है, यानी meta map और file नए result के साथ-साथ carry होते रहते हैं।

!!! note "नोट"

    Meta map को results के साथ जोड़े रखने का यह पैटर्न बाद में channels में results को join करना आसान बनाता है।
    तुम data को सही तरीके से associate करने के लिए channels में items के क्रम पर भरोसा नहीं कर सकते।
    इसके बजाय, तुम्हें keys का उपयोग करना होगा।
    Meta maps इस उद्देश्य के लिए एक आदर्श संरचना प्रदान करते हैं।

    इस use case को [Splitting & Grouping](../splitting_and_grouping/index.md) side quest में विस्तार से explore किया गया है।

### 2.2. Process outputs से मेटाडेटा बढ़ाना

Language prediction स्वयं फ़ाइल में मौजूद डेटा के बारे में मेटाडेटा है।
इसे एक अलग element के रूप में रखने के बजाय, चलो इसे meta map में वापस fold करते हैं।

#### 2.2.1. एक नया और expanded meta map बनाओ

हम Groovy `+` ऑपरेटर का उपयोग करके original को replace करने के लिए एक नया meta map बना सकते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="14" hl_lines="3-7"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएं
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएं
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

इस operation का मूल है `#!groovy meta + [lang: lang_id]`।

यह code essentially language code वाला एक temporary map बनाता है (`[lang: lang_id]`), फिर Groovy `+` ऑपरेटर का उपयोग करके इसे pre-existing मेटाडेटा वाले original `meta` map के साथ combine करता है, जिससे एक नया और expanded meta map बनता है।

अधिक विस्तृत स्पष्टीकरण के लिए, नीचे दिया गया box देखो।

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

    मान लो तुम फिर से `map1` से शुरू करते हो, लेकिन language prediction अपने map में नहीं है (कोई `map2` नहीं है)।
    इसके बजाय, यह `lang_id` नामक एक variable में है, और तुम जानते हो कि तुम इसकी value (`'fr'`) को `lang` key के साथ store करना चाहते हो।

    तुम वास्तव में निम्नलिखित कर सकते हो:

    ```groovy
    new_map = map1 + [lang: lang_id]
    ```

    यहाँ, `[lang: lang_id]` on the fly एक नया unnamed map बनाता है, और `map1 + ` `map1` को नए unnamed map के साथ merge करता है, जो पहले जैसी ही `new_map` सामग्री उत्पन्न करता है।

    अच्छा है, है ना?

    **अब इसे Nextflow `channel.map()` ऑपरेशन के संदर्भ में transpose करते हैं।**

    कोड बन जाता है:

    ```groovy
    .map { map1, lang_id ->
        map1 + [lang: lang_id]
    }
    ```

    यह निम्नलिखित करता है:

    - `#!groovy map1, lang_id ->` tuple में दो items लेता है
    - `#!groovy map1 + [lang: lang_id]` ऊपर बताए अनुसार नया map बनाता है

    आउटपुट हमारे उदाहरण में `new_map` जैसी ही सामग्री वाला एक single unnamed map है।
    इसलिए हमने effectively transform किया है:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    को:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    उम्मीद है कि तुम देख सकते हो कि अगर हम `map1` को `meta` में बदलते हैं, तो यह मूल रूप से वह सब है जो हमें अपने वर्कफ़्लो में meta map में language prediction जोड़ने के लिए चाहिए।

    सिवाय एक चीज़ के!

    हमारे वर्कफ़्लो के मामले में, **हमें tuple में `file` object की उपस्थिति का भी ध्यान रखना होगा**, जो `meta, file, lang_id` से बना है।

    तो यहाँ कोड बन जाएगा:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    अगर तुम्हें यह समझने में कठिनाई हो रही है कि `map` ऑपरेशन में `file` क्यों इधर-उधर होता दिखता है, तो कल्पना करो कि `#!groovy [meta + [lang: lang_id], file]` के बजाय वह लाइन `[new_map, file]` पढ़ती है।
    इससे यह स्पष्ट होना चाहिए कि हम बस `file` को tuple में दूसरे स्थान पर उसकी मूल जगह पर छोड़ रहे हैं। हमने बस `new_info` value ली और उसे पहले स्थान पर मौजूद map में fold कर दिया।

    **और यह हमें वापस `tuple val(meta), path(file)` चैनल संरचना पर लाता है!**

#### 2.2.2. वर्कफ़्लो चलाओ

एक बार जब तुम्हें विश्वास हो जाए कि तुम समझ गए हो यह कोड क्या कर रहा है, तो वर्कफ़्लो चलाओ और देखो कि यह काम किया या नहीं:

```bash
nextflow run main.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
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

!!! tip "Meta map से keys हटाना"

    तुम Groovy [`subMap`](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)) method का उपयोग करके meta map से एक key हटा सकते हो, जो केवल तुम्हारे द्वारा specify की गई keys वाला एक नया map return करती है:

    ```groovy
    meta.subMap(['id', 'character'])  // केवल 'id' और 'character' वाला map return करता है
    ```

    यह तब उपयोगी है जब downstream process या module को meta map में accumulated सभी fields की ज़रूरत नहीं होती।

### 2.3. Conditionals का उपयोग करके भाषा समूह असाइन करना

Meta map में language prediction के साथ, हम उससे और मेटाडेटा derive कर सकते हैं।
हमारे dataset की भाषाएँ दो families में आती हैं: Germanic (अंग्रेज़ी, जर्मन) और Romance (फ्रेंच, स्पेनिश, इतालवी)।
एक `lang_group` field जोड़ने से वह classification downstream उपलब्ध हो जाएगी।

#### 2.3.1. Conditional logic के साथ एक `map` operation जोड़ो

हम language family असाइन करने के लिए conditional logic के साथ एक दूसरे `map` operation का उपयोग करेंगे:

```groovy
.map { meta, file ->

    // lang_group define करने वाला conditional logic यहाँ जाता है

    [meta + [lang_group: lang_group], file]
}
```

यहाँ वह logic है जिसे apply करना है:

- Default के रूप में `lang_group = 'unknown'` से शुरू करो।
- अगर `meta.lang` `'de'` या `'en'` है, तो `lang_group` को `'germanic'` में set करो।
- Else if `meta.lang` `['fr', 'es', 'it']` में है, तो `lang_group` को `'romance'` में set करो।

!!! tip "सुझाव"

    तुम map operation के भीतर `meta.lang` से `lang` की value तक पहुँच सकते हो।

वर्कफ़्लो में निम्नलिखित बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="14" hl_lines="7-19 21"
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

        ch_languages.view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="14" hl_lines="7"
        // प्रत्येक अभिवादन की भाषा पहचानने के लिए langid चलाएं
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

मुख्य बिंदु:

- `def lang_group = "unknown"` variable को एक safe default के साथ initialize करता है।
- `if / else if` structure दो language families को handle करती है; बाकी सब `'unknown'` रहता है।
- `#!groovy .set { ch_languages }` resulting channel को अगले चरण में उपयोग के लिए एक नाम देता है।

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

#### 2.3.2. वर्कफ़्लो चलाओ:

वर्कफ़्लो चलाओ ताकि verify कर सको कि यह काम करता है:

```bash
nextflow run main.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Meta map में अब चार fields हैं: `id`, `character`, `lang`, और `lang_group`।
Channel structure अभी भी `[meta, file]` है।

### 2.4. आउटपुट को नाम देने और व्यवस्थित करने के लिए मेटाडेटा का उपयोग करना

Meta map में `lang` और `lang_group` अब उपलब्ध होने के साथ, हम उनका उपयोग output file names में language code जोड़ने और उन्हें language family के अनुसार subdirectories में व्यवस्थित करने के लिए कर सकते हैं।

इसके लिए तीन बदलाव चाहिए: `COWPY` प्रोसेस को update करना ताकि वह अपना output rename करे और `meta` को emit करे, `COWPY` call को `ch_languages` पर चलाने के लिए update करना, और output block को subdirectory path specify करने के लिए update करना।

#### 2.4.1. `COWPY` प्रोसेस update करो

Meta map से language code का उपयोग करके output file को rename करो, और output में `meta` जोड़ो ताकि output block subdirectory routing के लिए `lang_group` access कर सके:

=== "बाद में"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        tuple val(meta), path("${meta.lang}-${input_file}")

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
        """
    ```

=== "पहले"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    ```

यह दिखाता है कि हम input definition को बिल्कुल modify किए बिना प्रोसेस के behavior को customize करने के लिए अन्य metadata fields का लाभ कैसे उठा सकते हैं।

#### 2.4.2. `ch_languages` पर चलाने के लिए `COWPY` call update करो

`COWPY(ch_datasheet)` को `COWPY(ch_languages)` से बदलो:

=== "बाद में"

    ```groovy title="main.nf" linenums="32" hl_lines="3"
        .set { ch_languages }

        COWPY(ch_languages)

        publish:
        cowpy_art = COWPY.out
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="32" hl_lines="3 5"
        .set { ch_languages }

        ch_languages.view()

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    }
    ```

हम `ch_languages.view()` line भी हटा देते हैं क्योंकि अब हमें channel contents inspect करने की ज़रूरत नहीं है।

#### 2.4.3. Output block update करो

`output {}` block में एक `path` closure जोड़ो ताकि प्रत्येक फ़ाइल को उसके language group subdirectory में route किया जा सके:

=== "बाद में"

    ```groovy title="main.nf" linenums="40" hl_lines="3"
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="40" hl_lines="2 3"
    output {
        cowpy_art {
        }
    }
    ```

यह दिखाता है कि हम outputs को बड़ी flexibility के साथ व्यवस्थित करने के लिए मेटाडेटा का उपयोग कैसे कर सकते हैं।

#### 2.4.4. पूरी पाइपलाइन चलाओ

पिछले results delete करो और पूरी पाइपलाइन चलाओ:

```bash
rm -r results
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (14)
    [5d/dffd4e] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [e7/317c18] COWPY (7)             [100%] 7 of 7 ✔
    ```

Results डायरेक्टरी अब language family के अनुसार व्यवस्थित है, प्रत्येक फ़ाइल का नाम उसकी detected language के अनुसार है:

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

`output {}` block में `path` closure प्रत्येक `[meta, file]` tuple receive करता है और subdirectory name के रूप में `meta.lang_group` return करता है।
File name खुद वह है जो प्रोसेस output करती है (`#!groovy "${meta.lang}-${input_file}"`)।
मेटाडेटा के दोनों pieces (language code और language group) इस section में बनाए गए enriched meta map से आते हैं।

### सारांश

इस section में, तुमने सीखा:

- **Process outputs से meta map को कैसे बढ़ाएं:** `#!groovy meta + [key: value]` से नई keys जोड़ने से `[meta, file]` channel structure intact रहती है जबकि मेटाडेटा enrich होता है।
- **मेटाडेटा से मेटाडेटा कैसे derive करें:** `map` operation के अंदर conditional logic मौजूदा fields से नए fields compute कर सकती है।
- **Output organization के लिए मेटाडेटा का उपयोग कैसे करें:** `output {}` block में `path` closure meta map से पढ़कर files को subdirectories में route कर सकता है।

---

## 3. Robustness संबंधी विचार

जब metadata values process behavior को drive करती हैं, तो missing या incomplete data ऐसी समस्याएँ पैदा कर सकता है जिन्हें diagnose करना मुश्किल होता है।
यहाँ बताया गया है कि क्या expect करें और इसे कैसे handle करें।

### 3.1. जब एक required metadata field missing हो तो क्या होता है

`COWPY` प्रोसेस के लिए valid result produce करने के लिए `character` value आवश्यक है।
Failure mode इस बात पर निर्भर करता है कि datasheet में column मौजूद है लेकिन खाली है, या बिल्कुल absent है।

#### 3.1.1. Column मौजूद है लेकिन एक value खाली है

मान लो datasheet में एक entry का `character` field blank है:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

Datasheet parse होने पर सभी entries के लिए `character` key बनाई जाती है, लेकिन `sampleA` के लिए `meta.character` एक empty string होगी।
जब Nextflow command में `#!groovy ${meta.character}` substitute करता है, तो `COWPY` टूल को `-c` के लिए एक empty argument मिलता है और वह fail हो जाता है:

??? failure "कमांड आउटपुट"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > fr-bonjour.txt

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

Error message (`expected one argument`) empty `-c` flag की ओर इशारा करता है।
Work directory की `.command.sh` फ़ाइल जाँचने से confirm होता है कि command empty value के साथ चली।

#### 3.1.2. Datasheet में column बिल्कुल मौजूद नहीं है

अगर `character` column पूरी तरह absent है:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

Meta map में `character` key कभी नहीं बनाई जाती।
जब process script `#!groovy ${meta.character}` evaluate करती है, तो missing key `null` return करती है, और Nextflow literally string `null` को command में substitute कर देता है:

??? failure "कमांड आउटपुट"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > de-guten_tag.txt

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

Executed command में `cowpy -c null` diagnostic clue है।

### 3.2. Missing metadata को handle करने की strategies

Workflows को missing metadata के खिलाफ अधिक robust बनाने के लिए दो complementary approaches हैं।

**1. Input validation**

सबसे reliable solution यह है कि कोई भी processing शुरू होने से पहले datasheet को validate किया जाए, ताकि समस्याएँ एक clear error message के साथ जल्दी पकड़ी जाएँ बजाय run के बीच में एक cryptic process failure के रूप में सामने आने के।
[Hello nf-core](../../hello_nf-core/05_input_validation.md) training में nf-schema plugin का उपयोग करके input validation जोड़ने का तरीका cover किया गया है। <!-- TODO (future) pending a proper Validation side quest -->

**2. Required values के लिए explicit process inputs**

अगर तुम चाहते हो कि process interface खुद यह communicate करे कि कोई particular value mandatory है, तो उसे meta map से एक explicit input के रूप में extract करने पर विचार करो:

=== "Process definition"

    ```groovy title="modules/cowpy.nf" linenums="6"
    input:
    tuple val(meta), val(character), path(input_file)
    ```

=== "Workflow call"

    ```groovy title="main.nf"
    COWPY(ch_languages.map { meta, file -> [meta, meta.character, file] })
    ```

यह approach `character` को process contract का एक visible, required हिस्सा बनाती है।
Module पढ़ने वाला कोई भी व्यक्ति तुरंत देख सकता है कि एक character value provide करनी होगी।
अगर field absent है, तो workflow process के चलने से पहले ही channel level पर clearly fail हो जाता है।

यह एक उपयोगी design principle को उजागर करता है:

**Meta map का उपयोग optional या descriptive जानकारी के लिए करो; required values को explicit inputs के रूप में extract करो।**

Meta map channel structures को clean और stable रखता है, लेकिन उन values के लिए जो genuinely किसी प्रोसेस के लिए required हैं, उन्हें named inputs के रूप में surface करने से clarity बेहतर होती है और module को अन्य contexts में सही तरीके से उपयोग करना आसान हो जाता है।

### सारांश

इस section में, तुमने देखा:

- **Missing metadata कैसे manifest होता है:** एक empty field एक empty argument produce करता है; एक absent field `null` produce करता है जो literally command में substitute हो जाता है।
- **दो complementary strategies:** समस्याओं को जल्दी पकड़ने के लिए input validation, और requirements clearly communicate करने के लिए explicit process inputs।

---

## सारांश

इस side quest में, तुमने Nextflow वर्कफ़्लो में मेटाडेटा के साथ प्रभावी ढंग से काम करने का तरीका explore किया।

"Meta map + data file" tuple pattern Nextflow में एक core convention है, जो मेटाडेटा को individual values के रूप में pass करने की तुलना में कई फायदे प्रदान करता है:

- Datasheet evolve होने पर channel structure stable रहती है
- Process behavior को field names hard-code किए बिना प्रत्येक sample के लिए customize किया जा सकता है
- Naming, grouping, और outputs organize करने के लिए मेटाडेटा पूरी पाइपलाइन में उपलब्ध रहता है
- इस interface के लिए लिखे गए modules interchangeable हैं, nf-core modules सहित

### मुख्य पैटर्न

1.  **मेटाडेटा पढ़ना और संरचित करना:** CSV datasheet parse करो और meta map बनाओ।

    ```groovy
    channel.fromPath('datasheet.csv')
        .splitCsv(header: true)
        .map { row ->
            [ [id: row.id, character: row.character], row.recording ]
        }
    ```

2.  **वर्कफ़्लो के दौरान मेटाडेटा expand करना:** Process outputs या derived logic से नई keys जोड़ो।

    ```groovy
    // Process output से
    .map { meta, file, lang ->
        [ meta + [lang: lang], file ]
    }

    // Conditional logic से
    .map { meta, file ->
        def lang_group = "unknown"
        if (meta.lang in ["de", "en"]) { lang_group = "germanic" }
        else if (meta.lang in ["fr", "es", "it"]) { lang_group = "romance" }
        [ meta + [lang_group: lang_group], file ]
    }
    ```

3.  **Process के अंदर मेटाडेटा का उपयोग करना:** Script block में dot notation के ज़रिए किसी भी field को access करो।

    ```groovy
    cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
    ```

4.  **Metadata value के अनुसार outputs व्यवस्थित करना:** `output {}` block में `path` closure का उपयोग करो।

    ```groovy
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

### अतिरिक्त संसाधन

- [map ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#map)
- [multiMap ऑपरेटर](https://www.nextflow.io/docs/latest/reference/operator.html#multimap)
- [stdout output qualifier](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## आगे क्या है?

[Side Quests के मेनू](../index.md) पर वापस जाओ या सूची में अगले विषय पर जाने के लिए पृष्ठ के नीचे दाईं ओर बटन पर क्लिक करो।
