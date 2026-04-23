# विभाजन और समूहीकरण

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow डेटा के साथ लचीले तरीके से काम करने के लिए शक्तिशाली टूल प्रदान करता है। एक महत्वपूर्ण क्षमता है डेटा को अलग-अलग streams में विभाजित करना और फिर संबंधित items को वापस एक साथ समूहित करना। यह बायोइनफॉर्मेटिक्स वर्कफ़्लो में विशेष रूप से मूल्यवान है जहाँ तुम्हें विश्लेषण के लिए परिणामों को मिलाने से पहले विभिन्न प्रकार के नमूनों को अलग-अलग प्रोसेस करना होता है।

इसे डाक छाँटने जैसा समझो: तुम पत्रों को गंतव्य के अनुसार अलग करते हो, प्रत्येक ढेर को अलग तरह से प्रोसेस करते हो, फिर एक ही व्यक्ति के पास जाने वाले items को फिर से मिला देते हो। Nextflow वैज्ञानिक डेटा के साथ यह काम करने के लिए विशेष ऑपरेटर का उपयोग करता है। इस दृष्टिकोण को distributed computing और बायोइनफॉर्मेटिक्स वर्कफ़्लो में **scatter/gather** pattern के नाम से भी जाना जाता है।

Nextflow का चैनल सिस्टम इस लचीलेपन के केंद्र में है। चैनल तुम्हारे वर्कफ़्लो के विभिन्न हिस्सों को जोड़ते हैं, जिससे डेटा तुम्हारे विश्लेषण से होकर प्रवाहित होता है। तुम एक ही डेटा स्रोत से कई चैनल बना सकते हो, प्रत्येक चैनल को अलग तरह से प्रोसेस कर सकते हो, और जरूरत पड़ने पर चैनलों को वापस एक साथ मिला सकते हो। यह दृष्टिकोण तुम्हें ऐसे वर्कफ़्लो डिज़ाइन करने देता है जो जटिल बायोइनफॉर्मेटिक्स विश्लेषणों के शाखाओं और मिलन बिंदुओं को स्वाभाविक रूप से दर्शाते हैं।

### सीखने के लक्ष्य

इस side quest में, तुम Nextflow के चैनल ऑपरेटर का उपयोग करके डेटा को विभाजित और समूहित करना सीखोगे।
हम नमूना जानकारी और संबंधित डेटा फ़ाइलों वाली एक CSV फ़ाइल से शुरू करेंगे, फिर इस डेटा को manipulate और पुनर्व्यवस्थित करेंगे।

इस side quest के अंत तक, तुम निम्नलिखित तकनीकों का उपयोग करके डेटा streams को प्रभावी ढंग से अलग और मिला सकोगे:

- `splitCsv` का उपयोग करके फ़ाइलों से डेटा पढ़ना
- `filter` और `map` से डेटा को filter और transform करना
- `join` और `groupTuple` का उपयोग करके संबंधित डेटा को मिलाना
- समानांतर प्रोसेसिंग के लिए `combine` से डेटा combinations बनाना
- `subMap` और deduplication strategies का उपयोग करके डेटा संरचना को optimize करना
- चैनल संरचनाओं को manipulate करने में मदद के लिए named closures के साथ पुन: उपयोग योग्य functions बनाना

ये skills तुम्हें ऐसे वर्कफ़्लो बनाने में मदद करेंगी जो कई इनपुट फ़ाइलों और विभिन्न प्रकार के डेटा को कुशलतापूर्वक संभाल सकें, साथ ही clean और maintainable code structure बनाए रखें।

### पूर्वापेक्षाएँ

इस side quest को शुरू करने से पहले, तुम्हें:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष beginner's course पूरा करना चाहिए।
- बुनियादी Nextflow concepts और mechanisms (प्रोसेस, चैनल, ऑपरेटर, फ़ाइलों के साथ काम करना, meta data) का उपयोग करने में सहज होना चाहिए।

**वैकल्पिक:** हम पहले [Metadata in workflows](../metadata/) side quest पूरा करने की सलाह देते हैं।
वह `splitCsv` के साथ CSV फ़ाइलें पढ़ने और meta maps बनाने की मूल बातें कवर करता है, जिनका हम यहाँ बहुत उपयोग करेंगे।

---

## 0. शुरू करना

#### Training codespace खोलो

अगर तुमने अभी तक नहीं किया है, तो [Environment Setup](../envsetup/index.md) में बताए अनुसार training environment खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Project डायरेक्टरी में जाओ

चलो उस डायरेक्टरी में जाते हैं जहाँ इस ट्यूटोरियल की फ़ाइलें हैं।

```bash
cd side-quests/splitting_and_grouping
```

तुम VSCode को इस डायरेक्टरी पर focus करने के लिए सेट कर सकते हो:

```bash
code .
```

#### सामग्री की समीक्षा करो

तुम्हें एक main workflow फ़ाइल और एक `data` डायरेक्टरी मिलेगी जिसमें `samplesheet.csv` नाम की एक samplesheet है।

```console title="Directory contents"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

Samplesheet में विभिन्न मरीजों के नमूनों की जानकारी है, जिसमें patient ID, sample repeat number, type (normal या tumor), और काल्पनिक डेटा फ़ाइलों के paths शामिल हैं (जो वास्तव में मौजूद नहीं हैं, लेकिन हम मान लेंगे कि हैं)।

```console title="samplesheet.csv"
id,repeat,type,bam
patientA,1,normal,patientA_rep1_normal.bam
patientA,1,tumor,patientA_rep1_tumor.bam
patientA,2,normal,patientA_rep2_normal.bam
patientA,2,tumor,patientA_rep2_tumor.bam
patientB,1,normal,patientB_rep1_normal.bam
patientB,1,tumor,patientB_rep1_tumor.bam
patientC,1,normal,patientC_rep1_normal.bam
patientC,1,tumor,patientC_rep1_tumor.bam
```

इस samplesheet में तीन मरीजों (A, B, C) के आठ नमूने सूचीबद्ध हैं।

प्रत्येक मरीज के लिए, हमारे पास `tumor` type (आमतौर पर tumor biopsies से) या `normal` type (स्वस्थ ऊतक या रक्त से लिए गए) के नमूने हैं।
अगर तुम cancer analysis से परिचित नहीं हो, तो बस यह जानो कि यह एक experimental model से मेल खाता है जो contrastive analyses करने के लिए paired tumor/normal नमूनों का उपयोग करता है।

विशेष रूप से patient A के लिए, हमारे पास दो sets के technical replicates (repeats) हैं।

!!! note "नोट"

    अगर तुम इस experimental design से परिचित नहीं हो तो चिंता मत करो, यह इस ट्यूटोरियल को समझने के लिए महत्वपूर्ण नहीं है।

#### असाइनमेंट की समीक्षा करो

तुम्हारी चुनौती एक Nextflow वर्कफ़्लो लिखना है जो:

1. CSV फ़ाइल से नमूना डेटा **पढ़े** और इसे meta maps के साथ structure करे
2. type (normal vs tumor) के आधार पर नमूनों को अलग-अलग चैनलों में **अलग करे**
3. patient ID और replicate number के आधार पर matched tumor/normal pairs को **join करे**
4. समानांतर प्रोसेसिंग के लिए नमूनों को genomic intervals में **वितरित करे**
5. downstream analysis के लिए संबंधित नमूनों को वापस **समूहित करे**

यह एक सामान्य बायोइनफॉर्मेटिक्स pattern है जहाँ तुम्हें स्वतंत्र प्रोसेसिंग के लिए डेटा को विभाजित करना होता है, फिर तुलनात्मक विश्लेषण के लिए संबंधित items को फिर से मिलाना होता है।

#### तैयारी की जाँच सूची

क्या तुम शुरू करने के लिए तैयार हो?

- [ ] मैं इस course के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता/समझती हूँ
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी working directory उचित रूप से सेट की है
- [ ] मैं असाइनमेंट समझता/समझती हूँ

अगर तुम सभी boxes check कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. नमूना डेटा पढ़ना

### 1.1. `splitCsv` से नमूना डेटा पढ़ना और meta maps बनाना

चलो `splitCsv` से नमूना डेटा पढ़कर और इसे meta map pattern में व्यवस्थित करके शुरू करते हैं। `main.nf` में, तुम देखोगे कि हमने पहले से ही वर्कफ़्लो शुरू कर दिया है।

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note "नोट"

    इस ट्यूटोरियल में, हम सभी चैनल variables के लिए `ch_` prefix का उपयोग करेंगे ताकि यह स्पष्ट रूप से indicate हो कि वे Nextflow चैनल हैं।

अगर तुमने [Metadata in workflows](../metadata/) side quest पूरा किया है, तो तुम इस pattern को पहचानोगे। हम CSV पढ़ने के लिए `splitCsv` का उपयोग करेंगे और तुरंत file paths से metadata को अलग करने के लिए meta map के साथ डेटा को structure करेंगे।

!!! info "जानकारी"

    इस प्रशिक्षण में हम `map` नाम की दो अलग-अलग concepts से मिलेंगे:

    - **Data structure**: Groovy map (अन्य भाषाओं में dictionaries/hashes के समकक्ष) जो key-value pairs store करता है
    - **Channel operator**: `.map()` ऑपरेटर जो चैनल में items को transform करता है

    हम context में स्पष्ट करेंगे कि हमारा मतलब कौन सा है, लेकिन Nextflow के साथ काम करते समय यह अंतर समझना महत्वपूर्ण है।

`main.nf` में ये बदलाव करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="2-6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

यह `splitCsv` operation (headers के साथ CSV पढ़ना) और `map` operation (डेटा को `[meta, file]` tuples के रूप में structure करना) को एक ही step में जोड़ता है। वह बदलाव करो और pipeline चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

अब हमारे पास एक चैनल है जहाँ प्रत्येक item एक `[meta, file]` tuple है - file paths से अलग किया गया metadata। यह संरचना हमें metadata fields के आधार पर अपने workload को विभाजित और समूहित करने देती है।

---

## 2. डेटा को filter और transform करना

### 2.1. `filter` से डेटा filter करना

हम किसी condition के आधार पर डेटा filter करने के लिए [`filter` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#filter) का उपयोग कर सकते हैं। मान लो हम केवल normal नमूनों को प्रोसेस करना चाहते हैं। हम `type` field के आधार पर डेटा filter करके ऐसा कर सकते हैं। चलो इसे `view` ऑपरेटर से पहले insert करते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

filtered result देखने के लिए वर्कफ़्लो फिर से चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

हमने सफलतापूर्वक डेटा को केवल normal नमूनों तक filter किया है। चलो देखते हैं यह कैसे काम करता है।

`filter` ऑपरेटर एक closure लेता है जो चैनल के प्रत्येक element पर apply होता है। अगर closure `true` return करता है, तो element शामिल होता है; अगर `false` return करता है, तो element बाहर हो जाता है।

हमारे case में, हम केवल वे नमूने रखना चाहते हैं जहाँ `meta.type == 'normal'`। Closure tuple `meta,file` का उपयोग प्रत्येक नमूने को refer करने के लिए करता है, `meta.type` से sample type access करता है, और जाँचता है कि यह `'normal'` के बराबर है या नहीं।

यह ऊपर introduce किए गए single closure से accomplish होता है:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. अलग filtered चैनल बनाना

अभी हम filter को सीधे CSV से बनाए गए चैनल पर apply कर रहे हैं, लेकिन हम इसे एक से अधिक तरीकों से filter करना चाहते हैं, इसलिए चलो logic को फिर से लिखते हैं ताकि normal नमूनों के लिए एक अलग filtered चैनल बनाया जा सके:

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="6 8"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

परिणाम देखने के लिए pipeline चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

हमने सफलतापूर्वक डेटा filter किया और normal नमूनों के लिए एक अलग चैनल बनाया।

चलो tumor नमूनों के लिए भी एक filtered चैनल बनाते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="3-8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="3 4"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

हमने normal और tumor नमूनों को दो अलग-अलग चैनलों में अलग किया है, और `view()` को दिए गए closure का उपयोग करके आउटपुट में उन्हें अलग-अलग label किया है: `ch_tumor_samples.view{'Tumor sample: ' + it}`।

### सारांश

इस section में, तुमने सीखा:

- **डेटा filtering**: `filter` से डेटा कैसे filter करें
- **डेटा विभाजन**: किसी condition के आधार पर डेटा को अलग-अलग चैनलों में कैसे विभाजित करें
- **डेटा देखना**: अलग-अलग चैनलों से डेटा print करने और आउटपुट label करने के लिए `view` का उपयोग कैसे करें

हमने अब normal और tumor नमूनों को दो अलग-अलग चैनलों में अलग कर दिया है। अगले section में, हम `id` field पर normal और tumor नमूनों को join करेंगे।

---

## 3. Identifiers द्वारा चैनलों को join करना

पिछले section में, हमने normal और tumor नमूनों को दो अलग-अलग चैनलों में अलग किया। इन्हें उनके type के आधार पर specific processes या वर्कफ़्लो का उपयोग करके स्वतंत्र रूप से प्रोसेस किया जा सकता है। लेकिन क्या होता है जब हम एक ही मरीज के normal और tumor नमूनों की तुलना करना चाहते हैं? इस बिंदु पर, हमें उन्हें वापस एक साथ join करना होगा और यह सुनिश्चित करना होगा कि नमूने उनके `id` field के आधार पर match हों।

Nextflow में चैनलों को combine करने के कई तरीके हैं, लेकिन इस case में सबसे उपयुक्त ऑपरेटर [`join`](https://www.nextflow.io/docs/latest/operator.html#join) है। अगर तुम SQL से परिचित हो, तो यह `JOIN` operation की तरह काम करता है, जहाँ हम join करने की key और join का type specify करते हैं।

### 3.1. Patient ID के आधार पर combine करने के लिए `map` और `join` का उपयोग करना

अगर हम [`join`](https://www.nextflow.io/docs/latest/operator.html#join) documentation देखें, तो हम देख सकते हैं कि default रूप से यह प्रत्येक tuple के पहले item के आधार पर दो चैनलों को join करता है।

#### 3.1.1. Data structure की जाँच करना

अगर console output अभी भी उपलब्ध नहीं है, तो चलो pipeline चलाते हैं और अपनी data structure जाँचते हैं और देखते हैं कि `id` field पर join करने के लिए हमें इसे कैसे modify करना होगा।

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

हम देख सकते हैं कि `id` field प्रत्येक meta map का पहला element है। `join` काम करे इसके लिए, हमें प्रत्येक tuple में `id` field को isolate करना चाहिए। उसके बाद, हम दो चैनलों को combine करने के लिए `join` ऑपरेटर का उपयोग कर सकते हैं।

#### 3.1.2. `id` field को isolate करना

`id` field को isolate करने के लिए, हम [`map` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#map) का उपयोग करके `id` field को पहले element के रूप में एक नया tuple बना सकते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mad_lagrange] DSL2 - revision: 9940b3f23d

    Tumor sample: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Tumor sample: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    Normal sample: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

यह थोड़ा subtle हो सकता है, लेकिन तुम देख सकते हो कि प्रत्येक tuple का पहला element `id` field है।

#### 3.1.3. दो चैनलों को combine करना

अब हम `id` field के आधार पर दो चैनलों को combine करने के लिए `join` ऑपरेटर का उपयोग कर सकते हैं।

एक बार फिर, हम joined outputs print करने के लिए `view` का उपयोग करेंगे।

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_joined_samples = ch_normal_samples
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="7-10"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_wiles] DSL2 - revision: 3bc1979889

    [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

यह थोड़ा पढ़ने में मुश्किल है क्योंकि यह बहुत चौड़ा है, लेकिन तुम देख सकते हो कि नमूने `id` field द्वारा join हो गए हैं। प्रत्येक tuple का format अब इस प्रकार है:

- `id`: नमूना ID
- `normal_meta_map`: Normal नमूने का meta data जिसमें type, replicate और bam फ़ाइल का path शामिल है
- `normal_sample_file`: Normal नमूने की फ़ाइल
- `tumor_meta_map`: Tumor नमूने का meta data जिसमें type, replicate और bam फ़ाइल का path शामिल है
- `tumor_sample`: Tumor नमूना जिसमें type, replicate और bam फ़ाइल का path शामिल है

!!! warning "चेतावनी"

    `join` ऑपरेटर किसी भी un-matched tuples को discard कर देगा। इस उदाहरण में, हमने सुनिश्चित किया कि सभी नमूने tumor और normal के लिए matched थे, लेकिन अगर यह सच नहीं है तो तुम्हें unmatched tuples रखने के लिए `remainder: true` parameter का उपयोग करना होगा। अधिक जानकारी के लिए [documentation](https://www.nextflow.io/docs/latest/operator.html#join) देखो।

तो अब तुम जानते हो कि tuple में किसी field को isolate करने के लिए `map` का उपयोग कैसे करें, और पहले field के आधार पर tuples को combine करने के लिए `join` का उपयोग कैसे करें।
इस ज्ञान के साथ, हम एक shared field के आधार पर चैनलों को सफलतापूर्वक combine कर सकते हैं।

अगले section में, हम उस situation पर विचार करेंगे जहाँ तुम multiple fields पर join करना चाहते हो।

### 3.2. Multiple fields पर join करना

sampleA के लिए हमारे पास 2 replicates हैं, लेकिन sampleB और sampleC के लिए केवल 1। इस case में हम `id` field का उपयोग करके उन्हें प्रभावी ढंग से join कर सके, लेकिन क्या होता अगर वे sync से बाहर होते? हम अलग-अलग replicates के normal और tumor नमूनों को mix up कर सकते थे!

इससे बचने के लिए, हम multiple fields पर join कर सकते हैं। इसे achieve करने के वास्तव में कई तरीके हैं लेकिन हम एक नई joining key बनाने पर focus करेंगे जिसमें नमूना `id` और `replicate` number दोनों शामिल हों।

चलो एक नई joining key बनाकर शुरू करते हैं। हम यह उसी तरह कर सकते हैं जैसे पहले किया था, [`map` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#map) का उपयोग करके `id` और `repeat` fields को पहले element के रूप में एक नया tuple बनाते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

अब हमें दिखना चाहिए कि join `id` और `repeat` दोनों fields का उपयोग करके हो रहा है। वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

    [[patientA, 1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[patientA, 2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[patientB, 1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[patientC, 1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

ध्यान दो कि प्रत्येक joined result के पहले element के रूप में दो elements (`id` और `repeat` fields) का एक tuple है। यह दर्शाता है कि complex items को joining key के रूप में कैसे उपयोग किया जा सकता है, जो एक ही conditions के नमूनों के बीच काफी intricate matching को सक्षम बनाता है।

अगर तुम अलग-अलग keys पर join करने के अधिक तरीके explore करना चाहते हो, तो additional options और examples के लिए [join operator documentation](https://www.nextflow.io/docs/latest/operator.html#join) देखो।

### 3.3. नई joining key बनाने के लिए `subMap` का उपयोग करना

पिछला approach हमारी joining key से field names खो देता है - `id` और `repeat` fields केवल values की एक list बन जाती हैं। बाद में access के लिए field names retain करने के लिए, हम [`subMap` method](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>) का उपयोग कर सकते हैं।

`subMap` method एक map से केवल specified key-value pairs extract करती है। यहाँ हम अपनी joining key बनाने के लिए केवल `id` और `repeat` fields extract करेंगे।

=== "बाद में"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [reverent_wing] DSL2 - revision: 847016c3b7

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

अब हमारे पास एक नई joining key है जो न केवल `id` और `repeat` fields शामिल करती है बल्कि field names भी retain करती है ताकि हम उन्हें बाद में नाम से access कर सकें, जैसे `meta.id` और `meta.repeat`।

### 3.4. map में named closure का उपयोग करना

duplication से बचने और errors कम करने के लिए, हम एक named closure का उपयोग कर सकते हैं। एक named closure हमें एक reusable function बनाने देता है जिसे हम कई जगहों पर call कर सकते हैं।

ऐसा करने के लिए, पहले closure को एक नए variable के रूप में define करो:

=== "बाद में"

    ```groovy title="main.nf" linenums="2" hl_lines="7"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }

        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }

        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

हमने map transformation को एक named variable के रूप में define किया है जिसे हम reuse कर सकते हैं।

ध्यान दो कि हम `file()` का उपयोग करके file path को एक Path object में भी convert करते हैं ताकि इस चैनल को receive करने वाला कोई भी process फ़ाइल को सही तरीके से handle कर सके (अधिक जानकारी के लिए [Working with files](../working_with_files/) देखो)।

चलो अपने वर्कफ़्लो में closure implement करते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "पहले"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

!!! note "नोट"

    `map` ऑपरेटर ने closure pass करने के लिए `{ }` से `( )` में switch किया है। यह इसलिए है क्योंकि `map` ऑपरेटर एक argument के रूप में closure expect करता है और `{ }` का उपयोग anonymous closure define करने के लिए किया जाता है। Named closure call करते समय, `( )` syntax का उपयोग करो।

यह जाँचने के लिए कि सब कुछ अभी भी काम कर रहा है, वर्कफ़्लो एक बार और चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Named closure का उपयोग करने से हम कई जगहों पर एक ही transformation reuse कर सकते हैं, जिससे errors का जोखिम कम होता है और code अधिक readable और maintainable बनता है।

### 3.5. डेटा की duplication कम करना

हमारे वर्कफ़्लो में बहुत सारा duplicated डेटा है। Joined samples के प्रत्येक item में `id` और `repeat` fields repeat होते हैं। चूँकि यह जानकारी grouping key में पहले से उपलब्ध है, हम इस redundancy से बच सकते हैं। याद दिलाने के लिए, हमारी current data structure इस प्रकार दिखती है:

```groovy
[
  [
    "id": "sampleC",
    "repeat": "1",
  ],
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "normal",
  ],
  "sampleC_rep1_normal.bam"
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
  ],
  "sampleC_rep1_tumor.bam"
]
```

चूँकि `id` और `repeat` fields grouping key में उपलब्ध हैं, चलो duplication से बचने के लिए उन्हें प्रत्येक चैनल item के बाकी हिस्से से हटा देते हैं। हम केवल `type` field के साथ एक नया map बनाने के लिए `subMap` method का उपयोग कर सकते हैं। यह approach हमें redundancy को eliminate करते हुए सभी आवश्यक जानकारी बनाए रखने देती है।

=== "बाद में"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

अब closure एक tuple return करता है जहाँ पहले element में `id` और `repeat` fields हैं, और दूसरे element में केवल `type` field है। यह grouping key में एक बार `id` और `repeat` जानकारी store करके redundancy को eliminate करता है, साथ ही सभी आवश्यक जानकारी बनाए रखता है।

देखने के लिए वर्कफ़्लो चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    [[id:patientA, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_tumor.bam]
    ```

हम देख सकते हैं कि हम grouping key में `id` और `repeat` fields केवल एक बार बताते हैं और नमूना डेटा में `type` field है। हमने कोई जानकारी नहीं खोई लेकिन हम अपने चैनल की सामग्री को अधिक संक्षिप्त बनाने में सफल रहे।

### 3.6. Redundant जानकारी हटाना

हमने ऊपर duplicated जानकारी हटाई, लेकिन हमारे चैनलों में अभी भी कुछ अन्य redundant जानकारी है।

शुरुआत में, हमने `filter` का उपयोग करके normal और tumor नमूनों को अलग किया, फिर उन्हें `id` और `repeat` keys के आधार पर join किया। `join` ऑपरेटर tuples को merge करने के क्रम को preserve करता है, इसलिए हमारे case में, left side पर normal नमूनों और right side पर tumor नमूनों के साथ, resulting चैनल इस structure को maintain करता है: `id, <normal elements>, <tumor elements>`।

चूँकि हम अपने चैनल में प्रत्येक element की position जानते हैं, हम `[type:normal]` और `[type:tumor]` metadata को drop करके structure को और simplify कर सकते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

परिणाम देखने के लिए फिर से चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_leavitt] DSL2 - revision: a2303895bd

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

### सारांश

इस section में, तुमने सीखा:

- **Tuples को Manipulate करना**: Tuple में किसी field को isolate करने के लिए `map` का उपयोग कैसे करें
- **Tuples को Join करना**: पहले field के आधार पर tuples को combine करने के लिए `join` का उपयोग कैसे करें
- **Joining Keys बनाना**: नई joining key बनाने के लिए `subMap` का उपयोग कैसे करें
- **Named Closures**: map में named closure का उपयोग कैसे करें
- **Multiple Field Joining**: अधिक precise matching के लिए multiple fields पर join कैसे करें
- **Data Structure Optimization**: Redundant जानकारी हटाकर चैनल structure को कैसे streamline करें

अब तुम्हारे पास एक ऐसा वर्कफ़्लो है जो एक samplesheet को split कर सकता है, normal और tumor नमूनों को filter कर सकता है, उन्हें sample ID और replicate number के आधार पर join कर सकता है, फिर परिणाम print कर सकता है।

यह बायोइनफॉर्मेटिक्स वर्कफ़्लो में एक सामान्य pattern है जहाँ तुम्हें स्वतंत्र रूप से प्रोसेस करने के बाद नमूनों या अन्य प्रकार के डेटा को match करना होता है, इसलिए यह एक उपयोगी skill है। अगले section में, हम एक नमूने को कई बार repeat करने पर ध्यान देंगे।

## 4. Samples को intervals में फैलाना

बायोइनफॉर्मेटिक्स वर्कफ़्लो में एक महत्वपूर्ण pattern है genomic regions में analysis वितरित करना। उदाहरण के लिए, variant calling को genome को intervals (जैसे chromosomes या छोटे regions) में विभाजित करके parallelize किया जा सकता है। यह parallelization strategy multiple cores या nodes में computational load वितरित करके pipeline efficiency में काफी सुधार करती है, जिससे overall execution time कम होता है।

निम्नलिखित section में, हम दिखाएंगे कि अपने नमूना डेटा को multiple genomic intervals में कैसे वितरित करें। हम प्रत्येक नमूने को हर interval के साथ pair करेंगे, जिससे अलग-अलग genomic regions की parallel processing संभव होगी। यह हमारे dataset size को intervals की संख्या से गुणा कर देगा, जिससे कई independent analysis units बनेंगी जिन्हें बाद में वापस एक साथ लाया जा सकता है।

### 4.1. `combine` का उपयोग करके samples को intervals में फैलाना

चलो intervals का एक चैनल बनाकर शुरू करते हैं। जीवन को सरल रखने के लिए, हम केवल 3 intervals का उपयोग करेंगे जिन्हें हम manually define करेंगे। एक real वर्कफ़्लो में, तुम इन्हें एक file input से पढ़ सकते हो या यहाँ तक कि बहुत सारी interval files के साथ एक चैनल बना सकते हो।

=== "बाद में"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

याद रखो, हम प्रत्येक नमूने को प्रत्येक interval के लिए repeat करना चाहते हैं। इसे कभी-कभी samples और intervals का Cartesian product कहा जाता है। हम [`combine` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#combine) का उपयोग करके यह achieve कर सकते हैं। यह channel 1 से प्रत्येक item लेगा और channel 2 के प्रत्येक item के लिए इसे repeat करेगा। चलो अपने वर्कफ़्लो में combine ऑपरेटर जोड़ते हैं:

=== "बाद में"

    ```groovy title="main.nf" linenums="18" hl_lines="3-5"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')

        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="18"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

अब चलाते हैं और देखते हैं क्या होता है:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mighty_tesla] DSL2 - revision: ae013ab70b

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr1]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr2]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr3]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr1]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr2]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr3]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr1]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr2]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr3]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr1]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr2]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr3]
    ```

सफलता! हमने हमारी 3 interval list में प्रत्येक single interval के लिए हर नमूने को repeat किया है। हमने effectively अपने चैनल में items की संख्या तीन गुना कर दी है।

हालाँकि यह पढ़ने में थोड़ा मुश्किल है, इसलिए अगले section में हम इसे व्यवस्थित करेंगे।

### 4.2. चैनल को व्यवस्थित करना

हम अपने नमूना डेटा को tidy और refactor करने के लिए `map` ऑपरेटर का उपयोग कर सकते हैं ताकि इसे समझना आसान हो। चलो intervals string को पहले element पर joining map में move करते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="20" hl_lines="3-9"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="20"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

चलो step by step देखते हैं कि यह map operation क्या करती है।

पहले, हम code को अधिक readable बनाने के लिए named parameters का उपयोग करते हैं। `grouping_key`, `normal`, `tumor` और `interval` नामों का उपयोग करके, हम tuple के elements को index के बजाय नाम से refer कर सकते हैं:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

अगला, हम `grouping_key` को `interval` field के साथ combine करते हैं। `grouping_key` एक map है जिसमें `id` और `repeat` fields हैं। हम Groovy के map addition (`+`) का उपयोग करके `interval` के साथ एक नया map बनाते हैं:

```groovy
                grouping_key + [interval: interval],
```

अंत में, हम इसे तीन elements के tuple के रूप में return करते हैं: combined metadata map, normal नमूना फ़ाइल, और tumor नमूना फ़ाइल:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

चलो फिर से चलाते हैं और चैनल की सामग्री जाँचते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

    [[id:patientA, repeat:1, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

अपने डेटा को सही structure में coerce करने के लिए `map` का उपयोग करना tricky हो सकता है, लेकिन effective data manipulation के लिए यह महत्वपूर्ण है।

अब हमारे पास सभी genomic intervals में repeat किए गए हर नमूने हैं, जिससे कई independent analysis units बनती हैं जिन्हें parallel में प्रोसेस किया जा सकता है। लेकिन क्या होगा अगर हम संबंधित नमूनों को वापस एक साथ लाना चाहते हैं? अगले section में, हम सीखेंगे कि common attributes share करने वाले नमूनों को कैसे group करें।

### सारांश

इस section में, तुमने सीखा:

- **Samples को intervals में फैलाना**: Samples को intervals पर repeat करने के लिए `combine` का उपयोग कैसे करें
- **Cartesian products बनाना**: Samples और intervals के सभी combinations कैसे generate करें
- **चैनल structure व्यवस्थित करना**: बेहतर readability के लिए डेटा restructure करने के लिए `map` का उपयोग कैसे करें
- **Parallel processing की तैयारी**: Distributed analysis के लिए डेटा कैसे setup करें

## 5. `groupTuple` का उपयोग करके samples को aggregate करना

पिछले sections में, हमने सीखा कि input फ़ाइल से डेटा कैसे split करें और specific fields (हमारे case में normal और tumor नमूने) के आधार पर filter करें। लेकिन यह केवल एक प्रकार के joining को cover करता है। क्या होगा अगर हम किसी specific attribute के आधार पर नमूनों को group करना चाहते हैं? उदाहरण के लिए, matched normal-tumor pairs join करने के बजाय, हम "sampleA" के सभी नमूनों को उनके type की परवाह किए बिना एक साथ प्रोसेस करना चाहते हैं। यह pattern बायोइनफॉर्मेटिक्स वर्कफ़्लो में सामान्य है जहाँ तुम efficiency के कारणों से अंत में परिणामों की तुलना या combine करने से पहले संबंधित नमूनों को अलग-अलग प्रोसेस करना चाहते हो।

Nextflow में इसके लिए built-in methods हैं, जिनमें से मुख्य एक हम देखेंगे वह है `groupTuple`।

चलो उन सभी नमूनों को group करके शुरू करते हैं जिनके `id` और `interval` fields समान हैं, यह उस analysis के लिए typical होगा जहाँ हम technical replicates को group करना चाहते हैं लेकिन meaningfully different नमूनों को अलग रखना चाहते हैं।

ऐसा करने के लिए, हमें अपने grouping variables को अलग करना चाहिए ताकि हम उन्हें isolation में उपयोग कर सकें।

पहला step पिछले section में किए गए के समान है। हमें अपने grouping variable को tuple के पहले element के रूप में isolate करना होगा। याद रखो, हमारा पहला element वर्तमान में `id`, `repeat` और `interval` fields का एक map है:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

हम map से अपने `id` और `interval` fields को isolate करने के लिए पहले से `subMap` method का reuse कर सकते हैं। पहले की तरह, हम प्रत्येक नमूने के tuple के पहले element पर `subMap` method apply करने के लिए `map` ऑपरेटर का उपयोग करेंगे।

=== "बाद में"

    ```groovy title="main.nf" linenums="20" hl_lines="11-19"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }

        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="20" hl_lines="10"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

चलो फिर से चलाते हैं और चैनल की सामग्री जाँचते हैं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [hopeful_brenner] DSL2 - revision: 7f4f7fea76

    [[id:patientA, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

हम देख सकते हैं कि हमने सफलतापूर्वक `id` और `interval` fields को isolate किया है, लेकिन अभी तक नमूनों को group नहीं किया है।

!!! note "नोट"

    हम यहाँ `replicate` field को discard कर रहे हैं। यह इसलिए है क्योंकि हमें इसे आगे downstream processing के लिए नहीं चाहिए। इस ट्यूटोरियल को पूरा करने के बाद, देखो कि क्या तुम बाद के grouping को प्रभावित किए बिना इसे शामिल कर सकते हो!

चलो अब [`groupTuple` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#grouptuple) का उपयोग करके इस नए grouping element द्वारा नमूनों को group करते हैं।

=== "बाद में"

    ```groovy title="main.nf" linenums="30" hl_lines="9"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .groupTuple()
              .view()
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="30"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

बस इतना ही है! हमने केवल एक line of code जोड़ी। चलो देखते हैं जब हम इसे चलाते हैं तो क्या होता है:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_jang] DSL2 - revision: a1bee1c55d

    [[id:patientA, interval:chr1], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr2], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr3], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientB, interval:chr1], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr2], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr3], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientC, interval:chr1], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr2], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr3], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    ```

ध्यान दो कि हमारे डेटा की structure बदल गई है और प्रत्येक चैनल element के भीतर फ़ाइलें अब `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]` जैसे tuples में हैं। यह इसलिए है क्योंकि जब हम `groupTuple` का उपयोग करते हैं, तो Nextflow एक group के प्रत्येक नमूने की single files को combine करता है। Downstream में डेटा handle करने की कोशिश करते समय यह याद रखना महत्वपूर्ण है।

!!! note "नोट"

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) groupTuple का विपरीत है। यह चैनल में items को unpack करता है और उन्हें flatten करता है। `transpose` जोड़ने की कोशिश करो और ऊपर किए गए grouping को undo करो!

### सारांश

इस section में, तुमने सीखा:

- **संबंधित नमूनों को group करना**: Common attributes के आधार पर नमूनों को aggregate करने के लिए `groupTuple` का उपयोग कैसे करें
- **Grouping keys को isolate करना**: Grouping के लिए specific fields extract करने के लिए `subMap` का उपयोग कैसे करें
- **Grouped data structures को handle करना**: `groupTuple` द्वारा बनाई गई nested structure के साथ कैसे काम करें
- **Technical replicate handling**: एक ही experimental conditions share करने वाले नमूनों को कैसे group करें

---

## सारांश

इस side quest में, तुमने चैनलों का उपयोग करके डेटा को split और group करना सीखा।

Pipeline से होकर बहते समय डेटा को modify करके, तुम loops या while statements का उपयोग किए बिना एक scalable pipeline बना सकते हो, जो अधिक traditional approaches की तुलना में कई फायदे प्रदान करता है:

- हम बिना किसी additional code के जितने चाहें उतने कम या ज्यादा inputs तक scale कर सकते हैं
- हम iteration के बजाय pipeline से होकर डेटा के flow को handle करने पर focus करते हैं
- हम जितना जरूरी हो उतना complex या simple हो सकते हैं
- Pipeline अधिक declarative बन जाती है, यह focus करते हुए कि क्या होना चाहिए न कि कैसे होना चाहिए
- Nextflow हमारे लिए independent operations को parallel में चलाकर execution optimize करेगा

इन चैनल operations में महारत हासिल करने से तुम flexible, scalable pipelines बना सकोगे जो loops या iterative programming का सहारा लिए बिना complex data relationships को handle करती हैं, जिससे Nextflow execution optimize कर सके और independent operations को automatically parallelize कर सके।

### मुख्य patterns

1.  **Structured input data बनाना:** meta maps के साथ CSV फ़ाइल से शुरू करना ([Metadata in workflows](../metadata/) के patterns पर आधारित)

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **डेटा को अलग चैनलों में split करना:** हमने `type` field के आधार पर डेटा को independent streams में divide करने के लिए `filter` का उपयोग किया

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **Matched नमूनों को join करना:** हमने `id` और `repeat` fields के आधार पर संबंधित नमूनों को recombine करने के लिए `join` का उपयोग किया

    - Key (tuple के पहले element) द्वारा दो चैनलों को join करना

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - Joining key extract करना और इस value द्वारा join करना

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - subMap का उपयोग करके multiple fields पर join करना

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **Intervals में वितरित करना:** हमने parallel processing के लिए genomic intervals के साथ नमूनों के Cartesian products बनाने के लिए `combine` का उपयोग किया।

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Grouping keys द्वारा aggregate करना:** हमने प्रत्येक tuple के पहले element द्वारा group करने के लिए `groupTuple` का उपयोग किया, जिससे `id` और `interval` fields share करने वाले नमूने collect हुए और technical replicates merge हुए।

    ```groovy
    channel.groupTuple()
    ```

6.  **Data structure को optimize करना:** हमने specific fields extract करने के लिए `subMap` का उपयोग किया और transformations को reusable बनाने के लिए named closure बनाया।

    - Map से specific fields extract करना

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - Reusable transformations के लिए named closure का उपयोग करना

    ```groovy
    getSampleIdAndReplicate = { meta, file -> [meta.subMap(['id', 'repeat']), file] }
    channel.map(getSampleIdAndReplicate)
    ```

### अतिरिक्त संसाधन

- [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [join](https://www.nextflow.io/docs/latest/operator.html#join)
- [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
- [combine](https://www.nextflow.io/docs/latest/operator.html#combine)

---

## आगे क्या है?

[Side Quests के menu](../) पर वापस जाओ या list में अगले topic पर जाने के लिए page के नीचे दाईं ओर button पर click करो।
