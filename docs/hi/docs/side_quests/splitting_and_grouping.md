# विभाजन और समूहीकरण

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow डेटा के साथ लचीले ढंग से काम करने के लिए शक्तिशाली उपकरण प्रदान करता है। एक प्रमुख क्षमता है डेटा को विभिन्न स्ट्रीम में विभाजित करना और फिर संबंधित आइटम को वापस एक साथ समूहित करना। यह बायोइन्फॉर्मेटिक्स वर्कफ़्लो में विशेष रूप से मूल्यवान है जहाँ आपको विश्लेषण के लिए परिणामों को संयोजित करने से पहले विभिन्न प्रकार के नमूनों को अलग-अलग प्रोसेस करने की आवश्यकता होती है।

इसे मेल को सॉर्ट करने की तरह समझें: आप गंतव्य के आधार पर पत्रों को अलग करते हैं, प्रत्येक ढेर को अलग तरीके से प्रोसेस करते हैं, फिर एक ही व्यक्ति के पास जाने वाली वस्तुओं को फिर से जोड़ते हैं। Nextflow वैज्ञानिक डेटा के साथ ऐसा करने के लिए विशेष ऑपरेटर का उपयोग करता है। इस दृष्टिकोण को डिस्ट्रिब्यूटेड कंप्यूटिंग और बायोइन्फॉर्मेटिक्स वर्कफ़्लो में **scatter/gather** पैटर्न के रूप में भी जाना जाता है।

Nextflow का channel सिस्टम इस लचीलेपन का केंद्र है। Channel आपके वर्कफ़्लो के विभिन्न भागों को जोड़ते हैं, जिससे डेटा आपके विश्लेषण के माध्यम से प्रवाहित हो सकता है। आप एक ही डेटा स्रोत से कई channel बना सकते हैं, प्रत्येक channel को अलग तरीके से प्रोसेस कर सकते हैं, और फिर आवश्यकता पड़ने पर channel को वापस एक साथ मर्ज कर सकते हैं। यह दृष्टिकोण आपको ऐसे वर्कफ़्लो डिज़ाइन करने देता है जो स्वाभाविक रूप से जटिल बायोइन्फॉर्मेटिक्स विश्लेषणों के शाखाओं और अभिसरण पथों को दर्शाते हैं।

### सीखने के लक्ष्य

इस side quest में, आप Nextflow के channel ऑपरेटर का उपयोग करके डेटा को विभाजित और समूहित करना सीखेंगे।
हम एक CSV फ़ाइल से शुरू करेंगे जिसमें नमूना जानकारी और संबंधित डेटा फ़ाइलें होंगी, फिर इस डेटा को मैनिपुलेट और पुनर्गठित करेंगे।

इस side quest के अंत तक, आप निम्नलिखित तकनीकों का उपयोग करके डेटा स्ट्रीम को प्रभावी ढंग से अलग और संयोजित करने में सक्षम होंगे:

- `splitCsv` का उपयोग करके फ़ाइलों से डेटा पढ़ना
- `filter` और `map` के साथ डेटा को फ़िल्टर और ट्रांसफ़ॉर्म करना
- `join` और `groupTuple` का उपयोग करके संबंधित डेटा को संयोजित करना
- समानांतर प्रोसेसिंग के लिए `combine` के साथ डेटा संयोजन बनाना
- `subMap` और डुप्लीकेशन रणनीतियों का उपयोग करके डेटा संरचना को अनुकूलित करना
- channel संरचनाओं को मैनिपुलेट करने में मदद करने के लिए named closures के साथ पुन: प्रयोज्य फ़ंक्शन बनाना

ये कौशल आपको ऐसे वर्कफ़्लो बनाने में मदद करेंगे जो कई इनपुट फ़ाइलों और विभिन्न प्रकार के डेटा को कुशलता से संभाल सकते हैं, जबकि स्वच्छ, रखरखाव योग्य कोड संरचना बनाए रखते हैं।

### पूर्वापेक्षाएँ

इस side quest को शुरू करने से पहले, आपको चाहिए:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष शुरुआती पाठ्यक्रम पूरा किया हो।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (processes, channels, operators, फ़ाइलों के साथ काम करना, meta data) का उपयोग करने में सहज हों

**वैकल्पिक:** हम पहले [Metadata in workflows](./metadata.md) side quest पूरा करने की सलाह देते हैं।
यह `splitCsv` के साथ CSV फ़ाइलों को पढ़ने और meta maps बनाने की बुनियादी बातों को कवर करता है, जिसका हम यहाँ भारी उपयोग करेंगे।

---

## 0. शुरुआत करें

#### प्रशिक्षण codespace खोलें

यदि आपने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित अनुसार प्रशिक्षण वातावरण खोलना सुनिश्चित करें।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाएं

चलिए उस डायरेक्टरी में चलते हैं जहाँ इस ट्यूटोरियल के लिए फ़ाइलें स्थित हैं।

```bash
cd side-quests/splitting_and_grouping
```

आप VSCode को इस डायरेक्टरी पर फ़ोकस करने के लिए सेट कर सकते हैं:

```bash
code .
```

#### सामग्री की समीक्षा करें

आपको एक मुख्य वर्कफ़्लो फ़ाइल और एक `data` डायरेक्टरी मिलेगी जिसमें `samplesheet.csv` नाम की एक samplesheet है।

```console title="डायरेक्टरी सामग्री"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

Samplesheet में विभिन्न रोगियों के नमूनों के बारे में जानकारी होती है, जिसमें patient ID, sample repeat number, type (normal या tumor), और काल्पनिक डेटा फ़ाइलों के पथ शामिल हैं (जो वास्तव में मौजूद नहीं हैं, लेकिन हम मान लेंगे कि वे हैं)।

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

यह samplesheet तीन रोगियों (A, B, C) से आठ नमूनों की सूची बनाती है।

प्रत्येक रोगी के लिए, हमारे पास ऐसे नमूने हैं जो `tumor` प्रकार के हैं (आमतौर पर ट्यूमर बायोप्सी से उत्पन्न) या `normal` (स्वस्थ ऊतक या रक्त से लिए गए)।
यदि आप कैंसर विश्लेषण से परिचित नहीं हैं, तो बस जान लें कि यह एक प्रायोगिक मॉडल से मेल खाता है जो विपरीत विश्लेषण करने के लिए युग्मित tumor/normal नमूनों का उपयोग करता है।

विशेष रूप से रोगी A के लिए, हमारे पास तकनीकी प्रतिकृतियों (repeats) के दो सेट हैं।

!!! note

    यदि आप इस प्रायोगिक डिज़ाइन से परिचित नहीं हैं तो चिंता न करें, यह इस ट्यूटोरियल को समझने के लिए महत्वपूर्ण नहीं है।

#### असाइनमेंट की समीक्षा करें

आपकी चुनौती एक Nextflow वर्कफ़्लो लिखना है जो:

1. CSV फ़ाइल से नमूना डेटा **पढ़ेगा** और इसे meta maps के साथ संरचित करेगा
2. प्रकार (normal बनाम tumor) के आधार पर नमूनों को विभिन्न channel में **अलग** करेगा
3. patient ID और replicate number द्वारा मिलान किए गए tumor/normal जोड़ों को **जोड़ेगा**
4. समानांतर प्रोसेसिंग के लिए genomic intervals में नमूनों को **वितरित** करेगा
5. डाउनस्ट्रीम विश्लेषण के लिए संबंधित नमूनों को वापस एक साथ **समूहित** करेगा

यह एक सामान्य बायोइन्फॉर्मेटिक्स पैटर्न का प्रतिनिधित्व करता है जहाँ आपको स्वतंत्र प्रोसेसिंग के लिए डेटा को विभाजित करने की आवश्यकता होती है, फिर तुलनात्मक विश्लेषण के लिए संबंधित आइटम को फिर से संयोजित करना होता है।

#### तैयारी चेकलिस्ट

क्या आप गोता लगाने के लिए तैयार हैं?

- [ ] मैं इस पाठ्यक्रम के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता हूँ
- [ ] मेरा codespace चालू है
- [ ] मैंने अपनी कार्य डायरेक्टरी उचित रूप से सेट कर ली है
- [ ] मैं असाइनमेंट को समझता हूँ

यदि आप सभी बॉक्स को चेक कर सकते हैं, तो आप जाने के लिए तैयार हैं।

---

## 1. नमूना डेटा पढ़ें

### 1.1. `splitCsv` के साथ नमूना डेटा पढ़ें और meta maps बनाएं

चलिए `splitCsv` के साथ नमूना डेटा को पढ़कर और इसे meta map पैटर्न में व्यवस्थित करके शुरुआत करते हैं। `main.nf` में, आप देखेंगे कि हमने पहले से ही वर्कफ़्लो शुरू किया है।

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note

    इस ट्यूटोरियल के दौरान, हम सभी channel वेरिएबल के लिए `ch_` prefix का उपयोग करेंगे ताकि स्पष्ट रूप से संकेत दिया जा सके कि वे Nextflow channels हैं।

यदि आपने [Metadata in workflows](./metadata.md) side quest पूरा किया है, तो आप इस पैटर्न को पहचानेंगे। हम CSV को पढ़ने के लिए `splitCsv` का उपयोग करेंगे और तुरंत डेटा को meta map के साथ संरचित करेंगे ताकि मेटाडेटा को फ़ाइल पथों से अलग किया जा सके।

!!! info

    इस प्रशिक्षण में हम दो अलग-अलग अवधारणाओं का सामना करेंगे जिन्हें `map` कहा जाता है:

    - **डेटा संरचना**: Groovy map (अन्य भाषाओं में dictionaries/hashes के बराबर) जो key-value pairs को स्टोर करता है
    - **Channel ऑपरेटर**: `.map()` ऑपरेटर जो channel में आइटम को ट्रांसफ़ॉर्म करता है

    हम संदर्भ में स्पष्ट करेंगे कि हमारा क्या मतलब है, लेकिन Nextflow के साथ काम करते समय यह भेद समझना महत्वपूर्ण है।

`main.nf` में ये परिवर्तन लागू करें:

=== "बाद"

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

यह `splitCsv` ऑपरेशन (headers के साथ CSV पढ़ना) और `map` ऑपरेशन (डेटा को `[meta, file]` tuples के रूप में संरचित करना) को एक चरण में जोड़ता है। उस परिवर्तन को लागू करें और pipeline चलाएं:

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

अब हमारे पास एक channel है जहाँ प्रत्येक आइटम एक `[meta, file]` tuple है - मेटाडेटा फ़ाइल पथों से अलग है। यह संरचना हमें मेटाडेटा फ़ील्ड के आधार पर अपने workload को विभाजित और समूहित करने की अनुमति देती है।

---

## 2. डेटा को फ़िल्टर और ट्रांसफ़ॉर्म करें

### 2.1. `filter` के साथ डेटा को फ़िल्टर करें

हम एक शर्त के आधार पर डेटा को फ़िल्टर करने के लिए [`filter` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#filter) का उपयोग कर सकते हैं। मान लें कि हम केवल normal नमूनों को प्रोसेस करना चाहते हैं। हम `type` फ़ील्ड के आधार पर डेटा को फ़िल्टर करके ऐसा कर सकते हैं। चलिए इसे `view` ऑपरेटर से पहले डालते हैं।

=== "बाद"

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

फ़िल्टर किए गए परिणाम देखने के लिए वर्कफ़्लो को फिर से चलाएं:

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

हमने सफलतापूर्वक डेटा को केवल normal नमूनों को शामिल करने के लिए फ़िल्टर किया है। चलिए संक्षेप में देखें कि यह कैसे काम करता है।

`filter` ऑपरेटर एक closure लेता है जो channel में प्रत्येक तत्व पर लागू होता है। यदि closure `true` लौटाता है, तो तत्व शामिल है; यदि यह `false` लौटाता है, तो तत्व बाहर रखा जाता है।

हमारे मामले में, हम केवल उन नमूनों को रखना चाहते हैं जहाँ `meta.type == 'normal'` है। Closure tuple `meta,file` का उपयोग प्रत्येक नमूने को संदर्भित करने के लिए करता है, `meta.type` के साथ नमूना प्रकार को एक्सेस करता है, और जांचता है कि यह `'normal'` के बराबर है या नहीं।

यह ऊपर प्रस्तुत single closure के साथ पूरा होता है:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. अलग फ़िल्टर किए गए channel बनाएं

वर्तमान में हम CSV से सीधे बनाए गए channel पर फ़िल्टर लागू कर रहे हैं, लेकिन हम इसे एक से अधिक तरीकों से फ़िल्टर करना चाहते हैं, तो चलिए normal नमूनों के लिए एक अलग फ़िल्टर किया गया channel बनाने के लिए तर्क को फिर से लिखते हैं:

=== "बाद"

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

परिणाम देखने के लिए pipeline चलाएं:

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

हमने सफलतापूर्वक डेटा को फ़िल्टर किया है और normal नमूनों के लिए एक अलग channel बनाया है।

चलिए tumor नमूनों के लिए भी एक फ़िल्टर किया गया channel बनाते हैं:

=== "बाद"

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

हमने normal और tumor नमूनों को दो अलग-अलग channel में अलग कर दिया है, और आउटपुट में उन्हें अलग तरह से लेबल करने के लिए `view()` को प्रदान किए गए एक closure का उपयोग किया है: `ch_tumor_samples.view{'Tumor sample: ' + it}`।

### निष्कर्ष

इस खंड में, आपने सीखा:

- **डेटा फ़िल्टर करना**: `filter` के साथ डेटा को कैसे फ़िल्टर करें
- **डेटा विभाजित करना**: किसी शर्त के आधार पर डेटा को विभिन्न channel में कैसे विभाजित करें
- **डेटा देखना**: डेटा प्रिंट करने और विभिन्न channel से आउटपुट को लेबल करने के लिए `view` का उपयोग कैसे करें

हमने अब normal और tumor नमूनों को दो अलग-अलग channel में अलग कर दिया है। अगला, हम `id` फ़ील्ड पर normal और tumor नमूनों को जोड़ेंगे।

---

## 3. पहचानकर्ताओं द्वारा channel जोड़ना

पिछले खंड में, हमने normal और tumor नमूनों को दो अलग-अलग channel में अलग किया। इन्हें उनके प्रकार के आधार पर विशिष्ट processes या workflows का उपयोग करके स्वतंत्र रूप से प्रोसेस किया जा सकता है। लेकिन क्या होता है जब हम एक ही रोगी से normal और tumor नमूनों की तुलना करना चाहते हैं? इस बिंदु पर, हमें उन्हें वापस एक साथ जोड़ने की आवश्यकता है, यह सुनिश्चित करते हुए कि नमूनों को उनके `id` फ़ील्ड के आधार पर मैच किया जाए।

Nextflow में channel को संयोजित करने के कई तरीके शामिल हैं, लेकिन इस मामले में सबसे उपयुक्त ऑपरेटर [`join`](https://www.nextflow.io/docs/latest/operator.html#join) है। यदि आप SQL से परिचित हैं, तो यह `JOIN` ऑपरेशन की तरह कार्य करता है, जहाँ हम join करने के लिए key और join के प्रकार को निर्दिष्ट करते हैं।

### 3.1. patient ID के आधार पर संयोजित करने के लिए `map` और `join` का उपयोग करें

यदि हम [`join`](https://www.nextflow.io/docs/latest/operator.html#join) दस्तावेज़ीकरण की जांच करते हैं, तो हम देख सकते हैं कि डिफ़ॉल्ट रूप से यह प्रत्येक tuple के पहले आइटम के आधार पर दो channel को join करता है।

#### 3.1.1. डेटा संरचना की जांच करें

यदि आपके पास अभी भी console आउटपुट उपलब्ध नहीं है, तो चलिए हमारी डेटा संरचना की जांच करने और देखने के लिए pipeline चलाते हैं कि हमें `id` फ़ील्ड पर join करने के लिए इसे कैसे संशोधित करने की आवश्यकता है।

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

हम देख सकते हैं कि `id` फ़ील्ड प्रत्येक meta map में पहला तत्व है। `join` के काम करने के लिए, हमें प्रत्येक tuple में `id` फ़ील्ड को अलग करना चाहिए। उसके बाद, हम दोनों channel को संयोजित करने के लिए बस `join` ऑपरेटर का उपयोग कर सकते हैं।

#### 3.1.2. `id` फ़ील्ड को अलग करें

`id` फ़ील्ड को अलग करने के लिए, हम [`map` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#map) का उपयोग करके एक नया tuple बना सकते हैं जिसमें `id` फ़ील्ड पहले तत्व के रूप में है।

=== "बाद"

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

यह सूक्ष्म हो सकता है, लेकिन आपको प्रत्येक tuple में पहला तत्व `id` फ़ील्ड होना चाहिए।

#### 3.1.3. दोनों channel को संयोजित करें

अब हम `id` फ़ील्ड के आधार पर दोनों channel को संयोजित करने के लिए `join` ऑपरेटर का उपयोग कर सकते हैं।

एक बार फिर, हम joined आउटपुट प्रिंट करने के लिए `view` का उपयोग करेंगे।

=== "बाद"

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

यह बताना थोड़ा मुश्किल है क्योंकि यह इतना चौड़ा है, लेकिन आपको देखने में सक्षम होना चाहिए कि नमूनों को `id` फ़ील्ड द्वारा joined किया गया है। प्रत्येक tuple में अब यह प्रारूप है:

- `id`: नमूना ID
- `normal_meta_map`: normal नमूना मेटा डेटा जिसमें type, replicate और bam फ़ाइल का पथ शामिल है
- `normal_sample_file`: normal नमूना फ़ाइल
- `tumor_meta_map`: tumor नमूना मेटा डेटा जिसमें type, replicate और bam फ़ाइल का पथ शामिल है
- `tumor_sample`: tumor नमूना जिसमें type, replicate और bam फ़ाइल का पथ शामिल है

!!! warning

    `join` ऑपरेटर किसी भी बेमेल tuple को त्याग देगा। इस उदाहरण में, हमने यह सुनिश्चित किया कि सभी नमूने tumor और normal के लिए मेल खाते थे, लेकिन यदि यह सच नहीं है तो आपको बेमेल tuple रखने के लिए पैरामीटर `remainder: true` का उपयोग करना होगा। अधिक विवरण के लिए [दस्तावेज़ीकरण](https://www.nextflow.io/docs/latest/operator.html#join) देखें।

तो अब आप जानते हैं कि tuple में किसी फ़ील्ड को अलग करने के लिए `map` का उपयोग कैसे करें, और पहले फ़ील्ड के आधार पर tuple को संयोजित करने के लिए `join` का उपयोग कैसे करें।
इस ज्ञान के साथ, हम साझा फ़ील्ड के आधार पर channel को सफलतापूर्वक संयोजित कर सकते हैं।

अगला, हम उस स्थिति पर विचार करेंगे जहाँ आप कई फ़ील्ड पर join करना चाहते हैं।

### 3.2. कई फ़ील्ड पर join करें

हमारे पास sampleA के लिए 2 replicates हैं, लेकिन sampleB और sampleC के लिए केवल 1। इस मामले में हम `id` फ़ील्ड का उपयोग करके उन्हें प्रभावी ढंग से join करने में सक्षम थे, लेकिन क्या होगा यदि वे sync से बाहर हों? हम विभिन्न replicates से normal और tumor नमूनों को मिला सकते हैं!

इससे बचने के लिए, हम कई फ़ील्ड पर join कर सकते हैं। वास्तव में इसे प्राप्त करने के कई तरीके हैं लेकिन हम एक नई joining key बनाने पर ध्यान केंद्रित करने जा रहे हैं जिसमें नमूना `id` और `replicate` number दोनों शामिल हैं।

चलिए एक नई joining key बनाने से शुरू करते हैं। हम पहले की तरह [`map` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#map) का उपयोग करके ऐसा कर सकते हैं ताकि एक नया tuple बनाया जा सके जिसमें `id` और `repeat` फ़ील्ड पहले तत्व के रूप में हों।

=== "बाद"

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

अब हमें join होता दिखना चाहिए लेकिन `id` और `repeat` दोनों फ़ील्ड का उपयोग करते हुए। वर्कफ़्लो चलाएं:

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

ध्यान दें कि हमारे पास प्रत्येक joined परिणाम के पहले तत्व के रूप में दो तत्वों (`id` और `repeat` फ़ील्ड) का एक tuple है। यह दर्शाता है कि जटिल आइटम को joining key के रूप में कैसे उपयोग किया जा सकता है, जो एक ही स्थितियों से नमूनों के बीच काफी जटिल मिलान को सक्षम बनाता है।

यदि आप विभिन्न keys पर join करने के अधिक तरीकों का पता लगाना चाहते हैं, तो अतिरिक्त विकल्पों और उदाहरणों के लिए [join ऑपरेटर दस्तावेज़ीकरण](https://www.nextflow.io/docs/latest/operator.html#join) देखें।

### 3.3. नई joining key बनाने के लिए `subMap` का उपयोग करें

पिछला दृष्टिकोण हमारी joining key से फ़ील्ड नामों को खो देता है - `id` और `repeat` फ़ील्ड केवल मूल्यों की एक list बन जाते हैं। बाद में पहुँच के लिए फ़ील्ड नामों को बनाए रखने के लिए, हम [`subMap` method](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>) का उपयोग कर सकते हैं।

`subMap` method एक map से केवल निर्दिष्ट key-value pairs को extract करता है। यहाँ हम अपनी joining key बनाने के लिए केवल `id` और `repeat` फ़ील्ड को extract करेंगे।

=== "बाद"

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

अब हमारे पास एक नई joining key है जो न केवल `id` और `repeat` फ़ील्ड शामिल करती है, बल्कि फ़ील्ड नामों को भी बनाए रखती है ताकि हम उन्हें बाद में नाम से एक्सेस कर सकें, जैसे `meta.id` और `meta.repeat`।

### 3.4. map में named closure का उपयोग करें

डुप्लीकेशन से बचने और त्रुटियों को कम करने के लिए, हम एक named closure का उपयोग कर सकते हैं। एक named closure हमें एक पुन: प्रयोज्य फ़ंक्शन बनाने की अनुमति देता है जिसे हम कई स्थानों पर call कर सकते हैं।

ऐसा करने के लिए, पहले हम closure को एक नए वेरिएबल के रूप में परिभाषित करते हैं:

=== "बाद"

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

हमने map transformation को एक named वेरिएबल के रूप में परिभाषित किया है जिसे हम पुन: उपयोग कर सकते हैं।

ध्यान दें कि हम फ़ाइल पथ को `file()` का उपयोग करके एक Path object में भी convert करते हैं ताकि इस channel को प्राप्त करने वाली कोई भी process फ़ाइल को सही ढंग से handle कर सके (अधिक जानकारी के लिए [Working with files](./working_with_files.md) देखें)।

चलिए अपने वर्कफ़्लो में closure को implement करते हैं:

=== "बाद"

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

!!! note

    `map` ऑपरेटर closure को एक argument के रूप में पास करने के लिए `{ }` के बजाय `( )` का उपयोग करने के लिए स्विच हो गया है। ऐसा इसलिए है क्योंकि `map` ऑपरेटर एक closure को argument के रूप में अपेक्षा करता है और `{ }` का उपयोग anonymous closure को परिभाषित करने के लिए किया जाता है। named closure को call करते समय, `( )` syntax का उपयोग करें।

यह जांचने के लिए कि सबकुछ अभी भी काम कर रहा है, वर्कफ़्लो को एक बार और चलाएं:

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

named closure का उपयोग करने से हम कई स्थानों पर एक ही transformation को पुन: उपयोग कर सकते हैं, जिससे त्रुटियों का जोखिम कम होता है और कोड अधिक पठनीय और रखरखाव योग्य बनता है।

### 3.5. डेटा के डुप्लीकेशन को कम करें

हमारे वर्कफ़्लो में बहुत सारा डुप्लिकेट डेटा है। joined नमूनों में प्रत्येक आइटम `id` और `repeat` फ़ील्ड को दोहराता है। चूँकि यह जानकारी पहले से grouping key में उपलब्ध है, इसलिए हम इस अतिरेक से बच सकते हैं। एक reminder के रूप में, हमारी वर्तमान डेटा संरचना इस तरह दिखती है:

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

चूँकि `id` और `repeat` फ़ील्ड grouping key में उपलब्ध हैं, चलिए डुप्लीकेशन से बचने के लिए प्रत्येक channel आइटम के बाकी हिस्सों से उन्हें हटाते हैं। हम केवल `type` फ़ील्ड के साथ एक नया map बनाने के लिए `subMap` method का उपयोग कर सकते हैं। यह दृष्टिकोण हमें हमारे डेटा संरचना में अतिरेक को समाप्त करते हुए सभी आवश्यक जानकारी बनाए रखने की अनुमति देता है।

=== "बाद"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

अब closure एक tuple लौटाता है जहाँ पहले तत्व में `id` और `repeat` फ़ील्ड हैं, और दूसरे तत्व में केवल `type` फ़ील्ड है। यह grouping key
