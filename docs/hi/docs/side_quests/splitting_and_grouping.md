# डेटा को विभाजित और समूहित करना

Nextflow डेटा के साथ लचीले तरीके से काम करने के लिए शक्तिशाली टूल प्रदान करता है। एक प्रमुख क्षमता है डेटा को विभिन्न स्ट्रीम में विभाजित करना और फिर संबंधित आइटम को वापस एक साथ समूहित करना। यह बायोइनफॉर्मेटिक्स वर्कफ़्लो में विशेष रूप से मूल्यवान है जहाँ तुम्हें विभिन्न प्रकार के नमूनों को अलग-अलग प्रोसेस करने की आवश्यकता होती है और फिर विश्लेषण के लिए परिणामों को संयोजित करना होता है।

इसे डाक छाँटने की तरह समझो: तुम गंतव्य के अनुसार पत्रों को अलग करते हो, प्रत्येक ढेर को अलग तरीके से प्रोसेस करते हो, फिर एक ही व्यक्ति के पास जाने वाली वस्तुओं को फिर से जोड़ते हो। Nextflow वैज्ञानिक डेटा के साथ यह करने के लिए विशेष ऑपरेटर का उपयोग करता है। इस दृष्टिकोण को वितरित कंप्यूटिंग और बायोइनफॉर्मेटिक्स वर्कफ़्लो में आमतौर पर **scatter/gather** पैटर्न के रूप में भी जाना जाता है।

Nextflow का channel सिस्टम इस लचीलेपन के केंद्र में है। Channel तुम्हारे वर्कफ़्लो के विभिन्न भागों को जोड़ते हैं, जिससे डेटा तुम्हारे विश्लेषण के माध्यम से प्रवाहित हो सकता है। तुम एक ही डेटा स्रोत से कई channel बना सकते हो, प्रत्येक channel को अलग तरीके से प्रोसेस कर सकते हो, और फिर आवश्यकता पड़ने पर channel को वापस एक साथ मर्ज कर सकते हो। यह दृष्टिकोण तुम्हें ऐसे वर्कफ़्लो डिज़ाइन करने देता है जो स्वाभाविक रूप से जटिल बायोइनफॉर्मेटिक्स विश्लेषण के शाखाओं और अभिसरण पथों को दर्शाते हैं।

### सीखने के लक्ष्य

इस साइड क्वेस्ट में, तुम Nextflow के channel ऑपरेटर का उपयोग करके डेटा को विभाजित और समूहित करना सीखोगे।
हम एक CSV फ़ाइल से शुरू करेंगे जिसमें नमूना जानकारी और संबंधित डेटा फ़ाइलें होंगी, फिर इस डेटा को मैनिपुलेट और पुनर्गठित करेंगे।

इस साइड क्वेस्ट के अंत तक, तुम निम्नलिखित तकनीकों का उपयोग करके डेटा स्ट्रीम को प्रभावी ढंग से अलग और संयोजित करने में सक्षम होगे:

- `splitCsv` का उपयोग करके फ़ाइलों से डेटा पढ़ना
- `filter` और `map` के साथ डेटा को फ़िल्टर और ट्रांसफ़ॉर्म करना
- `join` और `groupTuple` का उपयोग करके संबंधित डेटा को संयोजित करना
- समानांतर प्रोसेसिंग के लिए `combine` के साथ डेटा संयोजन बनाना
- `subMap` और डिडुप्लीकेशन रणनीतियों का उपयोग करके डेटा संरचना को अनुकूलित करना
- channel संरचनाओं को मैनिपुलेट करने में मदद के लिए नामित closures के साथ पुन: प्रयोज्य फ़ंक्शन बनाना

ये कौशल तुम्हें ऐसे वर्कफ़्लो बनाने में मदद करेंगे जो कई इनपुट फ़ाइलों और विभिन्न प्रकार के डेटा को कुशलता से संभाल सकते हैं, साथ ही स्वच्छ, रखरखाव योग्य कोड संरचना बनाए रख सकते हैं।

### पूर्वापेक्षाएँ

इस साइड क्वेस्ट को शुरू करने से पहले, तुम्हें:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष शुरुआती पाठ्यक्रम पूरा कर लेना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (processes, channels, operators, फ़ाइलों के साथ काम करना, meta data) का उपयोग करने में सहज होना चाहिए

**वैकल्पिक:** हम पहले [Metadata in workflows](./metadata.md) साइड क्वेस्ट पूरा करने की सलाह देते हैं।
यह `splitCsv` के साथ CSV फ़ाइलें पढ़ने और meta maps बनाने की बुनियादी बातों को कवर करता है, जिसका हम यहाँ भारी उपयोग करेंगे।

---

## 0. शुरू करना

#### ट्रेनिंग codespace खोलें

यदि तुमने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित अनुसार ट्रेनिंग एनवायरनमेंट खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाओ

चलो उस डायरेक्टरी में चलते हैं जहाँ इस ट्यूटोरियल के लिए फ़ाइलें स्थित हैं।

```bash
cd side-quests/splitting_and_grouping
```

तुम VSCode को इस डायरेक्टरी पर फ़ोकस करने के लिए सेट कर सकते हो:

```bash
code .
```

#### सामग्री की समीक्षा करो

तुम्हें एक मुख्य वर्कफ़्लो फ़ाइल और एक `data` डायरेक्टरी मिलेगी जिसमें `samplesheet.csv` नाम की एक samplesheet है।

```console title="Directory contents"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

Samplesheet में विभिन्न रोगियों के नमूनों के बारे में जानकारी है, जिसमें रोगी ID, नमूना रिपीट नंबर, प्रकार (normal या tumor), और काल्पनिक डेटा फ़ाइलों के पथ शामिल हैं (जो वास्तव में मौजूद नहीं हैं, लेकिन हम मान लेंगे कि वे हैं)।

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

यह samplesheet तीन रोगियों (A, B, C) के आठ नमूनों को सूचीबद्ध करती है।

प्रत्येक रोगी के लिए, हमारे पास `tumor` प्रकार के नमूने हैं (आमतौर पर ट्यूमर बायोप्सी से उत्पन्न) या `normal` (स्वस्थ ऊतक या रक्त से लिए गए)।
यदि तुम कैंसर विश्लेषण से परिचित नहीं हो, तो बस यह जान लो कि यह एक प्रयोगात्मक मॉडल से मेल खाता है जो तुलनात्मक विश्लेषण करने के लिए युग्मित tumor/normal नमूनों का उपयोग करता है।

रोगी A के लिए विशेष रूप से, हमारे पास तकनीकी प्रतिकृतियों (repeats) के दो सेट हैं।

!!! note

    चिंता मत करो यदि तुम इस प्रयोगात्मक डिज़ाइन से परिचित नहीं हो, यह इस ट्यूटोरियल को समझने के लिए महत्वपूर्ण नहीं है।

#### असाइनमेंट की समीक्षा करो

तुम्हारी चुनौती एक Nextflow वर्कफ़्लो लिखना है जो:

1. CSV फ़ाइल से नमूना डेटा **पढ़ेगा** और इसे meta maps के साथ संरचित करेगा
2. प्रकार (normal बनाम tumor) के आधार पर नमूनों को विभिन्न channel में **अलग** करेगा
3. रोगी ID और रिपीट नंबर द्वारा मिलान किए गए tumor/normal जोड़ों को **जोड़ेगा**
4. समानांतर प्रोसेसिंग के लिए जीनोमिक अंतरालों में नमूनों को **वितरित** करेगा
5. डाउनस्ट्रीम विश्लेषण के लिए संबंधित नमूनों को वापस एक साथ **समूहित** करेगा

यह एक सामान्य बायोइनफॉर्मेटिक्स पैटर्न का प्रतिनिधित्व करता है जहाँ तुम्हें स्वतंत्र प्रोसेसिंग के लिए डेटा को विभाजित करने की आवश्यकता होती है, फिर तुलनात्मक विश्लेषण के लिए संबंधित आइटम को फिर से संयोजित करना होता है।

#### तैयारी चेकलिस्ट

लगता है कि तुम शुरू करने के लिए तैयार हो?

- [ ] मैं इस पाठ्यक्रम के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता हूँ
- [ ] मेरा codespace चालू है और चल रहा है
- [ ] मैंने अपनी वर्किंग डायरेक्टरी उचित रूप से सेट कर ली है
- [ ] मैं असाइनमेंट को समझता हूँ

यदि तुम सभी बॉक्स चेक कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. नमूना डेटा पढ़ना

### 1.1. `splitCsv` के साथ नमूना डेटा पढ़ें और meta maps बनाएँ

चलो `splitCsv` के साथ नमूना डेटा पढ़कर और इसे meta map पैटर्न में व्यवस्थित करके शुरू करते हैं। `main.nf` में, तुम देखोगे कि हमने पहले से ही वर्कफ़्लो शुरू कर दिया है।

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note

    इस ट्यूटोरियल के दौरान, हम सभी channel वेरिएबल के लिए `ch_` प्रीफ़िक्स का उपयोग करेंगे ताकि स्पष्ट रूप से संकेत मिल सके कि वे Nextflow channels हैं।

यदि तुमने [Metadata in workflows](./metadata.md) साइड क्वेस्ट पूरा किया है, तो तुम इस पैटर्न को पहचानोगे। हम CSV पढ़ने के लिए `splitCsv` का उपयोग करेंगे और तुरंत डेटा को meta map के साथ संरचित करेंगे ताकि मेटाडेटा को फ़ाइल पथों से अलग किया जा सके।

!!! info

    हम इस ट्रेनिंग में `map` नामक दो अलग-अलग अवधारणाओं का सामना करेंगे:

    - **डेटा संरचना**: Groovy map (अन्य भाषाओं में dictionaries/hashes के समकक्ष) जो key-value जोड़े संग्रहीत करता है
    - **Channel ऑपरेटर**: `.map()` ऑपरेटर जो channel में आइटम को ट्रांसफ़ॉर्म करता है

    हम संदर्भ में स्पष्ट करेंगे कि हमारा क्या मतलब है, लेकिन Nextflow के साथ काम करते समय यह अंतर समझना महत्वपूर्ण है।

`main.nf` में ये परिवर्तन लागू करो:

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

यह `splitCsv` ऑपरेशन (हेडर के साथ CSV पढ़ना) और `map` ऑपरेशन (डेटा को `[meta, file]` tuples के रूप में संरचित करना) को एक चरण में संयोजित करता है। वह परिवर्तन लागू करो और पाइपलाइन चलाओ:

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

अब हमारे पास एक channel है जहाँ प्रत्येक आइटम एक `[meta, file]` tuple है - मेटाडेटा फ़ाइल पथों से अलग है। यह संरचना हमें मेटाडेटा फ़ील्ड के आधार पर अपने वर्कलोड को विभाजित और समूहित करने की अनुमति देती है।

---

## 2. डेटा को फ़िल्टर और ट्रांसफ़ॉर्म करना

### 2.1. `filter` के साथ डेटा फ़िल्टर करें

हम किसी शर्त के आधार पर डेटा को फ़िल्टर करने के लिए [`filter` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#filter) का उपयोग कर सकते हैं। मान लो हम केवल normal नमूनों को प्रोसेस करना चाहते हैं। हम `type` फ़ील्ड के आधार पर डेटा को फ़िल्टर करके ऐसा कर सकते हैं। चलो इसे `view` ऑपरेटर से पहले डालते हैं।

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

फ़िल्टर किए गए परिणाम को देखने के लिए वर्कफ़्लो फिर से चलाओ:

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

हमने सफलतापूर्वक डेटा को केवल normal नमूनों को शामिल करने के लिए फ़िल्टर किया है। चलो संक्षेप में देखें कि यह कैसे काम करता है।

`filter` ऑपरेटर एक closure लेता है जो channel में प्रत्येक तत्व पर लागू होता है। यदि closure `true` लौटाता है, तो तत्व शामिल किया जाता है; यदि यह `false` लौटाता है, तो तत्व बाहर रखा जाता है।

हमारे मामले में, हम केवल उन नमूनों को रखना चाहते हैं जहाँ `meta.type == 'normal'`। Closure प्रत्येक नमूने को संदर्भित करने के लिए tuple `meta,file` का उपयोग करता है, `meta.type` के साथ नमूना प्रकार तक पहुँचता है, और जाँचता है कि यह `'normal'` के बराबर है या नहीं।

यह हमारे द्वारा ऊपर पेश किए गए एकल closure के साथ पूरा किया जाता है:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. अलग फ़िल्टर किए गए channel बनाएँ

वर्तमान में हम CSV से सीधे बनाए गए channel पर फ़िल्टर लागू कर रहे हैं, लेकिन हम इसे एक से अधिक तरीकों से फ़िल्टर करना चाहते हैं, इसलिए चलो normal नमूनों के लिए एक अलग फ़िल्टर किया गया channel बनाने के लिए लॉजिक को फिर से लिखते हैं:

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

परिणाम देखने के लिए पाइपलाइन चलाओ:

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

चलो tumor नमूनों के लिए भी एक फ़िल्टर किया गया channel बनाते हैं:

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

हमने normal और tumor नमूनों को दो अलग-अलग channel में अलग कर दिया है, और आउटपुट में उन्हें अलग तरीके से लेबल करने के लिए `view()` को दिए गए closure का उपयोग किया है: `ch_tumor_samples.view{'Tumor sample: ' + it}`।

### सारांश

इस सेक्शन में, तुमने सीखा:

- **डेटा फ़िल्टर करना**: `filter` के साथ डेटा को कैसे फ़िल्टर करें
- **डेटा विभाजित करना**: किसी शर्त के आधार पर डेटा को विभिन्न channel में कैसे विभाजित करें
- **डेटा देखना**: डेटा प्रिंट करने और विभिन्न channel से आउटपुट को लेबल करने के लिए `view` का उपयोग कैसे करें

हमने अब normal और tumor नमूनों को दो अलग-अलग channel में अलग कर दिया है। अगला, हम `id` फ़ील्ड पर normal और tumor नमूनों को जोड़ेंगे।

---

## 3. पहचानकर्ताओं द्वारा channel जोड़ना

पिछले सेक्शन में, हमने normal और tumor नमूनों को दो अलग-अलग channel में अलग किया। इन्हें उनके प्रकार के आधार पर विशिष्ट processes या वर्कफ़्लो का उपयोग करके स्वतंत्र रूप से प्रोसेस किया जा सकता है। लेकिन क्या होता है जब हम एक ही रोगी के normal और tumor नमूनों की तुलना करना चाहते हैं? इस बिंदु पर, हमें उन्हें वापस एक साथ जोड़ने की आवश्यकता है और यह सुनिश्चित करना होगा कि नमूनों को उनके `id` फ़ील्ड के आधार पर मिलान किया जाए।

Nextflow में channel को संयोजित करने के कई तरीके शामिल हैं, लेकिन इस मामले में सबसे उपयुक्त ऑपरेटर [`join`](https://www.nextflow.io/docs/latest/operator.html#join) है। यदि तुम SQL से परिचित हो, तो यह `JOIN` ऑपरेशन की तरह काम करता है, जहाँ हम join करने के लिए key और join के प्रकार को निर्दिष्ट करते हैं।

### 3.1. रोगी ID के आधार पर संयोजित करने के लिए `map` और `join` का उपयोग करें

#### 3.1.1. डेटा संरचना की जाँच करें

यदि तुम्हारे पास अभी भी कंसोल आउटपुट उपलब्ध नहीं है, तो चलो हमारी डेटा संरचना की जाँच करने के लिए पाइपलाइन चलाते हैं और देखते हैं कि हमें `id` फ़ील्ड पर join करने के लिए इसे कैसे संशोधित करने की आवश्यकता है।

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

`id` फ़ील्ड को अलग करने के लिए, हम [`map` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#map) का उपयोग करके `id` फ़ील्ड को पहले तत्व के रूप में एक नया tuple बना सकते हैं।

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

यह सूक्ष्म हो सकता है, लेकिन तुम्हें प्रत्येक tuple में पहला तत्व `id` फ़ील्ड देखने में सक्षम होना चाहिए।

#### 3.1.3. दोनों channel को संयोजित करें

अब हम `id` फ़ील्ड के आधार पर दोनों channel को संयोजित करने के लिए `join` ऑपरेटर का उपयोग कर सकते हैं।

एक बार फिर, हम जुड़े हुए आउटपुट को प्रिंट करने के लिए `view` का उपयोग करेंगे।

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

यह बताना थोड़ा मुश्किल है क्योंकि यह इतना चौड़ा है, लेकिन तुम्हें देखने में सक्षम होना चाहिए कि नमूनों को `id` फ़ील्ड द्वारा जोड़ा गया है। प्रत्येक tuple में अब यह प्रारूप है:

- `id`: नमूना ID
- `normal_meta_map`: normal नमूना मेटा डेटा जिसमें type, replicate और bam फ़ाइल का पथ शामिल है
- `normal_sample_file`: normal नमूना फ़ाइल
- `tumor_meta_map`: tumor नमूना मेटा डेटा जिसमें type, replicate और bam फ़ाइल का पथ शामिल है
- `tumor_sample`: tumor नमूना जिसमें type, replicate और bam फ़ाइल का पथ शामिल है

!!! warning

    `join` ऑपरेटर किसी भी अन-मैच किए गए tuples को त्याग देगा। इस उदाहरण में, हमने सुनिश्चित किया कि सभी नमूने tumor और normal के लिए मिलान किए गए थे लेकिन यदि यह सच नहीं है तो तुम्हें अन-मैच किए गए tuples को रखने के लिए पैरामीटर `remainder: true` का उपयोग करना होगा। अधिक विवरण के लिए [डॉक्यूमेंटेशन](https://www.nextflow.io/docs/latest/operator.html#join) देखो।

तो अब तुम जानते हो कि tuple में एक फ़ील्ड को अलग करने के लिए `map` का उपयोग कैसे करें, और पहले फ़ील्ड के आधार पर tuples को संयोजित करने के लिए `join` का उपयोग कैसे करें।
इस ज्ञान के साथ, हम एक साझा फ़ील्ड के आधार पर channel को सफलतापूर्वक संयोजित कर सकते हैं।

अगला, हम उस स्थिति पर विचार करेंगे जहाँ तुम कई फ़ील्ड पर join करना चाहते हो।

### 3.2. कई फ़ील्ड पर join करें

हमारे पास sampleA के लिए 2 replicates हैं, लेकिन sampleB और sampleC के लिए केवल 1। इस मामले में हम `id` फ़ील्ड का उपयोग करके उन्हें प्रभावी ढंग से join करने में सक्षम थे, लेकिन क्या होगा यदि वे सिंक से बाहर हों? हम विभिन्न replicates से normal और tumor नमूनों को मिला सकते हैं!

इससे बचने के लिए, हम कई फ़ील्ड पर join कर सकते हैं। वास्तव में इसे प्राप्त करने के कई तरीके हैं लेकिन हम एक नई joining key बनाने पर ध्यान केंद्रित करने जा रहे हैं जिसमें नमूना `id` और `replicate` नंबर दोनों शामिल हैं।

चलो एक नई joining key बनाकर शुरू करते हैं। हम पहले की तरह [`map` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#map) का उपयोग करके `id` और `repeat` फ़ील्ड को पहले तत्व के रूप में एक नया tuple बना सकते हैं।

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

अब हमें देखना चाहिए कि join हो रहा है लेकिन `id` और `repeat` दोनों फ़ील्ड का उपयोग कर रहा है। वर्कफ़्लो चलाओ:

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

ध्यान दो कि हमारे पास प्रत्येक जुड़े हुए परिणाम के पहले तत्व के रूप में दो तत्वों (`id` और `repeat` फ़ील्ड) का एक tuple है। यह दर्शाता है कि जटिल आइटम को joining key के रूप में कैसे उपयोग किया जा सकता है, जो एक ही स्थिति से नमूनों के बीच काफी जटिल मिलान को सक्षम करता है।

यदि तुम विभिन्न keys पर join करने के अधिक तरीकों का पता लगाना चाहते हो, तो अतिरिक्त विकल्पों और उदाहरणों के लिए [join ऑपरेटर डॉक्यूमेंटेशन](https://www.nextflow.io/docs/latest/operator.html#join) देखो।

### 3.3. नई joining key बनाने के लिए `subMap` का उपयोग करें

पिछला दृष्टिकोण हमारी joining key से फ़ील्ड नाम खो देता है - `id` और `repeat` फ़ील्ड केवल मानों की एक सूची बन जाते हैं। बाद में एक्सेस के लिए फ़ील्ड नाम बनाए रखने के लिए, हम [`subMap` मेथड](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>) का उपयोग कर सकते हैं।

`subMap` मेथड एक map से केवल निर्दिष्ट key-value जोड़े निकालता है। यहाँ हम हमारी joining key बनाने के लिए केवल `id` और `repeat` फ़ील्ड निकालेंगे।

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

अब हमारे पास एक नई joining key है जो न केवल `id` और `repeat` फ़ील्ड शामिल करती है बल्कि फ़ील्ड नाम भी बनाए रखती है ताकि हम उन्हें बाद में नाम से एक्सेस कर सकें, जैसे `meta.id` और `meta.repeat`।

### 3.4. map में नामित closure का उपयोग करें

डुप्लीकेशन से बचने और त्रुटियों को कम करने के लिए, हम एक नामित closure का उपयोग कर सकते हैं। एक नामित closure हमें एक पुन: प्रयोज्य फ़ंक्शन बनाने की अनुमति देता है जिसे हम कई स्थानों पर कॉल कर सकते हैं।

ऐसा करने के लिए, पहले हम closure को एक नए वेरिएबल के रूप में परिभाषित करते हैं:

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

हमने map ट्रांसफ़ॉर्मेशन को एक नामित वेरिएबल के रूप में परिभाषित किया है जिसे हम पुन: उपयोग कर सकते हैं।

ध्यान दो कि हम फ़ाइल पथ को `file()` का उपयोग करके Path ऑब्जेक्ट में भी परिवर्तित करते हैं ताकि इस channel को प्राप्त करने वाली कोई भी process फ़ाइल को सही ढंग से संभाल सके (अधिक जानकारी के लिए [Working with files](./working_with_files.md) देखें)।

चलो अपने वर्कफ़्लो में closure को लागू करते हैं:

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

!!! note

    `map` ऑपरेटर `{ }` का उपयोग करने से `( )` का उपयोग करने में बदल गया है ताकि closure को एक आर्गुमेंट के रूप में पास किया जा सके। ऐसा इसलिए है क्योंकि `map` ऑपरेटर एक closure को एक आर्गुमेंट के रूप में अपेक्षा करता है और `{ }` का उपयोग एक anonymous closure को परिभाषित करने के लिए किया जाता है। नामित closure को कॉल करते समय, `( )` सिंटैक्स का उपयोग करो।

सब कुछ अभी भी काम कर रहा है यह जाँचने के लिए वर्कफ़्लो एक बार फिर चलाओ:

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

नामित closure का उपयोग करने से हमें कई स्थानों पर एक ही ट्रांसफ़ॉर्मेशन का पुन: उपयोग करने की अनुमति मिलती है, त्रुटियों के जोखिम को कम करता है और कोड को अधिक पठनीय और रखरखाव योग्य बनाता है।

### 3.5. डेटा के डुप्लीकेशन को कम करें

हमारे वर्कफ़्लो में बहुत सारा डुप्लीकेट डेटा है। जुड़े हुए नमूनों में प्रत्येक आइटम `id` और `repeat` फ़ील्ड को दोहराता है। चूँकि यह जानकारी पहले से ही grouping key में उपलब्ध है, हम इस अतिरेक से बच सकते हैं। एक अनुस्मारक के रूप में, हमारी वर्तमान डेटा संरचना इस तरह दिखती है:

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

चूँकि `id` और `repeat` फ़ील्ड grouping key में उपलब्ध हैं, चलो डुप्लीकेशन से बचने के लिए उन्हें प्रत्येक channel आइटम के बाकी हिस्सों से हटा देते हैं। हम केवल `type` फ़ील्ड के साथ एक नया map बनाने के लिए `subMap` मेथड का उपयोग करके ऐसा कर सकते हैं। यह दृष्टिकोण हमें हमारे डेटा संरचना में अतिरेक को समाप्त करते हुए सभी आवश्यक जानकारी बनाए रखने की अनुमति देता है।

=== "बाद में"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

अब closure एक tuple लौटाता है जहाँ पहला तत्व `id` और `repeat` फ़ील्ड शामिल करता है, और दूसरा तत्व केवल `type` फ़ील्ड शामिल करता है। यह grouping key में `id` और `repeat` जानकारी को एक बार संग्रहीत करके अतिरेक को समाप्त करता है, जबकि सभी आवश्यक जानकारी बनाए रखता है।

यह कैसा दिखता है यह देखने के लिए वर्कफ़्लो चलाओ:

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

हम देख सकते हैं कि हम grouping key में केवल एक बार `id` और `repeat` फ़ील्ड बताते हैं और हमारे पास नमूना डेटा में `type` फ़ील्ड है। हमने कोई जानकारी नहीं खोई है लेकिन हम अपने channel सामग्री को अधिक संक्षिप्त बनाने में कामयाब रहे।

### 3.6. अनावश्यक जानकारी हटाएँ

हमने ऊपर डुप्लीकेट जानकारी हटा दी, लेकिन हमारे channel में अभी भी कुछ अन्य अनावश्यक जानकारी है।

शुरुआत में, हमने `filter` का उपयोग करके normal और tumor नमूनों को अलग किया, फिर उन्हें `id` और `repeat` keys के आधार पर जोड़ा। `join` ऑपरेटर उस क्रम को संरक्षित करता है जिसमें tuples को मर्ज किया जाता है, इसलिए हमारे मामले में, बाईं ओर normal नमूनों और दाईं ओर tumor नमूनों के साथ, परिणामी channel इस संरचना को बनाए रखता है: `id, <normal elements>, <tumor elements>`।

चूँकि हम अपने channel में प्रत्येक तत्व की स्थिति जानते हैं, हम `[type:normal]` और `[type:tumor]` मेटाडेटा को छोड़कर संरचना को और सरल बना सकते हैं।

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

इस सेक्शन में, तुमने सीखा:

- **Tuples को मैनिपुलेट करना**: tuple में एक फ़ील्ड को अलग करने के लिए `map` का उपयोग कैसे करें
- **Tuples को जोड़ना**: पहले फ़ील्ड के आधार पर tuples को संयोजित करने के लिए `join` का उपयोग कैसे करें
- **Joining Keys बनाना**: नई joining key बनाने के लिए `subMap` का उपयोग कैसे करें
- **नामित Closures**: map में नामित closure का उपयोग कैसे करें
- **कई फ़ील्ड पर Joining**: अधिक सटीक मिलान के लिए कई फ़ील्ड पर कैसे join करें
- **डेटा संरचना अनुकूलन**: अनावश्यक जानकारी हटाकर channel संरचना को कैसे सुव्यवस्थित करें

अब तुम्हारे पास एक वर्कफ़्लो है जो एक samplesheet को विभाजित कर सकता है, normal और tumor नमूनों को फ़िल्टर कर सकता है, उन्हें नमूना ID और रिपीट नंबर द्वारा एक साथ जोड़ सकता है, फिर परिणाम प्रिंट कर सकता है।

यह बायोइनफॉर्मेटिक्स वर्कफ़्लो में एक सामान्य पैटर्न है जहाँ तुम्हें स्वतंत्र रूप से प्रोसेसिंग के बाद नमूनों या अन्य प्रकार के डेटा को मिलाने की आवश्यकता होती है, इसलिए यह एक उपयोगी कौशल है। अगला, हम एक नमूने को कई बार दोहराने पर नज़र डालेंगे।

## 4. नमूनों को अंतरालों में फैलाना

बायोइनफॉर्मेटिक्स वर्कफ़्लो में एक प्रमुख पैटर्न जीनोमिक क्षेत्रों में विश्लेषण वितरित करना है। उदाहरण के लिए, variant calling को जीनोम को अंतरालों (जैसे chromosomes या छोटे क्षेत्रों) में विभाजित करके समानांतर किया जा सकता है। यह समानांतरण रणनीति कई cores या nodes में कम्प्यूटेशनल लोड वितरित करके पाइपलाइन दक्षता में काफी सुधार करती है, कुल निष्पादन समय को कम करती है।

निम्नलिखित सेक्शन में, हम प्रदर्शित करेंगे कि हमारे नमूना डेटा को कई जीनोमिक अंतरालों में कैसे वितरित किया जाए। हम प्रत्येक नमूने को हर अंतराल के साथ जोड़ेंगे, जिससे विभिन्न जीनोमिक क्षेत्रों की समानांतर प्रोसेसिंग की अनुमति मिलेगी। यह हमारे डेटासेट के आकार को अंतरालों की संख्या से गुणा करेगा, कई स्वतंत्र विश्लेषण इकाइयाँ बनाएगा जिन्हें बाद में वापस एक साथ लाया जा सकता है।

### 4.1. `combine` का उपयोग करके अंतरालों पर नमूनों को फैलाएँ

चलो अंतरालों का एक channel बनाकर शुरू करते हैं। जीवन को सरल रखने के लिए, हम केवल 3 अंतरालों का उपयोग करेंगे जिन्हें हम मैन्युअल रूप से परिभाषित करेंगे। एक वास्तविक वर्कफ़्लो में, तुम इन्हें एक फ़ाइल इनपुट से पढ़ सकते हो या यहाँ तक कि बहुत सारी अंतराल फ़ाइलों के साथ एक channel बना सकते हो।

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

अब याद रखो, हम प्रत्येक अंतराल के लिए प्रत्येक नमूने को दोहराना चाहते हैं। इसे कभी-कभी नमूनों और अंतरालों के Cartesian product के रूप में संदर्भित किया जाता है। हम [`combine` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#combine) का उपयोग करके इसे प्राप्त कर सकते हैं। यह channel 1 से प्रत्येक आइटम लेगा और इसे channel 2 में प्रत्येक आइटम के लिए दोहराएगा। चलो अपने वर्कफ़्लो में एक combine ऑपरेटर जोड़ते हैं:

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

अब चलो इसे चलाते हैं और देखते हैं कि क्या होता है:

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

सफलता! हमने हमारी 3 अंतराल सूची में प्रत्येक एकल अंतराल के लिए हर नमूने को दोहराया है। हमने प्रभावी रूप से अपने channel में आइटम की संख्या को तीन गुना कर दिया है।

हालाँकि इसे पढ़ना थोड़ा मुश्किल है, इसलिए अगले सेक्शन में हम इसे साफ करेंगे।

### 4.2. channel को व्यवस्थित करें

हम अपने नमूना डेटा को साफ और रीफैक्टर करने के लिए `map` ऑपरेटर का उपयोग कर सकते हैं ताकि इसे समझना आसान हो। चलो intervals string को पहले तत्व पर joining map में ले जाते हैं।

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

चलो तोड़ते हैं कि यह map ऑपरेशन चरण दर चरण क्या करता है।

पहले, हम कोड को अधिक पठनीय बनाने के लिए नामित पैरामीटर का उपयोग करते हैं। `grouping_key`, `normal`, `tumor` और `interval` नामों का उपयोग करके, हम tuple में तत्वों को इंडेक्स के बजाय नाम से संदर्भित कर सकते हैं:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

अगला, हम `grouping_key` को `interval` फ़ील्ड के साथ संयोजित करते हैं। `grouping_key` एक map है जिसमें `id` और `repeat` फ़ील्ड शामिल हैं। हम `interval` के साथ एक नया map बनाते हैं और उन्हें Groovy के map addition (`+`) का उपयोग करके मर्ज करते हैं:

```groovy
                grouping_key + [interval: interval],
```

अंत में, हम इसे तीन तत्वों के साथ एक tuple के रूप में लौटाते हैं: संयुक्त मेटाडेटा map, normal नमूना फ़ाइल, और tumor नमूना फ़ाइल:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

चलो इसे फिर से चलाते हैं और channel सामग्री की जाँच करते हैं:

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

अपने डेटा को सही संरचना में बाध्य करने के लिए `map` का उपयोग करना मुश्किल हो सकता है, लेकिन यह प्रभावी डेटा मैनिपुलेशन के लिए महत्वपूर्ण है।

अब हमारे पास सभी जीनोमिक अंतरालों में दोहराया गया हर नमूना है, कई स्वतंत्र विश्लेषण इकाइयाँ बनाते हुए जिन्हें समानांतर में प्रोसेस किया जा सकता है। लेकिन क्या होगा यदि हम संबंधित नमूनों को वापस एक साथ लाना चाहते हैं? अगले सेक्शन में, हम सीखेंगे कि उन नमूनों को कैसे समूहित करें जो सामान्य विशेषताओं को साझा करते हैं।

### सारांश

इस सेक्शन में, तुमने सीखा:

- **अंतरालों पर नमूनों को फैलाना**: अंतरालों पर नमूनों को दोहराने के लिए `combine` का उपयोग कैसे करें
- **Cartesian products बनाना**: नमूनों और अंतरालों के सभी संयोजन कैसे उत्पन्न करें
- **Channel संरचना को व्यवस्थित करना**: बेहतर पठनीयता के लिए डेटा को पुनर्संरचित करने के लिए `map` का उपयोग कैसे करें
- **समानांतर प्रोसेसिंग की तैयारी**: वितरित विश्लेषण के लिए डेटा कैसे सेट अप करें

## 5. `groupTuple` का उपयोग करके नमूनों को एकत्रित करना

पिछले सेक्शन में, हमने सीखा कि इनपुट फ़ाइल से डेटा को कैसे विभाजित करें और विशिष्ट फ़ील्ड (हमारे मामले में normal और tumor नमूने) द्वारा फ़िल्टर कैसे करें। लेकिन यह केवल एक प्रकार के joining को कवर करता है। क्या होगा यदि हम किसी विशिष्ट विशेषता द्वारा नमूनों को समूहित करना चाहते हैं? उदाहरण के लिए, मिलान किए गए normal-tumor जोड़ों को जोड़ने के बजाय, हम "sampleA" से सभी नमूनों को उनके प्रकार की परवाह किए बिना एक साथ प्रोसेस करना चाह सकते हैं। यह पैटर्न बायोइनफॉर्मेटिक्स वर्कफ़्लो में आम है जहाँ तुम अंत में परिणामों की तुलना या संयोजन करने से पहले दक्षता कारणों से संबंधित नमूनों को अलग से प्रोसेस करना चाह सकते हो।

Nextflow में ऐसा करने के लिए अंतर्निहित तरीके शामिल हैं, मुख्य जिसे हम देखेंगे वह `groupTuple` है।

चलो हमारे सभी नमूनों को समूहित करके शुरू करते हैं जिनमें समान `id` और `interval` फ़ील्ड हैं, यह एक विश्लेषण का विशिष्ट होगा जहाँ हम तकनीकी प्रतिकृतियों को समूहित करना चाहते हैं लेकिन सार्थक रूप से अलग नमूनों को अलग रखना चाहते हैं।

ऐसा करने के लिए, हमें अपने grouping वेरिएबल को अलग करना चाहिए ताकि हम उन्हें अलगाव में उपयोग कर सकें।

पहला चरण पिछले सेक्शन में हमने जो किया उसके समान है। हमें अपने grouping वेरिएबल को tuple के पहले तत्व के रूप में अलग करना होगा। याद रखो, हमारा पहला तत्व वर्तमान में `id`, `repeat` और `interval` फ़ील्ड का एक map है:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

हम map से अपने `id` और `interval` फ़ील्ड को अलग करने के लिए पहले से `subMap` मेथड का पुन: उपयोग कर सकते हैं। पहले की तरह, हम प्रत्येक नमूने के लिए tuple के पहले तत्व पर `subMap` मेथड लागू करने के लिए `map` ऑपरेटर का उपयोग करेंगे।

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

चलो इसे फिर से चलाते हैं और channel सामग्री की जाँच करते हैं:

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

हम देख सकते हैं कि हमने सफलतापूर्वक `id` और `interval` फ़ील्ड को अलग कर लिया है, लेकिन अभी तक नमूनों को समूहित नहीं किया है।

!!! note

    हम यहाँ `replicate` फ़ील्ड को त्याग रहे हैं। ऐसा इसलिए है क्योंकि हमें इसकी आगे डाउनस्ट्रीम प्रोसेसिंग के लिए आवश्यकता नहीं है। इस ट्यूटोरियल को पूरा करने के बाद, देखो कि क्या तुम बाद के grouping को प्रभावित किए बिना इसे शामिल कर सकते हो!

चलो अब इस नए grouping तत्व द्वारा नमूनों को समूहित करते हैं, [`groupTuple` ऑपरेटर](https://www.nextflow.io/docs/latest/operator.html#grouptuple) का उपयोग करके।

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

बस इतना ही! हमने कोड की केवल एक पंक्ति जोड़ी। चलो देखते हैं कि जब हम इसे चलाते हैं तो क्या होता है:

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

ध्यान दो कि हमारा डेटा संरचना बदल गया है और प्रत्येक channel तत्व के भीतर फ़ाइलें अब `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]` जैसे tuples में निहित हैं। ऐसा इसलिए है क्योंकि जब हम `groupTuple` का उपयोग करते हैं, तो Nextflow एक समूह के प्रत्येक नमूने के लिए एकल फ़ाइलों को संयोजित करता है। डाउनस्ट्रीम में डेटा को संभालने की कोशिश करते समय यह याद रखना महत्वपूर्ण है।

!!! note

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) `groupTuple` के विपरीत है। यह channel में आइटम को अनपैक करता है और उन्हें समतल करता है। `transpose` जोड़ने की कोशिश करो और हमने ऊपर किए गए grouping को पूर्ववत करो!

### सारांश

इस सेक्शन में, तुमने सीखा:

- **संबंधित नमूनों को समूहित करना**: सामान्य विशेषताओं द्वारा नमूनों को एकत्रित करने के लिए `groupTuple` का उपयोग कैसे करें
- **Grouping keys को अलग करना**: grouping के लिए विशिष्ट फ़ील्ड निकालने के लिए `subMap` का उपयोग कैसे करें
- **समूहित डेटा संरचनाओं को संभालना**: `groupTuple` द्वारा बनाई गई नेस्टेड संरचना के साथ कैसे काम करें
- **तकनीकी प्रतिकृति संभालना**: उन नमूनों को कैसे समूहित करें जो समान प्रयोगात्मक स्थितियों को साझा करते हैं

---

## सारांश

इस साइड क्वेस्ट में, तुमने सीखा कि channel का उपयोग करके डेटा को कैसे विभाजित और समूहित करें।

पाइपलाइन के माध्यम से प्रवाहित होने पर डेटा को संशोधित करके, तुम loops या while statements का उपयोग किए बिना एक स्केलेबल पाइपलाइन का निर्माण कर सकते हो, जो अधिक पारंपरिक दृष्टिकोणों पर कई लाभ प्रदान करता है:

- हम बिना किसी अतिरिक्त कोड के जितने चाहें उतने या कम इनपुट तक स्केल कर सकते हैं
- हम iteration के बजाय पाइपलाइन के माध्यम से डेटा के प्रवाह को संभालने पर ध्यान केंद्रित करते हैं
- हम आवश्यकतानुसार जटिल या सरल हो सकते हैं
- पाइपलाइन अधिक घोषणात्मक हो जाती है, यह कैसे होना चाहिए के बजाय क्या होना चाहिए पर ध्यान केंद्रित करती है
- Nextflow स्वतंत्र ऑपरेशन को समानांतर में चलाकर हमारे लिए निष्पादन को अनुकूलित करेगा

इन channel ऑपरेशन में महारत हासिल करना तुम्हें लचीले, स्केलेबल पाइपलाइन बनाने में सक्षम करेगा जो loops या iterative प्रोग्रामिंग का सहारा लिए बिना जटिल डेटा संबंधों को संभालते हैं, Nextflow को निष्पादन को अनुकूलित करने और स्वतंत्र ऑपरेशन को स्वचालित रूप से समानांतर करने की अनुमति देते हैं।

### प्रमुख पैटर्न

1.  **संरचित इनपुट डेटा बनाना:** meta maps के साथ CSV फ़ाइल से शुरू करना ([Metadata in workflows](./metadata.md) से पैटर्न पर निर्माण)

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **डेटा को अलग channel में विभाजित करना:** हमने `type` फ़ील्ड के आधार पर डेटा को स्वतंत्र स्ट्रीम में विभाजित करने के लिए `filter` का उपयोग किया

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **मिलान किए गए नमूनों को जोड़ना:** हमने `id` और `repeat` फ़ील्ड के आधार पर संबंधित नमूनों को फिर से संयोजित करने के लिए `join` का उपयोग किया

    - Key द्वारा दो channel join करें (tuple का पहला तत्व)

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - Joining key निकालें और इस मान द्वारा join करें

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - subMap का उपयोग करके कई फ़ील्ड पर join करें

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **अंतरालों में वितरण:** हमने समानांतर प्रोसेसिंग के लिए जीनोमिक अंतरालों के साथ नमूनों के Cartesian products बनाने के लिए `combine` का उपयोग किया।

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Grouping keys द्वारा एकत्रीकरण:** हमने प्रत्येक tuple में पहले तत्व द्वारा समूहित करने के लिए `groupTuple` का उपयोग किया, जिससे `id` और `interval` फ़ील्ड साझा करने वाले नमूनों को एकत्र किया और तकनीकी प्रतिकृतियों को मर्ज किया।

    ```groovy
    channel.groupTuple()
    ```

6.  **डेटा संरचना को अनुकूलित करना:** हमने विशिष्ट फ़ील्ड निकालने के लिए `subMap` का उपयोग किया और ट्रांसफ़ॉर्मेशन को पुन: प्रयोज्य बनाने के लिए एक नामित closure बनाया।

    - Map से विशिष्ट फ़ील्ड निकालें

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - पुन: प्रयोज्य ट्रांसफ़ॉर्मेशन के लिए नामित closure का उपयोग करें

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

[साइड क्वेस्ट के मेनू](./index.md) पर वापस जाओ या सूची में अगले विषय पर जाने के लिए पृष्ठ के निचले दाएं कोने में बटन पर क्लिक करो।
