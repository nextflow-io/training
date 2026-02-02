# भाग 1: एक डेमो पाइपलाइन चलाएं

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core प्रशिक्षण पाठ्यक्रम के इस पहले भाग में, हम आपको दिखाएंगे कि एक nf-core pipeline कैसे खोजें और आज़माएं, कोड कैसे व्यवस्थित है यह कैसे समझें, और यह पहचानें कि यह plain Nextflow कोड से कैसे अलग है जैसा कि [Hello Nextflow](../hello_nextflow/index.md) में दिखाया गया है।

हम nf-core/demo नामक एक pipeline का उपयोग करने जा रहे हैं जिसे nf-core प्रोजेक्ट अपने pipelines की सूची के हिस्से के रूप में कोड संरचना और टूल संचालन के प्रदर्शन के लिए बनाए रखता है।

सुनिश्चित करें कि आपकी वर्किंग डायरेक्टरी `hello-nf-core/` पर सेट है जैसा कि [Getting started](./00_orientation.md) पेज पर निर्देश दिया गया है।

---

## 1. nf-core/demo pipeline खोजें और प्राप्त करें

आइए [nf-co.re](https://nf-co.re) पर प्रोजेक्ट वेबसाइट पर nf-core/demo pipeline को खोजने से शुरू करें, जो सभी जानकारी को केंद्रीकृत करती है जैसे: सामान्य दस्तावेज़ीकरण और सहायता लेख, प्रत्येक pipeline के लिए दस्तावेज़ीकरण, ब्लॉग पोस्ट, इवेंट की घोषणाएं इत्यादि।

### 1.1. वेबसाइट पर pipeline खोजें

अपने वेब ब्राउज़र में, [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) पर जाएं और सर्च बार में `demo` टाइप करें।

![search results](./img/search-results.png)

Pipeline दस्तावेज़ीकरण पेज तक पहुँचने के लिए pipeline के नाम, `demo`, पर क्लिक करें।

प्रत्येक रिलीज़ की गई pipeline का एक समर्पित पेज होता है जिसमें निम्नलिखित दस्तावेज़ीकरण अनुभाग शामिल होते हैं:

- **Introduction:** Pipeline का परिचय और अवलोकन
- **Usage:** Pipeline को execute करने के तरीके का विवरण
- **Parameters:** विवरण के साथ समूहीकृत pipeline पैरामीटर
- **Output:** अपेक्षित आउटपुट फ़ाइलों का विवरण और उदाहरण
- **Results:** पूर्ण test dataset से उत्पन्न उदाहरण आउटपुट फ़ाइलें
- **Releases & Statistics:** Pipeline संस्करण इतिहास और आंकड़े

जब भी आप एक नई pipeline को अपनाने पर विचार कर रहे हों, तो आपको पहले pipeline दस्तावेज़ीकरण को ध्यान से पढ़ना चाहिए ताकि यह समझ सकें कि यह क्या करती है और इसे चलाने का प्रयास करने से पहले इसे कैसे कॉन्फ़िगर किया जाना चाहिए।

अभी देखें और पता लगाएं:

- Pipeline कौन से टूल चलाएगी (टैब देखें: `Introduction`)
- Pipeline कौन से इनपुट और पैरामीटर स्वीकार करती है या आवश्यक है (टैब देखें: `Parameters`)
- Pipeline द्वारा उत्पादित आउटपुट क्या हैं (टैब देखें: `Output`)

#### 1.1.1. Pipeline अवलोकन

`Introduction` टैब pipeline का एक अवलोकन प्रदान करता है, जिसमें एक दृश्य प्रतिनिधित्व (जिसे subway map कहा जाता है) और pipeline के हिस्से के रूप में चलाए जाने वाले टूल की एक सूची शामिल है।

![pipeline subway map](./img/nf-core-demo-subway-cropped.png)

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

#### 1.1.2. उदाहरण कमांड लाइन

दस्तावेज़ीकरण एक उदाहरण इनपुट फ़ाइल (जिसकी आगे चर्चा की जाएगी) और एक उदाहरण कमांड लाइन भी प्रदान करता है।

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

आप देखेंगे कि उदाहरण कमांड एक workflow फ़ाइल निर्दिष्ट नहीं करती है, केवल pipeline रिपॉज़िटरी का संदर्भ, `nf-core/demo`।

इस तरह से शुरू करने पर, Nextflow मान लेगा कि कोड एक निश्चित तरीके से व्यवस्थित है।
आइए कोड प्राप्त करें ताकि हम इस संरचना की जांच कर सकें।

### 1.2. Pipeline कोड प्राप्त करें

एक बार जब हम यह निर्धारित कर लेते हैं कि pipeline हमारे उद्देश्यों के लिए उपयुक्त प्रतीत होती है, तो आइए इसे आज़माएं।
सौभाग्य से Nextflow सही ढंग से प्रारूपित रिपॉज़िटरी से pipelines को मैन्युअल रूप से कुछ भी डाउनलोड किए बिना प्राप्त करना आसान बनाता है।

आइए टर्मिनल पर वापस जाएं और निम्नलिखित चलाएं:

```bash
nextflow pull nf-core/demo
```

??? success "कमांड आउटपुट"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
    ```

Nextflow pipeline कोड का एक `pull` करता है, मतलब यह पूरी रिपॉज़िटरी को आपकी लोकल ड्राइव पर डाउनलोड करता है।

स्पष्ट करने के लिए, आप GitHub में उचित रूप से सेट की गई किसी भी Nextflow pipeline के साथ ऐसा कर सकते हैं, न केवल nf-core pipelines के साथ।
हालाँकि nf-core Nextflow pipelines का सबसे बड़ा ओपन-सोर्स संग्रह है।

आप Nextflow से इस तरह से प्राप्त की गई pipelines की सूची दे सकते हैं:

```bash
nextflow list
```

??? success "कमांड आउटपुट"

    ```console
    nf-core/demo
    ```

आप देखेंगे कि फ़ाइलें आपकी वर्तमान कार्य डायरेक्टरी में नहीं हैं।
डिफ़ॉल्ट रूप से, Nextflow उन्हें `$NXF_HOME/assets` में सहेजता है।

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="डायरेक्टरी सामग्री"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo

2 directories, 0 files
```

!!! note

    यदि आप हमारे प्रशिक्षण वातावरण का उपयोग नहीं कर रहे हैं तो आपके सिस्टम पर पूर्ण पथ भिन्न हो सकता है।

Nextflow डाउनलोड किए गए स्रोत कोड को जानबूझकर 'बाहर' रखता है इस सिद्धांत पर कि इन pipelines का उपयोग उस कोड की तुलना में अधिक लाइब्रेरी की तरह किया जाना चाहिए जिसके साथ आप सीधे इंटरैक्ट करेंगे।

हालांकि, इस प्रशिक्षण के उद्देश्यों के लिए, हम इधर-उधर देखना और देखना चाहते हैं कि वहाँ क्या है।
तो इसे आसान बनाने के लिए, आइए अपनी वर्तमान कार्य डायरेक्टरी से उस स्थान पर एक symbolic link बनाएं।

```bash
ln -s $NXF_HOME/assets pipelines
```

यह एक शॉर्टकट बनाता है जो हमने अभी डाउनलोड किए गए कोड को एक्सप्लोर करना आसान बनाता है।

```bash
tree -L 2 pipelines
```

```console title="डायरेक्टरी सामग्री"
pipelines
└── nf-core
    └── demo

2 directories, 0 files
```

अब हम आवश्यकतानुसार स्रोत कोड में अधिक आसानी से झाँक सकते हैं।

लेकिन पहले, आइए अपनी पहली nf-core pipeline चलाने का प्रयास करें!

### मुख्य बात

अब आप जानते हैं कि nf-core वेबसाइट के माध्यम से एक pipeline कैसे खोजें और स्रोत कोड की एक स्थानीय प्रति कैसे प्राप्त करें।

### आगे क्या है?

सीखें कि न्यूनतम प्रयास के साथ एक nf-core pipeline कैसे आज़माएं।

---

## 2. अपने test profile के साथ pipeline आज़माएं

सुविधाजनक रूप से, प्रत्येक nf-core pipeline एक test profile के साथ आती है।
यह [nf-core/test-datasets](https://github.com/nf-core/test-datasets) रिपॉज़िटरी में होस्ट किए गए एक छोटे test dataset का उपयोग करके चलने के लिए pipeline के लिए कॉन्फ़िगरेशन सेटिंग्स का एक न्यूनतम सेट है।
यह छोटे पैमाने पर एक pipeline को जल्दी से आज़माने का एक शानदार तरीका है।

!!! note

    Nextflow का configuration profile सिस्टम आपको विभिन्न कंटेनर इंजन या execution वातावरण के बीच आसानी से स्विच करने की अनुमति देता है।
    अधिक विवरण के लिए, [Hello Nextflow भाग 6: Configuration](../hello_nextflow/06_hello_config.md) देखें।

### 2.1. test profile की जांच करें

Pipeline के test profile को चलाने से पहले यह जांचना अच्छा अभ्यास है कि यह क्या निर्दिष्ट करता है।
`nf-core/demo` के लिए `test` profile कॉन्फ़िगरेशन फ़ाइल `conf/test.config` में रहता है और नीचे दिखाया गया है।

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

process {
    resourceLimits = [
        cpus: 4,
        memory: '4.GB',
        time: '1.h'
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // इनपुट डेटा
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

}
```

आप तुरंत देखेंगे कि शीर्ष पर कमेंट ब्लॉक में एक उपयोग उदाहरण शामिल है जो दिखाता है कि इस test profile के साथ pipeline को कैसे चलाना है।

```groovy title="conf/test.config" linenums="7"
Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

हमें केवल वही प्रदान करने की आवश्यकता है जो उदाहरण कमांड में कैरेट के बीच दिखाया गया है: `<docker/singularity>` और `<OUTDIR>`।

याद दिलाने के लिए, `<docker/singularity>` कंटेनर सिस्टम की पसंद को संदर्भित करता है। सभी nf-core pipelines को reproducibility सुनिश्चित करने और सॉफ़्टवेयर इंस्टॉलेशन समस्याओं को खत्म करने के लिए कंटेनर (Docker, Singularity, आदि) के साथ उपयोग करने योग्य होने के लिए डिज़ाइन किया गया है।
इसलिए हमें यह निर्दिष्ट करने की आवश्यकता होगी कि हम pipeline का परीक्षण करने के लिए Docker या Singularity का उपयोग करना चाहते हैं या नहीं।

`--outdir <OUTDIR>` भाग उस डायरेक्टरी को संदर्भित करता है जहाँ Nextflow pipeline के आउटपुट लिखेगा।
हमें इसके लिए एक नाम प्रदान करने की आवश्यकता है, जिसे हम बस बना सकते हैं।
यदि यह पहले से मौजूद नहीं है, तो Nextflow इसे runtime पर हमारे लिए बना देगा।

कमेंट ब्लॉक के बाद के अनुभाग की ओर बढ़ते हुए, test profile हमें दिखाता है कि परीक्षण के लिए क्या पूर्व-कॉन्फ़िगर किया गया है: सबसे विशेष रूप से, `input` पैरामीटर पहले से ही एक test dataset की ओर इशारा करने के लिए सेट है, इसलिए हमें अपना डेटा प्रदान करने की आवश्यकता नहीं है।
यदि आप पूर्व-कॉन्फ़िगर किए गए इनपुट के लिंक का अनुसरण करते हैं, तो आप देखेंगे कि यह एक csv फ़ाइल है जिसमें कई प्रायोगिक नमूनों के लिए नमूना पहचानकर्ता और फ़ाइल पथ हैं।

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

इसे samplesheet कहा जाता है, और यह nf-core pipelines के लिए इनपुट का सबसे आम रूप है।

!!! note

    यदि आप डेटा प्रारूप और प्रकारों से परिचित नहीं हैं तो चिंता न करें, यह आगे के लिए महत्वपूर्ण नहीं है।

तो यह पुष्टि करता है कि हमारे पास pipeline को आज़माने के लिए आवश्यक सब कुछ है।

### 2.2. Pipeline चलाएं

आइए कंटेनर सिस्टम के लिए Docker का उपयोग करने और आउटपुट डायरेक्टरी के रूप में `demo-results` का निर्णय लें, और हम test कमांड चलाने के लिए तैयार हैं:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: db7f526ce1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.0.2
    ------------------------------------------------------
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-57-41

    Core Nextflow options
      revision                  : master
      runName                   : magical_pauling
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/demo
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/.nextflow/assets/nf-core/demo/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md


    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

यदि आपका आउटपुट उससे मेल खाता है, तो बधाई हो! आपने अभी अपनी पहली nf-core pipeline चलाई है।

आप देखेंगे कि जब आप एक बुनियादी Nextflow pipeline चलाते हैं तो कंसोल आउटपुट बहुत अधिक है।
एक हेडर है जिसमें pipeline के संस्करण, इनपुट और आउटपुट, और कॉन्फ़िगरेशन के कुछ तत्वों का सारांश शामिल है।

!!! note

    आपका आउटपुट अलग-अलग timestamps, execution नाम और फ़ाइल पथ दिखाएगा, लेकिन समग्र संरचना और process execution समान होनी चाहिए।

Execution आउटपुट की ओर बढ़ते हुए, आइए उन लाइनों पर एक नज़र डालें जो हमें बताती हैं कि कौन से processes चलाए गए थे:

```console
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
```

यह हमें बताता है कि तीन processes चलाई गईं, जो nf-core वेबसाइट पर pipeline दस्तावेज़ीकरण पेज में दिखाए गए तीन टूल से संबंधित हैं: FASTQC, SEQTK_TRIM और MULTIQC।

पूर्ण process नाम जैसा कि यहाँ दिखाया गया है, जैसे `NFCORE_DEMO:DEMO:MULTIQC`, परिचयात्मक Hello Nextflow सामग्री में आपने जो देखा होगा उससे लंबे हैं।
इनमें उनके parent workflows के नाम शामिल हैं और pipeline कोड की modularity को दर्शाते हैं।
हम इसके बारे में थोड़ी देर में अधिक विस्तार से बात करेंगे।

### 2.3. Pipeline के आउटपुट की जांच करें

अंत में, आइए pipeline द्वारा उत्पादित `demo-results` डायरेक्टरी पर एक नज़र डालें।

```bash
tree -L 2 demo-results
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    demo-results
    ├── fastqc
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── fq
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── multiqc
    │   ├── multiqc_data
    │   ├── multiqc_plots
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2025-11-21_04-57-41.html
        ├── execution_timeline_2025-11-21_04-57-41.html
        ├── execution_trace_2025-11-21_04-57-41.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2025-11-21_04-57-46.json
        └── pipeline_dag_2025-11-21_04-57-41.html
    ```

यह बहुत अधिक लग सकता है।
`nf-core/demo` pipeline के आउटपुट के बारे में अधिक जानने के लिए, इसका [दस्तावेज़ीकरण पेज](https://nf-co.re/demo/1.0.2/docs/output/) देखें।

इस स्तर पर, यह देखना महत्वपूर्ण है कि परिणाम module द्वारा व्यवस्थित हैं, और इसके अतिरिक्त `pipeline_info` नामक एक डायरेक्टरी है जिसमें pipeline execution के बारे में विभिन्न timestamped रिपोर्ट हैं।

उदाहरण के लिए, `execution_timeline_*` फ़ाइल आपको दिखाती है कि कौन से processes चलाए गए, किस क्रम में और उन्हें चलाने में कितना समय लगा:

![execution timeline report](./img/execution_timeline.png)

!!! note

    यहाँ कार्य समानांतर में नहीं चलाए गए क्योंकि हम Github Codespaces में एक minimalist मशीन पर चल रहे हैं।
    इन्हें समानांतर में चलते हुए देखने के लिए, अपने codespace के CPU allocation और test configuration में resource limits को बढ़ाने का प्रयास करें।

ये रिपोर्ट सभी nf-core pipelines के लिए स्वचालित रूप से उत्पन्न होती हैं।

### मुख्य बात

आप जानते हैं कि अपने अंतर्निहित test profile का उपयोग करके एक nf-core pipeline कैसे चलाएं और इसके आउटपुट कहाँ खोजें।

### आगे क्या है?

जानें कि pipeline कोड कैसे व्यवस्थित है।

---

## 3. Pipeline कोड संरचना की जांच करें

अब जब हमने उपयोगकर्ताओं के रूप में pipeline को सफलतापूर्वक चलाया है, तो आइए अपने दृष्टिकोण को बदलकर देखें कि nf-core pipelines आंतरिक रूप से कैसे संरचित हैं।

nf-core प्रोजेक्ट pipelines के संरचित होने के तरीके, और कोड को कैसे व्यवस्थित, कॉन्फ़िगर और दस्तावेज़ित किया जाता है, के लिए मजबूत दिशानिर्देश लागू करता है।
यह सब कैसे व्यवस्थित है यह समझना अपने स्वयं के nf-core-compatible pipelines विकसित करने की दिशा में पहला कदम है, जिसे हम इस पाठ्यक्रम के भाग 2 में निपटाएंगे।

आइए `pipelines` symlink का उपयोग करके `nf-core/demo` रिपॉज़िटरी में pipeline कोड कैसे व्यवस्थित है, इस पर एक नज़र डालें जिसे हमने पहले बनाया था।

आप या तो `tree` का उपयोग कर सकते हैं या `nf-core/demo` डायरेक्टरी को खोजने और खोलने के लिए फ़ाइल एक्सप्लोरर का उपयोग कर सकते हैं।

```bash
tree -L 1 pipelines/nf-core/demo
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    pipelines/nf-core/demo
    ├── assets
    ├── CHANGELOG.md
    ├── CITATIONS.md
    ├── CODE_OF_CONDUCT.md
    ├── conf
    ├── docs
    ├── LICENSE
    ├── main.nf
    ├── modules
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── nf-test.config
    ├── README.md
    ├── ro-crate-metadata.json
    ├── subworkflows
    ├── tests
    ├── tower.yml
    └── workflows
    ```

वहाँ बहुत कुछ चल रहा है, इसलिए हम इसे चरण दर चरण निपटाएंगे।

सबसे पहले, आइए ध्यान दें कि शीर्ष स्तर पर, आप सारांश जानकारी के साथ एक README फ़ाइल पा सकते हैं, साथ ही साथ सहायक फ़ाइलें जो प्रोजेक्ट की जानकारी जैसे लाइसेंसिंग, योगदान दिशानिर्देश, उद्धरण और आचार संहिता को सारांशित करती हैं।
विस्तृत pipeline दस्तावेज़ीकरण `docs` डायरेक्टरी में स्थित है।
यह सभी सामग्री nf-core वेबसाइट पर वेब पेज उत्पन्न करने के लिए प्रोग्रामेटिक रूप से उपयोग की जाती है, इसलिए वे हमेशा कोड के साथ अद्यतित रहते हैं।

अब, बाकी के लिए, हम अपने अन्वेषण को तीन चरणों में विभाजित करने जा रहे हैं:

1. Pipeline कोड घटक (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Pipeline कॉन्फ़िगरेशन
3. इनपुट और सत्यापन

आइए pipeline कोड घटकों से शुरू करें।
हम फ़ाइल पदानुक्रम और संरचनात्मक संगठन पर ध्यान केंद्रित करने जा रहे हैं, बजाय व्यक्तिगत फ़ाइलों के भीतर कोड में गोता लगाने के।

### 3.1. Pipeline कोड घटक

मानक nf-core pipeline कोड संगठन एक modular संरचना का अनुसरण करता है जो कोड पुन: उपयोग को अधिकतम करने के लिए डिज़ाइन किया गया है, जैसा कि [Hello Modules](../hello_nextflow/04_hello_modules.md), [Hello Nextflow](../hello_nextflow/index.md) पाठ्यक्रम के भाग 4 में पेश किया गया है, हालांकि सच्चे nf-core फैशन में, यह थोड़ी अतिरिक्त जटिलता के साथ लागू किया गया है।
विशेष रूप से, nf-core pipelines subworkflows का प्रचुर उपयोग करते हैं, अर्थात् workflow scripts जो एक parent workflow द्वारा import की जाती हैं।

यह थोड़ा अमूर्त लग सकता है, तो आइए देखें कि `nf-core/demo` pipeline में इसका व्यवहार में कैसे उपयोग किया जाता है।

!!! note

    हम इन modular घटकों को कैसे जोड़ा जाता है इसके लिए वास्तविक कोड पर नहीं जाएंगे, क्योंकि subworkflows के उपयोग से जुड़ी कुछ अतिरिक्त जटिलता है जो भ्रमित करने वाली हो सकती है, और उसे समझना प्रशिक्षण के इस चरण में आवश्यक नहीं है।
    फिलहाल, हम समग्र संगठन और तर्क पर ध्यान केंद्रित करने जा रहे हैं।

#### 3.1.1. सामान्य अवलोकन

यहाँ `nf-core/demo` pipeline के लिए प्रासंगिक कोड घटकों के बीच संबंध कैसा दिखता है:

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/nf-core_demo_code_organization.svg"
</figure>

एक तथाकथित _entrypoint_ स्क्रिप्ट है जिसे `main.nf` कहा जाता है, जो दो प्रकार के nested workflows के लिए एक wrapper के रूप में कार्य करता है: वास्तविक विश्लेषण तर्क वाली workflow, जो `workflows/` के तहत स्थित है और `demo.nf` कहलाती है, और housekeeping workflows का एक सेट `subworkflows/` के तहत स्थित है।
`demo.nf` workflow `modules/` के तहत स्थित **modules** को कॉल करती है; इनमें **processes** हैं जो वास्तविक विश्लेषण चरण करेंगी।

!!! note

    Subworkflows housekeeping कार्यों तक सीमित नहीं हैं, और वे process modules का उपयोग कर सकते हैं।

    यहाँ दिखाई गई `nf-core/demo` pipeline spectrum पर सरल पक्ष पर होती है, लेकिन अन्य nf-core pipelines (जैसे `nf-core/rnaseq`) subworkflows का उपयोग करती हैं जो वास्तविक विश्लेषण में शामिल हैं।

अब, आइए इन घटकों की बारी-बारी से समीक्षा करें।

#### 3.1.2. Entrypoint स्क्रिप्ट: `main.nf`

`main.nf` स्क्रिप्ट वह entrypoint है जहाँ से Nextflow शुरू होता है जब हम `nextflow run nf-core/demo` execute करते हैं।
इसका मतलब है कि जब आप pipeline चलाने के लिए `nextflow run nf-core/demo` चलाते हैं, तो Nextflow स्वचालित रूप से `main.nf` स्क्रिप्ट ढूंढता है और execute करता है।
यह इस पारंपरिक नामकरण और संरचना का पालन करने वाली किसी भी Nextflow pipeline के लिए काम करता है, न केवल nf-core pipelines के लिए।

एक entrypoint स्क्रिप्ट का उपयोग करने से वास्तविक विश्लेषण स्क्रिप्ट के चलने से पहले और बाद में मानकीकृत 'housekeeping' subworkflows को चलाना आसान हो जाता है।
हम वास्तविक विश्लेषण workflow और इसके modules की समीक्षा करने के बाद उन पर जाएंगे।

#### 3.1.3. विश्लेषण स्क्रिप्ट: `workflows/demo.nf`

`workflows/demo.nf` workflow वह जगह है जहाँ pipeline का केंद्रीय तर्क संग्रहीत है।
यह सामान्य Nextflow workflow की तरह संरचित है, सिवाय इसके कि इसे एक parent workflow से कॉल करने के लिए डिज़ाइन किया गया है, जिसके लिए कुछ अतिरिक्त सुविधाओं की आवश्यकता होती है।
हम इस पाठ्यक्रम के अगले भाग में प्रासंगिक अंतरों को कवर करेंगे, जब हम Hello Nextflow से साधारण Hello pipeline को nf-core-compatible रूप में परिवर्तित करने का सामना करेंगे।

`demo.nf` workflow `modules/` के तहत स्थित **modules** को कॉल करती है, जिसकी हम अगली बार समीक्षा करेंगे।

!!! note

    कुछ nf-core विश्लेषण workflows निम्न-स्तरीय subworkflows को कॉल करके nesting के अतिरिक्त स्तर प्रदर्शित करती हैं।
    यह ज्यादातर दो या अधिक modules को wrap करने के लिए उपयोग किया जाता है जो आमतौर पर एक साथ आसानी से पुन: उपयोग योग्य pipeline segments में उपयोग किए जाते हैं।
    आप nf-core वेबसाइट पर उपलब्ध [nf-core subworkflows](https://nf-co.re/subworkflows/) ब्राउज़ करके कुछ उदाहरण देख सकते हैं।

    जब विश्लेषण स्क्रिप्ट subworkflows का उपयोग करती है, तो वे `subworkflows/` डायरेक्टरी के तहत संग्रहीत होती हैं।

#### 3.1.4. Modules

Modules वह जगह हैं जहाँ process कोड रहता है, जैसा कि [Hello Nextflow प्रशिक्षण पाठ्यक्रम के भाग 4](../hello_nextflow/04_hello_modules.md) में वर्णित है।

nf-core प्रोजेक्ट में, modules को एक बहु-स्तरीय nested संरचना का उपयोग करके व्यवस्थित किया जाता है जो उनके मूल और उनकी सामग्री दोनों को दर्शाता है।
शीर्ष स्तर पर, modules को या तो `nf-core` या `local` (nf-core प्रोजेक्ट का हिस्सा नहीं) के रूप में विभेदित किया जाता है, और फिर उन टूल के नाम पर एक डायरेक्टरी में रखा जाता है जिन्हें वे wrap करते हैं।
यदि टूल एक toolkit से संबंधित है (अर्थात् एक पैकेज जिसमें कई टूल हैं) तो toolkit के नाम पर एक intermediate डायरेक्टरी स्तर है।

आप इसे `nf-core/demo` pipeline modules में व्यवहार में लागू देख सकते हैं:

```bash
tree -L 3 pipelines/nf-core/demo/modules
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    pipelines/nf-core/demo/modules
    └── nf-core
        ├── fastqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── multiqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── seqtk
            └── trim

    7 directories, 6 files
    ```

यहाँ आप देखते हैं कि `fastqc` और `multiqc` modules `nf-core` modules के भीतर शीर्ष स्तर पर बैठते हैं, जबकि `trim` module उस toolkit के तहत बैठता है जिससे यह संबंधित है, `seqtk`।
इस मामले में कोई `local` modules नहीं हैं।

Process का वर्णन करने वाली module कोड फ़ाइल हमेशा `main.nf` कहलाती है, और tests और `.yml` फ़ाइलों के साथ होती है जिन्हें हम अभी के लिए अनदेखा करेंगे।

एक साथ लिया जाए तो, entrypoint workflow, विश्लेषण workflow और modules pipeline के 'दिलचस्प' भागों को चलाने के लिए पर्याप्त हैं।
हालाँकि, हम जानते हैं कि वहाँ housekeeping subworkflows भी हैं, तो आइए अब उन्हें देखें।

#### 3.1.5. Housekeeping subworkflows

Modules की तरह, subworkflows को `local` और `nf-core` डायरेक्टरी में विभेदित किया जाता है, और प्रत्येक subworkflow की अपनी nested डायरेक्टरी संरचना है जिसमें अपनी `main.nf` स्क्रिप्ट, tests और `.yml` फ़ाइल है।

```bash
tree -L 3 pipelines/nf-core/demo/subworkflows
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    pipelines/nf-core/demo/subworkflows
    ├── local
    │   └── utils_nfcore_demo_pipeline
    │       └── main.nf
    └── nf-core
        ├── utils_nextflow_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── utils_nfcore_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── utils_nfschema_plugin
            ├── main.nf
            ├── meta.yml
            └── tests

    9 directories, 7 files
    ```

जैसा कि ऊपर बताया गया है, `nf-core/demo` pipeline में कोई विश्लेषण-विशिष्ट subworkflows शामिल नहीं हैं, इसलिए यहाँ हम जो सभी subworkflows देखते हैं वे तथाकथित 'housekeeping' या 'utility' workflows हैं, जैसा कि उनके नामों में `utils_` prefix द्वारा दर्शाया गया है।
ये subworkflows वे हैं जो कंसोल आउटपुट में fancy nf-core header उत्पन्न करती हैं, अन्य सहायक कार्यों के बीच।

!!! tip

    उनके नामकरण पैटर्न के अलावा, एक और संकेत है कि ये subworkflows कोई वास्तव में विश्लेषण-संबंधित कार्य नहीं करती हैं वह यह है कि वे किसी भी processes को बिल्कुल नहीं कॉल करती हैं।

यह `nf-core/demo` pipeline का गठन करने वाले मुख्य कोड घटकों के राउंड-अप को पूरा करता है।
अब आइए शेष तत्वों पर एक नज़र डालें जिनके बारे में आपको development में गोता लगाने से पहले थोड़ा पता होना चाहिए: pipeline कॉन्फ़िगरेशन और इनपुट सत्यापन।

### 3.2. Pipeline कॉन्फ़िगरेशन

आपने पहले सीखा है कि Nextflow pipeline execution को कॉन्फ़िगर करने के लिए कई विकल्प प्रदान करता है, चाहे वह इनपुट और पैरामीटर के संदर्भ में हो, कंप्यूटिंग संसाधन हों, और orchestration के अन्य पहलू हों।
nf-core प्रोजेक्ट pipeline कॉन्फ़िगरेशन के लिए अत्यधिक मानकीकृत दिशानिर्देश लागू करता है जो Nextflow के लचीले customization विकल्पों पर एक ऐसे तरीके से निर्माण करने का लक्ष्य रखते हैं जो pipelines में अधिक consistency और maintainability प्रदान करता है।

केंद्रीय कॉन्फ़िगरेशन फ़ाइल `nextflow.config` का उपयोग पैरामीटर और अन्य कॉन्फ़िगरेशन विकल्पों के लिए डिफ़ॉल्ट मान सेट करने के लिए किया जाता है।
इनमें से अधिकांश कॉन्फ़िगरेशन विकल्प डिफ़ॉल्ट रूप से लागू होते हैं जबकि अन्य (उदाहरण के लिए, सॉफ़्टवेयर dependency profiles) वैकल्पिक profiles के रूप में शामिल हैं।

कई अतिरिक्त कॉन्फ़िगरेशन फ़ाइलें हैं जो `conf` फ़ोल्डर में संग्रहीत हैं और जिन्हें डिफ़ॉल्ट रूप से या वैकल्पिक रूप से profiles के रूप में कॉन्फ़िगरेशन में जोड़ा जा सकता है:

- `base.config`: एक 'blank slate' config फ़ाइल, अधिकांश high-performance computing वातावरण पर सामान्य उपयोग के लिए उपयुक्त। यह संसाधन उपयोग के व्यापक bins को परिभाषित करती है, उदाहरण के लिए, जो modules पर लागू करने के लिए सुविधाजनक हैं।
- `modules.config`: अतिरिक्त module निर्देश और arguments।
- `test.config`: न्यूनतम test डेटा के साथ pipeline चलाने के लिए एक profile, जिसका उपयोग हमने demo pipeline चलाते समय किया था।
- `test_full.config`: पूर्ण आकार के test dataset के साथ pipeline चलाने के लिए एक profile।

हम पाठ्यक्रम में बाद में इनमें से कुछ फ़ाइलों को छुएंगे।

### 3.3. इनपुट और सत्यापन

जैसा कि हमने पहले बताया, जब हमने `nf-core/demo` pipeline के test profile की जांच की, यह इनपुट के रूप में एक samplesheet लेने के लिए डिज़ाइन की गई है जिसमें फ़ाइल पथ और नमूना पहचानकर्ता होते हैं।
फ़ाइल पथ `nf-core/test-datasets` रिपॉज़िटरी में स्थित वास्तविक डेटा से जुड़े थे।

एक उदाहरण samplesheet `assets` डायरेक्टरी के तहत भी प्रदान की गई है, हालांकि इसमें पथ वास्तविक नहीं हैं।

```csv title="assets/samplesheet.csv" linenums="1"
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_SINGLE_END,/path/to/fastq/files/AEG588A4_S4_L003_R1_001.fastq.gz,

```

यह विशेष samplesheet काफी सरल है, लेकिन कुछ pipelines samplesheets पर चलती हैं जो अधिक जटिल हैं, प्राथमिक इनपुट से जुड़े बहुत अधिक metadata के साथ।

दुर्भाग्य से, क्योंकि इन फ़ाइलों को आंख से जांचना मुश्किल हो सकता है, इनपुट डेटा की अनुचित formatting pipeline विफलताओं का एक बहुत ही सामान्य स्रोत है।
एक संबंधित समस्या तब होती है जब पैरामीटर गलत तरीके से प्रदान किए जाते हैं।

इन समस्याओं का समाधान सभी इनपुट फ़ाइलों पर स्वचालित सत्यापन जांच चलाना है ताकि यह सुनिश्चित किया जा सके कि उनमें अपेक्षित प्रकार की जानकारी है, सही ढंग से प्रारूपित है, और पैरामीटर पर यह सुनिश्चित करने के लिए कि वे अपेक्षित प्रकार के हैं।
इसे इनपुट सत्यापन कहा जाता है, और आदर्श रूप से pipeline चलाने की कोशिश करने से _पहले_ किया जाना चाहिए, बजाय यह पता लगाने के लिए pipeline के विफल होने की प्रतीक्षा करने के कि इनपुट के साथ कोई समस्या थी।

कॉन्फ़िगरेशन की तरह, nf-core प्रोजेक्ट इनपुट सत्यापन के बारे में बहुत opinionated है, और [nf-schema plugin](https://nextflow-io.github.io/nf-schema/latest/) के उपयोग की सिफारिश करता है, एक Nextflow plugin जो Nextflow pipelines के लिए व्यापक सत्यापन क्षमताएं प्रदान करता है।

हम इस पाठ्यक्रम के भाग 5 में इस विषय को अधिक विस्तार से कवर करेंगे।
अभी के लिए, बस जागरूक रहें कि उस उद्देश्य के लिए दो JSON फ़ाइलें प्रदान की गई हैं, `nextflow_schema.json` और `assets/schema_input.json`।

`nextflow_schema.json` एक फ़ाइल है जिसका उपयोग pipeline पैरामीटर के बारे में जानकारी संग्रहीत करने के लिए किया जाता है जिसमें प्रकार, विवरण और machine readable format में सहायता text शामिल है।
इसका उपयोग विभिन्न उद्देश्यों के लिए किया जाता है, जिसमें स्वचालित पैरामीटर सत्यापन, सहायता text generation, और UI इंटरफेस में interactive पैरामीटर फॉर्म rendering शामिल है।

`schema_input.json` एक फ़ाइल है जिसका उपयोग इनपुट samplesheet संरचना को परिभाषित करने के लिए किया जाता है।
प्रत्येक column का एक प्रकार, पैटर्न, विवरण और machine readable format में सहायता text हो सकता है।
स्कीमा का उपयोग विभिन्न उद्देश्यों के लिए किया जाता है, जिसमें स्वचालित सत्यापन और सहायक त्रुटि संदेश प्रदान करना शामिल है।

### मुख्य बात

आप जानते हैं कि nf-core pipeline के मुख्य घटक क्या हैं और कोड कैसे व्यवस्थित है; कॉन्फ़िगरेशन के मुख्य तत्व कहाँ स्थित हैं; और आप जागरूक हैं कि इनपुट सत्यापन किसके लिए है।

### आगे क्या है?

एक ब्रेक लें! वह बहुत कुछ था। जब आप तरोताज़ा महसूस करें और तैयार हों, तो एक nf-core compatible pipeline लिखने के लिए जो आपने सीखा है उसे लागू करने के लिए अगले अनुभाग पर जाएं।

!!! tip

    यदि आप अगले भाग पर जाने से पहले subworkflows के साथ workflows बनाने का तरीका सीखना चाहते हैं, तो [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Side Quest देखें।
