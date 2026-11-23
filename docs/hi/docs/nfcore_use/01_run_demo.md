# भाग 1: एक डेमो पाइपलाइन चलाएं

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Use nf-core प्रशिक्षण कोर्स के इस पहले भाग में, हम तुम्हें दिखाएंगे कि nf-core पाइपलाइन कैसे खोजें और इसके built-in test profile का उपयोग करके इसे कैसे आज़माएं।

हम `nf-core/demo` नामक एक पाइपलाइन का उपयोग करने जा रहे हैं जिसे nf-core प्रोजेक्ट द्वारा डेमोंस्ट्रेशन और प्रशिक्षण उद्देश्यों के लिए अपनी पाइपलाइनों की सूची के हिस्से के रूप में बनाए रखा जाता है।

सुनिश्चित करो कि तुम्हारी working directory `nfcore-use/` पर सेट है जैसा कि [Getting started](./00_orientation.md) पेज पर बताया गया है।

---

## 1. nf-core/demo पाइपलाइन खोजें और प्राप्त करें

आइए [nf-co.re](https://nf-co.re) पर प्रोजेक्ट वेबसाइट पर nf-core/demo पाइपलाइन को खोजकर शुरू करें, जो सभी जानकारी को केंद्रीकृत करती है जैसे: सामान्य दस्तावेज़ीकरण और सहायता लेख, प्रत्येक पाइपलाइन के लिए दस्तावेज़ीकरण, ब्लॉग पोस्ट, इवेंट घोषणाएं आदि।

### 1.1. वेबसाइट पर पाइपलाइन खोजें

अपने वेब ब्राउज़र में, [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) पर जाओ और सर्च बार में `demo` टाइप करो।

![search results](./img/search-results.png)

पाइपलाइन डॉक्यूमेंटेशन पेज तक पहुंचने के लिए पाइपलाइन के नाम, `demo`, पर क्लिक करो।

प्रत्येक released पाइपलाइन का एक dedicated पेज होता है जिसमें निम्नलिखित डॉक्यूमेंटेशन सेक्शन शामिल हैं:

- **Introduction:** पाइपलाइन का परिचय और अवलोकन
- **Usage:** पाइपलाइन को execute करने के तरीकों का विवरण
- **Parameters:** विवरण के साथ समूहीकृत पाइपलाइन पैरामीटर
- **Output:** अपेक्षित आउटपुट फ़ाइलों का विवरण और उदाहरण
- **Results:** पूर्ण test dataset से उत्पन्न उदाहरण आउटपुट फ़ाइलें
- **Releases & Statistics:** पाइपलाइन संस्करण इतिहास और आंकड़े

जब भी तुम किसी नई पाइपलाइन को अपनाने पर विचार कर रहे हो, तो इसे चलाने का प्रयास करने से पहले पाइपलाइन डॉक्यूमेंटेशन को ध्यान से पढ़ना चाहिए ताकि यह समझ सको कि यह क्या करती है और इसे कैसे कॉन्फ़िगर किया जाना चाहिए।

अभी एक नज़र डालो और देखो कि क्या तुम पता लगा सकते हो:

- पाइपलाइन कौन से टूल चलाएगी (टैब देखो: `Introduction`)
- पाइपलाइन कौन से इनपुट और पैरामीटर स्वीकार करती है या आवश्यक है (टैब देखो: `Parameters`)
- पाइपलाइन द्वारा उत्पादित आउटपुट क्या हैं (टैब देखो: `Output`)

#### 1.1.1. पाइपलाइन का अवलोकन

`Introduction` टैब पाइपलाइन का अवलोकन प्रदान करता है, जिसमें एक दृश्य प्रतिनिधित्व (जिसे subway map कहा जाता है) और पाइपलाइन के हिस्से के रूप में चलाए जाने वाले टूल की सूची शामिल है।

![pipeline subway map](./img/nf-core-demo-subway-cropped.png)

1. Read QC ([FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter और quality trimming ([SEQTK_TRIM](https://github.com/lh3/seqtk))
3. Raw reads के लिए QC प्रस्तुत करें ([MULTIQC](http://multiqc.info/))
4. एक गाय से हल्का-फुल्का टेक्स्ट संदेश उत्पन्न करें ([COWPY](https://github.com/jeffbuttars/cowpy))

#### 1.1.2. उदाहरण command line

डॉक्यूमेंटेशन एक उदाहरण इनपुट फ़ाइल (नीचे और चर्चा की गई) और एक उदाहरण command line भी प्रदान करता है।

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

तुम देखोगे कि उदाहरण कमांड एक workflow फ़ाइल निर्दिष्ट नहीं करता, बस पाइपलाइन रिपॉजिटरी का संदर्भ, `nf-core/demo`।

इस तरह invoke करने पर, Nextflow मान लेगा कि कोड एक निश्चित तरीके से व्यवस्थित है।
आइए कोड प्राप्त करें ताकि हम इस संरचना की जांच कर सकें।

### 1.2. पाइपलाइन कोड प्राप्त करें

एक बार जब हमने निर्धारित कर लिया कि पाइपलाइन हमारे उद्देश्यों के लिए उपयुक्त लगती है, तो आइए इसे आज़माएं।
सौभाग्य से Nextflow सही तरीके से formatted रिपॉजिटरी से पाइपलाइन प्राप्त करना आसान बनाता है बिना कुछ भी मैन्युअल रूप से डाउनलोड किए।

#### 1.2.1. `nextflow pull` का उपयोग करें

आइए टर्मिनल पर वापस जाएं और निम्नलिखित चलाएं:

```bash
nextflow pull nf-core/demo
```

??? success "कमांड आउटपुट"

    ```console
    Checking nf-core/demo ...
     downloaded from https://github.com/nf-core/demo.git - revision: 32893afef8 [master]
    ```

Nextflow पाइपलाइन कोड का `pull` करता है, यानी यह पूरी रिपॉजिटरी को तुम्हारी local drive पर डाउनलोड करता है।

स्पष्ट करने के लिए, तुम यह किसी भी Nextflow पाइपलाइन के साथ कर सकते हो जो GitHub में उचित रूप से सेट की गई है, न केवल nf-core पाइपलाइन के साथ।
हालांकि nf-core Nextflow पाइपलाइनों का सबसे बड़ा open-source संग्रह है।

#### 1.2.2. `nextflow list` का उपयोग करें

तुम Nextflow से उन पाइपलाइनों की सूची दे सकते हो जो तुमने इस तरह प्राप्त की हैं:

```bash
nextflow list
```

??? success "कमांड आउटपुट"

    ```console
    nf-core/demo
    ```

तुम कुछ अन्य पाइपलाइन pull करके देख सकते हो कि जब तुम्हारे पास एक से अधिक हों तो वे कैसे listed होती हैं।

#### 1.2.3. पता लगाएं कि पाइपलाइन कहां डाउनलोड हुई

तुम देखोगे कि फ़ाइलें तुम्हारी current work directory में नहीं हैं।
डिफ़ॉल्ट रूप से, Nextflow pulled पाइपलाइनों को `$NXF_HOME/assets` के अंतर्गत सहेजता है।

किसी विशिष्ट पाइपलाइन का पता लगाने के लिए, Nextflow से सीधे पूछो:

```bash
nextflow info nf-core/demo
```

??? success "कमांड आउटपुट"

    ```console
     project name: nf-core/demo
     repository  : https://github.com/nf-core/demo
     local path  : /workspaces/.nextflow/assets/.repos/nf-core/demo
     main script : main.nf
     description : An nf-core demo pipeline
     revisions   :
       TEMPLATE
       bumper
       dev
       fix-nxfversion
       manually-merge-3_0_2
     > master (default)
       nf-core-template-merge-2.13.2.dev0
       nf-core-template-merge-2.14.0
       nf-core-template-merge-2.14.1
       nf-core-template-merge-3.0.0
       nf-core-template-merge-3.0.1
       nf-core-template-merge-3.0.2
       nf-core-template-merge-3.1.0
       nf-core-template-merge-3.1.2
       nf-core-template-merge-3.2.0
       nf-core-template-merge-3.2.1
       nf-core-template-merge-3.3.1
       nf-core-template-merge-3.3.2
       nf-core-template-merge-4.0.0
       nf-core-template-merge-4.0.3
       nf-core-template-merge-4.1.0
       nf-core-template-merge-4.1.0-2
       patch
       1.0.0 [t]
       1.0.1 [t]
       1.0.2 [t]
       1.1.0 [t]
     > 1.2.0 [t]
    ```

!!! info "जानकारी"

    यदि तुम हमारे प्रशिक्षण वातावरण का उपयोग नहीं कर रहे हो तो तुम्हारे सिस्टम पर पूरा path अलग हो सकता है।

Nextflow डाउनलोड किए गए source code को जानबूझकर 'रास्ते से हटाकर' रखता है इस सिद्धांत पर कि इन पाइपलाइनों को libraries की तरह उपयोग किया जाना चाहिए न कि ऐसे कोड की तरह जिसके साथ तुम सीधे interact करोगे।

अंदर से, Nextflow प्रत्येक pulled पाइपलाइन को `$NXF_HOME/assets/.repos/` के अंतर्गत एक git रिपॉजिटरी के रूप में संग्रहीत करता है, और प्रत्येक revision के लिए कोड को `clones/<commit>/` उपडायरेक्टरी में checkout करता है।
क्योंकि `.repos` एक hidden डायरेक्टरी है, एक plain `tree -L 2 $NXF_HOME/assets/` खाली दिखेगी।

#### 1.2.4. source code तक आसानी से पहुंचने के लिए एक symlink बनाएं

हम कोड को विस्तार से नहीं देखेंगे, लेकिन आइए समग्र संगठन कैसा दिखता है इसका अंदाज़ा लगाने के लिए एक त्वरित झलक लें।

पाइपलाइन source code को browse करना आसान बनाने के लिए, pipeline के checked-out copy की ओर इशारा करने वाला एक symbolic link बनाओ:

```bash
mkdir -p pipelines/nf-core
ln -s "$(echo $NXF_HOME/assets/.repos/nf-core/demo/clones/*/)" pipelines/nf-core/demo
```

यह एक shortcut बनाता है ताकि तुम `tree -L 2 pipelines/nf-core/demo` से कोड explore कर सको या फ़ाइलें सीधे खोल सको।

#### 1.2.5. कोड संगठन का अवलोकन

तुम `tree` का उपयोग कर सकते हो या `nf-core/demo` डायरेक्टरी खोजने और खोलने के लिए file explorer का उपयोग कर सकते हो।

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

    7 directories, 12 files
    ```

जैसा कि तुम देख सकते हो, वहां बहुत कुछ चल रहा है, जिसके बारे में तुम्हें ज़्यादातर चिंता करने की ज़रूरत नहीं है।

संक्षेप में, ध्यान दें कि top level पर, तुम सारांश जानकारी के साथ एक README फ़ाइल पा सकते हो, साथ ही accessory फ़ाइलें जो licensing, contribution guidelines, citation और code of conduct जैसी प्रोजेक्ट जानकारी का सारांश देती हैं।
विस्तृत पाइपलाइन डॉक्यूमेंटेशन `docs` डायरेक्टरी में स्थित है।
यह सभी सामग्री nf-core वेबसाइट पर वेब पेज को programmatically उत्पन्न करने के लिए उपयोग की जाती है, इसलिए वे हमेशा कोड के साथ अप टू डेट रहती हैं।

बाकी के लिए, हम कोड फ़ाइलों के तीन functional समूहों को अलग कर सकते हैं:

1. पाइपलाइन कोड components (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. पाइपलाइन configuration
3. पाइपलाइन पैरामीटर / इनपुट और validation

हम इस कोर्स के इस भाग में पाइपलाइन कोड components पर नहीं जाएंगे, लेकिन हम configuration और validation के उन तत्वों को छूएंगे जो nf-core पाइपलाइनों के end user के रूप में तुम्हारे लिए प्रासंगिक होने की संभावना है।

!!! tip "सुझाव"

    तुम GitHub पर किसी भी nf-core पाइपलाइन का source code भी browse कर सकते हो, जैसे [github.com/nf-core/demo](https://github.com/nf-core/demo)।
    हर nf-core पाइपलाइन एक ही directory layout का पालन करती है, इसलिए एक बार जब तुम संरचना जान लो, तो तुम किसी भी पाइपलाइन के लिए configuration फ़ाइलें, modules और workflows उसी तरह खोज सकते हो।

अभी के लिए, पाइपलाइन चलाने की ओर बढ़ते हैं!

### सारांश

अब तुम जानते हो कि nf-core वेबसाइट के माध्यम से एक पाइपलाइन कैसे खोजें और source code की एक local copy कैसे प्राप्त करें।

### आगे क्या है?

जानो कि कम से कम प्रयास के साथ nf-core पाइपलाइन को कैसे आज़माएं।

---

## 2. test profile के साथ पाइपलाइन आज़माएं

सुविधाजनक रूप से, हर nf-core पाइपलाइन एक test profile के साथ आती है।
यह पाइपलाइन के लिए configuration settings का एक minimal सेट है जो [nf-core/test-datasets](https://github.com/nf-core/test-datasets) रिपॉजिटरी में hosted एक छोटे test dataset का उपयोग करके चलती है।
यह छोटे पैमाने पर पाइपलाइन को जल्दी आज़माने का एक शानदार तरीका है।

!!! tip "सुझाव"

    Nextflow का configuration profile सिस्टम तुम्हें विभिन्न container engines या execution environments के बीच आसानी से switch करने देता है।
    अधिक जानकारी के लिए, देखो [Hello Nextflow Part 6: Configuration](../hello_nextflow/06_hello_config.md)।

### 2.1. test profile की जांच करें

पाइपलाइन चलाने से पहले यह जांचना अच्छा अभ्यास है कि पाइपलाइन का test profile क्या निर्दिष्ट करता है।
`nf-core/demo` के लिए `test` profile configuration फ़ाइल `conf/test.config` में रहता है।
तुम इसे pipeline source के अंदर locally खोज सकते हो जो `nextflow pull` ने डाउनलोड किया, section 1.2.4 में बनाए गए `pipelines` symlink के माध्यम से:

```bash
code pipelines/nf-core/demo/conf/test.config
```

उस फ़ाइल की सामग्री यहां है:

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    न्यूनतम परीक्षण चलाने के लिए Nextflow config फ़ाइल
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    एक तेज़ और सरल पाइपलाइन परीक्षण चलाने के लिए आवश्यक इनपुट फ़ाइलें और सब कुछ परिभाषित करता है।

    इस प्रकार उपयोग करें:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

process {
    resourceLimits = [
        cpus: 2,
        memory: '4.GB',
        time: '1.h',
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // इनपुट डेटा
    input                      = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
}
```

तुम तुरंत देखोगे कि शीर्ष पर comment block में एक usage उदाहरण शामिल है जो दिखाता है कि इस test profile के साथ पाइपलाइन कैसे चलाएं।

```groovy title="conf/test.config" linenums="7"
    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

हमें केवल वही चीज़ें supply करनी हैं जो उदाहरण कमांड में carets के बीच दिखाई गई हैं: `<docker/singularity>` और `<OUTDIR>`।

याद दिलाने के लिए, `<docker/singularity>` container system की पसंद को संदर्भित करता है। सभी nf-core पाइपलाइनें reproducibility सुनिश्चित करने और software installation की समस्याओं को दूर करने के लिए containers (Docker, Singularity, आदि) के साथ उपयोग करने योग्य होने के लिए डिज़ाइन की गई हैं।
इसलिए हमें यह निर्दिष्ट करना होगा कि हम पाइपलाइन का परीक्षण करने के लिए Docker या Singularity का उपयोग करना चाहते हैं।

`--outdir <OUTDIR>` भाग उस डायरेक्टरी को संदर्भित करता है जहां Nextflow पाइपलाइन के आउटपुट लिखेगा।
हमें इसके लिए एक नाम प्रदान करना होगा, जो हम बस बना सकते हैं।
यदि यह पहले से मौजूद नहीं है, तो Nextflow runtime पर हमारे लिए इसे बना देगा।

comment block के बाद के section पर आगे बढ़ते हुए, test profile हमें दिखाता है कि परीक्षण के लिए क्या pre-configured किया गया है: सबसे उल्लेखनीय रूप से, `input` पैरामीटर पहले से ही एक test dataset की ओर इशारा करने के लिए सेट है, इसलिए हमें अपना डेटा प्रदान करने की आवश्यकता नहीं है।
यदि तुम pre-configured इनपुट के link का अनुसरण करते हो, तो तुम देखोगे कि यह एक csv फ़ाइल है जिसमें कई experimental नमूनों के लिए sample identifiers और file paths हैं।

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

इसे samplesheet कहा जाता है, और यह nf-core पाइपलाइनों के इनपुट का सबसे सामान्य रूप है।
चिंता मत करो अगर तुम data formats और types से परिचित नहीं हो, यह आगे के लिए महत्वपूर्ण नहीं है।

अब हमारे पास पाइपलाइन आज़माने के लिए सब कुछ है।

### 2.2. पाइपलाइन चलाएं

जैसा कि ऊपर उल्लेख किया गया है, हम उदाहरण testing कमांड को लगभग जैसा है वैसा उपयोग कर सकते हैं; हमें बस यह निर्दिष्ट करना है कि कौन सी software packaging का उपयोग करना है, और output डायरेक्टरी का नाम क्या रखना है।
यहां हम container system के लिए Docker और क्रमशः `demo-results` का उपयोग करेंगे।

इसके साथ, हम test कमांड चला सकते हैं:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Downloading plugin nf-schema@2.7.2
    Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_21-31-35

    Core Nextflow options
      revision                  : master
      runName                   : cranky_curry
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (8)
    [ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

यदि तुम्हारा आउटपुट उससे मेल खाता है, तो बधाई हो! तुमने अभी अपनी पहली nf-core पाइपलाइन चलाई है।

तुम देखोगे कि जब तुम एक basic Nextflow पाइपलाइन चलाते हो उससे कहीं अधिक console आउटपुट है।
एक header है जिसमें पाइपलाइन के संस्करण, इनपुट और आउटपुट का सारांश और configuration के कुछ तत्व शामिल हैं।

!!! info "जानकारी"

    तुम्हारा आउटपुट अलग timestamps, execution names और file paths दिखाएगा, लेकिन समग्र संरचना और process execution समान होनी चाहिए।

आउटपुट के शीर्ष के पास की लाइन पर ध्यान दो:

```console
Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]
```

यह तुम्हें बताता है कि पाइपलाइन का कौन सा revision उपयोग किया गया था।
क्योंकि हमने कोई संस्करण निर्दिष्ट नहीं किया, Nextflow ने `master` पर latest commit का उपयोग किया।
Reproducible runs के लिए, तुम्हें `-r` flag का उपयोग करके एक specific release pin करनी चाहिए:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir demo-results
```

यह सुनिश्चित करता है कि नए commits या releases की परवाह किए बिना हर बार एक ही पाइपलाइन कोड का उपयोग किया जाए।
इस प्रशिक्षण के लिए हम सरलता के लिए `-r` को छोड़ देते हैं, लेकिन production में तुम्हें हमेशा इसे निर्दिष्ट करना चाहिए।

execution आउटपुट पर आगे बढ़ते हुए, आइए उन लाइनों पर एक नज़र डालें जो हमें बताती हैं कि कौन से processes चलाए गए:

```console
executor >  local (8)
[ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
[b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
[ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
[09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

यह हमें बताता है कि चार processes चलाए गए, जो nf-core वेबसाइट पर पाइपलाइन डॉक्यूमेंटेशन पेज में दिखाए गए चार टूल के अनुरूप हैं: `FASTQC`, `SEQTK_TRIM`, `MULTIQC` और `COWPY`।

यहां दिखाए गए पूर्ण process names, जैसे `NFCORE_DEMO:DEMO:MULTIQC`, introductory Hello Nextflow material में तुमने जो देखा होगा उससे लंबे हैं।
इनमें उनके parent workflows के नाम शामिल हैं और पाइपलाइन कोड की modularity को दर्शाते हैं।
यदि तुम खुद nf-core-style पाइपलाइन develop करना सीखना चाहते हो, तो [Build with nf-core](../nfcore_build/index.md) कोर्स देखो।

### 2.3. पाइपलाइन के आउटपुट की जांच करें

अंत में, आइए पाइपलाइन द्वारा उत्पादित `demo-results` डायरेक्टरी पर एक नज़र डालें।

```bash
tree -L 2 demo-results
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    demo-results
    ├── cowpy
    │   └── cowpy.txt
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
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2026-07-03_21-31-35.html
        ├── execution_timeline_2026-07-03_21-31-35.html
        ├── execution_trace_2026-07-03_21-31-35.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2026-07-03_21-31-43.json
        └── pipeline_dag_2026-07-03_21-31-35.html

    12 directories, 8 files
    ```

यह बहुत कुछ लग सकता है।
`nf-core/demo` पाइपलाइन के आउटपुट के बारे में अधिक जानने के लिए, इसका [डॉक्यूमेंटेशन पेज](https://nf-co.re/demo/1.2.0/docs/output/) देखो।

इस चरण में, जो महत्वपूर्ण है वह यह है कि परिणाम module के अनुसार व्यवस्थित हैं, और इसके अतिरिक्त `pipeline_info` नामक एक डायरेक्टरी है जिसमें पाइपलाइन execution के बारे में विभिन्न timestamped reports हैं।

उदाहरण के लिए, `execution_timeline_*` फ़ाइल तुम्हें दिखाती है कि कौन से processes चलाए गए, किस क्रम में और उन्हें चलने में कितना समय लगा:

![execution timeline report](./img/execution_timeline.png)

!!! info "जानकारी"

    यहां कार्य parallel में नहीं चलाए गए क्योंकि हम Github Codespaces में एक minimalist machine पर चल रहे हैं।
    इन्हें parallel में चलते देखने के लिए, अपने codespace की CPU allocation और test configuration में resource limits बढ़ाने की कोशिश करो।

ये reports सभी nf-core पाइपलाइनों के लिए automatically उत्पन्न होती हैं।

### सारांश

तुम जानते हो कि nf-core पाइपलाइन को उसके built-in test profile का उपयोग करके कैसे चलाएं और उसके आउटपुट कहां खोजें।

### आगे क्या है?

[भाग 2](./02_configure_execution.md) पर जाओ, जहां तुम सीखोगे कि पाइपलाइन execution को कैसे configure करें।

---

## सारांश

इस भाग में तुमने सीखा:

- nf-core पाइपलाइन खोजना और प्राप्त करना और उसकी कोड संरचना की जांच करना
- पाइपलाइन को उसके built-in test profile का उपयोग करके चलाना
