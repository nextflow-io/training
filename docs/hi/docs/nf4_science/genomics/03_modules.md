# भाग 3: कोड को मॉड्यूल में ले जाना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

इस कोर्स के पहले भाग में, आपने एक variant calling पाइपलाइन बनाई जो पूरी तरह से linear थी और प्रत्येक sample के डेटा को दूसरों से स्वतंत्र रूप से प्रोसेस करती थी।

दूसरे भाग में, हमने आपको दिखाया कि GATK के साथ joint variant calling को implement करने के लिए channels और channel ऑपरेटर का उपयोग कैसे करें, जो भाग 1 की पाइपलाइन पर आधारित था।

इस भाग में, हम आपको दिखाएंगे कि उस workflow के कोड को मॉड्यूल में कैसे convert करें। इस प्रशिक्षण के इस भाग को फॉलो करने के लिए, आपको भाग 1 और भाग 2 के साथ-साथ [Hello Modules](../../../hello_nextflow/hello_modules.md) को पूरा करना चाहिए, जो मॉड्यूल की मूल बातें कवर करता है।

---

## 0. Warmup

जब हमने अपने workflow को develop करना शुरू किया, तो हमने सब कुछ एक ही कोड फ़ाइल में रखा।
अब हमारे कोड को **modularize** करने का समय आ गया है, _यानी_ process definitions को मॉड्यूल में extract करना।

हम भाग 2 के समान workflow से शुरू करने जा रहे हैं, जिसे हमने आपके लिए `genomics-3.nf` फ़ाइल में प्रदान किया है।

!!! note "नोट"

     सुनिश्चित करें कि आप सही working डायरेक्टरी में हैं:
     `cd /workspaces/training/nf4-science/genomics`

शुरुआती बिंदु को verify करने के लिए workflow को चलाएं:

```bash
nextflow run genomics-3.nf -resume
```

```console title="आउटपुट"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [serene_borg] DSL2 - revision: 0cbebb67a1

executor >  local (7)
[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (1) | 3 of 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

अब आपकी project डायरेक्टरी के अंदर एक `work` डायरेक्टरी और एक `results_genomics` डायरेक्टरी होगी।

### सारांश

आप अपने workflow को modularize करना शुरू करने के लिए तैयार हैं।

### आगे क्या है?

Genomics workflow की processes को मॉड्यूल में ले जाएं।

---

## 1. Processes को मॉड्यूल में ले जाएं

जैसा कि आपने [Hello Modules](../../../hello_nextflow/hello_modules.md) में सीखा, आप केवल process definition को अपनी खुद की फ़ाइल में copy करके, किसी भी डायरेक्टरी में, एक मॉड्यूल बना सकते हैं, और आप उस फ़ाइल को कुछ भी नाम दे सकते हैं।

उन कारणों से जो बाद में स्पष्ट हो जाएंगे (विशेष रूप से जब हम testing पर आएंगे), इस प्रशिक्षण में हम फ़ाइल को `main.nf` नाम देने की convention का पालन करेंगे, और इसे tool kit और कमांड के नाम पर आधारित डायरेक्टरी structure में रखेंगे।

### 1.1. `SAMTOOLS_INDEX` process के लिए एक मॉड्यूल बनाएं

`SAMTOOLS_INDEX` process के मामले में, 'samtools' toolkit है और 'index' कमांड है। इसलिए, हम एक डायरेक्टरी structure `modules/samtools/index` बनाएंगे और उस डायरेक्टरी के अंदर `main.nf` फ़ाइल में `SAMTOOLS_INDEX` process definition रखेंगे।

```bash
mkdir -p modules/samtools/index
touch modules/samtools/index/main.nf
```

`main.nf` फ़ाइल खोलें और इसमें `SAMTOOLS_INDEX` process definition को copy करें।

```groovy title="modules/samtools/index/main.nf" linenums="1"
/*
 * BAM इंडेक्स फ़ाइल बनाएं
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    tuple path(input_bam), path("${input_bam}.bai")

    script:
    """
    samtools index '$input_bam'
    """
}
```

फिर, `genomics-3.nf` से `SAMTOOLS_INDEX` process definition को हटाएं, और अगली process definition से पहले मॉड्यूल के लिए एक import declaration जोड़ें, इस तरह:

=== "बाद में"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1 2"
    // मॉड्यूल include करें
    include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'

    /*
     * GATK HaplotypeCaller के साथ variants को call करें
     */
    process GATK_HAPLOTYPECALLER {
    ```

=== "पहले"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1"
    /*
     * GATK HaplotypeCaller के साथ variants को call करें
     */
    process GATK_HAPLOTYPECALLER {
    ```

अब आप workflow को फिर से चला सकते हैं, और यह पहले की तरह ही काम करना चाहिए। यदि आप `-resume` flag प्रदान करते हैं, तो कोई नया कार्य भी चलाने की आवश्यकता नहीं होनी चाहिए:

```bash
nextflow run genomics-3.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-3.nf` [sleepy_snyder] DSL2 - revision: aa68d06c43

    [0f/71b55e] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [f1/18971b] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    [0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
    ```

### 1.2. `GATK_HAPLOTYPECALLER` और `GATK_JOINTGENOTYPING` processes के लिए मॉड्यूल बनाएं

शेष processes के लिए वही steps दोहराएं।
प्रत्येक process के लिए:

1. डायरेक्टरी structure बनाएं (`modules/gatk/haplotypecaller/` और `modules/gatk/jointgenotyping/`)
2. Process definition वाली एक `main.nf` फ़ाइल बनाएं
3. `genomics-3.nf` से process definition हटाएं
4. मॉड्यूल के लिए एक import declaration जोड़ें

एक बार जब आप काम कर लें, तो यह चलाकर जांचें कि आपकी modules डायरेक्टरी structure सही है:

```bash
tree modules/
```

??? abstract "डायरेक्टरी की सामग्री"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf

    5 directories, 3 files
    ```

पैरामीटर section के बाद, main workflow फ़ाइल में आपके पास इस तरह कुछ होना चाहिए:

```
include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'

workflow {
```

### सारांश

आपने genomics workflow को उदाहरण के रूप में लेकर एक workflow को modularize करने का अभ्यास किया है।

### आगे क्या है?

Modularized workflow को test करें।

---

## 2. Modularized workflow को test करें

सब कुछ अभी भी काम करता है यह verify करने के लिए modularized workflow को चलाएं।

```bash
nextflow run genomics-3.nf -resume
```

```console title="आउटपुट"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [astonishing_venter] DSL2 - revision: ca27264c13

[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
```

सब कुछ अभी भी काम करता है, जिसमें पाइपलाइन की resumability भी शामिल है।
परिणाम `results_genomics` डायरेक्टरी में publish होते रहते हैं।

```console title="डायरेक्टरी की सामग्री"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
```

### सारांश

आपने एक workflow को modularize किया है और verify किया है कि यह अभी भी पहले की तरह ही काम करता है।

### आगे क्या है?

आपने जो सीखा है उसकी समीक्षा करें और testing की ओर देखें।

---

## 3. सारांश

आपने workflow को modularize किया है, और पाइपलाइन के काम करने के तरीके में कुछ भी नहीं बदला है।
यह जानबूझकर किया गया है: आपने कोड को restructure किया है बिना इसके function को प्रभावित किए।

मॉड्यूल में केवल process logic होता है, जो उन्हें साफ और reusable बनाता है।
main script नियंत्रित करती है कि क्या publish होता है और कहां, जबकि मॉड्यूल अपने computational कार्य पर केंद्रित रहते हैं।

आपने उन चीजों के लिए एक नींव रखी है जो आपके कोड को maintain करना आसान बना देंगी।
उदाहरण के लिए, अब आप nf-test framework का उपयोग करके अपनी पाइपलाइन में tests जोड़ सकते हैं।
यही वह है जो हम इस कोर्स के अगले भाग में देखेंगे।
