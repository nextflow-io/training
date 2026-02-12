# भाग 3: बहु-नमूना युग्मित-अंत कार्यान्वयन

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

पहले, तुमने प्रति-नमूना वेरिएंट कॉलिंग पाइपलाइन बनाई थी जो हर नमूने के डेटा को स्वतंत्र रूप से प्रोसेस करती थी।
इस पाठ्यक्रम के इस भाग में, हम अपने सरल workflow को अगले स्तर पर ले जाने वाले हैं और इसे एक शक्तिशाली बैच ऑटोमेशन टूल में बदलने वाले हैं ताकि यह मनमाने संख्या के नमूनों को संभाल सके।
और जब हम यह कर रहे हैं, हम इसे युग्मित-अंत डेटा की अपेक्षा करने के लिए भी अपडेट करने वाले हैं, जो नए अध्ययनों में अधिक सामान्य है।

??? info "इस खंड से कैसे शुरू करें"

    पाठ्यक्रम का यह खंड मानता है कि तुमने [भाग 1: विधि अवलोकन](./01_method.md), [भाग 2: एकल-नमूना कार्यान्वयन](./02_single-sample.md) पूरा कर लिया है और तुम्हारे पास भरी हुई मॉड्यूल फ़ाइलों के साथ एक कार्यशील `rnaseq.nf` पाइपलाइन है।

    यदि तुमने भाग 2 पूरा नहीं किया या इस भाग के लिए नए सिरे से शुरू करना चाहते हो, तो तुम भाग 2 के समाधान को अपने शुरुआती बिंदु के रूप में उपयोग कर सकते हो।
    `nf4-science/rnaseq/` डायरेक्टरी के अंदर से ये कमांड चलाओ:

    ```bash
    cp solutions/part2/rnaseq-2.nf rnaseq.nf
    cp solutions/part2/modules/fastqc.nf modules/
    cp solutions/part2/modules/trim_galore.nf modules/
    cp solutions/part2/modules/hisat2_align.nf modules/
    cp solutions/part2/nextflow.config .
    ```

    यह तुम्हें एक पूर्ण एकल-नमूना प्रोसेसिंग workflow देता है।
    तुम परीक्षण कर सकते हो कि यह सफलतापूर्वक चलता है:

    ```bash
    nextflow run rnaseq.nf -profile test
    ```

## असाइनमेंट

पाठ्यक्रम के इस भाग में, हम workflow को निम्नलिखित करने के लिए विस्तारित करने वाले हैं:

1. CSV samplesheet से नमूना जानकारी पढ़ें
2. सभी नमूनों पर प्रति-नमूना QC, ट्रिमिंग और संरेखण समानांतर में चलाएं
3. सभी QC रिपोर्ट को एक व्यापक MultiQC रिपोर्ट में एकत्रित करें

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-03.svg"
</figure>

यह [भाग 1: विधि अवलोकन](./01_method.md#2-multi-sample-qc-aggregation) के दूसरे खंड से चरणों को स्वचालित करता है, जहां तुमने इन कमांड को उनके कंटेनरों में मैन्युअल रूप से चलाया था।

## पाठ योजना

हमने इसे तीन चरणों में विभाजित किया है:

1. **workflow को कई इनपुट नमूने स्वीकार करने योग्य बनाएं।**
   यह एकल फ़ाइल पथ से CSV samplesheet में स्विच करना, इसे `splitCsv()` से पार्स करना और सभी मौजूदा processes को कई नमूनों पर चलाना कवर करता है।
2. **व्यापक QC रिपोर्ट जनरेशन जोड़ें।**
   यह नमूनों में आउटपुट को एकत्रित करने के लिए `collect()` ऑपरेटर पेश करता है, और एक संयुक्त रिपोर्ट उत्पन्न करने के लिए MultiQC process जोड़ता है।
3. **युग्मित-अंत RNAseq डेटा पर स्विच करें।**
   यह युग्मित-अंत इनपुट (टपल का उपयोग करके) के लिए processes को अनुकूलित करना, युग्मित-अंत मॉड्यूल बनाना और एक अलग टेस्ट प्रोफ़ाइल सेट करना कवर करता है।

यह [भाग 1: विधि अवलोकन](./01_method.md) में वर्णित विधि को लागू करता है (बहु-नमूना उपयोग मामले को कवर करने वाला दूसरा खंड) और भाग 2 द्वारा उत्पादित workflow पर सीधे निर्माण करता है।

!!! tip "सुझाव"

     सुनिश्चित करो कि तुम सही कार्य डायरेक्टरी में हो:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. workflow को कई इनपुट नमूने स्वीकार करने योग्य बनाएं

कई नमूनों पर चलाने के लिए, हमें इनपुट को प्रबंधित करने के तरीके को बदलना होगा: एकल फ़ाइल पथ प्रदान करने के बजाय, हम CSV फ़ाइल से नमूना जानकारी पढ़ेंगे।

हम `data/` डायरेक्टरी में नमूना IDs और FASTQ फ़ाइल पथों वाली एक CSV फ़ाइल प्रदान करते हैं।

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

यह CSV फ़ाइल एक हेडर लाइन शामिल करती है जो कॉलम को नाम देती है।

ध्यान दें कि यह अभी भी single-end रीड डेटा है।

!!! warning "चेतावनी"

    CSV में फ़ाइल पथ पूर्ण पथ हैं जो तुम्हारे वातावरण से मेल खाने चाहिए।
    यदि तुम इसे हमारे द्वारा प्रदान किए गए प्रशिक्षण वातावरण में नहीं चला रहे हो, तो तुम्हें अपने सिस्टम से मेल खाने के लिए पथों को अपडेट करना होगा।

### 1.1. टेस्ट प्रोफ़ाइल में प्राथमिक इनपुट को फ़ाइल पथों के CSV में बदलें

सबसे पहले, हमें एकल FASTQ पथ के बजाय CSV फ़ाइल पथ प्रदान करने के लिए `nextflow.config` में टेस्ट प्रोफ़ाइल को अपडेट करना होगा।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

इसके बाद, हमें इस CSV से पढ़ने के लिए channel निर्माण को अपडेट करना होगा।

### 1.2. CSV इनपुट को पार्स करने के लिए channel फ़ैक्टरी को अपडेट करें

हमें केवल फ़ाइल पथ के बजाय फ़ाइल की सामग्री को channel में लोड करना होगा।

हम वही पैटर्न उपयोग कर सकते हैं जो हमने [Hello Nextflow के भाग 2](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file) में उपयोग किया था: फ़ाइल को पार्स करने के लिए [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) ऑपरेटर लागू करना, फिर हर पंक्ति से FASTQ फ़ाइल पथ निकालने के लिए `map` ऑपरेशन।

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="4-7"
    workflow {

        main:
        // CSV फ़ाइल की सामग्री से इनपुट channel बनाएं
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // फ़ाइल पथ से इनपुट channel बनाएं
        read_ch = channel.fromPath(params.input)
    ```

Hello Nextflow पाठ्यक्रम में जो तुमने देखा उसकी तुलना में एक नई बात यह है कि इस CSV में एक हेडर लाइन है, इसलिए हम `splitCsv()` कॉल में `#!groovy header: true` जोड़ते हैं।
यह हमें `map` ऑपरेशन में कॉलम को नाम से संदर्भित करने की अनुमति देता है: `#!groovy row.fastq_path` हर पंक्ति के `fastq_path` कॉलम से फ़ाइल पथ निकालता है।

इनपुट हैंडलिंग अपडेट हो गई है और workflow परीक्षण के लिए तैयार है।

### 1.3. workflow को चलाएं

workflow अब CSV फ़ाइल से नमूना जानकारी पढ़ता है और सभी नमूनों को समानांतर में प्रोसेस करता है।

```bash
nextflow run rnaseq.nf -profile test
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [golden_curry] DSL2 - revision: 2a5ba5be1e

    executor >  local (18)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6 ✔
    [cc/16859f] TRIM_GALORE (6)  [100%] 6 of 6 ✔
    [68/4c27b5] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
    ```

इस बार हर चरण 6 बार चलता है, CSV फ़ाइल में हर नमूने के लिए एक बार।

workflow को कई फ़ाइलों पर चलाने के लिए बस इतना ही काफी था।
Nextflow हमारे लिए सभी समानांतरता को संभालता है।

### मुख्य बातें

तुम जानते हो कि एकल-फ़ाइल इनपुट से CSV-आधारित बहु-नमूना इनपुट में कैसे स्विच करें जिसे Nextflow समानांतर में प्रोसेस करता है।

### आगे क्या है?

एक QC रिपोर्ट एकत्रीकरण चरण जोड़ें जो सभी नमूनों से मेट्रिक्स को संयोजित करता है।

---

## 2. प्री-प्रोसेसिंग QC मेट्रिक्स को एकल MultiQC रिपोर्ट में एकत्रित करें

यह सब बहुत सारी QC रिपोर्ट उत्पन्न करता है, और हम व्यक्तिगत रिपोर्टों को खोदना नहीं चाहते।
यह MultiQC रिपोर्ट एकत्रीकरण चरण लगाने के लिए एकदम सही बिंदु है।

[भाग 1](01_method.md) से `multiqc` कमांड याद करो:

```bash
multiqc . -n <output_name>.html
```

कमांड वर्तमान डायरेक्टरी को मान्यता प्राप्त QC आउटपुट फ़ाइलों के लिए स्कैन करती है और उन्हें एक एकल HTML रिपोर्ट में एकत्रित करती है।
कंटेनर URI `community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c` था।

हमें एक अतिरिक्त पैरामीटर सेट करना होगा, इनपुट तैयार करना होगा, process लिखना होगा, इसे वायर करना होगा और आउटपुट हैंडलिंग को अपडेट करना होगा।

### 2.1. इनपुट सेट करें

MultiQC process को एक रिपोर्ट नाम पैरामीटर और पिछले चरणों से सभी QC-संबंधित आउटपुट एक साथ बंडल किए हुए चाहिए।

#### 2.1.1. एक `report_id` पैरामीटर जोड़ें

आउटपुट रिपोर्ट को नाम देने के लिए एक पैरामीटर जोड़ें।

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="8-9"
    params {
        // प्राथमिक इनपुट
        input: Path

        // संदर्भ जीनोम आर्काइव
        hisat2_index_zip: Path

        // रिपोर्ट ID
        report_id: String
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // प्राथमिक इनपुट
        input: Path

        // संदर्भ जीनोम आर्काइव
        hisat2_index_zip: Path
    }
    ```

टेस्ट प्रोफ़ाइल में रिपोर्ट ID डिफ़ॉल्ट जोड़ें:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

इसके बाद, हमें MultiQC process के लिए इनपुट तैयार करने होंगे।

#### 2.1.2. पिछले चरणों से QC आउटपुट एकत्र और संयोजित करें

हमें `MULTIQC` process को पिछले चरणों से सभी QC-संबंधित आउटपुट एक साथ बंडल करके देने की आवश्यकता है।

इसके लिए, हम `.mix()` ऑपरेटर का उपयोग करते हैं, जो कई channels को एक में एकत्रित करता है।
हम `channel.empty()` से शुरू करते हैं और सभी आउटपुट channels को मिलाते हैं जिन्हें हम संयोजित करना चाहते हैं।
यह सीधे एक आउटपुट channel पर `.mix()` चेन करने से अधिक स्वच्छ है, क्योंकि यह सभी इनपुट को समान रूप से व्यवहार करता है।

हमारे workflow में, एकत्रित करने के लिए QC-संबंधित आउटपुट हैं:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

हम उन्हें एक एकल channel में मिलाते हैं, फिर सभी नमूनों में रिपोर्ट को एक एकल सूची में एकत्रित करने के लिए `.collect()` का उपयोग करते हैं।

`HISAT2_ALIGN` कॉल के बाद workflow बॉडी में ये लाइनें जोड़ें:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="38" hl_lines="4-13"
        // संदर्भ जीनोम के साथ संरेखण
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

        // व्यापक QC रिपोर्ट जनरेशन
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
        multiqc_files_list = multiqc_files_ch.collect()
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="38"
        // संदर्भ जीनोम के साथ संरेखण
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

मध्यवर्ती वेरिएबल का उपयोग हर चरण को स्पष्ट बनाता है: `multiqc_files_ch` में सभी व्यक्तिगत QC फ़ाइलें एक channel में मिली हुई हैं, और `multiqc_files_list` MultiQC को पास करने के लिए तैयार एकत्रित बंडल है।

### 2.2. QC एकत्रीकरण process लिखें और इसे workflow में कॉल करें

पहले की तरह, हमें process परिभाषा भरनी होगी, मॉड्यूल आयात करना होगा और process कॉल जोड़ना होगा।

#### 2.2.1. QC एकत्रीकरण process के लिए मॉड्यूल भरें

`modules/multiqc.nf` खोलें और process परिभाषा की रूपरेखा की जांच करें।

आगे बढ़ो और ऊपर दी गई जानकारी का उपयोग करके process परिभाषा को स्वयं भरो, फिर नीचे "बाद में" टैब में समाधान के विरुद्ध अपने काम की जांच करो।

=== "पहले"

    ```groovy title="modules/multiqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * MultiQC के साथ QC रिपोर्ट एकत्रित करें
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "बाद में"

    ```groovy title="modules/multiqc.nf" linenums="1" hl_lines="8 11 12 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * MultiQC के साथ QC रिपोर्ट एकत्रित करें
     */
    process MULTIQC {

        container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"

        input:
        path '*'
        val output_name

        output:
        path "${output_name}.html", emit: report
        path "${output_name}_data", emit: data

        script:
        """
        multiqc . -n ${output_name}.html
        """
    }
    ```

यह process QC फ़ाइलों के लिए इनपुट क्वालिफ़ायर के रूप में `#!groovy path '*'` का उपयोग करती है।
`'*'` वाइल्डकार्ड Nextflow को विशिष्ट नामों की आवश्यकता के बिना सभी एकत्रित फ़ाइलों को कार्य डायरेक्टरी में स्टेज करने के लिए कहता है।
`val output_name` इनपुट एक स्ट्रिंग है जो रिपोर्ट फ़ाइलनाम को नियंत्रित करती है।

`multiqc .` कमांड वर्तमान डायरेक्टरी (जहां सभी स्टेज की गई QC फ़ाइलें हैं) को स्कैन करती है और रिपोर्ट जनरेट करती है।

एक बार जब तुमने यह पूरा कर लिया, तो process उपयोग के लिए तैयार है।

#### 2.2.2. मॉड्यूल शामिल करें

`rnaseq.nf` में आयात स्टेटमेंट जोड़ें:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="5"
    // मॉड्यूल INCLUDE स्टेटमेंट
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    include { MULTIQC } from './modules/multiqc.nf'
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="3"
    // मॉड्यूल INCLUDE स्टेटमेंट
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

अब workflow में process कॉल जोड़ें।

#### 2.2.3. process कॉल जोड़ें

एकत्रित QC फ़ाइलों और रिपोर्ट ID को `MULTIQC` process को पास करें:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="49" hl_lines="2"
        multiqc_files_list = multiqc_files_ch.collect()
        MULTIQC(multiqc_files_list, params.report_id)
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="49"
        multiqc_files_list = multiqc_files_ch.collect()
    ```

MultiQC process अब workflow में वायर हो गई है।

### 2.3. आउटपुट हैंडलिंग को अपडेट करें

हमें MultiQC आउटपुट को publish घोषणा में जोड़ना होगा और कॉन्फ़िगर करना होगा कि वे कहां जाएं।

#### 2.3.1. MultiQC आउटपुट के लिए publish टारगेट जोड़ें

`publish:` सेक्शन में MultiQC आउटपुट जोड़ें:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="52" hl_lines="9-10"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

इसके बाद, हमें Nextflow को बताना होगा कि इन आउटपुट को कहां रखना है।

#### 2.3.2. नए आउटपुट टारगेट कॉन्फ़िगर करें

`output {}` ब्लॉक में MultiQC टारगेट के लिए एंट्री जोड़ें, उन्हें `multiqc/` सबडायरेक्टरी में प्रकाशित करते हुए:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="82" hl_lines="4-9"
        align_log {
            path 'align'
        }
        multiqc_report {
            path 'multiqc'
        }
        multiqc_data {
            path 'multiqc'
        }
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="82"
        align_log {
            path 'align'
        }
    }
    ```

आउटपुट कॉन्फ़िगरेशन पूर्ण है।

### 2.4. workflow को चलाएं

हम `-resume` का उपयोग करते हैं ताकि पिछले प्रोसेसिंग चरण कैश हो जाएं और केवल नया MultiQC चरण चले।

```bash
nextflow run rnaseq.nf -profile test -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [modest_pare] DSL2 - revision: fc724d3b49

    executor >  local (1)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6, cached: 6 ✔
    [2c/8d8e1e] TRIM_GALORE (5)  [100%] 6 of 6, cached: 6 ✔
    [a4/7f9c44] HISAT2_ALIGN (6) [100%] 6 of 6, cached: 6 ✔
    [56/e1f102] MULTIQC          [100%] 1 of 1 ✔
    ```

कैश किए गए process कॉलों के बाद MULTIQC के लिए एक एकल कॉल जोड़ा गया है।

तुम results डायरेक्टरी में MultiQC आउटपुट पा सकते हो।

```bash
tree -L 2 results/multiqc
```

```console title="आउटपुट"
results/multiqc
├── all_single-end_data
│   ├── cutadapt_filtered_reads_plot.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Counts.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt
│   ├── fastqc_adapter_content_plot.txt
│   ├── fastqc_overrepresented_sequences_plot.txt
│   ├── fastqc_per_base_n_content_plot.txt
│   ├── fastqc_per_base_sequence_quality_plot.txt
│   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   ├── fastqc_per_sequence_quality_scores_plot.txt
│   ├── fastqc_sequence_counts_plot.txt
│   ├── fastqc_sequence_duplication_levels_plot.txt
│   ├── fastqc_sequence_length_distribution_plot.txt
│   ├── fastqc-status-check-heatmap.txt
│   ├── fastqc_top_overrepresented_sequences_table.txt
│   ├── hisat2_se_plot.txt
│   ├── multiqc_citations.txt
│   ├── multiqc_cutadapt.txt
│   ├── multiqc_data.json
│   ├── multiqc_fastqc.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc_hisat2.txt
│   ├── multiqc.log
│   ├── multiqc_software_versions.txt
│   └── multiqc_sources.txt
└── all_single-end.html
```

वह अंतिम `all_single-end.html` फ़ाइल पूर्ण एकत्रित रिपोर्ट है, जो एक आसान ब्राउज़ करने योग्य HTML फ़ाइल में सुविधाजनक रूप से पैक की गई है।

### मुख्य बातें

तुम जानते हो कि कई channels से आउटपुट कैसे एकत्र करें, उन्हें `.mix()` और `.collect()` से कैसे बंडल करें, और उन्हें एकत्रीकरण process को कैसे पास करें।

### आगे क्या है?

युग्मित-अंत RNAseq डेटा को संभालने के लिए workflow को अनुकूलित करें।

---

## 3. युग्मित-अंत RNAseq डेटा प्रोसेसिंग सक्षम करें

अभी हमारा workflow केवल single-end RNAseq डेटा को संभाल सकता है।
युग्मित-अंत RNAseq डेटा देखना तेजी से सामान्य हो रहा है, इसलिए हम इसे संभालने में सक्षम होना चाहते हैं।

workflow को डेटा प्रकार से पूरी तरह अज्ञेयवादी बनाने के लिए थोड़ा अधिक उन्नत Nextflow भाषा सुविधाओं का उपयोग करने की आवश्यकता होगी, इसलिए हम यहां ऐसा नहीं करने वाले हैं, लेकिन हम यह प्रदर्शित करने के लिए एक युग्मित-अंत प्रोसेसिंग संस्करण बना सकते हैं कि क्या अनुकूलित करने की आवश्यकता है।

### 3.1. workflow की एक प्रति बनाएं और इनपुट अपडेट करें

हम workflow फ़ाइल की एक प्रति बनाकर शुरू करते हैं और इसे युग्मित-अंत डेटा के लिए अपडेट करते हैं।

#### 3.1.1. workflow फ़ाइल की प्रति बनाएं

युग्मित-अंत संस्करण के लिए शुरुआती बिंदु के रूप में उपयोग करने के लिए workflow फ़ाइल की एक प्रति बनाएं।

```bash
cp rnaseq.nf rnaseq_pe.nf
```

अब नई फ़ाइल में पैरामीटर और इनपुट हैंडलिंग को अपडेट करें।

#### 3.1.2. एक युग्मित-अंत टेस्ट प्रोफ़ाइल जोड़ें

हम `data/` डायरेक्टरी में नमूना IDs और युग्मित FASTQ फ़ाइल पथों वाली दूसरी CSV फ़ाइल प्रदान करते हैं।

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

`nextflow.config` में एक `test_pe` प्रोफ़ाइल जोड़ें जो इस फ़ाइल की ओर इंगित करती है और युग्मित-अंत रिपोर्ट ID का उपयोग करती है।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="9-13"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
        test_pe {
            params.input = "${projectDir}/data/paired-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_paired-end"
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

युग्मित-अंत डेटा के लिए टेस्ट प्रोफ़ाइल तैयार है।

#### 3.1.3. channel फ़ैक्टरी को अपडेट करें

`.map()` ऑपरेटर को दोनों FASTQ फ़ाइल पथ पकड़ने और उन्हें एक सूची के रूप में वापस करने की आवश्यकता है।

=== "बाद में"

    ```groovy title="rnaseq_pe.nf" linenums="25" hl_lines="4"
        // CSV फ़ाइल की सामग्री से इनपुट channel बनाएं
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
    ```

=== "पहले"

    ```groovy title="rnaseq_pe.nf" linenums="25"
        // CSV फ़ाइल की सामग्री से इनपुट channel बनाएं
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

इनपुट हैंडलिंग युग्मित-अंत डेटा के लिए कॉन्फ़िगर है।

### 3.2. युग्मित-अंत डेटा के लिए FASTQC मॉड्यूल को अनुकूलित करें

युग्मित-अंत संस्करण बनाने के लिए मॉड्यूल की प्रति बनाएं:

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

FASTQC process इनपुट को बदलने की आवश्यकता नहीं है — जब Nextflow को दो फ़ाइलों की सूची मिलती है, तो यह दोनों को स्टेज करता है और `reads` दोनों फ़ाइलनामों में विस्तारित होता है।
केवल आवश्यक परिवर्तन आउटपुट ब्लॉक में है: चूंकि हमें अब प्रति नमूना दो FastQC रिपोर्ट मिलती हैं, हम `simpleName`-आधारित पैटर्न से वाइल्डकार्ड में स्विच करते हैं।

=== "बाद में"

    ```groovy title="modules/fastqc_pe.nf" linenums="10" hl_lines="2 3"
        output:
        path "*_fastqc.zip", emit: zip
        path "*_fastqc.html", emit: html
    ```

=== "पहले"

    ```groovy title="modules/fastqc_pe.nf" linenums="10"
        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html
    ```

यह process को इस तरह से सामान्यीकृत करता है कि यह single-end या युग्मित-अंत डेटा को संभालने में सक्षम हो जाता है।

युग्मित-अंत संस्करण का उपयोग करने के लिए `rnaseq_pe.nf` में आयात को अपडेट करें:

=== "बाद में"

    ```groovy title="rnaseq_pe.nf" linenums="4" hl_lines="1"
    include { FASTQC } from './modules/fastqc_pe.nf'
    ```

=== "पहले"

    ```groovy title="rnaseq_pe.nf" linenums="4"
    include { FASTQC } from './modules/fastqc.nf'
    ```

FASTQC मॉड्यूल और इसका आयात युग्मित-अंत डेटा के लिए अपडेट हो गया है।

### 3.3. युग्मित-अंत डेटा के लिए TRIM_GALORE मॉड्यूल को अनुकूलित करें

युग्मित-अंत संस्करण बनाने के लिए मॉड्यूल की प्रति बनाएं:

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

इस मॉड्यूल को अधिक महत्वपूर्ण परिवर्तनों की आवश्यकता है:

- इनपुट एकल पथ से दो पथों के टपल में बदलता है
- कमांड `--paired` फ्लैग जोड़ती है और दोनों रीड फ़ाइलें लेती है
- आउटपुट Trim Galore के युग्मित-अंत नामकरण परंपराओं को प्रतिबिंबित करने के लिए बदलता है, हर रीड फ़ाइल के लिए अलग FastQC रिपोर्ट उत्पन्न करता है

=== "बाद में"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8" hl_lines="2 5 7 8 12"
        input:
        tuple path(read1), path(read2)

        output:
        tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
        path "*_trimming_report.txt", emit: trimming_reports
        path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
        path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

        script:
        """
        trim_galore --fastqc --paired ${read1} ${read2}
        """
    ```

=== "पहले"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8"
        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    ```

`rnaseq_pe.nf` में आयात को अपडेट करें:

=== "बाद में"

    ```groovy title="rnaseq_pe.nf" linenums="5" hl_lines="1"
    include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
    ```

=== "पहले"

    ```groovy title="rnaseq_pe.nf" linenums="5"
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

TRIM_GALORE मॉड्यूल और इसका आयात युग्मित-अंत डेटा के लिए अपडेट हो गया है।

### 3.4. युग्मित-अंत डेटा के लिए HISAT2_ALIGN मॉड्यूल को अनुकूलित करें

युग्मित-अंत संस्करण बनाने के लिए मॉड्यूल की प्रति बनाएं:

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

इस मॉड्यूल को समान परिवर्तनों की आवश्यकता है:

- इनपुट एकल पथ से दो पथों के टपल में बदलता है
- HISAT2 कमांड `-U` (unpaired) से `-1` और `-2` (paired) रीड तर्कों में बदलती है
- `reads.simpleName` के सभी उपयोग `read1.simpleName` में बदलते हैं क्योंकि हम अब जोड़ी के एक विशिष्ट सदस्य को संदर्भित करते हैं

=== "बाद में"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8" hl_lines="2 6 7 12 13 14"
        input:
        tuple path(read1), path(read2)
        path index_zip

        output:
        path "${read1.simpleName}.bam", emit: bam
        path "${read1.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
            --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
            samtools view -bS -o ${read1.simpleName}.bam
        """
    ```

=== "पहले"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8"
        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    ```

`rnaseq_pe.nf` में आयात को अपडेट करें:

=== "बाद में"

    ```groovy title="rnaseq_pe.nf" linenums="6" hl_lines="1"
    include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
    ```

=== "पहले"

    ```groovy title="rnaseq_pe.nf" linenums="6"
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

HISAT2_ALIGN मॉड्यूल और इसका आयात युग्मित-अंत डेटा के लिए अपडेट हो गया है।

### 3.5. युग्मित-अंत आउटपुट के लिए MultiQC एकत्रीकरण को अपडेट करें

युग्मित-अंत `TRIM_GALORE` process अब एक के बजाय दो अलग FastQC रिपोर्ट channels (`fastqc_reports_1` और `fastqc_reports_2`) उत्पन्न करती है।
दोनों को शामिल करने के लिए `rnaseq_pe.nf` में `.mix()` ब्लॉक को अपडेट करें:

=== "बाद में"

    ```groovy title="rnaseq_pe.nf" linenums="40" hl_lines="5 6"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            HISAT2_ALIGN.out.log,
        )
    ```

=== "पहले"

    ```groovy title="rnaseq_pe.nf" linenums="40"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
    ```

MultiQC एकत्रीकरण अब युग्मित-अंत FastQC रिपोर्ट के दोनों सेट शामिल करता है।

### 3.6. युग्मित-अंत आउटपुट के लिए आउटपुट हैंडलिंग को अपडेट करें

`publish:` सेक्शन और `output {}` ब्लॉक को भी युग्मित-अंत `TRIM_GALORE` process से दो अलग FastQC रिपोर्ट channels को प्रतिबिंबित करने की आवश्यकता है।

`rnaseq_pe.nf` में `publish:` सेक्शन को अपडेट करें:

=== "बाद में"

    ```groovy title="rnaseq_pe.nf" linenums="52" hl_lines="6-7"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc_1 = TRIM_GALORE.out.fastqc_reports_1
        trimming_fastqc_2 = TRIM_GALORE.out.fastqc_reports_2
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "पहले"

    ```groovy title="rnaseq_pe.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

`output {}` ब्लॉक में संबंधित एंट्री को अपडेट करें:

=== "बाद में"

    ```groovy title="rnaseq_pe.nf" linenums="77" hl_lines="4-9"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc_1 {
            path 'trimming'
        }
        trimming_fastqc_2 {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

=== "पहले"

    ```groovy title="rnaseq_pe.nf" linenums="77"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

युग्मित-अंत workflow अब पूरी तरह से अपडेट है और चलाने के लिए तैयार है।

### 3.7. workflow को चलाएं

हम `-resume` का उपयोग नहीं करते क्योंकि यह कैश नहीं होगा, और पहले की तुलना में दोगुना डेटा प्रोसेस करना है, लेकिन फिर भी यह एक मिनट से कम में पूरा होना चाहिए।

```bash
nextflow run rnaseq_pe.nf -profile test_pe
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq_pe.nf` [reverent_kare] DSL2 - revision: 9c376cc219

    executor >  local (19)
    [c5/cbde15] FASTQC (5)       [100%] 6 of 6 ✔
    [e4/fa2784] TRIM_GALORE (5)  [100%] 6 of 6 ✔
    [3a/e23049] HISAT2_ALIGN (5) [100%] 6 of 6 ✔
    [e6/a3ccd9] MULTIQC          [100%] 1 of 1 ✔
    ```

अब हमारे पास हमारे workflow के दो थोड़े भिन्न संस्करण हैं, एक single-end रीड डेटा के लिए और एक युग्मित-अंत डेटा के लिए।
अगला तार्किक कदम workflow को किसी भी डेटा प्रकार को तुरंत स्वीकार करने योग्य बनाना होगा, जो इस पाठ्यक्रम के दायरे से बाहर है, लेकिन हम इसे एक फॉलो-अप में संभाल सकते हैं।

---

### मुख्य बातें

तुम जानते हो कि एकल-नमूना workflow को कई नमूनों की प्रोसेसिंग को समानांतरित करने के लिए कैसे अनुकूलित करें, एक व्यापक QC रिपोर्ट उत्पन्न करें और युग्मित-अंत रीड डेटा का उपयोग करने के लिए workflow को कैसे अनुकूलित करें।

### आगे क्या है?

अपनी पीठ थपथपाओ! तुमने Nextflow For RNAseq पाठ्यक्रम पूरा कर लिया है।

अंतिम [पाठ्यक्रम सारांश](./next_steps.md) पर जाओ ताकि तुमने क्या सीखा उसकी समीक्षा करो और पता लगाओ कि आगे क्या आता है।
