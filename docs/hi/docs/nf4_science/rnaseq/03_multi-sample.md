# भाग 3: बहु-नमूना युग्मित-अंत कार्यान्वयन

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

इस पाठ्यक्रम के इस अंतिम भाग में, हम अपने सरल workflow को अगले स्तर पर ले जाने वाले हैं और इसे एक शक्तिशाली बैच ऑटोमेशन टूल में बदलने वाले हैं ताकि यह मनमाने संख्या के नमूनों को संभाल सके।
और जब हम यह कर रहे हैं, हम इसे युग्मित-अंत डेटा की अपेक्षा करने के लिए भी स्विच करने वाले हैं, जो नए अध्ययनों में अधिक सामान्य है।

हम यह तीन चरणों में करेंगे:

1. workflow को कई इनपुट नमूने स्वीकार करने योग्य बनाएं और निष्पादन को समानांतरित करें
2. व्यापक QC रिपोर्ट जनरेशन जोड़ें
3. युग्मित-अंत RNAseq डेटा पर स्विच करें

---

## 1. workflow को कई इनपुट नमूने स्वीकार करने योग्य बनाएं और निष्पादन को समानांतरित करें

हमें इनपुट को प्रबंधित करने के तरीके को बदलना होगा।

### 1.1. प्राथमिक इनपुट को एकल फ़ाइल के बजाय फ़ाइल पथों के CSV में बदलें

हम `data/` डायरेक्टरी में नमूना IDs और FASTQ फ़ाइल पथों वाली एक CSV फ़ाइल प्रदान करते हैं।
यह CSV फ़ाइल एक हेडर लाइन शामिल करती है।
ध्यान दें कि FASTQ फ़ाइल पथ पूर्ण पथ हैं।

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

आइए प्राथमिक इनपुट पैरामीटर का नाम बदलकर `input_csv` करें और डिफ़ॉल्ट को `single-end.csv` फ़ाइल के पथ में बदलें।

```groovy title="rnaseq.nf" linenums="13"
params {
    // प्राथमिक इनपुट
    input_csv: Path = "data/single-end.csv"

    // संदर्भ जीनोम आर्काइव
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 1.2. इनपुट channel फ़ैक्टरी को CSV को इनपुट के रूप में संभालने के लिए अपडेट करें

हम फ़ाइल के सामग्री को केवल फ़ाइल पथ के बजाय channel में लोड करना चाहते हैं, इसलिए हम CSV प्रारूप को पार्स करने के लिए `.splitCsv()` ऑपरेटर का उपयोग करते हैं, फिर जानकारी के विशिष्ट टुकड़े को पकड़ने के लिए `.map()` ऑपरेटर का उपयोग करते हैं (FASTQ फ़ाइल पथ)।

```groovy title="rnaseq.nf" linenums="16"
    // CSV फ़ाइल की सामग्री से इनपुट channel बनाएं
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }
```

### 1.3. workflow को चलाएं और परीक्षण करें कि यह काम करता है

```bash
nextflow run rnaseq.nf
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

इस बार हम देखते हैं कि हर चरण 6 बार चलता है, हर उस 6 डेटा फ़ाइलों पर जो हमने प्रदान की थीं।

बस इतना ही काफी था workflow को कई फ़ाइलों पर चलाने के लिए!
Nextflow हमारे लिए सभी समानांतरता को संभालता है।

---

## 2. प्री-प्रोसेसिंग QC मेट्रिक्स को एकल MultiQC रिपोर्ट में एकत्रित करें

यह सब बहुत सारी QC रिपोर्ट उत्पन्न करता है, और हम व्यक्तिगत रिपोर्टों को खोदना नहीं चाहते।
यह MultiQC रिपोर्ट एकत्रीकरण चरण लगाने के लिए एकदम सही बिंदु है!

### 2.1. QC एकत्रीकरण process के लिए एक मॉड्यूल बनाएं

आइए `MULTIQC` process को रखने के लिए `modules/multiqc.nf` नामक एक मॉड्यूल फ़ाइल बनाएं:

```bash
touch modules/multiqc.nf
```

कोड एडिटर में फ़ाइल खोलें और निम्नलिखित कोड को इसमें कॉपी करें:

```groovy title="modules/multiqc.nf" linenums="1"
#!/usr/bin/env nextflow

process MULTIQC {

    container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"
    publishDir "results/multiqc", mode: 'symlink'

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

### 2.2. मॉड्यूल को workflow फ़ाइल में आयात करें

`rnaseq.nf` फ़ाइल में स्टेटमेंट `include { MULTIQC } from './modules/multiqc.nf'` जोड़ें:

```groovy title="rnaseq.nf" linenums="3"
// मॉड्यूल INCLUDE स्टेटमेंट
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
```

### 2.3. एक `report_id` पैरामीटर जोड़ें और इसे एक उचित डिफ़ॉल्ट दें

```groovy title="rnaseq.nf" linenums="9"
params {
    // प्राथमिक इनपुट
    input_csv: Path = "data/single-end.csv"

    // संदर्भ जीनोम आर्काइव
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // रिपोर्ट ID
    report_id: String = "all_single-end"
}
```

### 2.4. पिछले चरणों के आउटपुट पर process को कॉल करें

हमें `MULTIQC` process को पिछले चरणों से सभी QC-संबंधित आउटपुट देने की आवश्यकता है।

इसके लिए, हम `.mix()` ऑपरेटर का उपयोग करने वाले हैं, जो कई channels को एक में एकत्रित करता है।

यदि हमारे पास A, B, C और D नामक चार processes थे जिनमें प्रत्येक में एक सरल `.out` channel था, तो सिंटैक्स कुछ इस तरह दिखेगा: `A.out.mix( B.out, C.out, D.out )`। जैसा कि आप देख सकते हैं, आप इसे पहले channel पर लागू करते हैं जिसे आप संयोजित करना चाहते हैं (कोई फर्क नहीं पड़ता कि कौन सा) और बाकी सभी को, कॉमा से अलग करके, कोष्ठक में जोड़ते हैं।

हमारे workflow के मामले में, हमारे पास एकत्रित करने के लिए निम्नलिखित आउटपुट हैं:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

तो सिंटैक्स उदाहरण बन जाता है:

```groovy title="MULTIQC कॉल में .mix() लागू करना"
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        )
```

यह प्रति नमूना QC रिपोर्ट एकत्र करेगा।
लेकिन चूंकि हम उन्हें सभी नमूनों में एकत्रित करना चाहते हैं, हमें सभी नमूनों के लिए रिपोर्ट को `MULTIQC` के एक कॉल में खींचने के लिए `collect()` ऑपरेटर जोड़ना होगा।
और हमें इसे `report_id` पैरामीटर भी देना होगा।

यह हमें निम्नलिखित देता है:

```groovy title="पूर्ण MULTIQC कॉल" linenums="33"
    // व्यापक QC रिपोर्ट जनरेशन
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

पूर्ण workflow ब्लॉक के संदर्भ में, यह इस तरह दिखता है:

```groovy title="rnaseq.nf" linenums="18"
workflow {
    // CSV फ़ाइल की सामग्री से इनपुट channel बनाएं
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }

    /// प्रारंभिक गुणवत्ता नियंत्रण
    FASTQC(read_ch)

    // एडॉप्टर ट्रिमिंग और पोस्ट-ट्रिमिंग QC
    TRIM_GALORE(read_ch)

    // संदर्भ जीनोम के साथ संरेखण
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))

    // व्यापक QC रिपोर्ट जनरेशन
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
}
```

### 2.5. workflow को चलाएं और परीक्षण करें कि यह काम करता है

```bash
nextflow run rnaseq.nf -resume
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

इस बार हम कैश किए गए process कॉलों के बाद MULTIQC के लिए एक एकल कॉल देखते हैं:

आप `TRIM_GALORE` process में `publishDir` निर्देश द्वारा निर्दिष्ट `results/trimming` के अंतर्गत आउटपुट पा सकते हैं।

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

---

## 3. युग्मित-अंत RNAseq डेटा प्रोसेसिंग सक्षम करें

अभी हमारा workflow केवल single-end RNAseq डेटा को संभाल सकता है।
युग्मित-अंत RNAseq डेटा देखना तेजी से सामान्य हो रहा है, इसलिए हम इसे संभालने में सक्षम होना चाहते हैं।

workflow को डेटा प्रकार से पूरी तरह अज्ञेयवादी बनाने के लिए थोड़ा अधिक उन्नत Nextflow भाषा सुविधाओं का उपयोग करने की आवश्यकता होगी, इसलिए हम यहां ऐसा नहीं करने वाले हैं, लेकिन हम यह प्रदर्शित करने के लिए एक युग्मित-अंत प्रोसेसिंग संस्करण बना सकते हैं कि क्या अनुकूलित करने की आवश्यकता है।

### 3.1. workflow की एक प्रति बनाएं जिसे `rnaseq_pe.nf` कहा जाता है

```bash
cp rnaseq.nf rnaseq_pe.nf
```

### 3.2. डिफ़ॉल्ट `input_csv` को युग्मित-अंत डेटा की ओर इंगित करने के लिए संशोधित करें

हम `data/` डायरेक्टरी में नमूना IDs और युग्मित FASTQ फ़ाइल पथों वाली दूसरी CSV फ़ाइल प्रदान करते हैं

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

आइए `input_csv` डिफ़ॉल्ट को `paired-end.csv` फ़ाइल का पथ बनाएं।

```groovy title="rnaseq_pe.nf" linenums="15"
params {
    // प्राथमिक इनपुट
    input_csv: Path = "data/paired-end.csv"

    // संदर्भ जीनोम आर्काइव
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // रिपोर्ट ID
    report_id: String = "all_single-end"
}
```

### 3.3. channel फ़ैक्टरी को अपडेट करें

हमें `.map()` ऑपरेटर को दोनों FASTQ फ़ाइल पथ अब पकड़ने के लिए बताना होगा।

तो `row -> file(row.fastq_path)` बन जाता है `row -> [file(row.fastq_1), file(row.fastq_2)]`

```groovy title="rnaseq_pe.nf" linenums="19"
    // CSV फ़ाइल की सामग्री से इनपुट channel बनाएं
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
```

### 3.4. FASTQC process का एक युग्मित-अंत संस्करण बनाएं

आइए मॉड्यूल की एक प्रति बनाएं ताकि हम दोनों संस्करणों को हाथ में रख सकें।

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

कोड एडिटर में नई `fastqc_pe.nf` मॉड्यूल फ़ाइल खोलें और निम्नलिखित कोड परिवर्तन करें:

- `script` ब्लॉक (पंक्ति 17) में `fastqc $reads` को `fastqc ${reads}` में बदलें ताकि `reads` इनपुट अनपैक हो जाए, क्योंकि यह अब एकल पथ के बजाय दो पथों का एक टपल है।
- आउटपुट फ़ाइलों को व्यक्तिगत रूप से संभालने से बचने के लिए `${reads.simpleName}` को वाइल्डकार्ड (`*`) से बदलें।

```groovy title="modules/fastqc_pe.nf" linenums="8"
    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads}
    """
```

तकनीकी रूप से यह `FASTQC` process को इस तरह से सामान्यीकृत करता है कि यह single-end या युग्मित-अंत RNAseq डेटा को संभालने में सक्षम हो जाता है।

अंत में, मॉड्यूल के युग्मित-अंत संस्करण का उपयोग करने के लिए मॉड्यूल आयात स्टेटमेंट को अपडेट करें।

```groovy title="rnaseq_pe.nf" linenums="4"
include { FASTQC } from './modules/fastqc_pe.nf'
```

### 3.5. TRIM_GALORE process का एक युग्मित-अंत संस्करण बनाएं

मॉड्यूल की एक प्रति बनाएं ताकि हम दोनों संस्करणों को हाथ में रख सकें।

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

कोड एडिटर में नई `trim_galore_pe.nf` मॉड्यूल फ़ाइल खोलें और निम्नलिखित कोड परिवर्तन करें:

- इनपुट घोषणा को `path reads` से `tuple path(read1), path(read2)` में बदलें
- `script` ब्लॉक में कमांड को अपडेट करें, `$reads` को `--paired ${read1} ${read2}` से बदलें
- जोड़ी गई फ़ाइलों और विभिन्न नामकरण परंपराओं को प्रतिबिंबित करने के लिए आउटपुट घोषणाओं को अपडेट करें, सब कुछ सूचीबद्ध करने से बचने के लिए वाइल्डकार्ड का उपयोग करें।

```groovy title="modules/trim_galore_pe.nf" linenums="8"
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

अंत में, मॉड्यूल के युग्मित-अंत संस्करण का उपयोग करने के लिए मॉड्यूल आयात स्टेटमेंट को अपडेट करें।

```groovy title="rnaseq_pe.nf" linenums="5"
include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
```

### 3.6. MULTIQC process के कॉल को TRIM_GALORE से दो रिपोर्ट अपेक्षित करने के लिए अपडेट करें

`TRIM_GALORE` process अब एक अतिरिक्त आउटपुट channel उत्पन्न करती है, इसलिए हमें इसे MultiQC को फीड करना होगा।

`TRIM_GALORE.out.fastqc_reports,` को `TRIM_GALORE.out.fastqc_reports_1,` और `TRIM_GALORE.out.fastqc_reports_2,` से बदलें:

```groovy title="rnaseq_pe.nf" linenums="33"
    // व्यापक QC रिपोर्ट जनरेशन
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports_1,
        TRIM_GALORE.out.fastqc_reports_2,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

जब हम MultiQC पर हैं, तो आइए `report_id` पैरामीटर डिफ़ॉल्ट को `"all_single-end"` से `"all_paired-end"` में भी अपडेट करें।

```groovy title="rnaseq_pe.nf" linenums="9"
params {
    // प्राथमिक इनपुट
    input_csv: Path = "data/paired-end.csv"

    // संदर्भ जीनोम आर्काइव
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // रिपोर्ट ID
    report_id: String = "all_paired-end"
}
```

### 3.7. HISAT2_ALIGN process का एक युग्मित-अंत संस्करण बनाएं

मॉड्यूल की एक प्रति बनाएं ताकि हम दोनों संस्करणों को हाथ में रख सकें।

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

कोड एडिटर में नई `hisat2_align_pe.nf` मॉड्यूल फ़ाइल खोलें और निम्नलिखित कोड परिवर्तन करें:

- इनपुट घोषणा को `path reads` से `tuple path(read1), path(read2)` में बदलें
- `script` ब्लॉक में कमांड को अपडेट करें, `-U $reads` को `-1 ${read1} -2 ${read2}` से बदलें
- `script` ब्लॉक में कमांड के साथ-साथ आउटपुट घोषणाओं में `${reads.simpleName}` के सभी उदाहरणों को `${read1.simpleName}` से बदलें।

```groovy title="modules/hisat2_align_pe.nf" linenums="8"
    input:
    tuple path(read1), path(read2)
    path index_zip

    output:
    path "${read1.simpleName}.bam", emit: bam
    path "${read1.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
        --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
        samtools view -bS -o ${read1.simpleName}.bam
    """
```

अंत में, मॉड्यूल के युग्मित-अंत संस्करण का उपयोग करने के लिए मॉड्यूल आयात स्टेटमेंट को अपडेट करें।

```groovy title="rnaseq_pe.nf" linenums="5"
include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
```

### 3.8. workflow को चलाएं और परीक्षण करें कि यह काम करता है

हम `-resume` का उपयोग नहीं करते क्योंकि यह कैश नहीं होगा, और पहले की तुलना में दोगुना डेटा प्रोसेस करना है, लेकिन फिर भी यह एक मिनट से कम में पूरा होना चाहिए।

```bash
nextflow run rnaseq_pe.nf
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

और बस इतना ही! अब हमारे पास हमारे workflow के दो थोड़े भिन्न संस्करण हैं, एक single-end रीड डेटा के लिए और एक युग्मित-अंत डेटा के लिए।
अगला तार्किक कदम workflow को किसी भी डेटा प्रकार को तुरंत स्वीकार करने योग्य बनाना होगा, जो इस पाठ्यक्रम के दायरे से बाहर है, लेकिन हम इसे एक फॉलो-अप में संभाल सकते हैं।

---

### मुख्य बातें

आप जानते हैं कि एकल-नमूना workflow को कई नमूनों की प्रोसेसिंग को समानांतरित करने के लिए कैसे अनुकूलित करें, एक व्यापक QC रिपोर्ट उत्पन्न करें और यदि आवश्यक हो तो युग्मित-अंत रीड डेटा का उपयोग करने के लिए workflow को अनुकूलित करें।

### आगे क्या है?

बधाई हो, आपने Nextflow For RNAseq मिनी-पाठ्यक्रम पूरा कर लिया है! अपनी सफलता का जश्न मनाएं और एक अच्छा आराम लें!

इसके बाद, हम आपसे इस प्रशिक्षण पाठ्यक्रम के साथ अपने अनुभव के बारे में एक बहुत छोटा सर्वेक्षण पूरा करने के लिए कहते हैं, फिर हम आपको आगे के प्रशिक्षण संसाधनों और सहायक लिंक के साथ एक पेज पर ले जाएंगे।
