# भाग 3: मल्टी-सैंपल पेयर्ड-एंड इम्प्लीमेंटेशन

इस कोर्स के अंतिम भाग में, हम अपने सरल workflow को अगले स्तर पर ले जाने वाले हैं और इसे एक शक्तिशाली बैच ऑटोमेशन टूल में बदलेंगे जो किसी भी संख्या में सैंपल को हैंडल कर सके।
और जब हम यह कर रहे हैं, तो हम इसे पेयर्ड-एंड डेटा के लिए भी स्विच करेंगे, जो नए अध्ययनों में अधिक सामान्य है।

हम यह तीन चरणों में करेंगे:

1. Workflow को कई इनपुट सैंपल स्वीकार करने और निष्पादन को समानांतर बनाने के लिए तैयार करना
2. व्यापक QC रिपोर्ट जनरेशन जोड़ना
3. पेयर्ड-एंड RNAseq डेटा पर स्विच करना

---

## 1. Workflow को कई इनपुट सैंपल स्वीकार करने और निष्पादन को समानांतर बनाने के लिए तैयार करना

हमें इनपुट को मैनेज करने के तरीके को बदलना होगा।

### 1.1. प्राथमिक इनपुट को एक सिंगल फ़ाइल के बजाय फ़ाइल पाथ की CSV में बदलें

हम `data/` डायरेक्टरी में सैंपल ID और FASTQ फ़ाइल पाथ वाली एक CSV फ़ाइल प्रदान करते हैं।
इस CSV फ़ाइल में एक हेडर लाइन शामिल है।
ध्यान दें कि FASTQ फ़ाइल पाथ एब्सोल्यूट पाथ हैं।

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

चलो प्राथमिक इनपुट पैरामीटर का नाम बदलकर `input_csv` करते हैं और डिफ़ॉल्ट को `single-end.csv` फ़ाइल के पाथ में बदलते हैं।

```groovy title="rnaseq.nf" linenums="13"
params {
    // Primary input
    input_csv: Path = "data/single-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 1.2. इनपुट चैनल फ़ैक्टरी को CSV को इनपुट के रूप में हैंडल करने के लिए अपडेट करें

हम फ़ाइल की सामग्री को चैनल में लोड करना चाहेंगे, सिर्फ फ़ाइल पाथ के बजाय, इसलिए हम CSV फ़ॉर्मेट को पार्स करने के लिए `.splitCsv()` ऑपरेटर का उपयोग करते हैं, फिर `.map()` ऑपरेटर का उपयोग करके वह विशिष्ट जानकारी प्राप्त करते हैं जो हम चाहते हैं (FASTQ फ़ाइल पाथ)।

```groovy title="rnaseq.nf" linenums="16"
    // Create input channel from the contents of a CSV file
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }
```

### 1.3. Workflow को चलाएं और टेस्ट करें कि यह काम करता है

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

इस बार हम देखते हैं कि प्रत्येक स्टेप 6 बार चलता है, हमारे द्वारा प्रदान की गई 6 डेटा फ़ाइलों में से प्रत्येक पर।

Workflow को कई फ़ाइलों पर चलाने के लिए बस इतना ही करना था!
Nextflow हमारे लिए सभी समानांतरता को हैंडल करता है।

---

## 2. प्री-प्रोसेसिंग QC मेट्रिक्स को एक सिंगल MultiQC रिपोर्ट में एकत्रित करें

यह सब बहुत सारी QC रिपोर्ट उत्पन्न करता है, और हम व्यक्तिगत रिपोर्ट को खोदना नहीं चाहते।
यह MultiQC रिपोर्ट एग्रीगेशन स्टेप डालने का सही बिंदु है!

### 2.1. QC एग्रीगेशन प्रोसेस के लिए एक मॉड्यूल बनाएं

चलो `MULTIQC` प्रोसेस को रखने के लिए `modules/multiqc.nf` नामक एक मॉड्यूल फ़ाइल बनाते हैं:

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

### 2.2. मॉड्यूल को workflow फ़ाइल में इम्पोर्ट करें

`rnaseq.nf` फ़ाइल में स्टेटमेंट `include { MULTIQC } from './modules/multiqc.nf'` जोड़ें:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
```

### 2.3. एक `report_id` पैरामीटर जोड़ें और इसे एक उचित डिफ़ॉल्ट दें

```groovy title="rnaseq.nf" linenums="9"
params {
    // Primary input
    input_csv: Path = "data/single-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_single-end"
}
```

### 2.4. पिछले स्टेप्स के आउटपुट पर प्रोसेस को कॉल करें

हमें `MULTIQC` प्रोसेस को पिछले स्टेप्स से सभी QC-संबंधित आउटपुट देने होंगे।

इसके लिए, हम `.mix()` ऑपरेटर का उपयोग करने जा रहे हैं, जो कई चैनलों को एक में एकत्रित करता है।

यदि हमारे पास A, B, C और D नामक चार प्रोसेस हैं, प्रत्येक में एक सरल `.out` चैनल है, तो सिंटैक्स इस तरह दिखेगा: `A.out.mix( B.out, C.out, D.out )`। जैसा कि आप देख सकते हैं, आप इसे उन चैनलों में से पहले पर लागू करते हैं जिन्हें आप कंबाइन करना चाहते हैं (कोई फर्क नहीं पड़ता कि कौन सा) और बाकी सभी को कॉमा से अलग करके, उसके बाद के कोष्ठक में जोड़ते हैं।

हमारे workflow के मामले में, हमारे पास एकत्रित करने के लिए निम्नलिखित आउटपुट हैं:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

तो सिंटैक्स उदाहरण बन जाता है:

```groovy title="Applying .mix() in the MULTIQC call"
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        )
```

यह प्रति सैंपल QC रिपोर्ट एकत्र करेगा।
लेकिन चूंकि हम उन्हें सभी सैंपल में एकत्रित करना चाहते हैं, इसलिए हमें सभी सैंपल के लिए रिपोर्ट को `MULTIQC` की एक सिंगल कॉल में खींचने के लिए `collect()` ऑपरेटर जोड़ने की आवश्यकता है।
और हमें इसे `report_id` पैरामीटर भी देना होगा।

यह हमें निम्नलिखित देता है:

```groovy title="The completed MULTIQC call" linenums="33"
    // Comprehensive QC report generation
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
    // Create input channel from the contents of a CSV file
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }

    /// Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))

    // Comprehensive QC report generation
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

### 2.5. Workflow को चलाएं और टेस्ट करें कि यह काम करता है

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

इस बार हम कैश्ड प्रोसेस कॉल के बाद MULTIQC की एक सिंगल कॉल देखते हैं:

आप `TRIM_GALORE` प्रोसेस में `publishDir` डायरेक्टिव द्वारा निर्दिष्ट `results/trimming` के तहत आउटपुट पा सकते हैं।

```bash
tree -L 2 results/multiqc
```

```console title="Output"
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

वह अंतिम `all_single-end.html` फ़ाइल पूर्ण एकत्रित रिपोर्ट है, जो सुविधाजनक रूप से एक आसान ब्राउज़ करने योग्य HTML फ़ाइल में पैक की गई है।

---

## 3. पेयर्ड-एंड RNAseq डेटा की प्रोसेसिंग को सक्षम करें

अभी हमारा workflow केवल सिंगल-एंड RNAseq डेटा को हैंडल कर सकता है।
पेयर्ड-एंड RNAseq डेटा देखना तेजी से सामान्य हो रहा है, इसलिए हम इसे हैंडल करने में सक्षम होना चाहते हैं।

Workflow को डेटा टाइप से पूरी तरह से अज्ञेयवादी बनाने के लिए थोड़ी अधिक उन्नत Nextflow भाषा सुविधाओं का उपयोग करने की आवश्यकता होगी, इसलिए हम यहां ऐसा नहीं करने जा रहे हैं, लेकिन हम यह प्रदर्शित करने के लिए एक पेयर्ड-एंड प्रोसेसिंग संस्करण बना सकते हैं कि क्या अनुकूलित करने की आवश्यकता है।

### 3.1. Workflow की एक कॉपी बनाएं जिसे `rnaseq_pe.nf` कहा जाता है

```bash
cp rnaseq.nf rnaseq_pe.nf
```

### 3.2. डिफ़ॉल्ट `input_csv` को पेयर्ड-एंड डेटा की ओर इशारा करने के लिए संशोधित करें

हम `data/` डायरेक्टरी में सैंपल ID और पेयर्ड FASTQ फ़ाइल पाथ वाली एक दूसरी CSV फ़ाइल प्रदान करते हैं

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

चलो `input_csv` डिफ़ॉल्ट को `paired-end.csv` फ़ाइल के पाथ में बदलते हैं।

```groovy title="rnaseq_pe.nf" linenums="15"
params {
    // Primary input
    input_csv: Path = "data/paired-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_single-end"
}
```

### 3.3. चैनल फ़ैक्टरी को अपडेट करें

हमें `.map()` ऑपरेटर को दोनों FASTQ फ़ाइल पाथ को ग्रैब करने के लिए कहना होगा।

तो `row -> file(row.fastq_path)` बन जाता है `row -> [file(row.fastq_1), file(row.fastq_2)]`

```groovy title="rnaseq_pe.nf" linenums="19"
    // Create input channel from the contents of a CSV file
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
```

### 3.4. FASTQC प्रोसेस का एक पेयर्ड-एंड संस्करण बनाएं

चलो मॉड्यूल की एक कॉपी बनाते हैं ताकि हम दोनों संस्करण हाथ में रख सकें।

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

कोड एडिटर में नई `fastqc_pe.nf` मॉड्यूल फ़ाइल खोलें और निम्नलिखित कोड परिवर्तन करें:

- `script` ब्लॉक (लाइन 17) में `fastqc $reads` को `fastqc ${reads}` में बदलें ताकि `reads` इनपुट अनपैक हो जाए, क्योंकि यह अब एक सिंगल पाथ के बजाय दो पाथ का टपल है।
- आउटपुट फ़ाइलों को व्यक्तिगत रूप से हैंडल करने से बचने के लिए `${reads.simpleName}` को वाइल्डकार्ड (`*`) से बदलें।

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

तकनीकी रूप से यह `FASTQC` प्रोसेस को इस तरह से सामान्यीकृत करता है जो इसे सिंगल-एंड या पेयर्ड-एंड RNAseq डेटा को हैंडल करने में सक्षम बनाता है।

अंत में, मॉड्यूल के पेयर्ड-एंड संस्करण का उपयोग करने के लिए मॉड्यूल इम्पोर्ट स्टेटमेंट को अपडेट करें।

```groovy title="rnaseq_pe.nf" linenums="4"
include { FASTQC } from './modules/fastqc_pe.nf'
```

### 3.5. TRIM_GALORE प्रोसेस का एक पेयर्ड-एंड संस्करण बनाएं

मॉड्यूल की एक कॉपी बनाएं ताकि हम दोनों संस्करण हाथ में रख सकें।

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

कोड एडिटर में नई `trim_galore_pe.nf` मॉड्यूल फ़ाइल खोलें और निम्नलिखित कोड परिवर्तन करें:

- इनपुट डिक्लेरेशन को `path reads` से `tuple path(read1), path(read2)` में बदलें
- `script` ब्लॉक में कमांड को अपडेट करें, `$reads` को `--paired ${read1} ${read2}` से बदलें
- जोड़ी गई फ़ाइलों और विभिन्न नामकरण परंपराओं को प्रतिबिंबित करने के लिए आउटपुट डिक्लेरेशन को अपडेट करें, सब कुछ सूचीबद्ध करने से बचने के लिए वाइल्डकार्ड का उपयोग करें।

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

अंत में, मॉड्यूल के पेयर्ड-एंड संस्करण का उपयोग करने के लिए मॉड्यूल इम्पोर्ट स्टेटमेंट को अपडेट करें।

```groovy title="rnaseq_pe.nf" linenums="5"
include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
```

### 3.6. MULTIQC प्रोसेस की कॉल को TRIM_GALORE से दो रिपोर्ट की अपेक्षा करने के लिए अपडेट करें

`TRIM_GALORE` प्रोसेस अब एक अतिरिक्त आउटपुट चैनल उत्पन्न करता है, इसलिए हमें इसे MultiQC को फीड करना होगा।

`TRIM_GALORE.out.fastqc_reports,` को `TRIM_GALORE.out.fastqc_reports_1,` प्लस `TRIM_GALORE.out.fastqc_reports_2,` से बदलें:

```groovy title="rnaseq_pe.nf" linenums="33"
    // Comprehensive QC report generation
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

जब हम MultiQC पर हैं, तो चलो `report_id` पैरामीटर डिफ़ॉल्ट को `"all_single-end"` से `"all_paired-end"` में भी अपडेट करते हैं।

```groovy title="rnaseq_pe.nf" linenums="9"
params {
    // Primary input
    input_csv: Path = "data/paired-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_paired-end"
}
```

### 3.7. HISAT2_ALIGN प्रोसेस का एक पेयर्ड-एंड संस्करण बनाएं

मॉड्यूल की एक कॉपी बनाएं ताकि हम दोनों संस्करण हाथ में रख सकें।

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

कोड एडिटर में नई `hisat2_align_pe.nf` मॉड्यूल फ़ाइल खोलें और निम्नलिखित कोड परिवर्तन करें:

- इनपुट डिक्लेरेशन को `path reads` से `tuple path(read1), path(read2)` में बदलें
- `script` ब्लॉक में कमांड को अपडेट करें, `-U $reads` को `-1 ${read1} -2 ${read2}` से बदलें
- `script` ब्लॉक में कमांड के साथ-साथ आउटपुट डिक्लेरेशन में `${reads.simpleName}` के सभी इंस्टेंस को `${read1.simpleName}` से बदलें।

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

अंत में, मॉड्यूल के पेयर्ड-एंड संस्करण का उपयोग करने के लिए मॉड्यूल इम्पोर्ट स्टेटमेंट को अपडेट करें।

```groovy title="rnaseq_pe.nf" linenums="5"
include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
```

### 3.8. Workflow को चलाएं और टेस्ट करें कि यह काम करता है

हम `-resume` का उपयोग नहीं करते क्योंकि यह कैश नहीं होगा, और पहले की तुलना में दोगुना डेटा प्रोसेस करना है, लेकिन यह अभी भी एक मिनट से कम में पूरा होना चाहिए।

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

और बस! अब हमारे पास अपने workflow के दो थोड़े अलग संस्करण हैं, एक सिंगल-एंड रीड डेटा के लिए और एक पेयर्ड-एंड डेटा के लिए।
अगला तार्किक कदम workflow को किसी भी डेटा टाइप को फ्लाई पर स्वीकार करने के लिए बनाना होगा, जो इस कोर्स के दायरे से बाहर है, लेकिन हम इसे एक फॉलो-अप में निपटा सकते हैं।

---

### सारांश

तुम जानते हो कि एक सिंगल-सैंपल workflow को कई सैंपल की प्रोसेसिंग को समानांतर बनाने के लिए कैसे अनुकूलित किया जाए, एक व्यापक QC रिपोर्ट कैसे जनरेट की जाए और यदि आवश्यक हो तो workflow को पेयर्ड-एंड रीड डेटा का उपयोग करने के लिए कैसे अनुकूलित किया जाए।

### आगे क्या है?

बधाई हो, तुमने Nextflow For RNAseq मिनी-कोर्स पूरा कर लिया है! अपनी सफलता का जश्न मनाओ और एक अच्छा आराम लो!

इसके बाद, हम तुमसे इस प्रशिक्षण कोर्स के बारे में अपने अनुभव के बारे में एक बहुत छोटा सर्वेक्षण पूरा करने के लिए कहते हैं, फिर हम तुम्हें आगे के प्रशिक्षण संसाधनों और सहायक लिंक के साथ एक पेज पर ले जाएंगे।
