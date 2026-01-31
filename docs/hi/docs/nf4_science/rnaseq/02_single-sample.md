# भाग 2: एकल-नमूना कार्यान्वयन

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

कोर्स के इस भाग में, हम सबसे सरल संभव workflow लिखने जा रहे हैं जो भाग 1 में चलाए गए सभी commands को wrap करता है ताकि उन्हें चलाना automate हो सके, और हम एक बार में सिर्फ एक नमूना प्रोसेस करने का लक्ष्य रखेंगे।

हम यह तीन चरणों में करेंगे:

1. एक single-stage workflow लिखें जो शुरुआती QC step चलाता है
2. Adapter trimming और post-trimming QC जोड़ें
3. Reference genome के लिए alignment जोड़ें

!!! warning "पूर्वापेक्षा"

    इस पाठ को शुरू करने से पहले आपको कोर्स के भाग 1 को पूरा करना होगा।
    विशेष रूप से, sections 2.1-3 पर काम करने से genome index फ़ाइल (`data/genome_index.tar.gz`) बनती है जो इस पाठ में alignment step के लिए आवश्यक है।

---

## 1. एक single-stage workflow लिखें जो शुरुआती QC चलाता है

आइए एक सरल workflow लिखकर शुरू करें जो single-end RNAseq reads वाली FASTQ फ़ाइल पर FastQC tool चलाता है।

हम आपको एक workflow फ़ाइल, `rnaseq.nf`, प्रदान करते हैं, जो workflow के मुख्य भागों की रूपरेखा देती है।

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Module INCLUDE statements

/*
 * Pipeline parameters
 */

// प्राथमिक इनपुट

workflow {

    // इनपुट channel बनाएं

    // Processes को call करें

}
```

ध्यान रखें यह workflow code सही है लेकिन यह functional नहीं है; इसका उद्देश्य सिर्फ एक skeleton के रूप में काम करना है जिसका उपयोग आप वास्तविक workflow लिखने के लिए करेंगे।

### 1.1. Modules स्टोर करने के लिए एक डायरेक्टरी बनाएं

हम प्रत्येक process के लिए standalone modules बनाएंगे ताकि उन्हें manage करना और reuse करना आसान हो, तो चलिए उन्हें स्टोर करने के लिए एक डायरेक्टरी बनाते हैं।

```bash
mkdir modules
```

### 1.2. QC metrics collection process के लिए एक मॉड्यूल बनाएं

आइए `FASTQC` process को रखने के लिए `modules/fastqc.nf` नाम की एक मॉड्यूल फ़ाइल बनाएं:

```bash
touch modules/fastqc.nf
```

फ़ाइल को code editor में खोलें और निम्नलिखित code को इसमें कॉपी करें:

```groovy title="modules/fastqc.nf" linenums="1"
#!/usr/bin/env nextflow

process FASTQC {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/fastqc", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_fastqc.zip", emit: zip
    path "${reads.simpleName}_fastqc.html", emit: html

    script:
    """
    fastqc $reads
    """
}
```

आपको इस training series के भाग 1 और भाग 2 में जो सीखा है उससे सभी टुकड़े पहचानने चाहिए; एकमात्र उल्लेखनीय बदलाव यह है कि इस बार हम `publishDir` directive के लिए `mode: symlink` का उपयोग कर रहे हैं, और हम `publishDir` को define करने के लिए एक पैरामीटर का उपयोग कर रहे हैं।

!!! note

    भले ही हम यहां जो data फ़ाइलें उपयोग कर रहे हैं वे बहुत छोटी हैं, genomics में वे बहुत बड़ी हो सकती हैं। teaching वातावरण में प्रदर्शन के उद्देश्य से, हम अनावश्यक फ़ाइल copies से बचने के लिए 'symlink' publishing mode का उपयोग कर रहे हैं। आपको अपने final workflows में ऐसा नहीं करना चाहिए, क्योंकि जब आप अपनी `work` डायरेक्टरी को साफ करेंगे तो आप results खो देंगे।

### 1.3. Workflow फ़ाइल में मॉड्यूल को import करें

`rnaseq.nf` फ़ाइल में statement `include { FASTQC } from './modules/fastqc.nf'` जोड़ें:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
```

### 1.4. एक इनपुट declaration जोड़ें

एक default value के साथ एक इनपुट पैरामीटर declare करें:

```groovy title="rnaseq.nf" linenums="10"
params {
    // प्राथमिक इनपुट
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
}
```

### 1.5. Workflow block में एक इनपुट channel बनाएं

इनपुट channel बनाने के लिए एक basic `.fromPath()` channel factory का उपयोग करें:

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // एक फ़ाइल path से इनपुट channel बनाएं
    read_ch = channel.fromPath(params.reads)

    // Processes को call करें

}
```

### 1.6. इनपुट channel पर `FASTQC` process को call करें

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // एक फ़ाइल path से इनपुट channel बनाएं
    read_ch = channel.fromPath(params.reads)

    // प्रारंभिक quality control
    FASTQC(read_ch)

}
```

### 1.7. Workflow को चलाएं यह जांचने के लिए कि यह काम करता है

हम command line से इनपुट specify करने के लिए `--reads` पैरामीटर का उपयोग कर सकते हैं, लेकिन development के दौरान हम lazy हो सकते हैं और बस हमारे द्वारा set किए गए test default का उपयोग कर सकते हैं।

```bash
nextflow run rnaseq.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    ```

यह बहुत जल्दी चलना चाहिए अगर आपने भाग 1 पर काम किया है और पहले से ही कंटेनर को pull कर लिया है।
यदि आपने इसे छोड़ दिया है, तो Nextflow आपके लिए कंटेनर को pull करेगा; इसके होने के लिए आपको कुछ भी करने की आवश्यकता नहीं है, लेकिन आपको एक मिनट तक प्रतीक्षा करने की आवश्यकता हो सकती है।

आप `FASTQC` process द्वारा `publishDir` directive द्वारा निर्दिष्ट के अनुसार `results/fastqc` के तहत आउटपुट पा सकते हैं।

```bash
ls results/fastqc
```

```console title="आउटपुट"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

## 2. Adapter trimming और post-trimming quality control जोड़ें

हम Trim_Galore wrapper का उपयोग करने जा रहे हैं, जो trimming के लिए Cutadapt और post-trimming quality control के लिए FastQC को bundle करता है।

### 2.1. Trimming और QC process के लिए एक मॉड्यूल बनाएं

आइए `TRIM_GALORE` process को रखने के लिए `modules/trim_galore.nf` नाम की एक मॉड्यूल फ़ाइल बनाएं:

```bash
touch modules/trim_galore.nf
```

फ़ाइल को code editor में खोलें और निम्नलिखित code को इसमें कॉपी करें:

```groovy title="modules/trim_galore.nf" linenums="1"
#!/usr/bin/env nextflow

process TRIM_GALORE {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/trimming", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
    path "${reads}_trimming_report.txt", emit: trimming_reports
    path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

    script:
    """
    trim_galore --fastqc $reads
    """
}
```

### 2.2. Workflow फ़ाइल में मॉड्यूल को import करें

`rnaseq.nf` फ़ाइल में statement `include { TRIM_GALORE } from './modules/trim_galore.nf'` जोड़ें:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. इनपुट channel पर process को call करें

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // एक फ़ाइल path से इनपुट channel बनाएं
    read_ch = channel.fromPath(params.reads)

    // प्रारंभिक quality control
    FASTQC(read_ch)

    // Adapter trimming और post-trimming QC
    TRIM_GALORE(read_ch)
}
```

### 2.4. Workflow को चलाएं यह जांचने के लिए कि यह काम करता है

```bash
nextflow run rnaseq.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    [c2/e4a9bb] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    ```

यह भी बहुत जल्दी चलना चाहिए, क्योंकि हम इतनी छोटी इनपुट फ़ाइल पर चला रहे हैं।

आप `TRIM_GALORE` process द्वारा `publishDir` directive द्वारा निर्दिष्ट के अनुसार `results/trimming` के तहत आउटपुट पा सकते हैं।

```bash
ls results/trimming
```

```console title="आउटपुट"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. Reads को reference genome से align करें

अंत में हम Hisat2 का उपयोग करके genome alignment step चला सकते हैं, जो FastQC-style quality control metrics भी emit करेगा।

### 3.1. HiSat2 process के लिए एक मॉड्यूल बनाएं

आइए `HISAT2_ALIGN` process को रखने के लिए `modules/hisat2_align.nf` नाम की एक मॉड्यूल फ़ाइल बनाएं:

```bash
touch modules/hisat2_align.nf
```

फ़ाइल को code editor में खोलें और निम्नलिखित code को इसमें कॉपी करें:

```groovy title="modules/hisat2_align.nf" linenums="1"
#!/usr/bin/env nextflow

process HISAT2_ALIGN {

    container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"
    publishDir "results/align", mode: 'symlink'

    input:
    path reads
    path index_zip

    output:
    path "${reads.simpleName}.bam", emit: bam
    path "${reads.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -U $reads \
        --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
        samtools view -bS -o ${reads.simpleName}.bam
    """
}
```

### 3.2. Workflow फ़ाइल में मॉड्यूल को import करें

`rnaseq.nf` फ़ाइल में statement `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` जोड़ें:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. Genome index प्रदान करने के लिए एक पैरामीटर declaration जोड़ें

एक default value के साथ एक इनपुट पैरामीटर declare करें:

```groovy title="rnaseq.nf" linenums="8"
params {
    // प्राथमिक इनपुट
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

    // Reference genome archive (संदर्भ जीनोम आर्काइव)
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 3.4. `TRIM_GALORE` द्वारा आउटपुट किए गए trimmed reads पर `HISAT2_ALIGN` process को call करें

Trimmed reads पिछले step द्वारा आउटपुट किए गए `TRIM_GALORE.out.trimmed_reads` channel में हैं।

इसके अलावा, हम Hisat2 tool को gzipped genome index tarball प्रदान करने के लिए `file (params.hisat2_index_zip)` का उपयोग करते हैं।

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // एक फ़ाइल path से इनपुट channel बनाएं
    read_ch = channel.fromPath(params.reads)

    // प्रारंभिक quality control
    FASTQC(read_ch)

    // Adapter trimming और post-trimming QC
    TRIM_GALORE(read_ch)

    // Reference genome के साथ alignment
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.5. Workflow को चलाएं यह जांचने के लिए कि यह काम करता है

```bash
nextflow run rnaseq.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [extravagant_khorana] DSL2 - revision: 701b41bd16

    executor >  local (3)
    [e4/d15ad4] FASTQC (1)       [100%] 1 of 1 ✔
    [c6/12b2be] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    [c6/7a9f13] HISAT2_ALIGN (1) [100%] 1 of 1 ✔
    ```

आप `HISAT2_ALIGN` process द्वारा `publishDir` directive द्वारा निर्दिष्ट के अनुसार `results/align` के तहत आउटपुट पा सकते हैं।

```bash
ls results/align
```

```console title="आउटपुट"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

यह प्रत्येक नमूने पर लागू करने के लिए आवश्यक basic processing को पूरा करता है।

_हम भाग 2 में MultiQC report aggregation जोड़ेंगे, जब हम workflow को एक बार में कई नमूने स्वीकार करने योग्य बना लेंगे।_

---

### सारांश

आप जानते हैं कि single-end RNAseq नमूनों को व्यक्तिगत रूप से प्रोसेस करने के लिए सभी मुख्य steps को कैसे wrap करें।

### आगे क्या है?

जानें कि कैसे workflow को कई नमूनों को parallel में प्रोसेस करने के लिए modify करें, सभी नमूनों के लिए सभी steps में QC reports को aggregate करें, और paired-end RNAseq data पर workflow चलाना enable करें।
