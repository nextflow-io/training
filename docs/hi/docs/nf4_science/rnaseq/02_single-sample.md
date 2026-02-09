# भाग 2: सिंगल-सैंपल इम्प्लीमेंटेशन

इस कोर्स के इस भाग में, हम सबसे सरल संभव वर्कफ़्लो लिखने जा रहे हैं जो भाग 1 में चलाए गए सभी कमांड्स को wrap करता है ताकि उन्हें automate किया जा सके, और हम एक बार में सिर्फ एक सैंपल को प्रोसेस करने का लक्ष्य रखेंगे।

हम इसे तीन चरणों में करेंगे:

1. एक single-stage वर्कफ़्लो लिखें जो शुरुआती QC स्टेप को चलाता है
2. Adapter trimming और post-trimming QC जोड़ें
3. Reference genome के साथ alignment जोड़ें

!!! warning "पूर्वापेक्षा"

    इस पाठ को शुरू करने से पहले तुम्हें कोर्स के भाग 1 को पूरा करना होगा।
    विशेष रूप से, सेक्शन 2.1-3 में काम करने से genome index फ़ाइल (`data/genome_index.tar.gz`) बनती है जो इस पाठ में alignment स्टेप के लिए आवश्यक है।

---

## 1. एक single-stage वर्कफ़्लो लिखें जो शुरुआती QC चलाता है

चलो एक सरल वर्कफ़्लो लिखकर शुरू करते हैं जो single-end RNAseq reads वाली FASTQ फ़ाइल पर FastQC टूल चलाता है।

हम तुम्हें एक वर्कफ़्लो फ़ाइल, `rnaseq.nf`, प्रदान करते हैं, जो वर्कफ़्लो के मुख्य भागों को outline करती है।

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Module INCLUDE statements

/*
 * Pipeline parameters
 */

// Primary input

workflow {

    // Create input channel

    // Call processes

}
```

ध्यान रखो कि यह वर्कफ़्लो कोड सही है लेकिन यह functional नहीं है; इसका उद्देश्य सिर्फ एक skeleton के रूप में काम करना है जिसे तुम वास्तविक वर्कफ़्लो लिखने के लिए उपयोग करोगे।

### 1.1. मॉड्यूल्स को स्टोर करने के लिए एक डायरेक्टरी बनाएं

हम प्रत्येक प्रोसेस के लिए standalone मॉड्यूल्स बनाएंगे ताकि उन्हें manage और reuse करना आसान हो, तो चलो उन्हें स्टोर करने के लिए एक डायरेक्टरी बनाते हैं।

```bash
mkdir modules
```

### 1.2. QC मेट्रिक्स कलेक्शन प्रोसेस के लिए एक मॉड्यूल बनाएं

चलो `FASTQC` प्रोसेस को रखने के लिए `modules/fastqc.nf` नाम की एक मॉड्यूल फ़ाइल बनाते हैं:

```bash
touch modules/fastqc.nf
```

कोड एडिटर में फ़ाइल खोलो और निम्नलिखित कोड को इसमें कॉपी करो:

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

तुम्हें इस ट्रेनिंग सीरीज़ के भाग 1 और भाग 2 में जो सीखा था उसके सभी हिस्से पहचानने चाहिए; एकमात्र उल्लेखनीय बदलाव यह है कि इस बार हम `publishDir` directive के लिए `mode: symlink` का उपयोग कर रहे हैं, और हम `publishDir` को define करने के लिए एक पैरामीटर का उपयोग कर रहे हैं।

!!! note

    भले ही हम यहां जिन डेटा फ़ाइलों का उपयोग कर रहे हैं वे बहुत छोटी हैं, genomics में वे बहुत बड़ी हो सकती हैं। शिक्षण वातावरण में प्रदर्शन के उद्देश्यों के लिए, हम अनावश्यक फ़ाइल कॉपी से बचने के लिए 'symlink' publishing mode का उपयोग कर रहे हैं। तुम्हें अपने अंतिम वर्कफ़्लो में ऐसा नहीं करना चाहिए, क्योंकि जब तुम अपनी `work` डायरेक्टरी को साफ करोगे तो तुम परिणाम खो दोगे।

### 1.3. मॉड्यूल को वर्कफ़्लो फ़ाइल में import करें

`rnaseq.nf` फ़ाइल में statement `include { FASTQC } from './modules/fastqc.nf'` जोड़ो:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
```

### 1.4. एक इनपुट declaration जोड़ें

एक डिफ़ॉल्ट वैल्यू के साथ एक इनपुट पैरामीटर declare करो:

```groovy title="rnaseq.nf" linenums="10"
params {
    // Primary input
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
}
```

### 1.5. वर्कफ़्लो ब्लॉक में एक इनपुट चैनल बनाएं

इनपुट चैनल बनाने के लिए एक बेसिक `.fromPath()` channel factory का उपयोग करो:

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Call processes

}
```

### 1.6. इनपुट चैनल पर `FASTQC` प्रोसेस को कॉल करें

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

}
```

### 1.7. वर्कफ़्लो को चलाकर टेस्ट करें कि यह काम करता है

हम कमांड लाइन से इनपुट specify करने के लिए `--reads` पैरामीटर का उपयोग कर सकते हैं, लेकिन development के दौरान हम आलसी हो सकते हैं और बस हमने जो टेस्ट डिफ़ॉल्ट सेट किया है उसका उपयोग कर सकते हैं।

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

यह बहुत जल्दी चलना चाहिए अगर तुमने भाग 1 में काम किया है और पहले से ही कंटेनर को pull कर लिया है।
अगर तुमने इसे छोड़ दिया है, तो Nextflow तुम्हारे लिए कंटेनर को pull करेगा; तुम्हें इसके होने के लिए कुछ भी करने की ज़रूरत नहीं है, लेकिन तुम्हें एक मिनट तक इंतज़ार करना पड़ सकता है।

तुम `results/fastqc` के अंदर आउटपुट पा सकते हो जैसा कि `FASTQC` प्रोसेस में `publishDir` directive द्वारा निर्दिष्ट किया गया है।

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

## 2. Adapter trimming और post-trimming quality control जोड़ें

हम Trim_Galore wrapper का उपयोग करने जा रहे हैं, जो trimming के लिए Cutadapt और post-trimming quality control के लिए FastQC को bundle करता है।

### 2.1. Trimming और QC प्रोसेस के लिए एक मॉड्यूल बनाएं

चलो `TRIM_GALORE` प्रोसेस को रखने के लिए `modules/trim_galore.nf` नाम की एक मॉड्यूल फ़ाइल बनाते हैं:

```bash
touch modules/trim_galore.nf
```

कोड एडिटर में फ़ाइल खोलो और निम्नलिखित कोड को इसमें कॉपी करो:

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

### 2.2. मॉड्यूल को वर्कफ़्लो फ़ाइल में import करें

`rnaseq.nf` फ़ाइल में statement `include { TRIM_GALORE } from './modules/trim_galore.nf'` जोड़ो:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. इनपुट चैनल पर प्रोसेस को कॉल करें

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)
}
```

### 2.4. वर्कफ़्लो को चलाकर टेस्ट करें कि यह काम करता है

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

तुम `results/trimming` के अंदर आउटपुट पा सकते हो जैसा कि `TRIM_GALORE` प्रोसेस में `publishDir` directive द्वारा निर्दिष्ट किया गया है।

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. Reads को reference genome के साथ align करें

अंत में हम Hisat2 का उपयोग करके genome alignment स्टेप चला सकते हैं, जो FastQC-style quality control मेट्रिक्स भी emit करेगा।

### 3.1. HiSat2 प्रोसेस के लिए एक मॉड्यूल बनाएं

चलो `HISAT2_ALIGN` प्रोसेस को रखने के लिए `modules/hisat2_align.nf` नाम की एक मॉड्यूल फ़ाइल बनाते हैं:

```bash
touch modules/hisat2_align.nf
```

कोड एडिटर में फ़ाइल खोलो और निम्नलिखित कोड को इसमें कॉपी करो:

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

### 3.2. मॉड्यूल को वर्कफ़्लो फ़ाइल में import करें

`rnaseq.nf` फ़ाइल में statement `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` जोड़ो:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. Genome index प्रदान करने के लिए एक पैरामीटर declaration जोड़ें

एक डिफ़ॉल्ट वैल्यू के साथ एक इनपुट पैरामीटर declare करो:

```groovy title="rnaseq.nf" linenums="8"
params {
    // Primary input
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 3.4. `TRIM_GALORE` द्वारा आउटपुट किए गए trimmed reads पर `HISAT2_ALIGN` प्रोसेस को कॉल करें

Trimmed reads पिछले स्टेप द्वारा आउटपुट किए गए `TRIM_GALORE.out.trimmed_reads` चैनल में हैं।

इसके अलावा, हम Hisat2 टूल को gzipped genome index tarball प्रदान करने के लिए `file (params.hisat2_index_zip)` का उपयोग करते हैं।

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.5. वर्कफ़्लो को चलाकर टेस्ट करें कि यह काम करता है

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

तुम `results/align` के अंदर आउटपुट पा सकते हो जैसा कि `HISAT2_ALIGN` प्रोसेस में `publishDir` directive द्वारा निर्दिष्ट किया गया है।

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

यह प्रत्येक सैंपल पर लागू करने के लिए आवश्यक बेसिक प्रोसेसिंग को पूरा करता है।

_हम भाग 2 में MultiQC रिपोर्ट aggregation जोड़ेंगे, जब हम वर्कफ़्लो को एक बार में कई सैंपल स्वीकार करने योग्य बना लेंगे।_

---

### सारांश

तुम जानते हो कि single-end RNAseq सैंपल को व्यक्तिगत रूप से प्रोसेस करने के लिए सभी मुख्य स्टेप्स को कैसे wrap करना है।

### आगे क्या है?

सीखो कि वर्कफ़्लो को कैसे modify करें ताकि कई सैंपल को parallel में प्रोसेस किया जा सके, सभी सैंपल के लिए सभी स्टेप्स में QC रिपोर्ट्स को aggregate किया जा सके, और paired-end RNAseq डेटा पर वर्कफ़्लो चलाना enable किया जा सके।
