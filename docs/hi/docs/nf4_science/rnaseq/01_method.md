# भाग 1: विधि का अवलोकन और मैनुअल परीक्षण

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

बल्क RNAseq डेटा को प्रोसेस और विश्लेषण करने के लिए कई वैध विधियाँ हैं।
इस कोर्स के लिए, हम [Babraham Institute](https://www.babraham.ac.uk/) में डॉ. साइमन एंड्रूज और लौरा बिगिंस द्वारा [यहाँ](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) वर्णित विधि का अनुसरण कर रहे हैं।

हमारा लक्ष्य एक workflow विकसित करना है जो निम्नलिखित प्रोसेसिंग चरणों को लागू करता है: बल्क RNAseq नमूने में reads पर प्रारंभिक गुणवत्ता नियंत्रण चलाना, reads से adapter अनुक्रम ट्रिम करना, reads को reference genome के साथ संरेखित करना, और एक व्यापक गुणवत्ता नियंत्रण (QC) रिपोर्ट तैयार करना।

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/preprocess.svg"
</figure>

- **FASTQC:** ट्रिमिंग से पहले FastQC का उपयोग करके read डेटा पर QC करें
- **TRIM_GALORE:** Trim Galore (Cutadapt और FastQC को बंडल करता है) का उपयोग करके adapter अनुक्रम ट्रिम करें और ट्रिमिंग के बाद QC करें
- **HISAT2_ALIGN:** Hisat2 का उपयोग करके reads को reference genome के साथ संरेखित करें
- **MULTIQC:** MultiQC का उपयोग करके एक व्यापक QC रिपोर्ट तैयार करें

हालांकि, इससे पहले कि हम कोई भी workflow कोड लिखना शुरू करें, हम कुछ टेस्ट डेटा पर commands को मैनुअल रूप से आज़माने जा रहे हैं।
हमें जिन tools की आवश्यकता है वे GitHub Codespaces वातावरण में इंस्टॉल नहीं हैं, इसलिए हम उन्हें containers के माध्यम से उपयोग करेंगे ([Hello Containers](../../hello_nextflow/05_hello_containers.md) देखें)।

!!! note "नोट"

     सुनिश्चित करें कि आप `nf4-science/rnaseq` डायरेक्टरी में हैं। जब आप `pwd` टाइप करते हैं तो दिखाए गए path का अंतिम भाग `rnaseq` होना चाहिए।

---

## 1. प्रारंभिक QC और adapter ट्रिमिंग

हम एक container image पुल करने जा रहे हैं जिसमें `fastqc` और `trim_galore` दोनों इंस्टॉल हैं, इसे इंटरैक्टिव रूप से स्पिन करेंगे और उदाहरण डेटा फ़ाइलों में से एक पर ट्रिमिंग और QC commands चलाएंगे।

### 1.1. कंटेनर पुल करें

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

यह आपको निम्नलिखित console आउटपुट देता है क्योंकि सिस्टम image डाउनलोड करता है:

??? success "कमांड आउटपुट"

    ```console
    0.6.10--1bf8ca4e1967cd18: Pulling from library/trim-galore
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    32ec762be2d0: Pull complete
    d2cb90387285: Pull complete
    Digest: sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907
    Status: Downloaded newer image for community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    ```

### 1.2. कंटेनर को इंटरैक्टिव रूप से स्पिन करें

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

<!--
??? success "कमांड आउटपुट"

    ```console

    ```
-->

आपका prompt कुछ इस तरह `(base) root@b645838b3314:/tmp#` में बदल जाएगा, जो इंगित करता है कि आप अब कंटेनर के अंदर हैं।

कमांड का `-v ./data:/data` भाग हमें कंटेनर के अंदर से `data/` डायरेक्टरी की सामग्री तक पहुँचने में सक्षम बनाएगा।

```bash
ls /data/reads
```

??? success "कमांड आउटपुट"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gzO
    ```

### 1.3. पहली `fastqc` कमांड चलाएं

आइए read डेटा पर गुणवत्ता नियंत्रण मैट्रिक्स एकत्र करने के लिए `fastqc` चलाएं।

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "कमांड आउटपुट"

    ```console
    application/gzip
    Started analysis of ENCSR000COQ1_1.fastq.gz
    Approx 5% complete for ENCSR000COQ1_1.fastq.gz
    Approx 10% complete for ENCSR000COQ1_1.fastq.gz
    Approx 15% complete for ENCSR000COQ1_1.fastq.gz
    Approx 20% complete for ENCSR000COQ1_1.fastq.gz
    Approx 25% complete for ENCSR000COQ1_1.fastq.gz
    Approx 30% complete for ENCSR000COQ1_1.fastq.gz
    Approx 35% complete for ENCSR000COQ1_1.fastq.gz
    Approx 40% complete for ENCSR000COQ1_1.fastq.gz
    Approx 45% complete for ENCSR000COQ1_1.fastq.gz
    Approx 50% complete for ENCSR000COQ1_1.fastq.gz
    Approx 55% complete for ENCSR000COQ1_1.fastq.gz
    Approx 60% complete for ENCSR000COQ1_1.fastq.gz
    Approx 65% complete for ENCSR000COQ1_1.fastq.gz
    Approx 70% complete for ENCSR000COQ1_1.fastq.gz
    Approx 75% complete for ENCSR000COQ1_1.fastq.gz
    Approx 80% complete for ENCSR000COQ1_1.fastq.gz
    Approx 85% complete for ENCSR000COQ1_1.fastq.gz
    Approx 90% complete for ENCSR000COQ1_1.fastq.gz
    Approx 95% complete for ENCSR000COQ1_1.fastq.gz
    Analysis complete for ENCSR000COQ1_1.fastq.gz
    ```

यह बहुत जल्दी चलना चाहिए।
आप मूल डेटा के समान डायरेक्टरी में आउटपुट फ़ाइलें पा सकते हैं:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

<!-- switch to tree -->

```console title="आउटपुट"
/data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
```

### 1.4. `trim_galore` के साथ adapter अनुक्रम ट्रिम करें

अब आइए `trim_galore` चलाएं, जो Cutadapt और FastQC को बंडल करता है, adapter अनुक्रमों को ट्रिम करने और पोस्ट-ट्रिमिंग QC मैट्रिक्स एकत्र करने के लिए।

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

`--fastqc` flag कमांड को ट्रिमिंग पूर्ण होने के बाद स्वचालित रूप से एक QC संग्रह चरण चलाने का कारण बनता है।

_आउटपुट बहुत verbose है इसलिए निम्नलिखित संक्षिप्त है।_

??? success "कमांड आउटपुट"

    ```console
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    <...>

    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

आप working डायरेक्टरी में आउटपुट फ़ाइलें पा सकते हैं:

```bash
ls ENCSR000COQ1_1*
```

```console title="आउटपुट"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
```

### 1.5. आउटपुट फ़ाइलों को कंटेनर के बाहर फाइलसिस्टम में ले जाएं

कंटेनर के अंदर जो कुछ भी रहता है वह भविष्य के काम के लिए दुर्गम होगा इसलिए आइए इन्हें एक नई डायरेक्टरी में ले जाएं।

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

### 1.6. कंटेनर से बाहर निकलें

```bash
exit
```

---

## 2. reads को reference genome के साथ संरेखित करें

हम एक container image पुल करने जा रहे हैं जिसमें `hisat2` इंस्टॉल है, इसे इंटरैक्टिव रूप से स्पिन करेंगे और RNAseq डेटा को reference genome के साथ संरेखित करने के लिए संरेखण कमांड चलाएंगे।

### 2.1. `hisat2` कंटेनर पुल करें

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "कमांड आउटपुट"

    ```console
    Unable to find image 'community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e' locally
    5e49f68a37dc010e: Pulling from library/hisat2_samtools
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    e74ed5dd390b: Pull complete
    abfcf0185e51: Pull complete
    Digest: sha256:29d8e1a3172a2bdde7be813f7ebec22d331388194a7c0de872b4ccca4bed8f45
    Status: Downloaded newer image for community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
    ```

### 2.2. `hisat2` कंटेनर को इंटरैक्टिव रूप से स्पिन करें

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

कमांड पहले के समान है, संबंधित container URI को स्वैप कर दिया गया है।

### 2.3. Hisat2 genome index फ़ाइलें बनाएं

Hisat2 को genome संदर्भ को बहुत विशिष्ट प्रारूप में प्रदान किए जाने की आवश्यकता होती है, और यह सिर्फ हमारे द्वारा प्रदान की गई `genome.fa` FASTA फ़ाइल का उपभोग नहीं कर सकता है, इसलिए हम प्रासंगिक संसाधन बनाने के लिए इस अवसर का लाभ उठाने जा रहे हैं।

```bash
hisat2-build /data/genome.fa genome_index
```

आउटपुट बहुत verbose है इसलिए निम्नलिखित संक्षिप्त है:

<!-- TODO: switch to full output -->

??? success "कमांड आउटपुट"

    ```console
    Settings:
      Output files: "genome_index.*.ht2"
    <...>
    Total time for call to driver() for forward index: 00:00:16
    ```

यह कई genome index फ़ाइलें बनाता है, जिन्हें आप working डायरेक्टरी में पा सकते हैं।

```bash
ls genome_index.*
```

```console title="आउटपुट"
genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
```

हम इनका उपयोग एक क्षण में करेंगे, लेकिन पहले आइए इन genome index फ़ाइलों के साथ एक gzipped tarball बनाएं; हमें बाद में इनकी आवश्यकता होगी और इन्हें generate करना आम तौर पर कुछ ऐसा नहीं है जो हम workflow के हिस्से के रूप में करना चाहते हैं।

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

यह genome index फ़ाइलों वाला एक `genome_index.tar.gz` tarball हमारे फ़ाइल सिस्टम पर `data/` डायरेक्टरी में स्टोर करता है, जो इस कोर्स के भाग 2 में काम आएगा।

### 2.4. `hisat2` कमांड चलाएं

अब हम संरेखण कमांड चला सकते हैं, जो `hisat2` के साथ संरेखण चरण करता है फिर आउटपुट को BAM फ़ाइल के रूप में लिखने के लिए `samtools` को पाइप करता है।

read डेटा इनपुट `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` फ़ाइल है जिसे हमने पिछले चरण में `trim_galore` के साथ उत्पन्न किया था।

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "कमांड आउटपुट"

    ```console
    HISAT2 summary stats:
            Total reads: 27816
                    Aligned 0 time: 1550 (5.57%)
                    Aligned 1 time: 25410 (91.35%)
                    Aligned >1 times: 856 (3.08%)
            Overall alignment rate: 94.43%
    ```

यह लगभग तुरंत चलता है क्योंकि यह एक बहुत छोटी टेस्ट फ़ाइल है।
वास्तविक स्केल पर यह बहुत अधिक समय ले सकता है।

एक बार फिर आप working डायरेक्टरी में आउटपुट फ़ाइलें पा सकते हैं:

```bash
ls ENCSR000COQ1_1*
```

```console title="आउटपुट"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

### 2.5. आउटपुट फ़ाइलों को कंटेनर के बाहर फाइलसिस्टम में ले जाएं

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

### 2.6. कंटेनर से बाहर निकलें

```bash
exit
```

---

## 3. एक व्यापक QC रिपोर्ट तैयार करें

हम एक container image पुल करने जा रहे हैं जिसमें `multiqc` इंस्टॉल है, इसे इंटरैक्टिव रूप से स्पिन करेंगे और before/after FastQC रिपोर्ट फ़ाइलों पर एक रिपोर्ट generation कमांड चलाएंगे।

### 3.1. `multiqc` कंटेनर पुल करें

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "कमांड आउटपुट"

    ```console
    ad8f247edb55897c: Pulling from library/pip_multiqc
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    3f229294c69a: Pull complete
    5a5ad47fd84c: Pull complete
    Digest: sha256:0ebb1d9605395a7df49ad0eb366b21f46afd96a5090376b0d8941cf5294a895a
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

### 3.2. `multiqc` कंटेनर को इंटरैक्टिव रूप से स्पिन करें

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

### 3.3. `multiqc` कमांड चलाएं

```bash
multiqc /data/reads /data/trimmed /data/aligned -n ENCSR000COQ1_1_QC
```

??? success "कमांड आउटपुट"

    ```console

    /// MultiQC 🔍 v1.27.1

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 20/20
                hisat2 | Found 1 reports
              cutadapt | Found 1 reports
                fastqc | Found 1 reports
        write_results | Data        : ENCSR000COQ1_1_QC_data
        write_results | Report      : ENCSR000COQ1_1_QC.html
              multiqc | MultiQC complete
    ```

MultiQC संगत QC रिपोर्ट के लिए डायरेक्टरी में खोज करने में सक्षम है और जो कुछ भी मिलता है उसे एकत्र करेगा।

यहाँ हम देखते हैं कि टूल ने हमारे द्वारा उत्पन्न तीनों QC रिपोर्ट पाईं: प्रारंभिक QC जो हमने `fastqc` के साथ की, पोस्ट-ट्रिमिंग रिपोर्ट `cutadapt` से (`trim_galore` के माध्यम से बनाई गई) और `hisat2` द्वारा उत्पादित पोस्ट-संरेखण QC।

आउटपुट फ़ाइलें एक बार फिर working डायरेक्टरी में हैं:

```bash
ls ENCSR000COQ1_1_QC*
```

```console title="आउटपुट"
ENCSR000COQ1_1_QC.html

ENCSR000COQ1_1_QC_data:
cutadapt_filtered_reads_plot.txt                     fastqc_top_overrepresented_sequences_table.txt
cutadapt_trimmed_sequences_plot_3_Counts.txt         hisat2_se_plot.txt
cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc.log
fastqc-status-check-heatmap.txt                      multiqc_citations.txt
fastqc_adapter_content_plot.txt                      multiqc_cutadapt.txt
fastqc_per_base_n_content_plot.txt                   multiqc_data.json
fastqc_per_base_sequence_quality_plot.txt            multiqc_fastqc.txt
fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_general_stats.txt
fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_hisat2.txt
fastqc_per_sequence_quality_scores_plot.txt          multiqc_software_versions.txt
fastqc_sequence_counts_plot.txt                      multiqc_sources.txt
fastqc_sequence_duplication_levels_plot.txt
```

### 3.4. आउटपुट फ़ाइलों को कंटेनर के बाहर फाइलसिस्टम में ले जाएं

```bash
mkdir /data/final_qc
mv ENCSR000COQ1_1_QC** /data/final_qc
```

### 3.5. कंटेनर से बाहर निकलें

```bash
exit
```

---

### निष्कर्ष

आपने संबंधित containers में सभी व्यक्तिगत commands को इंटरैक्टिव रूप से परीक्षण किया है।

### आगे क्या है?

सीखें कि उन्हीं commands को एक बहु-चरणीय workflow में कैसे लपेटा जाए जो काम को निष्पादित करने के लिए containers का उपयोग करता है।
