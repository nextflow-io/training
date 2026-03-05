# भाग 1: विधि का अवलोकन

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

बल्क RNAseq डेटा को प्रोसेस और विश्लेषण करने के लिए कई वैध विधियाँ हैं।
इस कोर्स के लिए, हम [Babraham Institute](https://www.babraham.ac.uk/) में डॉ. साइमन एंड्रूज और लौरा बिगिंस द्वारा [यहाँ](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) वर्णित विधि का अनुसरण कर रहे हैं।

हमारा लक्ष्य एक workflow विकसित करना है जो निम्नलिखित प्रोसेसिंग चरणों को लागू करता है: बल्क RNAseq नमूने में reads पर प्रारंभिक गुणवत्ता नियंत्रण चलाना, reads से adapter अनुक्रम ट्रिम करना, reads को reference genome के साथ संरेखित करना, और एक व्यापक गुणवत्ता नियंत्रण (QC) रिपोर्ट तैयार करना।

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

- **FASTQC:** ट्रिमिंग से पहले FastQC का उपयोग करके read डेटा पर QC करें
- **TRIM_GALORE:** Trim Galore (Cutadapt और FastQC को बंडल करता है) का उपयोग करके adapter अनुक्रम ट्रिम करें और ट्रिमिंग के बाद QC करें
- **HISAT2_ALIGN:** Hisat2 का उपयोग करके reads को reference genome के साथ संरेखित करें
- **MULTIQC:** MultiQC का उपयोग करके एक व्यापक QC रिपोर्ट तैयार करें

### विधियाँ

हम तुम्हें दिखाने जा रहे हैं कि इन प्रोसेसिंग चरणों को दो चरणों में कैसे लागू किया जाए।
पहले हम **सिंगल-सैंपल प्रोसेसिंग** से शुरू करेंगे जो एक नमूने पर QC, ट्रिमिंग और संरेखण टूल चलाता है।
फिर हम **मल्टी-सैंपल प्रोसेसिंग** तक विस्तार करेंगे जो कई नमूनों पर समान टूल चलाता है और एक एकत्रित गुणवत्ता नियंत्रण रिपोर्ट तैयार करता है।

इससे पहले कि हम किसी भी दृष्टिकोण के लिए कोई भी workflow कोड लिखना शुरू करें, हम कुछ टेस्ट डेटा पर commands को मैनुअल रूप से आज़माने जा रहे हैं।

### डेटासेट

हम निम्नलिखित डेटा और संबंधित संसाधन प्रदान करते हैं:

- **RNAseq डेटा** (`reads/`): छह नमूनों से FASTQ फ़ाइलें, फ़ाइल आकार कम रखने के लिए एक छोटे क्षेत्र तक सीमित। प्रत्येक नमूने में paired-end reads हैं (प्रति नमूना दो फ़ाइलें), हालांकि हम केवल single-end reads के साथ काम करना शुरू करते हैं।
- **एक reference genome** (`genome.fa`): मानव क्रोमोसोम 20 का एक छोटा क्षेत्र (hg19/b37 से)।
- **CSV samplesheets** (`single-end.csv` और `paired-end.csv`): उदाहरण डेटा फ़ाइलों के IDs और paths को सूचीबद्ध करने वाली फ़ाइलें।

### सॉफ़्टवेयर

शामिल चार मुख्य टूल हैं गुणवत्ता नियंत्रण मैट्रिक्स संग्रह के लिए [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), adapter ट्रिमिंग के लिए [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) (पोस्ट-ट्रिमिंग QC के लिए Cutadapt और FastQC को बंडल करता है), reference genome के साथ spliced संरेखण के लिए [HISAT2](http://daehwankimlab.github.io/hisat2/), और एकत्रित QC रिपोर्ट generation के लिए [MultiQC](https://multiqc.info/)।

ये टूल GitHub Codespaces वातावरण में इंस्टॉल नहीं हैं, इसलिए हम उन्हें Seqera Containers सेवा के माध्यम से प्राप्त containers के माध्यम से उपयोग करेंगे ([Hello Containers](../../hello_nextflow/05_hello_containers.md) देखें)।

!!! tip "सुझाव"

     सुनिश्चित करें कि तुम `nf4-science/rnaseq` डायरेक्टरी में हो। जब तुम `pwd` टाइप करते हो तो दिखाए गए path का अंतिम भाग `rnaseq` होना चाहिए।

---

## 1. सिंगल-सैंपल प्रोसेसिंग

इस सेक्शन में हम उन commands का परीक्षण करते हैं जो एक single RNAseq नमूने को प्रोसेस करते हैं: गुणवत्ता नियंत्रण, adapter ट्रिमिंग, और reference genome के साथ संरेखण।
ये वे commands हैं जिन्हें हम इस कोर्स के भाग 2 में एक Nextflow workflow में लपेटेंगे।

1. FastQC का उपयोग करके FASTQ फ़ाइल पर प्रारंभिक QC चलाएं
2. Trim Galore का उपयोग करके adapter अनुक्रम ट्रिम करें और पोस्ट-ट्रिमिंग QC चलाएं
3. HISAT2 का उपयोग करके ट्रिम किए गए reads को reference genome के साथ संरेखित करें

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-02.svg"
</figure>

हम केवल एक नमूने पर इन commands का परीक्षण करके शुरू करते हैं।

### 1.1. QC और adapter ट्रिमिंग

पहले, हम उदाहरण डेटा फ़ाइलों में से एक पर QC और ट्रिमिंग commands चलाना चाहते हैं।

#### 1.1.1. कंटेनर पुल करें

आइए एक container image पुल करें जिसमें `fastqc` और `trim_galore` दोनों इंस्टॉल हैं:

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

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

यदि तुमने पहले इस image को डाउनलोड नहीं किया है, तो इसे पूरा होने में एक मिनट लग सकता है।
एक बार यह हो जाने के बाद, तुम्हारे पास container image की एक लोकल कॉपी है।

#### 1.1.2. कंटेनर को इंटरैक्टिव रूप से स्पिन करें

कंटेनर को इंटरैक्टिव रूप से चलाने के लिए, `-it` flags के साथ `docker run` का उपयोग करें।
`-v ./data:/data` विकल्प हमारी लोकल `data/` डायरेक्टरी को माउंट करता है ताकि हम कंटेनर के अंदर से इनपुट फ़ाइलों तक पहुँच सकें।

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "कमांड आउटपुट"

    ```console
    (base) root@b645838b3314:/tmp#
    ```

तुम्हारा prompt कुछ इस तरह `(base) root@b645838b3314:/tmp#` में बदल जाएगा, जो इंगित करता है कि तुम अब कंटेनर के अंदर हो।

सत्यापित करें कि तुम `/data/reads` के अंतर्गत sequence डेटा फ़ाइलें देख सकते हो:

```bash
ls /data/reads
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gz
    ```

इसके साथ, तुम अपनी पहली कमांड आज़माने के लिए तैयार हो।

#### 1.1.3. FastQC कमांड चलाएं

ऊपर संदर्भित विधि हमें एक single फ़ाइल पर QC चलाने के लिए command line देती है।
हमें केवल इनपुट फ़ाइल प्रदान करने की आवश्यकता है; टूल स्वचालित रूप से मूल डेटा के समान डायरेक्टरी में आउटपुट फ़ाइलें generate करेगा।

एक डेटा फ़ाइल पर `fastqc` कमांड चलाएं:

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
तुम मूल डेटा के समान डायरेक्टरी में आउटपुट फ़ाइलें पा सकते हो:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    /data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
    ```

तुम्हें एक HTML रिपोर्ट और QC मैट्रिक्स वाला एक ZIP archive दिखना चाहिए।
यह पहले चरण के परीक्षण को पूरा करता है।

#### 1.1.4. Trim Galore के साथ adapter अनुक्रम ट्रिम करें

अब आइए `trim_galore` चलाएं, जो Cutadapt और FastQC को बंडल करता है, adapter अनुक्रमों को ट्रिम करने और पोस्ट-ट्रिमिंग QC मैट्रिक्स एकत्र करने के लिए।
जैसा कि ऊपर बताया गया है, सॉफ़्टवेयर उसी कंटेनर में शामिल है, इसलिए वहाँ कोई बदलाव की आवश्यकता नहीं है।

कमांड सीधी है; हमें बस ट्रिमिंग पूर्ण होने के बाद स्वचालित रूप से एक QC संग्रह चरण चलाने के लिए `--fastqc` flag जोड़ने की आवश्यकता है।

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "कमांड आउटपुट"

    ```console hl_lines="54 55 56 58 59 60"
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)



    AUTO-DETECTING ADAPTER TYPE
    ===========================
    Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> /data/reads/ENCSR000COQ1_1.fastq.gz <<)

    Found perfect matches for the following adapter sequences:
    Adapter type	Count	Sequence	Sequences analysed	Percentage
    Illumina	9	AGATCGGAAGAGC	27816	0.03
    smallRNA	0	TGGAATTCTCGG	27816	0.00
    Nextera	0	CTGTCTCTTATA	27816	0.00
    Using Illumina adapter for trimming (count: 9). Second best hit was smallRNA (count: 0)

    Writing report to 'ENCSR000COQ1_1.fastq.gz_trimming_report.txt'

    SUMMARISING RUN PARAMETERS
    ==========================
    Input filename: /data/reads/ENCSR000COQ1_1.fastq.gz
    Trimming mode: single-end
    Trim Galore version: 0.6.10
    Cutadapt version: 4.9
    Number of cores used for trimming: 1
    Quality Phred score cutoff: 20
    Quality encoding type selected: ASCII+33
    Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    Maximum trimming error rate: 0.1 (default)
    Minimum required adapter overlap (stringency): 1 bp
    Minimum required sequence length before a sequence gets removed: 20 bp
    Running FastQC on the data once trimming has completed
    Output file(s) will be GZIP compressed

    Cutadapt seems to be fairly up-to-date (version 4.9). Setting -j 1
    Writing final adapter and quality trimmed output to ENCSR000COQ1_1_trimmed.fq.gz


      >>> Now performing quality (cutoff '-q 20') and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file /data/reads/ENCSR000COQ1_1.fastq.gz <<<
    This is cutadapt 4.9 with Python 3.12.7
    Command line parameters: -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC /data/reads/ENCSR000COQ1_1.fastq.gz
    Processing single-end reads on 1 core ...
    Finished in 0.373 s (13.399 µs/read; 4.48 M reads/minute).

    === Summary ===

    Total reads processed:                  27,816
    Reads with adapters:                     9,173 (33.0%)
    Reads written (passing filters):        27,816 (100.0%)

    Total basepairs processed:     2,114,016 bp
    Quality-trimmed:                       0 bp (0.0%)
    Total written (filtered):      2,100,697 bp (99.4%)

    === Adapter 1 ===

    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 9173 times

    Minimum overlap: 1
    No. of allowed errors:
    1-9 bp: 0; 10-13 bp: 1

    Bases preceding removed adapters:
      A: 27.4%
      C: 37.4%
      G: 20.9%
      T: 14.3%
      none/other: 0.0%

    Overview of removed sequences
    length	count	expect	max.err	error counts
    1	6229	6954.0	0	6229
    2	2221	1738.5	0	2221
    3	581	434.6	0	581
    4	88	108.7	0	88
    5	33	27.2	0	33
    6	2	6.8	0	2
    7	1	1.7	0	1
    9	1	0.1	0	1
    10	2	0.0	1	2
    12	1	0.0	1	0 1
    14	4	0.0	1	3 1
    16	1	0.0	1	1
    19	1	0.0	1	1
    22	1	0.0	1	1
    29	4	0.0	1	0 4
    33	3	0.0	1	3

    RUN STATISTICS FOR INPUT FILE: /data/reads/ENCSR000COQ1_1.fastq.gz
    =============================================
    27816 sequences processed in total
    Sequences removed because they became shorter than the length cutoff of 20 bp:	0 (0.0%)


      >>> Now running FastQC on the data <<<

    application/gzip
    Started analysis of ENCSR000COQ1_1_trimmed.fq.gz
    Approx 5% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 10% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 15% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 20% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 25% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 30% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 35% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 40% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 45% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 50% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 55% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 60% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 65% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 70% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 75% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 80% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 85% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 90% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 95% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

आउटपुट बहुत verbose है, इसलिए हमने ऊपर दिए गए उदाहरण में सबसे प्रासंगिक लाइनों को हाइलाइट किया है।
तुम working डायरेक्टरी में आउटपुट फ़ाइलें पा सकते हो:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
    ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

इसमें ट्रिम किए गए reads, ट्रिमिंग रिपोर्ट और पोस्ट-ट्रिमिंग QC फ़ाइलें शामिल हैं।

#### 1.1.5. आउटपुट फ़ाइलों को ले जाएं

कंटेनर के अंदर जो कुछ भी रहता है वह भविष्य के काम के लिए दुर्गम होगा, इसलिए हमें इन फ़ाइलों को माउंट किए गए फाइलसिस्टम पर एक डायरेक्टरी में ले जाने की आवश्यकता है।

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    └── ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

फ़ाइलें अब तुम्हारे सामान्य फाइलसिस्टम में accessible हैं।

#### 1.1.6. कंटेनर से बाहर निकलें

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करें।

```bash
exit
```

तुम्हारा prompt सामान्य हो जाना चाहिए; यह पहले दो चरणों के परीक्षण को पूरा करता है।

### 1.2. reads को reference genome के साथ संरेखित करें

अगला, हम ट्रिम किए गए RNAseq reads को reference genome के साथ संरेखित करने के लिए संरेखण कमांड चलाना चाहते हैं।

#### 1.2.1. कंटेनर पुल करें

आइए एक container image पुल करें जिसमें `hisat2` और `samtools` इंस्टॉल हैं:

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

तुम देखोगे कि कुछ layers `Already exists` दिखाते हैं क्योंकि वे पहले हमने जो Trim Galore container image पुल की थी उसके साथ साझा किए गए हैं।
परिणामस्वरूप, यह pull पहले वाले की तुलना में तेज़ होना चाहिए।

#### 1.2.2. कंटेनर को इंटरैक्टिव रूप से स्पिन करें

कंटेनर को इंटरैक्टिव रूप से स्पिन करें, पहले की तरह समान दृष्टिकोण का उपयोग करते हुए संबंधित container URI को स्वैप करके।

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

तुम्हारा prompt फिर से बदल जाएगा यह इंगित करने के लिए कि तुम कंटेनर के अंदर हो।

#### 1.2.3. genome index फ़ाइलें बनाएं

HISAT2 को genome संदर्भ को बहुत विशिष्ट प्रारूप में प्रदान किए जाने की आवश्यकता होती है, और यह सिर्फ हमारे द्वारा प्रदान की गई `genome.fa` FASTA फ़ाइल का उपभोग नहीं कर सकता है, इसलिए हम प्रासंगिक संसाधन बनाने के लिए इस अवसर का लाभ उठाने जा रहे हैं।

```bash
hisat2-build /data/genome.fa genome_index
```

??? success "कमांड आउटपुट"

    ```console hl_lines="1 2 218"
    Settings:
      Output files: "genome_index.*.ht2"
      Line rate: 6 (line is 64 bytes)
      Lines per side: 1 (side is 64 bytes)
      Offset rate: 4 (one in 16)
      FTable chars: 10
      Strings: unpacked
      Local offset rate: 3 (one in 8)
      Local fTable chars: 6
      Local sequence length: 57344
      Local sequence overlap between two consecutive indexes: 1024
      Endianness: little
      Actual local endianness: little
      Sanity checking: disabled
      Assertions: disabled
      Random seed: 0
      Sizeofs: void*:8, int:4, long:8, size_t:8
    Input files DNA, FASTA:
      /data/genome.fa
    Reading reference sizes
      Time reading reference sizes: 00:00:00
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:00
      Time to read SNPs and splice sites: 00:00:00
    Using parameters --bmax 6542727 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 6542727 --dcv 1024
    Constructing suffix-array element generator
    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:00:01
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:00
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:00:00
      Sanity-checking and returning
    Building samples
    Reserving space for 12 sample suffixes
    Generating random suffixes
    QSorting 12 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 12 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 1, merged 7; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 4.98493e+06 (target: 6542726)
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering GFM loop
    Getting block 1 of 7
      Reserving size (6542727) for bucket 1
      Calculating Z arrays for bucket 1
      Entering block accumulator loop for bucket 1:
      bucket 1: 10%
      bucket 1: 20%
      bucket 1: 30%
      bucket 1: 40%
      bucket 1: 50%
      bucket 1: 60%
      bucket 1: 70%
      bucket 1: 80%
      bucket 1: 90%
      bucket 1: 100%
      Sorting block of length 3540952 for bucket 1
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 3540953 for bucket 1
    Getting block 2 of 7
      Reserving size (6542727) for bucket 2
      Calculating Z arrays for bucket 2
      Entering block accumulator loop for bucket 2:
      bucket 2: 10%
      bucket 2: 20%
      bucket 2: 30%
      bucket 2: 40%
      bucket 2: 50%
      bucket 2: 60%
      bucket 2: 70%
      bucket 2: 80%
      bucket 2: 90%
      bucket 2: 100%
      Sorting block of length 6195795 for bucket 2
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6195796 for bucket 2
    Getting block 3 of 7
      Reserving size (6542727) for bucket 3
      Calculating Z arrays for bucket 3
      Entering block accumulator loop for bucket 3:
      bucket 3: 10%
      bucket 3: 20%
      bucket 3: 30%
      bucket 3: 40%
      bucket 3: 50%
      bucket 3: 60%
      bucket 3: 70%
      bucket 3: 80%
      bucket 3: 90%
      bucket 3: 100%
      Sorting block of length 6199288 for bucket 3
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6199289 for bucket 3
    Getting block 4 of 7
      Reserving size (6542727) for bucket 4
      Calculating Z arrays for bucket 4
      Entering block accumulator loop for bucket 4:
      bucket 4: 10%
      bucket 4: 20%
      bucket 4: 30%
      bucket 4: 40%
      bucket 4: 50%
      bucket 4: 60%
      bucket 4: 70%
      bucket 4: 80%
      bucket 4: 90%
      bucket 4: 100%
      Sorting block of length 6454986 for bucket 4
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 6454987 for bucket 4
    Getting block 5 of 7
      Reserving size (6542727) for bucket 5
      Calculating Z arrays for bucket 5
      Entering block accumulator loop for bucket 5:
      bucket 5: 10%
      bucket 5: 20%
      bucket 5: 30%
      bucket 5: 40%
      bucket 5: 50%
      bucket 5: 60%
      bucket 5: 70%
      bucket 5: 80%
      bucket 5: 90%
      bucket 5: 100%
      Sorting block of length 3493181 for bucket 5
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3493182 for bucket 5
    Getting block 6 of 7
      Reserving size (6542727) for bucket 6
      Calculating Z arrays for bucket 6
      Entering block accumulator loop for bucket 6:
      bucket 6: 10%
      bucket 6: 20%
      bucket 6: 30%
      bucket 6: 40%
      bucket 6: 50%
      bucket 6: 60%
      bucket 6: 70%
      bucket 6: 80%
      bucket 6: 90%
      bucket 6: 100%
      Sorting block of length 5875908 for bucket 6
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 5875909 for bucket 6
    Getting block 7 of 7
      Reserving size (6542727) for bucket 7
      Calculating Z arrays for bucket 7
      Entering block accumulator loop for bucket 7:
      bucket 7: 10%
      bucket 7: 20%
      bucket 7: 30%
      bucket 7: 40%
      bucket 7: 50%
      bucket 7: 60%
      bucket 7: 70%
      bucket 7: 80%
      bucket 7: 90%
      bucket 7: 100%
      Sorting block of length 3134429 for bucket 7
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3134430 for bucket 7
    Exited GFM loop
    fchr[A]: 0
    fchr[C]: 9094775
    fchr[G]: 17470759
    fchr[T]: 25839994
    fchr[$]: 34894545
    Exiting GFM::buildToDisk()
    Returning from initFromVector
    Wrote 15826295 bytes to primary GFM file: genome_index.1.ht2
    Wrote 8723644 bytes to secondary GFM file: genome_index.2.ht2
    Re-opening _in1 and _in2 as input streams
    Returning from GFM constructor
    Returning from initFromVector
    Wrote 15353415 bytes to primary GFM file: genome_index.5.ht2
    Wrote 8883598 bytes to secondary GFM file: genome_index.6.ht2
    Re-opening _in5 and _in5 as input streams
    Returning from HGFM constructor
    Headers:
        len: 34894545
        gbwtLen: 34894546
        nodes: 34894546
        sz: 8723637
        gbwtSz: 8723637
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 0
        eftabSz: 0
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 2180910
        offsSz: 8723640
        lineSz: 64
        sideSz: 64
        sideGbwtSz: 48
        sideGbwtLen: 192
        numSides: 181743
        numLines: 181743
        gbwtTotLen: 11631552
        gbwtTotSz: 11631552
        reverse: 0
        linearFM: Yes
    Total time for call to driver() for forward index: 00:00:12
    ```

आउटपुट बहुत verbose है, इसलिए हमने ऊपर दिए गए उदाहरण में कुछ प्रासंगिक लाइनों को हाइलाइट किया है।

यह कई genome index फ़ाइलें बनाता है, जिन्हें तुम working डायरेक्टरी में पा सकते हो।

```bash
ls genome_index.*
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
    genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
    ```

हमें बाद में इन फ़ाइलों की आवश्यकता होगी, और इन्हें generate करना आम तौर पर कुछ ऐसा नहीं है जो हम workflow के हिस्से के रूप में करना चाहते हैं, इसलिए हम genome index फ़ाइलों वाला एक gzipped tarball generate करने जा रहे हैं जिसे हम आवश्यकतानुसार आसानी से पास कर सकते हैं।

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

??? success "कमांड आउटपुट"

    ```console
    genome_index.1.ht2
    genome_index.2.ht2
    genome_index.3.ht2
    genome_index.4.ht2
    genome_index.5.ht2
    genome_index.6.ht2
    genome_index.7.ht2
    genome_index.8.ht2
    ```

हम परिणामी `genome_index.tar.gz` tarball को कुछ मिनटों में हमारे फाइलसिस्टम पर `data/` डायरेक्टरी में ले जाएंगे।
यह इस कोर्स के भाग 2 में काम आएगा।

#### 1.2.4. संरेखण कमांड चलाएं

अब हम संरेखण कमांड चला सकते हैं, जो `hisat2` के साथ संरेखण चरण करता है फिर आउटपुट को BAM फ़ाइल के रूप में लिखने के लिए `samtools` को पाइप करता है।

read डेटा इनपुट `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` फ़ाइल है जिसे हमने पिछले चरण में `trim_galore` के साथ generate किया था।

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 27816
    		Aligned 0 time: 1550 (5.57%)
    		Aligned 1 time: 25410 (91.35%)
    		Aligned >1 times: 856 (3.08%)
    	Overall alignment rate: 94.43%
    ```

यह लगभग तुरंत चलता है क्योंकि यह एक बहुत छोटी टेस्ट फ़ाइल है।
पूर्ण स्केल पर यह बहुत अधिक समय ले सकता है।

एक बार फिर तुम working डायरेक्टरी में आउटपुट फ़ाइलें पा सकते हो:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "डायरेक्टरी सामग्री"

    ```console title="आउटपुट"
    ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
    ```

संरेखण ने एक BAM फ़ाइल और संरेखण statistics के साथ एक log फ़ाइल उत्पन्न की।

#### 1.2.5. आउटपुट फ़ाइलों को ले जाएं

पहले की तरह, आउटपुट फ़ाइलों को माउंट किए गए फाइलसिस्टम पर एक डायरेक्टरी में ले जाएं ताकि वे कंटेनर से बाहर निकलने के बाद accessible रहें।

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

इसके साथ, हमारे पास वह सब कुछ है जो हमें चाहिए।

#### 1.2.6. कंटेनर से बाहर निकलें

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करें।

```bash
exit
```

तुम्हारा prompt सामान्य हो जाना चाहिए।
यह सिंगल-सैंपल प्रोसेसिंग परीक्षण रन को समाप्त करता है।

!!! example "इसे एक workflow के रूप में लिखें!"

    यदि तुम इस विश्लेषण को Nextflow workflow के रूप में लागू करना शुरू करना चाहते हो तो तुम तुरंत [भाग 2](./02_single-sample.md) पर जा सकते हो।
    तुम्हें बस भाग 3 पर जाने से पहले परीक्षण के दूसरे दौर को पूरा करने के लिए वापस आना होगा।

---

## 2. मल्टी-सैंपल QC एकत्रीकरण

हमने अभी जो commands का परीक्षण किया वे एक समय में एक नमूने को प्रोसेस करते हैं।
व्यवहार में, हमें आम तौर पर कई नमूनों को प्रोसेस करने और फिर समग्र डेटासेट की गुणवत्ता का मूल्यांकन करने के लिए उन सभी में QC परिणामों को एकत्र करने की आवश्यकता होती है।

[MultiQC](https://multiqc.info/) एक टूल है जो कई सामान्य bioinformatics टूल से QC रिपोर्ट के लिए डायरेक्टरी में खोज करता है और उन्हें एक single व्यापक HTML रिपोर्ट में एकत्र करता है।
यह FastQC, Cutadapt (Trim Galore के माध्यम से) और HISAT2 से आउटपुट को पहचान सकता है, कई अन्य के बीच।

यहाँ हम समान per-sample टूल के माध्यम से दो अतिरिक्त नमूनों को प्रोसेस करते हैं, फिर सभी तीन नमूनों में QC रिपोर्ट को एकत्र करने के लिए MultiQC का उपयोग करते हैं।
ये वे commands हैं जिन्हें हम इस कोर्स के भाग 3 में एक Nextflow workflow में लपेटेंगे।

1. Trim Galore का उपयोग करके अतिरिक्त नमूनों पर QC और ट्रिमिंग चलाएं
2. HISAT2 का उपयोग करके अतिरिक्त नमूनों पर संरेखण चलाएं
3. MultiQC का उपयोग करके सभी QC रिपोर्ट को एक व्यापक रिपोर्ट में एकत्र करें

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

### 2.1. अतिरिक्त नमूनों पर QC और ट्रिम करें

per-sample QC और ट्रिमिंग commands वही हैं जो हमने सेक्शन 1.1 में चलाई थीं।
हमने पहले ही container image पुल कर ली है, इसलिए हम इसे सीधे स्पिन कर सकते हैं।

#### 2.1.1. कंटेनर स्पिन करें

हमने पहले ही सेक्शन 1.1 में इस container image को पुल कर लिया था, इसलिए हम इसे सीधे स्पिन कर सकते हैं:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

तुम्हारा prompt बदल जाता है यह इंगित करने के लिए कि तुम कंटेनर के अंदर हो।

#### 2.1.2. अतिरिक्त नमूनों पर QC और ट्रिमिंग चलाएं

दो और नमूनों पर FastQC और Trim Galore चलाएं, एक के बाद एक।

```bash
fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
trim_galore --fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

एक बार यह पूरा हो जाने पर, तुम्हारे पास working डायरेक्टरी में दोनों नमूनों के लिए Trim Galore आउटपुट फ़ाइलें होनी चाहिए।

#### 2.1.3. आउटपुट फ़ाइलों को ले जाएं

Trim Galore आउटपुट फ़ाइलों को उसी डायरेक्टरी में ले जाएं जिसका हमने सेक्शन 1 में उपयोग किया था।

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/trimmed
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    ├── ENCSR000COQ1_1_trimmed_fastqc.zip
    ├── ENCSR000COQ2_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ2_1_trimmed.fq.gz
    ├── ENCSR000COQ2_1_trimmed_fastqc.html
    ├── ENCSR000COQ2_1_trimmed_fastqc.zip
    ├── ENCSR000COR1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COR1_1_trimmed.fq.gz
    ├── ENCSR000COR1_1_trimmed_fastqc.html
    └── ENCSR000COR1_1_trimmed_fastqc.zip
    ```

फ़ाइलें अब तुम्हारे सामान्य फाइलसिस्टम में accessible हैं।

#### 2.1.4. कंटेनर से बाहर निकलें

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करें।

```bash
exit
```

तुम्हारा prompt सामान्य हो जाना चाहिए।

### 2.2. अतिरिक्त नमूनों को संरेखित करें

संरेखण commands वही हैं जो हमने सेक्शन 1.2 में चलाई थीं।
हमें पहले सहेजे गए tarball से genome index को extract करने की आवश्यकता है, क्योंकि मूल index फ़ाइलें एक कंटेनर के अंदर बनाई गई थीं जो अब मौजूद नहीं है।

#### 2.2.1. कंटेनर स्पिन करें

हमने पहले ही सेक्शन 1.2 में इस container image को पुल कर लिया था, इसलिए हम इसे सीधे स्पिन कर सकते हैं:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

तुम्हारा prompt बदल जाता है यह इंगित करने के लिए कि तुम कंटेनर के अंदर हो।

#### 2.2.2. genome index extract करें

हमने माउंट किए गए फाइलसिस्टम में सहेजे गए tarball से genome index फ़ाइलों को extract करें:

```bash
tar -xzf /data/genome_index.tar.gz
```

यह working डायरेक्टरी में `genome_index.*` फ़ाइलों को restore करता है।

#### 2.2.3. अतिरिक्त नमूनों पर संरेखण चलाएं

दो नए ट्रिम किए गए नमूनों पर HISAT2 संरेखण चलाएं, एक के बाद एक।

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ2_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ2_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ2_1_trimmed.bam
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 18736
    		Aligned 0 time: 1531 (8.17%)
    		Aligned 1 time: 16726 (89.27%)
    		Aligned >1 times: 479 (2.56%)
    	Overall alignment rate: 91.83%
    ```

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COR1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COR1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COR1_1_trimmed.bam
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 38056
    		Aligned 0 time: 2311 (6.07%)
    		Aligned 1 time: 33289 (87.47%)
    		Aligned >1 times: 2456 (6.45%)
    	Overall alignment rate: 93.93%
    ```

एक बार यह पूरा हो जाने पर, तुम्हारे पास working डायरेक्टरी में दोनों नमूनों के लिए BAM और log फ़ाइलें होनी चाहिए।

#### 2.2.4. आउटपुट फ़ाइलों को ले जाएं

संरेखण आउटपुट फ़ाइलों को उसी डायरेक्टरी में ले जाएं जिसका हमने सेक्शन 1 में उपयोग किया था।

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/aligned
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    /data/aligned
    ├── ENCSR000COQ1_1_trimmed.bam
    ├── ENCSR000COQ1_1_trimmed.hisat2.log
    ├── ENCSR000COQ2_1_trimmed.bam
    ├── ENCSR000COQ2_1_trimmed.hisat2.log
    ├── ENCSR000COR1_1_trimmed.bam
    └── ENCSR000COR1_1_trimmed.hisat2.log
    ```

फ़ाइलें अब तुम्हारे सामान्य फाइलसिस्टम में accessible हैं।

#### 2.2.5. कंटेनर से बाहर निकलें

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करें।

```bash
exit
```

तुम्हारा prompt सामान्य हो जाना चाहिए।

### 2.3. एक व्यापक QC रिपोर्ट तैयार करें

अब जब हमारे पास तीन नमूनों के लिए QC, ट्रिमिंग और संरेखण आउटपुट है, हम उन्हें एक single रिपोर्ट में एकत्र करने के लिए MultiQC का उपयोग कर सकते हैं।
MultiQC compatible QC रिपोर्ट के लिए डायरेक्टरी में खोज करता है और जो कुछ भी मिलता है उसे एकत्र करता है।

#### 2.3.1. कंटेनर पुल करें

आइए एक container image पुल करें जिसमें `multiqc` इंस्टॉल है:

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "कमांड आउटपुट"

    ```console
    a3c26f6199d64b7c: Pulling from library/pip_multiqc
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
    2ed162b168e8: Pull complete
    ca06fe148f21: Pull complete
    Digest: sha256:af0e9de56896805aa2a065f7650362956f4213d99e95314f6fec472c6a3bf091
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

तुम देखोगे कि कुछ layers `Already exists` दिखाते हैं क्योंकि वे पहले हमने जो container images पुल की थीं उनके साथ साझा किए गए हैं।
परिणामस्वरूप, यह pull पिछले वाले की तुलना में तेज़ होना चाहिए।

#### 2.3.2. कंटेनर को इंटरैक्टिव रूप से स्पिन करें

डेटा डायरेक्टरी माउंट के साथ कंटेनर को इंटरैक्टिव रूप से स्पिन करें, पहले की तरह।

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

तुम्हारा prompt बदल जाएगा यह इंगित करने के लिए कि तुम कंटेनर के अंदर हो।

#### 2.3.3. MultiQC कमांड चलाएं

`multiqc` चलाएं, इसे उन डायरेक्टरी की ओर इंगित करते हुए जहाँ हमने सभी तीन नमूनों के लिए QC-संबंधित आउटपुट फ़ाइलें स्टोर कीं।
`-n` flag आउटपुट रिपोर्ट का नाम सेट करता है।

```bash
multiqc /data/reads /data/trimmed /data/aligned -n all_samples_QC
```

??? success "कमांड आउटपुट"

    ```console hl_lines="8 9 10 11 12"

    /// MultiQC 🔍 v1.32

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 36/36
               hisat2 | Found 3 reports
             cutadapt | Found 3 reports
               fastqc | Found 3 reports
        write_results | Data        : all_samples_QC_data
        write_results | Report      : all_samples_QC.html
              multiqc | MultiQC complete
    ```

यहाँ हम देखते हैं कि टूल ने सभी तीन नमूनों के लिए QC रिपोर्ट पाईं: `fastqc` से प्रारंभिक QC, `cutadapt` से पोस्ट-ट्रिमिंग रिपोर्ट (`trim_galore` के माध्यम से) और `hisat2` से संरेखण summaries।

आउटपुट फ़ाइलें working डायरेक्टरी में हैं:

```bash
ls all_samples_QC*
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    all_samples_QC.html

    all_samples_QC_data:
    cutadapt_filtered_reads_plot.txt                     multiqc.log
    cutadapt_trimmed_sequences_plot_3_Counts.txt         multiqc.parquet
    cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc_citations.txt
    fastqc-status-check-heatmap.txt                      multiqc_cutadapt.txt
    fastqc_adapter_content_plot.txt                      multiqc_data.json
    fastqc_overrepresented_sequences_plot.txt            multiqc_fastqc.txt
    fastqc_per_base_n_content_plot.txt                   multiqc_general_stats.txt
    fastqc_per_base_sequence_quality_plot.txt            multiqc_hisat2.txt
    fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_software_versions.txt
    fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_sources.txt
    fastqc_per_sequence_quality_scores_plot.txt
    fastqc_sequence_counts_plot.txt
    fastqc_sequence_duplication_levels_plot.txt
    fastqc_top_overrepresented_sequences_table.txt
    hisat2_se_plot.txt
    llms-full.txt
    ```

मुख्य आउटपुट `all_samples_QC.html` रिपोर्ट है, जो underlying मैट्रिक्स वाली एक डेटा डायरेक्टरी के साथ है।

#### 2.3.4. आउटपुट फ़ाइलों को ले जाएं

रिपोर्ट और इसकी डेटा डायरेक्टरी को माउंट किए गए फाइलसिस्टम में ले जाएं।

```bash
mkdir /data/multiqc
mv all_samples_QC* /data/multiqc
```

फ़ाइलें अब तुम्हारे सामान्य फाइलसिस्टम में accessible हैं।

#### 2.3.5. कंटेनर से बाहर निकलें

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करें।

```bash
exit
```

तुम्हारा prompt सामान्य हो जाना चाहिए।
यह सभी RNAseq प्रोसेसिंग commands के परीक्षण को समाप्त करता है।

---

### निष्कर्ष

तुम जानते हो कि FastQC, Trim Galore, HISAT2 और MultiQC commands को उनके संबंधित containers में कैसे चलाया जाए, जिसमें कई नमूनों को प्रोसेस करना और QC रिपोर्ट को एकत्र करना शामिल है।

### आगे क्या है?

एक ब्रेक लें, फिर [भाग 2](./02_single-sample.md) पर जाएं यह सीखने के लिए कि उन्हीं commands को workflows में कैसे लपेटा जाए जो काम को निष्पादित करने के लिए containers का उपयोग करते हैं।
