# भाग 1: विधि का अवलोकन

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

वेरिएंट कॉलिंग एक जीनोमिक विश्लेषण विधि है जिसका उद्देश्य संदर्भ जीनोम की तुलना में जीनोम अनुक्रम में भिन्नता की पहचान करना है।
यहां हम होल-जीनोम सीक्वेंसिंग डेटा में शॉर्ट जर्मलाइन वेरिएंट, _यानी_ SNPs और indels, को कॉल करने के लिए डिज़ाइन किए गए टूल्स और विधियों का उपयोग करने जा रहे हैं।

![GATK pipeline](img/gatk-pipeline.png)

एक पूर्ण वेरिएंट कॉलिंग पाइपलाइन में आमतौर पर कई चरण शामिल होते हैं, जिनमें रेफरेंस से मैपिंग (कभी-कभी जीनोम संरेखण के रूप में संदर्भित) और वेरिएंट फ़िल्टरिंग और प्राथमिकता शामिल हैं।
सरलता के लिए, इस कोर्स में हम सिर्फ वेरिएंट कॉलिंग भाग पर ध्यान केंद्रित करने जा रहे हैं।

### विधियां

हम तुम्हें होल-जीनोम सीक्वेंसिंग नमूनों पर वेरिएंट कॉलिंग लागू करने के दो तरीके दिखाने जा रहे हैं ताकि जर्मलाइन SNPs और indels की पहचान की जा सके।
पहले हम एक सरल **प्रति-नमूना दृष्टिकोण** से शुरुआत करेंगे जो प्रत्येक नमूने से स्वतंत्र रूप से वेरिएंट को कॉल करता है।
फिर हम तुम्हें एक अधिक परिष्कृत **संयुक्त कॉलिंग दृष्टिकोण** दिखाएंगे जो कई नमूनों का एक साथ विश्लेषण करता है, और अधिक सटीक और जानकारीपूर्ण परिणाम उत्पन्न करता है।

इससे पहले कि हम किसी भी दृष्टिकोण के लिए कोई वर्कफ़्लो कोड लिखना शुरू करें, हम कुछ परीक्षण डेटा पर मैनुअल रूप से कमांड्स को आज़माने जा रहे हैं।

### डेटासेट

हम निम्नलिखित डेटा और संबंधित संसाधन प्रदान करते हैं:

- **एक रेफरेंस जीनोम** जो मानव क्रोमोसोम 20 (hg19/b37 से) के एक छोटे से क्षेत्र और इसकी सहायक फ़ाइलों (इंडेक्स और अनुक्रम शब्दकोश) से बना है।
- **तीन होल जीनोम सीक्वेंसिंग नमूने** जो एक परिवार ट्रायो (मां, पिता और बेटा) से संबंधित हैं, जिन्हें क्रोमोसोम 20 पर डेटा के एक छोटे से हिस्से तक सीमित कर दिया गया है ताकि फ़ाइल का आकार छोटा रहे।
  यह Illumina शॉर्ट-रीड सीक्वेंसिंग डेटा है जिसे पहले से ही रेफरेंस जीनोम से मैप किया जा चुका है, [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) फॉर्मेट में प्रदान किया गया है (Binary Alignment Map, SAM का संकुचित संस्करण, Sequence Alignment Map)।
- **जीनोमिक अंतराल की एक सूची**, यानी जीनोम पर निर्देशांक जहां हमारे नमूनों में वेरिएंट कॉल करने के लिए उपयुक्त डेटा है, BED फॉर्मेट में प्रदान किया गया है।

### सॉफ्टवेयर

शामिल दो मुख्य टूल हैं [Samtools](https://www.htslib.org/), अनुक्रम संरेखण फ़ाइलों में हेरफेर करने के लिए व्यापक रूप से उपयोग किया जाने वाला टूलकिट, और [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), Broad Institute में विकसित वेरिएंट खोज के लिए टूल्स का एक सेट।

ये टूल्स GitHub Codespaces वातावरण में इंस्टॉल नहीं हैं, इसलिए हम उन्हें Seqera Containers सेवा के माध्यम से प्राप्त कंटेनरों के माध्यम से उपयोग करेंगे ([Hello Containers](../../hello_nextflow/05_hello_containers.md) देखें)।

!!! tip "सुझाव"

    सुनिश्चित करो कि तुम `nf4-science/genomics` डायरेक्टरी में हो ताकि जब तुम `pwd` टाइप करो तो पथ का अंतिम भाग `genomics` दिखे।

---

## 1. प्रति-नमूना वेरिएंट कॉलिंग

प्रति-नमूना वेरिएंट कॉलिंग प्रत्येक नमूने को स्वतंत्र रूप से प्रोसेस करती है: वेरिएंट कॉलर एक समय में एक नमूने के लिए सीक्वेंसिंग डेटा की जांच करता है और उन पोजीशन्स की पहचान करता है जहां नमूना रेफरेंस से अलग है।

इस सेक्शन में हम प्रति-नमूना वेरिएंट कॉलिंग दृष्टिकोण के दो कमांड्स का परीक्षण करते हैं: Samtools के साथ BAM फ़ाइल को इंडेक्सिंग करना और GATK HaplotypeCaller के साथ वेरिएंट कॉल करना।
ये वे कमांड्स हैं जिन्हें हम इस कोर्स के भाग 2 में Nextflow वर्कफ़्लो में wrap करेंगे।

1. [Samtools](https://www.htslib.org/) का उपयोग करके BAM इनपुट फ़ाइल के लिए एक इंडेक्स फ़ाइल जेनरेट करें
2. इंडेक्स की गई BAM फ़ाइल पर GATK HaplotypeCaller चलाएं ताकि VCF (Variant Call Format) में प्रति-नमूना वेरिएंट कॉल जेनरेट हो सकें

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

हम सिर्फ एक नमूने पर दो कमांड्स का परीक्षण करके शुरुआत करते हैं।

### 1.1. Samtools के साथ BAM इनपुट फ़ाइल को इंडेक्स करें

इंडेक्स फ़ाइलें बायोइन्फॉर्मेटिक्स फ़ाइल फॉर्मेट की एक सामान्य विशेषता हैं; इनमें मुख्य फ़ाइल की संरचना के बारे में जानकारी होती है जो GATK जैसे टूल्स को पूरी फ़ाइल पढ़े बिना डेटा के एक सबसेट तक पहुंचने की अनुमति देती है।
यह महत्वपूर्ण है क्योंकि ये फ़ाइलें कितनी बड़ी हो सकती हैं।

BAM फ़ाइलें अक्सर इंडेक्स के बिना प्रदान की जाती हैं, इसलिए कई विश्लेषण वर्कफ़्लो में पहला चरण `samtools index` का उपयोग करके एक जेनरेट करना है।

हम एक Samtools कंटेनर को पुल करने जा रहे हैं, इसे इंटरएक्टिव रूप से स्पिन अप करेंगे और BAM फ़ाइलों में से एक पर `samtools index` कमांड चलाएंगे।

#### 1.1.1. Samtools कंटेनर को पुल करें

Samtools कंटेनर इमेज डाउनलोड करने के लिए `docker pull` कमांड चलाएं:

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

??? success "कमांड आउटपुट"

    ```console
    1.20--b5dfbd93de237464: Pulling from library/samtools
    6360b3717211: Pull complete
    2ec3f7ad9b3c: Pull complete
    7716ca300600: Pull complete
    4f4fb700ef54: Pull complete
    8c61d418774c: Pull complete
    03dae77ff45c: Pull complete
    aab7f787139d: Pull complete
    4f4fb700ef54: Pull complete
    837d55536720: Pull complete
    897362c12ca7: Pull complete
    3893cbe24e91: Pull complete
    d1b61e94977b: Pull complete
    c72ff66fb90f: Pull complete
    0e0388f29b6d: Pull complete
    Digest: sha256:bbfc45b4f228975bde86cba95e303dd94ecf2fdacea5bfb2e2f34b0d7b141e41
    Status: Downloaded newer image for community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    ```

अगर तुमने यह इमेज पहले डाउनलोड नहीं की है, तो इसे पूरा होने में एक मिनट लग सकता है।
एक बार यह हो जाने के बाद, तुम्हारे पास कंटेनर इमेज की एक लोकल कॉपी है।

#### 1.1.2. Samtools कंटेनर को इंटरएक्टिव रूप से स्पिन अप करें

कंटेनर को इंटरएक्टिव रूप से चलाने के लिए, `-it` फ्लैग के साथ `docker run` का उपयोग करो।
`-v ./data:/data` विकल्प लोकल `data` डायरेक्टरी को कंटेनर में माउंट करता है ताकि टूल्स इनपुट फ़ाइलों तक पहुंच सकें।

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

??? success "कमांड आउटपुट"

    ```console
    (base) root@1409896f77b1:/tmp#
    ```

तुम्हारा प्रॉम्प्ट `(base) root@a1b2c3d4e5f6:/tmp#` जैसे कुछ में बदल जाता है, जो दर्शाता है कि तुम अब कंटेनर के अंदर हो।

सत्यापित करो कि तुम `/data/bam` के तहत सीक्वेंस डेटा फ़ाइलें देख सकते हो:

```bash
ls /data/bam
```

??? success "कमांड आउटपुट"

    ```console
    reads_father.bam  reads_mother.bam  reads_mother.bam.bai  reads_son.bam
    ```

इसके साथ, तुम अपनी पहली कमांड आज़माने के लिए तैयार हो।

#### 1.1.3. इंडेक्सिंग कमांड चलाएं

[Samtools डॉक्यूमेंटेशन](https://www.htslib.org/doc/samtools-index.html) हमें BAM फ़ाइल को इंडेक्स करने के लिए चलाने वाली कमांड लाइन देता है।
हमें केवल इनपुट फ़ाइल प्रदान करने की आवश्यकता है; टूल स्वचालित रूप से इनपुट फ़ाइलनाम में `.bai` जोड़कर आउटपुट के लिए एक नाम जेनरेट करेगा।

एक डेटा फ़ाइल पर `samtools index` कमांड चलाएं:

```bash
samtools index /data/bam/reads_mother.bam
```

कमांड टर्मिनल में कोई आउटपुट उत्पन्न नहीं करती है, लेकिन तुम्हें अब मूल BAM इनपुट फ़ाइल के समान डायरेक्टरी में `reads_mother.bam.bai` नामक एक फ़ाइल दिखनी चाहिए।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

यह पहले चरण के परीक्षण को पूरा करता है।

#### 1.1.4. Samtools कंटेनर से बाहर निकलें

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करो।

```bash
exit
```

तुम्हारा प्रॉम्प्ट अब वापस वैसा होना चाहिए जैसा कंटेनर शुरू करने से पहले था।

### 1.2. GATK HaplotypeCaller के साथ वेरिएंट कॉल करें

हम BAM फ़ाइल पर `gatk HaplotypeCaller` कमांड चलाना चाहते हैं जिसे हमने अभी इंडेक्स किया है।

#### 1.2.1. GATK कंटेनर को पुल करें

पहले, आओ GATK कंटेनर इमेज डाउनलोड करने के लिए `docker pull` कमांड चलाएं:

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

??? success "कमांड आउटपुट"

    कुछ लेयर्स `Already exists` दिखाती हैं क्योंकि वे Samtools कंटेनर इमेज के साथ शेयर की गई हैं जिसे हमने पहले पुल किया था।

    ```console
    4.5.0.0--730ee8817e436867: Pulling from library/gatk4
    6360b3717211: Already exists
    2ec3f7ad9b3c: Already exists
    7716ca300600: Already exists
    4f4fb700ef54: Already exists
    8c61d418774c: Already exists
    03dae77ff45c: Already exists
    aab7f787139d: Already exists
    4f4fb700ef54: Already exists
    837d55536720: Already exists
    897362c12ca7: Already exists
    3893cbe24e91: Already exists
    d1b61e94977b: Already exists
    e5c558f54708: Pull complete
    087cce32d294: Pull complete
    Digest: sha256:e33413b9100f834fcc62fd5bc9edc1e881e820aafa606e09301eac2303d8724b
    Status: Downloaded newer image for community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    ```

यह पहले पुल की तुलना में तेज़ होना चाहिए क्योंकि दोनों कंटेनर इमेजेस अपनी अधिकांश लेयर्स शेयर करती हैं।

#### 1.2.2. GATK कंटेनर को इंटरएक्टिव रूप से स्पिन अप करें

GATK कंटेनर को डेटा डायरेक्टरी माउंट के साथ इंटरएक्टिव रूप से स्पिन अप करो, जैसे हमने Samtools के लिए किया था।

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

तुम्हारा प्रॉम्प्ट बदल जाता है जो दर्शाता है कि तुम अब GATK कंटेनर के अंदर हो।

#### 1.2.3. वेरिएंट कॉलिंग कमांड चलाएं

[GATK डॉक्यूमेंटेशन](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) हमें BAM फ़ाइल पर वेरिएंट कॉलिंग करने के लिए चलाने वाली कमांड लाइन देता है।

हमें BAM इनपुट फ़ाइल (`-I`) के साथ-साथ रेफरेंस जीनोम (`-R`), आउटपुट फ़ाइल के लिए एक नाम (`-O`) और विश्लेषण करने के लिए जीनोमिक अंतराल की एक सूची (`-L`) प्रदान करने की आवश्यकता है।

हालांकि, हमें इंडेक्स फ़ाइल का पथ निर्दिष्ट करने की आवश्यकता नहीं है; टूल स्वचालित रूप से स्थापित नामकरण और सह-स्थान सम्मेलन के आधार पर इसे उसी डायरेक्टरी में देखेगा।
यही रेफरेंस जीनोम की सहायक फ़ाइलों (इंडेक्स और अनुक्रम शब्दकोश फ़ाइलें, `*.fai` और `*.dict`) पर भी लागू होता है।

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O /data/vcf/reads_mother.vcf \
        -L /data/ref/intervals.bed
```

??? success "कमांड आउटपुट"

    ```console hl_lines="37 51 56 57"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.vcf -L /data/ref/intervals.bed
    00:27:50.687 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:27:50.854 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.858 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:27:50.858 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:27:50.858 INFO  HaplotypeCaller - Executing as root@a1fe8ff42d07 on Linux v6.10.14-linuxkit amd64
    00:27:50.858 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:27:50.859 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:27:50 AM GMT
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.861 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:27:50.861 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:27:50.861 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:27:50.863 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:27:50.864 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:27:50.864 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:27:50.864 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:27:50.864 INFO  HaplotypeCaller - Requester pays: disabled
    00:27:50.865 INFO  HaplotypeCaller - Initializing engine
    00:27:50.991 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:27:51.016 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:27:51.029 INFO  HaplotypeCaller - Done initializing engine
    00:27:51.040 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:27:51.042 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:27:51.042 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:27:51.046 INFO  HaplotypeCallerEngine - Disabling physical phasing, which is supported only for reference-model confidence output
    00:27:51.063 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:27:51.085 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:27:51.086 INFO  IntelPairHmm - Available threads: 10
    00:27:51.086 INFO  IntelPairHmm - Requested threads: 4
    00:27:51.086 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:27:51.128 INFO  ProgressMeter - Starting traversal
    00:27:51.136 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:27:51.882 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:27:52.969 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:27:52.971 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35           1145.7
    00:27:52.971 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:27:52.976 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003346916
    00:27:52.976 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.045731709
    00:27:52.977 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:27:52.981 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:27:52 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.04 minutes.
    Runtime.totalMemory()=203423744
    ```

लॉग आउटपुट बहुत verbose है, इसलिए हमने ऊपर दिए गए उदाहरण में सबसे प्रासंगिक लाइनों को हाइलाइट किया है।

आउटपुट फ़ाइलें, `reads_mother.vcf` और इसकी इंडेक्स फ़ाइल, `reads_mother.vcf.idx`, कंटेनर में तुम्हारी कार्यशील डायरेक्टरी के अंदर बनाई गई हैं।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    conda.yml  hsperfdata_root  reads_mother.vcf  reads_mother.vcf.idx
    ```

VCF फ़ाइल में वेरिएंट कॉल होते हैं, जैसा कि हम एक मिनट में देखेंगे, और इंडेक्स फ़ाइल का वही कार्य है जो BAM इंडेक्स फ़ाइल का है, टूल्स को संपूर्ण फ़ाइल लोड किए बिना डेटा के सबसेट को खोजने और पुनः प्राप्त करने की अनुमति देना।

चूंकि VCF एक टेक्स्ट फॉर्मेट है और यह एक छोटी परीक्षण फ़ाइल है, तुम इसे खोलने और सामग्री देखने के लिए `cat reads_mother.vcf` चला सकते हो।
अगर तुम फ़ाइल की शुरुआत तक स्क्रॉल करते हो, तो तुम्हें मेटाडेटा की कई लाइनों से बना एक हेडर मिलेगा, उसके बाद वेरिएंट कॉल की एक सूची होगी, प्रति लाइन एक।

??? abstract "फ़ाइल की सामग्री (संक्षिप्त)"

    ```console title="reads_mother.vcf" linenums="1" hl_lines="26"
    ##fileformat=VCFv4.2
    ##FILTER=<ID=LowQual,Description="Low quality">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
    ##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller --output reads_mother.vcf --intervals /data/ref/intervals.bed --input /data/bam/reads_mother.bam --reference /data/ref/ref.fasta [abridged]",Version="4.5.0.0",Date="February 11, 2026 at 4:23:43 PM GMT">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    ##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
    ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
    ##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
    ##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
    ##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
    ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
    ##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
    ##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
    ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
    ##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
    ##contig=<ID=20_10037292_10066351,length=29059>
    ##source=HaplotypeCaller
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_mother
    20_10037292_10066351    3480    .       C       CT      503.03  .       AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179     GT:AD:DP:GQ:PL  1/1:0,18:18:54:517,54,0
    20_10037292_10066351    3520    .       AT      A       609.03  .       AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693     GT:AD:DP:GQ:PL  1/1:0,18:18:54:623,54,0
    20_10037292_10066351    3529    .       T       A       155.64  .       AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034       GT:AD:DP:GQ:PL  0/1:12,8:20:99:163,0,328
    20_10037292_10066351    4012    .       C       T       1398.06 .       AC=2;AF=1.00;AN=2;DP=44;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=32.51;SOR=0.739     GT:AD:DP:GQ:PL  1/1:0,43:43:99:1412,129,0
    20_10037292_10066351    4409    .       A       ATATG   710.03  .       AC=2;AF=1.00;AN=2;DP=31;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=30.87;SOR=0.784     GT:AD:DP:GQ:PL  1/1:0,23:23:69:724,69,0
    20_10037292_10066351    5027    .       C       T       784.06  .       AC=2;AF=1.00;AN=2;DP=27;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=30.16;SOR=0.693     GT:AD:DP:GQ:PL  1/1:0,26:26:77:798,77,0
    20_10037292_10066351    5469    .       A       G       1297.06 .       AC=2;AF=1.00;AN=2;DP=42;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=30.88;SOR=1.005     GT:AD:DP:GQ:PL  1/1:0,42:42:99:1311,126,0
    20_10037292_10066351    7557    .       A       G       935.06  .       AC=2;AF=1.00;AN=2;DP=36;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.50;SOR=0.693     GT:AD:DP:GQ:PL  1/1:0,34:34:99:949,100,0
    20_10037292_10066351    7786    .       G       T       1043.06 .       AC=2;AF=1.00;AN=2;DP=35;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=30.68;SOR=0.941     GT:AD:DP:GQ:PL  1/1:0,34:34:99:1057,102,0
    20_10037292_10066351    8350    .       G       C       1162.06 .       AC=2;AF=1.00;AN=2;DP=39;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=29.80;SOR=1.096     GT:AD:DP:GQ:PL  1/1:0,39:39:99:1176,115,0
    20_10037292_10066351    8886    .       AAGAAAGAAAG     A       1268.03 .       AC=2;AF=1.00;AN=2;DP=34;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=25.36;SOR=1.071     GT:AD:DP:GQ:PL  1/1:0,29:29:88:1282,88,0
    20_10037292_10066351    13536   .       T       C       437.64  .       AC=1;AF=0.500;AN=2;BaseQRankSum=1.454;DP=45;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=9.95;ReadPosRankSum=-1.613;SOR=0.818        GT:AD:DP:GQ:PL  0/1:26,18:44:99:445,0,672
    20_10037292_10066351    14156   .       T       C       183.64  .       AC=1;AF=0.500;AN=2;BaseQRankSum=0.703;DP=20;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=9.18;ReadPosRankSum=-0.193;SOR=1.034        GT:AD:DP:GQ:PL  0/1:12,8:20:99:191,0,319
    ```

ऊपर दिए गए उदाहरण आउटपुट में, हमने अंतिम हेडर लाइन को हाइलाइट किया है, जो उसके बाद आने वाले टेबुलर डेटा के लिए कॉलम के नाम देती है।
डेटा की प्रत्येक लाइन नमूने के सीक्वेंसिंग डेटा में पहचाने गए संभावित वेरिएंट का वर्णन करती है। VCF फॉर्मेट की व्याख्या करने के लिए मार्गदर्शन के लिए, [यह उपयोगी लेख](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/) देखें।

#### 1.2.4. आउटपुट फ़ाइलों को स्थानांतरित करें

कंटेनर के अंदर जो कुछ भी रहता है वह भविष्य के काम के लिए दुर्गम होगा।
BAM इंडेक्स फ़ाइल सीधे माउंट की गई फ़ाइलसिस्टम पर `/data/bam` डायरेक्टरी में बनाई गई थी, लेकिन VCF फ़ाइल और इसकी इंडेक्स नहीं, इसलिए हमें उन दोनों को मैनुअल रूप से स्थानांतरित करने की आवश्यकता है।

```bash
mkdir /data/vcf
mv reads_mother.vcf* /data/vcf
```

??? abstract "डायरेक्टरी सामग्री"

    ```console hl_lines="5 13-15"
    data/
    ├── bam
    │   ├── reads_father.bam
    │   ├── reads_mother.bam
    │   ├── reads_mother.bam.bai
    │   └── reads_son.bam
    ├── ref
    │   ├── intervals.bed
    │   ├── ref.dict
    │   ├── ref.fasta
    │   └── ref.fasta.fai
    ├── samplesheet.csv
    └── vcf
        ├── reads_mother.vcf
        └── reads_mother.vcf.idx
    ```

एक बार यह हो जाने के बाद, सभी फ़ाइलें अब तुम्हारी सामान्य फ़ाइलसिस्टम में सुलभ हैं।

#### 1.2.5. GATK कंटेनर से बाहर निकलें

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करो।

```bash
exit
```

तुम्हारा प्रॉम्प्ट वापस सामान्य होना चाहिए।
यह प्रति-नमूना वेरिएंट कॉलिंग परीक्षण समाप्त करता है।

!!! example "इसे वर्कफ़्लो के रूप में लिखो!"

    अगर तुम इस विश्लेषण को Nextflow वर्कफ़्लो के रूप में लागू करना शुरू करना चाहते हो तो तुम तुरंत [भाग 2](./02_per_sample_variant_calling.md) पर जा सकते हो।
    तुम्हें बस भाग 3 पर जाने से पहले परीक्षण के दूसरे दौर को पूरा करने के लिए वापस आना होगा।

---

## 2. कोहॉर्ट पर संयुक्त कॉलिंग

जिस वेरिएंट कॉलिंग दृष्टिकोण का हमने अभी उपयोग किया, वह प्रति नमूना वेरिएंट कॉल जेनरेट करता है।
यह प्रत्येक नमूने से अलग-अलग वेरिएंट को देखने के लिए ठीक है, लेकिन यह सीमित जानकारी देता है।
यह अक्सर अधिक दिलचस्प होता है कि वेरिएंट कॉल कई नमूनों में कैसे भिन्न होते हैं।
GATK इस उद्देश्य के लिए संयुक्त वेरिएंट कॉलिंग नामक एक वैकल्पिक विधि प्रदान करता है।

संयुक्त वेरिएंट कॉलिंग में प्रत्येक नमूने के लिए GVCF (Genomic VCF के लिए) नामक एक विशेष प्रकार का वेरिएंट आउटपुट जेनरेट करना शामिल है, फिर सभी नमूनों से GVCF डेटा को जोड़ना और एक 'संयुक्त जीनोटाइपिंग' सांख्यिकीय विश्लेषण चलाना।

![Joint analysis](img/joint-calling.png)

एक नमूने के GVCF में जो विशेष है वह यह है कि इसमें जीनोम के लक्षित क्षेत्र में सभी पोजीशन्स के बारे में अनुक्रम डेटा सांख्यिकी को सारांशित करने वाले रिकॉर्ड होते हैं, न कि केवल उन पोजीशन्स के जहां प्रोग्राम को भिन्नता के प्रमाण मिले।
यह संयुक्त जीनोटाइपिंग गणना के लिए महत्वपूर्ण है ([आगे पढ़ना](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants))।

GVCF को GATK HaplotypeCaller द्वारा उत्पादित किया जाता है, वही टूल जिसे हमने अभी परीक्षण किया, एक अतिरिक्त पैरामीटर (`-ERC GVCF`) के साथ।
GVCFs को जोड़ना GATK GenomicsDBImport के साथ किया जाता है, जो प्रति-नमूना कॉल को GenomicsDB डेटा स्टोर में जोड़ता है।
GenomicsDB डेटा स्टोर एक प्रकार का डेटाबेस फॉर्मेट है जो वेरिएंट जानकारी के लिए एक मध्यवर्ती भंडारण के रूप में कार्य करता है।

वास्तविक 'संयुक्त जीनोटाइपिंग' विश्लेषण फिर GATK GenotypeGVCFs के साथ किया जाता है।

यहां हम GVCFs जेनरेट करने और संयुक्त जीनोटाइपिंग चलाने के लिए आवश्यक कमांड्स का परीक्षण करते हैं।
ये वे कमांड्स हैं जिन्हें हम इस कोर्स के भाग 3 में Nextflow वर्कफ़्लो में wrap करेंगे।

1. Samtools का उपयोग करके प्रत्येक BAM इनपुट फ़ाइल के लिए एक इंडेक्स फ़ाइल जेनरेट करें
2. प्रति-नमूना जीनोमिक वेरिएंट कॉल का GVCF जेनरेट करने के लिए प्रत्येक BAM इनपुट फ़ाइल पर GATK HaplotypeCaller चलाएं
3. सभी GVCFs को एकत्र करें और उन्हें GenomicsDB डेटा स्टोर में जोड़ें
4. कोहॉर्ट-स्तर VCF उत्पन्न करने के लिए संयुक्त GVCF डेटा स्टोर पर संयुक्त जीनोटाइपिंग चलाएं

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

हमें अब इन सभी कमांड्स का परीक्षण करने की आवश्यकता है, सभी तीन BAM फ़ाइलों को इंडेक्सिंग करने से शुरुआत करते हुए।

### 2.1. सभी तीन नमूनों के लिए BAM फ़ाइलों को इंडेक्स करें

ऊपर पहले सेक्शन में, हमने केवल एक BAM फ़ाइल को इंडेक्स किया था।
अब हमें सभी तीन नमूनों को इंडेक्स करने की आवश्यकता है ताकि GATK HaplotypeCaller उन्हें प्रोसेस कर सके।

#### 2.1.1. Samtools कंटेनर को इंटरएक्टिव रूप से स्पिन अप करें

हमने पहले ही Samtools कंटेनर इमेज को पुल कर लिया है, इसलिए हम इसे सीधे स्पिन अप कर सकते हैं:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

तुम्हारा प्रॉम्प्ट बदलता है जो दर्शाता है कि तुम कंटेनर के अंदर हो, डेटा डायरेक्टरी पहले की तरह माउंट है।

#### 2.1.2. सभी तीन नमूनों पर इंडेक्सिंग कमांड चलाएं

तीनों BAM फ़ाइलों में से प्रत्येक पर इंडेक्सिंग कमांड चलाएं:

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

यह संबंधित BAM फ़ाइलों के समान डायरेक्टरी में इंडेक्स फ़ाइलें उत्पन्न करना चाहिए।

#### 2.1.3. Samtools कंटेनर से बाहर निकलें

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करो।

```bash
exit
```

तुम्हारा प्रॉम्प्ट वापस सामान्य होना चाहिए।

### 2.2. सभी तीन नमूनों के लिए GVCFs जेनरेट करें

संयुक्त जीनोटाइपिंग चरण चलाने के लिए, हमें सभी तीन नमूनों के लिए GVCFs की आवश्यकता है।

#### 2.2.1. GATK कंटेनर को इंटरएक्टिव रूप से स्पिन अप करें

हमने पहले ही GATK कंटेनर इमेज को पुल कर लिया है, इसलिए हम इसे सीधे स्पिन अप कर सकते हैं:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

तुम्हारा प्रॉम्प्ट बदलता है जो दर्शाता है कि तुम GATK कंटेनर के अंदर हो।

#### 2.2.2. GVCF विकल्प के साथ वेरिएंट कॉलिंग कमांड चलाएं

जीनोमिक VCF (GVCF) उत्पन्न करने के लिए, हम बेस कमांड में `-ERC GVCF` विकल्प जोड़ते हैं, जो HaplotypeCaller के GVCF मोड को चालू करता है।

हम आउटपुट फ़ाइल के लिए फ़ाइल एक्सटेंशन को `.vcf` से `.g.vcf` में भी बदल देते हैं।
यह तकनीकी रूप से आवश्यकता नहीं है, लेकिन यह एक दृढ़ता से अनुशंसित सम्मेलन है।

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

??? success "कमांड आउटपुट"

    ```console hl_lines="39 53 58 59"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.g.vcf -L /data/ref/intervals.bed -ERC GVCF
    16:51:00.620 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    16:51:00.749 INFO  HaplotypeCaller - ------------------------------------------------------------
    16:51:00.751 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    16:51:00.751 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    16:51:00.751 INFO  HaplotypeCaller - Executing as root@be1a0302f6c7 on Linux v6.8.0-1030-azure amd64
    16:51:00.751 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    16:51:00.752 INFO  HaplotypeCaller - Start Date/Time: February 11, 2026 at 4:51:00 PM GMT
    16:51:00.752 INFO  HaplotypeCaller - ------------------------------------------------------------
    16:51:00.752 INFO  HaplotypeCaller - ------------------------------------------------------------
    16:51:00.752 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    16:51:00.753 INFO  HaplotypeCaller - Picard Version: 3.1.1
    16:51:00.753 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    16:51:00.753 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    16:51:00.753 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    16:51:00.753 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    16:51:00.754 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    16:51:00.754 INFO  HaplotypeCaller - Deflater: IntelDeflater
    16:51:00.754 INFO  HaplotypeCaller - Inflater: IntelInflater
    16:51:00.754 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    16:51:00.754 INFO  HaplotypeCaller - Requester pays: disabled
    16:51:00.755 INFO  HaplotypeCaller - Initializing engine
    16:51:00.893 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    16:51:00.905 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    16:51:00.910 INFO  HaplotypeCaller - Done initializing engine
    16:51:00.912 INFO  HaplotypeCallerEngine - Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled
    16:51:00.917 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    16:51:00.919 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    16:51:00.919 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    16:51:00.923 INFO  HaplotypeCallerEngine - Standard Emitting and Calling confidence set to -0.0 for reference-model confidence output
    16:51:00.923 INFO  HaplotypeCallerEngine - All sites annotated with PLs forced to true for reference-model confidence output
    16:51:00.933 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    16:51:00.945 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    16:51:00.945 INFO  IntelPairHmm - Available threads: 4
    16:51:00.945 INFO  IntelPairHmm - Requested threads: 4
    16:51:00.945 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    16:51:00.984 INFO  ProgressMeter - Starting traversal
    16:51:00.985 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    16:51:01.452 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    16:51:02.358 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    16:51:02.359 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35           1529.5
    16:51:02.359 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    16:51:02.361 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.0022800000000000003
    16:51:02.361 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.061637120000000004
    16:51:02.361 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    16:51:02.362 INFO  HaplotypeCaller - Shutting down engine
    [February 11, 2026 at 4:51:02 PM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.03 minutes.
    Runtime.totalMemory()=257949696
    ```

यह कंटेनर में वर्तमान कार्यशील डायरेक्टरी में GVCF आउटपुट फ़ाइल `reads_mother.g.vcf` बनाता है, साथ ही इसकी इंडेक्स फ़ाइल, `reads_mother.g.vcf.idx`।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    conda.yml  hsperfdata_root  reads_mother.g.vcf  reads_mother.g.vcf.idx
    ```

अगर तुम फ़ाइल की सामग्री की पहली 200 लाइनें देखने के लिए `head -200 reads_mother.g.vcf` चलाते हो, तो तुम देखोगे कि यह सेक्शन 1 में जेनरेट किए गए समकक्ष VCF की तुलना में बहुत लंबा है, और अधिकांश लाइनें VCF में देखी गई चीज़ों से काफी अलग दिखती हैं।

??? abstract "फ़ाइल की सामग्री (संक्षिप्त)"

    ```console title="reads_mother.g.vcf" linenums="1" hl_lines="92 175 191 195"
    ##fileformat=VCFv4.2
    ##ALT=<ID=NON_REF,Description="Represents any possible alternative allele not already represented at this location by REF and ALT">
    ##FILTER=<ID=LowQual,Description="Low quality">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
    ##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles">
    ##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
    ##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">
    ##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
    ##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller --emit-ref-confidence GVCF --output reads_mother.g.vcf --intervals /data/ref/intervals.bed --input /data/bam/reads_mother.bam --reference /data/ref/ref.fasta [abridged]",Version="4.5.0.0",Date="February 11, 2026 at 4:51:00 PM GMT">
    ##GVCFBlock0-1=minGQ=0(inclusive),maxGQ=1(exclusive)
    ##GVCFBlock1-2=minGQ=1(inclusive),maxGQ=2(exclusive)
    ##GVCFBlock10-11=minGQ=10(inclusive),maxGQ=11(exclusive)
    ##GVCFBlock11-12=minGQ=11(inclusive),maxGQ=12(exclusive)
    ##GVCFBlock12-13=minGQ=12(inclusive),maxGQ=13(exclusive)
    ##GVCFBlock13-14=minGQ=13(inclusive),maxGQ=14(exclusive)
    ##GVCFBlock14-15=minGQ=14(inclusive),maxGQ=15(exclusive)
    ##GVCFBlock15-16=minGQ=15(inclusive),maxGQ=16(exclusive)
    ##GVCFBlock16-17=minGQ=16(inclusive),maxGQ=17(exclusive)
    ##GVCFBlock17-18=minGQ=17(inclusive),maxGQ=18(exclusive)
    ##GVCFBlock18-19=minGQ=18(inclusive),maxGQ=19(exclusive)
    ##GVCFBlock19-20=minGQ=19(inclusive),maxGQ=20(exclusive)
    ##GVCFBlock2-3=minGQ=2(inclusive),maxGQ=3(exclusive)
    ##GVCFBlock20-21=minGQ=20(inclusive),maxGQ=21(exclusive)
    ##GVCFBlock21-22=minGQ=21(inclusive),maxGQ=22(exclusive)
    ##GVCFBlock22-23=minGQ=22(inclusive),maxGQ=23(exclusive)
    ##GVCFBlock23-24=minGQ=23(inclusive),maxGQ=24(exclusive)
    ##GVCFBlock24-25=minGQ=24(inclusive),maxGQ=25(exclusive)
    ##GVCFBlock25-26=minGQ=25(inclusive),maxGQ=26(exclusive)
    ##GVCFBlock26-27=minGQ=26(inclusive),maxGQ=27(exclusive)
    ##GVCFBlock27-28=minGQ=27(inclusive),maxGQ=28(exclusive)
    ##GVCFBlock28-29=minGQ=28(inclusive),maxGQ=29(exclusive)
    ##GVCFBlock29-30=minGQ=29(inclusive),maxGQ=30(exclusive)
    ##GVCFBlock3-4=minGQ=3(inclusive),maxGQ=4(exclusive)
    ##GVCFBlock30-31=minGQ=30(inclusive),maxGQ=31(exclusive)
    ##GVCFBlock31-32=minGQ=31(inclusive),maxGQ=32(exclusive)
    ##GVCFBlock32-33=minGQ=32(inclusive),maxGQ=33(exclusive)
    ##GVCFBlock33-34=minGQ=33(inclusive),maxGQ=34(exclusive)
    ##GVCFBlock34-35=minGQ=34(inclusive),maxGQ=35(exclusive)
    ##GVCFBlock35-36=minGQ=35(inclusive),maxGQ=36(exclusive)
    ##GVCFBlock36-37=minGQ=36(inclusive),maxGQ=37(exclusive)
    ##GVCFBlock37-38=minGQ=37(inclusive),maxGQ=38(exclusive)
    ##GVCFBlock38-39=minGQ=38(inclusive),maxGQ=39(exclusive)
    ##GVCFBlock39-40=minGQ=39(inclusive),maxGQ=40(exclusive)
    ##GVCFBlock4-5=minGQ=4(inclusive),maxGQ=5(exclusive)
    ##GVCFBlock40-41=minGQ=40(inclusive),maxGQ=41(exclusive)
    ##GVCFBlock41-42=minGQ=41(inclusive),maxGQ=42(exclusive)
    ##GVCFBlock42-43=minGQ=42(inclusive),maxGQ=43(exclusive)
    ##GVCFBlock43-44=minGQ=43(inclusive),maxGQ=44(exclusive)
    ##GVCFBlock44-45=minGQ=44(inclusive),maxGQ=45(exclusive)
    ##GVCFBlock45-46=minGQ=45(inclusive),maxGQ=46(exclusive)
    ##GVCFBlock46-47=minGQ=46(inclusive),maxGQ=47(exclusive)
    ##GVCFBlock47-48=minGQ=47(inclusive),maxGQ=48(exclusive)
    ##GVCFBlock48-49=minGQ=48(inclusive),maxGQ=49(exclusive)
    ##GVCFBlock49-50=minGQ=49(inclusive),maxGQ=50(exclusive)
    ##GVCFBlock5-6=minGQ=5(inclusive),maxGQ=6(exclusive)
    ##GVCFBlock50-51=minGQ=50(inclusive),maxGQ=51(exclusive)
    ##GVCFBlock51-52=minGQ=51(inclusive),maxGQ=52(exclusive)
    ##GVCFBlock52-53=minGQ=52(inclusive),maxGQ=53(exclusive)
    ##GVCFBlock53-54=minGQ=53(inclusive),maxGQ=54(exclusive)
    ##GVCFBlock54-55=minGQ=54(inclusive),maxGQ=55(exclusive)
    ##GVCFBlock55-56=minGQ=55(inclusive),maxGQ=56(exclusive)
    ##GVCFBlock56-57=minGQ=56(inclusive),maxGQ=57(exclusive)
    ##GVCFBlock57-58=minGQ=57(inclusive),maxGQ=58(exclusive)
    ##GVCFBlock58-59=minGQ=58(inclusive),maxGQ=59(exclusive)
    ##GVCFBlock59-60=minGQ=59(inclusive),maxGQ=60(exclusive)
    ##GVCFBlock6-7=minGQ=6(inclusive),maxGQ=7(exclusive)
    ##GVCFBlock60-70=minGQ=60(inclusive),maxGQ=70(exclusive)
    ##GVCFBlock7-8=minGQ=7(inclusive),maxGQ=8(exclusive)
    ##GVCFBlock70-80=minGQ=70(inclusive),maxGQ=80(exclusive)
    ##GVCFBlock8-9=minGQ=8(inclusive),maxGQ=9(exclusive)
    ##GVCFBlock80-90=minGQ=80(inclusive),maxGQ=90(exclusive)
    ##GVCFBlock9-10=minGQ=9(inclusive),maxGQ=10(exclusive)
    ##GVCFBlock90-99=minGQ=90(inclusive),maxGQ=99(exclusive)
    ##GVCFBlock99-100=minGQ=99(inclusive),maxGQ=100(exclusive)
    ##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    ##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
    ##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
    ##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
    ##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
    ##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
    ##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
    ##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
    ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
    ##contig=<ID=20_10037292_10066351,length=29059>
    ##source=HaplotypeCaller
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3278	GT:DP:GQ:MIN_DP:PL	0/0:38:99:37:0,102,1530
    20_10037292_10066351	3279	.	A	<NON_REF>	.	.	END=3279	GT:DP:GQ:MIN_DP:PL	0/0:37:81:37:0,81,1084
    20_10037292_10066351	3280	.	A	<NON_REF>	.	.	END=3281	GT:DP:GQ:MIN_DP:PL	0/0:37:99:37:0,99,1485
    20_10037292_10066351	3282	.	T	<NON_REF>	.	.	END=3282	GT:DP:GQ:MIN_DP:PL	0/0:36:96:36:0,96,1440
    20_10037292_10066351	3283	.	T	<NON_REF>	.	.	END=3283	GT:DP:GQ:MIN_DP:PL	0/0:35:87:35:0,87,1305
    20_10037292_10066351	3284	.	T	<NON_REF>	.	.	END=3284	GT:DP:GQ:MIN_DP:PL	0/0:37:90:37:0,90,1350
    20_10037292_10066351	3285	.	T	<NON_REF>	.	.	END=3293	GT:DP:GQ:MIN_DP:PL	0/0:35:81:31:0,81,1215
    20_10037292_10066351	3294	.	A	<NON_REF>	.	.	END=3302	GT:DP:GQ:MIN_DP:PL	0/0:33:90:30:0,90,970
    20_10037292_10066351	3303	.	G	<NON_REF>	.	.	END=3304	GT:DP:GQ:MIN_DP:PL	0/0:33:72:32:0,72,963
    20_10037292_10066351	3305	.	C	<NON_REF>	.	.	END=3309	GT:DP:GQ:MIN_DP:PL	0/0:34:99:33:0,99,1053
    20_10037292_10066351	3310	.	A	<NON_REF>	.	.	END=3319	GT:DP:GQ:MIN_DP:PL	0/0:35:90:35:0,90,1086
    20_10037292_10066351	3320	.	A	<NON_REF>	.	.	END=3320	GT:DP:GQ:MIN_DP:PL	0/0:33:68:33:0,68,959
    20_10037292_10066351	3321	.	A	<NON_REF>	.	.	END=3322	GT:DP:GQ:MIN_DP:PL	0/0:32:84:32:0,84,1260
    20_10037292_10066351	3323	.	G	<NON_REF>	.	.	END=3323	GT:DP:GQ:MIN_DP:PL	0/0:32:79:32:0,79,953
    20_10037292_10066351	3324	.	T	<NON_REF>	.	.	END=3325	GT:DP:GQ:MIN_DP:PL	0/0:32:81:32:0,81,1215
    20_10037292_10066351	3326	.	G	<NON_REF>	.	.	END=3326	GT:DP:GQ:MIN_DP:PL	0/0:31:60:31:0,60,873
    20_10037292_10066351	3327	.	C	<NON_REF>	.	.	END=3328	GT:DP:GQ:MIN_DP:PL	0/0:30:78:30:0,78,1170
    20_10037292_10066351	3329	.	T	<NON_REF>	.	.	END=3329	GT:DP:GQ:MIN_DP:PL	0/0:31:81:31:0,81,1215
    20_10037292_10066351	3330	.	G	<NON_REF>	.	.	END=3330	GT:DP:GQ:MIN_DP:PL	0/0:31:76:31:0,76,949
    20_10037292_10066351	3331	.	T	<NON_REF>	.	.	END=3332	GT:DP:GQ:MIN_DP:PL	0/0:30:81:29:0,81,1215
    20_10037292_10066351	3333	.	A	<NON_REF>	.	.	END=3335	GT:DP:GQ:MIN_DP:PL	0/0:30:72:30:0,72,892
    20_10037292_10066351	3336	.	T	<NON_REF>	.	.	END=3337	GT:DP:GQ:MIN_DP:PL	0/0:30:84:30:0,84,1260
    20_10037292_10066351	3338	.	C	<NON_REF>	.	.	END=3338	GT:DP:GQ:MIN_DP:PL	0/0:30:59:30:0,59,851
    20_10037292_10066351	3339	.	C	<NON_REF>	.	.	END=3339	GT:DP:GQ:MIN_DP:PL	0/0:30:84:30:0,84,1260
    20_10037292_10066351	3340	.	T	<NON_REF>	.	.	END=3340	GT:DP:GQ:MIN_DP:PL	0/0:30:77:30:0,77,888
    20_10037292_10066351	3341	.	A	<NON_REF>	.	.	END=3343	GT:DP:GQ:MIN_DP:PL	0/0:30:84:28:0,84,910
    20_10037292_10066351	3344	.	T	<NON_REF>	.	.	END=3344	GT:DP:GQ:MIN_DP:PL	0/0:29:73:29:0,73,832
    20_10037292_10066351	3345	.	T	<NON_REF>	.	.	END=3348	GT:DP:GQ:MIN_DP:PL	0/0:29:87:29:0,87,891
    20_10037292_10066351	3349	.	A	<NON_REF>	.	.	END=3349	GT:DP:GQ:MIN_DP:PL	0/0:29:72:29:0,72,904
    20_10037292_10066351	3350	.	G	<NON_REF>	.	.	END=3350	GT:DP:GQ:MIN_DP:PL	0/0:29:87:29:0,87,910
    20_10037292_10066351	3351	.	T	<NON_REF>	.	.	END=3352	GT:DP:GQ:MIN_DP:PL	0/0:30:90:30:0,90,975
    20_10037292_10066351	3353	.	G	<NON_REF>	.	.	END=3354	GT:DP:GQ:MIN_DP:PL	0/0:31:72:30:0,72,846
    20_10037292_10066351	3355	.	A	<NON_REF>	.	.	END=3355	GT:DP:GQ:MIN_DP:PL	0/0:31:93:31:0,93,978
    20_10037292_10066351	3356	.	T	<NON_REF>	.	.	END=3357	GT:DP:GQ:MIN_DP:PL	0/0:31:67:31:0,67,916
    20_10037292_10066351	3358	.	A	<NON_REF>	.	.	END=3363	GT:DP:GQ:MIN_DP:PL	0/0:31:90:31:0,90,1017
    20_10037292_10066351	3364	.	G	<NON_REF>	.	.	END=3364	GT:DP:GQ:MIN_DP:PL	0/0:32:82:32:0,82,947
    20_10037292_10066351	3365	.	A	<NON_REF>	.	.	END=3365	GT:DP:GQ:MIN_DP:PL	0/0:32:90:32:0,90,1350
    20_10037292_10066351	3366	.	C	<NON_REF>	.	.	END=3366	GT:DP:GQ:MIN_DP:PL	0/0:32:79:32:0,79,963
    20_10037292_10066351	3367	.	T	<NON_REF>	.	.	END=3369	GT:DP:GQ:MIN_DP:PL	0/0:32:90:31:0,90,1350
    20_10037292_10066351	3370	.	C	<NON_REF>	.	.	END=3370	GT:DP:GQ:MIN_DP:PL	0/0:32:46:32:0,46,903
    20_10037292_10066351	3371	.	A	<NON_REF>	.	.	END=3371	GT:DP:GQ:MIN_DP:PL	0/0:32:90:32:0,90,1350
    20_10037292_10066351	3372	.	A	<NON_REF>	.	.	END=3372	GT:DP:GQ:MIN_DP:PL	0/0:32:80:32:0,80,905
    20_10037292_10066351	3373	.	A	<NON_REF>	.	.	END=3374	GT:DP:GQ:MIN_DP:PL	0/0:32:90:32:0,90,1350
    20_10037292_10066351	3375	.	G	<NON_REF>	.	.	END=3375	GT:DP:GQ:MIN_DP:PL	0/0:31:76:31:0,76,922
    20_10037292_10066351	3376	.	C	<NON_REF>	.	.	END=3376	GT:DP:GQ:MIN_DP:PL	0/0:33:93:33:0,93,1395
    20_10037292_10066351	3377	.	A	<NON_REF>	.	.	END=3381	GT:DP:GQ:MIN_DP:PL	0/0:32:84:31:0,84,1260
    20_10037292_10066351	3382	.	A	<NON_REF>	.	.	END=3385	GT:DP:GQ:MIN_DP:PL	0/0:33:90:33:0,90,1350
    20_10037292_10066351	3386	.	A	<NON_REF>	.	.	END=3387	GT:DP:GQ:MIN_DP:PL	0/0:34:84:33:0,84,964
    20_10037292_10066351	3388	.	A	<NON_REF>	.	.	END=3397	GT:DP:GQ:MIN_DP:PL	0/0:32:90:31:0,90,1350
    20_10037292_10066351	3398	.	A	<NON_REF>	.	.	END=3398	GT:DP:GQ:MIN_DP:PL	0/0:31:75:31:0,75,920
    20_10037292_10066351	3399	.	T	<NON_REF>	.	.	END=3399	GT:DP:GQ:MIN_DP:PL	0/0:31:87:31:0,87,1305
    20_10037292_10066351	3400	.	T	<NON_REF>	.	.	END=3400	GT:DP:GQ:MIN_DP:PL	0/0:32:90:32:0,90,1350
    20_10037292_10066351	3401	.	T	<NON_REF>	.	.	END=3402	GT:DP:GQ:MIN_DP:PL	0/0:32:87:31:0,87,1305
    20_10037292_10066351	3403	.	T	<NON_REF>	.	.	END=3403	GT:DP:GQ:MIN_DP:PL	0/0:32:90:32:0,90,1350
    20_10037292_10066351	3404	.	C	<NON_REF>	.	.	END=3405	GT:DP:GQ:MIN_DP:PL	0/0:32:80:32:0,80,944
    20_10037292_10066351	3406	.	G	<NON_REF>	.	.	END=3407	GT:DP:GQ:MIN_DP:PL	0/0:31:72:30:0,72,859
    20_10037292_10066351	3408	.	G	<NON_REF>	.	.	END=3408	GT:DP:GQ:MIN_DP:PL	0/0:32:62:32:0,62,890
    20_10037292_10066351	3409	.	T	<NON_REF>	.	.	END=3418	GT:DP:GQ:MIN_DP:PL	0/0:33:81:30:0,81,1215
    20_10037292_10066351	3419	.	T	<NON_REF>	.	.	END=3419	GT:DP:GQ:MIN_DP:PL	0/0:29:74:29:0,74,827
    20_10037292_10066351	3420	.	T	<NON_REF>	.	.	END=3421	GT:DP:GQ:MIN_DP:PL	0/0:30:84:29:0,84,1260
    20_10037292_10066351	3422	.	T	<NON_REF>	.	.	END=3430	GT:DP:GQ:MIN_DP:PL	0/0:29:75:28:0,75,1125
    20_10037292_10066351	3431	.	T	<NON_REF>	.	.	END=3439	GT:DP:GQ:MIN_DP:PL	0/0:28:81:28:0,81,1215
    20_10037292_10066351	3440	.	A	<NON_REF>	.	.	END=3440	GT:DP:GQ:MIN_DP:PL	0/0:28:70:28:0,70,782
    20_10037292_10066351	3441	.	T	<NON_REF>	.	.	END=3442	GT:DP:GQ:MIN_DP:PL	0/0:28:81:28:0,81,1215
    20_10037292_10066351	3443	.	T	<NON_REF>	.	.	END=3443	GT:DP:GQ:MIN_DP:PL	0/0:28:78:28:0,78,1170
    20_10037292_10066351	3444	.	T	<NON_REF>	.	.	END=3445	GT:DP:GQ:MIN_DP:PL	0/0:28:64:28:0,64,722
    20_10037292_10066351	3446	.	G	<NON_REF>	.	.	END=3446	GT:DP:GQ:MIN_DP:PL	0/0:28:78:28:0,78,1170
    20_10037292_10066351	3447	.	C	<NON_REF>	.	.	END=3447	GT:DP:GQ:MIN_DP:PL	0/0:29:53:29:0,53,694
    20_10037292_10066351	3448	.	T	<NON_REF>	.	.	END=3449	GT:DP:GQ:MIN_DP:PL	0/0:31:76:30:0,76,827
    20_10037292_10066351	3450	.	A	<NON_REF>	.	.	END=3450	GT:DP:GQ:MIN_DP:PL	0/0:31:87:31:0,87,1305
    20_10037292_10066351	3451	.	C	<NON_REF>	.	.	END=3452	GT:DP:GQ:MIN_DP:PL	0/0:31:74:31:0,74,715
    20_10037292_10066351	3453	.	T	<NON_REF>	.	.	END=3455	GT:DP:GQ:MIN_DP:PL	0/0:31:84:31:0,84,1260
    20_10037292_10066351	3456	.	A	<NON_REF>	.	.	END=3456	GT:DP:GQ:MIN_DP:PL	0/0:31:23:31:0,23,766
    20_10037292_10066351	3457	.	T	<NON_REF>	.	.	END=3460	GT:DP:GQ:MIN_DP:PL	0/0:31:90:31:0,90,1350
    20_10037292_10066351	3461	.	C	<NON_REF>	.	.	END=3461	GT:DP:GQ:MIN_DP:PL	0/0:30:89:30:0,89,873
    20_10037292_10066351	3462	.	T	<NON_REF>	.	.	END=3462	GT:DP:GQ:MIN_DP:PL	0/0:31:90:31:0,90,1350
    20_10037292_10066351	3463	.	G	<NON_REF>	.	.	END=3463	GT:DP:GQ:MIN_DP:PL	0/0:31:44:31:0,44,739
    20_10037292_10066351	3464	.	T	<NON_REF>	.	.	END=3468	GT:DP:GQ:MIN_DP:PL	0/0:32:90:32:0,90,1350
    20_10037292_10066351	3469	.	C	<NON_REF>	.	.	END=3469	GT:DP:GQ:MIN_DP:PL	0/0:32:79:32:0,79,816
    20_10037292_10066351	3470	.	T	<NON_REF>	.	.	END=3470	GT:DP:GQ:MIN_DP:PL	0/0:31:84:31:0,84,1260
    20_10037292_10066351	3471	.	T	<NON_REF>	.	.	END=3478	GT:DP:GQ:MIN_DP:PL	0/0:32:75:32:0,75,1125
    20_10037292_10066351	3479	.	T	<NON_REF>	.	.	END=3479	GT:DP:GQ:MIN_DP:PL	0/0:34:36:34:0,36,906
    20_10037292_10066351	3480	.	C	CT,<NON_REF>	503.03	.	DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23	GT:AD:DP:GQ:PL:SB	1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
    20_10037292_10066351	3481	.	T	<NON_REF>	.	.	END=3481	GT:DP:GQ:MIN_DP:PL	0/0:21:51:21:0,51,765
    20_10037292_10066351	3482	.	T	<NON_REF>	.	.	END=3482	GT:DP:GQ:MIN_DP:PL	0/0:21:54:21:0,54,810
    20_10037292_10066351	3483	.	T	<NON_REF>	.	.	END=3487	GT:DP:GQ:MIN_DP:PL	0/0:20:51:19:0,51,765
    20_10037292_10066351	3488	.	T	<NON_REF>	.	.	END=3488	GT:DP:GQ:MIN_DP:PL	0/0:19:42:19:0,42,571
    20_10037292_10066351	3489	.	A	<NON_REF>	.	.	END=3489	GT:DP:GQ:MIN_DP:PL	0/0:17:51:17:0,51,521
    20_10037292_10066351	3490	.	C	<NON_REF>	.	.	END=3490	GT:DP:GQ:MIN_DP:PL	0/0:17:35:17:0,35,431
    20_10037292_10066351	3491	.	A	<NON_REF>	.	.	END=3495	GT:DP:GQ:MIN_DP:PL	0/0:17:48:17:0,48,720
    20_10037292_10066351	3496	.	A	<NON_REF>	.	.	END=3498	GT:DP:GQ:MIN_DP:PL	0/0:17:51:17:0,51,473
    20_10037292_10066351	3499	.	C	<NON_REF>	.	.	END=3499	GT:DP:GQ:MIN_DP:PL	0/0:16:48:16:0,48,428
    20_10037292_10066351	3500	.	G	<NON_REF>	.	.	END=3500	GT:DP:GQ:MIN_DP:PL	0/0:16:31:16:0,31,379
    20_10037292_10066351	3501	.	T	<NON_REF>	.	.	END=3501	GT:DP:GQ:MIN_DP:PL	0/0:17:48:17:0,48,720
    20_10037292_10066351	3502	.	A	<NON_REF>	.	.	END=3503	GT:DP:GQ:MIN_DP:PL	0/0:19:54:18:0,54,550
    20_10037292_10066351	3504	.	T	<NON_REF>	.	.	END=3504	GT:DP:GQ:MIN_DP:PL	0/0:19:48:19:0,48,720
    20_10037292_10066351	3505	.	A	<NON_REF>	.	.	END=3506	GT:DP:GQ:MIN_DP:PL	0/0:20:51:20:0,51,765
    20_10037292_10066351	3507	.	T	<NON_REF>	.	.	END=3519	GT:DP:GQ:MIN_DP:PL	0/0:19:54:18:0,54,501
    20_10037292_10066351	3520	.	AT	A,<NON_REF>	609.03	.	DP=18;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=64800,18	GT:AD:DP:GQ:PL:SB	1/1:0,18,0:18:54:623,54,0,623,54,623:0,0,9,9
    20_10037292_10066351	3522	.	T	<NON_REF>	.	.	END=3525	GT:DP:GQ:MIN_DP:PL	0/0:18:54:18:0,54,550
    20_10037292_10066351	3526	.	T	<NON_REF>	.	.	END=3527	GT:DP:GQ:MIN_DP:PL	0/0:19:57:19:0,57,607
    20_10037292_10066351	3528	.	T	<NON_REF>	.	.	END=3528	GT:DP:GQ:MIN_DP:PL	0/0:19:54:19:0,54,810
    20_10037292_10066351	3529	.	T	A,<NON_REF>	155.64	.	BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;MLEAC=1,0;MLEAF=0.500,0.00;MQRankSum=0.000;RAW_MQandDP=75600,21;ReadPosRankSum=-1.158	GT:AD:DP:GQ:PL:SB	0/1:12,8,0:20:99:163,0,328,199,352,551:5,7,5,3
    20_10037292_10066351	3530	.	A	<NON_REF>	.	.	END=3530	GT:DP:GQ:MIN_DP:PL	0/0:33:64:33:0,64,941
    20_10037292_10066351	3531	.	A	<NON_REF>	.	.	END=3533	GT:DP:GQ:MIN_DP:PL	0/0:33:81:33:0,81,1215
    20_10037292_10066351	3534	.	A	<NON_REF>	.	.	END=3534	GT:DP:GQ:MIN_DP:PL	0/0:33:78:33:0,78,1170
    20_10037292_10066351	3535	.	A	<NON_REF>	.	.	END=3536	GT:DP:GQ:MIN_DP:PL	0/0:33:68:33:0,68,891
    20_10037292_10066351	3537	.	A	<NON_REF>	.	.	END=3546	GT:DP:GQ:MIN_DP:PL	0/0:29:72:26:0,72,1080
    ```

हमने एक बार फिर अंतिम हेडर लाइन को हाइलाइट किया है, साथ ही फ़ाइल में पहले तीन 'उचित' वेरिएंट कॉल को भी।

तुम देखोगे कि वेरिएंट कॉल लाइनें कई गैर-वेरिएंट लाइनों के बीच बिखरी हुई हैं, जो गैर-वेरिएंट क्षेत्रों का प्रतिनिधित्व करती हैं जहां वेरिएंट कॉलर को भिन्नता का कोई प्रमाण नहीं मिला।
जैसा कि ऊपर संक्षेप में उल्लेख किया गया है, यह वेरिएंट कॉलिंग के GVCF मोड में विशेष है: वेरिएंट कॉलर कुछ सांख्यिकी को कैप्चर करता है जो भिन्नता की अनुपस्थिति में इसके विश्वास के स्तर का वर्णन करती हैं।
यह दो बहुत अलग मामलों के बीच अंतर करना संभव बनाता है: (1) अच्छी गुणवत्ता वाला डेटा है जो दिखाता है कि नमूना होमोजाइगस-रेफरेंस है, और (2) किसी भी तरह से निर्धारण करने के लिए पर्याप्त अच्छा डेटा उपलब्ध नहीं है।

इस तरह के GVCF में, आमतौर पर ऐसी कई गैर-वेरिएंट लाइनें होती हैं, जिनके बीच कम संख्या में वेरिएंट रिकॉर्ड बिखरे होते हैं।

#### 2.2.3. अन्य दो नमूनों पर प्रक्रिया को दोहराएं

अब आओ बाकी दो नमूनों के लिए GVCFs जेनरेट करने के लिए नीचे दी गई कमांड्स चलाएं, एक के बाद एक।

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

??? success "कमांड आउटपुट"

    ```console
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_father.bam -O reads_father.g.vcf -L /data/ref/intervals.bed -ERC GVCF
    17:28:30.677 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    17:28:30.801 INFO  HaplotypeCaller - ------------------------------------------------------------
    17:28:30.803 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    17:28:30.804 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    17:28:30.804 INFO  HaplotypeCaller - Executing as root@be1a0302f6c7 on Linux v6.8.0-1030-azure amd64
    17:28:30.804 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    17:28:30.804 INFO  HaplotypeCaller - Start Date/Time: February 11, 2026 at 5:28:30 PM GMT
    17:28:30.804 INFO  HaplotypeCaller - ------------------------------------------------------------
    17:28:30.804 INFO  HaplotypeCaller - ------------------------------------------------------------
    17:28:30.805 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    17:28:30.805 INFO  HaplotypeCaller - Picard Version: 3.1.1
    17:28:30.805 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    17:28:30.806 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    17:28:30.806 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    17:28:30.806 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    17:28:30.806 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    17:28:30.806 INFO  HaplotypeCaller - Deflater: IntelDeflater
    17:28:30.807 INFO  HaplotypeCaller - Inflater: IntelInflater
    17:28:30.807 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    17:28:30.807 INFO  HaplotypeCaller - Requester pays: disabled
    17:28:30.807 INFO  HaplotypeCaller - Initializing engine
    17:28:30.933 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    17:28:30.946 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    17:28:30.951 INFO  HaplotypeCaller - Done initializing engine
    17:28:30.953 INFO  HaplotypeCallerEngine - Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled
    17:28:30.957 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    17:28:30.959 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    17:28:30.960 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    17:28:30.963 INFO  HaplotypeCallerEngine - Standard Emitting and Calling confidence set to -0.0 for reference-model confidence output
    17:28:30.963 INFO  HaplotypeCallerEngine - All sites annotated with PLs forced to true for reference-model confidence output
    17:28:30.972 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    17:28:30.987 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    17:28:30.987 INFO  IntelPairHmm - Available threads: 4
    17:28:30.987 INFO  IntelPairHmm - Requested threads: 4
    17:28:30.987 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    17:28:31.034 INFO  ProgressMeter - Starting traversal
    17:28:31.034 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    17:28:31.570 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    17:28:32.865 INFO  HaplotypeCaller - 9 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    9 total reads filtered out of 2064 reads processed
    17:28:32.866 INFO  ProgressMeter - 20_10037292_10066351:13338              0.0                    38           1245.2
    17:28:32.866 INFO  ProgressMeter - Traversal complete. Processed 38 total regions in 0.0 minutes.
    17:28:32.868 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.0035923200000000004
    17:28:32.868 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.10765202500000001
    17:28:32.868 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.03 sec
    17:28:32.869 INFO  HaplotypeCaller - Shutting down engine
    [February 11, 2026 at 5:28:32 PM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.04 minutes.
    Runtime.totalMemory()=299892736
    ```

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

??? success "कमांड आउटपुट"

    ```console
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_son.bam -O reads_son.g.vcf -L /data/ref/intervals.bed -ERC GVCF
    17:30:10.017 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    17:30:10.156 INFO  HaplotypeCaller - ------------------------------------------------------------
    17:30:10.159 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    17:30:10.159 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    17:30:10.159 INFO  HaplotypeCaller - Executing as root@be1a0302f6c7 on Linux v6.8.0-1030-azure amd64
    17:30:10.159 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    17:30:10.159 INFO  HaplotypeCaller - Start Date/Time: February 11, 2026 at 5:30:09 PM GMT
    17:30:10.159 INFO  HaplotypeCaller - ------------------------------------------------------------
    17:30:10.160 INFO  HaplotypeCaller - ------------------------------------------------------------
    17:30:10.160 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    17:30:10.160 INFO  HaplotypeCaller - Picard Version: 3.1.1
    17:30:10.161 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    17:30:10.161 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    17:30:10.161 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    17:30:10.161 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    17:30:10.161 INFO  HaplotypeCaller -HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    17:30:10.161 INFO  HaplotypeCaller - Deflater: IntelDeflater
    17:30:10.162 INFO  HaplotypeCaller - Inflater: IntelInflater
    17:30:10.162 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    17:30:10.162 INFO  HaplotypeCaller - Requester pays: disabled
    17:30:10.162 INFO  HaplotypeCaller - Initializing engine
    17:30:10.277 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    17:30:10.290 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    17:30:10.296 INFO  HaplotypeCaller - Done initializing engine
    17:30:10.298 INFO  HaplotypeCallerEngine - Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled
    17:30:10.302 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    17:30:10.303 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    17:30:10.304 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    17:30:10.307 INFO  HaplotypeCallerEngine - Standard Emitting and Calling confidence set to -0.0 for reference-model confidence output
    17:30:10.307 INFO  HaplotypeCallerEngine - All sites annotated with PLs forced to true for reference-model confidence output
    17:30:10.315 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    17:30:10.328 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    17:30:10.329 INFO  IntelPairHmm - Available threads: 4
    17:30:10.329 INFO  IntelPairHmm - Requested threads: 4
    17:30:10.329 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    17:30:10.368 INFO  ProgressMeter - Starting traversal
    17:30:10.369 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    17:30:10.875 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    17:30:11.980 INFO  HaplotypeCaller - 14 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    14 total reads filtered out of 1981 reads processed
    17:30:11.981 INFO  ProgressMeter - 20_10037292_10066351:13223              0.0                    35           1302.7
    17:30:11.981 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    17:30:11.983 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.0034843710000000004
    17:30:11.983 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.048108363
    17:30:11.983 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    17:30:11.984 INFO  HaplotypeCaller - Shutting down engine
    [February 11, 2026 at 5:30:11 PM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.03 minutes.
    Runtime.totalMemory()=226492416
    ```

एक बार यह पूरा हो जाने के बाद, तुम्हारे पास अपनी वर्तमान डायरेक्टरी में `.g.vcf` में समाप्त होने वाली तीन फ़ाइलें (प्रति नमूना एक) और उनकी संबंधित इंडेक्स फ़ाइलें `.g.vcf.idx` में समाप्त होनी चाहिए।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    conda.yml        reads_father.g.vcf      reads_mother.g.vcf      reads_son.g.vcf
    hsperfdata_root  reads_father.g.vcf.idx  reads_mother.g.vcf.idx  reads_son.g.vcf.idx
    ```

इस बिंदु पर, हमने अपने प्रत्येक इनपुट नमूने के लिए GVCF मोड में वेरिएंट कॉल किए हैं।
अब संयुक्त कॉलिंग पर आगे बढ़ने का समय है।

लेकिन कंटेनर से बाहर मत निकलो!
हम अगले चरण में उसी का उपयोग करने जा रहे हैं।

### 2.3. संयुक्त जीनोटाइपिंग चलाएं

अब जब हमारे पास सभी GVCFs हैं, तो हम नमूनों के कोहॉर्ट के लिए वेरिएंट कॉल जेनरेट करने के लिए संयुक्त जीनोटाइपिंग दृष्टिकोण को आज़मा सकते हैं।
यह एक दो-चरणीय विधि है जिसमें सभी GVCFs से डेटा को डेटा स्टोर में जोड़ना शामिल है, फिर संयुक्त-कॉल किए गए वेरिएंट के अंतिम VCF को जेनरेट करने के लिए वास्तविक संयुक्त जीनोटाइपिंग विश्लेषण चलाना।

#### 2.3.1. सभी प्रति-नमूना GVCFs को जोड़ें

यह पहला चरण एक अन्य GATK टूल का उपयोग करता है, जिसे GenomicsDBImport कहा जाता है, सभी GVCFs से डेटा को GenomicsDB डेटा स्टोर में जोड़ने के लिए।
GenomicsDB डेटा स्टोर एक प्रकार का डेटाबेस फॉर्मेट है जो वेरिएंट जानकारी के लिए एक मध्यवर्ती भंडारण के रूप में कार्य करता है।

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

??? success "कमांड आउटपुट"

    ```console hl_lines="33 36 37 39 40"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenomicsDBImport -V reads_mother.g.vcf -V reads_father.g.vcf -V reads_son.g.vcf -L /data/ref/intervals.bed --genomicsdb-workspace-path family_trio_gdb
    17:37:07.569 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    17:37:07.699 INFO  GenomicsDBImport - ------------------------------------------------------------
    17:37:07.702 INFO  GenomicsDBImport - The Genome Analysis Toolkit (GATK) v4.5.0.0
    17:37:07.702 INFO  GenomicsDBImport - For support and documentation go to https://software.broadinstitute.org/gatk/
    17:37:07.703 INFO  GenomicsDBImport - Executing as root@be1a0302f6c7 on Linux v6.8.0-1030-azure amd64
    17:37:07.703 INFO  GenomicsDBImport - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    17:37:07.704 INFO  GenomicsDBImport - Start Date/Time: February 11, 2026 at 5:37:07 PM GMT
    17:37:07.704 INFO  GenomicsDBImport - ------------------------------------------------------------
    17:37:07.704 INFO  GenomicsDBImport - ------------------------------------------------------------
    17:37:07.706 INFO  GenomicsDBImport - HTSJDK Version: 4.1.0
    17:37:07.706 INFO  GenomicsDBImport - Picard Version: 3.1.1
    17:37:07.707 INFO  GenomicsDBImport - Built for Spark Version: 3.5.0
    17:37:07.709 INFO  GenomicsDBImport - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    17:37:07.709 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    17:37:07.709 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    17:37:07.710 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    17:37:07.710 INFO  GenomicsDBImport - Deflater: IntelDeflater
    17:37:07.711 INFO  GenomicsDBImport - Inflater: IntelInflater
    17:37:07.711 INFO  GenomicsDBImport - GCS max retries/reopens: 20
    17:37:07.711 INFO  GenomicsDBImport - Requester pays: disabled
    17:37:07.712 INFO  GenomicsDBImport - Initializing engine
    17:37:07.883 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    17:37:07.886 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    17:37:07.889 INFO  GenomicsDBImport - Done initializing engine
    17:37:08.560 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    17:37:08.561 INFO  GenomicsDBImport - Vid Map JSON file will be written to /tmp/family_trio_gdb/vidmap.json
    17:37:08.561 INFO  GenomicsDBImport - Callset Map JSON file will be written to /tmp/family_trio_gdb/callset.json
    17:37:08.561 INFO  GenomicsDBImport - Complete VCF Header will be written to /tmp/family_trio_gdb/vcfheader.vcf
    17:37:08.561 INFO  GenomicsDBImport - Importing to workspace - /tmp/family_trio_gdb
    17:37:08.878 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    17:37:09.359 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    17:37:09.487 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    17:37:09.591 INFO  GenomicsDBImport - Done importing batch 1/1
    17:37:09.592 INFO  GenomicsDBImport - Import completed!
    17:37:09.592 INFO  GenomicsDBImport - Shutting down engine
    [February 11, 2026 at 5:37:09 PM GMT] org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport done. Elapsed time: 0.03 minutes.
    Runtime.totalMemory()=113246208
    Tool returned:
    true
    ```

इस चरण का आउटपुट प्रभावी रूप से एक डायरेक्टरी है जिसमें कई अलग-अलग फ़ाइलों के रूप में संयुक्त वेरिएंट डेटा रखने वाली आगे नेस्टेड डायरेक्टरियों का एक सेट होता है।
तुम इसके चारों ओर घूम सकते हो लेकिन तुम जल्दी ही देखोगे कि यह डेटा स्टोर फॉर्मेट मनुष्यों द्वारा सीधे पढ़े जाने के लिए नहीं है।

!!! tip "सुझाव"

    GATK में ऐसे टूल शामिल हैं जो आवश्यकतानुसार डेटा स्टोर से वेरिएंट कॉल डेटा का निरीक्षण और निष्कर्षण संभव बनाते हैं।

#### 2.3.2. वास्तविक संयुक्त जीनोटाइपिंग विश्लेषण चलाएं

यह दूसरा चरण एक और GATK टूल का उपयोग करता है, जिसे GenotypeGVCFs कहा जाता है, कोहॉर्ट में सभी नमूनों में उपलब्ध डेटा के आलोक में वेरिएंट सांख्यिकी और व्यक्तिगत जीनोटाइप्स की पुनर्गणना करने के लिए।

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

??? success "कमांड आउटपुट"

    ```console hl_lines="30 35 37 38"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenotypeGVCFs -R /data/ref/ref.fasta -V gendb://family_trio_gdb -O family_trio.vcf
    17:38:45.084 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    17:38:45.217 INFO  GenotypeGVCFs - ------------------------------------------------------------
    17:38:45.220 INFO  GenotypeGVCFs - The Genome Analysis Toolkit (GATK) v4.5.0.0
    17:38:45.220 INFO  GenotypeGVCFs - For support and documentation go to https://software.broadinstitute.org/gatk/
    17:38:45.220 INFO  GenotypeGVCFs - Executing as root@be1a0302f6c7 on Linux v6.8.0-1030-azure amd64
    17:38:45.220 INFO  GenotypeGVCFs - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    17:38:45.221 INFO  GenotypeGVCFs - Start Date/Time: February 11, 2026 at 5:38:45 PM GMT
    17:38:45.221 INFO  GenotypeGVCFs - ------------------------------------------------------------
    17:38:45.221 INFO  GenotypeGVCFs - ------------------------------------------------------------
    17:38:45.221 INFO  GenotypeGVCFs - HTSJDK Version: 4.1.0
    17:38:45.222 INFO  GenotypeGVCFs - Picard Version: 3.1.1
    17:38:45.222 INFO  GenotypeGVCFs - Built for Spark Version: 3.5.0
    17:38:45.222 INFO  GenotypeGVCFs - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    17:38:45.222 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    17:38:45.222 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    17:38:45.222 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    17:38:45.223 INFO  GenotypeGVCFs - Deflater: IntelDeflater
    17:38:45.223 INFO  GenotypeGVCFs - Inflater: IntelInflater
    17:38:45.223 INFO  GenotypeGVCFs - GCS max retries/reopens: 20
    17:38:45.223 INFO  GenotypeGVCFs - Requester pays: disabled
    17:38:45.223 INFO  GenotypeGVCFs - Initializing engine
    17:38:45.544 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    17:38:45.561 INFO  NativeGenomicsDB - pid=221 tid=222 No valid combination operation found for INFO field InbreedingCoeff  - the field will NOT be part of INFO fields in the generated VCF records
    17:38:45.561 INFO  NativeGenomicsDB - pid=221 tid=222 No valid combination operation found for INFO field MLEAC  - the field will NOT be part of INFO fields in the generated VCF records
    17:38:45.561 INFO  NativeGenomicsDB - pid=221 tid=222 No valid combination operation found for INFO field MLEAF  - the field will NOT be part of INFO fields in the generated VCF records
    17:38:45.577 INFO  GenotypeGVCFs - Done initializing engine
    17:38:45.615 INFO  ProgressMeter - Starting traversal
    17:38:45.615 INFO  ProgressMeter -        Current Locus  Elapsed Minutes    Variants Processed  Variants/Minute
    17:38:45.903 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    GENOMICSDB_TIMER,GenomicsDB iterator next() timer,Wall-clock time(s),0.07757032800000006,Cpu time(s),0.07253379200000037
    17:38:46.421 INFO  ProgressMeter - 20_10037292_10066351:13953              0.0                  3390         252357.3
    17:38:46.422 INFO  ProgressMeter - Traversal complete. Processed 3390 total variants in 0.0 minutes.
    17:38:46.423 INFO  GenotypeGVCFs - Shutting down engine
    [February 11, 2026 at 5:38:46 PM GMT] org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs done. Elapsed time: 0.02 minutes.
    Runtime.totalMemory()=203423744
    ```

यह कंटेनर में वर्तमान कार्यशील डायरेक्टरी में VCF आउटपुट फ़ाइल `family_trio.vcf` बनाता है, साथ ही इसकी इंडेक्स, `family_trio.vcf.idx`।
यह एक और उचित रूप से छोटी फ़ाइल है, इसलिए तुम फ़ाइल की सामग्री देखने के लिए `cat family_trio.vcf` चला सकते हो, और पहली कुछ वेरिएंट लाइनों को खोजने के लिए नीचे स्क्रॉल कर सकते हो।

??? abstract "फ़ाइल की सामग्री (संक्षिप्त)"

    ```console title="family_trio.vcf" linenums="1" hl_lines="39"
    ##fileformat=VCFv4.2
    ##ALT=<ID=NON_REF,Description="Represents any possible alternative allele not already represented at this location by REF and ALT">
    ##FILTER=<ID=LowQual,Description="Low quality">
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
    ##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles">
    ##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
    ##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">
    ##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
    ##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
    ##GATKCommandLine=<ID=GenomicsDBImport,CommandLine="GenomicsDBImport --genomicsdb-workspace-path family_trio_gdb --variant reads_mother.g.vcf --variant reads_father.g.vcf --variant reads_son.g.vcf --intervals /data/ref/intervals.bed [abridged]",Version="4.5.0.0",Date="February 11, 2026 at 5:37:07 PM GMT">
    ##GATKCommandLine=<ID=GenotypeGVCFs,CommandLine="GenotypeGVCFs --output family_trio.vcf --variant gendb://family_trio_gdb --reference /data/ref/ref.fasta --include-non-variant-sites false [abridged]",Version="4.5.0.0",Date="February 11, 2026 at 5:38:45 PM GMT">
    ##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller --emit-ref-confidence GVCF --output reads_mother.g.vcf --intervals /data/ref/intervals.bed --input /data/bam/reads_mother.bam --reference /data/ref/ref.fasta [abridged]",Version="4.5.0.0",Date="February 11, 2026 at 4:51:00 PM GMT">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    ##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
    ##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
    ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
    ##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
    ##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
    ##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
    ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
    ##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
    ##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
    ##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
    ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
    ##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
    ##contig=<ID=20_10037292_10066351,length=29059>
    ##source=GenomicsDBImport
    ##source=GenotypeGVCFs
    ##source=HaplotypeCaller
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
    20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
    20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730  GT:AD:DP:GQ:PL  0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
    20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
    20_10037292_10066351    4012    .       C       T       3950.73 .       AC=6;AF=1.00;AN=6;DP=127;ExcessHet=0.0000;FS=0.000;MLEAC=6;MLEAF=1.00;MQ=60.00;QD=31.86;SOR=0.725    GT:AD:DP:GQ:PL  1/1:0,46:46:99:1446,137,0   1/1:0,43:43:99:1412,129,0        1/1:0,35:35:99:1106,105,0
    20_10037292_10066351    4409    .       A       ATATG   2478.69 .       AC=6;AF=1.00;AN=6;DP=96;ExcessHet=0.0000;FS=0.000;MLEAC=6;MLEAF=1.00;MQ=60.00;QD=33.95;SOR=0.963     GT:AD:DP:GQ:PL  1/1:0,28:28:90:969,90,0 1/1:0,21:21:69:724,69,0      1/1:0,24:24:72:799,72,0
    20_10037292_10066351    4464    .       T       TA      620.25  .       AC=1;AF=0.167;AN=6;BaseQRankSum=0.108;DP=102;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=19.38;ReadPosRankSum=1.27;SOR=0.892 GT:AD:DP:GQ:PGT:PID:PL:PS       0|1:15,17:32:99:0|1:4464_T_TA:629,0,554:4464    0/0:30,0:30:78:.:.:0,78,1170 0/0:39,0:39:99:.:.:0,108,1286
    20_10037292_10066351    4465    .       T       TA      620.25  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-2.250e-01;DP=101;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=19.38;ReadPosRankSum=0.910;SOR=0.892   GT:AD:DP:GQ:PGT:PID:PL:PS       0|1:15,17:32:99:0|1:4464_T_TA:629,0,554:4464    0/0:30,0:30:78:.:.:0,78,1170 0/0:39,0:39:99:.:.:0,108,1286
    20_10037292_10066351    5027    .       C       T       3339.73 .       AC=6;AF=1.00;AN=6;DP=108;ExcessHet=0.0000;FS=0.000;MLEAC=6;MLEAF=1.00;MQ=60.00;QD=31.51;SOR=0.731    GT:AD:DP:GQ:PL  1/1:0,36:36:99:1164,108,0   1/1:0,26:26:77:798,77,0  1/1:0,44:44:99:1391,132,0
    20_10037292_10066351    5469    .       A       G       2725.93 .       AC=5;AF=0.833;AN=6;BaseQRankSum=-3.665e+00;DP=113;ExcessHet=0.0000;FS=6.914;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=24.34;ReadPosRankSum=1.50;SOR=0.320    GT:AD:DP:GQ:PL  0/1:18,23:41:99:553,0,486       1/1:0,42:42:99:1311,126,0       1/1:0,29:29:86:876,86,0
    20_10037292_10066351    7557    .       A       G       2257.93 .       AC=5;AF=0.833;AN=6;BaseQRankSum=-1.362e+00;DP=111;ExcessHet=0.0000;FS=3.400;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.50;ReadPosRankSum=1.11;SOR=0.566    GT:AD:DP:GQ:PL  0/1:19,15:34:99:313,0,493       1/1:0,34:34:99:949,100,0        1/1:0,37:37:99:1010,108,0
    20_10037292_10066351    7786    .       G       T       3503.73 .       AC=6;AF=1.00;AN=6;DP=114;ExcessHet=0.0000;FS=0.000;MLEAC=6;MLEAF=1.00;MQ=60.00;QD=31.28;SOR=0.970    GT:AD:DP:GQ:PL  1/1:0,34:34:99:1066,102,0   1/1:0,34:34:99:1057,102,0        1/1:0,44:44:99:1394,132,0
    20_10037292_10066351    8350    .       G       C       2663.93 .       AC=5;AF=0.833;AN=6;BaseQRankSum=-1.608e+00;DP=106;ExcessHet=0.0000;FS=5.378;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=25.37;ReadPosRankSum=-1.870e-01;SOR=0.950      GT:AD:DP:GQ:PL  0/1:16,14:30:99:356,0,430       1/1:0,39:39:99:1176,115,0       1/1:0,36:36:99:1146,108,0
    20_10037292_10066351    8886    .       AAGAAAGAAAG     A       3037.69 .       AC=6;AF=1.00;AN=6;DP=89;ExcessHet=0.0000;FS=0.000;MLEAC=6;MLEAF=1.00;MQ=60.00;QD=25.36;SOR=2.269     GT:AD:DP:GQ:PL  1/1:0,18:18:55:804,55,0      1/1:0,29:29:88:1282,88,0        1/1:0,22:22:67:965,67,0
    20_10037292_10066351    9536    .       T       C       1089.95 .       AC=3;AF=0.500;AN=6;BaseQRankSum=-5.640e-01;DP=82;ExcessHet=0.0000;FS=12.258;MLEAC=3;MLEAF=0.500;MQ=60.00;MQRankSum=0.00;QD=20.57;ReadPosRankSum=0.860;SOR=0.373   GT:AD:DP:GQ:PL  1/1:0,32:32:95:950,95,0 0/0:29,0:29:81:0,81,1215        0/1:14,7:21:99:156,0,353
    20_10037292_10066351    13375   .       C       T       724.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=0.171;DP=121;ExcessHet=0.0000;FS=7.398;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=12.71;ReadPosRankSum=0.415;SOR=1.688        GT:AD:DP:GQ:PL  0/1:28,29:57:99:733,0,679       0/0:29,0:29:81:0,81,1215        0/0:34,0:34:99:0,99,1485
    20_10037292_10066351    13536   .       T       C       1025.16 .       AC=2;AF=0.333;AN=6;BaseQRankSum=1.63;DP=118;ExcessHet=0.9691;FS=1.719;MLEAC=2;MLEAF=0.333;MQ=60.00;MQRankSum=0.00;QD=11.65;ReadPosRankSum=-2.000e-01;SOR=0.904    GT:AD:DP:GQ:PL  0/1:21,23:44:99:591,0,526       0/1:26,18:44:99:445,0,672       0/0:29,0:29:84:0,84,1260
    20_10037292_10066351    14156   .       T       C       438.16  .       AC=2;AF=0.333;AN=6;BaseQRankSum=3.20;DP=96;ExcessHet=0.9691;FS=2.381;MLEAC=2;MLEAF=0.333;MQ=60.00;MQRankSum=0.00;QD=7.82;ReadPosRankSum=1.13;SOR=0.592    GT:AD:DP:GQ:PL  0/1:25,11:36:99:258,0,676       0/1:12,8:20:99:191,0,319        0/0:38,0:38:99:0,99,1117
    20_10037292_10066351    14403   .       G       A       144.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=2.63;DP=116;ExcessHet=0.0000;FS=1.435;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=3.52;ReadPosRankSum=0.252;SOR=0.802  GT:AD:DP:GQ:PL  0/1:32,9:41:99:153,0,821        0/0:37,0:37:99:0,109,1169       0/0:37,0:37:99:0,99,1113
    ```

हमने एक बार फिर अंतिम हेडर लाइन को हाइलाइट किया है, जो वेरिएंट कॉल डेटा की शुरुआत को चिह्नित करती है।

यह पहले जेनरेट किए गए VCF के समान दिखता है, सिवाय इस बार हमारे पास सभी तीन नमूनों के लिए जीनोटाइप-स्तरीय जानकारी है।
फ़ाइल में अंतिम तीन कॉलम नमूनों के लिए जीनोटाइप ब्लॉक हैं, जो हाइलाइट की गई हेडर लाइन में दिखाए अनुसार उनके ID फ़ील्ड के वर्णानुक्रम में सूचीबद्ध हैं।

अगर हम बहुत पहले वेरिएंट के लिए हमारे परीक्षण परिवार ट्रायो के लिए कॉल किए गए जीनोटाइप्स को देखते हैं, तो हम देखते हैं कि पिता हेटेरोजाइगस-वेरिएंट (`0/1`) है, और मां और बेटा दोनों होमोजाइगस-वेरिएंट (`1/1`) हैं।

अंततः यही वह जानकारी है जिसे हम डेटासेट से निकालना चाहते हैं!

#### 2.3.3. आउटपुट फ़ाइलों को स्थानांतरित करें

जैसा कि पहले उल्लेख किया गया है, कंटेनर के अंदर जो कुछ भी रहता है वह भविष्य के काम के लिए दुर्गम होगा।
कंटेनर से बाहर निकलने से पहले, हम GVCF फ़ाइलों, अंतिम मल्टी-सैंपल VCF और उनकी सभी इंडेक्स फ़ाइलों को कंटेनर के बाहर फ़ाइलसिस्टम में मैनुअल रूप से स्थानांतरित करने जा रहे हैं।
इस तरह, हमारे पास तुलना करने के लिए कुछ होगा जब हम इस सभी काम को स्वचालित करने के लिए अपना वर्कफ़्लो बनाएंगे।

```bash
mv *.vcf* /data/vcf
```

??? abstract "डायरेक्टरी की सामग्री" hl_lines="14-19 22-23"

    ```console
    data
    ├── bam
    │   ├── reads_father.bam
    │   ├── reads_father.bam.bai
    │   ├── reads_mother.bam
    │   ├── reads_mother.bam.bai
    │   ├── reads_son.bam
    │   └── reads_son.bam.bai
    ├── ref
    │   ├── intervals.bed
    │   ├── ref.dict
    │   ├── ref.fasta
    │   └── ref.fasta.fai
    ├── samplesheet.csv
    └── vcf
        ├── family_trio.vcf
        ├── family_trio.vcf.idx
        ├── reads_father.g.vcf
        ├── reads_father.g.vcf.idx
        ├── reads_mother.g.vcf
        ├── reads_mother.g.vcf.idx
        ├── reads_mother.vcf
        ├── reads_mother.vcf.idx
        ├── reads_son.g.vcf
        └── reads_son.g.vcf.idx
    ```

एक बार यह हो जाने के बाद, सभी फ़ाइलें अब तुम्हारी सामान्य फ़ाइलसिस्टम में सुलभ हैं।

#### 2.3.4. GATK कंटेनर से बाहर निकलें

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करो।

```bash
exit
```

तुम्हारा प्रॉम्प्ट वापस सामान्य होना चाहिए।
यह संयुक्त वेरिएंट कॉलिंग कमांड्स के मैनुअल परीक्षण को समाप्त करता है।

---

### सारांश

तुम जानते हो कि Samtools इंडेक्सिंग और GATK वेरिएंट कॉलिंग कमांड्स को उनके संबंधित कंटेनरों में कैसे परीक्षण करना है, जिसमें GVCFs जेनरेट करना और कई नमूनों पर संयुक्त जीनोटाइपिंग चलाना शामिल है।

### आगे क्या है?

एक ब्रेक लो, फिर [भाग 2](./02_per_sample_variant_calling.md) पर जाओ ताकि सीख सको कि उन्हीं कमांड्स को वर्कफ़्लो में कैसे wrap करें जो काम को निष्पादित करने के लिए कंटेनरों का उपयोग करते हैं।
