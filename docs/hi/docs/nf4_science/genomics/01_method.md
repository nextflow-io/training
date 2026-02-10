# भाग 1: विधि का अवलोकन और मैनुअल परीक्षण

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता से अनुवाद किया गया - [और जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

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

ये टूल्स GitHub Codespaces वातावरण में इंस्टॉल नहीं हैं, इसलिए हम उन्हें कंटेनरों के माध्यम से उपयोग करेंगे ([Hello Containers](../../hello_nextflow/05_hello_containers.md) देखें)।

!!! note

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

तुम्हारा प्रॉम्प्ट `(base) root@a1b2c3d4e5f6:/tmp#` जैसे कुछ में बदल जाता है, जो दर्शाता है कि तुम अब कंटेनर के अंदर हो।
डेटा फ़ाइलें `/data` के तहत सुलभ हैं।

#### 1.1.3. इंडेक्सिंग कमांड चलाएं

[Samtools डॉक्यूमेंटेशन](https://www.htslib.org/doc/samtools-index.html) हमें BAM फ़ाइल को इंडेक्स करने के लिए चलाने वाली कमांड लाइन देता है।

हमें केवल इनपुट फ़ाइल प्रदान करने की आवश्यकता है; टूल स्वचालित रूप से इनपुट फ़ाइलनाम में `.bai` जोड़कर आउटपुट के लिए एक नाम जेनरेट करेगा।

```bash
samtools index /data/bam/reads_mother.bam
```

??? abstract "डायरेक्टरी की सामग्री"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

तुम्हें अब मूल BAM इनपुट फ़ाइल के समान डायरेक्टरी में `reads_mother.bam.bai` नामक एक फ़ाइल दिखनी चाहिए।

#### 1.1.4. Samtools कंटेनर से बाहर निकलें

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करो।

```bash
exit
```

तुम्हारा प्रॉम्प्ट अब वापस वैसा होना चाहिए जैसा कंटेनर शुरू करने से पहले था।

### 1.2. GATK HaplotypeCaller के साथ वेरिएंट कॉल करें

हम एक GATK कंटेनर को पुल करने जा रहे हैं, इसे इंटरएक्टिव रूप से स्पिन अप करेंगे और BAM फ़ाइल पर `gatk HaplotypeCaller` कमांड चलाएंगे जिसे हमने अभी इंडेक्स किया है।

#### 1.2.1. GATK कंटेनर को पुल करें

GATK कंटेनर इमेज डाउनलोड करने के लिए `docker pull` कमांड चलाएं:

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
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

??? success "कमांड आउटपुट"

    टूल verbose लॉगिंग आउटपुट उत्पन्न करता है। हाइलाइट की गई लाइनें सफल पूर्णता की पुष्टि करती हैं।

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

आउटपुट फ़ाइल `reads_mother.vcf` कंटेनर में तुम्हारी कार्यशील डायरेक्टरी के अंदर बनाई गई है, इसलिए तुम्हें यह VS Code फ़ाइल एक्सप्लोरर में तब तक नहीं दिखेगी जब तक तुम आउटपुट फ़ाइल पथ नहीं बदलते।
हालांकि, यह एक छोटी परीक्षण फ़ाइल है, इसलिए तुम इसे खोलने और सामग्री देखने के लिए `cat` कर सकते हो।
अगर तुम फ़ाइल की शुरुआत तक स्क्रॉल करते हो, तो तुम्हें मेटाडेटा की कई लाइनों से बना एक हेडर मिलेगा, उसके बाद वेरिएंट कॉल की एक सूची होगी, प्रति लाइन एक।

??? abstract "फ़ाइल की सामग्री"

    ```console title="reads_mother.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

प्रत्येक लाइन नमूने के सीक्वेंसिंग डेटा में पहचाने गए संभावित वेरिएंट का वर्णन करती है। VCF फॉर्मेट की व्याख्या करने के लिए मार्गदर्शन के लिए, [यह उपयोगी लेख](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/) देखें।

आउटपुट VCF फ़ाइल के साथ `reads_mother.vcf.idx` नामक एक इंडेक्स फ़ाइल है जो स्वचालित रूप से GATK द्वारा बनाई गई थी।
इसका वही कार्य है जो BAM इंडेक्स फ़ाइल का है, टूल्स को संपूर्ण फ़ाइल लोड किए बिना डेटा के सबसेट को खोजने और पुनः प्राप्त करने की अनुमति देना।

#### 1.2.4. GATK कंटेनर से बाहर निकलें

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करो।

```bash
exit
```

तुम्हारा प्रॉम्प्ट वापस सामान्य होना चाहिए।
यह प्रति-नमूना वेरिएंट कॉलिंग परीक्षण समाप्त करता है।

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
GVCFs को जोड़ना GATK GenomicsDBImport के साथ किया जाता है, जो प्रति-नमूना कॉल को डेटा स्टोर (डेटाबेस के समान) में जोड़ता है।
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

??? abstract "डायरेक्टरी की सामग्री"

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
    00:28:03.593 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:03.765 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.768 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:03.768 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:03.768 INFO  HaplotypeCaller - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:03.768 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:03.769 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:28:03 AM GMT
    00:28:03.769 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.770 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.772 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:28:03.773 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:28:03.773 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:03.774 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:03.774 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:28:03.774 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:28:03.775 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:28:03.775 INFO  HaplotypeCaller - Requester pays: disabled
    00:28:03.776 INFO  HaplotypeCaller - Initializing engine
    00:28:03.896 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:03.919 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:03.934 INFO  HaplotypeCaller - Done initializing engine
    00:28:03.935 INFO  HaplotypeCallerEngine - Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled
    00:28:03.943 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:28:03.945 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:28:03.946 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:28:03.955 INFO  HaplotypeCallerEngine - Standard Emitting and Calling confidence set to -0.0 for reference-model confidence output
    00:28:03.956 INFO  HaplotypeCallerEngine - All sites annotated with PLs forced to true for reference-model confidence output
    00:28:03.972 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:28:03.993 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:28:03.994 INFO  IntelPairHmm - Available threads: 10
    00:28:03.994 INFO  IntelPairHmm - Requested threads: 4
    00:28:03.994 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:28:04.044 INFO  ProgressMeter - Starting traversal
    00:28:04.070 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:28:04.874 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:28:06.535 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:28:06.537 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35            851.6
    00:28:06.538 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:28:06.543 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003648749
    00:28:06.544 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.031498916
    00:28:06.544 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:28:06.547 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:28:06 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.05 minutes.
    Runtime.totalMemory()=281018368
    ```

यह कंटेनर में वर्तमान कार्यशील डायरेक्टरी में GVCF आउटपुट फ़ाइल `reads_mother.g.vcf` बनाता है।

अगर तुम सामग्री देखने के लिए इसे `cat` करते हो, तो तुम्हें दिखेगा कि यह सेक्शन 1 में जेनरेट किए गए समकक्ष VCF की तुलना में बहुत लंबा है। तुम फ़ाइल की शुरुआत तक स्क्रॉल भी नहीं कर सकते, और अधिकांश लाइनें VCF में देखी गई चीज़ों से काफी अलग दिखती हैं।

??? abstract "फ़ाइल की सामग्री"

    ```console title="reads_mother.g.vcf" linenums="1674"
    20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
    20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
    20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
    ```

ये गैर-वेरिएंट क्षेत्रों का प्रतिनिधित्व करते हैं जहां वेरिएंट कॉलर को भिन्नता का कोई प्रमाण नहीं मिला, इसलिए इसने कुछ सांख्यिकी को कैप्चर किया जो भिन्नता की अनुपस्थिति में इसके विश्वास के स्तर का वर्णन करती हैं।
यह दो बहुत अलग मामलों के बीच अंतर करना संभव बनाता है: (1) अच्छी गुणवत्ता वाला डेटा है जो दिखाता है कि नमूना होमोजाइगस-रेफरेंस है, और (2) किसी भी तरह से निर्धारण करने के लिए पर्याप्त अच्छा डेटा उपलब्ध नहीं है।

एक GVCF में, आमतौर पर ऐसी कई गैर-वेरिएंट लाइनें होती हैं, जिनके बीच कम संख्या में वेरिएंट रिकॉर्ड बिखरे होते हैं।
एक वास्तविक वेरिएंट कॉल खोजने के लिए GVCF पर `head -176` चलाकर फ़ाइल की पहली 176 लाइनें लोड करने का प्रयास करो।

??? abstract "फ़ाइल की सामग्री"

    ```console title="reads_mother.g.vcf" linenums="174"
    20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
    20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
    20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
    ```

दूसरी लाइन फ़ाइल में पहला वेरिएंट रिकॉर्ड दिखाती है, जो VCF फ़ाइल में पहले वेरिएंट से मेल खाती है जिसे हमने पहले देखा था।

मूल VCF की तरह ही, आउटपुट GVCF फ़ाइल भी एक इंडेक्स फ़ाइल के साथ है, जिसे `reads_mother.g.vcf.idx` कहा जाता है।

#### 2.2.3. अन्य दो नमूनों पर प्रक्रिया को दोहराएं

बाकी दो नमूनों के लिए GVCFs जेनरेट करने के लिए नीचे दी गई कमांड्स चलाएं, एक के बाद एक।

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

एक बार यह पूरा हो जाने के बाद, तुम्हारे पास अपनी वर्तमान डायरेक्टरी में `.g.vcf` में समाप्त होने वाली तीन फ़ाइलें (प्रति नमूना एक) और उनकी संबंधित इंडेक्स फ़ाइलें `.g.vcf.idx` में समाप्त होनी चाहिए।

लेकिन कंटेनर से बाहर मत निकलो!
हम अगले चरण में उसी कंटेनर का उपयोग करने जा रहे हैं।

### 2.3. संयुक्त जीनोटाइपिंग चलाएं

अब जब हमारे पास सभी GVCFs हैं, तो हम नमूनों के कोहॉर्ट के लिए वेरिएंट कॉल जेनरेट करने के लिए संयुक्त जीनोटाइपिंग दृष्टिकोण को आज़मा सकते हैं।
यह एक दो-चरणीय विधि है जिसमें सभी GVCFs से डेटा को डेटा स्टोर में जोड़ना शामिल है, फिर संयुक्त-कॉल किए गए वेरिएंट के अंतिम VCF को जेनरेट करने के लिए वास्तविक संयुक्त जीनोटाइपिंग विश्लेषण चलाना।

#### 2.3.1. सभी प्रति-नमूना GVCFs को जोड़ें

यह पहला चरण एक अन्य GATK टूल का उपयोग करता है, जिसे GenomicsDBImport कहा जाता है, सभी GVCFs से डेटा को GenomicsDB डेटा स्टोर में जोड़ने के लिए।

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
    00:28:20.772 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:20.914 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.917 INFO  GenomicsDBImport - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:20.917 INFO  GenomicsDBImport - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:20.917 INFO  GenomicsDBImport - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:20.917 INFO  GenomicsDBImport - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:20.918 INFO  GenomicsDBImport - Start Date/Time: February 8, 2026 at 12:28:20 AM GMT
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.920 INFO  GenomicsDBImport - HTSJDK Version: 4.1.0
    00:28:20.921 INFO  GenomicsDBImport - Picard Version: 3.1.1
    00:28:20.921 INFO  GenomicsDBImport - Built for Spark Version: 3.5.0
    00:28:20.922 INFO  GenomicsDBImport - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:20.923 INFO  GenomicsDBImport - Deflater: IntelDeflater
    00:28:20.924 INFO  GenomicsDBImport - Inflater: IntelInflater
    00:28:20.924 INFO  GenomicsDBImport - GCS max retries/reopens: 20
    00:28:20.924 INFO  GenomicsDBImport - Requester pays: disabled
    00:28:20.925 INFO  GenomicsDBImport - Initializing engine
    00:28:21.144 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:21.152 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:21.157 INFO  GenomicsDBImport - Done initializing engine
    00:28:21.287 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:21.290 INFO  GenomicsDBImport - Vid Map JSON file will be written to /tmp/family_trio_gdb/vidmap.json
    00:28:21.290 INFO  GenomicsDBImport - Callset Map JSON file will be written to /tmp/family_trio_gdb/callset.json
    00:28:21.291 INFO  GenomicsDBImport - Complete VCF Header will be written to /tmp/family_trio_gdb/vcfheader.vcf
    00:28:21.291 INFO  GenomicsDBImport - Importing to workspace - /tmp/family_trio_gdb
    00:28:21.453 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.757 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.859 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.979 INFO  GenomicsDBImport - Done importing batch 1/1
    00:28:21.988 INFO  GenomicsDBImport - Import completed!
    00:28:21.988 INFO  GenomicsDBImport - Shutting down engine
    [February 8, 2026 at 12:28:21 AM GMT] org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport done. Elapsed time: 0.02 minutes.
    Runtime.totalMemory()=305135616
    ```

इस चरण का आउटपुट प्रभावी रूप से एक डायरेक्टरी है जिसमें कई अलग-अलग फ़ाइलों के रूप में संयुक्त वेरिएंट डेटा रखने वाली आगे नेस्टेड डायरेक्टरियों का एक सेट होता है।
तुम इसके चारों ओर घूम सकते हो लेकिन तुम जल्दी ही देखोगे कि यह डेटा स्टोर फॉर्मेट मनुष्यों द्वारा सीधे पढ़े जाने के लिए नहीं है।

!!! note

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
    00:28:24.625 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:24.798 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.801 INFO  GenotypeGVCFs - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:24.801 INFO  GenotypeGVCFs - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:24.801 INFO  GenotypeGVCFs - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:24.801 INFO  GenotypeGVCFs - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:24.802 INFO  GenotypeGVCFs - Start Date/Time: February 8, 2026 at 12:28:24 AM GMT
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.804 INFO  GenotypeGVCFs - HTSJDK Version: 4.1.0
    00:28:24.804 INFO  GenotypeGVCFs - Picard Version: 3.1.1
    00:28:24.804 INFO  GenotypeGVCFs - Built for Spark Version: 3.5.0
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:24.806 INFO  GenotypeGVCFs - Deflater: IntelDeflater
    00:28:24.806 INFO  GenotypeGVCFs - Inflater: IntelInflater
    00:28:24.807 INFO  GenotypeGVCFs - GCS max retries/reopens: 20
    00:28:24.807 INFO  GenotypeGVCFs - Requester pays: disabled
    00:28:24.808 INFO  GenotypeGVCFs - Initializing engine
    00:28:25.023 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:25.081 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field InbreedingCoeff  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAC  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAF  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.109 INFO  GenotypeGVCFs - Done initializing engine
    00:28:25.184 INFO  ProgressMeter - Starting traversal
    00:28:25.187 INFO  ProgressMeter -        Current Locus  Elapsed Minutes    Variants Processed  Variants/Minute
    00:28:25.446 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    GENOMICSDB_TIMER,GenomicsDB iterator next() timer,Wall-clock time(s),0.15034835899999904,Cpu time(s),0.1355218420000006
    00:28:26.189 INFO  ProgressMeter - 20_10037292_10066351:13953              0.0                  3390         202994.0
    00:28:26.190 INFO  ProgressMeter - Traversal complete. Processed 3390 total variants in 0.0 minutes.
    00:28:26.194 INFO  GenotypeGVCFs - Shutting down engine
    [February 8, 2026 at 12:28:26 AM GMT] org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs done. Elapsed time: 0.03 minutes.
    Runtime.totalMemory()=296747008
    ```

यह कंटेनर में वर्तमान कार्यशील डायरेक्टरी में VCF आउटपुट फ़ाइल `family_trio.vcf` बनाता है।
यह एक और उचित रूप से छोटी फ़ाइल है इसलिए तुम इसकी सामग्री देखने के लिए इस फ़ाइल को `cat` कर सकते हो, और पहली कुछ वेरिएंट लाइनों को खोजने के लिए ऊपर स्क्रॉल कर सकते हो।

??? abstract "फ़ाइल की सामग्री"

    ```console title="family_trio.vcf" linenums="40"
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
    20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
    20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
    20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
    ```

यह पहले जेनरेट किए गए VCF के समान दिखता है, सिवाय इस बार हमारे पास सभी तीन नमूनों के लिए जीनोटाइप-स्तरीय जानकारी है।
फ़ाइल में अंतिम तीन कॉलम नमूनों के लिए जीनोटाइप ब्लॉक हैं, जो वर्णानुक्रम में सूचीबद्ध हैं।

अगर हम बहुत पहले वेरिएंट के लिए हमारे परीक्षण परिवार ट्रायो के लिए कॉल किए गए जीनोटाइप्स को देखते हैं, तो हम देखते हैं कि पिता हेटेरोजाइगस-वेरिएंट (`0/1`) है, और मां और बेटा दोनों होमोजाइगस-वेरिएंट (`1/1`) हैं।

अंततः यही वह जानकारी है जिसे हम डेटासेट से निकालना चाहते हैं!

#### 2.3.3. GATK कंटेनर से बाहर निकलें

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करो।

```bash
exit
```

तुम्हारा प्रॉम्प्ट वापस सामान्य होना चाहिए।
यह वेरिएंट कॉलिंग कमांड्स के मैनुअल परीक्षण को समाप्त करता है।

---

### सारांश

तुम जानते हो कि Samtools इंडेक्सिंग और GATK वेरिएंट कॉलिंग कमांड्स को उनके संबंधित कंटेनरों में कैसे परीक्षण करना है, जिसमें GVCFs जेनरेट करना और कई नमूनों पर संयुक्त जीनोटाइपिंग चलाना शामिल है।

### आगे क्या है?

सीखो कि उन्हीं कमांड्स को वर्कफ़्लो में कैसे wrap करें जो काम को निष्पादित करने के लिए कंटेनरों का उपयोग करते हैं।
