# भाग 1: प्रति-नमूना वेरिएंट कॉलिंग

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

इस कोर्स के पहले भाग में, हम आपको एक सरल वेरिएंट कॉलिंग पाइपलाइन बनाना सिखाएंगे जो व्यक्तिगत सीक्वेंसिंग नमूनों पर GATK वेरिएंट कॉलिंग लागू करती है।

### विधि का अवलोकन

वेरिएंट कॉलिंग एक जीनोमिक विश्लेषण विधि है जिसका उद्देश्य एक संदर्भ जीनोम के सापेक्ष एक जीनोम अनुक्रम में भिन्नताओं की पहचान करना है।
यहाँ हम छोटे वेरिएंट कॉल करने के लिए डिज़ाइन किए गए टूल और विधियों का उपयोग करने जा रहे हैं, _यानी_ SNPs और indels।

![GATK pipeline](img/gatk-pipeline.png)

एक पूर्ण वेरिएंट कॉलिंग पाइपलाइन में आम तौर पर कई चरण शामिल होते हैं, जिनमें संदर्भ के साथ मैपिंग (कभी-कभी जीनोम संरेखण के रूप में संदर्भित) और वेरिएंट फ़िल्टरिंग और प्राथमिकता शामिल है।
सरलता के लिए, इस कोर्स के इस भाग में हम केवल वेरिएंट कॉलिंग भाग पर ध्यान केंद्रित करने जा रहे हैं।

### डेटासेट

हम निम्नलिखित डेटा और संबंधित संसाधन प्रदान करते हैं:

- **एक संदर्भ जीनोम** जिसमें मानव क्रोमोसोम 20 (hg19/b37 से) का एक छोटा क्षेत्र और इसकी सहायक फ़ाइलें (इंडेक्स और अनुक्रम शब्दकोश) शामिल हैं।
- **तीन संपूर्ण जीनोम सीक्वेंसिंग नमूने** जो एक परिवार के तीन सदस्यों (माँ, पिता और बेटा) से संबंधित हैं, जिन्हें फ़ाइल आकार को छोटा रखने के लिए क्रोमोसोम 20 पर डेटा के एक छोटे से हिस्से तक सीमित किया गया है।
  यह Illumina शॉर्ट-रीड सीक्वेंसिंग डेटा है जो पहले ही संदर्भ जीनोम से मैप किया जा चुका है, [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) प्रारूप (Binary Alignment Map, SAM का एक संपीड़ित संस्करण, Sequence Alignment Map) में प्रदान किया गया है।
- **जीनोमिक अंतरालों की एक सूची**, यानी जीनोम पर निर्देशांक जहाँ हमारे नमूनों में वेरिएंट कॉल करने के लिए उपयुक्त डेटा है, BED प्रारूप में प्रदान किया गया है।

### वर्कफ़्लो

इस कोर्स के इस भाग में, हम एक वर्कफ़्लो विकसित करने जा रहे हैं जो निम्नलिखित करती है:

1. [Samtools](https://www.htslib.org/) का उपयोग करके प्रत्येक BAM इनपुट फ़ाइल के लिए एक इंडेक्स फ़ाइल बनाएं
2. प्रत्येक BAM इनपुट फ़ाइल पर GATK HaplotypeCaller चलाएं ताकि VCF (Variant Call Format) में प्रति-नमूना वेरिएंट कॉल उत्पन्न हो सकें

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

!!! note

    इंडेक्स फ़ाइलें बायोइन्फॉर्मेटिक्स फ़ाइल प्रारूपों की एक सामान्य विशेषता हैं; उनमें मुख्य फ़ाइल की संरचना के बारे में जानकारी होती है जो GATK जैसे टूल को पूरी फ़ाइल को पढ़े बिना डेटा के एक सबसेट तक पहुँचने की अनुमति देती है।
    यह महत्वपूर्ण है क्योंकि ये फ़ाइलें कितनी बड़ी हो सकती हैं।

---

## 0. वार्मअप: Samtools और GATK कमांड को इंटरैक्टिव रूप से टेस्ट करें

पहले हम कमांड को मैन्युअल रूप से आज़माना चाहते हैं इससे पहले कि हम उन्हें वर्कफ़्लो में लपेटने का प्रयास करें।
हमें जिन टूल की आवश्यकता है (Samtools और GATK) GitHub Codespaces वातावरण में इंस्टॉल नहीं हैं, इसलिए हम उन्हें कंटेनर के माध्यम से उपयोग करेंगे ([Hello Containers](../../hello_nextflow/05_hello_containers.md) देखें)।

!!! note

     सुनिश्चित करें कि आप `nf4-science/genomics` डायरेक्टरी में हैं ताकि जब आप `pwd` टाइप करें तो दिखाए गए पथ का अंतिम भाग `genomics` हो।

### 0.1. Samtools के साथ BAM इनपुट फ़ाइल को इंडेक्स करें

हम एक Samtools कंटेनर को पुल करने जा रहे हैं, इसे इंटरैक्टिव रूप से स्पिन करें और BAM फ़ाइलों में से एक पर `samtools index` कमांड चलाएं।

#### 0.1.1. Samtools कंटेनर को पुल करें

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Command output"

    ```console

    ```
-->

#### 0.1.2. Samtools कंटेनर को इंटरैक्टिव रूप से स्पिन करें

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Command output"

    ```console

    ```
-->

#### 0.1.3. इंडेक्सिंग कमांड चलाएं

[Samtools डॉक्यूमेंटेशन](https://www.htslib.org/doc/samtools-index.html) हमें BAM फ़ाइल को इंडेक्स करने के लिए चलाने के लिए कमांड लाइन देता है।

हमें केवल इनपुट फ़ाइल प्रदान करनी होगी; टूल स्वचालित रूप से इनपुट फ़ाइलनाम में `.bai` जोड़कर आउटपुट के लिए एक नाम उत्पन्न करेगा।

```bash
samtools index /data/bam/reads_mother.bam
```

यह तुरंत पूरा होना चाहिए, और अब आपको मूल BAM इनपुट फ़ाइल के समान डायरेक्टरी में `reads_mother.bam.bai` नामक एक फ़ाइल दिखनी चाहिए।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

#### 0.1.4. Samtools कंटेनर से बाहर निकलें

```bash
exit
```

### 0.2. GATK HaplotypeCaller के साथ वेरिएंट कॉल करें

हम एक GATK कंटेनर को पुल करने जा रहे हैं, इसे इंटरैक्टिव रूप से स्पिन करें और BAM फ़ाइल पर `gatk HaplotypeCaller` कमांड चलाएं जिसे हमने अभी इंडेक्स किया है।

#### 0.2.1. GATK कंटेनर को पुल करें

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Command output"

    ```console

    ```
-->

#### 0.2.2. GATK कंटेनर को इंटरैक्टिव रूप से स्पिन करें

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Command output"

    ```console

    ```
-->

#### 0.2.3. वेरिएंट कॉलिंग कमांड चलाएं

[GATK डॉक्यूमेंटेशन](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) हमें BAM फ़ाइल पर वेरिएंट कॉलिंग करने के लिए चलाने के लिए कमांड लाइन देता है।

हमें BAM इनपुट फ़ाइल (`-I`) के साथ-साथ संदर्भ जीनोम (`-R`), आउटपुट फ़ाइल के लिए एक नाम (`-O`) और विश्लेषण करने के लिए जीनोमिक अंतरालों की एक सूची (`-L`) प्रदान करनी होगी।

हालाँकि, हमें इंडेक्स फ़ाइल के पथ को निर्दिष्ट करने की आवश्यकता नहीं है; टूल स्वचालित रूप से इसे समान डायरेक्टरी में ढूंढेगा, स्थापित नामकरण और सह-स्थान सम्मेलन के आधार पर।
यही बात संदर्भ जीनोम की सहायक फ़ाइलों (इंडेक्स और अनुक्रम शब्दकोश फ़ाइलें, `*.fai` और `*.dict`) पर भी लागू होती है।

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

<!--
??? success "Command output"

    ```console

    ```
-->

आउटपुट फ़ाइल `reads_mother.vcf` कंटेनर में आपकी वर्किंग डायरेक्टरी के अंदर बनाई गई है, इसलिए आप इसे VS Code फ़ाइल एक्सप्लोरर में तब तक नहीं देखेंगे जब तक आप आउटपुट फ़ाइल पथ नहीं बदलते।
हालाँकि, यह एक छोटी परीक्षण फ़ाइल है, इसलिए आप इसे खोलने और सामग्री देखने के लिए `cat` कर सकते हैं।
यदि आप फ़ाइल की शुरुआत तक स्क्रॉल करते हैं, तो आपको मेटाडेटा की कई पंक्तियों से बना एक हेडर मिलेगा, जिसके बाद वेरिएंट कॉल की एक सूची होगी, प्रति पंक्ति एक।

```console title="reads_mother.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

प्रत्येक पंक्ति नमूने के सीक्वेंसिंग डेटा में पहचाने गए संभावित वेरिएंट का वर्णन करती है। VCF प्रारूप की व्याख्या करने के लिए मार्गदर्शन के लिए, [यह उपयोगी लेख](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/) देखें।

आउटपुट VCF फ़ाइल के साथ `reads_mother.vcf.idx` नामक एक इंडेक्स फ़ाइल है जो स्वचालित रूप से GATK द्वारा बनाई गई थी।
इसका वही कार्य है जो BAM इंडेक्स फ़ाइल का है, टूल को पूरी फ़ाइल लोड किए बिना डेटा के सबसेट को खोजने और पुनर्प्राप्त करने की अनुमति देना।

#### 0.2.4. GATK कंटेनर से बाहर निकलें

```bash
exit
```

### निष्कर्ष

आप जानते हैं कि Samtools इंडेक्सिंग और GATK वेरिएंट कॉलिंग कमांड को उनके संबंधित कंटेनर में कैसे टेस्ट करें।

### आगे क्या है?

सीखें कि उन्हीं कमांड को दो-चरणीय वर्कफ़्लो में कैसे लपेटें जो कार्य निष्पादित करने के लिए कंटेनर का उपयोग करती है।

---

## 1. एक सिंगल-स्टेज वर्कफ़्लो लिखें जो BAM फ़ाइल पर Samtools index चलाती है

हम आपको एक वर्कफ़्लो फ़ाइल प्रदान करते हैं, `genomics-1.nf`, जो वर्कफ़्लो के मुख्य भागों की रूपरेखा देती है।
यह कार्यात्मक नहीं है; इसका उद्देश्य केवल एक कंकाल के रूप में काम करना है जिसका उपयोग आप वास्तविक वर्कफ़्लो लिखने के लिए करेंगे।

### 1.1. इंडेक्सिंग प्रोसेस को परिभाषित करें

आइए एक प्रोसेस लिखें, जिसे हम `SAMTOOLS_INDEX` कहेंगे, जो इंडेक्सिंग ऑपरेशन का वर्णन करती है।

```groovy title="genomics-1.nf" linenums="9"
/*
 * BAM इंडेक्स फ़ाइल जेनरेट करें
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    path "${input_bam}.bai"

    script:
    """
    samtools index '$input_bam'
    """
}
```

आपको इस प्रशिक्षण श्रृंखला के भाग 1 और भाग 2 में जो कुछ सीखा है, उससे सभी टुकड़ों को पहचानना चाहिए।

यह प्रोसेस हमें `input_bam` इनपुट के माध्यम से एक फ़ाइल पथ पास करने की आवश्यकता होगी, तो आइए इसे अगले सेट करें।

### 1.2. एक इनपुट पैरामीटर घोषणा जोड़ें

फ़ाइल के शीर्ष पर, `Pipeline parameters` सेक्शन के तहत, हम `reads_bam` नामक एक CLI पैरामीटर घोषित करते हैं और इसे एक डिफ़ॉल्ट मान देते हैं।
इस तरह, हम आलसी हो सकते हैं और पाइपलाइन को लॉन्च करने के लिए कमांड टाइप करते समय इनपुट निर्दिष्ट नहीं करते (विकास उद्देश्यों के लिए)।

```groovy title="genomics-1.nf" linenums="3"
/*
 * Pipeline पैरामीटर
 */
params {
    // प्राथमिक इनपुट
    reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
}
```

अब हमारे पास एक प्रोसेस तैयार है, साथ ही इसे चलाने के लिए एक इनपुट देने के लिए एक पैरामीटर है, तो आइए उन चीजों को एक साथ जोड़ें।

!!! note

    `${projectDir}` एक बिल्ट-इन Nextflow वेरिएबल है जो उस डायरेक्टरी की ओर इशारा करता है जहाँ वर्तमान Nextflow वर्कफ़्लो स्क्रिप्ट (`genomics-1.nf`) स्थित है।

    यह वर्कफ़्लो रिपॉजिटरी में शामिल फ़ाइलों, डेटा डायरेक्टरी और अन्य संसाधनों को संदर्भित करना आसान बनाता है बिना पूर्ण पथों को हार्डकोड किए।

### 1.3. SAMTOOLS_INDEX चलाने के लिए वर्कफ़्लो ब्लॉक जोड़ें

`workflow` ब्लॉक में, हमें `SAMTOOLS_INDEX` प्रोसेस को इनपुट देने के लिए एक **channel** सेट करना होगा; फिर हम उस चैनल की सामग्री पर चलाने के लिए प्रोसेस को स्वयं कॉल कर सकते हैं।

```groovy title="genomics-1.nf" linenums="24"
workflow {

    main:
    // इनपुट चैनल बनाएं (CLI पैरामीटर के माध्यम से एकल फ़ाइल)
    reads_ch = channel.fromPath(params.reads_bam)

    // इनपुट BAM फ़ाइल के लिए इंडेक्स फ़ाइल बनाएं
    SAMTOOLS_INDEX(reads_ch)

    publish:
    bam_index = SAMTOOLS_INDEX.out
}
```

वर्कफ़्लो ब्लॉक में दो सेक्शन हैं:

- `main:` में चैनल ऑपरेशन और प्रोसेस कॉल शामिल हैं
- `publish:` घोषित करता है कि कौन से आउटपुट प्रकाशित किए जाने चाहिए, उन्हें नामित लक्ष्यों को सौंपते हुए

आप देखेंगे कि हम वही `.fromPath` चैनल फैक्ट्री का उपयोग कर रहे हैं जो हमने [Hello Channels](../../hello_nextflow/02_hello_channels.md) में उपयोग किया था।
वास्तव में, हम कुछ बहुत समान कर रहे हैं।
अंतर यह है कि हम Nextflow को बता रहे हैं कि फ़ाइल पथ को ही चैनल में इनपुट तत्व के रूप में लोड करें, इसकी सामग्री को पढ़ने के बजाय।

### 1.4. यह परिभाषित करने के लिए एक आउटपुट ब्लॉक जोड़ें कि परिणाम कहाँ प्रकाशित किए जाते हैं

वर्कफ़्लो ब्लॉक के बाद, हम एक `output` ब्लॉक जोड़ते हैं जो निर्दिष्ट करता है कि वर्कफ़्लो आउटपुट कहाँ प्रकाशित करें।

```groovy title="genomics-1.nf" linenums="37"
output {
    bam_index {
        path '.'
    }
}
```

`publish:` सेक्शन से प्रत्येक नामित लक्ष्य (जैसे `bam_index`) को अपना ब्लॉक मिलता है जहाँ आप बेस आउटपुट डायरेक्टरी के सापेक्ष आउटपुट पथ को कॉन्फ़िगर कर सकते हैं।

!!! note

    भले ही हम यहाँ जिन डेटा फ़ाइलों का उपयोग कर रहे हैं वे बहुत छोटी हैं, जीनोमिक्स में वे बहुत बड़ी हो सकती हैं।
    डिफ़ॉल्ट रूप से, Nextflow प्रकाशन डायरेक्टरी में आउटपुट फ़ाइलों के लिए सिम्बोलिक लिंक बनाता है, जो अनावश्यक फ़ाइल प्रतियों से बचता है।
    आप `mode` विकल्प का उपयोग करके इस व्यवहार को बदल सकते हैं (जैसे, `mode 'copy'`) वास्तविक प्रतियाँ बनाने के लिए।
    ध्यान रखें कि जब आप अपनी `work` डायरेक्टरी को साफ़ करते हैं तो सिमलिंक टूट जाएंगे, इसलिए प्रोडक्शन वर्कफ़्लो के लिए आप `mode 'copy'` का उपयोग करना चाह सकते हैं।

### 1.5. आउटपुट डायरेक्टरी को कॉन्फ़िगर करें

बेस आउटपुट डायरेक्टरी `outputDir` कॉन्फ़िग विकल्प के माध्यम से सेट की जाती है। इसे `nextflow.config` में जोड़ें:

=== "बाद में"

    ```groovy title="nextflow.config" hl_lines="2"
    docker.enabled = true
    outputDir = 'results_genomics'
    ```

=== "पहले"

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

### 1.6. वर्कफ़्लो चलाएं यह सत्यापित करने के लिए कि इंडेक्सिंग चरण काम करता है

चलो वर्कफ़्लो चलाएं! एक अनुस्मारक के रूप में, हमें कमांड लाइन में इनपुट निर्दिष्ट करने की आवश्यकता नहीं है क्योंकि हमने इनपुट पैरामीटर घोषित करते समय इनपुट के लिए एक डिफ़ॉल्ट मान सेट किया था।

```bash
nextflow run genomics-1.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

आप वर्क डायरेक्टरी या परिणाम डायरेक्टरी में देखकर जाँच सकते हैं कि इंडेक्स फ़ाइल सही ढंग से जेनरेट की गई है।

??? abstract "वर्क डायरेक्टरी सामग्री"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "परिणाम डायरेक्टरी सामग्री"

    ```console
    results_genomics/
    └── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ```

यह रहा!

### निष्कर्ष

आप जानते हैं कि एक जीनोमिक्स टूल को सिंगल-स्टेप Nextflow वर्कफ़्लो में कैसे लपेटें और इसे कंटेनर का उपयोग करके चलाएं।

### आगे क्या है?

एक दूसरा चरण जोड़ें जो पहले के आउटपुट का उपभोग करता है।

---

## 2. इंडेक्स की गई BAM फ़ाइल पर GATK HaplotypeCaller चलाने के लिए दूसरा प्रोसेस जोड़ें

अब जब हमारे पास हमारी इनपुट फ़ाइल के लिए एक इंडेक्स है, तो हम वेरिएंट कॉलिंग चरण को सेट करने के लिए आगे बढ़ सकते हैं, जो वर्कफ़्लो का दिलचस्प हिस्सा है।

### 2.1. वेरिएंट कॉलिंग प्रोसेस को परिभाषित करें

आइए एक प्रोसेस लिखें, जिसे हम `GATK_HAPLOTYPECALLER` कहेंगे, जो वेरिएंट कॉलिंग ऑपरेशन का वर्णन करती है।

```groovy title="genomics-1.nf" linenums="44"
/*
 * GATK HaplotypeCaller के साथ वेरिएंट कॉल करें
 */
process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    input:
    path input_bam
    path input_bam_index
    path ref_fasta
    path ref_index
    path ref_dict
    path interval_list

    output:
    path "${input_bam}.vcf"     , emit: vcf
    path "${input_bam}.vcf.idx" , emit: idx

    script:
    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.vcf \
        -L ${interval_list}
    """
}
```

आप देखेंगे कि हमने यहाँ कुछ नया सिंटैक्स पेश किया है (`emit:`) हमारे प्रत्येक आउटपुट चैनल को विशिष्ट रूप से नाम देने के लिए, और इसके कारण जल्द ही स्पष्ट हो जाएंगे।

यह कमांड काफी अधिक इनपुट लेता है, क्योंकि GATK को एक सरल इंडेक्सिंग कार्य की तुलना में विश्लेषण करने के लिए अधिक जानकारी की आवश्यकता होती है।
लेकिन आप ध्यान देंगे कि GATK कमांड में सूचीबद्ध की तुलना में इनपुट ब्लॉक में और भी अधिक इनपुट परिभाषित हैं। ऐसा क्यों है?

!!! note

    GATK BAM इंडेक्स फ़ाइल और संदर्भ जीनोम की सहायक फ़ाइलों को देखना जानता है क्योंकि यह उन फ़ाइलों के आसपास के सम्मेलनों के बारे में जानता है।
    हालाँकि, Nextflow को डोमेन-अज्ञेयवादी होने के लिए डिज़ाइन किया गया है और यह बायोइन्फॉर्मेटिक्स फ़ाइल प्रारूप आवश्यकताओं के बारे में कुछ नहीं जानता है।

हमें Nextflow को स्पष्ट रूप से बताना होगा कि उसे रनटाइम पर वर्किंग डायरेक्टरी में उन फ़ाइलों को स्टेज करना होगा; अन्यथा यह ऐसा नहीं करेगा, और GATK (सही ढंग से) इंडेक्स फ़ाइलों के गायब होने के बारे में एक त्रुटि फेंकेगा।

इसी तरह, हमें आउटपुट VCF की इंडेक्स फ़ाइल (`"${input_bam}.vcf.idx"` फ़ाइल) को स्पष्ट रूप से सूचीबद्ध करना होगा ताकि Nextflow को पता चले कि यदि बाद के चरणों में इसकी आवश्यकता हो तो उस फ़ाइल का ट्रैक रखना होगा।

### 2.2. सहायक इनपुट के लिए परिभाषाएं जोड़ें

चूंकि हमारी नई प्रोसेस को कुछ अतिरिक्त फ़ाइलें प्रदान करने की अपेक्षा है, हम `Pipeline parameters` सेक्शन के तहत उनके लिए कुछ CLI पैरामीटर सेट करते हैं, कुछ डिफ़ॉल्ट मानों के साथ (पहले जैसे ही कारण)।

```groovy title="genomics-1.nf" linenums="8"
    // सहायक फ़ाइलें
    reference: Path = "${projectDir}/data/ref/ref.fasta"
    reference_index: Path = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict: Path = "${projectDir}/data/ref/ref.dict"
    intervals: Path = "${projectDir}/data/ref/intervals.bed"
```

### 2.3. सहायक फ़ाइल पथों को रखने के लिए वेरिएबल बनाएं

जबकि मुख्य डेटा इनपुट चैनलों के माध्यम से गतिशील रूप से स्ट्रीम किए जाते हैं, सहायक फ़ाइलों को संभालने के लिए दो दृष्टिकोण हैं। अनुशंसित दृष्टिकोण स्पष्ट चैनल बनाना है, जो डेटा प्रवाह को स्पष्ट और अधिक सुसंगत बनाता है। वैकल्पिक रूप से, सरल मामलों के लिए वेरिएबल बनाने के लिए file() फ़ंक्शन का उपयोग किया जा सकता है, विशेष रूप से जब आपको कई प्रोसेस में एक ही फ़ाइल को संदर्भित करने की आवश्यकता हो - हालाँकि ध्यान रखें कि यह अभी भी चैनल बनाता है अस्पष्ट रूप से। <!-- TODO: स्पष्ट करें: क्या यह अभी भी टाइप किए गए इनपुट के साथ आवश्यक है? -->

इसे वर्कफ़्लो ब्लॉक में जोड़ें (`reads_ch` निर्माण के बाद, `main:` सेक्शन के अंदर):

```groovy title="genomics-1.nf" linenums="79"
    // सहायक फ़ाइलों (संदर्भ और अंतराल) के लिए फ़ाइल पथ लोड करें
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)
```

यह सहायक फ़ाइल पथों को किसी भी प्रोसेस को इनपुट के रूप में प्रदान करने के लिए उपलब्ध कराएगा जिन्हें उनकी आवश्यकता है।

### 2.4. GATK_HAPLOTYPECALLER चलाने के लिए वर्कफ़्लो ब्लॉक में एक कॉल जोड़ें

अब जब हमने अपनी दूसरी प्रोसेस सेट कर ली है और सभी इनपुट और सहायक फ़ाइलें तैयार और उपलब्ध हैं, तो हम वर्कफ़्लो बॉडी में `GATK_HAPLOTYPECALLER` प्रोसेस में एक कॉल जोड़ सकते हैं।

```groovy title="genomics-1.nf" linenums="88"
    // इंडेक्स की गई BAM फ़ाइल से वेरिएंट कॉल करें
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )
```

आपको इस प्रशिक्षण श्रृंखला के भाग 1 से `*.out` सिंटैक्स को पहचानना चाहिए; हम Nextflow को बता रहे हैं कि `SAMTOOLS_INDEX` द्वारा आउटपुट किए गए चैनल को लें और उसे `GATK_HAPLOTYPECALLER` प्रोसेस कॉल में प्लग करें।

!!! note

    आप देखेंगे कि इनपुट प्रोसेस की कॉल में उसी क्रम में प्रदान किए जाते हैं जैसे वे प्रोसेस के इनपुट ब्लॉक में सूचीबद्ध हैं।
    Nextflow में, इनपुट स्थितीय हैं, जिसका अर्थ है कि आपको _उसी क्रम का पालन करना होगा_; और बेशक तत्वों की संख्या समान होनी चाहिए।

### 2.5. प्रकाशन सेक्शन और आउटपुट ब्लॉक को अपडेट करें

हमें VCF आउटपुट को शामिल करने के लिए `publish:` सेक्शन को अपडेट करने की आवश्यकता है, और `output` ब्लॉक में संबंधित लक्ष्य जोड़ने की आवश्यकता है।

```groovy title="genomics-1.nf" linenums="99"
    publish:
    bam_index = SAMTOOLS_INDEX.out
    vcf = GATK_HAPLOTYPECALLER.out.vcf
    vcf_idx = GATK_HAPLOTYPECALLER.out.idx
}

output {
    bam_index {
        path '.'
    }
    vcf {
        path '.'
    }
    vcf_idx {
        path '.'
    }
}
```

### 2.6. वर्कफ़्लो चलाएं यह सत्यापित करने के लिए कि वेरिएंट कॉलिंग चरण काम करता है

चलो विस्तारित वर्कफ़्लो को `-resume` के साथ चलाएं ताकि हमें इंडेक्सिंग चरण को फिर से चलाने की आवश्यकता न हो।

```bash
nextflow run genomics-1.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

अब यदि हम कंसोल आउटपुट को देखते हैं, तो हम दो प्रोसेस सूचीबद्ध देखते हैं।

पहली प्रोसेस को कैशिंग के कारण छोड़ दिया गया था, जैसा कि अपेक्षित था, जबकि दूसरी प्रोसेस चलाई गई क्योंकि यह बिल्कुल नई है।

आप परिणाम डायरेक्टरी में आउटपुट फ़ाइलें पाएंगे (वर्क डायरेक्टरी के सिम्बोलिक लिंक के रूप में)।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results_genomics/
    ├── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */cf/36f756*/reads_mother.bam.vcf
    └── reads_mother.bam.vcf.idx -> */cf/36f756*/reads_mother.bam.vcf.idx
    ```

यदि आप VCF फ़ाइल खोलते हैं, तो आपको उसी सामग्री को देखना चाहिए जो आपने GATK कमांड को सीधे कंटेनर में चलाकर जेनरेट की गई फ़ाइल में की थी।

```console title="reads_mother.bam.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

यह वह आउटपुट है जिसे हम अपने अध्ययन में प्रत्येक नमूने के लिए उत्पन्न करने की परवाह करते हैं।

### निष्कर्ष

आप जानते हैं कि एक बहुत ही बुनियादी दो-चरणीय वर्कफ़्लो कैसे बनाएं जो वास्तविक विश्लेषण कार्य करती है और जीनोमिक्स फ़ाइल प्रारूप की विशिष्टताओं जैसे सहायक फ़ाइलों से निपटने में सक्षम है।

### आगे क्या है?

वर्कफ़्लो को नमूनों के एक बैच को थोक में संभालने के लिए अनुकूलित करें।

---

## 3. वर्कफ़्लो को नमूनों के बैच पर चलाने के लिए अनुकूलित करें

यह सब ठीक है कि एक वर्कफ़्लो है जो एक एकल नमूने पर प्रसंस्करण को स्वचालित कर सकती है, लेकिन यदि आपके पास 1000 नमूने हैं तो क्या होगा?
क्या आपको एक bash स्क्रिप्ट लिखने की आवश्यकता है जो आपके सभी नमूनों के माध्यम से लूप करती है?

नहीं, भगवान का शुक्र है! बस कोड में एक मामूली बदलाव करें और Nextflow आपके लिए भी इसे संभाल लेगा।

### 3.1. इनपुट पैरामीटर घोषणा को तीन नमूनों को सूचीबद्ध करने वाली एक सरणी में बदलें

आइए इनपुट BAM फ़ाइल घोषणा में उस डिफ़ॉल्ट फ़ाइल पथ को हमारे तीन परीक्षण नमूनों के लिए फ़ाइल पथों को सूचीबद्ध करने वाली एक सरणी में बदलें, `Pipeline parameters` सेक्शन के तहत।

=== "बाद में"

    ```groovy title="genomics-1.nf" linenums="7"
    // प्राथमिक इनपुट (तीन नमूनों की सरणी)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

=== "पहले"

    ```groovy title="genomics-1.nf" linenums="7"
        // प्राथमिक इनपुट
        reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
    ```

!!! note

    टाइप किए गए पैरामीटर घोषणाओं (जैसे `reads_bam: Path`) का उपयोग करते समय, आप एक सरणी मान नहीं सौंप सकते।
    सरणियों के लिए, प्रकार एनोटेशन को छोड़ दें।

और वास्तव में यह सब हमें करने की आवश्यकता है, क्योंकि चैनल फैक्ट्री जो हम वर्कफ़्लो बॉडी में उपयोग करते हैं (`.fromPath`) इनपुट चैनल में लोड करने के लिए कई फ़ाइल पथों को स्वीकार करने में उतना ही खुश है जितना यह एक को लोड करने में था।

!!! note

    सामान्य रूप से, आप नमूनों की सूची को अपनी वर्कफ़्लो फ़ाइल में हार्डकोड नहीं करना चाहेंगे, लेकिन हम यहां चीजों को सरल रखने के लिए ऐसा कर रहे हैं।
    हम इस प्रशिक्षण श्रृंखला में बाद में इनपुट को संभालने के लिए अधिक सुरुचिपूर्ण तरीके प्रस्तुत करेंगे।

### 3.2. वर्कफ़्लो चलाएं यह सत्यापित करने के लिए कि यह सभी तीन नमूनों पर चलती है

आइए अब वर्कफ़्लो चलाने का प्रयास करें जब प्लंबिंग सभी तीन परीक्षण नमूनों पर चलाने के लिए सेट की गई है।

```bash
nextflow run genomics-1.nf -resume
```

मज़ेदार बात: यह _काम कर सकता है_, OR यह _विफल हो सकता है_। उदाहरण के लिए, यहाँ एक रन है जो सफल हुआ:

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

यदि आपका वर्कफ़्लो रन सफल रहा, तो इसे तब तक फिर से चलाएं जब तक आपको इस तरह की त्रुटि न मिले:

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

    executor >  local (4)
    [01/eea165] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 1 ✔
    [a5/fa9fd0] GATK_HAPLOTYPECALLER (3) | 1 of 3, cached: 1
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Process `GATK_HAPLOTYPECALLER (2)` terminated with an error exit status (2)

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_father.bam         -O reads_father.bam.vcf         -L intervals.bed

    Command exit status:
      2

    Command error:
      ...
      A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
      ...
    ```

यदि आप GATK कमांड त्रुटि आउटपुट को देखते हैं, तो इस तरह की एक पंक्ति होगी:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

खैर, यह अजीब है, यह देखते हुए कि हमने वर्कफ़्लो के पहले चरण में BAM फ़ाइलों को स्पष्ट रूप से इंडेक्स किया था। क्या प्लंबिंग में कुछ गलत हो सकता है?

#### 3.2.1. प्रासंगिक कॉल के लिए वर्क डायरेक्टरी की जाँच करें

आइए कंसोल आउटपुट में सूचीबद्ध विफल `GATK_HAPLOTYPECALLER` प्रोसेस कॉल के लिए वर्क डायरेक्टरी के अंदर देखें।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    work/a5/fa9fd0994b6beede5fb9ea073596c2
    ├── intervals.bed -> /workspaces/training/nf4-science/genomics/data/ref/intervals.bed
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/01/eea16597bd6e810fb4cf89e60f8c2d/reads_father.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    ├── reads_son.bam.vcf
    ├── reads_son.bam.vcf.idx
    ├── ref.dict -> /workspaces/training/nf4-science/genomics/data/ref/ref.dict
    ├── ref.fasta -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta
    └── ref.fasta.fai -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta.fai
    ```

इस डायरेक्टरी में सूचीबद्ध BAM फ़ाइल और BAM इंडेक्स के नामों पर विशेष ध्यान दें: `reads_son.bam` और `reads_father.bam.bai`।

यह क्या है? Nextflow ने इस प्रोसेस कॉल की वर्क डायरेक्टरी में एक इंडेक्स फ़ाइल स्टेज की है, लेकिन यह गलत है। यह कैसे हुआ?

#### 3.2.2. चैनल सामग्री का निरीक्षण करने के लिए [view() ऑपरेटर](https://www.nextflow.io/docs/latest/reference/operator.html#view) का उपयोग करें

`GATK_HAPLOTYPER` प्रोसेस कॉल से पहले वर्कफ़्लो बॉडी में ये दो पंक्तियाँ जोड़ें:

```groovy title="genomics-1.nf" linenums="84"
    // अस्थायी निदान
    reads_ch.view()
    SAMTOOLS_INDEX.out.view()
```

फिर वर्कफ़्लो कमांड को फिर से चलाएं।

```bash
nextflow run genomics-1.nf
```

एक बार फिर, यह सफल हो सकता है या विफल हो सकता है। यहाँ एक सफल रन है:

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [fervent_pasteur] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    [a2/dbd8d5] GATK_HAPLOTYPECALLER (3) | 3 of 3 ✔
    ```

और यहाँ एक विफल है:

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2
