# भाग 2: प्रति-नमूना वेरिएंट कॉलिंग

भाग 1 में, तुमने Samtools और GATK कमांड्स को उनके संबंधित कंटेनर्स में मैन्युअली टेस्ट किया।
अब हम उन्हीं कमांड्स को एक Nextflow वर्कफ़्लो में रैप करने जा रहे हैं।

## असाइनमेंट

इस कोर्स के इस भाग में, हम एक वर्कफ़्लो विकसित करने जा रहे हैं जो निम्नलिखित करता है:

1. [Samtools](https://www.htslib.org/) का उपयोग करके प्रत्येक BAM इनपुट फ़ाइल के लिए एक इंडेक्स फ़ाइल जेनरेट करना
2. प्रत्येक BAM इनपुट फ़ाइल पर GATK HaplotypeCaller चलाना ताकि VCF (Variant Call Format) में प्रति-नमूना वेरिएंट कॉल्स जेनरेट हो सकें

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

यह भाग 1 के स्टेप्स को दोहराता है, जहाँ तुमने इन कमांड्स को उनके कंटेनर्स में मैन्युअली चलाया था।

शुरुआती बिंदु के रूप में, हम तुम्हें एक वर्कफ़्लो फ़ाइल, `genomics.nf`, प्रदान करते हैं, जो वर्कफ़्लो के मुख्य भागों की रूपरेखा देती है, साथ ही दो मॉड्यूल फ़ाइलें, samtools_index.nf और gatk_haplotypecaller.nf, जो मॉड्यूल्स की संरचना की रूपरेखा देती हैं।
ये फ़ाइलें फंक्शनल नहीं हैं; इनका उद्देश्य सिर्फ स्कैफोल्ड के रूप में काम करना है जिन्हें तुम कोड के दिलचस्प भागों से भर सकते हो।

## पाठ योजना

विकास प्रक्रिया को अधिक शैक्षिक बनाने के लिए, हमने इसे चार चरणों में विभाजित किया है:

1. **एक सिंगल-स्टेज वर्कफ़्लो लिखना जो BAM फ़ाइल पर Samtools index चलाता है।**
   यह एक मॉड्यूल बनाना, उसे इम्पोर्ट करना, और वर्कफ़्लो में उसे कॉल करना कवर करता है।
2. **GATK HaplotypeCaller को इंडेक्स्ड BAM फ़ाइल पर चलाने के लिए दूसरा प्रोसेस जोड़ना।**
   यह प्रोसेस आउटपुट्स को इनपुट्स से चेन करना और सहायक फ़ाइलों को हैंडल करना पेश करता है।
3. **वर्कफ़्लो को नमूनों के बैच पर चलाने के लिए अनुकूलित करना।**
   यह समानांतर निष्पादन को कवर करता है और संबंधित फ़ाइलों को एक साथ रखने के लिए टपल्स पेश करता है।
4. **वर्कफ़्लो को एक टेक्स्ट फ़ाइल स्वीकार करने योग्य बनाना जिसमें इनपुट फ़ाइलों का बैच हो।**
   यह बल्क में इनपुट प्रदान करने के लिए एक सामान्य पैटर्न प्रदर्शित करता है।

प्रत्येक चरण वर्कफ़्लो विकास के एक विशिष्ट पहलू पर केंद्रित है।

---

## 1. एक सिंगल-स्टेज वर्कफ़्लो लिखना जो BAM फ़ाइल पर Samtools index चलाता है

यह पहला चरण मूल बातों पर केंद्रित है: एक BAM फ़ाइल लोड करना और उसके लिए एक इंडेक्स जेनरेट करना।

[भाग 1](01_method.md) से `samtools index` कमांड को याद करो:

```bash
samtools index '<input_bam>'
```

यह कमांड इनपुट के रूप में एक BAM फ़ाइल लेता है और उसके साथ एक `.bai` इंडेक्स फ़ाइल प्रोड्यूस करता है।
कंटेनर URI था `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464`।

हम इस जानकारी को लेने जा रहे हैं और इसे तीन चरणों में Nextflow में रैप करने जा रहे हैं:

1. इनपुट सेट अप करना
2. इंडेक्सिंग प्रोसेस लिखना और वर्कफ़्लो में उसे कॉल करना
3. आउटपुट हैंडलिंग कॉन्फ़िगर करना

### 1.1. इनपुट सेट अप करना

हमें एक इनपुट पैरामीटर डिक्लेयर करना होगा, एक सुविधाजनक डिफ़ॉल्ट वैल्यू प्रदान करने के लिए एक टेस्ट प्रोफ़ाइल बनाना होगा, और एक इनपुट चैनल बनाना होगा।

#### 1.1.1. इनपुट पैरामीटर डिक्लेरेशन जोड़ना

मुख्य वर्कफ़्लो फ़ाइल `genomics.nf` में, `Pipeline parameters` सेक्शन के अंतर्गत, `reads_bam` नाम का एक CLI पैरामीटर डिक्लेयर करो।

=== "बाद में"

    ```groovy title="genomics.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        reads_bam: Path
    }
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

यह CLI पैरामीटर सेट अप करता है, लेकिन हम विकास के दौरान हर बार वर्कफ़्लो चलाते समय फ़ाइल पाथ टाइप नहीं करना चाहते।
डिफ़ॉल्ट वैल्यू प्रदान करने के लिए कई विकल्प हैं; यहाँ हम एक टेस्ट प्रोफ़ाइल का उपयोग करते हैं।

#### 1.1.2. `nextflow.config` में डिफ़ॉल्ट वैल्यू के साथ टेस्ट प्रोफ़ाइल बनाना

एक टेस्ट प्रोफ़ाइल कमांड लाइन पर इनपुट निर्दिष्ट किए बिना वर्कफ़्लो को आज़माने के लिए सुविधाजनक डिफ़ॉल्ट वैल्यूज़ प्रदान करती है।
यह Nextflow इकोसिस्टम में एक सामान्य परंपरा है (अधिक विवरण के लिए [Hello Config](../../hello_nextflow/06_hello_config.md) देखो)।

`nextflow.config` में एक `profiles` ब्लॉक जोड़ो जिसमें एक `test` प्रोफ़ाइल हो जो `reads_bam` पैरामीटर को टेस्ट BAM फ़ाइलों में से एक पर सेट करे।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

यहाँ, हम `${projectDir}` का उपयोग कर रहे हैं, एक बिल्ट-इन Nextflow वेरिएबल जो उस डायरेक्टरी की ओर इशारा करता है जहाँ वर्कफ़्लो स्क्रिप्ट स्थित है।
यह एब्सोल्यूट पाथ्स को हार्डकोड किए बिना डेटा फ़ाइलों और अन्य संसाधनों को संदर्भित करना आसान बनाता है।

#### 1.1.3. इनपुट चैनल सेट अप करना

वर्कफ़्लो ब्लॉक में, `.fromPath` चैनल फ़ैक्टरी का उपयोग करके पैरामीटर वैल्यू से एक इनपुट चैनल बनाओ ([Hello Channels](../../hello_nextflow/02_hello_channels.md) में उपयोग किए गए अनुसार)।

=== "बाद में"

    ```groovy title="genomics.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="13"
    workflow {

        main:
        // Create input channel
    ```

अब हमें इस इनपुट पर इंडेक्सिंग चलाने के लिए प्रोसेस बनाना होगा।

### 1.2. इंडेक्सिंग प्रोसेस लिखना और वर्कफ़्लो में उसे कॉल करना

हमें मॉड्यूल फ़ाइल में प्रोसेस डेफिनिशन लिखनी होगी, include स्टेटमेंट का उपयोग करके इसे वर्कफ़्लो में इम्पोर्ट करना होगा, और इनपुट पर इसे कॉल करना होगा।

#### 1.2.1. इंडेक्सिंग प्रोसेस के लिए मॉड्यूल भरना

`modules/samtools_index.nf` खोलो और प्रोसेस डेफिनिशन की रूपरेखा की जांच करो।
तुम्हें मुख्य संरचनात्मक तत्वों को पहचानना चाहिए; यदि नहीं, तो रिफ्रेशर के लिए [Hello Nextflow](../../hello_nextflow/01_hello_world.md) पढ़ने पर विचार करो।

आगे बढ़ो और ऊपर दी गई जानकारी का उपयोग करके प्रोसेस डेफिनिशन को खुद भरो, फिर नीचे "बाद में" टैब में समाधान के विरुद्ध अपने काम की जांच करो।

=== "पहले"

    ```groovy title="modules/samtools_index.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Generate BAM index file
     */
    process SAMTOOLS_INDEX {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "बाद में"

    ```groovy title="modules/samtools_index.nf" linenums="1" hl_lines="8 11 14 18"
    #!/usr/bin/env nextflow

    /*
     * Generate BAM index file
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

एक बार जब तुम इसे पूरा कर लो, तो प्रोसेस पूर्ण हो जाता है।
वर्कफ़्लो में इसका उपयोग करने के लिए, तुम्हें मॉड्यूल इम्पोर्ट करना होगा और एक प्रोसेस कॉल जोड़नी होगी।

#### 1.2.2. मॉड्यूल इन्क्लूड करना

`genomics.nf` में, प्रोसेस को वर्कफ़्लो के लिए उपलब्ध कराने के लिए एक `include` स्टेटमेंट जोड़ो:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    ```

प्रोसेस अब वर्कफ़्लो स्कोप में उपलब्ध है।

#### 1.2.3. इनपुट पर इंडेक्सिंग प्रोसेस कॉल करना

अब, वर्कफ़्लो ब्लॉक में `SAMTOOLS_INDEX` के लिए एक कॉल जोड़ें, इनपुट चैनल को आर्ग्युमेंट के रूप में पास करते हुए।

=== "बाद में"

    ```groovy title="genomics.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="14"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Call processes
    ```

वर्कफ़्लो अब इनपुट लोड करता है और उस पर इंडेक्सिंग प्रोसेस चलाता है।
अगला, हमें कॉन्फ़िगर करना होगा कि आउटपुट कैसे पब्लिश किया जाए।

### 1.3. आउटपुट हैंडलिंग कॉन्फ़िगर करना

हमें यह डिक्लेयर करना होगा कि कौन से प्रोसेस आउटपुट्स पब्लिश करने हैं और निर्दिष्ट करना होगा कि वे कहाँ जाने चाहिए।

#### 1.3.1. `publish:` सेक्शन में आउटपुट डिक्लेयर करना

वर्कफ़्लो ब्लॉक के अंदर `publish:` सेक्शन डिक्लेयर करता है कि कौन से प्रोसेस आउटपुट्स पब्लिश किए जाने चाहिए।
`SAMTOOLS_INDEX` के आउटपुट को `bam_index` नामक टार्गेट को असाइन करो।

=== "बाद में"

    ```groovy title="genomics.nf" linenums="22" hl_lines="2"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="22"
        publish:
        // Declare outputs to publish
    }
    ```

अब हमें Nextflow को बताना होगा कि पब्लिश किए गए आउटपुट को कहाँ रखना है।

#### 1.3.2. `output {}` ब्लॉक में आउटपुट टार्गेट कॉन्फ़िगर करना

`output {}` ब्लॉक वर्कफ़्लो के बाहर बैठता है और निर्दिष्ट करता है कि प्रत्येक नामित टार्गेट कहाँ पब्लिश किया जाता है।
चलो `bam_index` के लिए एक टार्गेट जोड़ें जो `bam/` सबडायरेक्टरी में पब्लिश करता है।

=== "बाद में"

    ```groovy title="genomics.nf" linenums="26" hl_lines="2-4"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="26"
    output {
        // Configure publish targets
    }
    ```

!!! note

    डिफ़ॉल्ट रूप से, Nextflow आउटपुट फ़ाइलों को सिम्बोलिक लिंक्स के रूप में पब्लिश करता है, जो अनावश्यक डुप्लीकेशन से बचता है।
    भले ही हम यहाँ जिन डेटा फ़ाइलों का उपयोग कर रहे हैं वे बहुत छोटी हैं, जीनोमिक्स में वे बहुत बड़ी हो सकती हैं।
    जब तुम अपनी `work` डायरेक्टरी को साफ करते हो तो सिम्लिंक्स टूट जाएंगे, इसलिए प्रोडक्शन वर्कफ़्लो के लिए तुम डिफ़ॉल्ट पब्लिश मोड को `'copy'` में ओवरराइड करना चाह सकते हो।

### 1.4. वर्कफ़्लो चलाना

इस बिंदु पर, हमारे पास एक वन-स्टेप इंडेक्सिंग वर्कफ़्लो है जो पूरी तरह से फंक्शनल होना चाहिए। चलो टेस्ट करें कि यह काम करता है!

हम इसे `-profile test` के साथ चला सकते हैं ताकि टेस्ट प्रोफ़ाइल में सेट अप की गई डिफ़ॉल्ट वैल्यू का उपयोग किया जा सके और कमांड लाइन पर पाथ लिखने से बचा जा सके।

```bash
nextflow run genomics.nf -profile test
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

तुम जांच सकते हो कि इंडेक्स फ़ाइल सही तरीके से जेनरेट हुई है या नहीं, work डायरेक्टरी या results डायरेक्टरी में देखकर।

??? abstract "Work डायरेक्टरी कंटेंट्स"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Results डायरेक्टरी कंटेंट्स"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

वहाँ है!

### सारांश

तुम जानते हो कि एक प्रोसेस युक्त मॉड्यूल कैसे बनाएं, इसे वर्कफ़्लो में इम्पोर्ट करें, इसे इनपुट चैनल के साथ कॉल करें, और परिणाम पब्लिश करें।

### आगे क्या है?

एक दूसरा स्टेप जोड़ो जो इंडेक्सिंग प्रोसेस के आउटपुट को लेता है और वेरिएंट कॉलिंग चलाने के लिए इसका उपयोग करता है।

---

## 2. GATK HaplotypeCaller को इंडेक्स्ड BAM फ़ाइल पर चलाने के लिए दूसरा प्रोसेस जोड़ना

अब जब हमारे पास अपनी इनपुट फ़ाइल के लिए एक इंडेक्स है, तो हम वेरिएंट कॉलिंग स्टेप सेट अप करने की ओर बढ़ सकते हैं।

[भाग 1](01_method.md) से `gatk HaplotypeCaller` कमांड को याद करो:

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

यह कमांड एक BAM फ़ाइल (`-I`), एक रेफरेंस जीनोम (`-R`), और एक intervals फ़ाइल (`-L`) लेता है, और एक VCF फ़ाइल (`-O`) उसके इंडेक्स के साथ प्रोड्यूस करता है।
टूल BAM इंडेक्स, रेफरेंस इंडेक्स, और रेफरेंस डिक्शनरी को उनकी संबंधित फ़ाइलों के साथ को-लोकेटेड होने की भी अपेक्षा करता है।
कंटेनर URI था `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`।

हम पहले की तरह ही तीन चरणों का पालन करते हैं:

1. इनपुट्स सेट अप करना
2. वेरिएंट कॉलिंग प्रोसेस लिखना और वर्कफ़्लो में उसे कॉल करना
3. आउटपुट हैंडलिंग कॉन्फ़िगर करना

### 2.1. इनपुट्स सेट अप करना

वेरिएंट कॉलिंग स्टेप को कई अतिरिक्त इनपुट फ़ाइलों की आवश्यकता होती है।
हमें उनके लिए पैरामीटर्स डिक्लेयर करने होंगे, टेस्ट प्रोफ़ाइल में डिफ़ॉल्ट वैल्यूज़ जोड़नी होंगी, और उन्हें लोड करने के लिए वेरिएबल्स बनाने होंगे।

#### 2.1.1. सहायक इनपुट्स के लिए पैरामीटर डिक्लेरेशन जोड़ना

चूंकि हमारा नया प्रोसेस कुछ अतिरिक्त फ़ाइलों की अपेक्षा करता है, `genomics.nf` में `Pipeline parameters` सेक्शन के अंतर्गत उनके लिए पैरामीटर डिक्लेरेशन जोड़ो:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="9" hl_lines="5-9"
    params {
        // Primary input
        reads_bam: Path

        // Accessory files
        reference: Path
        reference_index: Path
        reference_dict: Path
        intervals: Path
    }
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="9"
    params {
        // Primary input
        reads_bam: Path
    }
    ```

पहले की तरह, हम इनलाइन के बजाय टेस्ट प्रोफ़ाइल के माध्यम से डिफ़ॉल्ट वैल्यूज़ प्रदान करते हैं।

#### 2.1.2. टेस्ट प्रोफ़ाइल में सहायक फ़ाइल डिफ़ॉल्ट्स जोड़ना

जैसे हमने सेक्शन 1.1.2 में `reads_bam` के लिए किया था, `nextflow.config` में टेस्ट प्रोफ़ाइल में सहायक फ़ाइलों के लिए डिफ़ॉल्ट वैल्यूज़ जोड़ो:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
    }
    ```

अब हमें वेरिएबल्स बनाने होंगे जो वर्कफ़्लो में उपयोग के लिए इन फ़ाइल पाथ्स को लोड करें।

#### 2.1.3. सहायक फ़ाइलों के लिए वेरिएबल्स बनाना

वर्कफ़्लो ब्लॉक के अंदर सहायक फ़ाइल पाथ्स के लिए वेरिएबल्स जोड़ो:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="21" hl_lines="7-11"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Load the file paths for the accessory files (reference and intervals)
        ref_file        = file(params.reference)
        ref_index_file  = file(params.reference_index)
        ref_dict_file   = file(params.reference_dict)
        intervals_file  = file(params.intervals)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="21"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

`file()` सिंटैक्स Nextflow को स्पष्ट रूप से बताता है कि इन इनपुट्स को फ़ाइल पाथ्स के रूप में हैंडल करे।
तुम इसके बारे में Side Quest [Working with files](../../side_quests/working_with_files.md) में अधिक जान सकते हो।

### 2.2. वेरिएंट कॉलिंग प्रोसेस लिखना और वर्कफ़्लो में उसे कॉल करना

हमें मॉड्यूल फ़ाइल में प्रोसेस डेफिनिशन लिखनी होगी, include स्टेटमेंट का उपयोग करके इसे वर्कफ़्लो में इम्पोर्ट करना होगा, और इनपुट रीड्स प्लस इंडेक्सिंग स्टेप के आउटपुट और सहायक फ़ाइलों पर इसे कॉल करना होगा।

#### 2.2.1. वेरिएंट कॉलिंग प्रोसेस के लिए मॉड्यूल भरना

`modules/gatk_haplotypecaller.nf` खोलो और प्रोसेस डेफिनिशन की रूपरेखा की जांच करो।

आगे बढ़ो और ऊपर दी गई जानकारी का उपयोग करके प्रोसेस डेफिनिशन को खुद भरो, फिर नीचे "बाद में" टैब में समाधान के विरुद्ध अपने काम की जांच करो।

=== "पहले"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Call variants with GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "बाद में"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1" hl_lines="8 11-16 19-20 24-28"
    #!/usr/bin/env nextflow

    /*
     * Call variants with GATK HaplotypeCaller
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

तुम देखोगे कि इस प्रोसेस में GATK कमांड की आवश्यकता से अधिक इनपुट्स हैं।
GATK नामकरण परंपराओं के आधार पर BAM इंडेक्स फ़ाइल और रेफरेंस जीनोम की सहायक फ़ाइलों को खोजना जानता है, लेकिन Nextflow डोमेन-अज्ञेयवादी है और इन परंपराओं के बारे में नहीं जानता।
हमें उन्हें स्पष्ट रूप से सूचीबद्ध करना होगा ताकि Nextflow उन्हें रनटाइम पर वर्किंग डायरेक्टरी में स्टेज करे; अन्यथा GATK लापता फ़ाइलों के बारे में एक एरर थ्रो करेगा।

इसी तरह, हम आउटपुट VCF की इंडेक्स फ़ाइल (`"${input_bam}.vcf.idx"`) को स्पष्ट रूप से सूचीबद्ध करते हैं ताकि Nextflow बाद के स्टेप्स के लिए इसका ट्रैक रखे।
हम प्रत्येक आउटपुट चैनल को एक नाम असाइन करने के लिए `emit:` सिंटैक्स का उपयोग करते हैं, जो तब उपयोगी हो जाएगा जब हम आउटपुट्स को publish ब्लॉक में वायर करेंगे।

एक बार जब तुम इसे पूरा कर लो, तो प्रोसेस पूर्ण हो जाता है।
वर्कफ़्लो में इसका उपयोग करने के लिए, तुम्हें मॉड्यूल इम्पोर्ट करना होगा और एक प्रोसेस कॉल जोड़नी होगी।

#### 2.2.2. नया मॉड्यूल इम्पोर्ट करना

नए मॉड्यूल को इम्पोर्ट करने के लिए `genomics.nf` को अपडेट करो:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

प्रोसेस अब वर्कफ़्लो स्कोप में उपलब्ध है।

#### 2.2.3. प्रोसेस कॉल जोड़ना

वर्कफ़्लो बॉडी में, `main:` के अंतर्गत प्रोसेस कॉल जोड़ो:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="33" hl_lines="4-12"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
            ref_file,
            ref_index_file,
            ref_dict_file,
            intervals_file
        )
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="33"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

तुम्हें Hello Nextflow ट्रेनिंग सीरीज़ से `*.out` सिंटैक्स को पहचानना चाहिए; हम Nextflow को बता रहे हैं कि `SAMTOOLS_INDEX` द्वारा आउटपुट किए गए चैनल को लें और उसे `GATK_HAPLOTYPECALLER` प्रोसेस कॉल में प्लग करें।

!!! note

    ध्यान दें कि इनपुट्स प्रोसेस की कॉल में उसी क्रम में प्रदान किए जाते हैं जैसे वे प्रोसेस के input ब्लॉक में सूचीबद्ध हैं।
    Nextflow में, इनपुट्स पोज़िशनल हैं, जिसका अर्थ है कि तुम्हें _उसी क्रम का पालन करना चाहिए_; और निश्चित रूप से तत्वों की संख्या समान होनी चाहिए।

### 2.3. आउटपुट हैंडलिंग कॉन्फ़िगर करना

हमें publish डिक्लेरेशन में नए आउटपुट्स जोड़ने होंगे और कॉन्फ़िगर करना होगा कि वे कहाँ जाते हैं।

#### 2.3.1. वेरिएंट कॉलिंग आउटपुट्स के लिए publish टार्गेट्स जोड़ना

`publish:` सेक्शन में VCF और इंडेक्स आउटपुट्स जोड़ो:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="45" hl_lines="3-4"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="45"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

अब हमें Nextflow को बताना होगा कि नए आउटपुट्स को कहाँ रखना है।

#### 2.3.2. नए आउटपुट टार्गेट्स कॉन्फ़िगर करना

`output {}` ब्लॉक में `vcf` और `vcf_idx` टार्गेट्स के लिए एंट्रीज़ जोड़ो, दोनों को `vcf/` सबडायरेक्टरी में पब्लिश करते हुए:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="51" hl_lines="5-10"
    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="49"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

VCF और उसका इंडेक्स अलग-अलग टार्गेट्स के रूप में पब्लिश किए जाते हैं जो दोनों `vcf/` सबडायरेक्टरी में जाते हैं।

### 2.4. वर्कफ़्लो चलाना

विस्तारित वर्कफ़्लो चलाओ, इस बार `-resume` जोड़ते हुए ताकि हमें इंडेक्सिंग स्टेप फिर से न चलाना पड़े।

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

अब अगर हम कंसोल आउटपुट को देखें, तो हम दो प्रोसेसेज़ सूचीबद्ध देखते हैं।

पहला प्रोसेस कैशिंग के कारण स्किप किया गया था, जैसी अपेक्षा थी, जबकि दूसरा प्रोसेस चलाया गया क्योंकि यह बिल्कुल नया है।

तुम results डायरेक्टरी में आउटपुट फ़ाइलें पाओगे (work डायरेक्टरी के सिम्बोलिक लिंक्स के रूप में)।

??? abstract "डायरेक्टरी कंटेंट्स"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

अगर तुम VCF फ़ाइल खोलते हो, तो तुम्हें वही कंटेंट्स दिखने चाहिए जो उस फ़ाइल में थे जो तुमने कंटेनर में सीधे GATK कमांड चलाकर जेनरेट की थी।

??? abstract "फ़ाइल कंटेंट्स"

    ```console title="reads_mother.bam.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

यह वह आउटपुट है जिसे हम अपने अध्ययन में प्रत्येक नमूने के लिए जेनरेट करने की परवाह करते हैं।

### सारांश

तुम जानते हो कि दो-स्टेप मॉड्यूलर वर्कफ़्लो कैसे बनाएं जो वास्तविक विश्लेषण कार्य करता है और जीनोमिक्स फ़ाइल फ़ॉर्मेट की विशिष्टताओं जैसे सहायक फ़ाइलों से निपटने में सक्षम है।

### आगे क्या है?

वर्कफ़्लो को बल्क में कई नमूनों को हैंडल करने योग्य बनाओ।

---

## 3. वर्कफ़्लो को नमूनों के बैच पर चलाने के लिए अनुकूलित करना

एक वर्कफ़्लो होना जो एक सिंगल नमूने पर प्रोसेसिंग को ऑटोमेट कर सके, यह सब अच्छा है, लेकिन क्या होगा अगर तुम्हारे पास 1000 नमूने हैं?
क्या तुम्हें एक bash स्क्रिप्ट लिखनी होगी जो तुम्हारे सभी नमूनों के माध्यम से लूप करे?

नहीं, भगवान का शुक्र है! बस कोड में एक छोटा सा बदलाव करो और Nextflow तुम्हारे लिए वह भी हैंडल कर लेगा।

### 3.1. तीन नमूनों को सूचीबद्ध करने के लिए इनपुट अपडेट करना

कई नमूनों पर चलाने के लिए, टेस्ट प्रोफ़ाइल को एक सिंगल के बजाय फ़ाइल पाथ्स की एक ऐरे प्रदान करने के लिए अपडेट करो।
यह मल्टी-सैंपल एक्ज़ीक्यूशन को टेस्ट करने का एक त्वरित तरीका है; अगले स्टेप में हम इनपुट्स की एक फ़ाइल का उपयोग करके अधिक स्केलेबल दृष्टिकोण पर स्विच करेंगे।

पहले, पैरामीटर डिक्लेरेशन में टाइप एनोटेशन को कमेंट आउट करो, क्योंकि ऐरे टाइप्ड डिक्लेरेशन का उपयोग नहीं कर सकते:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (array of three samples)
        reads_bam //: Path
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input
        reads_bam: Path
    ```

फिर सभी तीन नमूनों को सूचीबद्ध करने के लिए टेस्ट प्रोफ़ाइल को अपडेट करो:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2-6"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

वर्कफ़्लो बॉडी में चैनल फ़ैक्टरी (`.fromPath`) एक सिंगल के साथ-साथ कई फ़ाइल पाथ्स को भी स्वीकार करती है, इसलिए कोई अन्य बदलाव की आवश्यकता नहीं है।

### 3.2. वर्कफ़्लो चलाना

अब वर्कफ़्लो चलाने की कोशिश करो जब प्लंबिंग सभी तीन टेस्ट नमूनों पर चलाने के लिए सेट अप है।

```bash
nextflow run genomics.nf -profile test -resume
```

मज़ेदार बात: यह _काम कर सकता है_, या यह _फेल हो सकता है_। उदाहरण के लिए, यहाँ एक रन है जो सफल रहा:

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

अगर तुम्हारा वर्कफ़्लो रन सफल रहा, तो इसे फिर से चलाओ जब तक तुम्हें इस तरह का एरर न मिले:

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

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

अगर तुम GATK कमांड एरर आउटपुट को देखो, तो इस तरह की एक लाइन होगी:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

खैर, यह अजीब है, यह देखते हुए कि हमने वर्कफ़्लो के पहले स्टेप में स्पष्ट रूप से BAM फ़ाइलों को इंडेक्स किया था। क्या प्लंबिंग में कुछ गड़बड़ हो सकती है?

### 3.3. समस्या का निवारण करना

हम work डायरेक्टरीज़ का निरीक्षण करेंगे और यह पता लगाने के लिए `view()` ऑपरेटर का उपयोग करेंगे कि क्या गलत हुआ।

#### 3.3.1. संबंधित कॉल्स के लिए work डायरेक्टरीज़ की जांच करना

कंसोल आउटपुट में सूचीबद्ध फेल हुए `GATK_HAPLOTYPECALLER` प्रोसेस कॉल के लिए work डायरेक्टरी के अंदर देखो।

??? abstract "डायरेक्टरी कंटेंट्स"

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

इस डायरेक्टरी में सूचीबद्ध BAM फ़ाइल और BAM इंडेक्स के नामों पर विशेष ध्यान दो: `reads_son.bam` और `reads_father.bam.bai`।

यह क्या है? Nextflow ने इस प्रोसेस कॉल की work डायरेक्टरी में एक इंडेक्स फ़ाइल स्टेज की है, लेकिन यह गलत है। यह कैसे हो सकता है?

#### 3.3.2. चैनल कंटेंट्स का निरीक्षण करने के लिए [view() ऑपरेटर](https://www.nextflow.io/docs/latest/reference/operator.html#view) का उपयोग करना

`GATK_HAPLOTYPECALLER` प्रोसेस कॉल से पहले वर्कफ़्लो बॉडी में ये दो लाइनें जोड़ो ताकि चैनल के कंटेंट्स को देख सको:

=== "बाद में"

    ```groovy title="genomics.nf" hl_lines="3-5"
        SAMTOOLS_INDEX(reads_ch)

        // temporary diagnostics
        reads_ch.view()
        SAMTOOLS_INDEX.out.view()

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

=== "पहले"

    ```groovy title="genomics.nf"
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

फिर वर्कफ़्लो कमांड फिर से चलाओ।

```bash
nextflow run genomics.nf -profile test
```

एक बार फिर, यह सफल हो सकता है या फेल हो सकता है। यहाँ फेल हुए रन के लिए दो `.view()` कॉल्स का आउटपुट है:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

पहली तीन लाइनें इनपुट चैनल से संबंधित हैं और दूसरी, आउटपुट चैनल से।
तुम देख सकते हो कि तीन नमूनों के लिए BAM फ़ाइलें और इंडेक्स फ़ाइलें एक ही क्रम में सूचीबद्ध नहीं हैं!

!!! note

    जब तुम कई तत्वों वाले चैनल पर Nextflow प्रोसेस कॉल करते हो, तो Nextflow निष्पादन को यथासंभव समानांतर करने की कोशिश करेगा, और जो भी क्रम में उपलब्ध हो जाएं उस क्रम में आउटपुट एकत्र करेगा।
    परिणाम यह है कि संबंधित आउटपुट्स मूल इनपुट्स के दिए गए क्रम से अलग क्रम में एकत्र किए जा सकते हैं।

जैसा कि वर्तमान में लिखा गया है, हमारी वर्कफ़्लो स्क्रिप्ट मानती है कि इंडेक्स फ़ाइलें इंडेक्सिंग स्टेप से उसी mother/father/son क्रम में सूचीबद्ध होकर बाहर आएंगी जैसे इनपुट्स दिए गए थे।
लेकिन यह गारंटी नहीं है कि ऐसा होगा, यही कारण है कि कभी-कभी (हालांकि हमेशा नहीं) गलत फ़ाइलें दूसरे स्टेप में जोड़ी जाती हैं।

इसे ठीक करने के लिए, हमें यह सुनिश्चित करना होगा कि BAM फ़ाइलें और उनकी इंडेक्स फ़ाइलें चैनलों के माध्यम से एक साथ यात्रा करें।

!!! tip

    वर्कफ़्लो कोड में `view()` स्टेटमेंट्स कुछ नहीं करते, इसलिए उन्हें छोड़ना कोई समस्या नहीं है।
    हालांकि वे तुम्हारे कंसोल आउटपुट को अव्यवस्थित करेंगे, इसलिए हम सुझाव देते हैं कि जब तुम समस्या का निवारण कर लो तो उन्हें हटा दो।

### 3.4. इंडेक्स फ़ाइलों को सही तरीके से हैंडल करने के लिए वर्कफ़्लो अपडेट करना

फिक्स यह है कि प्रत्येक BAM फ़ाइल को उसके इंडेक्स के साथ एक टपल में बंडल करें, फिर डाउनस्ट्रीम प्रोसेस और वर्कफ़्लो प्लंबिंग को मैच करने के लिए अपडेट करें।

#### 3.4.1. SAMTOOLS_INDEX मॉड्यूल के आउटपुट को टपल में बदलना

यह सुनिश्चित करने का सबसे सरल तरीका कि एक BAM फ़ाइल और उसका इंडेक्स निकटता से जुड़े रहें, उन्हें इंडेक्स टास्क से बाहर आने वाले टपल में एक साथ पैकेज करना है।

!!! note

    एक **tuple** तत्वों की एक परिमित, क्रमबद्ध सूची है जो आमतौर पर एक फ़ंक्शन से कई वैल्यूज़ रिटर्न करने के लिए उपयोग की जाती है। टपल्स विशेष रूप से उपयोगी हैं जब प्रोसेसेज़ के बीच कई इनपुट्स या आउटपुट्स को उनके संबंध और क्रम को संरक्षित करते हुए पास करना हो।

`modules/samtools_index.nf` में आउटपुट को BAM फ़ाइल शामिल करने के लिए अपडेट करो:

=== "बाद में"

    ```groovy title="modules/samtools_index.nf" linenums="14" hl_lines="2"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "पहले"

    ```groovy title="modules/samtools_index.nf" linenums="14"
        output:
        path "${input_bam}.bai"
    ```

इस तरह, प्रत्येक इंडेक्स फ़ाइल अपनी मूल BAM फ़ाइल के साथ कसकर जुड़ी होगी, और इंडेक्सिंग स्टेप का समग्र आउटपुट फ़ाइलों के जोड़े युक्त एक सिंगल चैनल होगा।

#### 3.4.2. GATK_HAPLOTYPECALLER मॉड्यूल के इनपुट को टपल स्वीकार करने के लिए बदलना

चूंकि हमने पहले प्रोसेस के आउटपुट के 'शेप' को बदल दिया है, हमें दूसरे प्रोसेस की इनपुट डेफिनिशन को मैच करने के लिए अपडेट करना होगा।

`modules/gatk_haplotypecaller.nf` अपडेट करो:

=== "बाद में"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11" hl_lines="2"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "पहले"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11"
        input:
        path input_bam
        path input_bam_index
    ```

अब हमें प्रोसेस कॉल और publish टार्गेट्स में नई टपल संरचना को प्रतिबिंबित करने के लिए वर्कफ़्लो को अपडेट करना होगा।

#### 3.4.3. वर्कफ़्लो में GATK_HAPLOTYPECALLER की कॉल अपडेट करना

हमें अब `GATK_HAPLOTYPECALLER` प्रोसेस को मूल `reads_ch` प्रदान करने की आवश्यकता नहीं है, क्योंकि BAM फ़ाइल अब `SAMTOOLS_INDEX` द्वारा आउटपुट किए गए चैनल में बंडल है।

`genomics.nf` में कॉल अपडेट करो:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="42" hl_lines="2"
        GATK_HAPLOTYPECALLER(
            SAMTOOLS_INDEX.out,
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="42"
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
    ```

अंत में, हमें नई आउटपुट संरचना को प्रतिबिंबित करने के लिए publish टार्गेट्स को अपडेट करना होगा।

#### 3.4.4. इंडेक्स्ड BAM आउटपुट के लिए publish टार्गेट अपडेट करना

चूंकि SAMTOOLS_INDEX आउटपुट अब BAM फ़ाइल और उसके इंडेक्स दोनों युक्त एक टपल है, publish टार्गेट को `bam_index` से `indexed_bam` में रीनेम करो ताकि इसकी कंटेंट्स को बेहतर तरीके से प्रतिबिंबित किया जा सके:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="46" hl_lines="2 8"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="46"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

इन बदलावों के साथ, BAM और उसका इंडेक्स एक साथ यात्रा करने की गारंटी है, इसलिए जोड़ी हमेशा सही होगी।

### 3.5. सही किए गए वर्कफ़्लो को चलाना

यह सुनिश्चित करने के लिए वर्कफ़्लो फिर से चलाओ कि यह आगे विश्वसनीय रूप से काम करेगा।

```bash
nextflow run genomics.nf -profile test
```

इस बार (और हर बार) सब कुछ सही तरीके से चलना चाहिए:

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

results डायरेक्टरी में अब प्रत्येक नमूने के लिए BAM और BAI फ़ाइलें (टपल से) शामिल हैं, साथ ही VCF आउटपुट्स:

??? abstract "Results डायरेक्टरी कंटेंट्स"

    ```console
    results/
    ├── bam/
    │   ├── reads_father.bam -> ...
    │   ├── reads_father.bam.bai -> ...
    │   ├── reads_mother.bam -> ...
    │   ├── reads_mother.bam.bai -> ...
    │   ├── reads_son.bam -> ...
    │   └── reads_son.bam.bai -> ...
    └── vcf/
        ├── reads_father.bam.vcf -> ...
        ├── reads_father.bam.vcf.idx -> ...
        ├── reads_mother.bam.vcf -> ...
        ├── reads_mother.bam.vcf.idx -> ...
        ├── reads_son.bam.vcf -> ...
        └── reads_son.bam.vcf.idx -> ...
    ```

संबंधित फ़ाइलों को टपल्स में बंडल करके, हमने सुनिश्चित किया कि सही फ़ाइलें हमेशा वर्कफ़्लो के माध्यम से एक साथ यात्रा करती हैं।
वर्कफ़्लो अब किसी भी संख्या में नमूनों को विश्वसनीय रूप से प्रोसेस करता है, लेकिन उन्हें config में व्यक्तिगत रूप से सूचीबद्ध करना बहुत स्केलेबल नहीं है।
अगले स्टेप में, हम एक फ़ाइल से इनपुट्स पढ़ने पर स्विच करेंगे।

### सारांश

तुम जानते हो कि अपने वर्कफ़्लो को कई नमूनों पर (स्वतंत्र रूप से) कैसे चलाएं।

### आगे क्या है?

बल्क में नमूनों को हैंडल करना आसान बनाओ।

---

## 4. वर्कफ़्लो को इनपुट फ़ाइलों के बैच युक्त टेक्स्ट फ़ाइल स्वीकार करने योग्य बनाना

वर्कफ़्लो को कई डेटा इनपुट फ़ाइलें प्रदान करने का एक बहुत सामान्य तरीका फ़ाइल पाथ्स युक्त टेक्स्ट फ़ाइल के साथ करना है।
यह एक सरल टेक्स्ट फ़ाइल हो सकती है जो प्रति लाइन एक फ़ाइल पाथ सूचीबद्ध करती है और कुछ नहीं, या फ़ाइल में अतिरिक्त मेटाडेटा हो सकता है, जिस स्थिति में इसे अक्सर samplesheet कहा जाता है।

यहाँ हम तुम्हें दिखाने जा रहे हैं कि सरल केस कैसे करें।

### 4.1. इनपुट फ़ाइल पाथ्स सूचीबद्ध करने वाली प्रदान की गई टेक्स्ट फ़ाइल की जांच करना

हमने पहले से ही इनपुट फ़ाइल पाथ्स सूचीबद्ध करने वाली एक टेक्स्ट फ़ाइल बनाई है, जिसे `sample_bams.txt` कहा जाता है, जिसे तुम `data/` डायरेक्टरी में पा सकते हो।

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

जैसा कि तुम देख सकते हो, हमने प्रति लाइन एक फ़ाइल पाथ सूचीबद्ध किया है, और वे एब्सोल्यूट पाथ्स हैं।

!!! note

    हम यहाँ जिन फ़ाइलों का उपयोग कर रहे हैं वे सिर्फ तुम्हारे GitHub Codespaces के लोकल फ़ाइलसिस्टम पर हैं, लेकिन हम क्लाउड स्टोरेज में फ़ाइलों की ओर भी इशारा कर सकते हैं।
    अगर तुम प्रदान किए गए Codespaces एनवायरनमेंट का उपयोग नहीं कर रहे हो, तो तुम्हें अपने लोकल सेटअप से मैच करने के लिए फ़ाइल पाथ्स को अनुकूलित करने की आवश्यकता हो सकती है।

### 4.2. पैरामीटर और टेस्ट प्रोफ़ाइल अपडेट करना

व्यक्तिगत नमूनों को सूचीबद्ध करने के बजाय `reads_bam` पैरामीटर को `sample_bams.txt` फ़ाइल की ओर इशारा करने के लिए स्विच करो।

params ब्लॉक में टाइप एनोटेशन रिस्टोर करो (क्योंकि यह फिर से एक सिंगल पाथ है):

=== "बाद में"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (file of input files, one per line)
        reads_bam: Path
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input (array of three samples)
        reads_bam
    ```

फिर टेक्स्ट फ़ाइल की ओर इशारा करने के लिए टेस्ट प्रोफ़ाइल को अपडेट करो:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2"
    test {
        params.reads_bam = "${projectDir}/data/sample_bams.txt"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

फ़ाइलों की सूची अब कोड में बिल्कुल नहीं रहती, जो सही दिशा में एक बड़ा कदम है।

### 4.3. फ़ाइल से लाइनें पढ़ने के लिए चैनल फ़ैक्टरी अपडेट करना

वर्तमान में, हमारी इनपुट चैनल फ़ैक्टरी हमें दी गई किसी भी फ़ाइल को डेटा इनपुट के रूप में ट्रीट करती है जिसे हम इंडेक्सिंग प्रोसेस को फीड करना चाहते हैं।
चूंकि हम अब इसे एक फ़ाइल दे रहे हैं जो इनपुट फ़ाइल पाथ्स सूचीबद्ध करती है, हमें फ़ाइल को पार्स करने और इसमें शामिल फ़ाइल पाथ्स को डेटा इनपुट के रूप में ट्रीट करने के लिए इसके व्यवहार को बदलने की आवश्यकता है।

हम उसी पैटर्न का उपयोग करके ऐसा कर सकते हैं जिसका हमने [Hello Nextflow के भाग 2](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file) में उपयोग किया था: फ़ाइल को पार्स करने के लिए [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) ऑपरेटर लागू करना, फिर प्रत्येक लाइन के पहले फ़ील्ड को चुनने के लिए एक `map` ऑपरेशन।

=== "बाद में"

    ```groovy title="genomics.nf" linenums="24" hl_lines="1-4"
        // Create input channel from a CSV file listing input file paths
        reads_ch = Channel.fromPath(params.reads_bam)
                .splitCsv()
                .map { line -> file(line[0]) }
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="24"
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

तकनीकी रूप से हम इसे [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) ऑपरेटर का उपयोग करके अधिक सरलता से कर सकते हैं, क्योंकि हमारी इनपुट फ़ाइल वर्तमान में केवल फ़ाइल पाथ्स शामिल करती है।
हालांकि, अधिक बहुमुखी `splitCsv` ऑपरेटर (`map` द्वारा पूरक) का उपयोग करके, हम अपने वर्कफ़्लो को भविष्य-प्रूफ कर सकते हैं यदि हम फ़ाइल पाथ्स युक्त फ़ाइल में मेटाडेटा जोड़ने का निर्णय लेते हैं।

!!! tip

    अगर तुम्हें विश्वास नहीं है कि तुम समझते हो कि ऑपरेटर्स यहाँ क्या कर रहे हैं, तो यह `.view()` ऑपरेटर का उपयोग करके यह देखने का एक और शानदार अवसर है कि उन्हें लागू करने से पहले और बाद में चैनल कंटेंट्स कैसे दिखते हैं।

### 4.4. वर्कफ़्लो चलाना

वर्कफ़्लो एक बार और चलाओ। इससे पहले जैसा ही परिणाम मिलना चाहिए, सही?

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

हाँ! वास्तव में, Nextflow सही तरीके से पता लगाता है कि प्रोसेस कॉल्स बिल्कुल समान हैं, और सब कुछ फिर से चलाने की जहमत भी नहीं उठाता, क्योंकि हम `-resume` के साथ चला रहे थे।

और बस! हमारे सरल वेरिएंट कॉलिंग वर्कफ़्लो में वे सभी बुनियादी फीचर्स हैं जो हम चाहते थे।

### सारांश

तुम जानते हो कि GATK का उपयोग करके BAM फ़ाइल को इंडेक्स करने और प्रति-नमूना वेरिएंट कॉलिंग लागू करने के लिए मल्टी-स्टेप मॉड्यूलर वर्कफ़्लो कैसे बनाएं।

अधिक सामान्य रूप से, तुमने सीखा है कि जीनोमिक्स फ़ाइल फ़ॉर्मेट्स और टूल आवश्यकताओं की विशिष्टताओं को ध्यान में रखते हुए वास्तविक काम करने वाली एक सरल जीनोमिक्स पाइपलाइन बनाने के लिए आवश्यक Nextflow कंपोनेंट्स और लॉजिक का उपयोग कैसे करें।

### आगे क्या है?

अपनी सफलता का जश्न मनाओ और एक अतिरिक्त लंबा ब्रेक लो!

इस कोर्स के अगले भाग में, तुम सीखोगे कि डेटा पर joint वेरिएंट कॉलिंग लागू करने के लिए इस सरल प्रति-नमूना वेरिएंट कॉलिंग वर्कफ़्लो को कैसे ट्रांसफ़ॉर्म करें।
