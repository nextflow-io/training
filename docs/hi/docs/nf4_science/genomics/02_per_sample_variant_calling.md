# भाग 2: प्रति-नमूना वेरिएंट कॉलिंग

भाग 1 में, तुमने Samtools और GATK कमांड को उनके संबंधित कंटेनर में मैन्युअल रूप से टेस्ट किया।
अब हम उन्हीं कमांड को एक Nextflow workflow में रैप करने जा रहे हैं।

## असाइनमेंट

इस कोर्स के इस भाग में, हम एक workflow विकसित करने जा रहे हैं जो निम्नलिखित करता है:

1. [Samtools](https://www.htslib.org/) का उपयोग करके प्रत्येक BAM इनपुट फ़ाइल के लिए एक इंडेक्स फ़ाइल जनरेट करना
2. प्रत्येक BAM इनपुट फ़ाइल पर GATK HaplotypeCaller चलाना ताकि VCF (Variant Call Format) में प्रति-नमूना वेरिएंट कॉल जनरेट हो सकें

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

यह भाग 1 के चरणों को दोहराता है, जहां तुमने इन कमांड को उनके कंटेनर में मैन्युअल रूप से चलाया था।

शुरुआती बिंदु के रूप में, हम तुम्हें एक workflow फ़ाइल, `genomics.nf`, प्रदान करते हैं, जो workflow के मुख्य भागों को रेखांकित करती है, साथ ही दो मॉड्यूल फ़ाइलें, samtools_index.nf और gatk_haplotypecaller.nf, जो मॉड्यूल की संरचना को रेखांकित करती हैं।
ये फ़ाइलें फ़ंक्शनल नहीं हैं; उनका उद्देश्य केवल स्कैफ़ोल्ड के रूप में काम करना है जिन्हें तुम कोड के दिलचस्प भागों से भर सकते हो।

## पाठ योजना

विकास प्रक्रिया को अधिक शैक्षिक बनाने के लिए, हमने इसे चार चरणों में विभाजित किया है:

1. **एक सिंगल-स्टेज workflow लिखो जो एक BAM फ़ाइल पर Samtools index चलाए।**
   यह एक मॉड्यूल बनाना, उसे इम्पोर्ट करना, और workflow में कॉल करना कवर करता है।
2. **इंडेक्स की गई BAM फ़ाइल पर GATK HaplotypeCaller चलाने के लिए एक दूसरी प्रोसेस जोड़ो।**
   यह प्रोसेस आउटपुट को इनपुट में चेन करना और एक्सेसरी फ़ाइलों को हैंडल करना परिचित कराता है।
3. **Workflow को नमूनों के बैच पर चलाने के लिए अनुकूलित करो।**
   यह समानांतर निष्पादन को कवर करता है और संबंधित फ़ाइलों को एक साथ रखने के लिए टपल परिचित कराता है।
4. **Workflow को एक टेक्स्ट फ़ाइल स्वीकार करने योग्य बनाओ जिसमें इनपुट फ़ाइलों का बैच हो।**
   यह थोक में इनपुट प्रदान करने के लिए एक सामान्य पैटर्न प्रदर्शित करता है।

प्रत्येक चरण workflow विकास के एक विशिष्ट पहलू पर केंद्रित है।

---

## 1. एक सिंगल-स्टेज workflow लिखो जो एक BAM फ़ाइल पर Samtools index चलाए

यह पहला चरण मूल बातों पर केंद्रित है: एक BAM फ़ाइल लोड करना और उसके लिए एक इंडेक्स जनरेट करना।

[भाग 1](01_method.md) से `samtools index` कमांड को याद करो:

```bash
samtools index '<input_bam>'
```

यह कमांड इनपुट के रूप में एक BAM फ़ाइल लेता है और उसके साथ एक `.bai` इंडेक्स फ़ाइल बनाता है।
कंटेनर URI था `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464`।

हम इस जानकारी को लेने और इसे तीन चरणों में Nextflow में रैप करने जा रहे हैं:

1. इनपुट सेटअप करो
2. इंडेक्सिंग प्रोसेस लिखो और workflow में इसे कॉल करो
3. आउटपुट हैंडलिंग को कॉन्फ़िगर करो

### 1.1. इनपुट सेटअप करो

हमें एक इनपुट पैरामीटर डिक्लेयर करना होगा, एक सुविधाजनक डिफ़ॉल्ट वैल्यू प्रदान करने के लिए एक टेस्ट प्रोफ़ाइल बनाना होगा, और एक इनपुट चैनल बनाना होगा।

#### 1.1.1. एक इनपुट पैरामीटर डिक्लेरेशन जोड़ो

मुख्य workflow फ़ाइल `genomics.nf` में, `Pipeline parameters` सेक्शन के तहत, `reads_bam` नाम का एक CLI पैरामीटर डिक्लेयर करो।

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

यह CLI पैरामीटर सेटअप करता है, लेकिन हम विकास के दौरान हर बार workflow चलाते समय फ़ाइल पथ टाइप नहीं करना चाहते।
डिफ़ॉल्ट वैल्यू प्रदान करने के लिए कई विकल्प हैं; यहां हम एक टेस्ट प्रोफ़ाइल का उपयोग करते हैं।

#### 1.1.2. `nextflow.config` में डिफ़ॉल्ट वैल्यू के साथ एक टेस्ट प्रोफ़ाइल बनाओ

एक टेस्ट प्रोफ़ाइल कमांड लाइन पर इनपुट निर्दिष्ट किए बिना workflow आज़माने के लिए सुविधाजनक डिफ़ॉल्ट वैल्यू प्रदान करती है।
यह Nextflow इकोसिस्टम में एक सामान्य परंपरा है (अधिक विवरण के लिए [Hello Config](../../hello_nextflow/06_hello_config.md) देखें)।

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

यहां, हम `${projectDir}` का उपयोग कर रहे हैं, एक बिल्ट-इन Nextflow वेरिएबल जो उस डायरेक्टरी की ओर इशारा करता है जहां workflow स्क्रिप्ट स्थित है।
यह एब्सोल्यूट पथ हार्डकोड किए बिना डेटा फ़ाइलों और अन्य संसाधनों को रेफ़रेंस करना आसान बनाता है।

#### 1.1.3. इनपुट चैनल सेटअप करो

workflow ब्लॉक में, `.fromPath` चैनल फ़ैक्ट्री का उपयोग करके पैरामीटर वैल्यू से एक इनपुट चैनल बनाओ (जैसा कि [Hello Channels](../../hello_nextflow/02_hello_channels.md) में उपयोग किया गया था)।

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

अब हमें इस इनपुट पर इंडेक्सिंग चलाने के लिए प्रोसेस बनाने की आवश्यकता है।

### 1.2. इंडेक्सिंग प्रोसेस लिखो और workflow में इसे कॉल करो

हमें मॉड्यूल फ़ाइल में प्रोसेस डेफिनिशन लिखनी होगी, include स्टेटमेंट का उपयोग करके इसे workflow में इम्पोर्ट करना होगा, और इनपुट पर इसे कॉल करना होगा।

#### 1.2.1. इंडेक्सिंग प्रोसेस के लिए मॉड्यूल भरो

`modules/samtools_index.nf` खोलो और प्रोसेस डेफिनिशन की रूपरेखा जांचो।
तुम्हें मुख्य संरचनात्मक तत्वों को पहचानना चाहिए; यदि नहीं, तो रिफ़्रेशर के लिए [Hello Nextflow](../../hello_nextflow/01_hello_world.md) पढ़ने पर विचार करो।

आगे बढ़ो और ऊपर दी गई जानकारी का उपयोग करके खुद प्रोसेस डेफिनिशन भरो, फिर नीचे "बाद में" टैब में समाधान के विरुद्ध अपना काम चेक करो।

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

एक बार जब तुम इसे पूरा कर लो, तो प्रोसेस पूर्ण हो जाती है।
इसे workflow में उपयोग करने के लिए, तुम्हें मॉड्यूल इम्पोर्ट करना होगा और एक प्रोसेस कॉल जोड़नी होगी।

#### 1.2.2. मॉड्यूल को इन्क्लूड करो

`genomics.nf` में, प्रोसेस को workflow के लिए उपलब्ध कराने के लिए एक `include` स्टेटमेंट जोड़ो:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    ```

अब प्रोसेस workflow स्कोप में उपलब्ध है।

#### 1.2.3. इनपुट पर इंडेक्सिंग प्रोसेस कॉल करो

अब, workflow ब्लॉक में `SAMTOOLS_INDEX` के लिए एक कॉल जोड़ें, इनपुट चैनल को आर्गुमेंट के रूप में पास करते हुए।

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

अब workflow इनपुट लोड करता है और उस पर इंडेक्सिंग प्रोसेस चलाता है।
आगे, हमें कॉन्फ़िगर करने की आवश्यकता है कि आउटपुट कैसे पब्लिश किया जाए।

### 1.3. आउटपुट हैंडलिंग को कॉन्फ़िगर करो

हमें यह डिक्लेयर करने की आवश्यकता है कि कौन से प्रोसेस आउटपुट पब्लिश करने हैं और निर्दिष्ट करना है कि वे कहां जाने चाहिए।

#### 1.3.1. `publish:` सेक्शन में एक आउटपुट डिक्लेयर करो

workflow ब्लॉक के अंदर `publish:` सेक्शन यह डिक्लेयर करता है कि कौन से प्रोसेस आउटपुट पब्लिश किए जाने चाहिए।
`SAMTOOLS_INDEX` के आउटपुट को `bam_index` नामक एक नामित टार्गेट को असाइन करो।

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

अब हमें Nextflow को यह बताना होगा कि पब्लिश किए गए आउटपुट को कहां रखना है।

#### 1.3.2. `output {}` ब्लॉक में आउटपुट टार्गेट को कॉन्फ़िगर करो

`output {}` ब्लॉक workflow के बाहर बैठता है और निर्दिष्ट करता है कि प्रत्येक नामित टार्गेट कहां पब्लिश किया जाता है।
`bam_index` के लिए एक टार्गेट जोड़ें जो `bam/` सबडायरेक्टरी में पब्लिश करे।

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

    डिफ़ॉल्ट रूप से, Nextflow आउटपुट फ़ाइलों को सिंबॉलिक लिंक के रूप में पब्लिश करता है, जो अनावश्यक डुप्लीकेशन से बचता है।
    भले ही हम यहां जिन डेटा फ़ाइलों का उपयोग कर रहे हैं वे बहुत छोटी हैं, genomics में वे बहुत बड़ी हो सकती हैं।
    जब तुम अपनी `work` डायरेक्टरी को साफ़ करते हो तो Symlinks टूट जाएंगे, इसलिए प्रोडक्शन workflow के लिए तुम डिफ़ॉल्ट पब्लिश मोड को `'copy'` में ओवरराइड करना चाह सकते हो।

### 1.4. Workflow चलाओ

इस बिंदु पर, हमारे पास एक वन-स्टेप इंडेक्सिंग workflow है जो पूरी तरह से फ़ंक्शनल होनी चाहिए। चलो परीक्षण करें कि यह काम करता है!

हम इसे `-profile test` के साथ चला सकते हैं ताकि टेस्ट प्रोफ़ाइल में सेटअप की गई डिफ़ॉल्ट वैल्यू का उपयोग किया जा सके और कमांड लाइन पर पथ लिखने से बचा जा सके।

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

तुम चेक कर सकते हो कि इंडेक्स फ़ाइल सही तरीके से जनरेट हुई है या नहीं, work डायरेक्टरी या results डायरेक्टरी में देखकर।

??? abstract "Work डायरेक्टरी कंटेंट"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Results डायरेक्टरी कंटेंट"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

यह रहा!

### सारांश

तुम जानते हो कि एक प्रोसेस वाला मॉड्यूल कैसे बनाया जाए, इसे workflow में इम्पोर्ट कैसे किया जाए, इनपुट चैनल के साथ इसे कॉल कैसे किया जाए, और परिणाम पब्लिश कैसे किए जाएं।

### आगे क्या है?

एक दूसरा चरण जोड़ो जो इंडेक्सिंग प्रोसेस के आउटपुट को लेता है और इसे वेरिएंट कॉलिंग चलाने के लिए उपयोग करता है।

---

## 2. इंडेक्स की गई BAM फ़ाइल पर GATK HaplotypeCaller चलाने के लिए एक दूसरी प्रोसेस जोड़ो

अब जब हमारे पास अपनी इनपुट फ़ाइल के लिए एक इंडेक्स है, तो हम वेरिएंट कॉलिंग चरण सेटअप करने के लिए आगे बढ़ सकते हैं।

[भाग 1](01_method.md) से `gatk HaplotypeCaller` कमांड को याद करो:

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

यह कमांड एक BAM फ़ाइल (`-I`), एक रेफ़रेंस जीनोम (`-R`), और एक intervals फ़ाइल (`-L`) लेता है, और एक VCF फ़ाइल (`-O`) उसके इंडेक्स के साथ बनाता है।
टूल BAM इंडेक्स, रेफ़रेंस इंडेक्स, और रेफ़रेंस डिक्शनरी को उनकी संबंधित फ़ाइलों के साथ कोलोकेट होने की उम्मीद करता है।
कंटेनर URI था `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`।

हम पहले की तरह ही तीन चरणों का पालन करते हैं:

1. इनपुट सेटअप करो
2. वेरिएंट कॉलिंग प्रोसेस लिखो और workflow में इसे कॉल करो
3. आउटपुट हैंडलिंग को कॉन्फ़िगर करो

### 2.1. इनपुट सेटअप करो

वेरिएंट कॉलिंग चरण को कई अतिरिक्त इनपुट फ़ाइलों की आवश्यकता होती है।
हमें उनके लिए पैरामीटर डिक्लेयर करने होंगे, टेस्ट प्रोफ़ाइल में डिफ़ॉल्ट वैल्यू जोड़नी होंगी, और उन्हें लोड करने के लिए वेरिएबल बनाने होंगे।

#### 2.1.1. एक्सेसरी इनपुट के लिए पैरामीटर डिक्लेरेशन जोड़ो

चूंकि हमारी नई प्रोसेस को कई अतिरिक्त फ़ाइलें प्रदान करने की उम्मीद है, `genomics.nf` में `Pipeline parameters` सेक्शन के तहत उनके लिए पैरामीटर डिक्लेरेशन जोड़ो:

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

पहले की तरह, हम inline के बजाय टेस्ट प्रोफ़ाइल के माध्यम से डिफ़ॉल्ट वैल्यू प्रदान करते हैं।

#### 2.1.2. टेस्ट प्रोफ़ाइल में एक्सेसरी फ़ाइल डिफ़ॉल्ट जोड़ो

जैसे हमने सेक्शन 1.1.2 में `reads_bam` के लिए किया था, `nextflow.config` में टेस्ट प्रोफ़ाइल में एक्सेसरी फ़ाइलों के लिए डिफ़ॉल्ट वैल्यू जोड़ो:

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

अब हमें वेरिएबल बनाने की आवश्यकता है जो workflow में उपयोग के लिए इन फ़ाइल पथों को लोड करें।

#### 2.1.3. एक्सेसरी फ़ाइलों के लिए वेरिएबल बनाओ

workflow ब्लॉक के अंदर एक्सेसरी फ़ाइल पथों के लिए वेरिएबल जोड़ो:

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

`file()` सिंटैक्स Nextflow को स्पष्ट रूप से इन इनपुट को फ़ाइल पथ के रूप में हैंडल करने के लिए कहता है।
तुम Side Quest [Working with files](../../side_quests/working_with_files.md) में इसके बारे में अधिक जान सकते हो।

### 2.2. वेरिएंट कॉलिंग प्रोसेस लिखो और workflow में इसे कॉल करो

हमें मॉड्यूल फ़ाइल में प्रोसेस डेफिनिशन लिखनी होगी, include स्टेटमेंट का उपयोग करके इसे workflow में इम्पोर्ट करना होगा, और इनपुट रीड्स प्लस इंडेक्सिंग चरण के आउटपुट और एक्सेसरी फ़ाइलों पर इसे कॉल करना होगा।

#### 2.2.1. वेरिएंट कॉलिंग प्रोसेस के लिए मॉड्यूल भरो

`modules/gatk_haplotypecaller.nf` खोलो और प्रोसेस डेफिनिशन की रूपरेखा जांचो।

आगे बढ़ो और ऊपर दी गई जानकारी का उपयोग करके खुद प्रोसेस डेफिनिशन भरो, फिर नीचे "बाद में" टैब में समाधान के विरुद्ध अपना काम चेक करो।

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

तुम देखोगे कि इस प्रोसेस में GATK कमांड की आवश्यकता से अधिक इनपुट हैं।
GATK नामकरण परंपराओं के आधार पर BAM इंडेक्स फ़ाइल और रेफ़रेंस जीनोम की एक्सेसरी फ़ाइलों को देखना जानता है, लेकिन Nextflow डोमेन-अज्ञेयवादी है और इन परंपराओं के बारे में नहीं जानता।
हमें उन्हें स्पष्ट रूप से सूचीबद्ध करने की आवश्यकता है ताकि Nextflow रनटाइम पर उन्हें working डायरेक्टरी में स्टेज करे; अन्यथा GATK लापता फ़ाइलों के बारे में एक त्रुटि फेंकेगा।

इसी तरह, हम आउटपुट VCF की इंडेक्स फ़ाइल (`"${input_bam}.vcf.idx"`) को स्पष्ट रूप से सूचीबद्ध करते हैं ताकि Nextflow बाद के चरणों के लिए इसका ट्रैक रखे।
हम प्रत्येक आउटपुट चैनल को एक नाम असाइन करने के लिए `emit:` सिंटैक्स का उपयोग करते हैं, जो तब उपयोगी होगा जब हम आउटपुट को publish ब्लॉक में वायर करेंगे।

एक बार जब तुम इसे पूरा कर लो, तो प्रोसेस पूर्ण हो जाती है।
इसे workflow में उपयोग करने के लिए, तुम्हें मॉड्यूल इम्पोर्ट करना होगा और एक प्रोसेस कॉल जोड़नी होगी।

#### 2.2.2. नए मॉड्यूल को इम्पोर्ट करो

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

अब प्रोसेस workflow स्कोप में उपलब्ध है।

#### 2.2.3. प्रोसेस कॉल जोड़ो

`main:` के तहत workflow body में प्रोसेस कॉल जोड़ो:

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

तुम्हें Hello Nextflow ट्रेनिंग सीरीज़ से `*.out` सिंटैक्स को पहचानना चाहिए; हम Nextflow को बता रहे हैं कि `SAMTOOLS_INDEX` द्वारा आउटपुट चैनल लें और उसे `GATK_HAPLOTYPECALLER` प्रोसेस कॉल में प्लग करें।

!!! note

    ध्यान दें कि प्रोसेस के कॉल में इनपुट उसी क्रम में प्रदान किए जाते हैं जिस क्रम में वे प्रोसेस के इनपुट ब्लॉक में सूचीबद्ध हैं।
    Nextflow में, इनपुट पोज़िशनल होते हैं, मतलब तुम्हें उसी क्रम का पालन करना _चाहिए_; और निश्चित रूप से तत्वों की संख्या समान होनी चाहिए।

### 2.3. आउटपुट हैंडलिंग को कॉन्फ़िगर करो

हमें नए आउटपुट को publish डिक्लेरेशन में जोड़ने और कॉन्फ़िगर करने की आवश्यकता है कि वे कहां जाते हैं।

#### 2.3.1. वेरिएंट कॉलिंग आउटपुट के लिए publish टार्गेट जोड़ो

`publish:` सेक्शन में VCF और इंडेक्स आउटपुट जोड़ो:

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

अब हमें Nextflow को यह बताना होगा कि नए आउटपुट को कहां रखना है।

#### 2.3.2. नए आउटपुट टार्गेट को कॉन्फ़िगर करो

`output {}` ब्लॉक में `vcf` और `vcf_idx` टार्गेट के लिए एंट्री जोड़ो, दोनों को `vcf/` सबडायरेक्टरी में पब्लिश करते हुए:

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

VCF और इसका इंडेक्स अलग-अलग टार्गेट के रूप में पब्लिश किए जाते हैं जो दोनों `vcf/` सबडायरेक्टरी में जाते हैं।

### 2.4. Workflow चलाओ

विस्तारित workflow चलाओ, इस बार `-resume` जोड़ते हुए ताकि हमें इंडेक्सिंग चरण फिर से चलाना न पड़े।

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

अब यदि हम कंसोल आउटपुट को देखें, तो हम दो प्रोसेस सूचीबद्ध देखते हैं।

पहली प्रोसेस को caching की बदौलत स्किप कर दिया गया, जैसी उम्मीद थी, जबकि दूसरी प्रोसेस को चलाया गया क्योंकि यह बिल्कुल नई है।

तुम results डायरेक्टरी में आउटपुट फ़ाइलें पाओगे (work डायरेक्टरी के सिंबॉलिक लिंक के रूप में)।

??? abstract "डायरेक्टरी कंटेंट"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

यदि तुम VCF फ़ाइल खोलते हो, तो तुम्हें वही कंटेंट देखना चाहिए जो उस फ़ाइल में था जिसे तुमने कंटेनर में सीधे GATK कमांड चलाकर जनरेट किया था।

??? abstract "फ़ाइल कंटेंट"

    ```console title="reads_mother.bam.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

यह वह आउटपुट है जिसे हम अपने अध्ययन में प्रत्येक नमूने के लिए जनरेट करना चाहते हैं।

### सारांश

तुम जानते हो कि दो-चरणीय मॉड्यूलर workflow कैसे बनाया जाए जो वास्तविक विश्लेषण कार्य करता है और genomics फ़ाइल फ़ॉर्मेट की विशेषताओं जैसे एक्सेसरी फ़ाइलों से निपटने में सक्षम है।

### आगे क्या है?

Workflow को थोक में कई नमूने हैंडल करने योग्य बनाओ।

---

## 3. Workflow को नमूनों के बैच पर चलाने के लिए अनुकूलित करो

एक workflow होना जो एक single नमूने पर प्रोसेसिंग को स्वचालित कर सके, सब ठीक है, लेकिन अगर तुम्हारे पास 1000 नमूने हैं तो क्या होगा?
क्या तुम्हें एक bash script लिखने की आवश्यकता है जो तुम्हारे सभी नमूनों के माध्यम से लूप करे?

नहीं, भगवान का शुक्र है! बस कोड में एक छोटा सा बदलाव करो और Nextflow तुम्हारे लिए भी इसे हैंडल करेगा।

### 3.1. तीन नमूनों को सूचीबद्ध करने के लिए इनपुट को अपडेट करो

कई नमूनों पर चलाने के लिए, टेस्ट प्रोफ़ाइल को एक के बजाय फ़ाइल पथों का एक ऐरे प्रदान करने के लिए अपडेट करो।
यह मल्टी-सैंपल निष्पादन का परीक्षण करने का एक त्वरित तरीका है; अगले चरण में हम इनपुट की एक फ़ाइल का उपयोग करके अधिक स्केलेबल दृष्टिकोण पर स्विच करेंगे।

सबसे पहले, पैरामीटर डिक्लेरेशन में टाइप एनोटेशन को कमेंट करो, क्योंकि ऐरे typed डिक्लेरेशन का उपयोग नहीं कर सकते:

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

फिर तीनों नमूनों को सूचीबद्ध करने के लिए टेस्ट प्रोफ़ाइल को अपडेट करो:

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

Workflow body में चैनल फ़ैक्ट्री (`.fromPath`) एकल की तरह ही कई फ़ाइल पथों को स्वीकार करती है, इसलिए कोई अन्य परिवर्तन की आवश्यकता नहीं है।

### 3.2. Workflow चलाओ

अब workflow चलाने की कोशिश करो जब प्लंबिंग सभी तीन टेस्ट नमूनों पर चलाने के लिए सेटअप है।

```bash
nextflow run genomics.nf -profile test -resume
```

मज़ेदार बात: यह _काम कर सकता है_, या यह _फेल हो सकता है_। उदाहरण के लिए, यहां एक रन है जो सफल रहा:

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

यदि तुम्हारा workflow रन सफल रहा, तो इसे तब तक फिर से चलाओ जब तक तुम्हें इस तरह की त्रुटि न मिले:

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

यदि तुम GATK कमांड त्रुटि आउटपुट को देखते हो, तो इस तरह की एक लाइन होगी:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

खैर, यह अजीब है, यह देखते हुए कि हमने workflow के पहले चरण में BAM फ़ाइलों को स्पष्ट रूप से इंडेक्स किया। क्या प्लंबिंग में कुछ गड़बड़ हो सकती है?

### 3.3. समस्या को ट्रबलशूट करो

हम work डायरेक्टरी का निरीक्षण करेंगे और यह पता लगाने के लिए `view()` ऑपरेटर का उपयोग करेंगे कि क्या गलत हुआ।

#### 3.3.1. संबंधित कॉल के लिए work डायरेक्टरी चेक करो

कंसोल आउटपुट में सूचीबद्ध फेल `GATK_HAPLOTYPECALLER` प्रोसेस कॉल के लिए work डायरेक्टरी के अंदर देखो।

??? abstract "डायरेक्टरी कंटेंट"

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

विशेष रूप से BAM फ़ाइल और BAM इंडेक्स के नामों पर ध्यान दो जो इस डायरेक्टरी में सूचीबद्ध हैं: `reads_son.bam` और `reads_father.bam.bai`।

क्या हुआ? Nextflow ने इस प्रोसेस कॉल की work डायरेक्टरी में एक इंडेक्स फ़ाइल स्टेज की है, लेकिन यह गलत है। यह कैसे हो सकता है?

#### 3.3.2. चैनल कंटेंट का निरीक्षण करने के लिए [view() ऑपरेटर](https://www.nextflow.io/docs/latest/reference/operator.html#view) का उपयोग करो

`GATK_HAPLOTYPECALLER` प्रोसेस कॉल से पहले workflow body में चैनल की कंटेंट देखने के लिए ये दो लाइनें जोड़ो:

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

फिर workflow कमांड फिर से चलाओ।

```bash
nextflow run genomics.nf -profile test
```

फिर से, यह सफल हो सकता है या फेल हो सकता है। यहां एक फेल रन के लिए दो `.view()` कॉल के आउटपुट जैसे दिखते हैं:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

पहली तीन लाइनें इनपुट चैनल से मेल खाती हैं और दूसरी, आउटपुट चैनल से।
तुम देख सकते हो कि तीन नमूनों के लिए BAM फ़ाइलें और इंडेक्स फ़ाइलें एक ही क्रम में सूचीबद्ध नहीं हैं!

!!! note

    जब तुम कई तत्वों वाले चैनल पर एक Nextflow प्रोसेस को कॉल करते हो, तो Nextflow निष्पादन को यथासंभव समानांतर करने की कोशिश करेगा, और जो भी क्रम में उपलब्ध होंगे उस क्रम में आउटपुट एकत्र करेगा।
    परिणाम यह है कि संबंधित आउटपुट मूल इनपुट की तुलना में एक अलग क्रम में एकत्र किए जा सकते हैं।

वर्तमान में लिखे गए अनुसार, हमारा workflow script मानता है कि इंडेक्स फ़ाइलें इंडेक्सिंग चरण से उसी mother/father/son क्रम में सूचीबद्ध होंगी जिस क्रम में इनपुट दिए गए थे।
लेकिन ऐसा होने की गारंटी नहीं है, यही कारण है कि कभी-कभी (हालांकि हमेशा नहीं) गलत फ़ाइलें दूसरे चरण में जोड़ी बन जाती हैं।

इसे ठीक करने के लिए, हमें यह सुनिश्चित करने की आवश्यकता है कि BAM फ़ाइलें और उनकी इंडेक्स फ़ाइलें चैनलों के माध्यम से एक साथ यात्रा करें।

!!! tip

    Workflow कोड में `.view()` स्टेटमेंट कुछ नहीं करते, इसलिए उन्हें छोड़ना कोई समस्या नहीं है।
    हालाँकि वे तुम्हारे कंसोल आउटपुट को अव्यवस्थित करेंगे, इसलिए हम समस्या को ट्रबलशूट करने के बाद उन्हें हटाने की सलाह देते हैं।

### 3.4. इंडेक्स फ़ाइलों को सही ढंग से हैंडल करने के लिए workflow को अपडेट करो

फिक्स यह है कि प्रत्येक BAM फ़ाइल को उसके इंडेक्स के साथ एक tuple में बंडल करें, फिर downstream प्रोसेस और workflow प्लंबिंग को मिलान के लिए अपडेट करें।

#### 3.4.1. SAMTOOLS_INDEX मॉड्यूल के आउटपुट को tuple में बदलो

यह सुनिश्चित करने का सबसे सरल तरीका कि BAM फ़ाइल और उसका इंडेक्स निकटता से जुड़े रहें, उन्हें इंडेक्स टास्क से बाहर आने वाले tuple में पैकेज करना है।

!!! note

    एक **tuple** तत्वों की एक सीमित, क्रमबद्ध सूची है जिसका उपयोग आमतौर पर एक फ़ंक्शन से कई वैल्यू लौटाने के लिए किया जाता है। Tuples उनके संबंध और क्रम को संरक्षित करते हुए प्रोसेस के बीच कई इनपुट या आउटपुट पास करने के लिए विशेष रूप से उपयोगी हैं।

BAM फ़ाइल को शामिल करने के लिए `modules/samtools_index.nf` में आउटपुट को अपडेट करो:

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

इस तरह, प्रत्येक इंडेक्स फ़ाइल को उसकी मूल BAM फ़ाइल के साथ कसकर जोड़ा जाएगा, और इंडेक्सिंग चरण का समग्र आउटपुट फ़ाइलों के जोड़े वाला एकल चैनल होगा।

#### 3.4.2. GATK_HAPLOTYPECALLER मॉड्यूल के इनपुट को tuple स्वीकार करने के लिए बदलो

चूंकि हमने पहली प्रोसेस के आउटपुट का 'शेप' बदल दिया है, हमें मिलान के लिए दूसरी प्रोसेस की इनपुट डेफिनिशन को अपडेट करने की आवश्यकता है।

`modules/gatk_haplotypecaller.nf` को अपडेट करो:

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

अब हमें प्रोसेस कॉल और publish टार्गेट में नए tuple संरचना को प्रतिबिंबित करने के लिए workflow को अपडेट करने की आवश्यकता है।

#### 3.4.3. Workflow में GATK_HAPLOTYPECALLER के कॉल को अपडेट करो

हमें अब `GATK_HAPLOTYPECALLER` प्रोसेस को मूल `reads_ch` प्रदान करने की आवश्यकता नहीं है, क्योंकि BAM फ़ाइल अब `SAMTOOLS_INDEX` द्वारा आउटपुट चैनल में बंडल है।

`genomics.nf` में कॉल को अपडेट करो:

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

अंत में, हमें नई आउटपुट संरचना को प्रतिबिंबित करने के लिए publish टार्गेट को अपडेट करने की आवश्यकता है।

#### 3.4.4. इंडेक्स की गई BAM आउटपुट के लिए publish टार्गेट को अपडेट करो

चूंकि SAMTOOLS_INDEX आउटपुट अब एक tuple है जिसमें BAM फ़ाइल और उसका इंडेक्स दोनों हैं, publish टार्गेट को `bam_index` से `indexed_bam` में बदल दो ताकि इसकी सामग्री को बेहतर ढंग से प्रतिबिंबित किया जा सके:

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

इन परिवर्तनों के साथ, BAM और उसका इंडेक्स एक साथ यात्रा करने की गारंटी है, इसलिए जोड़ी हमेशा सही होगी।

### 3.5. सही की गई workflow चलाओ

यह सुनिश्चित करने के लिए workflow को फिर से चलाओ कि यह आगे बढ़ने पर विश्वसनीय रूप से काम करेगा।

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

Results डायरेक्टरी में अब प्रत्येक नमूने के लिए BAM और BAI फ़ाइलें (tuple से) हैं, VCF आउटपुट के साथ:

??? abstract "Results डायरेक्टरी कंटेंट"

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

संबंधित फ़ाइलों को tuples में बंडल करके, हमने सुनिश्चित किया कि सही फ़ाइलें हमेशा workflow के माध्यम से एक साथ यात्रा करें।
Workflow अब किसी भी संख्या में नमूनों को विश्वसनीय रूप से प्रोसेस करता है, लेकिन उन्हें config में व्यक्तिगत रूप से सूचीबद्ध करना बहुत स्केलेबल नहीं है।
अगले चरण में, हम एक फ़ाइल से इनपुट पढ़ने पर स्विच करेंगे।

### सारांश

तुम जानते हो कि अपने workflow को कई नमूनों पर (स्वतंत्र रूप से) चलाने के लिए कैसे बनाया जाए।

### आगे क्या है?

थोक में नमूनों को हैंडल करना आसान बनाओ।

---

## 4. Workflow को एक टेक्स्ट फ़ाइल स्वीकार करने योग्य बनाओ जिसमें इनपुट फ़ाइलों का बैच हो

एक workflow को कई डेटा इनपुट फ़ाइलें प्रदान करने का एक बहुत सामान्य तरीका इसे फ़ाइल पथों वाली एक टेक्स्ट फ़ाइल के साथ करना है।
यह एक टेक्स्ट फ़ाइल जितना सरल हो सकता है जिसमें प्रति पंक्ति एक फ़ाइल पथ हो और कुछ नहीं, या फ़ाइल में अतिरिक्त मेटाडेटा हो सकता है, इस स्थिति में इसे अक्सर samplesheet कहा जाता है।

यहां हम तुम्हें दिखाने जा रहे हैं कि सरल केस कैसे करना है।

### 4.1. इनपुट फ़ाइल पथों को सूचीबद्ध करने वाली प्रदान की गई टेक्स्ट फ़ाइल जांचो

हमने पहले से ही इनपुट फ़ाइल पथों को सूचीबद्ध करने वाली एक टेक्स्ट फ़ाइल बनाई, जिसे `sample_bams.txt` कहा जाता है, जिसे तुम `data/` डायरेक्टरी में पा सकते हो।

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

जैसा कि तुम देख सकते हो, हमने प्रति पंक्ति एक फ़ाइल पथ सूचीबद्ध किया, और वे एब्सोल्यूट पथ हैं।

!!! note

    हम यहां जो फ़ाइलें उपयोग कर रहे हैं वे बस तुम्हारे GitHub Codespaces के लोकल फ़ाइल सिस्टम पर हैं, लेकिन हम क्लाउड स्टोरेज में फ़ाइलों की ओर भी इशारा कर सकते हैं।
    यदि तुम प्रदान किए गए Codespaces वातावरण का उपयोग नहीं कर रहे हो, तो तुम्हें अपने लोकल सेटअप से मिलान करने के लिए फ़ाइल पथों को अनुकूलित करने की आवश्यकता हो सकती है।

### 4.2. पैरामीटर और टेस्ट प्रोफ़ाइल को अपडेट करो

व्यक्तिगत नमूनों को सूचीबद्ध करने के बजाय `reads_bam` पैरामीटर को `sample_bams.txt` फ़ाइल की ओर इशारा करने के लिए स्विच करो।

params ब्लॉक में टाइप एनोटेशन को रिस्टोर करो (क्योंकि यह फिर से एकल पथ है):

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

फिर टेस्ट प्रोफ़ाइल को टेक्स्ट फ़ाइल की ओर इशारा करने के लिए अपडेट करो:

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

### 4.3. फ़ाइल से लाइनें पढ़ने के लिए चैनल फ़ैक्ट्री को अपडेट करो

वर्तमान में, हमारी इनपुट चैनल फ़ैक्ट्री हम जो भी फ़ाइलें देते हैं उन्हें डेटा इनपुट के रूप में ट्रीट करती है जिन्हें हम इंडेक्सिंग प्रोसेस को फीड करना चाहते हैं।
चूंकि हम अब इसे एक फ़ाइल दे रहे हैं जो इनपुट फ़ाइल पथों को सूचीबद्ध करती है, हमें फ़ाइल को पार्स करने और इसमें मौजूद फ़ाइल पथों को डेटा इनपुट के रूप में ट्रीट करने के लिए इसके व्यवहार को बदलने की आवश्यकता है।

हम उसी पैटर्न का उपयोग करके यह कर सकते हैं जिसका उपयोग हमने [Hello Nextflow के भाग 2](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file) में किया था: फ़ाइल को पार्स करने के लिए [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) ऑपरेटर को लागू करना, फिर प्रत्येक पंक्ति के पहले फ़ील्ड को चुनने के लिए एक `map` ऑपरेशन।

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

तकनीकी रूप से हम इसे [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) ऑपरेटर का उपयोग करके अधिक सरलता से कर सकते हैं, क्योंकि हमारी इनपुट फ़ाइल वर्तमान में केवल फ़ाइल पथ रखती है।
हालाँकि, अधिक बहुमुखी `splitCsv` ऑपरेटर का उपयोग करके (`map` द्वारा पूरक), हम अपने workflow को भविष्य के लिए तैयार कर सकते हैं यदि हम फ़ाइल पथ वाली फ़ाइल में मेटाडेटा जोड़ने का निर्णय लेते हैं।

!!! tip

    यदि तुम्हें भरोसा नहीं है कि तुम समझते हो कि ऑपरेटर यहां क्या कर रहे हैं, तो यह उन्हें लागू करने से पहले और बाद में चैनल कंटेंट कैसी दिखती है यह देखने के लिए `.view()` ऑपरेटर का उपयोग करने का एक और शानदार अवसर है।

### 4.4. Workflow चलाओ

Workflow को एक बार और चलाओ। इसे पहले जैसा ही परिणाम देना चाहिए, सही?

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

हाँ! वास्तव में, Nextflow सही तरीके से पता लगाता है कि प्रोसेस कॉल बिल्कुल समान हैं, और सब कुछ फिर से चलाने की ज़रूरत भी नहीं है, क्योंकि हम `-resume` के साथ चला रहे थे।

और बस इतना ही! हमारे सरल वेरिएंट कॉलिंग workflow में वे सभी बुनियादी सुविधाएँ हैं जो हम चाहते थे।

### सारांश

तुम जानते हो कि BAM फ़ाइल को इंडेक्स करने और GATK का उपयोग करके प्रति-नमूना वेरिएंट कॉलिंग लागू करने के लिए मल्टी-स्टेप मॉड्यूलर workflow कैसे बनाया जाए।

अधिक सामान्य रूप से, तुमने सीखा है कि वास्तविक काम करने वाली एक सरल genomics पाइपलाइन बनाने के लिए आवश्यक Nextflow घटकों और तर्क का उपयोग कैसे करें, genomics फ़ाइल फ़ॉर्मेट की विशेषताओं और टूल आवश्यकताओं को ध्यान में रखते हुए।

### आगे क्या है?

अपनी सफलता का जश्न मनाओ और एक अतिरिक्त लंबा ब्रेक लो!

इस कोर्स के अगले भाग में, तुम सीखोगे कि डेटा पर joint वेरिएंट कॉलिंग लागू करने के लिए इस सरल प्रति-नमूना वेरिएंट कॉलिंग workflow को कैसे रूपांतरित किया जाए।
