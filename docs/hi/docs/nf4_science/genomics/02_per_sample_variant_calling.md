# भाग 2: प्रति-नमूना वेरिएंट कॉलिंग

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

भाग 1 में, तुमने Samtools और GATK कमांड को उनके संबंधित कंटेनर में मैन्युअल रूप से टेस्ट किया।
अब हम उन्हीं कमांड को एक Nextflow workflow में wrap करने जा रहे हैं।

## असाइनमेंट

कोर्स के इस भाग में, हम एक workflow विकसित करने जा रहे हैं जो निम्नलिखित करता है:

1. [Samtools](https://www.htslib.org/) का उपयोग करके प्रत्येक BAM इनपुट फ़ाइल के लिए एक इंडेक्स फ़ाइल जेनरेट करना
2. प्रत्येक BAM इनपुट फ़ाइल पर GATK HaplotypeCaller चलाना ताकि VCF (Variant Call Format) में प्रति-नमूना वेरिएंट कॉल जेनरेट हो सकें

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

यह भाग 1 के चरणों को दोहराता है, जहाँ तुमने इन कमांड को उनके कंटेनर में मैन्युअल रूप से चलाया था।

शुरुआती बिंदु के रूप में, हम तुम्हें एक workflow फ़ाइल, `genomics.nf`, प्रदान करते हैं, जो workflow के मुख्य भागों की रूपरेखा देती है, साथ ही दो मॉड्यूल फ़ाइलें, samtools_index.nf और gatk_haplotypecaller.nf, जो मॉड्यूल की संरचना की रूपरेखा देती हैं।
ये फ़ाइलें कार्यात्मक नहीं हैं; उनका उद्देश्य केवल scaffolds के रूप में काम करना है जिन्हें तुम कोड के दिलचस्प भागों से भर सकते हो।

## पाठ योजना

विकास प्रक्रिया को अधिक शैक्षिक बनाने के लिए, हमने इसे चार चरणों में विभाजित किया है:

1. **एक single-stage workflow लिखो जो BAM फ़ाइल पर Samtools index चलाता है।**
   यह एक मॉड्यूल बनाना, उसे import करना, और workflow में उसे कॉल करना कवर करता है।
2. **indexed BAM फ़ाइल पर GATK HaplotypeCaller चलाने के लिए दूसरा प्रोसेस जोड़ो।**
   यह प्रोसेस आउटपुट को इनपुट से जोड़ना और accessory फ़ाइलों को हैंडल करना पेश करता है।
3. **workflow को नमूनों के batch पर चलाने के लिए अनुकूलित करो।**
   यह parallel execution को कवर करता है और संबंधित फ़ाइलों को एक साथ रखने के लिए tuples पेश करता है।
4. **workflow को एक टेक्स्ट फ़ाइल स्वीकार करने योग्य बनाओ जिसमें इनपुट फ़ाइलों का batch हो।**
   यह bulk में इनपुट प्रदान करने के लिए एक सामान्य पैटर्न प्रदर्शित करता है।

प्रत्येक चरण workflow विकास के एक विशिष्ट पहलू पर केंद्रित है।

---

## 1. एक single-stage workflow लिखो जो BAM फ़ाइल पर Samtools index चलाता है

यह पहला चरण मूल बातों पर केंद्रित है: एक BAM फ़ाइल लोड करना और उसके लिए एक इंडेक्स जेनरेट करना।

[भाग 1](01_method.md) से `samtools index` कमांड याद करो:

```bash
samtools index '<input_bam>'
```

कमांड इनपुट के रूप में एक BAM फ़ाइल लेता है और उसके साथ एक `.bai` इंडेक्स फ़ाइल उत्पन्न करता है।
कंटेनर URI था `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464`।

हम इस जानकारी को लेने जा रहे हैं और इसे तीन चरणों में Nextflow में wrap करने जा रहे हैं:

1. इनपुट सेट अप करो
2. indexing प्रोसेस लिखो और workflow में उसे कॉल करो
3. आउटपुट हैंडलिंग कॉन्फ़िगर करो

### 1.1. इनपुट सेट अप करो

हमें एक इनपुट पैरामीटर घोषित करना होगा, एक सुविधाजनक डिफ़ॉल्ट मान प्रदान करने के लिए एक test profile बनाना होगा, और एक इनपुट चैनल बनाना होगा।

#### 1.1.1. एक इनपुट पैरामीटर घोषणा जोड़ो

मुख्य workflow फ़ाइल `genomics.nf` में, `Pipeline parameters` सेक्शन के अंतर्गत, `reads_bam` नामक एक CLI पैरामीटर घोषित करो।

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

यह CLI पैरामीटर सेट अप करता है, लेकिन हम विकास के दौरान हर बार workflow चलाते समय फ़ाइल पथ टाइप नहीं करना चाहते।
डिफ़ॉल्ट मान प्रदान करने के लिए कई विकल्प हैं; यहाँ हम एक test profile का उपयोग करते हैं।

#### 1.1.2. `nextflow.config` में डिफ़ॉल्ट मान के साथ एक test profile बनाओ

एक test profile कमांड लाइन पर इनपुट निर्दिष्ट किए बिना workflow को आज़माने के लिए सुविधाजनक डिफ़ॉल्ट मान प्रदान करता है।
यह Nextflow ecosystem में एक सामान्य परंपरा है (अधिक विवरण के लिए [Hello Config](../../hello_nextflow/06_hello_config.md) देखो)।

`nextflow.config` में एक `profiles` ब्लॉक जोड़ो जिसमें एक `test` profile हो जो `reads_bam` पैरामीटर को test BAM फ़ाइलों में से एक पर सेट करता है।

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

यहाँ, हम `${projectDir}` का उपयोग कर रहे हैं, एक built-in Nextflow वेरिएबल जो उस डायरेक्टरी की ओर इशारा करता है जहाँ workflow स्क्रिप्ट स्थित है।
यह absolute paths को hardcode किए बिना डेटा फ़ाइलों और अन्य संसाधनों को संदर्भित करना आसान बनाता है।

#### 1.1.3. इनपुट चैनल सेट अप करो

workflow ब्लॉक में, `.fromPath` चैनल फ़ैक्टरी का उपयोग करके पैरामीटर मान से एक इनपुट चैनल बनाओ ([Hello Channels](../../hello_nextflow/02_hello_channels.md) में उपयोग किए गए अनुसार)।

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

अब हमें इस इनपुट पर indexing चलाने के लिए प्रोसेस बनाना होगा।

### 1.2. indexing प्रोसेस लिखो और workflow में उसे कॉल करो

हमें मॉड्यूल फ़ाइल में प्रोसेस परिभाषा लिखनी होगी, include statement का उपयोग करके इसे workflow में import करना होगा, और इनपुट पर इसे कॉल करना होगा।

#### 1.2.1. indexing प्रोसेस के लिए मॉड्यूल भरो

`modules/samtools_index.nf` खोलो और प्रोसेस परिभाषा की रूपरेखा की जाँच करो।
तुम्हें मुख्य संरचनात्मक तत्वों को पहचानना चाहिए; यदि नहीं, तो refresher के लिए [Hello Nextflow](../../hello_nextflow/01_hello_world.md) पढ़ने पर विचार करो।

आगे बढ़ो और ऊपर दी गई जानकारी का उपयोग करके प्रोसेस परिभाषा को स्वयं भरो, फिर नीचे "बाद में" टैब में समाधान के विरुद्ध अपने काम की जाँच करो।

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
इसे workflow में उपयोग करने के लिए, तुम्हें मॉड्यूल import करना होगा और एक प्रोसेस कॉल जोड़नी होगी।

#### 1.2.2. मॉड्यूल include करो

`genomics.nf` में, प्रोसेस को workflow के लिए उपलब्ध कराने के लिए एक `include` statement जोड़ो:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    ```

प्रोसेस अब workflow scope में उपलब्ध है।

#### 1.2.3. इनपुट पर indexing प्रोसेस कॉल करो

अब, workflow ब्लॉक में `SAMTOOLS_INDEX` के लिए एक कॉल जोड़ें, इनपुट चैनल को argument के रूप में पास करते हुए।

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

workflow अब इनपुट लोड करता है और उस पर indexing प्रोसेस चलाता है।
इसके बाद, हमें कॉन्फ़िगर करना होगा कि आउटपुट कैसे publish किया जाए।

### 1.3. आउटपुट हैंडलिंग कॉन्फ़िगर करो

हमें यह घोषित करना होगा कि कौन से प्रोसेस आउटपुट publish करने हैं और निर्दिष्ट करना होगा कि वे कहाँ जाने चाहिए।

#### 1.3.1. `publish:` सेक्शन में एक आउटपुट घोषित करो

workflow ब्लॉक के अंदर `publish:` सेक्शन घोषित करता है कि कौन से प्रोसेस आउटपुट publish किए जाने चाहिए।
`SAMTOOLS_INDEX` के आउटपुट को `bam_index` नामक एक named target को assign करो।

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

अब हमें Nextflow को बताना होगा कि published आउटपुट कहाँ रखना है।

#### 1.3.2. `output {}` ब्लॉक में आउटपुट target कॉन्फ़िगर करो

`output {}` ब्लॉक workflow के बाहर बैठता है और निर्दिष्ट करता है कि प्रत्येक named target कहाँ publish किया जाता है।
चलो `bam_index` के लिए एक target जोड़ें जो `bam/` subdirectory में publish करता है।

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

!!! note "नोट"

    डिफ़ॉल्ट रूप से, Nextflow आउटपुट फ़ाइलों को symbolic links के रूप में publish करता है, जो अनावश्यक duplication से बचता है।
    भले ही हम यहाँ जो डेटा फ़ाइलें उपयोग कर रहे हैं वे बहुत छोटी हैं, genomics में वे बहुत बड़ी हो सकती हैं।
    जब तुम अपनी `work` डायरेक्टरी को साफ़ करते हो तो Symlinks टूट जाएंगे, इसलिए production workflows के लिए तुम डिफ़ॉल्ट publish mode को `'copy'` में override करना चाह सकते हो।

### 1.4. workflow चलाओ

इस बिंदु पर, हमारे पास एक one-step indexing workflow है जो पूरी तरह से कार्यात्मक होना चाहिए। चलो परीक्षण करें कि यह काम करता है!

हम इसे `-profile test` के साथ चला सकते हैं ताकि test profile में सेट किए गए डिफ़ॉल्ट मान का उपयोग किया जा सके और कमांड लाइन पर पथ लिखने से बचा जा सके।

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

तुम जाँच सकते हो कि इंडेक्स फ़ाइल सही तरीके से जेनरेट हुई है या नहीं work डायरेक्टरी या results डायरेक्टरी में देखकर।

??? abstract "Work डायरेक्टरी सामग्री"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Results डायरेक्टरी सामग्री"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

वहाँ है!

### सारांश

तुम जानते हो कि एक प्रोसेस युक्त मॉड्यूल कैसे बनाएँ, इसे workflow में import करें, इसे इनपुट चैनल के साथ कॉल करें, और परिणाम publish करें।

### आगे क्या है?

एक दूसरा चरण जोड़ो जो indexing प्रोसेस के आउटपुट को लेता है और इसका उपयोग वेरिएंट कॉलिंग चलाने के लिए करता है।

---

## 2. indexed BAM फ़ाइल पर GATK HaplotypeCaller चलाने के लिए दूसरा प्रोसेस जोड़ो

अब जब हमारे पास अपनी इनपुट फ़ाइल के लिए एक इंडेक्स है, तो हम वेरिएंट कॉलिंग चरण सेट अप करने के लिए आगे बढ़ सकते हैं।

[भाग 1](01_method.md) से `gatk HaplotypeCaller` कमांड याद करो:

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

कमांड एक BAM फ़ाइल (`-I`), एक reference genome (`-R`), और एक intervals फ़ाइल (`-L`) लेता है, और अपने इंडेक्स के साथ एक VCF फ़ाइल (`-O`) उत्पन्न करता है।
टूल BAM इंडेक्स, reference इंडेक्स, और reference dictionary को उनकी संबंधित फ़ाइलों के साथ co-located होने की भी अपेक्षा करता है।
कंटेनर URI था `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`।

हम पहले की तरह तीन चरणों का पालन करते हैं:

1. इनपुट सेट अप करो
2. वेरिएंट कॉलिंग प्रोसेस लिखो और workflow में उसे कॉल करो
3. आउटपुट हैंडलिंग कॉन्फ़िगर करो

### 2.1. इनपुट सेट अप करो

वेरिएंट कॉलिंग चरण के लिए कई अतिरिक्त इनपुट फ़ाइलों की आवश्यकता होती है।
हमें उनके लिए पैरामीटर घोषित करने होंगे, test profile में डिफ़ॉल्ट मान जोड़ने होंगे, और उन्हें लोड करने के लिए वेरिएबल बनाने होंगे।

#### 2.1.1. accessory इनपुट के लिए पैरामीटर घोषणाएँ जोड़ो

चूंकि हमारा नया प्रोसेस कुछ अतिरिक्त फ़ाइलों की अपेक्षा करता है, `genomics.nf` में `Pipeline parameters` सेक्शन के अंतर्गत उनके लिए पैरामीटर घोषणाएँ जोड़ो:

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

पहले की तरह, हम inline के बजाय test profile के माध्यम से डिफ़ॉल्ट मान प्रदान करते हैं।

#### 2.1.2. test profile में accessory फ़ाइल डिफ़ॉल्ट जोड़ो

जैसा कि हमने सेक्शन 1.1.2 में `reads_bam` के लिए किया था, `nextflow.config` में test profile में accessory फ़ाइलों के लिए डिफ़ॉल्ट मान जोड़ो:

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

अब हमें वेरिएबल बनाने होंगे जो workflow में उपयोग के लिए इन फ़ाइल पथों को लोड करते हैं।

#### 2.1.3. accessory फ़ाइलों के लिए वेरिएबल बनाओ

workflow ब्लॉक के अंदर accessory फ़ाइल पथों के लिए वेरिएबल जोड़ो:

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

`file()` सिंटैक्स Nextflow को स्पष्ट रूप से बताता है कि इन इनपुट को फ़ाइल पथों के रूप में हैंडल करना है।
तुम Side Quest [Working with files](../../side_quests/working_with_files.md) में इसके बारे में अधिक जान सकते हो।

### 2.2. वेरिएंट कॉलिंग प्रोसेस लिखो और workflow में उसे कॉल करो

हमें मॉड्यूल फ़ाइल में प्रोसेस परिभाषा लिखनी होगी, include statement का उपयोग करके इसे workflow में import करना होगा, और इनपुट reads plus indexing चरण के आउटपुट और accessory फ़ाइलों पर इसे कॉल करना होगा।

#### 2.2.1. वेरिएंट कॉलिंग प्रोसेस के लिए मॉड्यूल भरो

`modules/gatk_haplotypecaller.nf` खोलो और प्रोसेस परिभाषा की रूपरेखा की जाँच करो।

आगे बढ़ो और ऊपर दी गई जानकारी का उपयोग करके प्रोसेस परिभाषा को स्वयं भरो, फिर नीचे "बाद में" टैब में समाधान के विरुद्ध अपने काम की जाँच करो।

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
GATK naming conventions के आधार पर BAM इंडेक्स फ़ाइल और reference genome की accessory फ़ाइलों को देखना जानता है, लेकिन Nextflow domain-agnostic है और इन conventions के बारे में नहीं जानता।
हमें उन्हें स्पष्ट रूप से सूचीबद्ध करना होगा ताकि Nextflow उन्हें runtime पर working डायरेक्टरी में stage करे; अन्यथा GATK missing फ़ाइलों के बारे में एक error throw करेगा।

इसी तरह, हम आउटपुट VCF की इंडेक्स फ़ाइल (`"${input_bam}.vcf.idx"`) को स्पष्ट रूप से सूचीबद्ध करते हैं ताकि Nextflow बाद के चरणों के लिए इसका ट्रैक रखे।
हम प्रत्येक आउटपुट चैनल को एक नाम assign करने के लिए `emit:` सिंटैक्स का उपयोग करते हैं, जो तब उपयोगी हो जाएगा जब हम आउटपुट को publish ब्लॉक में wire करेंगे।

एक बार जब तुम इसे पूरा कर लो, तो प्रोसेस पूर्ण हो जाता है।
इसे workflow में उपयोग करने के लिए, तुम्हें मॉड्यूल import करना होगा और एक प्रोसेस कॉल जोड़नी होगी।

#### 2.2.2. नया मॉड्यूल import करो

नए मॉड्यूल को import करने के लिए `genomics.nf` को अपडेट करो:

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

प्रोसेस अब workflow scope में उपलब्ध है।

#### 2.2.3. प्रोसेस कॉल जोड़ो

workflow body में, `main:` के अंतर्गत प्रोसेस कॉल जोड़ो:

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

तुम्हें Hello Nextflow प्रशिक्षण श्रृंखला से `*.out` सिंटैक्स को पहचानना चाहिए; हम Nextflow को बता रहे हैं कि `SAMTOOLS_INDEX` द्वारा आउटपुट किए गए चैनल को लें और उसे `GATK_HAPLOTYPECALLER` प्रोसेस कॉल में plug करें।

!!! note "नोट"

    ध्यान दें कि इनपुट प्रोसेस के कॉल में उसी क्रम में प्रदान किए जाते हैं जिस क्रम में वे प्रोसेस के input ब्लॉक में सूचीबद्ध हैं।
    Nextflow में, इनपुट positional हैं, जिसका अर्थ है कि तुम्हें उसी क्रम का पालन _करना चाहिए_; और निश्चित रूप से तत्वों की संख्या समान होनी चाहिए।

### 2.3. आउटपुट हैंडलिंग कॉन्फ़िगर करो

हमें publish घोषणा में नए आउटपुट जोड़ने होंगे और कॉन्फ़िगर करना होगा कि वे कहाँ जाते हैं।

#### 2.3.1. वेरिएंट कॉलिंग आउटपुट के लिए publish targets जोड़ो

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

अब हमें Nextflow को बताना होगा कि नए आउटपुट कहाँ रखने हैं।

#### 2.3.2. नए आउटपुट targets कॉन्फ़िगर करो

`output {}` ब्लॉक में `vcf` और `vcf_idx` targets के लिए entries जोड़ो, दोनों को `vcf/` subdirectory में publish करते हुए:

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

VCF और इसका इंडेक्स अलग targets के रूप में publish किए जाते हैं जो दोनों `vcf/` subdirectory में जाते हैं।

### 2.4. workflow चलाओ

विस्तारित workflow चलाओ, इस बार `-resume` जोड़ते हुए ताकि हमें indexing चरण को फिर से चलाना न पड़े।

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

अब यदि हम console आउटपुट को देखें, तो हम दो प्रोसेस सूचीबद्ध देखते हैं।

पहला प्रोसेस caching के कारण skip किया गया था, जैसा कि अपेक्षित था, जबकि दूसरा प्रोसेस चलाया गया क्योंकि यह बिल्कुल नया है।

तुम results डायरेक्टरी में आउटपुट फ़ाइलें पाओगे (work डायरेक्टरी के symbolic links के रूप में)।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

यदि तुम VCF फ़ाइल खोलते हो, तो तुम्हें उस फ़ाइल के समान सामग्री दिखनी चाहिए जो तुमने कंटेनर में सीधे GATK कमांड चलाकर जेनरेट की थी।

??? abstract "फ़ाइल सामग्री"

    ```console title="reads_mother.bam.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

यह वह आउटपुट है जिसे हम अपने अध्ययन में प्रत्येक नमूने के लिए जेनरेट करने की परवाह करते हैं।

### सारांश

तुम जानते हो कि दो-चरणीय modular workflow कैसे बनाएँ जो वास्तविक विश्लेषण कार्य करता है और genomics फ़ाइल format की विशिष्टताओं जैसे accessory फ़ाइलों से निपटने में सक्षम है।

### आगे क्या है?

workflow को bulk में कई नमूनों को हैंडल करने योग्य बनाओ।

---

## 3. workflow को नमूनों के batch पर चलाने के लिए अनुकूलित करो

एक workflow होना जो एक single नमूने पर processing को automate कर सकता है, यह सब अच्छा है, लेकिन यदि तुम्हारे पास 1000 नमूने हैं तो क्या होगा?
क्या तुम्हें एक bash स्क्रिप्ट लिखने की आवश्यकता है जो तुम्हारे सभी नमूनों के माध्यम से loop करती है?

नहीं, भगवान का शुक्र है! बस कोड में एक minor tweak करो और Nextflow तुम्हारे लिए भी इसे हैंडल करेगा।

### 3.1. तीन नमूनों को सूचीबद्ध करने के लिए इनपुट अपडेट करो

कई नमूनों पर चलाने के लिए, test profile को एक single के बजाय फ़ाइल पथों की एक array प्रदान करने के लिए अपडेट करो।
यह multi-sample execution का परीक्षण करने का एक त्वरित तरीका है; अगले चरण में हम इनपुट की एक फ़ाइल का उपयोग करके अधिक scalable दृष्टिकोण पर स्विच करेंगे।

पहले, पैरामीटर घोषणा में type annotation को comment out करो, क्योंकि arrays typed घोषणाओं का उपयोग नहीं कर सकते:

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

फिर सभी तीन नमूनों को सूचीबद्ध करने के लिए test profile को अपडेट करो:

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

workflow body में चैनल फ़ैक्टरी (`.fromPath`) एक single के साथ-साथ कई फ़ाइल पथों को स्वीकार करती है, इसलिए कोई अन्य परिवर्तन की आवश्यकता नहीं है।

### 3.2. workflow चलाओ

अब workflow चलाने की कोशिश करो जब plumbing सभी तीन test नमूनों पर चलाने के लिए सेट अप है।

```bash
nextflow run genomics.nf -profile test -resume
```

मज़ेदार बात: यह _काम कर सकता है_, या यह _fail हो सकता है_। उदाहरण के लिए, यहाँ एक run है जो सफल रहा:

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

यदि तुम्हारा workflow run सफल रहा, तो इसे तब तक फिर से चलाओ जब तक तुम्हें इस तरह की error न मिले:

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

यदि तुम GATK कमांड error आउटपुट को देखते हो, तो इस तरह की एक लाइन होगी:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

खैर, यह अजीब है, यह देखते हुए कि हमने workflow के पहले चरण में स्पष्ट रूप से BAM फ़ाइलों को indexed किया। क्या plumbing में कुछ गड़बड़ हो सकती है?

### 3.3. समस्या का निवारण करो

हम work डायरेक्टरी का निरीक्षण करेंगे और यह पता लगाने के लिए `view()` ऑपरेटर का उपयोग करेंगे कि क्या गलत हुआ।

#### 3.3.1. संबंधित कॉल के लिए work डायरेक्टरी की जाँच करो

console आउटपुट में सूचीबद्ध failed `GATK_HAPLOTYPECALLER` प्रोसेस कॉल के लिए work डायरेक्टरी के अंदर देखो।

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

इस डायरेक्टरी में सूचीबद्ध BAM फ़ाइल और BAM इंडेक्स के नामों पर विशेष ध्यान दो: `reads_son.bam` और `reads_father.bam.bai`।

क्या बात है? Nextflow ने इस प्रोसेस कॉल की work डायरेक्टरी में एक इंडेक्स फ़ाइल stage की है, लेकिन यह गलत है। यह कैसे हो सकता है?

#### 3.3.2. चैनल सामग्री का निरीक्षण करने के लिए [view() ऑपरेटर](https://www.nextflow.io/docs/latest/reference/operator.html#view) का उपयोग करो

चैनल की सामग्री को देखने के लिए `GATK_HAPLOTYPECALLER` प्रोसेस कॉल से पहले workflow body में ये दो लाइनें जोड़ो:

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

एक बार फिर, यह सफल या fail हो सकता है। यहाँ एक failed run के लिए दो `.view()` कॉल के आउटपुट की तरह दिखता है:

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

!!! note "नोट"

    जब तुम कई तत्वों वाले चैनल पर Nextflow प्रोसेस कॉल करते हो, तो Nextflow execution को यथासंभव parallelize करने की कोशिश करेगा, और जो भी क्रम में उपलब्ध हो जाएंगे उस क्रम में आउटपुट एकत्र करेगा।
    परिणाम यह है कि संबंधित आउटपुट मूल इनपुट के दिए गए क्रम से अलग क्रम में एकत्र किए जा सकते हैं।

जैसा कि वर्तमान में लिखा गया है, हमारी workflow स्क्रिप्ट मानती है कि इंडेक्स फ़ाइलें indexing चरण से उसी mother/father/son क्रम में सूचीबद्ध होंगी जैसे इनपुट दिए गए थे।
लेकिन यह गारंटी नहीं है कि ऐसा होगा, यही कारण है कि कभी-कभी (हालांकि हमेशा नहीं) दूसरे चरण में गलत फ़ाइलें जोड़ी जाती हैं।

इसे ठीक करने के लिए, हमें यह सुनिश्चित करना होगा कि BAM फ़ाइलें और उनकी इंडेक्स फ़ाइलें चैनलों के माध्यम से एक साथ यात्रा करें।

!!! tip "सुझाव"

    workflow कोड में `view()` statements कुछ नहीं करते, इसलिए उन्हें छोड़ना कोई समस्या नहीं है।
    हालाँकि वे तुम्हारे console आउटपुट को अव्यवस्थित करेंगे, इसलिए हम अनुशंसा करते हैं कि जब तुम समस्या का निवारण कर लो तो उन्हें हटा दो।

### 3.4. इंडेक्स फ़ाइलों को सही ढंग से हैंडल करने के लिए workflow अपडेट करो

fix यह है कि प्रत्येक BAM फ़ाइल को उसके इंडेक्स के साथ एक tuple में bundle करें, फिर downstream प्रोसेस और workflow plumbing को match करने के लिए अपडेट करें।

#### 3.4.1. SAMTOOLS_INDEX मॉड्यूल के आउटपुट को tuple में बदलो

यह सुनिश्चित करने का सबसे सरल तरीका है कि एक BAM फ़ाइल और उसका इंडेक्स निकटता से जुड़े रहें, उन्हें index task से बाहर आने वाले tuple में एक साथ package करना है।

!!! note "नोट"

    एक **tuple** तत्वों की एक finite, ordered list है जो आमतौर पर एक function से कई मान वापस करने के लिए उपयोग की जाती है। Tuples विशेष रूप से उपयोगी हैं जब प्रोसेस के बीच कई इनपुट या आउटपुट पास करते समय उनके association और order को संरक्षित करना हो।

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

इस तरह, प्रत्येक इंडेक्स फ़ाइल अपनी मूल BAM फ़ाइल के साथ कसकर जुड़ी होगी, और indexing चरण का समग्र आउटपुट फ़ाइलों के जोड़े युक्त एक single चैनल होगा।

#### 3.4.2. GATK_HAPLOTYPECALLER मॉड्यूल के इनपुट को tuple स्वीकार करने के लिए बदलो

चूंकि हमने पहले प्रोसेस के आउटपुट के 'shape' को बदल दिया है, हमें दूसरे प्रोसेस की इनपुट परिभाषा को match करने के लिए अपडेट करना होगा।

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

अब हमें प्रोसेस कॉल और publish targets में नई tuple संरचना को reflect करने के लिए workflow को अपडेट करना होगा।

#### 3.4.3. workflow में GATK_HAPLOTYPECALLER के कॉल को अपडेट करो

हमें अब `GATK_HAPLOTYPECALLER` प्रोसेस को मूल `reads_ch` प्रदान करने की आवश्यकता नहीं है, क्योंकि BAM फ़ाइल अब `SAMTOOLS_INDEX` द्वारा आउटपुट किए गए चैनल में bundled है।

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

अंत में, हमें नई आउटपुट संरचना को reflect करने के लिए publish targets को अपडेट करना होगा।

#### 3.4.4. indexed BAM आउटपुट के लिए publish target अपडेट करो

चूंकि SAMTOOLS_INDEX आउटपुट अब BAM फ़ाइल और उसके इंडेक्स दोनों युक्त एक tuple है, इसकी सामग्री को बेहतर ढंग से reflect करने के लिए publish target का नाम `bam_index` से `indexed_bam` में बदलो:

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

इन परिवर्तनों के साथ, BAM और उसका इंडेक्स एक साथ यात्रा करने की गारंटी है, इसलिए pairing हमेशा सही होगी।

### 3.5. सही किए गए workflow को चलाओ

आगे बढ़ने के लिए यह सुनिश्चित करने के लिए workflow को फिर से चलाओ कि यह विश्वसनीय रूप से काम करेगा।

```bash
nextflow run genomics.nf -profile test
```

इस बार (और हर बार) सब कुछ सही ढंग से चलना चाहिए:

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

results डायरेक्टरी में अब प्रत्येक नमूने के लिए BAM और BAI फ़ाइलें (tuple से) शामिल हैं, साथ ही VCF आउटपुट भी:

??? abstract "Results डायरेक्टरी सामग्री"

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

संबंधित फ़ाइलों को tuples में bundle करके, हमने सुनिश्चित किया कि सही फ़ाइलें हमेशा workflow के माध्यम से एक साथ यात्रा करती हैं।
workflow अब किसी भी संख्या में नमूनों को विश्वसनीय रूप से process करता है, लेकिन उन्हें config में व्यक्तिगत रूप से सूचीबद्ध करना बहुत scalable नहीं है।
अगले चरण में, हम एक फ़ाइल से इनपुट पढ़ने पर स्विच करेंगे।

### सारांश

तुम जानते हो कि अपने workflow को कई नमूनों पर (स्वतंत्र रूप से) चलाने योग्य कैसे बनाएँ।

### आगे क्या है?

bulk में नमूनों को हैंडल करना आसान बनाओ।

---

## 4. workflow को इनपुट फ़ाइलों के batch युक्त टेक्स्ट फ़ाइल स्वीकार करने योग्य बनाओ

workflow को कई डेटा इनपुट फ़ाइलें प्रदान करने का एक बहुत सामान्य तरीका फ़ाइल पथों युक्त टेक्स्ट फ़ाइल के साथ करना है।
यह एक टेक्स्ट फ़ाइल के रूप में सरल हो सकती है जो प्रति लाइन एक फ़ाइल पथ सूचीबद्ध करती है और कुछ नहीं, या फ़ाइल में अतिरिक्त metadata हो सकता है, जिस स्थिति में इसे अक्सर samplesheet कहा जाता है।

यहाँ हम तुम्हें दिखाने जा रहे हैं कि simple case कैसे करें।

### 4.1. इनपुट फ़ाइल पथों को सूचीबद्ध करने वाली प्रदान की गई टेक्स्ट फ़ाइल की जाँच करो

हमने पहले से ही इनपुट फ़ाइल पथों को सूचीबद्ध करने वाली एक टेक्स्ट फ़ाइल बनाई है, जिसे `sample_bams.txt` कहा जाता है, जिसे तुम `data/` डायरेक्टरी में पा सकते हो।

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

जैसा कि तुम देख सकते हो, हमने प्रति लाइन एक फ़ाइल पथ सूचीबद्ध किया है, और वे absolute paths हैं।

!!! note "नोट"

    हम यहाँ जो फ़ाइलें उपयोग कर रहे हैं वे केवल तुम्हारे GitHub Codespaces के local filesystem पर हैं, लेकिन हम cloud storage में फ़ाइलों की ओर भी इशारा कर सकते हैं।
    यदि तुम प्रदान किए गए Codespaces वातावरण का उपयोग नहीं कर रहे हो, तो तुम्हें अपने local setup से match करने के लिए फ़ाइल पथों को अनुकूलित करने की आवश्यकता हो सकती है।

### 4.2. पैरामीटर और test profile अपडेट करो

व्यक्तिगत नमूनों को सूचीबद्ध करने के बजाय `sample_bams.txt` फ़ाइल की ओर इशारा करने के लिए `reads_bam` पैरामीटर को switch करो।

params ब्लॉक में type annotation restore करो (क्योंकि यह फिर से एक single path है):

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

फिर टेक्स्ट फ़ाइल की ओर इशारा करने के लिए test profile को अपडेट करो:

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

फ़ाइलों की list अब कोड में बिल्कुल नहीं रहती, जो सही दिशा में एक बड़ा कदम है।

### 4.3. फ़ाइल से लाइनें पढ़ने के लिए चैनल फ़ैक्टरी अपडेट करो

वर्तमान में, हमारी इनपुट चैनल फ़ैक्टरी हमें दी गई किसी भी फ़ाइल को डेटा इनपुट के रूप में treat करती है जिसे हम indexing प्रोसेस को feed करना चाहते हैं।
चूंकि हम अब इसे एक फ़ाइल दे रहे हैं जो इनपुट फ़ाइल पथों को सूचीबद्ध करती है, हमें फ़ाइल को parse करने और इसमें शामिल फ़ाइल पथों को डेटा इनपुट के रूप में treat करने के लिए इसके व्यवहार को बदलने की आवश्यकता है।

हम उसी पैटर्न का उपयोग करके ऐसा कर सकते हैं जिसका हमने [Hello Nextflow के भाग 2](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file) में उपयोग किया था: फ़ाइल को parse करने के लिए [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) ऑपरेटर लागू करना, फिर प्रत्येक लाइन के पहले field को select करने के लिए एक `map` operation।

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

तकनीकी रूप से हम इसे [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) ऑपरेटर का उपयोग करके अधिक सरलता से कर सकते हैं, क्योंकि हमारी इनपुट फ़ाइल वर्तमान में केवल फ़ाइल पथों को contain करती है।
हालाँकि, अधिक versatile `splitCsv` ऑपरेटर (`map` द्वारा पूरक) का उपयोग करके, हम अपने workflow को future-proof कर सकते हैं यदि हम फ़ाइल पथों युक्त फ़ाइल में metadata जोड़ने का निर्णय लेते हैं।

!!! tip "सुझाव"

    यदि तुम्हें विश्वास नहीं है कि तुम समझते हो कि ऑपरेटर यहाँ क्या कर रहे हैं, तो यह उन्हें लागू करने से पहले और बाद में चैनल सामग्री कैसी दिखती है यह देखने के लिए `.view()` ऑपरेटर का उपयोग करने का एक और शानदार अवसर है।

### 4.4. workflow चलाओ

workflow को एक बार और चलाओ। इसे पहले जैसा ही परिणाम उत्पन्न करना चाहिए, है ना?

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

हाँ! वास्तव में, Nextflow सही ढंग से detect करता है कि प्रोसेस कॉल बिल्कुल समान हैं, और सब कुछ फिर से चलाने की परवाह भी नहीं करता, क्योंकि हम `-resume` के साथ चला रहे थे।

और बस! हमारे simple वेरिएंट कॉलिंग workflow में वे सभी बुनियादी features हैं जो हम चाहते थे।

### सारांश

तुम जानते हो कि BAM फ़ाइल को index करने और GATK का उपयोग करके प्रति-नमूना वेरिएंट कॉलिंग लागू करने के लिए multi-step modular workflow कैसे बनाएँ।

अधिक सामान्य रूप से, तुमने सीखा है कि genomics फ़ाइल formats और tool requirements की विशिष्टताओं को ध्यान में रखते हुए वास्तविक काम करने वाली एक simple genomics pipeline बनाने के लिए आवश्यक Nextflow components और logic का उपयोग कैसे करें।

### आगे क्या है?

अपनी सफलता का जश्न मनाओ और एक extra long break लो!

कोर्स के अगले भाग में, तुम सीखोगे कि डेटा पर joint वेरिएंट कॉलिंग लागू करने के लिए इस simple प्रति-नमूना वेरिएंट कॉलिंग workflow को कैसे transform करें।
