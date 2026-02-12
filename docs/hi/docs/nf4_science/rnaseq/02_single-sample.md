# भाग 2: एकल-नमूना कार्यान्वयन

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

कोर्स के इस भाग में, हम सबसे सरल संभव workflow लिखने जा रहे हैं जो भाग 1 में चलाए गए सभी commands को wrap करता है ताकि उन्हें चलाना automate हो सके, और हम एक बार में सिर्फ एक नमूना प्रोसेस करने का लक्ष्य रखेंगे।

!!! warning "पूर्वापेक्षा"

    इस पाठ को शुरू करने से पहले तुम्हें [भाग 1: विधि अवलोकन](./01_method.md) पर काम करना होगा।
    विशेष रूप से, section 1.2.3 पर काम करने से genome index फ़ाइल (`data/genome_index.tar.gz`) बनती है जो इस पाठ में alignment step के लिए आवश्यक है।

## असाइनमेंट

कोर्स के इस भाग में, हम एक workflow विकसित करने जा रहे हैं जो निम्नलिखित करता है:

1. इनपुट reads पर quality control (FastQC) चलाएं
2. Adapters को trim करें और post-trimming QC (Trim Galore) चलाएं
3. Trimmed reads को reference genome (HISAT2) से align करें

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-02.svg"
</figure>

यह [भाग 1: विधि अवलोकन](./01_method.md#1-single-sample-processing) के पहले section के steps को automate करता है, जहां तुमने इन commands को manually उनके containers में चलाया था।

शुरुआती बिंदु के रूप में, हम तुम्हें एक workflow फ़ाइल, `rnaseq.nf`, प्रदान करते हैं, जो workflow के मुख्य भागों की रूपरेखा देती है, साथ ही `modules/` डायरेक्टरी में चार मॉड्यूल फ़ाइलें (`fastqc.nf`, `trim_galore.nf`, `hisat2_align.nf`, और `multiqc.nf`) जो प्रत्येक process की संरचना की रूपरेखा देती हैं।

??? full-code "Scaffold फ़ाइलें"

    ```groovy title="rnaseq.nf"
    #!/usr/bin/env nextflow

    // Module INCLUDE statements

    /*
     * Pipeline parameters
     */

    // Primary input

    workflow {

        main:
        // Create input channel

        // Call processes

        publish:
        // Declare outputs to publish
    }

    output {
        // Configure publish targets
    }
    ```

    ```groovy title="modules/fastqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/trim_galore.nf"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/hisat2_align.nf"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/multiqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Aggregate QC reports with MultiQC
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

ये फ़ाइलें functional नहीं हैं; इनका उद्देश्य सिर्फ scaffolds के रूप में काम करना है जिन्हें तुम code के interesting भागों से भर सकते हो।

## पाठ योजना

Development process को अधिक शैक्षिक बनाने के लिए, हमने इसे तीन चरणों में विभाजित किया है:

1. **एक single-stage workflow लिखें जो शुरुआती QC step चलाता है।**
   यह CLI पैरामीटर सेट करना, इनपुट channel बनाना, process मॉड्यूल लिखना, और आउटपुट publishing को configure करना cover करता है।
2. **Adapter trimming और post-trimming QC जोड़ें।**
   यह एक process के आउटपुट को दूसरे के इनपुट से जोड़कर processes को chain करना introduce करता है।
3. **Reference genome के लिए alignment जोड़ें।**
   यह अतिरिक्त reference इनपुट को handle करना और compressed archives के साथ काम करना cover करता है।

प्रत्येक step workflow development के एक विशिष्ट पहलू पर focus करता है।

!!! tip "सुझाव"

     सुनिश्चित करो कि तुम सही working डायरेक्टरी में हो:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. एक single-stage workflow लिखें जो शुरुआती QC चलाता है

यह पहला step basics पर focus करता है: एक FASTQ फ़ाइल को load करना और उस पर quality control चलाना।

[भाग 1](01_method.md) से `fastqc` कमांड को याद करो:

```bash
fastqc <reads>
```

कमांड इनपुट के रूप में एक FASTQ फ़ाइल लेता है और एक quality control report को `.zip` archive और एक `.html` summary के रूप में produce करता है।
Container URI था `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`।

हम इस जानकारी को लेने जा रहे हैं और इसे Nextflow में तीन चरणों में wrap करेंगे:

1. इनपुट को सेट करें
2. QC process लिखें और इसे workflow में call करें
3. आउटपुट handling को configure करें

### 1.1. इनपुट को सेट करें

हमें एक इनपुट पैरामीटर declare करना होगा, एक सुविधाजनक default value प्रदान करने के लिए एक test profile बनाना होगा, और एक इनपुट channel बनाना होगा।

#### 1.1.1. एक इनपुट पैरामीटर declaration जोड़ें

`rnaseq.nf` में, `Pipeline parameters` section के तहत, `input` नाम का एक पैरामीटर `Path` type के साथ declare करो।

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        input: Path
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

यह CLI पैरामीटर को सेट करता है, लेकिन हम development के दौरान हर बार workflow चलाते समय फ़ाइल path को type नहीं करना चाहते।
Default value प्रदान करने के लिए कई विकल्प हैं; यहां हम एक test profile का उपयोग करते हैं।

#### 1.1.2. `nextflow.config` में default value के साथ एक test profile बनाएं

एक test profile command line पर इनपुट specify किए बिना workflow को try करने के लिए सुविधाजनक default values प्रदान करता है।
यह Nextflow ecosystem में एक सामान्य convention है (अधिक विवरण के लिए [Hello Config](../../hello_nextflow/06_hello_config.md) देखें)।

`nextflow.config` में एक `profiles` block जोड़ें जिसमें एक `test` profile हो जो `input` पैरामीटर को test FASTQ फ़ाइलों में से एक पर सेट करता है।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

यहां, हम `#!groovy ${projectDir}` का उपयोग कर रहे हैं, एक built-in Nextflow variable जो उस डायरेक्टरी को point करता है जहां workflow script स्थित है।
यह absolute paths को hardcode किए बिना data फ़ाइलों और अन्य resources को reference करना आसान बनाता है।

पैरामीटर के पास अब एक सुविधाजनक default है। अगला, हमें इससे एक channel बनाना होगा।

#### 1.1.3. इनपुट channel को सेट करें

Workflow block में, `.fromPath` channel factory का उपयोग करके पैरामीटर value से एक इनपुट channel बनाएं ([Hello Channels](../../hello_nextflow/02_hello_channels.md) में उपयोग किया गया)।

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // एक फ़ाइल path से इनपुट channel बनाएं
        read_ch = channel.fromPath(params.input)

        // Processes को call करें

        publish:
        // Declare outputs to publish
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="13"
    workflow {

        main:
        // Create input channel

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

अगला, हमें इस इनपुट पर QC चलाने के लिए process बनाना होगा।

### 1.2. QC process लिखें और इसे workflow में call करें

हमें मॉड्यूल फ़ाइल में process definition को भरना होगा, include statement का उपयोग करके इसे workflow में import करना होगा, और इसे इनपुट पर call करना होगा।

#### 1.2.1. QC process के लिए मॉड्यूल को भरें

`modules/fastqc.nf` खोलें और process definition की रूपरेखा को examine करें।
तुम्हें मुख्य structural elements को पहचानना चाहिए; यदि नहीं, तो refresher के लिए [Hello Nextflow](../../hello_nextflow/01_hello_world.md) पढ़ने पर विचार करो।

आगे बढ़ो और ऊपर दी गई जानकारी का उपयोग करके process definition को खुद भरो, फिर नीचे "बाद में" tab में solution के खिलाफ अपने काम की जांच करो।

=== "पहले"

    ```groovy title="modules/fastqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "बाद में"

    ```groovy title="modules/fastqc.nf" linenums="1" hl_lines="8 11 14 15 19"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html

        script:
        """
        fastqc ${reads}
        """
    }
    ```

`simpleName` accessor filename से सभी extensions को strip करता है, इसलिए `ENCSR000COQ1_1.fastq.gz` `ENCSR000COQ1_1` बन जाता है।
हम प्रत्येक आउटपुट channel को नाम assign करने के लिए `emit:` syntax का उपयोग करते हैं, जो outputs को publish block में wire करने के लिए उपयोगी होगा।

एक बार जब तुमने यह पूरा कर लिया, तो process पूरा हो गया है।
इसे workflow में उपयोग करने के लिए, तुम्हें मॉड्यूल को import करना होगा और एक process call जोड़ना होगा।

#### 1.2.2. मॉड्यूल को include करें

`rnaseq.nf` में, process को workflow के लिए उपलब्ध कराने के लिए एक `include` statement जोड़ें:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    ```

Process अब workflow scope में उपलब्ध है।

#### 1.2.3. इनपुट पर QC process को call करें

Workflow block में `FASTQC` के लिए एक call जोड़ें, इनपुट channel को argument के रूप में pass करते हुए।

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // एक फ़ाइल path से इनपुट channel बनाएं
        read_ch = channel.fromPath(params.input)

        // प्रारंभिक quality control
        FASTQC(read_ch)

        publish:
        // Declare outputs to publish
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="14"
    workflow {

        main:
        // एक फ़ाइल path से इनपुट channel बनाएं
        read_ch = channel.fromPath(params.input)

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

Workflow अब इनपुट को load करता है और उस पर QC process चलाता है।
अगला, हमें configure करना होगा कि आउटपुट कैसे publish किया जाता है।

### 1.3. आउटपुट handling को configure करें

हमें declare करना होगा कि कौन से process outputs को publish करना है और specify करना होगा कि उन्हें कहां जाना चाहिए।

#### 1.3.1. `publish:` section में outputs को declare करें

Workflow block के अंदर `publish:` section declare करता है कि कौन से process outputs को publish किया जाना चाहिए।
`FASTQC` के outputs को named targets को assign करें।

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="23" hl_lines="2-3"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="23"
        publish:
        // Declare outputs to publish
    }
    ```

अगला, हमें Nextflow को बताना होगा कि published outputs को कहां रखना है।

#### 1.3.2. `output {}` block में आउटपुट targets को configure करें

`output {}` block workflow के बाहर बैठता है और specify करता है कि प्रत्येक named target को कहां publish किया जाता है।
दोनों targets को `fastqc/` subdirectory में publish करने के लिए configure करें।

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="28" hl_lines="2-7"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="28"
    output {
        // Configure publish targets
    }
    ```

!!! note "नोट"

    डिफ़ॉल्ट रूप से, Nextflow आउटपुट फ़ाइलों को symbolic links के रूप में publish करता है, जो अनावश्यक duplication से बचता है।
    भले ही हम यहां जो data फ़ाइलें उपयोग कर रहे हैं वे बहुत छोटी हैं, genomics में वे बहुत बड़ी हो सकती हैं।
    Symlinks टूट जाएंगे जब तुम अपनी `work` डायरेक्टरी को साफ करोगे, इसलिए production workflows के लिए तुम default publish mode को `'copy'` में override करना चाह सकते हो।

### 1.4. Workflow को चलाएं

इस बिंदु पर, हमारे पास एक one-step QC workflow है जो पूरी तरह से functional होना चाहिए।

हम test profile में सेट किए गए default value का उपयोग करने के लिए `-profile test` के साथ चलाते हैं, command line पर path लिखने की आवश्यकता से बचते हुए।

```bash
nextflow run rnaseq.nf -profile test
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [mad_lorenz] DSL2 - revision: 5846a164d2

    executor >  local (1)
    [7b/8ee79e] FASTQC (1) | 1 of 1 ✔
    ```

यह बहुत जल्दी चलना चाहिए अगर तुमने भाग 1 पर काम किया है और पहले से ही कंटेनर को pull कर लिया है।
यदि तुमने इसे छोड़ दिया है, तो Nextflow तुम्हारे लिए कंटेनर को pull करेगा; इसके होने के लिए तुम्हें कुछ भी करने की आवश्यकता नहीं है, लेकिन तुम्हें एक मिनट तक प्रतीक्षा करने की आवश्यकता हो सकती है।

तुम results डायरेक्टरी में outputs की जांच कर सकते हो।

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

नमूने के लिए QC reports अब `fastqc/` subdirectory में publish हैं।

### सारांश

तुम जानते हो कि एक process युक्त मॉड्यूल कैसे बनाएं, इसे workflow में import करें, इसे इनपुट channel के साथ call करें, और workflow-level output block का उपयोग करके results को publish करें।

### आगे क्या है?

Workflow में दूसरे step के रूप में post-trimming QC के साथ adapter trimming जोड़ें।

---

## 2. Adapter trimming और post-trimming QC जोड़ें

अब जब हमारे पास शुरुआती QC है, तो हम अपने built-in post-trimming QC के साथ adapter trimming step जोड़ सकते हैं।

[भाग 1](01_method.md) से `trim_galore` कमांड को याद करो:

```bash
trim_galore --fastqc <reads>
```

कमांड एक FASTQ फ़ाइल से adapters को trim करता है और trimmed आउटपुट पर FastQC चलाता है।
यह trimmed reads, एक trimming report, और trimmed reads के लिए FastQC reports produce करता है।
Container URI था `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`।

हमें बस process definition लिखना है, इसे import करना है, इसे workflow में call करना है, और आउटपुट handling को update करना है।

### 2.1. Trimming process लिखें और इसे workflow में call करें

पहले की तरह, हमें process definition को भरना होगा, मॉड्यूल को import करना होगा, और process call जोड़ना होगा।

#### 2.1.1. Trimming process के लिए मॉड्यूल को भरें

`modules/trim_galore.nf` खोलें और process definition की रूपरेखा को examine करें।

आगे बढ़ो और ऊपर दी गई जानकारी का उपयोग करके process definition को खुद भरो, फिर नीचे "बाद में" tab में solution के खिलाफ अपने काम की जांच करो।

=== "पहले"

    ```groovy title="modules/trim_galore.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "बाद में"

    ```groovy title="modules/trim_galore.nf" linenums="1" hl_lines="8 11 14 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    }
    ```

इस process के तीन named outputs हैं: trimmed reads जो alignment step में feed होते हैं, trimming report, और post-trimming FastQC reports।
`--fastqc` flag Trim Galore को trimmed आउटपुट पर automatically FastQC चलाने के लिए कहता है।

#### 2.1.2. मॉड्यूल को include करें

नए मॉड्यूल को import करने के लिए `rnaseq.nf` को update करें:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    ```

अगला, हम workflow में process call जोड़ेंगे।

#### 2.1.3. इनपुट पर trimming process को call करें

Workflow block में process call जोड़ें:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="15" hl_lines="10-11"
    workflow {

        main:
        // एक फ़ाइल path से इनपुट channel बनाएं
        read_ch = channel.fromPath(params.input)

        // प्रारंभिक quality control
        FASTQC(read_ch)

        // Adapter trimming और post-trimming QC
        TRIM_GALORE(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="15"
    workflow {

        main:
        // एक फ़ाइल path से इनपुट channel बनाएं
        read_ch = channel.fromPath(params.input)

        // प्रारंभिक quality control
        FASTQC(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

Trimming process अब workflow में wired है।

### 2.2. आउटपुट handling को update करें

हमें trimming outputs को publish declaration में जोड़ना होगा और configure करना होगा कि वे कहां जाते हैं।

#### 2.2.1. Trimming outputs के लिए publish targets जोड़ें

`publish:` section में trimming outputs जोड़ें:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="27" hl_lines="4-6"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="27"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

अगला, हमें Nextflow को बताना होगा कि इन outputs को कहां रखना है।

#### 2.2.2. नए आउटपुट targets को configure करें

`output {}` block में trimming targets के लिए entries जोड़ें, उन्हें `trimming/` subdirectory में publish करते हुए:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="35" hl_lines="8-16"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="35"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

आउटपुट configuration पूरा हो गया है।

### 2.3. Workflow को चलाएं

Workflow में अब शुरुआती QC और adapter trimming दोनों शामिल हैं।

```bash
nextflow run rnaseq.nf -profile test
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [gloomy_becquerel] DSL2 - revision: bb11055736

    executor >  local (2)
    [f6/c8ef2e] FASTQC (1)      | 1 of 1 ✔
    [58/c58d8a] TRIM_GALORE (1) | 1 of 1 ✔
    ```

यह भी बहुत जल्दी चलना चाहिए, क्योंकि हम इतनी छोटी इनपुट फ़ाइल पर चला रहे हैं।

तुम results डायरेक्टरी में trimming outputs पा सकते हो।

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

Trimming outputs और post-trimming QC reports अब `trimming/` subdirectory में हैं।

### सारांश

तुम जानते हो कि एक दूसरा processing step कैसे जोड़ें जो same इनपुट पर independently चलता है, कई named outputs produce करते हुए।

### आगे क्या है?

Alignment step जोड़ें जो trimmed reads आउटपुट से chain off करता है।

---

## 3. Reference genome के लिए alignment जोड़ें

अंत में हम HISAT2 का उपयोग करके genome alignment step जोड़ सकते हैं।

[भाग 1](01_method.md) से alignment कमांड को याद करो:

```bash
hisat2 -x <genome_index> -U <reads> \
    --new-summary --summary-file <reads>.hisat2.log | \
    samtools view -bS -o <reads>.bam
```

कमांड reads को reference genome से align करता है और आउटपुट को BAM format में convert करता है।
इसे एक pre-built genome index archive की आवश्यकता होती है और एक BAM फ़ाइल और एक alignment summary log produce करता है।
Container URI था `community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e`।

इस process को एक अतिरिक्त इनपुट (genome index archive) की आवश्यकता है, इसलिए हमें पहले इसे सेट करना होगा, फिर process को लिखना और wire करना होगा।

### 3.1. इनपुट को सेट करें

हमें genome index archive के लिए एक पैरामीटर declare करना होगा।

#### 3.1.1. Genome index के लिए एक पैरामीटर जोड़ें

`rnaseq.nf` में genome index archive के लिए एक पैरामीटर declaration जोड़ें:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="5-6"
    params {
        // Primary input
        input: Path

        // Reference genome archive
        hisat2_index_zip: Path
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Primary input
        input: Path
    }
    ```

#### 3.1.2. Test profile में genome index default जोड़ें

जैसा कि हमने section 1.1.2 में `input` के लिए किया था, `nextflow.config` में test profile में genome index के लिए एक default value जोड़ें:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="6"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

पैरामीटर तैयार है; अब हम alignment process बना सकते हैं।

### 3.2. Alignment process लिखें और इसे workflow में call करें

पहले की तरह, हमें process definition को भरना होगा, मॉड्यूल को import करना होगा, और process call जोड़ना होगा।

#### 3.2.1. Alignment process के लिए मॉड्यूल को भरें

`modules/hisat2_align.nf` खोलें और process definition की रूपरेखा को examine करें।

आगे बढ़ो और ऊपर दी गई जानकारी का उपयोग करके process definition को खुद भरो, फिर नीचे "बाद में" tab में solution के खिलाफ अपने काम की जांच करो।

=== "पहले"

    ```groovy title="modules/hisat2_align.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "बाद में"

    ```groovy title="modules/hisat2_align.nf" linenums="1" hl_lines="8 11 12 15 16 20 21 22 23"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"

        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    }
    ```

यह process दो इनपुट लेता है: reads और genome index archive।
Script block पहले archive से index को extract करता है, फिर HISAT2 alignment को `samtools view` में pipe करता है ताकि आउटपुट को BAM format में convert किया जा सके।
`index_zip` पर `simpleName` accessor archive के base name (`genome_index`) को extract करता है ताकि इसे index prefix के रूप में उपयोग किया जा सके।

#### 3.2.2. मॉड्यूल को include करें

नए मॉड्यूल को import करने के लिए `rnaseq.nf` को update करें:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="4"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

अगला, हम workflow में process call जोड़ेंगे।

#### 3.2.3. Alignment process को call करें

Trimmed reads पिछले step द्वारा आउटपुट किए गए `TRIM_GALORE.out.trimmed_reads` channel में हैं।
हम genome index archive प्रदान करने के लिए `#!groovy file(params.hisat2_index_zip)` का उपयोग करते हैं।

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="14-15"
    workflow {

        main:
        // एक फ़ाइल path से इनपुट channel बनाएं
        read_ch = channel.fromPath(params.input)

        // प्रारंभिक quality control
        FASTQC(read_ch)

        // Adapter trimming और post-trimming QC
        TRIM_GALORE(read_ch)

        // Reference genome के साथ alignment
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // एक फ़ाइल path से इनपुट channel बनाएं
        read_ch = channel.fromPath(params.input)

        // प्रारंभिक quality control
        FASTQC(read_ch)

        // Adapter trimming और post-trimming QC
        TRIM_GALORE(read_ch)
    ```

Alignment process अब workflow में wired है।

### 3.3. आउटपुट handling को update करें

हमें alignment outputs को publish declaration में जोड़ना होगा और configure करना होगा कि वे कहां जाते हैं।

#### 3.3.1. Alignment outputs के लिए publish targets जोड़ें

`publish:` section में alignment outputs जोड़ें:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="34" hl_lines="7-8"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="34"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

अगला, हमें Nextflow को बताना होगा कि इन outputs को कहां रखना है।

#### 3.3.2. नए आउटपुट targets को configure करें

`output {}` block में alignment targets के लिए entries जोड़ें, उन्हें `align/` subdirectory में publish करते हुए:

=== "बाद में"

    ```groovy title="rnaseq.nf" linenums="44" hl_lines="17-22"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
        align_log {
            path 'align'
        }
    }
    ```

=== "पहले"

    ```groovy title="rnaseq.nf" linenums="44"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

आउटपुट configuration पूरा हो गया है।

### 3.4. Workflow को चलाएं

Workflow में अब तीनों processing steps शामिल हैं: QC, trimming, और alignment।

```bash
nextflow run rnaseq.nf -profile test
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [elated_stonebraker] DSL2 - revision: e8e57d0cdd

    executor >  local (3)
    [e8/fa29d6] FASTQC (1)       | 1 of 1 ✔
    [ca/ffdde2] TRIM_GALORE (1)  | 1 of 1 ✔
    [b6/1c6ca3] HISAT2_ALIGN (1) | 1 of 1 ✔
    ```

तुम results डायरेक्टरी में alignment outputs पा सकते हो।

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

यह प्रत्येक नमूने पर लागू करने के लिए आवश्यक basic processing को पूरा करता है।

_हम भाग 3 में MultiQC report aggregation जोड़ेंगे, जब हम workflow को एक बार में कई नमूने स्वीकार करने योग्य बना लेंगे।_

---

### सारांश

तुम जानते हो कि single-end RNAseq नमूनों को व्यक्तिगत रूप से प्रोसेस करने के लिए सभी मुख्य steps को कैसे wrap करें।

### आगे क्या है?

ब्रेक लो! यह बहुत कुछ था।

जब तुम तरोताजा महसूस करो, तो [भाग 3](./03_multi-sample.md) पर जाओ, जहां तुम सीखोगे कि कैसे workflow को कई नमूनों को parallel में प्रोसेस करने के लिए modify करें, सभी नमूनों के लिए सभी steps में QC reports को aggregate करें, और paired-end RNAseq data पर workflow चलाना enable करें।
