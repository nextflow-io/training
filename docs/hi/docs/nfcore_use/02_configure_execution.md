# भाग 2: पाइपलाइन निष्पादन को कॉन्फ़िगर करना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[भाग 1](./01_run_demo.md) में, तुमने nf-core/demo पाइपलाइन को उसके test profile का उपयोग करके खोजा और चलाया।
अब हम देखेंगे कि पाइपलाइन निष्पादन को कैसे कॉन्फ़िगर किया जाए: पैरामीटर सेट करना, वैलिडेशन को समझना, और रिसोर्स आवंटन तथा टूल आर्गुमेंट को कस्टमाइज़ करना।

जैसा कि [Hello Config](../hello_nextflow/06_hello_config.md) में समझाया गया है, हम चाहते हैं कि पाइपलाइन कोड को बदले बिना यह बदल सकें कि हमारी पाइपलाइन किस डेटा पर और कैसे चलेगी।
इसके लिए, Nextflow पाइपलाइन कॉन्फ़िगरेशन को नियंत्रित करने के कई तरीके सपोर्ट करता है, जो थोड़ा भारी लग सकता है।

nf-core प्रोजेक्ट कॉन्फ़िगरेशन तत्वों को व्यवस्थित करने के लिए नियम निर्धारित करता है, और शीर्ष स्तर पर दो प्रकार के कॉन्फ़िगरेशन को अलग करता है: **पाइपलाइन पैरामीटर** और सख्त अर्थ में **कॉन्फ़िगरेशन**।

- **पाइपलाइन पैरामीटर** (`params` सिस्टम के माध्यम से सेट किए जाते हैं) में आमतौर पर इनपुट फ़ाइलें, टूल व्यवहार फ्लैग और विश्लेषण पैरामीटर शामिल होते हैं।
- सख्त अर्थ में **कॉन्फ़िगरेशन** का मतलब है पाइपलाइन कैसे चलती है इसकी लॉजिस्टिक्स, यानी executor, कंप्यूट रिसोर्स आवंटन आदि।

<figure class="excalidraw">
    --8<-- "docs/en/docs/nfcore_use/img/params_vs_config.excalidraw.svg"
</figure>

चलो पहले पाइपलाइन पैरामीटर से शुरू करते हैं, फिर सख्त अर्थ में कॉन्फ़िगरेशन देखेंगे।

---

## 1. पाइपलाइन पैरामीटर

सभी nf-core पाइपलाइनों के लिए, तुम `--help` फ्लैग का उपयोग करके सीधे कमांड लाइन से पाइपलाइन पैरामीटर की पूरी सूची प्राप्त कर सकते हो, जो खुद एक पाइपलाइन पैरामीटर है।

### 1.1. `--help` से पैरामीटर की सूची प्राप्त करें

demo पाइपलाइन के लिए help कमांड चलाओ:

```bash
nextflow run nf-core/demo --help
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [adoring_meucci] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>


    Input/output options
      --input                       [string] Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string] The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
      --email                       [string] Email address for completion summary.
      --multiqc_title               [string] MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string] Name of iGenomes reference.
      --fasta                       [string] Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean] Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]          Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string] Display the help message.
      --help_full                   [boolean]         Display the full detailed help message.
      --show_hidden                 [boolean]         Display hidden parameters in the help message (only works when --help or --help_full are provided).
     !! Hiding 19 param(s), use the `--showHidden` parameter to show them !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

जैसा कि तुम देख सकते हो, आउटपुट पैरामीटर को श्रेणियों में समूहित करता है (Input/output options, Reference genome options, आदि) और प्रत्येक के लिए प्रकार और विवरण देता है।

यह वर्गीकरण एक schema फ़ाइल द्वारा निर्धारित होता है, जिसे नीचे और विस्तार से कवर किया गया है।
सामान्य Nextflow पाइपलाइनों में, `--help` तभी काम करता है जब डेवलपर ने इसे मैन्युअल रूप से लागू किया हो।

!!! tip "सुझाव"

    अतिरिक्त पैरामीटर देखने के लिए `--help --show_hidden` का उपयोग करो जो डिफ़ॉल्ट रूप से छिपे होते हैं, जैसे `--publish_dir_mode` या `--monochrome_logs`।

### 1.2. पैरामीटर मान सेट करें

जैसा कि [Hello Config](../hello_nextflow/06_hello_config.md) में कवर किया गया है, तुम कमांड लाइन पर `--param_name` के साथ पैरामीटर मान सेट कर सकते हो या एक YAML फ़ाइल में पैरामीटर का एक सेट एकत्र करके `-params-file` के साथ पास कर सकते हो।
दोनों तरीके nf-core पाइपलाइनों के साथ एक ही तरह काम करते हैं।

उदाहरण के लिए, trimming चरण को छोड़ने के लिए, हम boolean पैरामीटर `skip_trim` को `true` पर सेट करना चाहते हैं।
तुम्हारी working directory में `my_params.yml` नाम की एक params फ़ाइल दी गई है जिसमें वह मान पहले से सेट है:

```yaml title="my_params.yml"
skip_trim: true
```

इसे `-params-file` के साथ पास करो:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-notrim -params-file my_params.yml
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [focused_heisenberg] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results-notrim

    Process skipping options
      skip_trim                 : true

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_22-08-47

    Core Nextflow options
      revision                  : master
      runName                   : focused_heisenberg
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               | 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

`SEQTK_TRIM` प्रोसेस अब आउटपुट में नहीं दिखती।

!!! warning "पैरामीटर इनपुट के बारे में महत्वपूर्ण सीमाएं"

    **कमांड लाइन पर boolean पैरामीटर सेट करना**

    Nextflow संस्करण 26.04 से शुरू होकर, कमांड लाइन पर दिए गए सभी मान string के रूप में टाइप किए जाते हैं।
    `skip_trim` जैसे boolean पैरामीटर के लिए, इसे bare flag (`--skip_trim`) के रूप में या `--skip_trim true` के रूप में पास करने पर यह **string** `"true"` के रूप में मूल्यांकित होता है, जो schema वैलिडेशन में विफल हो जाता है:

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    किसी boolean पैरामीटर को वास्तविक `true`/`false` मान पर सेट करने के लिए, ऊपर दिखाए अनुसार `-params-file` का उपयोग करो, या इसे एक config फ़ाइल में सेट करो।
    String, integer और file-path पैरामीटर प्रभावित नहीं होते और अभी भी सीधे कमांड लाइन पर सेट किए जा सकते हैं।
    यह कोर्स boolean पैरामीटर के लिए इस पैटर्न का उपयोग करता है।

    **कस्टम कॉन्फ़िगरेशन फ़ाइलों का उपयोग करना**

    हालांकि तकनीकी रूप से `-c` के साथ पास की गई कस्टम कॉन्फ़िगरेशन फ़ाइल में पाइपलाइन पैरामीटर सेट करना संभव है, लेकिन Nextflow के कॉन्फ़िगरेशन प्राथमिकता नियमों के आधार पर यह पाइपलाइन के अपने `nextflow.config` में पहले से सेट डिफ़ॉल्ट को ओवरराइड नहीं कर सकता।
    कमांड लाइन पर `--param_name` या `-params-file` का उपयोग करना अधिक विश्वसनीय है, क्योंकि ये हमेशा प्राथमिकता लेते हैं।

    एक सामान्य नियम के रूप में: यदि यह `--help` आउटपुट में दिखता है, तो इसे config फ़ाइल के बजाय कमांड लाइन या params फ़ाइल के माध्यम से सेट करो।

### 1.3. पैरामीटर वैलिडेशन

एक रोचक तथ्य: `--help` कमांड सभी nf-core पाइपलाइनों के लिए काम करता है क्योंकि nf-core प्रोजेक्ट डेवलपर्स को सभी पाइपलाइन पैरामीटर को एक JSON schema फ़ाइल (`nextflow_schema.json`) में औपचारिक रूप से परिभाषित करने की आवश्यकता होती है।
यह schema प्रत्येक पैरामीटर का प्रकार, विवरण, डिफ़ॉल्ट मान और समूहीकरण रिकॉर्ड करता है।

`--help` आउटपुट को शक्ति देने के अलावा, schema फ़ाइल लॉन्च के समय स्वचालित वैलिडेशन भी सक्षम करती है।
इसका मतलब है कि Nextflow जांच कर सकता है कि तुम्हारे द्वारा पास किया गया हर पैरामीटर मौजूद है और उसे उचित मान दिया गया है (उचित प्रकार का, अनुमत मानों की सीमा के भीतर आदि)।

हम इसे [इनपुट वैलिडेशन अनुभाग](../nfcore_build/04_input_validation.md) में अधिक विस्तार से कवर करते हैं, लेकिन तुम demo पाइपलाइन को कुछ अमान्य पैरामीटर इनपुट देकर इसे पहले से ही काम करते देख सकते हो।

#### 1.3.1. अपरिचित पैरामीटर

एक ऐसा पैरामीटर पास करने की कोशिश करो जो मौजूद नहीं है:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --foobar "invalid"
```

कंसोल आउटपुट में एक चेतावनी शामिल है:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

पाइपलाइन फिर भी चलती है, लेकिन चेतावनी तुरंत तुम्हें सचेत करती है कि `--foobar` एक मान्यता प्राप्त पैरामीटर नहीं है।
यह तुम्हारा ध्यान non-breaking टाइपो की ओर आकर्षित करने के लिए है, जैसे `--outdir` के बजाय `--outDir` का उपयोग करना, जो तुम्हें समय और कंप्यूट बर्बाद करने से बचाने में मदद कर सकता है।

#### 1.3.2. अमान्य पैरामीटर मान

वैलिडेशन पैरामीटर **मानों** की भी जांच करता है।
`--skip_trim` पैरामीटर एक boolean फ्लैग है, इसलिए string मान पास करने पर पाइपलाइन तुरंत विफल हो जाती है:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

कोई भी प्रोसेस चलने से पहले पाइपलाइन रुक जाती है, जिससे तुम एक विफल या गलत निष्पादन से बच जाते हो।
जैसा कि [1.2](#12-set-parameter-values) में उल्लेख किया गया है, boolean पैरामीटर को कमांड लाइन पर पास करने के बजाय params फ़ाइल में वास्तविक `true`/`false` मान पर सेट किया जाना चाहिए, क्योंकि कमांड-लाइन मान string के रूप में टाइप किए जाते हैं।

### 1.4. इनपुट वैलिडेशन

वही वैलिडेशन लॉजिक इनपुट फ़ाइलों की वैधता जांचने के लिए भी उपयोग किया जा सकता है।
उदाहरण के लिए, यदि कोई पाइपलाइन अपने मुख्य डेटा इनपुट के रूप में एक samplesheet की अपेक्षा करती है (जो कि कई nf-core पाइपलाइनों का मामला है), तो डेवलपर एक इनपुट schema (पैरामीटर schema से अलग) प्रदान कर सकता है जो बताता है कि इनपुट फ़ाइल कैसे संरचित होनी चाहिए।

फिर, runtime पर, Nextflow जांच कर सकता है कि प्रदान की गई इनपुट फ़ाइल वैध है।

हम इसे [इनपुट वैलिडेशन अनुभाग](../nfcore_build/04_input_validation.md) में भी अधिक विस्तार से कवर करते हैं, लेकिन तुम demo पाइपलाइन को एक अमान्य इनपुट samplesheet देकर इसे पहले से ही काम करते देख सकते हो।

`nf-core/demo` पाइपलाइन `sample`, `fastq_1`, और `fastq_2` कॉलम वाली एक CSV फ़ाइल की अपेक्षा करती है।
यह एक schema फ़ाइल (`assets/schema_input.json`) में परिभाषित है जो अपेक्षित संरचना, कॉलम प्रकार और बाधाओं को निर्दिष्ट करती है।

??? abstract "इनपुट के लिए Schema फ़ाइल"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

Schema निर्दिष्ट करता है कि `sample` और `fastq_1` आवश्यक हैं, जबकि `fastq_2` वैकल्पिक है (paired-end और single-end दोनों डेटा को सपोर्ट करता है)।
फ़ाइल पथों को अस्तित्व और extension pattern के लिए वैलिडेट किया जाता है।

इसे प्रदर्शित करने के लिए, हम तुम्हारी working directory में `malformed_samplesheet.csv` नाम की एक खराब samplesheet प्रदान करते हैं:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

इस samplesheet में आवश्यक `fastq_1` कॉलम गायब है और `fastq_2` में एक गैर-मौजूद फ़ाइल पथ है।

`malformed_samplesheet.csv` को इनपुट के रूप में उपयोग करके demo पाइपलाइन चलाओ:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

जैसा कि तुम देख सकते हो, पाइपलाइन तुरंत विफल हो जाती है और **सभी** वैलिडेशन त्रुटियों को एक साथ रिपोर्ट करती है।
nf-schema पहली त्रुटि पर नहीं रुकता — यह हर समस्या को एकत्र करता है और उन्हें एक साथ सूचीबद्ध करता है, ताकि तुम एक बार में सब कुछ ठीक कर सको बजाय एक-एक करके समस्याएं खोजने के।

प्रत्येक त्रुटि उस सटीक entry और field की पहचान करती है जिसने समस्या उत्पन्न की, ताकि तुम अपनी samplesheet ठीक कर सको और फिर पाइपलाइन को इस विश्वास के साथ फिर से लॉन्च कर सको कि जब Nextflow वास्तव में फ़ाइल पथ तक पहुंचने की कोशिश करेगा तो यह बाद में विफल नहीं होगी।

डेवलपर्स के लिए, यह सब [Build with nf-core के भाग 4](../nfcore_build/04_input_validation.md) में अधिक विस्तार से कवर किया गया है।

### सारांश

तुम जानते हो कि `--help` के साथ पाइपलाइन के पैरामीटर की पूरी सूची कैसे प्राप्त करें, उन्हें कमांड लाइन या params फ़ाइल के माध्यम से कैसे सेट करें, और Nextflow पाइपलाइन के schema के विरुद्ध पैरामीटर मानों और इनपुट फ़ाइलों दोनों को कैसे वैलिडेट करता है।

### आगे क्या है?

दूसरे प्रकार के कॉन्फ़िगरेशन के बारे में जानो: पाइपलाइन कैसे चलती है, रिसोर्स आवंटन और टूल आर्गुमेंट को कवर करते हुए।

---

## 2. कॉन्फ़िगरेशन

सख्त अर्थ में कॉन्फ़िगरेशन नियंत्रित करती है कि पाइपलाइन **कैसे** चलती है: रिसोर्स आवंटन, टूल-विशिष्ट आर्गुमेंट, jobs कहां execute होते हैं, और कौन सा सॉफ़्टवेयर पैकेजिंग सिस्टम उपयोग करना है।

nf-core पाइपलाइनों में `nextflow.config` और `conf/` डायरेक्टरी में डिफ़ॉल्ट कॉन्फ़िगरेशन शामिल है।
कुछ भी ओवरराइड करने से पहले, यह जानना मददगार है कि डिफ़ॉल्ट कहां रहते हैं।

### 2.1. कॉन्फ़िगरेशन फ़ाइलें एक्सप्लोर करें

तुमने [भाग 1](./01_run_demo.md) में पहले ही देखा कि पाइपलाइन सोर्स कोड `$NXF_HOME/assets` के अंतर्गत रहता है।
[भाग 1](./01_run_demo.md) में बनाए गए `pipelines` symlink का उपयोग करके, उपलब्ध config फ़ाइलें देखने के लिए सूचीबद्ध करो:

```bash
ls pipelines/nf-core/demo/conf/
```

```console
base.config
containers_conda_lock_files_amd64.config
containers_conda_lock_files_arm64.config
containers_docker_amd64.config
containers_docker_arm64.config
containers_singularity_https_amd64.config
containers_singularity_https_arm64.config
containers_singularity_oras_amd64.config
containers_singularity_oras_arm64.config
igenomes.config
igenomes_ignored.config
modules.config
test.config
test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/nfcore_use/img/nfcore_config_files.excalidraw.svg"
</figure>

सबसे महत्वपूर्ण कॉन्फ़िगरेशन फ़ाइलें हैं:

- **`conf/base.config`**: रिसोर्स लेबल (`process_low`, `process_medium`, `process_high`) परिभाषित करता है जो प्रोसेस को CPUs, memory और time असाइन करते हैं। जब तुम देखते हो कि कोई प्रोसेस अपेक्षा से अधिक रिसोर्स उपयोग कर रही है, तो वे डिफ़ॉल्ट यहीं से आते हैं।
- **`conf/modules.config`**: प्रति-प्रोसेस टूल आर्गुमेंट (`ext.args`) और आउटपुट publishing सेटिंग (`publishDir`) सेट करता है। यह फ़ाइल खोलो यह देखने के लिए कि प्रत्येक टूल को डिफ़ॉल्ट रूप से कौन से आर्गुमेंट मिलते हैं।
- **`conf/test.config`**: [भाग 1](./01_run_demo.md) में उपयोग किया गया test profile, जो `resourceLimits` के माध्यम से रिसोर्स को सीमित करता है और एक test samplesheet सेट करता है। `-profile test` के साथ सक्रिय होता है।
  पूर्ण आकार के test dataset के साथ चलाने के लिए `conf/test_full.config` भी है, जो benchmarking के लिए उपयोगी है।

केंद्रीय `nextflow.config` उपरोक्त सभी को लोड करता है और सब कुछ के लिए उचित डिफ़ॉल्ट मान सेट करता है।

यदि तुम इन फ़ाइलों में निर्दिष्ट किसी भी सेटिंग को संशोधित करना चाहते हो, तो इनमें से किसी भी फ़ाइल को सीधे संशोधित मत करो।
इसके बजाय, अपनी खुद की config फ़ाइल बनाओ और इसे `-c` के साथ पास करो।
तुम्हारे द्वारा निर्दिष्ट मान उन अन्य फ़ाइलों में सेट किए गए डिफ़ॉल्ट मानों को ओवरराइड करेंगे।

चलो इसे व्यवहार में आज़माते हैं।

### 2.2. प्रोसेस रिसोर्स और टूल आर्गुमेंट कस्टमाइज़ करें

nf-core मॉड्यूल दो सामान्य प्रकार के कॉन्फ़िगरेशन ओवरराइड को सपोर्ट करते हैं: **रिसोर्स आवंटन** (CPUs, memory, time) और `ext.args` के माध्यम से **टूल आर्गुमेंट**।

कई कमांड-लाइन टूल में ऐसे आर्गुमेंट होते हैं जो पाइपलाइन पैरामीटर के रूप में उजागर करने के लिए पर्याप्त सामान्य नहीं होते।
`ext.args` कन्वेंशन तुम्हें इन आर्गुमेंट को एक config फ़ाइल के माध्यम से अंतर्निहित टूल को पास करने देता है।

तुम्हारी working directory में दी गई `custom.config` फ़ाइल दोनों ओवरराइड प्रदर्शित करती है:

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

पहला ब्लॉक `FASTQC` रिसोर्स आवंटन को ओवरराइड करता है।
डिफ़ॉल्ट रूप से, `FASTQC` `base.config` से `process_medium` लेबल का उपयोग करता है, जो 6 CPUs और 36 GB memory आवंटित करता है; यहां हम इसे 2 CPUs और 4 GB तक सीमित करते हैं।

दूसरा ब्लॉक `ext.args` के माध्यम से `SEQTK_TRIM` को एक अतिरिक्त आर्गुमेंट पास करता है।
`-b 5` फ्लैग `seqtk trimfq` को quality trimming के अलावा प्रत्येक read की शुरुआत से 5 bases trim करने के लिए कहता है।

इस config के साथ पाइपलाइन चलाओ:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-custom -c custom.config
```

??? success "कमांड आउटपुट"

    ```console
    executor >  local (8)
    [95/b32876] NFCORE_DEMO:DEMO:FASTQC (SAMPLE1_PE)     | 3 of 3 ✔
    [17/428668] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE1_PE) | 3 of 3 ✔
    [cf/85991a] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [3c/94a7a0] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

`-c` फ्लैग तुम्हारे config को पाइपलाइन के built-in कॉन्फ़िगरेशन के ऊपर जोड़ता है।

`ext.args` ओवरराइड प्रभावी हुआ या नहीं यह सत्यापित करने के लिए, run आउटपुट से `SEQTK_TRIM` work directory hash खोजो (जैसे `work/17/428668...`) और उसके अंदर `.command.sh` फ़ाइल जांचो:

```bash
cat work/17/428668/.command.sh
```

??? success "कमांड आउटपुट"

    ```console
    #!/usr/bin/env bash -e -u -o pipefail
    printf "%s\n" sample1_R1.fastq.gz sample1_R2.fastq.gz | while read f;
    do
        seqtk \
            trimfq \
            -b 5 \
            $f \
            | gzip --no-name > SAMPLE1_PE_$(basename $f)
    done
    ...
    ```

तुम्हें `seqtk trimfq` कमांड में `-b 5` दिखना चाहिए।

`ext.args` के बारे में एक महत्वपूर्ण बात: यदि किसी मॉड्यूल में पहले से एक डिफ़ॉल्ट मान सेट है, तो तुम्हारा मान उसमें जुड़ने के बजाय उसे **पूरी तरह से बदल** देगा।
उदाहरण के लिए, `FASTQC` में `conf/modules.config` में डिफ़ॉल्ट रूप से `ext.args = '--quiet'` सेट है:

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args   = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}",
        ]
    }
```

यदि तुम `FASTQC` के लिए `ext.args = '--kmers 8'` सेट करते हो, तो `--quiet` फ्लैग अब लागू नहीं होगा।
दोनों को रखने के लिए, `ext.args = '--quiet --kmers 8'` सेट करो।

`ext.args` को ओवरराइड करने से पहले तुम्हें हमेशा किसी मॉड्यूल के डिफ़ॉल्ट कॉन्फ़िगरेशन की जांच करनी चाहिए।

### सारांश

तुम जानते हो कि nf-core पाइपलाइन कॉन्फ़िगरेशन डिफ़ॉल्ट कहां रहते हैं, और एक कस्टम config फ़ाइल के साथ रिसोर्स आवंटन और टूल आर्गुमेंट को कैसे ओवरराइड करें।

### आगे क्या है?

[भाग 3](./03_run_production_pipeline.md) पर जाओ, जहां तुम एक वास्तविक production पाइपलाइन पर जो सीखा है उसे लागू करोगे।

---

## सारांश

इस भाग में तुमने सीखा:

- help प्राप्त करना, पैरामीटर सेट करना, और पैरामीटर तथा इनपुट वैलिडेशन को समझना
- कॉन्फ़िगरेशन फ़ाइलों के माध्यम से रिसोर्स आवंटन और टूल आर्गुमेंट को कस्टमाइज़ करना
