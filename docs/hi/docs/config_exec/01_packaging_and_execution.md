# भाग 1: अपने कंप्यूट वातावरण के अनुसार अनुकूलित करें

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


[Nextflow Run](../nextflow_run/index.md) में, तुमने एक पाइपलाइन के इनपुट, पैरामीटर और आउटपुट कॉन्फ़िगर किए।
यह कोर्स तस्वीर के दूसरे हिस्से को कवर करता है: वर्कफ़्लो कोड बदले बिना, पाइपलाइन के execution को किसी भी कंप्यूट वातावरण के अनुसार अनुकूलित करना।

!!! example "परिदृश्य"

    तुमने अपनी पाइपलाइन Docker का उपयोग करके अपने लैपटॉप पर develop और test की।
    अब तुम्हें इसे आगे देना है: एक सहयोगी के पास केवल Conda सेटअप है, और तुम्हारे संस्थान का HPC cluster अपने scheduler और अपनी resource limits के साथ jobs चलाने की उम्मीद करता है।
    इनमें से किसी के लिए भी पाइपलाइन को फिर से लिखने की ज़रूरत नहीं होनी चाहिए।

वही पाइपलाइन कोड इन सभी जगहों पर चल सकता है, क्योंकि इनमें से कुछ भी वर्कफ़्लो में बेक-इन नहीं है।
Software packaging, execution platform, और resource allocation सभी configuration के ज़रिए नियंत्रित होते हैं, कोड के ऊपर layered होते हैं — और यही यह कोर्स कवर करता है: कोड नहीं, config बदलकर उसी पाइपलाइन को नए वातावरण के अनुसार कैसे अनुकूलित करें।

---

## 1. Software packaging technology चुनें

[Nextflow Run](../nextflow_run/index.md) में, तुमने देखा कि `nextflow.config` में Docker के विकल्प के रूप में एक `conda` profile पहले से सेटअप था।
यहाँ तुम वही switch खुद बनाओगे, और देखोगे कि किसी process को Conda के साथ वास्तव में उपयोग करने योग्य बनाने के लिए क्या चाहिए।

### 1.1. Docker बंद करें और Conda चालू करें

`docker.enabled` को `false` पर switch करो और Conda को enable करने वाला एक निर्देश जोड़ो।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

इससे Nextflow किसी भी process के लिए Conda वातावरण बना और उपयोग कर सकता है जिसमें Conda package निर्दिष्ट हो।
`cowpy` process में अभी तक कोई नहीं है, तो चलो एक जोड़ते हैं, पूरी तरह config से।

### 1.2. Config के ज़रिए Conda package जोड़ें

एक `conda` निर्देश process definition में ही सेट किया जा सकता है, उसी तरह जैसे `container` पहले से `modules/cowpy.nf` में है, लेकिन यह ज़रूरी नहीं है: `withName` तुम्हें इसे config से सेट करने देता है, केवल `cowpy` process तक scoped।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="6" hl_lines="3-5"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
    }
    ```

यह पाइपलाइन कोड में पहले से मौजूद `container` निर्देश को replace नहीं करता, यह उस कोड को बिल्कुल छुए बिना उसके साथ एक विकल्प जोड़ता है।

!!! tip "सुझाव"

    [Seqera Containers](https://seqera.io/containers/) search किसी tool के लिए Conda package URI देखने का एक सुविधाजनक तरीका है, भले ही तुम उससे container बनाने की योजना न बना रहे हो।

### 1.3. यह verify करने के लिए वर्कफ़्लो चलाओ कि वह Conda उपयोग कर सकता है

```bash
nextflow run main.nf --batch conda
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_pike] revision: c3c85dec78

    executor >  local (8)
    [6d/d48030] sayHello (2)       | 3 of 3 ✔
    [f5/7a9d76] convertToUpper (1) | 3 of 3 ✔
    [1c/79b693] collectGreetings   | 1 of 1 ✔
    Creating env using conda: conda-forge::cowpy==1.1.5 [cache /workspaces/training/config-exec/work/conda/env-898314d566668b6587ad714ae06b8520]
    [bb/64b67c] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-conda-output.txt

      batch_report: full_pipeline/conda-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-conda-output.txt
    ```

यह Docker के साथ चलाने जैसा ही आउटपुट देता है, भले ही पर्दे के पीछे mechanics अलग हों: Nextflow container image pull करने की बजाय Conda package retrieve करता है और उससे एक वातावरण बनाता है।

!!! info "जानकारी"

    पहली बार एक नया Conda वातावरण बनाने में container pull करने से थोड़ा अधिक समय लग सकता है, लेकिन यहाँ उपयोग किया गया package छोटा है इसलिए यह जल्दी होना चाहिए।

अब इस कोर्स के बाकी हिस्से के लिए Docker पर वापस switch करो।

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

??? tip "Docker और Conda को मिलाना"

    क्योंकि ये settings प्रति process scoped हैं, तुम इन्हें मिला सकते हो: कुछ processes Docker उपयोग करती हैं, अन्य Conda, इस पर निर्भर करते हुए कि प्रत्येक tool के लिए क्या उपलब्ध है।
    अगर एक ही process के लिए `container` निर्देश (पाइपलाइन कोड में) और `conda` निर्देश (यहाँ, config से) दोनों सेट हों और दोनों packaging systems enabled हों, तो Nextflow containers को प्राथमिकता देता है।

### सारांश

तुम जानते हो कि किसी process को कौन सी software packaging technology उपयोग करनी चाहिए यह कैसे configure करें, और Docker और Conda के बीच कैसे switch करें।

### आगे क्या है?

जानो कि Nextflow तुम्हारे tasks को वास्तव में चलाने के लिए execution platform कैसे बदलें।

---

## 2. Execution platform चुनें

अब तक तुमने जो भी पाइपलाइन चलाई है वह local executor का उपयोग करती है: प्रत्येक कार्य Nextflow के समान machine पर चलता है।
Nextflow उपलब्ध CPUs और memory की जाँच करता है, और tasks को तब तक रोकता है जब तक पर्याप्त resources मुक्त न हो जाएं।

Local executor सुविधाजनक है, लेकिन यह एक machine से आगे scale नहीं होता।
Nextflow [कई अन्य execution backends](https://nextflow.io/docs/latest/executor.html) को support करता है, जिनमें HPC schedulers (Slurm, LSF, SGE, PBS, और अन्य) और cloud platforms (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes, और अधिक) शामिल हैं।

### 2.1. एक अलग backend को target करें

Executor `executor` नामक एक process निर्देश द्वारा सेट किया जाता है।
Default रूप से यह `local` है, इसलिए निम्नलिखित implied है:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

किसी अलग backend को target करने के लिए, निर्देश को अपने मनचाहे executor पर सेट करो।

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "चेतावनी"

    Training वातावरण किसी HPC cluster से connected नहीं है, इसलिए यह कुछ ऐसा नहीं है जो तुम यहाँ चला सकते हो।

### 2.2. Backend-specific syntax को abstract किया गया है

अधिकांश HPC platforms को job submissions में resource requests निर्दिष्ट करने की आवश्यकता होती है, जैसे CPUs, memory, और queue का नाम, अपनी खुद की syntax का उपयोग करके।
`my-science-work` नामक queue पर 8 CPUs और 4 GB RAM के लिए वही request scheduler के आधार पर बिल्कुल अलग दिखती है।

??? abstract "उदाहरण"

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=8
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=8
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Nextflow यह सब abstract कर देता है: तुम एक बार standardized properties जैसे `cpus`, `memory`, और `queue` निर्दिष्ट करते हो (पूरी सूची के लिए [process directives](https://nextflow.io/docs/latest/reference/process.html#process-directives) देखो), और Nextflow उन्हें runtime पर उचित backend-specific scripts में translate करता है।

### 2.3. देखो कि Nextflow वास्तव में क्या चलाता है

वह translation केवल एक config-file की सुविधा नहीं है: यह किसी ठोस चीज़ पर आधारित है जिसे तुम अभी inspect कर सकते हो, local executor के साथ भी।
[Nextflow Run, section 1.3](../nextflow_run/01_run_nextflow.md#13-explore-the-work-directory) में, तुमने `work/` के अंतर्गत एक task directory के अंदर देखा और `.command.sh` पाया, वह exact कमांड जो Nextflow ने चलाई।
उसी directory में एक फ़ाइल भी है जिसे तुमने अभी तक नहीं देखा: `.command.run`।

```bash
cat work/0a/0df4a1*/.command.run
```

??? success "कमांड आउटपुट (अंश)"

    ```console
    #!/bin/bash
    ### ---
    ### name: 'convertToUpper (3)'
    ### container: 'null'
    ### outputs:
    ### - 'UPPER-Bonjour-output.txt'
    ### ...
    set -e
    set -u
    ...
    nxf_launch() {
        /bin/bash -ue /workspaces/training/nextflow-run/work/0a/0df4a1028c2001758b1841cff92fc7/.command.sh
    }
    ...
    ```

`.command.run` वह असली script है जिसे Nextflow execution के लिए सौंपता है।
यह `.command.sh` को वह सब कुछ के साथ wrap करता है जो इसे वास्तव में चलाने के लिए चाहिए: वातावरण सेटअप, input/output staging, और Nextflow को परिणाम वापस report करना।
`local` executor के साथ, Nextflow बस इस script को उसी machine पर चलाता है।

यही वह है जो बदलता है जब तुम एक अलग `executor` सेट करते हो।
Slurm या PBS जैसे HPC scheduler के लिए, Nextflow उसी तरह का wrapper script generate करता है, उसमें scheduler-specific header जोड़ता है जो तुमने [2.2](#22-backend-specific-syntax-is-abstracted-away) में देखा (तुम्हारी `cpus`, `memory`, और `queue` settings से translate किया गया), और परिणाम को उस scheduler की अपनी submission कमांड को सौंपता है, उदाहरण के लिए Slurm के लिए `sbatch`।
वहाँ से, Nextflow local process को सीधे देखने की बजाय job status के लिए scheduler को poll करता है।
Cloud batch backends थोड़े अलग तरीके से काम करते हैं, क्योंकि वे submission कमांड की बजाय API calls द्वारा driven होते हैं, लेकिन वही underlying विचार लागू होता है: वही task script चलती है, केवल इसे कैसे launch और track किया जाता है यह बदलता है।

### सारांश

तुम जानते हो कि अलग-अलग compute infrastructure को target करने के लिए executor कैसे बदलें, कि Nextflow backend-specific submission syntax को abstract कर देता है, और जब कोई कार्य किसी अलग backend पर चलता है तो पर्दे के पीछे वास्तव में क्या होता है।

### आगे क्या है?

[भाग 2](./02_resources_and_retries.md) पर जाओ, जहाँ तुम सीखोगे कि compute resources को कैसे profile और allocate करें, और retries के साथ task failures को कैसे handle करें।

---

## सारांश

इस भाग में तुमने सीखा:

- Docker और Conda के बीच software packaging technology switch करना
- Process definition में `conda` निर्देश जोड़ना
- `executor` निर्देश के साथ execution platform बदलना
- Nextflow किसी कार्य के लिए वास्तव में क्या generate और चलाता है यह inspect करना, और यह executors के बीच कैसे बदलता है
