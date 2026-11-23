# भाग 3: एक प्रोडक्शन पाइपलाइन चलाना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[भाग 2](./02_configure_execution.md) में, तुमने nf-core/demo के लिए पैरामीटर सेट करना और कॉन्फ़िगरेशन कस्टमाइज़ करना सीखा।
अब हम जो तुमने सीखा है उसे एक असली प्रोडक्शन पाइपलाइन, nf-core/rnaseq, पर लागू करेंगे।

---

## 1. nf-core/rnaseq को Pull और Run करना

अब तक हमने `nf-core/demo` का उपयोग किया है, जो प्रशिक्षण के लिए बनाई गई एक न्यूनतम पाइपलाइन है।
अब हम एक असली प्रोडक्शन पाइपलाइन pull करेंगे और उसे उसके test profile के साथ चलाएंगे।

`nf-core/rnaseq` पाइपलाइन bulk RNA sequencing विश्लेषण के मुख्य चरण करती है: quality control, adapter trimming, read alignment, और gene-level quantification।
यह अब तक की सबसे अधिक उपयोग की जाने वाली nf-core पाइपलाइन है।

### 1.1. पाइपलाइन Pull करना

इसे डाउनलोड करने के लिए निम्नलिखित कमांड चलाओ।

```bash
nextflow pull nf-core/rnaseq
```

??? success "कमांड आउटपुट"

    ```console
    Checking nf-core/rnaseq ...
     downloaded from https://github.com/nf-core/rnaseq.git - revision: 1f03b53ef7 [master]
    ```

पाइपलाइन अब लोकल में cache हो गई है और चलाने के लिए तैयार है।

### 1.2. Test Profile चलाना

इसे test profile और Docker के साथ चलाओ:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir rnaseq-results
```

??? failure "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [suspicious_dijkstra] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    [-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT  -
    Plus 47 more processes waiting for tasks…

    Execution cancelled -- Finishing pending tasks before exit
    -[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (WT_REP2)'

    Caused by:
      Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB


    Command executed:

      fq lint \
          --disable-validator P001 \
          SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > WT_REP2.fq_lint.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/nfcore-use/work/xx/xxxxxxxxxxxxxxxxxxxxxx

    Container:
      quay.io/biocontainers/fq:0.12.0--h9ee0642_0

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting
    ```

उस error में मुख्य लाइन यह है:

```console
Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB
```

डिफ़ॉल्ट Codespaces मशीन में 8 GB RAM है, जो Docker Desktop के लिए भी सामान्य डिफ़ॉल्ट है।
पाइपलाइन `FQ_LINT` प्रोसेस के लिए 12 GB मांग रही है — जो मशीन दे सकती है उससे अधिक।

वह 12 GB `conf/base.config` में परिभाषित `process_low` resource label से आता है:

```groovy title="conf/base.config"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

एक विकल्प यह होगा कि बड़ी मशीन का उपयोग किया जाए, लेकिन परीक्षण के उद्देश्य से हम जो भी हार्डवेयर उपलब्ध हो उस पर चलाना चाहते हैं।
बेहतर तरीका यह है कि एक custom config फ़ाइल में resource defaults को override किया जाए।

### 1.3. Custom Configuration के साथ फिर से चलाना

हम तुम्हें एक custom config फ़ाइल प्रदान करते हैं जो label-based resource defaults को override करती है।

??? full-code "laptop.config"

    ```groovy title="laptop.config"
    process {
        withLabel: 'process_low' {
            cpus   = 2
            memory = 6.GB
        }
        withLabel: 'process_medium' {
            cpus   = 4
            memory = 6.GB
        }
        withLabel: 'process_high' {
            cpus   = 6
            memory = 6.GB
        }
        withLabel: 'process_high_memory' {
            memory = 6.GB
        }
    }
    ```

[भाग 2](./02_configure_execution.md) ने `withName:` को एक प्रोसेस को नाम से target करने के लिए पेश किया था।
यहाँ हम `withLabel:` का उपयोग करते हैं ताकि एक साथ उन सभी प्रोसेस को target किया जा सके जो एक label साझा करते हैं।

यह फ़ाइल पहले से तुम्हारी working directory में मौजूद है।
Overrides लागू करने के लिए इसे `-c` के साथ पास करो:

```bash
nextflow run nf-core/rnaseq -profile test,docker -c laptop.config --outdir rnaseq-results
```

??? success "कमांड आउटपुट (पाइपलाइन लॉन्च हो रही है)"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [romantic_faraday] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    executor >  local (7)
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FQ_LINT (RAP1_IAA_30M_REP1)   | 3 of 5
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQC (RAP1_IAA_30M_REP1)    | 2 of 5
    ...
    ```

पाइपलाइन अब चल रही है, और तुम कार्यों को एक-एक करके पूरा होते देख सकते हो।
इस न्यूनतम test dataset पर यह 15–20 मिनट में पूरी होगी, कुल मिलाकर 200 से अधिक कार्य निष्पादित करेगी।

असली RNA-seq प्रयोगों में आमतौर पर दर्जनों नमूने होते हैं और वे घंटों या दिनों तक चलते हैं।
Nextflow HPC schedulers (SLURM, PBS, LSF) और cloud platforms (AWS, Google Cloud, Azure) को सपोर्ट करता है, जो कई nodes में काम वितरित करके wall-clock time को काफी कम कर सकते हैं।
हालाँकि, उन वातावरणों को सेटअप करने में काफी जटिलता होती है।

Seqera Platform (Nextflow के निर्माताओं द्वारा विकसित) HPC या cloud infrastructure पर Nextflow पाइपलाइन लॉन्च करने के लिए एक web-based interface प्रदान करता है (या तो तुम्हारा अपना या तुम्हारे लिए प्रबंधित), जिसमें compute और data management क्षमताएं हैं जो बड़े पैमाने पर पाइपलाइन चलाने की प्रक्रिया को सरल बनाती हैं।

!!! tip "सुझाव"

    शैक्षणिक शोधकर्ता [Seqera academic program](https://seqera.io/academic-program/) के माध्यम से Seqera Platform को निःशुल्क उपयोग कर सकते हैं।

### सारांश

तुमने `nf-core/rnaseq` pull किया, देखा कि nf-core resource labels कैसे काम करते हैं, और एक custom config फ़ाइल के साथ उन्हें override करना सीखा।
इससे भी महत्वपूर्ण बात, तुमने देखा कि वास्तविक पैमाने के विश्लेषण के लिए लोकल execution एक शुरुआती बिंदु है, न कि अंतिम लक्ष्य।

### आगे क्या है?

तुमने nf-core पाइपलाइन चलाने की मूल बातें सीख ली हैं।
यहाँ से आगे कहाँ जाना है, इसके लिए [अगले चरण](next_steps.md) देखो।

---

## सारांश

इस भाग में तुमने सीखा:

- एक प्रोडक्शन-स्केल पाइपलाइन (nf-core/rnaseq) को pull और run करना, और उसके डिफ़ॉल्ट resource labels को override करना
