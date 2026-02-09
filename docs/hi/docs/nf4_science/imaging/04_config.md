# भाग 4: कॉन्फ़िगरेशन

भाग 1-3 में, हमने सीखा कि Nextflow कैसे चलाएं, nf-core पाइपलाइन कैसे चलाएं, और पैरामीटर फ़ाइलों और samplesheets के साथ इनपुट कैसे मैनेज करें।
अब हम जानेंगे कि **कॉन्फ़िगरेशन फ़ाइलों** और **प्रोफ़ाइलों** का उपयोग करके विभिन्न कंप्यूटिंग वातावरणों के लिए पाइपलाइनों को कैसे कॉन्फ़िगर करें।

## सीखने के उद्देश्य

इस भाग के अंत तक, तुम यह कर पाओगे:

- समझना कि Nextflow कई स्रोतों से कॉन्फ़िगरेशन कैसे resolve करता है
- कंटेनरों और टेस्टिंग के लिए nf-core बिल्ट-इन प्रोफ़ाइलों का उपयोग करना
- विभिन्न कंप्यूटिंग वातावरणों के लिए कस्टम प्रोफ़ाइल बनाना
- प्रोसेस लेबल का उपयोग करके रिसोर्स रिक्वेस्ट को कस्टमाइज़ करना
- सीमित वातावरणों में रिसोर्स लिमिट मैनेज करना
- `nextflow config` के साथ resolved कॉन्फ़िगरेशन को inspect करना

---

## 1. Nextflow कॉन्फ़िगरेशन को समझना

### 1.1. कॉन्फ़िगरेशन फ़ाइल क्या है?

Nextflow कॉन्फ़िगरेशन फ़ाइलों का उपयोग **वर्कफ़्लो लॉजिक** (क्या करना है) को **एक्ज़ीक्यूशन सेटिंग्स** (कैसे और कहाँ करना है) से अलग करने के लिए करता है।

कॉन्फ़िगरेशन फ़ाइलें नियंत्रित करती हैं:

- कंटेनर इंजन (Docker, Singularity, Conda)
- कंप्यूट रिसोर्स (CPUs, मेमोरी, समय)
- एक्ज़ीक्यूशन प्लेटफ़ॉर्म (लोकल, HPC, क्लाउड)
- पाइपलाइन पैरामीटर

### 1.2. कॉन्फ़िगरेशन precedence

Nextflow कई स्रोतों से कॉन्फ़िगरेशन लोड करता है, बाद के स्रोत पहले वाले को override करते हैं:

1. **पाइपलाइन config**: पाइपलाइन रिपॉज़िटरी में `nextflow.config`
2. **डायरेक्टरी config**: तुम्हारी वर्तमान वर्किंग डायरेक्टरी में `nextflow.config`
3. **यूज़र config**: `~/.nextflow/config`
4. **कमांड-लाइन**: सीधे पास किए गए पैरामीटर और ऑप्शन

यह लेयर्ड दृष्टिकोण तुम्हें पाइपलाइन में डिफ़ॉल्ट रखने, यूज़र-विशिष्ट सेटिंग्स के साथ override करने, और कमांड लाइन पर त्वरित समायोजन करने देता है।

### 1.3. हमारा वर्तमान कॉन्फ़िगरेशन

आओ देखें कि हम किस कॉन्फ़िगरेशन का उपयोग कर रहे हैं:

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

आओ भाग 2 से `docker.enabled = true` लाइन को कमेंट आउट करें या वापस बदलें, और पता लगाएं कि हम molkart में प्रोफ़ाइल का उपयोग करके वही परिणाम कैसे प्राप्त कर सकते हैं।

---

## 2. प्रोफ़ाइलों का उपयोग करना

### 2.1. प्रोफ़ाइल क्या हैं?

प्रोफ़ाइल कॉन्फ़िगरेशन के नामित सेट हैं जिन्हें `nextflow run` कमांड के माध्यम से `-profile` फ़्लैग के साथ सक्रिय किया जा सकता है।
ये config फ़ाइलों को एडिट किए बिना विभिन्न कंप्यूट परिदृश्यों के बीच स्विच करना आसान बनाते हैं।

सभी nf-core पाइपलाइनों में कई डिफ़ॉल्ट प्रोफ़ाइल होते हैं जिनका हम उपयोग कर सकते हैं।

### 2.2. बिल्ट-इन प्रोफ़ाइलों को inspect करना

आओ पाइपलाइन कोडबेस से जुड़ी `molkart/nextflow.config` फ़ाइल में उन्हें inspect करें:

```bash
code molkart/nextflow.config
```

`profiles` ब्लॉक खोजें:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
    docker {
        docker.enabled          = true
        singularity.enabled     = false
        conda.enabled           = false
    }
    singularity {
        singularity.enabled     = true
        singularity.autoMounts  = true
        docker.enabled          = false
    }
    conda {
        conda.enabled           = true
        docker.enabled          = false
        conda.channels          = ['conda-forge', 'bioconda']
    }
}
```

सामान्य कंटेनर प्रोफ़ाइल:

- `docker`: Docker कंटेनर का उपयोग करें (लोकल डेवलपमेंट के लिए सबसे आम)
- `singularity`: Singularity/Apptainer का उपयोग करें (HPC पर आम)
- `conda`: Conda वातावरण का उपयोग करें
- `apptainer`: Apptainer कंटेनर का उपयोग करें

### 2.3. nextflow.config के बजाय प्रोफ़ाइलों के साथ फिर से चलाना

अब जब हमने अपनी लोकल `nextflow.config` फ़ाइल में docker कॉन्फ़िगरेशन को डिसेबल कर दिया है और प्रोफ़ाइलों को समझ लिया है, तो आओ `-profile` फ़्लैग का उपयोग करके पाइपलाइन को फिर से चलाएं।

पहले भाग 3 में, हमने अपने कस्टम पैरामीटर के साथ एक `params.yaml` फ़ाइल बनाई थी।
अब हम इसे बिल्ट-इन Docker प्रोफ़ाइल के साथ जोड़ सकते हैं:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

आओ समझें कि प्रत्येक फ़्लैग क्या करता है:

- `-profile docker`: molkart के `nextflow.config` से Docker प्रोफ़ाइल को सक्रिय करता है, जो `docker.enabled = true` सेट करता है
- `-params-file params.yaml`: हमारी YAML फ़ाइल से सभी पाइपलाइन पैरामीटर लोड करता है
- `-resume`: पिछले रन से cached परिणामों का पुन: उपयोग करता है

क्योंकि हम `-resume` का उपयोग कर रहे हैं, Nextflow जांचेगा कि पिछले रन के बाद से कुछ बदला है या नहीं।
यदि पैरामीटर, इनपुट और कोड समान हैं, तो सभी कार्य cache से पुनर्प्राप्त किए जाएंगे और पाइपलाइन लगभग तुरंत पूरी हो जाएगी।

```console title="Output (excerpt)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

ध्यान दो कि सभी प्रोसेस `cached: 2` या `cached: 1` दिखाते हैं - कुछ भी फिर से execute नहीं हुआ!

### 2.4. टेस्ट प्रोफ़ाइल

टेस्ट प्रोफ़ाइल डिफ़ॉल्ट इनपुट पैरामीटर और डेटाफ़ाइलों को निर्दिष्ट करने के त्वरित तरीके प्रदान करते हैं ताकि तुम सत्यापित कर सको कि पाइपलाइन काम करती है।
nf-core पाइपलाइनों में हमेशा कम से कम दो टेस्ट प्रोफ़ाइल शामिल होंगे:

- `test`: त्वरित टेस्टिंग के लिए छोटा डेटासेट और तेज़ पैरामीटर
- `test_full`: बड़े डेटा के साथ अधिक व्यापक टेस्ट

आओ molkart में `test` प्रोफ़ाइल पर करीब से नज़र डालें जो `includeConfig` निर्देश का उपयोग करके शामिल किया गया है:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

इसका मतलब है कि जब भी हम `-profile test` के साथ पाइपलाइन चलाते हैं, Nextflow `conf/test.config` से कॉन्फ़िगरेशन लोड करेगा।

```groovy title="molkart/conf/test.config (excerpt)"
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    input = 'https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv'
    mindagap_tilesize = 90
    mindagap_boxsize = 7
    mindagap_loopnum = 100
    clahe_pyramid_tile = 368
    segmentation_method = "mesmer,cellpose,stardist"
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

ध्यान दो कि इस प्रोफ़ाइल में वही पैरामीटर हैं जो हमने पहले अपनी `params.yaml` फ़ाइल में उपयोग किए थे।

तुम कॉमा से अलग करके कई प्रोफ़ाइल सक्रिय कर सकते हो।
आओ इसका उपयोग करके अपनी पाइपलाइन को params फ़ाइल की आवश्यकता के बिना टेस्ट करें:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

यह जोड़ता है:

- `docker`: Docker कंटेनर सक्षम करें
- `test`: टेस्ट डेटासेट और पैरामीटर का उपयोग करें

प्रोफ़ाइल बाएं से दाएं लागू होते हैं, इसलिए बाद के प्रोफ़ाइल पहले वाले को override करते हैं यदि वे समान मान सेट करते हैं।

### सारांश

nf-core पाइपलाइनों में कंटेनरों, टेस्टिंग और विशेष वातावरणों के लिए बिल्ट-इन प्रोफ़ाइल होते हैं।
तुम अपनी आवश्यक कॉन्फ़िगरेशन बनाने के लिए कई प्रोफ़ाइलों को जोड़ सकते हो।

### आगे क्या है?

विभिन्न कंप्यूटिंग वातावरणों के लिए अपने स्वयं के कस्टम प्रोफ़ाइल बनाना सीखो।

---

## 3. कस्टम प्रोफ़ाइल बनाना

### 3.1. लोकल डेवलपमेंट और HPC पर एक्ज़ीक्यूशन के बीच स्विच करने के लिए प्रोफ़ाइल बनाएं

आओ दो परिदृश्यों के लिए कस्टम प्रोफ़ाइल बनाएं:

1. Docker के साथ लोकल डेवलपमेंट
2. Slurm शेड्यूलर और Singularity के साथ यूनिवर्सिटी HPC

अपनी `nextflow.config` में निम्नलिखित जोड़ें:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'standard_queue'
        singularity.cacheDir = '/shared/containers'
    }
}
```

अब तुम आसानी से वातावरणों के बीच स्विच कर सकते हो:

```bash
# लोकल डेवलपमेंट के लिए
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# HPC के लिए (जब उपलब्ध हो)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! note "नोट"

    हम इस प्रशिक्षण वातावरण में HPC प्रोफ़ाइल का टेस्ट नहीं कर सकते क्योंकि हमारे पास Slurm शेड्यूलर तक पहुंच नहीं है।
    लेकिन यह दिखाता है कि तुम इसे वास्तविक दुनिया के उपयोग के लिए कैसे कॉन्फ़िगर करोगे।

### 3.2. कॉन्फ़िगरेशन को inspect करने के लिए `nextflow config` का उपयोग करें

`nextflow config` कमांड पाइपलाइन चलाए बिना पूरी तरह से resolved कॉन्फ़िगरेशन दिखाता है।

डिफ़ॉल्ट कॉन्फ़िगरेशन देखें:

```bash
nextflow config ./molkart
```

किसी विशिष्ट प्रोफ़ाइल के साथ कॉन्फ़िगरेशन देखें:

```bash
nextflow config -profile local_dev ./molkart
```

यह इसके लिए बेहद उपयोगी है:

- कॉन्फ़िगरेशन समस्याओं को डीबग करना
- समझना कि वास्तव में कौन से मान उपयोग किए जाएंगे
- जांचना कि कई प्रोफ़ाइल कैसे इंटरैक्ट करते हैं

### सारांश

कस्टम प्रोफ़ाइल तुम्हें एक कमांड-लाइन फ़्लैग के साथ विभिन्न कंप्यूटिंग वातावरणों के बीच स्विच करने देते हैं।
चलाने से पहले resolved कॉन्फ़िगरेशन को inspect करने के लिए `nextflow config` का उपयोग करो।

### आगे क्या है?

nf-core के प्रोसेस लेबल सिस्टम का उपयोग करके व्यक्तिगत प्रोसेस के लिए रिसोर्स रिक्वेस्ट को कस्टमाइज़ करना सीखो।

---

## 4. रिसोर्स रिक्वेस्ट को कस्टमाइज़ करना

### 4.1. nf-core पाइपलाइनों में प्रोसेस लेबल को समझना

सरलता के लिए, nf-core पाइपलाइनें सभी पाइपलाइनों में रिसोर्स आवंटन को मानकीकृत करने के लिए [**प्रोसेस लेबल**](https://www.nextflow.io/docs/latest/reference/process.html#process-label) का उपयोग करती हैं।
प्रत्येक प्रोसेस को क्रमशः कम, मध्यम या उच्च कंप्यूट रिसोर्स आवश्यकताओं का वर्णन करने के लिए `process_low`, `process_medium`, या `process_high` जैसे लेबल के साथ टैग किया जाता है।
ये लेबल पाइपलाइन की `conf/` डायरेक्टरी में स्थित कॉन्फ़िगरेशन फ़ाइलों में से एक में विशिष्ट रिसोर्स रिक्वेस्ट में परिवर्तित हो जाते हैं।

```groovy title="molkart/conf/base.config (excerpt)"
process {
    cpus   = { 1      * task.attempt }
    memory = { 6.GB   * task.attempt }
    time   = { 4.h    * task.attempt }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_single {
        cpus   = { 1                   }
        memory = { 6.GB * task.attempt }
        time   = { 4.h  * task.attempt }
    }
    withLabel:process_low {
        cpus   = { 2     * task.attempt }
        memory = { 12.GB * task.attempt }
        time   = { 4.h   * task.attempt }
    }
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 8.h   * task.attempt }
    }
    withLabel:process_high {
        cpus   = { 12    * task.attempt }
        memory = { 72.GB * task.attempt }
        time   = { 16.h  * task.attempt }
    }
}
```

`task.attempt` गुणक पर ध्यान दो - यह बाद के कार्य पुनः प्रयासों को अधिक रिसोर्स का अनुरोध करने की अनुमति देता है, यदि पाइपलाइन `process.maxRetries > 1` के साथ सेट है।

### 4.2. विशिष्ट प्रोसेस के लिए रिसोर्स को override करना

बारीक नियंत्रण के लिए, नाम से व्यक्तिगत प्रोसेस को लक्षित करो:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

यदि हम उपरोक्त override के साथ इस पाइपलाइन को चलाने की कोशिश करते हैं, तो `CELLPOSE` प्रोसेस अपने लेबल द्वारा परिभाषित डिफ़ॉल्ट के बजाय 16 CPUs और 32 GB मेमोरी का अनुरोध करेगा।
यह हमारे वर्तमान वातावरण में पाइपलाइन को विफल कर देगा क्योंकि हमारे पास इतनी RAM उपलब्ध नहीं है।
हम अगले सेक्शन में सीखेंगे कि इस प्रकार की विफलताओं को कैसे रोकें।

!!! tip "सुझाव"

    प्रोसेस नाम खोजने के लिए, पाइपलाइन एक्ज़ीक्यूशन आउटपुट देखो या `.nextflow.log` जांचो।
    प्रोसेस नाम `WORKFLOW:SUBWORKFLOW:PROCESS` पैटर्न का पालन करते हैं।

### सारांश

nf-core पाइपलाइनें रिसोर्स आवंटन को मानकीकृत करने के लिए प्रोसेस लेबल का उपयोग करती हैं।
तुम लेबल द्वारा (कई प्रोसेस को प्रभावित करता है) या नाम द्वारा (एक विशिष्ट प्रोसेस को प्रभावित करता है) रिसोर्स को override कर सकते हो।

### आगे क्या है?

GitHub Codespaces जैसे सीमित वातावरणों में रिसोर्स लिमिट को मैनेज करना सीखो।

---

## 5. सीमित वातावरणों में रिसोर्स मैनेज करना

### 5.1. रिसोर्स लिमिट की समस्या

यदि हमने 16 CPUs और 32 GB मेमोरी का अनुरोध करने वाली प्रोसेस के साथ molkart चलाने की कोशिश की (जैसा कि सेक्शन 4.2 में दिखाया गया है), तो यह हमारे वर्तमान वातावरण में विफल हो जाएगा क्योंकि हमारे पास इतने रिसोर्स उपलब्ध नहीं हैं।
बड़े नोड्स वाले क्लस्टर वातावरण में, ऐसे अनुरोध शेड्यूलर को सबमिट किए जाएंगे।

GitHub Codespaces जैसे सीमित वातावरणों में, लिमिट के बिना, Nextflow उन प्रोसेस को चलाने से इनकार कर देगा जो उपलब्ध रिसोर्स से अधिक हैं।

### 5.2. रिसोर्स लिमिट सेट करना

`resourceLimits` निर्देश निर्दिष्ट मानों पर रिसोर्स रिक्वेस्ट को cap करता है:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

यह Nextflow को बताता है: "यदि कोई प्रोसेस 2 CPUs या 7 GB मेमोरी से अधिक का अनुरोध करता है, तो इसे इन लिमिट पर cap करें।"

### 5.3. कस्टम प्रोफ़ाइलों में रिसोर्स लिमिट जोड़ना

उपयुक्त लिमिट शामिल करने के लिए अपने कस्टम प्रोफ़ाइलों को अपडेट करो:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
        process.resourceLimits = [
            cpus: 2,
            memory: 7.GB
        ]
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'batch'
        process.resourceLimits = [
            cpus: 32,
            memory: 128.GB,
            time: 24.h
        ]
    }
}
```

!!! warning "चेतावनी"

    रिसोर्स लिमिट बहुत कम सेट करने से प्रोसेस विफल हो सकते हैं या धीरे चल सकते हैं।
    पाइपलाइन को कम मेमोरी-गहन एल्गोरिदम का उपयोग करने या छोटे chunks में डेटा प्रोसेस करने की आवश्यकता हो सकती है।

### सारांश

प्रोसेस रिसोर्स रिक्वेस्ट को cap करके रिसोर्स-सीमित वातावरणों में पाइपलाइनों को चलाने के लिए `resourceLimits` का उपयोग करो।
विभिन्न प्रोफ़ाइलों में उनके वातावरण के लिए उपयुक्त विभिन्न लिमिट हो सकती हैं।

### आगे क्या है?

तुमने कोर Nextflow for Bioimaging प्रशिक्षण पूरा कर लिया है!

---

## निष्कर्ष

अब तुम समझते हो कि विभिन्न कंप्यूटिंग वातावरणों के लिए Nextflow पाइपलाइनों को कैसे कॉन्फ़िगर करें।

तुमने जो मुख्य कौशल सीखे हैं:

- **कॉन्फ़िगरेशन precedence**: Nextflow कई स्रोतों से सेटिंग्स कैसे resolve करता है
- **nf-core प्रोफ़ाइल**: कंटेनरों, टेस्टिंग और उपयोगिताओं के लिए बिल्ट-इन प्रोफ़ाइलों का उपयोग करना
- **कस्टम प्रोफ़ाइल**: विभिन्न वातावरणों के लिए अपने स्वयं के प्रोफ़ाइल बनाना
- **प्रोसेस लेबल**: लेबल द्वारा रिसोर्स रिक्वेस्ट को समझना और override करना
- **रिसोर्स लिमिट**: `resourceLimits` के साथ सीमित वातावरणों को मैनेज करना
- **कॉन्फ़िगरेशन inspection**: सेटिंग्स को डीबग और सत्यापित करने के लिए `nextflow config` का उपयोग करना

ये कॉन्फ़िगरेशन कौशल किसी भी Nextflow पाइपलाइन में स्थानांतरणीय हैं और तुम्हें लोकल मशीनों, HPC क्लस्टरों और क्लाउड प्लेटफ़ॉर्मों पर वर्कफ़्लो को कुशलता से चलाने में मदद करेंगे।

### आगे क्या है?

Nextflow for Bioimaging कोर्स पूरा करने पर बधाई!

अगले कदम:

- फ़ीडबैक देने के लिए कोर्स सर्वे भरो
- वर्कफ़्लो विकसित करने के बारे में अधिक जानने के लिए [Hello Nextflow](../hello_nextflow/index.md) देखो
- nf-core टूलिंग में गहराई से जाने के लिए [Hello nf-core](../hello_nf-core/index.md) एक्सप्लोर करो
- [प्रशिक्षण संग्रह](../training_collections/index.md) में अन्य कोर्स ब्राउज़ करो
