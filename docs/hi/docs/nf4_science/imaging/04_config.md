# भाग 4: कॉन्फ़िगरेशन

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

भाग 1-3 में, हमने सीखा कि Nextflow को कैसे चलाएं, एक nf-core pipeline को कैसे चलाएं, और parameter files और samplesheets के साथ inputs को कैसे मैनेज करें।
अब हम **configuration files** और **profiles** का उपयोग करके विभिन्न computing environments के लिए pipelines को configure करने का तरीका जानेंगे।

## सीखने के उद्देश्य

इस भाग के अंत तक, आप सक्षम होंगे:

- समझना कि Nextflow कई स्रोतों से configuration को कैसे resolve करता है
- containers और testing के लिए nf-core built-in profiles का उपयोग करना
- विभिन्न computing environments के लिए custom profiles बनाना
- process labels का उपयोग करके resource requests को customize करना
- सीमित environments में resource limits को मैनेज करना
- `nextflow config` के साथ resolved configuration को inspect करना

---

## 1. Nextflow configuration को समझना

### 1.1. configuration file क्या है?

Nextflow configuration files का उपयोग **workflow logic** (क्या करना है) को **execution settings** (कैसे और कहाँ करना है) से अलग करने के लिए करता है।

Configuration files निम्नलिखित को नियंत्रित करती हैं:

- Container engines (Docker, Singularity, Conda)
- Compute resources (CPUs, memory, time)
- Execution platforms (local, HPC, cloud)
- Pipeline parameters

### 1.2. Configuration precedence

Nextflow कई स्रोतों से configuration लोड करता है, बाद के स्रोत पहले वाले को override करते हैं:

1. **Pipeline config**: pipeline repository में `nextflow.config`
2. **Directory config**: आपकी वर्तमान working directory में `nextflow.config`
3. **User config**: `~/.nextflow/config`
4. **Command-line**: सीधे पास किए गए parameters और options

यह layered approach आपको pipeline में defaults रखने, user-specific settings के साथ override करने, और command line पर त्वरित समायोजन करने देता है।

### 1.3. हमारा वर्तमान configuration

आइए उस configuration को देखें जिसका हम उपयोग कर रहे हैं:

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

आइए भाग 2 से `docker.enabled = true` लाइन को comment out करें या वापस बदलें, और पता लगाएं कि molkart में एक profile का उपयोग करके हम इसी परिणाम को कैसे प्राप्त कर सकते हैं।

---

## 2. Profiles का उपयोग करना

### 2.1. Profiles क्या हैं?

Profiles configuration के नामित सेट हैं जिन्हें `nextflow run` कमांड के माध्यम से `-profile` flag के साथ सक्रिय किया जा सकता है।
ये config files को edit किए बिना विभिन्न compute scenarios के बीच स्विच करना आसान बनाते हैं।

सभी nf-core pipelines कई default profiles के साथ आते हैं जिनका हम उपयोग कर सकते हैं।

### 2.2. Built-in profiles को inspect करना

आइए उन्हें pipeline codebase से जुड़ी `molkart/nextflow.config` फ़ाइल में inspect करें:

```bash
code molkart/nextflow.config
```

`profiles` block को खोजें:

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

सामान्य container profiles:

- `docker`: Docker containers का उपयोग करें (local development के लिए सबसे आम)
- `singularity`: Singularity/Apptainer का उपयोग करें (HPC पर आम)
- `conda`: Conda environments का उपयोग करें
- `apptainer`: Apptainer containers का उपयोग करें

### 2.3. nextflow.config के बजाय profiles के साथ फिर से चलाना

अब जब हमने अपनी local `nextflow.config` फ़ाइल में docker configuration को disable कर दिया है और profiles को समझ गए हैं, तो आइए `-profile` flag का उपयोग करके pipeline को फिर से चलाएं।

पहले भाग 3 में, हमने अपने custom parameters के साथ एक `params.yaml` फ़ाइल बनाई थी।
अब हम इसे built-in Docker profile के साथ combine कर सकते हैं:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

आइए समझें कि प्रत्येक flag क्या करता है:

- `-profile docker`: molkart के `nextflow.config` से Docker profile को सक्रिय करता है, जो `docker.enabled = true` set करता है
- `-params-file params.yaml`: हमारी YAML फ़ाइल से सभी pipeline parameters लोड करता है
- `-resume`: पिछले runs से cached results का पुनः उपयोग करता है

क्योंकि हम `-resume` का उपयोग कर रहे हैं, Nextflow जांच करेगा कि पिछली बार से कुछ बदला है या नहीं।
यदि parameters, inputs, और code समान हैं, तो सभी tasks cache से पुनर्प्राप्त किए जाएंगे और pipeline लगभग तुरंत पूरा हो जाएगा।

```console title="Output (excerpt)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

ध्यान दें कि सभी processes `cached: 2` या `cached: 1` दिखाते हैं - कुछ भी फिर से execute नहीं किया गया!

### 2.4. Test profiles

Test profiles default input parameters और datafiles को निर्दिष्ट करने के त्वरित तरीके प्रदान करते हैं ताकि आप सत्यापित कर सकें कि pipeline काम करता है।
nf-core pipelines में हमेशा कम से कम दो test profiles शामिल होंगे:

- `test`: त्वरित testing के लिए छोटा dataset और fast parameters
- `test_full`: बड़े data के साथ अधिक व्यापक test

आइए molkart में `test` profile पर करीब से नज़र डालें जो `includeConfig` निर्देश का उपयोग करके शामिल किया गया है:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

इसका मतलब है कि जब भी हम `-profile test` के साथ pipeline चलाते हैं, Nextflow `conf/test.config` से configuration लोड करेगा।

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

ध्यान दें कि इस profile में वही parameters हैं जो हमने पहले अपनी `params.yaml` फ़ाइल में उपयोग किए थे।

आप कॉमा से अलग करके कई profiles को सक्रिय कर सकते हैं।
आइए इसका उपयोग करके अपनी params फ़ाइल की आवश्यकता के बिना अपने pipeline को test करें:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

यह निम्नलिखित को combine करता है:

- `docker`: Docker containers को enable करें
- `test`: Test dataset और parameters का उपयोग करें

Profiles बाएं से दाएं applied होते हैं, इसलिए बाद के profiles पहले वाले को override करते हैं यदि वे समान values set करते हैं।

### निष्कर्ष

nf-core pipelines containers, testing, और विशेष environments के लिए built-in profiles के साथ आते हैं।
आप अपनी आवश्यकता के अनुसार configuration बनाने के लिए कई profiles को combine कर सकते हैं।

### आगे क्या है?

विभिन्न computing environments के लिए अपने custom profiles बनाना सीखें।

---

## 3. Custom profiles बनाना

### 3.1. Local development और HPC पर execution के बीच switching के लिए profiles बनाएं

आइए दो scenarios के लिए custom profiles बनाएं:

1. Docker के साथ local development
2. Slurm scheduler और Singularity के साथ University HPC

अपने `nextflow.config` में निम्नलिखित जोड़ें:

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

अब आप environments के बीच आसानी से switch कर सकते हैं:

```bash
# Local development के लिए
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# HPC के लिए (जब उपलब्ध हो)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! note "नोट"

    हम इस training environment में HPC profile का test नहीं कर सकते क्योंकि हमारे पास Slurm scheduler तक पहुंच नहीं है।
    लेकिन यह दिखाता है कि आप वास्तविक उपयोग के लिए इसे कैसे configure करेंगे।

### 3.2. Configuration को inspect करने के लिए `nextflow config` का उपयोग करें

`nextflow config` कमांड pipeline को चलाए बिना पूरी तरह से resolved configuration दिखाता है।

Default configuration देखें:

```bash
nextflow config ./molkart
```

किसी विशिष्ट profile के साथ configuration देखें:

```bash
nextflow config -profile local_dev ./molkart
```

यह निम्नलिखित के लिए अत्यंत उपयोगी है:

- Configuration issues को debug करना
- समझना कि वास्तव में कौन सी values का उपयोग किया जाएगा
- जांचना कि कई profiles कैसे interact करते हैं

### निष्कर्ष

Custom profiles आपको एक single command-line flag के साथ विभिन्न computing environments के बीच switch करने देते हैं।
चलाने से पहले resolved configuration को inspect करने के लिए `nextflow config` का उपयोग करें।

### आगे क्या है?

nf-core की process label system का उपयोग करके व्यक्तिगत processes के लिए resource requests को customize करना सीखें।

---

## 4. Resource requests को customize करना

### 4.1. nf-core pipelines में process labels को समझना

सरलता के लिए, nf-core pipelines सभी pipelines में resource allocation को standardize करने के लिए [**process labels**](https://www.nextflow.io/docs/latest/reference/process.html#process-label) का उपयोग करते हैं।
प्रत्येक process को low, medium, या high compute resource आवश्यकताओं का वर्णन करने के लिए `process_low`, `process_medium`, या `process_high` जैसे label के साथ tag किया जाता है।
ये labels pipeline की `conf/` directory में स्थित configuration files में से एक में विशिष्ट resource requests में परिवर्तित हो जाते हैं।

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

`task.attempt` multiplier पर ध्यान दें - यह subsequent task retries को अधिक resources request करने की अनुमति देता है, यदि pipeline `process.maxRetries > 1` के साथ set है।

### 4.2. विशिष्ट processes के लिए resources को override करना

fine-grained नियंत्रण के लिए, व्यक्तिगत processes को नाम से target करें:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

यदि हम उपरोक्त override के साथ इस pipeline को चलाने का प्रयास करें, तो `CELLPOSE` process अपने label द्वारा परिभाषित default के बजाय 16 CPUs और 32 GB memory का request करेगा।
यह हमारे वर्तमान environment में pipeline को fail कर देगा क्योंकि हमारे पास इतना RAM उपलब्ध नहीं है।
हम अगले section में इस प्रकार की failures को रोकने का तरीका सीखेंगे।

!!! tip "सुझाव"

    Process names खोजने के लिए, pipeline execution output देखें या `.nextflow.log` की जांच करें।
    Process names `WORKFLOW:SUBWORKFLOW:PROCESS` pattern का पालन करते हैं।

### निष्कर्ष

nf-core pipelines resource allocation को standardize करने के लिए process labels का उपयोग करते हैं।
आप label के द्वारा (कई processes को प्रभावित करता है) या नाम से (एक विशिष्ट process को प्रभावित करता है) resources को override कर सकते हैं।

### आगे क्या है?

GitHub Codespaces जैसे सीमित environments में resource limits को मैनेज करना सीखें।

---

## 5. सीमित environments में resources को मैनेज करना

### 5.1. Resource limits की समस्या

यदि हमने 16 CPUs और 32 GB memory का request करने वाली process के साथ molkart को चलाने का प्रयास किया (जैसा कि section 4.2 में दिखाया गया है), तो यह हमारे वर्तमान environment में fail हो जाएगा क्योंकि हमारे पास इतने resources उपलब्ध नहीं हैं।
बड़े nodes के साथ cluster environment में, ऐसे requests scheduler को submit किए जाएंगे।

GitHub Codespaces जैसे सीमित environments में, limits के बिना, Nextflow उन processes को चलाने से इनकार कर देगा जो उपलब्ध resources से अधिक हैं।

### 5.2. Resource limits सेट करना

`resourceLimits` निर्देश निर्दिष्ट values पर resource requests को cap करता है:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

यह Nextflow को बताता है: "यदि कोई process 2 CPUs या 7 GB memory से अधिक का request करता है, तो इसे इन limits पर cap करें।"

### 5.3. Custom profiles में resource limits जोड़ना

उपयुक्त limits शामिल करने के लिए अपने custom profiles को update करें:

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

    Resource limits को बहुत कम सेट करने से processes fail हो सकते हैं या धीरे चल सकते हैं।
    Pipeline को कम memory-intensive algorithms का उपयोग करने या छोटे chunks में data को process करने की आवश्यकता हो सकती है।

### निष्कर्ष

Process resource requests को cap करके resource-constrained environments में pipelines चलाने के लिए `resourceLimits` का उपयोग करें।
विभिन्न profiles में उनके environment के लिए उपयुक्त अलग-अलग limits हो सकते हैं।

### आगे क्या है?

आपने core Nextflow for Bioimaging training पूरा कर लिया है!

---

## निष्कर्ष

अब आप समझते हैं कि विभिन्न computing environments के लिए Nextflow pipelines को कैसे configure करें।

आपने जो मुख्य कौशल सीखे हैं:

- **Configuration precedence**: Nextflow कई स्रोतों से settings को कैसे resolve करता है
- **nf-core profiles**: Containers, testing, और utilities के लिए built-in profiles का उपयोग करना
- **Custom profiles**: विभिन्न environments के लिए अपने profiles बनाना
- **Process labels**: Label के द्वारा resource requests को समझना और override करना
- **Resource limits**: `resourceLimits` के साथ सीमित environments को मैनेज करना
- **Configuration inspection**: Settings को debug और verify करने के लिए `nextflow config` का उपयोग करना

ये configuration skills किसी भी Nextflow pipeline में transferable हैं और local machines, HPC clusters, और cloud platforms पर workflows को efficiently चलाने में आपकी मदद करेंगी।

### आगे क्या है?

Nextflow for Bioimaging course पूरा करने के लिए बधाई!

अगले कदम:

- Feedback प्रदान करने के लिए course survey भरें
- Workflows develop करने के बारे में अधिक जानने के लिए [Hello Nextflow](../hello_nextflow/index.md) देखें
- nf-core tooling में गहराई से जाने के लिए [Hello nf-core](../hello_nf-core/index.md) explore करें
- [Training collections](../training_collections/index.md) में अन्य courses देखें
