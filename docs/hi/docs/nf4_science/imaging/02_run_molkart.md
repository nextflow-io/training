# भाग 2: nf-core/molkart चलाना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

भाग 1 में, हमने Nextflow निष्पादन की मूल बातें समझने के लिए एक सरल Hello World वर्कफ़्लो चलाया।
अब हम एक वास्तविक बायोइमेजिंग पाइपलाइन चलाने जा रहे हैं: **nf-core/molkart**।

यह पाइपलाइन Resolve Bioscience से Molecular Cartography spatial transcriptomics डेटा को प्रोसेस करती है।
हालांकि, यहां आप जो Nextflow पैटर्न सीखेंगे वे किसी भी nf-core पाइपलाइन या प्रोडक्शन वर्कफ़्लो पर लागू होते हैं।

## 1. nf-core पाइपलाइनों को समझना

पाइपलाइन चलाने से पहले, आइए समझें कि nf-core क्या है और वर्कफ़्लो चलाने के लिए यह क्यों महत्वपूर्ण है।

### 1.1. nf-core क्या है?

[nf-core](https://nf-co.re/) उच्च-गुणवत्ता वाली Nextflow पाइपलाइनों का एक समुदाय-संचालित संग्रह है।
सभी nf-core पाइपलाइनें समान संरचना और परंपराओं का पालन करती हैं, जिसका मतलब है कि एक बार जब आप एक को चलाना सीख लेते हैं, तो आप उनमें से किसी को भी चला सकते हैं।

nf-core पाइपलाइनों की मुख्य विशेषताएं:

- **मानकीकृत संरचना**: सभी पाइपलाइनों में सुसंगत पैरामीटर नाम और उपयोग पैटर्न हैं
- **अंतर्निहित टेस्ट डेटा**: प्रत्येक पाइपलाइन में त्वरित सत्यापन के लिए टेस्ट प्रोफाइल शामिल हैं
- **व्यापक दस्तावेज़ीकरण**: विस्तृत उपयोग निर्देश और पैरामीटर विवरण
- **गुणवत्ता नियंत्रण**: MultiQC का उपयोग करके स्वचालित QC रिपोर्ट
- **कंटेनर समर्थन**: पुनरुत्पादकता के लिए पूर्व-निर्मित कंटेनर

!!! tip "nf-core के बारे में अधिक जानना चाहते हैं?"

    nf-core पाइपलाइन विकास के गहन परिचय के लिए, [Hello nf-core](../../hello_nf-core/index.md) प्रशिक्षण पाठ्यक्रम देखें।
    यह शुरुआत से nf-core पाइपलाइनों को बनाने और कस्टमाइज़ करने को कवर करता है।

### 1.2. molkart पाइपलाइन

![nf-core/molkart पाइपलाइन](img/molkart.png)

[nf-core/molkart](https://nf-co.re/molkart) पाइपलाइन कई चरणों में spatial transcriptomics इमेजिंग डेटा को प्रोसेस करती है:

1. **इमेज प्रीप्रोसेसिंग**: ग्रिड पैटर्न फिलिंग और वैकल्पिक कंट्रास्ट एन्हांसमेंट
2. **सेल सेगमेंटेशन**: कई एल्गोरिदम विकल्प (Cellpose, Mesmer, ilastik, Stardist)
3. **स्पॉट असाइनमेंट**: सेगमेंटेड सेल्स को ट्रांसक्रिप्ट स्पॉट असाइन करना
4. **गुणवत्ता नियंत्रण**: व्यापक QC रिपोर्ट जेनरेट करना

मुख्य आउटपुट हैं:

- सेल-बाय-ट्रांसक्रिप्ट काउंट टेबल
- सेगमेंटेशन मास्क
- MultiQC गुणवत्ता नियंत्रण रिपोर्ट

---

## 2. टेस्ट डेटा के साथ molkart चलाना

शुरू करने से पहले, आइए molkart रिपॉजिटरी को स्थानीय रूप से क्लोन करें ताकि हम इसके कोड का निरीक्षण कर सकें:

```bash
cd /workspaces/training/nf4-science/imaging
git clone --branch 1.2.0 --depth 1 https://github.com/nf-core/molkart
```

यह एक `molkart/` डायरेक्टरी बनाता है जिसमें संपूर्ण पाइपलाइन सोर्स कोड होता है।

!!! note "हम स्थानीय रूप से क्यों क्लोन कर रहे हैं?"

    आमतौर पर, आप nf-core पाइपलाइनों को सीधे GitHub से `nextflow run nf-core/molkart -r 1.2.0` का उपयोग करके चलाते हैं।
    Nextflow आपके लिए स्वचालित रूप से अनुरोधित पाइपलाइन संस्करण को `$HOME/.nextflow/assets/nf-core/molkart` पर डाउनलोड करता है और वहां से चलाता है।
    हालांकि, इस प्रशिक्षण के लिए, हम पाइपलाइन को एक अलग स्थानीय डायरेक्टरी में क्लोन कर रहे हैं ताकि हम कोड का अधिक आसानी से निरीक्षण कर सकें।

### 2.1. कंटेनर आवश्यकताओं को समझना

पूर्ण पाइपलाइन चलाने से पहले, आइए जानें कि nf-core पाइपलाइनों के लिए कंटेनर क्यों आवश्यक हैं।

आइए molkart टेस्ट कॉन्फ़िगरेशन से टेस्ट डेटासेट और पैरामीटर का उपयोग करके पाइपलाइन चलाने का प्रयास करें:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "mesmer,cellpose,stardist" \
  --outdir results
```

आइए इन पैरामीटर्स को समझें:

- `--input`: नमूना मेटाडेटा युक्त samplesheet का पथ
- `--mindagap_tilesize`, `--mindagap_boxsize`, `--mindagap_loopnum`: ग्रिड पैटर्न फिलिंग के लिए पैरामीटर
- `--clahe_pyramid_tile`: कंट्रास्ट एन्हांसमेंट के लिए कर्नेल साइज़
- `--segmentation_method`: सेल सेगमेंटेशन के लिए कौन सा एल्गोरिदम उपयोग करना है
- `--outdir`: परिणाम कहां सेव करने हैं

!!! Warning "यह कमांड विफल होगी - यह जानबूझकर है!"

    हम जानबूझकर इसे कंटेनर के बिना चला रहे हैं ताकि यह प्रदर्शित कर सकें कि उनकी आवश्यकता क्यों है।

कुछ क्षणों के बाद, आपको इस तरह की एक एरर दिखेगी:

??? failure "कमांड आउटपुट"

    ```console
    ERROR ~ Error executing process > 'NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)'

    Caused by:
      Process `NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)` terminated with an error exit status (127)

    Command executed:

      duplicate_finder.py \
          spots.txt \
          90

    Command exit status:
      127

    Command error:
      .command.sh: line 3: duplicate_finder.py: command not found
    ```

**यहां क्या हो रहा है?**

एरर `command not found` (exit status 127) का मतलब है कि Nextflow ने `duplicate_finder.py` चलाने का प्रयास किया लेकिन इसे आपके सिस्टम पर नहीं मिला।
ऐसा इसलिए है क्योंकि:

1. पाइपलाइन को विशेष बायोइंफॉर्मेटिक्स सॉफ़्टवेयर इंस्टॉल होने की उम्मीद है
2. ये टूल (जैसे `duplicate_finder.py`, `apply_clahe.dask.py`, आदि) मानक Linux वितरणों का हिस्सा नहीं हैं
3. कंटेनर के बिना, Nextflow आपकी स्थानीय मशीन पर सीधे कमांड चलाने का प्रयास करता है

**इन टूल्स को कहां से आना चाहिए?**

आइए एक प्रोसेस मॉड्यूल का निरीक्षण करें यह देखने के लिए कि यह अपनी सॉफ़्टवेयर आवश्यकताओं को कैसे घोषित करता है।

CLAHE प्रीप्रोसेसिंग मॉड्यूल खोलें:

```bash
code molkart/modules/local/clahe/main.nf
```

लाइन 5 को देखें - आपको दिखेगा:

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

यह लाइन Nextflow को बताती है: "इस प्रोसेस को चलाने के लिए, Docker इमेज `ghcr.io/schapirolabor/molkart-local:v0.0.4` का उपयोग करें, जिसमें सभी आवश्यक सॉफ़्टवेयर हैं।"

प्रत्येक प्रोसेस घोषित करता है कि कौन सी कंटेनर इमेज अपने आवश्यक टूल प्रदान करती है।
हालांकि, Nextflow इन कंटेनरों का उपयोग तभी करता है जब आप इसे बताते हैं!

**समाधान: कॉन्फ़िगरेशन में Docker सक्षम करें**

### 2.2. Docker कॉन्फ़िगर करें और पाइपलाइन लॉन्च करें

Docker सक्षम करने के लिए, हमें `nextflow.config` फ़ाइल में `docker.enabled` को `false` से `true` में बदलना होगा।

कॉन्फ़िग फ़ाइल खोलें:

```bash
code nextflow.config
```

`docker.enabled = false` को `docker.enabled = true` में बदलें:

```groovy
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}
```

अब समान कमांड के साथ पाइपलाइन फिर से चलाएं:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose,mesmer,stardist" \
  --outdir results
```

इस बार, Nextflow:

1. कॉन्फ़िग से `docker.enabled = true` सेटिंग पढ़ेगा
2. आवश्यक Docker इमेज पुल करेगा (केवल पहली बार)
3. प्रत्येक प्रोसेस को उसकी निर्दिष्ट कंटेनर के अंदर चलाएगा
4. सफलतापूर्वक निष्पादित होगा क्योंकि सभी टूल कंटेनरों के अंदर उपलब्ध हैं

!!! Tip "कंटेनर क्यों महत्वपूर्ण हैं"

    अधिकांश nf-core पाइपलाइनों को कंटेनराइज़ेशन (Docker, Singularity, Podman, आदि) की **आवश्यकता** होती है क्योंकि:

    - वे विशेष बायोइंफॉर्मेटिक्स सॉफ़्टवेयर का उपयोग करते हैं जो मानक वातावरणों में उपलब्ध नहीं होते
    - कंटेनर पुनरुत्पादकता सुनिश्चित करते हैं - बिल्कुल समान सॉफ़्टवेयर संस्करण हर जगह चलते हैं
    - आपको दर्जनों टूल और उनकी निर्भरताओं को मैन्युअल रूप से इंस्टॉल करने की आवश्यकता नहीं है

    Nextflow में कंटेनरों के बारे में अधिक विवरण के लिए, Hello Nextflow प्रशिक्षण से [Hello Containers](../../hello_nextflow/05_hello_containers.md) देखें।

### 2.3. निष्पादन की निगरानी करें

जैसे-जैसे पाइपलाइन चलती है, आपको इस तरह का आउटपुट दिखेगा:

??? success "कमांड आउटपुट"

    ```console
    Nextflow 25.04.8 is available - Please consider updating your version to it

    N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/molkart` [soggy_kalam] DSL2 - revision: 5e54b29cb3 [dev]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/molkart 1.2.0dev
    ------------------------------------------------------
    Segmentation methods and options
      segmentation_method       : mesmer,cellpose,stardist

    Image preprocessing
      mindagap_boxsize          : 7
      mindagap_loopnum          : 100
      clahe_kernel              : 25
      mindagap_tilesize         : 90
      clahe_pyramid_tile        : 368

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv
      outdir                    : results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-10-18_22-22-21

    Core Nextflow options
      revision                  : dev
      runName                   : soggy_kalam
      containerEngine           : docker
      launchDir                 : /workspaces/training/nf4-science/imaging
      workDir                   : /workspaces/training/nf4-science/imaging/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/molkart
      userName                  : root
      profile                   : docker,test
      configFiles               :

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.10650748

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/molkart/blob/master/CITATIONS.md

    executor >  local (22)
    [c1/da5009] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2 ✔
    [73/8f5e8a] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2 ✔
    [ec/8f84d5] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1 ✔
    [a2/99349b] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1 ✔
    [95/c9b4b1] NFCORE_MOLKART:MOLKART:DEEPCELL_MESMER (mem_only)          [100%] 1 of 1 ✔
    [d4/1ebd1e] NFCORE_MOLKART:MOLKART:STARDIST (mem_only)                 [100%] 1 of 1 ✔
    [3e/3c0736] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1 ✔
    [a0/415c6a] NFCORE_MOLKART:MOLKART:MASKFILTER (mem_only)               [100%] 3 of 3 ✔
    [14/a830c9] NFCORE_MOLKART:MOLKART:SPOT2CELL (mem_only)                [100%] 3 of 3 ✔
    [b5/391836] NFCORE_MOLKART:MOLKART:CREATE_ANNDATA (mem_only)           [100%] 3 of 3 ✔
    [77/aed558] NFCORE_MOLKART:MOLKART:MOLKARTQC (mem_only)                [100%] 3 of 3 ✔
    [e6/b81475] NFCORE_MOLKART:MOLKART:MULTIQC                             [100%] 1 of 1 ✔
    -[nf-core/molkart] Pipeline completed successfully-
    Completed at: 19-Oct-2025 22:23:01
    Duration    : 2m 52s
    CPU hours   : 0.1
    Succeeded   : 22
    ```

ध्यान दें कि यह आउटपुट हमारे Hello World उदाहरण से अधिक विस्तृत है क्योंकि पाइपलाइन द्वारा अपनाई गई nf-core परंपराओं के कारण:

- पाइपलाइन अपना संस्करण और लोगो दिखाती है
- कॉन्फ़िगरेशन पैरामीटर प्रदर्शित किए जाते हैं
- कई प्रोसेस समानांतर में चलते हैं (कई प्रोसेस लाइनों द्वारा संकेत किया गया)
- प्रोसेस नामों में पूर्ण मॉड्यूल पथ शामिल है (जैसे, `NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP`)

### 2.4. प्रोसेस निष्पादन को समझना

executor लाइन `executor > local (22)` आपको बताती है:

- **executor**: किस कंप्यूट वातावरण का उपयोग किया जा रहा है (`local` = आपकी मशीन)
- **(22)**: लॉन्च किए गए कुल कार्यों की संख्या

प्रत्येक प्रोसेस लाइन दिखाती है:

- **Hash** (`[1a/2b3c4d]`): work डायरेक्टरी आइडेंटिफ़ायर (पहले की तरह)
- **Process name**: पूर्ण मॉड्यूल पथ और प्रोसेस नाम
- **Input identifier**: कोष्ठक में नमूना नाम
- **Progress**: पूर्णता का प्रतिशत और गिनती (जैसे, `1 of 1 ✔`)

### निष्कर्ष

आप जानते हैं कि टेस्ट डेटा के साथ nf-core पाइपलाइन कैसे लॉन्च करें और इसके निष्पादन आउटपुट की व्याख्या कैसे करें।

### आगे क्या है?

जानें कि परिणाम कहां मिलेंगे और उनकी व्याख्या कैसे करें।

---

## 3. आउटपुट ढूंढें और जांचें

जब पाइपलाइन सफलतापूर्वक पूर्ण होती है, तो आपको एक पूर्णता संदेश और निष्पादन सारांश दिखेगा।

### 3.1. परिणाम डायरेक्टरी का पता लगाएं

डिफ़ॉल्ट रूप से, nf-core पाइपलाइनें `outdir` पैरामीटर द्वारा निर्दिष्ट डायरेक्टरी में आउटपुट लिखती हैं, जिसे हमने `results/` पर सेट किया था।

सामग्री की सूची बनाएं:

```bash
tree results/
```

आपको कई subdirectories दिखनी चाहिए:

```console title="results/"
results/
├── anndata/
├── clahe/
├── mindagap/
├── molkartqc/
├── multiqc/
├── pipeline_info/
├── segmentation/
├── spot2cell/
└── stack/
```

प्रत्येक subdirectory में पाइपलाइन के एक विशिष्ट चरण से आउटपुट होता है:

- **mindagap/**: MindaGap प्रीप्रोसेसिंग स्टेप से ग्रिड-फिल्ड इमेज
- **clahe/**: CLAHE प्रीप्रोसेसिंग से कंट्रास्ट-एन्हांस्ड इमेज
- **stack/**: सेगमेंटेशन के लिए बनाए गए मल्टी-चैनल इमेज स्टैक
- **segmentation/**: विभिन्न एल्गोरिदम से सेगमेंटेशन परिणाम (cellpose/, mesmer/, stardist/, filtered_masks/)
- **spot2cell/**: सेल-बाय-ट्रांसक्रिप्ट काउंट टेबल
- **anndata/**: सेल-बाय-ट्रांसक्रिप्ट मैट्रिसेस और स्थानिक निर्देशांक युक्त AnnData ऑब्जेक्ट
- **molkartqc/**: स्पॉट असाइनमेंट के लिए गुणवत्ता नियंत्रण मेट्रिक्स
- **multiqc/**: व्यापक गुणवत्ता नियंत्रण रिपोर्ट
- **pipeline_info/**: निष्पादन रिपोर्ट और लॉग

### 3.2. MultiQC रिपोर्ट जांचें

MultiQC रिपोर्ट एक व्यापक HTML फ़ाइल है जो सभी पाइपलाइन स्टेप्स से गुणवत्ता मेट्रिक्स एकत्रित करती है।

फ़ाइल ब्राउज़र में रिपोर्ट खोलें और फिर VS Code में सीधे इसे रेंडर देखने के लिए "Show Preview" बटन पर क्लिक करें।

रिपोर्ट में शामिल है:

- सभी नमूनों के लिए सामान्य सांख्यिकी
- प्रीप्रोसेसिंग मेट्रिक्स
- सेगमेंटेशन गुणवत्ता मेट्रिक्स
- पाई गई सेल्स और स्पॉट्स की संख्या

!!! Tip

    MultiQC रिपोर्ट आमतौर पर सभी nf-core पाइपलाइनों में शामिल होती हैं।
    वे हमेशा पाइपलाइन निष्पादन और डेटा गुणवत्ता का उच्च-स्तरीय अवलोकन प्रदान करती हैं।

### 3.3. सेल-बाय-ट्रांसक्रिप्ट टेबल जांचें

सबसे महत्वपूर्ण वैज्ञानिक आउटपुट सेल-बाय-ट्रांसक्रिप्ट काउंट टेबल है।
यह आपको बताता है कि प्रत्येक सेल में प्रत्येक ट्रांसक्रिप्ट कितने पाए गए।

spot2cell डायरेक्टरी पर नेविगेट करें:

```bash
ls results/spot2cell/
```

आपको इस तरह की फ़ाइलें मिलेंगी:

- `cellxgene_mem_only_cellpose.csv`: Cellpose सेगमेंटेशन का उपयोग करके सेल-बाय-ट्रांसक्रिप्ट टेबल
- `cellxgene_mem_only_mesmer.csv`: Mesmer सेगमेंटेशन का उपयोग करके सेल-बाय-ट्रांसक्रिप्ट टेबल
- `cellxgene_mem_only_stardist.csv`: Stardist सेगमेंटेशन का उपयोग करके सेल-बाय-ट्रांसक्रिप्ट टेबल

हमने इस टेस्ट डेटासेट में केवल 1 नमूना चलाया, लेकिन एक वास्तविक प्रयोग में हमारे पास प्रत्येक नमूने के लिए ये टेबल होंगी।
ध्यान दें कि Nextflow कैसे समानांतर में कई सेगमेंटेशन विधियों को प्रोसेस कर सकता है, जिससे परिणामों की तुलना करना आसान हो जाता है।

### 3.4. निष्पादन रिपोर्ट देखें

Nextflow स्वचालित रूप से कई निष्पादन रिपोर्ट जेनरेट करता है।

pipeline_info डायरेक्टरी जांचें:

```bash
ls results/pipeline_info/
```

मुख्य फ़ाइलें:

- **execution_report.html**: टाइमलाइन और संसाधन उपयोग विज़ुअलाइज़ेशन
- **execution_timeline.html**: प्रोसेस निष्पादन का Gantt चार्ट
- **execution_trace.txt**: विस्तृत कार्य निष्पादन मेट्रिक्स
- **pipeline_dag.html**: वर्कफ़्लो संरचना दिखाने वाला निर्देशित एसाइक्लिक ग्राफ़

संसाधन उपयोग देखने के लिए निष्पादन रिपोर्ट खोलें:

```bash
code results/pipeline_info/execution_report.html
```

यह दिखाता है:

- प्रत्येक प्रोसेस को कितना समय लगा
- CPU और मेमोरी उपयोग
- कौन से कार्य कैश्ड थे बनाम निष्पादित

!!! Tip

    ये रिपोर्ट संसाधन आवंटन को अनुकूलित करने और प्रदर्शन समस्याओं को ट्रबलशूट करने के लिए अविश्वसनीय रूप से उपयोगी हैं।

### निष्कर्ष

आप जानते हैं कि पाइपलाइन आउटपुट कैसे खोजें, गुणवत्ता नियंत्रण रिपोर्ट की जांच करें, और निष्पादन मेट्रिक्स तक पहुंचें।

### आगे क्या है?

work डायरेक्टरी के बारे में जानें और Nextflow intermediate फ़ाइलों को कैसे प्रबंधित करता है।

---

## 4. work डायरेक्टरी का अन्वेषण करें

हमारे Hello World उदाहरण की तरह, सभी वास्तविक कार्य `work/` डायरेक्टरी में होता है।

### 4.1. work डायरेक्टरी संरचना को समझना

work डायरेक्टरी में निष्पादित प्रत्येक कार्य के लिए एक subdirectory होती है।
22 कार्यों वाली इस पाइपलाइन के लिए, 22 work subdirectories होंगी।

work डायरेक्टरी की सूची बनाएं:

```bash
ls -d work/*/*/ | head -5
```

यह पहली 5 कार्य डायरेक्टरियां दिखाता है।

### 4.2. एक कार्य डायरेक्टरी का निरीक्षण करें

console आउटपुट से एक सेगमेंटेशन प्रोसेस हैश (जैसे, `[3m/4n5o6p]`) चुनें और अंदर देखें:

```bash
ls -la work/3m/4n5o6p*/
```

आपको दिखेगा:

- **.command.\* फ़ाइलें**: Nextflow निष्पादन स्क्रिप्ट और लॉग (पहले की तरह)
- **Staged इनपुट फ़ाइलें**: वास्तविक इनपुट फ़ाइलों के symlinks
- **आउटपुट फ़ाइलें**: सेगमेंटेशन मास्क, intermediate परिणाम, आदि।

Hello World से मुख्य अंतर:

- वास्तविक पाइपलाइनें बड़ी इनपुट फ़ाइलों को stage करती हैं (इमेज, संदर्भ डेटा)
- आउटपुट फ़ाइलें काफी बड़ी हो सकती हैं (सेगमेंटेशन मास्क, प्रोसेस्ड इमेज)
- प्रति कार्य कई इनपुट और आउटपुट फ़ाइलें

!!! Tip

    यदि कोई प्रोसेस विफल होती है, तो आप उसकी work डायरेक्टरी पर नेविगेट कर सकते हैं, एरर संदेशों के लिए `.command.err` की जांच कर सकते हैं, और समस्या को debug करने के लिए `.command.sh` को मैन्युअली फिर से चला भी सकते हैं।

### 4.3. work डायरेक्टरी की सफाई

कई पाइपलाइन रन पर work डायरेक्टरी काफी बड़ी हो सकती है।
जैसा कि हमने भाग 1 में सीखा, आप पुराने रन से work डायरेक्टरियों को हटाने के लिए `nextflow clean` का उपयोग कर सकते हैं।

हालांकि, बड़ी intermediate फ़ाइलों वाली nf-core पाइपलाइनों के लिए, नियमित रूप से साफ करना विशेष रूप से महत्वपूर्ण है।

### निष्कर्ष

आप समझते हैं कि nf-core पाइपलाइनें अपनी work डायरेक्टरियों को कैसे व्यवस्थित करती हैं और debugging के लिए व्यक्तिगत कार्यों का निरीक्षण कैसे करें।

### आगे क्या है?

Nextflow cache के बारे में जानें और विफल पाइपलाइन रन को कैसे resume करें।

---

## 5. पाइपलाइन रन को resume करें

Nextflow की सबसे शक्तिशाली विशेषताओं में से एक विफलता के बिंदु से पाइपलाइन को resume करने की क्षमता है।

### 5.1. cache तंत्र

जब आप `-resume` के साथ पाइपलाइन चलाते हैं, तो Nextflow:

1. प्रत्येक कार्य के लिए cache की जांच करता है
2. यदि इनपुट, कोड, और पैरामीटर समान हैं, तो cached परिणाम का पुन: उपयोग करता है
3. केवल उन कार्यों को फिर से चलाता है जो बदल गए या विफल हुए

यह लंबे समय तक चलने वाली पाइपलाइनों के लिए आवश्यक है जहां निष्पादन के अंत में विफलताएं हो सकती हैं।

### 5.2. molkart के साथ resume का प्रयास करें

समान कमांड फिर से चलाएं, लेकिन `-resume` जोड़ें:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results \
  -resume
```

आपको इस तरह का आउटपुट दिखना चाहिए: <!-- TODO: full output -->

```console
executor >  local (0)
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2, cached: 2 ✔
[7f/8g9h0i] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1, cached: 1 ✔
[9h/0i1j2k] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1, cached: 1 ✔
[2k/3l4m5n] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1, cached: 1 ✔
...
```

ध्यान दें प्रत्येक प्रोसेस के लिए `cached: 2` या `cached: 1` - कुछ भी फिर से निष्पादित नहीं किया गया!

### 5.3. resume कब उपयोगी है

Resume विशेष रूप से मूल्यवान है जब:

- संसाधन सीमाओं के कारण पाइपलाइन विफल हो जाती है (मेमोरी समाप्त, समय सीमा पार)
- आपको upstream स्टेप्स को फिर से चलाए बिना downstream प्रोसेस को संशोधित करने की आवश्यकता है
- डेटा डाउनलोड के दौरान आपका नेटवर्क कनेक्शन टूट जाता है
- आप गणना को फिर से किए बिना अतिरिक्त आउटपुट जोड़ना चाहते हैं

!!! Warning

    Resume केवल तभी काम करता है जब आपने इनपुट डेटा, पाइपलाइन कोड, या पैरामीटर को नहीं बदला हो।
    यदि आप इनमें से किसी को भी बदलते हैं, तो Nextflow सही ढंग से प्रभावित कार्यों को फिर से चलाएगा।

### निष्कर्ष

आप जानते हैं कि सफल कार्यों को दोहराए बिना पाइपलाइनों को कुशलतापूर्वक फिर से चलाने के लिए `-resume` का उपयोग कैसे करें।

### आगे क्या है?

अब जब आप टेस्ट डेटा के साथ nf-core/molkart चला सकते हैं, तो आप अपने स्वयं के डेटासेट के लिए इसे कॉन्फ़िगर करना सीखने के लिए तैयार हैं।
