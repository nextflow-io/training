# भाग 2: nf-core/molkart चलाना

भाग 1 में, हमने Nextflow execution की मूल बातें समझने के लिए एक सरल Hello World workflow चलाया था।
अब हम एक वास्तविक bioimaging pipeline चलाने जा रहे हैं: **nf-core/molkart**।

यह pipeline Resolve Bioscience से Molecular Cartography spatial transcriptomics डेटा को प्रोसेस करती है।
हालांकि, यहाँ तुम जो Nextflow patterns सीखोगे वे किसी भी nf-core pipeline या production workflow पर लागू होते हैं।

## 1. nf-core pipelines को समझना

Pipeline चलाने से पहले, आओ समझें कि nf-core क्या है और workflows चलाने के लिए यह क्यों महत्वपूर्ण है।

### 1.1. nf-core क्या है?

[nf-core](https://nf-co.re/) उच्च-गुणवत्ता वाली Nextflow pipelines का एक community-driven संग्रह है।
सभी nf-core pipelines एक ही संरचना और conventions का पालन करती हैं, जिसका मतलब है कि एक बार जब तुम एक को चलाना सीख जाते हो, तो तुम उनमें से किसी को भी चला सकते हो।

nf-core pipelines की मुख्य विशेषताएं:

- **मानकीकृत संरचना**: सभी pipelines में consistent parameter नाम और उपयोग patterns होते हैं
- **Built-in test डेटा**: हर pipeline में त्वरित validation के लिए test profiles शामिल हैं
- **व्यापक documentation**: विस्तृत उपयोग निर्देश और parameter विवरण
- **Quality control**: MultiQC का उपयोग करके automated QC रिपोर्ट
- **Container समर्थन**: reproducibility के लिए पहले से बने containers

!!! tip "nf-core के बारे में और जानना चाहते हो?"

    nf-core pipeline development के गहन परिचय के लिए, [Hello nf-core](../../hello_nf-core/index.md) प्रशिक्षण कोर्स देखो।
    यह शुरू से nf-core pipelines बनाने और customize करने को कवर करता है।

### 1.2. molkart pipeline

![nf-core/molkart pipeline](img/molkart.png)

[nf-core/molkart](https://nf-co.re/molkart) pipeline कई चरणों के माध्यम से spatial transcriptomics imaging डेटा को प्रोसेस करती है:

1. **Image preprocessing**: Grid pattern filling और optional contrast enhancement
2. **Cell segmentation**: कई algorithm विकल्प (Cellpose, Mesmer, ilastik, Stardist)
3. **Spot assignment**: Segmented cells को transcript spots assign करना
4. **Quality control**: व्यापक QC रिपोर्ट generate करना

मुख्य outputs हैं:

- Cell-by-transcript count tables
- Segmentation masks
- MultiQC quality control रिपोर्ट

---

## 2. Test डेटा के साथ molkart चलाना

शुरू करने से पहले, आओ molkart repository को locally clone करें ताकि हम इसके code को inspect कर सकें:

```bash
cd /workspaces/training/nf4-science/imaging
git clone --branch 1.2.0 --depth 1 https://github.com/nf-core/molkart
```

यह एक `molkart/` डायरेक्टरी बनाता है जिसमें पूरा pipeline source code होता है।

!!! note "हम locally क्यों clone कर रहे हैं?"

    आमतौर पर, तुम `nextflow run nf-core/molkart -r 1.2.0` का उपयोग करके सीधे GitHub से nf-core pipelines चलाओगे।
    Nextflow स्वचालित रूप से तुम्हारे लिए requested pipeline version को `$HOME/.nextflow/assets/nf-core/molkart` में download करता है और वहाँ से इसे चलाता है।
    हालांकि, इस प्रशिक्षण के लिए, हम pipeline को एक अलग local डायरेक्टरी में clone कर रहे हैं ताकि हम code को आसानी से inspect कर सकें।

### 2.1. Container requirements को समझना

पूरी pipeline चलाने से पहले, आओ सीखें कि nf-core pipelines के लिए containers क्यों आवश्यक हैं।

आओ molkart test configuration से test dataset और parameters का उपयोग करके pipeline चलाने की कोशिश करें:

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

आओ इन parameters को समझें:

- `--input`: Sample metadata युक्त samplesheet का path
- `--mindagap_tilesize`, `--mindagap_boxsize`, `--mindagap_loopnum`: Grid pattern filling के लिए parameters
- `--clahe_pyramid_tile`: Contrast enhancement के लिए kernel size
- `--segmentation_method`: Cell segmentation के लिए कौन सा algorithm(s) उपयोग करना है
- `--outdir`: परिणाम कहाँ save करने हैं

!!! Warning "यह कमांड fail होगी - यह जानबूझकर है!"

    हम जानबूझकर इसे containers के बिना चला रहे हैं ताकि यह प्रदर्शित कर सकें कि वे क्यों आवश्यक हैं।

कुछ क्षणों के बाद, तुम्हें इस तरह की एक error दिखाई देगी:

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

**यहाँ क्या हो रहा है?**

Error `command not found` (exit status 127) का मतलब है कि Nextflow ने `duplicate_finder.py` चलाने की कोशिश की लेकिन इसे तुम्हारे system पर नहीं मिला।
ऐसा इसलिए है क्योंकि:

1. Pipeline को विशेष bioinformatics software installed होने की उम्मीद है
2. ये tools (जैसे `duplicate_finder.py`, `apply_clahe.dask.py`, आदि) standard Linux distributions का हिस्सा नहीं हैं
3. Containers के बिना, Nextflow सीधे तुम्हारी local machine पर commands चलाने की कोशिश करता है

**ये tools कहाँ से आने चाहिए?**

आओ एक process module को inspect करें यह देखने के लिए कि यह अपनी software requirements को कैसे declare करता है।

CLAHE preprocessing module खोलो:

```bash
code molkart/modules/local/clahe/main.nf
```

Line 5 को देखो - तुम्हें दिखाई देगा:

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

यह line Nextflow को बताती है: "इस process को चलाने के लिए, Docker image `ghcr.io/schapirolabor/molkart-local:v0.0.4` का उपयोग करो, जिसमें सभी आवश्यक software हैं।"

प्रत्येक process declare करता है कि कौन सा container image इसके आवश्यक tools प्रदान करता है।
हालांकि, Nextflow इन containers का उपयोग तभी करता है जब तुम इसे बताते हो!

**समाधान: Configuration में Docker enable करना**

### 2.2. Docker configure करना और pipeline launch करना

Docker enable करने के लिए, हमें `nextflow.config` फ़ाइल में `docker.enabled` को `false` से `true` में बदलना होगा।

Config फ़ाइल खोलो:

```bash
code nextflow.config
```

`docker.enabled = false` को `docker.enabled = true` में बदलो:

```groovy
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}
```

अब उसी कमांड के साथ pipeline फिर से चलाओ:

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

1. Config से `docker.enabled = true` setting पढ़ेगा
2. आवश्यक Docker images pull करेगा (पहली बार केवल)
3. प्रत्येक process को उसके specified container के अंदर चलाएगा
4. सफलतापूर्वक execute होगा क्योंकि सभी tools containers के अंदर उपलब्ध हैं

!!! Tip "Containers क्यों महत्वपूर्ण हैं"

    अधिकांश nf-core pipelines को containerization (Docker, Singularity, Podman, आदि) की **आवश्यकता** होती है क्योंकि:

    - वे विशेष bioinformatics software का उपयोग करती हैं जो standard environments में उपलब्ध नहीं है
    - Containers reproducibility सुनिश्चित करते हैं - बिल्कुल वही software versions हर जगह चलते हैं
    - तुम्हें manually दर्जनों tools और उनकी dependencies install करने की आवश्यकता नहीं है

    Nextflow में containers के बारे में अधिक विवरण के लिए, Hello Nextflow प्रशिक्षण से [Hello Containers](../../hello_nextflow/05_hello_containers.md) देखो।

### 2.3. Execution monitor करना

जैसे ही pipeline चलती है, तुम्हें इस तरह का आउटपुट दिखाई देगा:

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

ध्यान दो कि यह आउटपुट हमारे Hello World उदाहरण से अधिक विस्तृत है क्योंकि pipeline nf-core conventions का पालन करती है:

- Pipeline अपना version और logo दिखाती है
- Configuration parameters प्रदर्शित किए जाते हैं
- कई processes parallel में चलती हैं (कई process lines द्वारा संकेतित)
- Process नाम पूरा module path शामिल करते हैं (जैसे, `NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP`)

### 2.4. Process execution को समझना

Executor line `executor > local (22)` तुम्हें बताती है:

- **executor**: कौन सा compute environment उपयोग किया जा रहा है (`local` = तुम्हारी machine)
- **(22)**: Launch किए गए tasks की कुल संख्या

प्रत्येक process line दिखाती है:

- **Hash** (`[1a/2b3c4d]`): Work डायरेक्टरी identifier (पहले की तरह)
- **Process नाम**: पूरा module path और process नाम
- **Input identifier**: Parentheses में sample नाम
- **Progress**: पूर्ण प्रतिशत और count (जैसे, `1 of 1 ✔`)

### सारांश

तुम जानते हो कि test डेटा के साथ nf-core pipeline कैसे launch करें और इसके execution आउटपुट को कैसे interpret करें।

### आगे क्या है?

सीखो कि परिणाम कहाँ मिलेंगे और उन्हें कैसे interpret करें।

---

## 3. Outputs खोजना और examine करना

जब pipeline सफलतापूर्वक पूरी होती है, तो तुम्हें एक completion संदेश और execution सारांश दिखाई देगा।

### 3.1. Results डायरेक्टरी locate करना

डिफ़ॉल्ट रूप से, nf-core pipelines outputs को `outdir` parameter द्वारा निर्दिष्ट डायरेक्टरी में लिखती हैं, जिसे हमने `results/` पर सेट किया था।

Contents को list करो:

```bash
tree results/
```

तुम्हें कई subdirectories दिखाई देनी चाहिए:

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

प्रत्येक subdirectory में pipeline के एक विशिष्ट चरण से outputs होते हैं:

- **mindagap/**: MindaGap preprocessing step से grid-filled images
- **clahe/**: CLAHE preprocessing से contrast-enhanced images
- **stack/**: Segmentation के लिए बनाए गए multi-channel image stacks
- **segmentation/**: विभिन्न algorithms से segmentation परिणाम (cellpose/, mesmer/, stardist/, filtered_masks/)
- **spot2cell/**: Cell-by-transcript count tables
- **anndata/**: Cell-by-transcript matrices और spatial coordinates युक्त AnnData objects
- **molkartqc/**: Spot assignment के लिए quality control metrics
- **multiqc/**: व्यापक quality control रिपोर्ट
- **pipeline_info/**: Execution रिपोर्ट और logs

### 3.2. MultiQC रिपोर्ट examine करना

MultiQC रिपोर्ट एक व्यापक HTML फ़ाइल है जो सभी pipeline steps से quality metrics को aggregate करती है।

File browser में रिपोर्ट खोलो और फिर इसे सीधे VS Code में rendered देखने के लिए "Show Preview" बटन पर क्लिक करो।

रिपोर्ट में शामिल हैं:

- सभी samples के लिए सामान्य statistics
- Preprocessing metrics
- Segmentation quality metrics
- Detected cells और spots की संख्या

!!! Tip

    MultiQC रिपोर्ट आमतौर पर सभी nf-core pipelines में शामिल होती हैं।
    वे हमेशा pipeline execution और डेटा quality का high-level overview प्रदान करती हैं।

### 3.3. Cell-by-transcript tables examine करना

सबसे महत्वपूर्ण scientific आउटपुट cell-by-transcript count table है।
यह तुम्हें बताता है कि प्रत्येक cell में प्रत्येक transcript कितने detect किए गए।

spot2cell डायरेक्टरी में navigate करो:

```bash
ls results/spot2cell/
```

तुम्हें इस तरह की फ़ाइलें मिलेंगी:

- `cellxgene_mem_only_cellpose.csv`: Cellpose segmentation का उपयोग करके cell-by-transcript table
- `cellxgene_mem_only_mesmer.csv`: Mesmer segmentation का उपयोग करके cell-by-transcript table
- `cellxgene_mem_only_stardist.csv`: Stardist segmentation का उपयोग करके cell-by-transcript table

हमने इस test dataset में केवल 1 sample चलाया, लेकिन एक वास्तविक experiment में हमारे पास प्रत्येक sample के लिए ये tables होंगे।
ध्यान दो कि Nextflow कैसे कई segmentation methods को parallel में प्रोसेस करने में सक्षम है, जिससे परिणामों की तुलना करना आसान हो जाता है।

### 3.4. Execution रिपोर्ट देखना

Nextflow स्वचालित रूप से कई execution रिपोर्ट generate करता है।

pipeline_info डायरेक्टरी check करो:

```bash
ls results/pipeline_info/
```

मुख्य फ़ाइलें:

- **execution_report.html**: Timeline और resource usage visualization
- **execution_timeline.html**: Process execution का Gantt chart
- **execution_trace.txt**: विस्तृत task execution metrics
- **pipeline_dag.html**: Workflow संरचना दिखाने वाला directed acyclic graph

Resource usage देखने के लिए execution रिपोर्ट खोलो:

```bash
code results/pipeline_info/execution_report.html
```

यह दिखाता है:

- प्रत्येक process को कितना समय लगा
- CPU और memory उपयोग
- कौन से tasks cached बनाम executed थे

!!! Tip

    ये रिपोर्ट resource allocation को optimize करने और performance issues को troubleshoot करने के लिए अविश्वसनीय रूप से उपयोगी हैं।

### सारांश

तुम जानते हो कि pipeline outputs को कैसे locate करें, quality control रिपोर्ट को कैसे examine करें, और execution metrics को कैसे access करें।

### आगे क्या है?

Work डायरेक्टरी के बारे में सीखो और Nextflow intermediate फ़ाइलों को कैसे manage करता है।

---

## 4. Work डायरेक्टरी explore करना

हमारे Hello World उदाहरण की तरह, सभी वास्तविक कार्य `work/` डायरेक्टरी में होता है।

### 4.1. Work डायरेक्टरी संरचना को समझना

Work डायरेक्टरी में प्रत्येक task के लिए एक subdirectory होती है जो execute किया गया था।
12 tasks वाली इस pipeline के लिए, 12 work subdirectories होंगी।

Work डायरेक्टरी को list करो:

```bash
ls -d work/*/*/ | head -5
```

यह पहली 5 task डायरेक्टरी दिखाता है।

### 4.2. एक task डायरेक्टरी inspect करना

Console आउटपुट से segmentation process hashes में से एक चुनो (जैसे, `[3m/4n5o6p]`) और अंदर देखो:

```bash
ls -la work/3m/4n5o6p*/
```

तुम्हें दिखाई देगा:

- **.command.\* फ़ाइलें**: Nextflow execution scripts और logs (पहले की तरह)
- **Staged input फ़ाइलें**: वास्तविक input फ़ाइलों के symlinks
- **Output फ़ाइलें**: Segmentation masks, intermediate परिणाम, आदि।

Hello World से मुख्य अंतर:

- वास्तविक pipelines बड़ी input फ़ाइलें stage करती हैं (images, reference डेटा)
- Output फ़ाइलें काफी बड़ी हो सकती हैं (segmentation masks, processed images)
- प्रति task कई input और output फ़ाइलें

!!! Tip

    यदि कोई process fail होती है, तो तुम इसकी work डायरेक्टरी में navigate कर सकते हो, error संदेशों के लिए `.command.err` examine कर सकते हो, और issue को debug करने के लिए `.command.sh` को manually भी फिर से चला सकते हो।

### 4.3. Work डायरेक्टरी cleanup

कई pipeline runs के बाद work डायरेक्टरी काफी बड़ी हो सकती है।
जैसा कि हमने भाग 1 में सीखा, तुम पुराने runs से work डायरेक्टरी हटाने के लिए `nextflow clean` का उपयोग कर सकते हो।

हालांकि, बड़ी intermediate फ़ाइलों वाली nf-core pipelines के लिए, नियमित रूप से cleanup करना विशेष रूप से महत्वपूर्ण है।

### सारांश

तुम समझते हो कि nf-core pipelines अपनी work डायरेक्टरी को कैसे organize करती हैं और debugging के लिए individual tasks को कैसे inspect करें।

### आगे क्या है?

Nextflow cache के बारे में सीखो और failed pipeline runs को कैसे resume करें।

---

## 5. एक pipeline run resume करना

Nextflow की सबसे शक्तिशाली विशेषताओं में से एक failure के बिंदु से pipeline को resume करने की क्षमता है।

### 5.1. Cache mechanism

जब तुम `-resume` के साथ pipeline चलाते हो, तो Nextflow:

1. प्रत्येक task के लिए cache check करता है
2. यदि inputs, code, और parameters समान हैं, तो cached परिणाम का पुन: उपयोग करता है
3. केवल उन tasks को फिर से चलाता है जो बदले या fail हुए

यह लंबे समय तक चलने वाली pipelines के लिए आवश्यक है जहाँ failures execution में देर से हो सकती हैं।

### 5.2. molkart के साथ resume try करना

उसी कमांड को फिर से चलाओ, लेकिन `-resume` जोड़ो:

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

तुम्हें इस तरह का आउटपुट दिखाई देना चाहिए: <!-- TODO: full output -->

```console
executor >  local (0)
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2, cached: 2 ✔
[7f/8g9h0i] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1, cached: 1 ✔
[9h/0i1j2k] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1, cached: 1 ✔
[2k/3l4m5n] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1, cached: 1 ✔
...
```

प्रत्येक process के लिए `cached: 2` या `cached: 1` ध्यान दो - कुछ भी फिर से execute नहीं किया गया!

### 5.3. Resume कब उपयोगी है

Resume विशेष रूप से तब valuable है जब:

- एक pipeline resource limits के कारण fail होती है (out of memory, time limit exceeded)
- तुम्हें upstream steps को फिर से चलाए बिना downstream processes को modify करने की आवश्यकता है
- डेटा download के दौरान तुम्हारा network connection drop हो जाता है
- तुम computation को फिर से किए बिना अतिरिक्त outputs जोड़ना चाहते हो

!!! Warning

    Resume केवल तभी काम करता है जब तुमने input डेटा, pipeline code, या parameters को नहीं बदला है।
    यदि तुम इनमें से कोई भी बदलते हो, तो Nextflow सही तरीके से प्रभावित tasks को फिर से चलाएगा।

### सारांश

तुम जानते हो कि सफल tasks को दोहराए बिना pipelines को efficiently फिर से चलाने के लिए `-resume` का उपयोग कैसे करें।

### आगे क्या है?

अब जब तुम test डेटा के साथ nf-core/molkart चला सकते हो, तो तुम अपने स्वयं के datasets के लिए इसे configure करना सीखने के लिए तैयार हो।
