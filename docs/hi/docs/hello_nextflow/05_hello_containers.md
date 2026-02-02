# भाग 5: Hello Containers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube channel पर [पूरी playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) देखें।

:green_book: वीडियो transcript [यहाँ](./transcripts/05_hello_containers.md) उपलब्ध है।
///
-->

इस training course के Parts 1-4 में, तुमने सीखा कि कुछ text process करने, multiple inputs होने पर execution parallelize करने, और further processing के लिए results collect करने में सक्षम simple workflow assemble करने के लिए Nextflow के basic building blocks कैसे use करें।

हालाँकि, तुम अपने environment में available basic UNIX tools तक limited थे।
Real-world tasks को अक्सर विभिन्न tools और packages की need होती है जो default रूप से included नहीं हैं।
Typically, तुम्हें इन tools को install करना होगा, उनकी dependencies manage करनी होंगी, और किसी भी conflicts को resolve करना होगा।

यह सब बहुत tedious और annoying है, इसलिए हम तुम्हें दिखाएंगे कि इस problem को बहुत अधिक conveniently solve करने के लिए **containers** कैसे use करें।

एक **container** एक lightweight, standalone, executable unit of software है जो container **image** से बनाई जाती है जिसमें application run करने के लिए आवश्यक सब कुछ शामिल है जिसमें code, system libraries और settings शामिल हैं।
जैसा तुम imagine कर सकते हो, यह तुम्हारी pipelines को अधिक reproducible बनाने में बहुत helpful होगा।

Note करो कि हम यह [Docker](https://www.docker.com/get-started/) का उपयोग करके teach करेंगे, लेकिन ध्यान रखो कि Nextflow [कई अन्य container technologies](https://www.nextflow.io/docs/latest/container.html#) को भी support करता है।

??? info "इस section से कैसे शुरू करें"

    Course का यह section मानता है कि तुमने [Hello Nextflow](./index.md) course के Parts 1-4 complete कर लिए हैं और एक complete working pipeline है।

    यदि तुम इस point से course शुरू कर रहे हो, तो तुम्हें solutions से `modules` directory copy करनी होगी:

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. Warmup: `hello-containers.nf` चलाएं

हम starting point के रूप में workflow script `hello-containers.nf` use करेंगे।

यह sure करने के लिए कि सब कुछ काम कर रहा है, कोई भी changes करने से पहले script को एक बार run करो:

```bash
nextflow run hello-containers.nf
```

यदि यह तुम्हारे लिए काम किया, तो तुम containers use करना सीखने के लिए ready हो।

---

## 1. Container 'manually' use करें

हम क्या करना चाहते हैं अपने workflow में एक step add करना जो execution के लिए container use करेगा।

हालाँकि, हम पहले कुछ basic concepts और operations cover करेंगे ताकि Nextflow में use करने से पहले containers क्या हैं इसकी तुम्हारी understanding solidify हो।

### 1.1. Container image pull करें

Container use करने के लिए, तुम आमतौर पर container registry से container image download या _pull_ करते हो, और फिर container instance बनाने के लिए container image run करते हो।

General syntax इस प्रकार है:

```bash title="Syntax"
docker pull '<container>'
```

एक example के रूप में, चलो एक container image pull करते हैं जिसमें [cowpy](https://github.com/jeffbuttars/cowpy) है, `cowsay` नामक tool का python implementation जो arbitrary text inputs को fun तरीके से display करने के लिए ASCII art generate करता है।

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

### 1.2. One-off command के रूप में `cowpy` run करने के लिए container use करें

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

### 1.3. `cowpy` interactively run करने के लिए container use करें

#### 1.3.1. Container spin up करें

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

#### 1.3.2. Desired tool command(s) run करें

```bash
cowpy "Hello Containers" -c tux
```

#### 1.3.3. Container से exit करें

```bash
exit
```

#### 1.3.4. Container में data mount करें

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

#### 1.3.5. Mounted data use करें

```bash
cat /my_project/data/greetings.csv | cowpy -c turkey
```

### सीख

तुम जानते हो कि container pull करना और इसे one-off या interactively run करना। तुम यह भी जानते हो कि अपने data को अपने container के भीतर से accessible कैसे बनाना।

### आगे क्या?

सीखो कि Nextflow processes के execution के लिए containers कैसे use करें।

---

## 2. Nextflow में containers use करें

Nextflow में processes को containers के अंदर run करने के लिए built-in support है।

### 2.1. `cowpy` module लिखें

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// cowpy के साथ ASCII art generate करें
process cowpy {

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
    """

}
```

### 2.2. Workflow में cowpy add करें

#### 2.2.1. `hello-containers.nf` में `cowpy` process import करें

```groovy title="hello-containers.nf" linenums="3" hl_lines="5"
// Modules को include करें
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'
```

#### 2.2.2. Workflow में `cowpy` process का call add करें

```groovy title="hello-containers.nf" linenums="31" hl_lines="2"
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // cowpy के साथ अभिवादनों का ASCII art generate करें
        cowpy(collectGreetings.out.outfile, params.character)
```

### 2.3. `cowpy` process run करने के लिए container use करें

#### 2.3.1. `cowpy` के लिए container specify करें

```groovy title="modules/cowpy.nf" linenums="4" hl_lines="3"
process cowpy {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
```

#### 2.3.2. `nextflow.config` file के माध्यम से Docker use enable करें

```console title="nextflow.config" linenums="1" hl_lines="1"
docker.enabled = true
```

#### 2.3.3. Docker enabled के साथ workflow run करें

```bash
nextflow run hello-containers.nf -resume
```

### सीख

तुम जानते हो कि Nextflow में processes run करने के लिए containers कैसे use करें।

### आगे क्या?

Break लो!

जब तुम ready हो, तो [**Part 6: Hello Config**](./06_hello_config.md) पर move करो यह सीखने के लिए कि अपने infrastructure fit करने के लिए अपने pipeline के execution को configure कैसे करें साथ ही inputs और parameters के configuration को manage करें।

यह बिल्कुल last part है, और फिर तुम इस course से done हो जाओगे!

---

## Quiz

<quiz>
Container क्या है?
- [ ] एक type की virtual machine
- [ ] एक file compression format
- [x] एक lightweight, standalone executable unit जिसमें application run करने के लिए आवश्यक सब कुछ शामिल है
- [ ] एक network protocol
</quiz>

<quiz>
Container image और container instance में क्या difference है?
- [ ] वे same thing हैं
- [x] एक image एक template है; एक instance उस image से बनाया गया running container है
- [ ] एक instance एक template है; एक image एक running container है
- [ ] Images Docker के लिए हैं; instances Singularity के लिए हैं
</quiz>

<quiz>
`docker run` command में `-v` flag क्या करता है?
- [ ] Verbose output enable करता है
- [ ] Container validate करता है
- [x] Host system से container में volume mount करता है
- [ ] Container का version specify करता है

और जानें: [1.3.4. Mount data into the container](#134-mount-data-into-the-container)
</quiz>

<quiz>
Containers use करते समय volumes mount करने की need क्यों होती है?
- [ ] Container performance improve करने के लिए
- [ ] Disk space save करने के लिए
- [x] क्योंकि containers default रूप से host filesystem से isolated हैं
- [ ] Networking enable करने के लिए

और जानें: [1.3.4. Mount data into the container](#134-mount-data-into-the-container)
</quiz>

<quiz>
Nextflow process के लिए container कैसे specify करते हो?
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

और जानें: [2.3.1. Specify a container for cowpy](#231-specify-a-container-for-cowpy)
</quiz>

<quiz>
कौन सी `nextflow.config` setting तुम्हारे workflow के लिए Docker enable करती है?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

और जानें: [2.3.2. Enable use of Docker via the `nextflow.config` file](#232-enable-use-of-docker-via-the-nextflowconfig-file)
</quiz>

<quiz>
Container में process run करते समय Nextflow automatically क्या handle करता है? (सभी लागू select करें)
- [x] यदि needed हो तो container image pull करना
- [x] Work directory mount करना
- [x] Container के अंदर process script run करना
- [x] Execution के बाद container instance clean up करना

और जानें: [2.3.4. Inspect how Nextflow launched the containerized task](#234-inspect-how-nextflow-launched-the-containerized-task)
</quiz>
