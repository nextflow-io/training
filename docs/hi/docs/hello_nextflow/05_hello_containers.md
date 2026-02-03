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

Note करो कि हम यह [Docker](https://www.docker.com/get-started/) का उपयोग करके teach करेंगे, लेकिन ध्यान रखो कि Nextflow [कई अन्य container technologies](https://nextflow.io/docs/latest/container.html) को भी support करता है।

??? info "इस section से कैसे शुरू करें"

    Course का यह section मानता है कि तुमने [Hello Nextflow](./index.md) course के Parts 1-4 complete कर लिए हैं और एक complete working pipeline है।

    यदि तुम इस point से course शुरू कर रहे हो, तो तुम्हें solutions से `modules` directory copy करनी होगी:

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. वार्मअप: `hello-containers.nf` चलाएं

हम starting point के रूप में workflow script `hello-containers.nf` use करेंगे।
यह Part 4 of this training course को पूरा करके बनाई गई script के बराबर है, सिवाय इसके कि हमने output destinations बदल दिए हैं:

```groovy title="hello-containers.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_containers'
        mode 'copy'
    }
    uppercased {
        path 'hello_containers'
        mode 'copy'
    }
    collected {
        path 'hello_containers'
        mode 'copy'
    }
    batch_report {
        path 'hello_containers'
        mode 'copy'
    }
}
```

यह sure करने के लिए कि सब कुछ काम कर रहा है, कोई भी changes करने से पहले script को एक बार run करो:

```bash
nextflow run hello-containers.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [5a/ec1fa1] sayHello (2) [100%] 3 of 3 ✔
    [30/32b5b8] convertToUpper (3) [100%] 3 of 3 ✔
    [d3/be01bc] collectGreetings [100%] 1 of 1 ✔

    ```

पहले की तरह, तुम्हें output files `output` block (`results/hello_containers/`) में specified directory में मिलेंगी।

??? abstract "डायरेक्टरी कंटेंट्स"

    ```console
    results/hello_containers/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

यदि यह तुम्हारे लिए काम किया, तो तुम containers use करना सीखने के लिए ready हो।

---

## 1. Container 'manually' use करें

हम क्या करना चाहते हैं अपने workflow में एक step add करना जो execution के लिए container use करेगा।

हालाँकि, हम पहले कुछ basic concepts और operations cover करेंगे ताकि Nextflow में use करने से पहले containers क्या हैं इसकी तुम्हारी understanding solidify हो।

### 1.1. Container image pull करें

Container use करने के लिए, तुम आमतौर पर container registry से container image download या _pull_ करते हो, और फिर container instance बनाने के लिए container image run करते हो।

General syntax इस प्रकार है:

```bash title="सिंटैक्स"
docker pull '<container>'
```

`docker pull` part container system को container repository से container image pull करने का instruction है।

`'<container>'` part container image का URI address है।

एक example के रूप में, चलो एक container image pull करते हैं जिसमें [cowpy](https://github.com/jeffbuttars/cowpy) है, `cowsay` नामक tool का python implementation जो arbitrary text inputs को fun तरीके से display करने के लिए ASCII art generate करता है।

```txt title="उदाहरण"
 ________________________
< Are we having fun yet? >
 ------------------------
    \                                  ___-------___
     \                             _-~~             ~~-_
      \                         _-~                    /~-_
             /^\__/^\         /~  \                   /    \
           /|  O|| O|        /      \_______________/        \
          | |___||__|      /       /                \          \
          |          \    /      /                    \          \
          |   (_______) /______/                        \_________ \
          |         / /         \                      /            \
           \         \^\\         \                  /               \     /
             \         ||           \______________/      _-_       //\__//
               \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
                 ~-----||====/~     |==================|       |/~~~~~
                  (_(__/  ./     /                    \_\      \.
                         (_(___/                         \_____)_)
```

विभिन्न repositories हैं जहाँ तुम published containers पा सकते हो।
हमने [Seqera Containers](https://seqera.io/containers/) service use की `cowpy` Conda package से यह Docker container image generate करने के लिए: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`।

पूरी pull command run करो:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "कमांड आउटपुट"

    ```console
    1.1.5--3db457ae1977a273: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    c23bdb422167: Pull complete
    e1686ff32a11: Pull complete
    Digest: sha256:1ebc0043e8cafa61203bf42d29fd05bd14e7b4298e5e8cf986504c15f5aa4160
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

यदि तुमने image पहले कभी download नहीं की है, तो complete होने में एक minute लग सकता है।
एक बार हो जाने पर, तुम्हारे पास container image की local copy होगी।

### 1.2. One-off command के रूप में `cowpy` run करने के लिए container use करें

एक बहुत common तरीका जो लोग containers use करते हैं वह है उन्हें directly run करना, _यानी_ non-interactively।
यह one-off commands run करने के लिए बहुत अच्छा है।

General syntax इस प्रकार है:

```bash title="सिंटैक्स"
docker run --rm '<container>' [tool command]
```

`docker run --rm '<container>'` part container system को container image से container instance spin up करने और उसमें command execute करने का instruction है।
`--rm` flag system को command complete होने के बाद container instance shut down करने के लिए बताता है।

`[tool command]` syntax tool पर depend करता है जो तुम use कर रहे हो और container कैसे set up है।
चलो बस `cowpy` से शुरू करते हैं।

पूरी तरह से assembled, container execution command इस तरह दिखती है; आगे बढ़ो और इसे run करो।

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

??? success "कमांड आउटपुट"

    ```console
    ______________________________________________________
    < Cowacter, eyes:default, tongue:False, thoughts:False >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

System ने container spin up किया, `cowpy` command अपने parameters के साथ run की, output console में भेजा और अंत में, container instance shut down कर दिया।

### 1.3. `cowpy` interactively run करने के लिए container use करें

तुम एक container को interactively भी run कर सकते हो, जो तुम्हें container के अंदर shell prompt देता है और तुम्हें command के साथ play करने की अनुमति देता है।

#### 1.3.1. Container spin up करें

Interactively run करने के लिए, हम बस `docker run` command में `-it` add करते हैं।
Optionally, हम shell specify कर सकते हैं जो हम container के अंदर use करना चाहते हैं _जैसे_ `/bin/bash` को command के अंत में जोड़कर।

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

ध्यान दो कि तुम्हारा prompt कुछ `(base) root@b645838b3314:/tmp#` जैसा बदल जाता है, जो indicate करता है कि तुम अब container के अंदर हो।

तुम यह verify कर सकते हो `ls /` run करके filesystem के root से directory contents list करने के लिए:

```bash
ls /
```

??? abstract "कमांड आउटपुट"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

हम यहाँ `tree` के बजाय `ls` use करते हैं क्योंकि `tree` utility इस container में available नहीं है।
तुम देख सकते हो कि container के अंदर की filesystem तुम्हारे host system की filesystem से अलग है।

जो हमने अभी किया उसकी एक limitation यह है कि container default रूप से host system से पूरी तरह से isolated है।
इसका मतलब है कि container host system पर किसी भी files को access नहीं कर सकता जब तक तुम explicitly इसे ऐसा करने की अनुमति नहीं दो।

हम तुम्हें एक minute में दिखाएंगे कि यह कैसे करना है।

#### 1.3.2. Desired tool command(s) run करें

अब जब तुम container के अंदर हो, तुम सीधे `cowpy` command run कर सकते हो और इसे कुछ parameters दे सकते हो।
उदाहरण के लिए, tool documentation कहता है कि हम character ('cowacter') को `-c` के साथ बदल सकते हैं।

```bash
cowpy "Hello Containers" -c tux
```

??? success "कमांड आउटपुट"

    ```console
    __________________
    < Hello Containers >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

अब output Linux penguin, Tux, दिखाता है, default cow के बजाय, क्योंकि हमने `-c tux` parameter specify किया।

क्योंकि तुम container के अंदर हो, तुम जितनी बार चाहो `cowpy` command run कर सकते हो, input parameters vary करते हुए, Docker commands की चिंता किए बिना।

!!! Tip

    एक अलग character pick करने के लिए '-c' flag use करो, जिसमें शामिल हैं:
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

यह neat है। और भी neater होगा अगर हम अपनी `greetings.csv` को इसमें input के रूप में feed कर सकें।
लेकिन चूंकि हमारे पास filesystem तक access नहीं है, हम नहीं कर सकते।

चलो इसे fix करें।

#### 1.3.3. Container से exit करें

Container से exit करने के लिए, तुम prompt पर `exit` type कर सकते हो या ++ctrl+d++ keyboard shortcut use कर सकते हो।

```bash
exit
```

तुम्हारा prompt अब वापस वह होना चाहिए जो container शुरू करने से पहले था।

#### 1.3.4. Container में data mount करें

जैसा earlier note किया गया, container default रूप से host system से isolated है।

Container को host filesystem access करने की अनुमति देने के लिए, तुम host system से container में एक **volume** **mount** कर सकते हो निम्नलिखित syntax use करके:

```bash title="सिंटैक्स"
-v <outside_path>:<inside_path>
```

हमारे case में `<outside_path>` current working directory होगी, इसलिए हम बस एक dot (`.`) use कर सकते हैं, और `<inside_path>` बस एक alias है जो हम बनाते हैं; चलो इसे `/my_project` कहते हैं (inside path absolute होना चाहिए)।

एक volume mount करने के लिए, हम paths replace करते हैं और docker run command में volume mounting argument add करते हैं इस प्रकार:

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

यह current working directory को एक volume के रूप में mount करता है जो container के अंदर `/my_project` के तहत accessible होगी।

तुम check कर सकते हो कि यह काम करता है `/my_project` की contents list करके:

```bash
ls /my_project
```

??? success "कमांड आउटपुट"

    ```console
    data               hello-config.nf      hello-modules.nf   hello-world.nf  nextflow.config  solutions         work
    hello-channels.nf  hello-containers.nf  hello-workflow.nf  modules         results          test-params.json
    ```

अब तुम container के अंदर से working directory की contents देख सकते हो, जिसमें `data/` के तहत `greetings.csv` file शामिल है।

इसने effectively container wall के through एक tunnel establish की जिसे तुम अपने filesystem के उस भाग को access करने के लिए use कर सकते हो।

#### 1.3.5. Mounted data use करें

अब जब हमने working directory को container में mount कर लिया है, हम `cowpy` command use कर सकते हैं `greetings.csv` file की contents display करने के लिए।

ऐसा करने के लिए, हम `cat /my_project/data/greetings.csv | ` use करेंगे CSV file की contents को `cowpy` command में pipe करने के लिए।

```bash
cat /my_project/data/greetings.csv | cowpy -c turkey
```

??? success "कमांड आउटपुट"

    ```console title="data/greetings.csv"
     ____________________
    / Hello,English,123  \
    | Bonjour,French,456 |
    \ Holà,Spanish,789   /
    --------------------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

यह desired ASCII art of a turkey produce करता है जो हमारे example greetings rattle कर रही है!
सिवाय यहाँ turkey full rows repeat कर रही है बजाय सिर्फ greetings के।
हम पहले से ही जानते हैं कि हमारा Nextflow workflow बेहतर job करेगा!

इस command के साथ play करने में feel free। जब तुम done हो, पहले की तरह container से exit करो:

```bash
exit
```

तुम खुद को अपने normal shell में वापस पाओगे।

### सीख

तुम जानते हो कि container pull करना और इसे one-off या interactively run करना। तुम यह भी जानते हो कि अपने data को अपने container के भीतर से accessible कैसे बनाना, जो तुम्हें अपने system पर कोई software install किए बिना real data पर किसी भी tool को try करने देता है।

### आगे क्या है?

सीखो कि Nextflow processes के execution के लिए containers कैसे use करें।

---

## 2. Nextflow में containers use करें

Nextflow में processes को containers के अंदर run करने के लिए built-in support है ताकि तुम ऐसे tools run कर सको जो तुम्हारे compute environment में installed नहीं हैं।
इसका मतलब है कि तुम अपनी processes run करने के लिए कोई भी container image use कर सकते हो, और Nextflow image pull करने, data mount करने, और इसके अंदर process run करने का ध्यान रखेगा।

इसे demonstrate करने के लिए, हम pipeline में एक `cowpy` step add करने जा रहे हैं जो हम develop कर रहे हैं, `collectGreetings` step के बाद।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

### 2.1. `cowpy` module लिखें

पहले, चलो `cowpy` process module बनाते हैं।

#### 2.1.1. नए module के लिए file stub बनाएं

Module के लिए एक empty file बनाओ जिसका नाम `cowpy.nf` है।

```bash
touch modules/cowpy.nf
```

यह हमें process code रखने के लिए एक जगह देता है।

#### 2.1.2. Module file में `cowpy` process code copy करें

हम अपनी `cowpy` process को पहले लिखी गई अन्य processes पर model कर सकते हैं।

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

Process एक `input_file` expect करती है जिसमें greetings हैं साथ ही एक `character` value।

Output एक नई text file होगी जिसमें `cowpy` tool द्वारा generated ASCII art होगी।

### 2.2. Workflow में cowpy add करें

अब हमें module import करना है और process call करना है।

#### 2.2.1. `hello-containers.nf` में `cowpy` process import करें

Workflow block के ऊपर import declaration insert करो और इसे appropriately fill करो।

=== "बाद में"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="5"
    // Modules को include करें
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

=== "पहले"

    ```groovy title="hello-containers.nf" linenums="3"
    // Modules को include करें
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    ```

अब `cowpy` module workflow में use करने के लिए available है।

#### 2.2.2. Workflow में `cowpy` process का call add करें

चलो `cowpy()` process को `collectGreetings()` process के output से connect करें, जो जैसा तुम याद कर सकते हो दो outputs produce करती है:

- `collectGreetings.out.outfile` output file contain करती है <--_जो हम चाहते हैं_
- `collectGreetings.out.report` report file contain करती है जिसमें greetings per batch की count है

Workflow block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-containers.nf" linenums="19" hl_lines="12-13"
        main:
        // CSV file से inputs के लिए channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक greeting emit करें
        sayHello(greeting_ch)
        // greeting को uppercase में convert करें
        convertToUpper(sayHello.out)
        // सभी greetings को एक file में collect करें
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // cowpy के साथ greetings का ASCII art generate करें
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "पहले"

    ```groovy title="hello-containers.nf" linenums="19"
        main:
        // CSV file से inputs के लिए channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक greeting emit करें
        sayHello(greeting_ch)
        // greeting को uppercase में convert करें
        convertToUpper(sayHello.out)
        // सभी greetings को एक file में collect करें
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

ध्यान दो कि हमने एक नया CLI parameter, `params.character`, declare किया ताकि specify कर सकें कि हम कौन सा character greetings कहता देखना चाहते हैं।

#### 2.2.3. `params` block में `character` parameter add करें

यह technically optional है लेकिन यह recommended practice है और यह character के लिए default value set करने का अवसर है।

=== "बाद में"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="7"
    /*
    * Pipeline parameters
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

=== "पहले"

    ```groovy title="hello-containers.nf" linenums="9"
    /*
    * Pipeline parameters
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

अब हम lazy हो सकते हैं और अपनी command lines में character parameter type करना skip कर सकते हैं।

#### 2.2.4. Workflow outputs update करें

हमें workflow outputs update करने की need है `cowpy` process की output publish करने के लिए।

##### 2.2.4.1. `publish:` section update करें

`Workflow block` में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-containers.nf" linenums="34" hl_lines="6"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    ```

=== "पहले"

    ```groovy title="hello-containers.nf" linenums="34"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

`cowpy` process केवल एक output produce करती है इसलिए हम इसे usual तरीके से refer कर सकते हैं `.out` append करके।

लेकिन अभी के लिए, चलो workflow-level outputs update करना finish करें।

##### 2.2.4.2. `output` block update करें

हमें final `cowpy_art` output को `output` block में add करने की need है। जब हम यह कर रहे हैं, चलो publishing destinations भी edit करें क्योंकि अब हमारी pipeline complete है और हम जानते हैं कि कौन से outputs हमें वाकई में चाहिए।

`Output` block में, निम्नलिखित code changes करो:

=== "बाद में"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15 18-21"
    output {
        first_output {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

=== "पहले"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15"
    output {
        first_output {
            path 'hello_containers'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers'
            mode 'copy'
        }
        collected {
            path 'hello_containers'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

अब published outputs थोड़ा अधिक organized होंगे।

#### 2.2.5. Workflow run करें

बस recap करने के लिए, यह वह है जिसका हम aim कर रहे हैं:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

क्या तुम सोचते हो यह काम करेगा?

चलो एक clean slate के लिए previous published outputs delete करें, और workflow को `-resume` flag के साथ run करें।

```bash
rm -r hello_containers/
nextflow run hello-containers.nf -resume
```

??? failure "कमांड आउटपुट (clarity के लिए edited)"

    ```console hl_lines="10 13 20-21 26-27"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [lonely_woese] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [9b/02e776] cowpy              [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (127)


    Command executed:

      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/hello-nextflow/work/9b/02e7761db848f82db3c3e59ff3a9b6

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ERROR ~ Cannot access first() element from an empty List

    -- Check '.nextflow.log' file for details
    ```

ओह नहीं, एक error है!
`error exit status (127)` द्वारा दिया गया error code का मतलब है कि हमने जो executable माँगा वह नहीं मिला।

यह make sense करता है, क्योंकि हम `cowpy` tool call कर रहे हैं लेकिन हमने actually अभी तक container specify नहीं किया है (ओह)।

### 2.3. `cowpy` process run करने के लिए container use करें

हमें एक container specify करना है और Nextflow को बताना है कि `cowpy()` process के लिए इसे use करे।

#### 2.3.1. `cowpy` के लिए container specify करें

हम same image use कर सकते हैं जो हम इस tutorial के first section में directly use कर रहे थे।

`cowpy.nf` module को edit करो process definition में `container` directive add करने के लिए इस प्रकार:

=== "बाद में"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

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

=== "पहले"

    ```groovy title="modules/cowpy.nf" linenums="4"
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

यह Nextflow को बताता है कि _यदि Docker का use enabled है_, तो process execute करने के लिए यहाँ specified container image use करना चाहिए।

#### 2.3.2. `nextflow.config` file के माध्यम से Docker use enable करें

ध्यान दो हमने कहा _'यदि Docker का use enabled है'_। Default रूप से, यह नहीं है, इसलिए हमें Nextflow को बताना होगा कि यह Docker use करने की अनुमति है।
इस उद्देश्य के लिए, हम इस course के next और last part (Part 6) के topic को थोड़ा anticipate करने जा रहे हैं, जो configuration cover करता है।

Workflow execution configure करने के लिए Nextflow जो main ways offer करता है उनमें से एक `nextflow.config` file use करना है।
जब ऐसी file current directory में present होती है, तो Nextflow automatically इसे load करेगा और कोई भी configuration apply करेगा जो इसमें contain है।

हमने एक `nextflow.config` file provide की एक single line of code के साथ जो explicitly Docker disable करती है: `docker.enabled = false`।

अब, चलो इसे `true` में switch करें Docker enable करने के लिए:

=== "बाद में"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "पहले"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = false
    ```

!!! tip

    Command-line से, per-run basis पर, `-with-docker <container>` parameter use करके Docker execution enable करना possible है।
    हालाँकि, वह केवल हमें entire workflow के लिए एक container specify करने की अनुमति देता है, जबकि जो approach हमने अभी तुम्हें दिखाया वह हमें per process एक अलग container specify करने की अनुमति देता है।
    यह modularity, code maintenance और reproducibility के लिए बेहतर है।

#### 2.3.3. Docker enabled के साथ workflow run करें

Workflow को `-resume` flag के साथ run करो:

```bash
nextflow run hello-containers.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [drunk_perlman] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [98/656c6c] cowpy              [100%] 1 of 1 ✔
    ```

इस बार यह वास्तव में work करता है!
हमेशा की तरह तुम corresponding results directory में workflow outputs पा सकते हो, हालाँकि इस बार वे थोड़े अधिक neatly organized हैं, केवल report और final output top level पर है, और सभी intermediate files एक subdirectory में out of the way shoved हैं।

??? abstract "डायरेक्टरी कंटेंट्स"

    ```console
    results/hello_containers/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

Final ASCII art output `results/hello_containers/` directory में है, `cowpy-COLLECTED-batch-output.txt` नाम के तहत।

??? abstract "फ़ाइल कंटेंट्स"

    ```console title="results/hello_containers/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

और वहाँ यह है, हमारी beautiful turkey greetings कहते हुए जैसा desired था।

#### 2.3.4. Inspect करें कि Nextflow ने containerized task कैसे launch की

इस section के final coda के रूप में, चलो `cowpy` process calls में से एक के लिए work subdirectory पर एक नज़र डालें ताकि थोड़ी और insight मिले कि Nextflow containers के साथ under the hood कैसे काम करता है।

अपने `nextflow run` command के output को check करो `cowpy` process के लिए work subdirectory का path खोजने के लिए।
ऊपर दिखाए गए run के लिए हमें जो मिला उसे देखते हुए, `cowpy` process के लिए console log line `[98/656c6c]` से शुरू होती है।
वह निम्नलिखित truncated directory path से correspond करती है: `work/98/656c6c`।

उस directory में, तुम्हें `.command.run` file मिलेगी जिसमें वे सभी commands contain हैं जो Nextflow ने pipeline execute करने के course में तुम्हारी ओर से run कीं।

??? abstract "फ़ाइल कंटेंट्स"

    ```console title="work/98/656c6c90cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/hello-nextflow/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY

    ```

यदि तुम इस file में `nxf_launch` search करो, तो तुम्हें कुछ ऐसा दिखना चाहिए:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
}
```

जैसा तुम देख सकते हो, Nextflow process call launch करने के लिए `docker run` command use कर रहा है।
यह corresponding work subdirectory को भी container में mount करता है, container के अंदर working directory accordingly set करता है, और `.command.sh` file में हमारी templated bash script run करता है।

वह सारा hard work जो हमें first section में manually करना पड़ा? Nextflow हमारे लिए behind the scenes करता है!

```txt
 _______________________
< Hurray for robots...! >
 -----------------------
                                   ,-----.
                                   |     |
                                ,--|     |-.
                         __,----|  |     | |
                       ,;::     |  `_____' |
                       `._______|    i^i   |
                                `----| |---'| .
                           ,-------._| |== ||//
                           |       |_|P`.  /'/
                           `-------' 'Y Y/'/'
                                     .==\ /_\
   ^__^                             /   /'|  `i
   (oo)\_______                   /'   /  |   |
   (__)\       )\/\             /'    /   |   `i
       ||----w |           ___,;`----'.___L_,-'`\__
       ||     ||          i_____;----\.____i""\____\
```

### सीख

तुम जानते हो कि Nextflow में processes run करने के लिए containers कैसे use करें।

### आगे क्या है?

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

और जानें: [1.3.4. Container में data mount करें](#134-container-में-data-mount-करें)
</quiz>

<quiz>
Containers use करते समय volumes mount करने की need क्यों होती है?
- [ ] Container performance improve करने के लिए
- [ ] Disk space save करने के लिए
- [x] क्योंकि containers default रूप से host filesystem से isolated हैं
- [ ] Networking enable करने के लिए

और जानें: [1.3.4. Container में data mount करें](#134-container-में-data-mount-करें)
</quiz>

<quiz>
Nextflow process के लिए container कैसे specify करते हो?
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

और जानें: [2.3.1. cowpy के लिए container specify करें](#231-cowpy-के-लिए-container-specify-करें)
</quiz>

<quiz>
कौन सी `nextflow.config` setting तुम्हारे workflow के लिए Docker enable करती है?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

और जानें: [2.3.2. `nextflow.config` file के माध्यम से Docker use enable करें](#232-nextflowconfig-file-के-माध्यम-से-docker-use-enable-करें)
</quiz>

<quiz>
Container में process run करते समय Nextflow automatically क्या handle करता है? (सभी लागू select करें)
- [x] यदि needed हो तो container image pull करना
- [x] Work directory mount करना
- [x] Container के अंदर process script run करना
- [x] Execution के बाद container instance clean up करना

और जानें: [2.3.4. Inspect करें कि Nextflow ने containerized task कैसे launch की](#234-inspect-करें-कि-nextflow-ने-containerized-task-कैसे-launch-की)
</quiz>
