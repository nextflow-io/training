# भाग 5: Hello Containers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [और जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [पूरी प्लेलिस्ट](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) Nextflow YouTube चैनल पर देखें।

:green_book: वीडियो ट्रांसक्रिप्ट [यहाँ](./transcripts/05_hello_containers.md) उपलब्ध है।
///

इस ट्रेनिंग कोर्स के भाग 1-4 में, तुमने Nextflow के बेसिक बिल्डिंग ब्लॉक्स का उपयोग करके एक सरल workflow बनाना सीखा जो कुछ टेक्स्ट को प्रोसेस कर सकता है, अगर कई इनपुट हों तो execution को parallelize कर सकता है, और आगे की प्रोसेसिंग के लिए परिणाम एकत्र कर सकता है।

हालांकि, तुम अपने environment में उपलब्ध बेसिक UNIX टूल्स तक सीमित थे।
वास्तविक दुनिया के कार्यों में अक्सर विभिन्न टूल्स और पैकेजों की आवश्यकता होती है जो डिफ़ॉल्ट रूप से शामिल नहीं होते।
आमतौर पर, तुम्हें इन टूल्स को इंस्टॉल करना होगा, उनकी dependencies को मैनेज करना होगा, और किसी भी conflicts को हल करना होगा।

यह सब बहुत थकाऊ और परेशान करने वाला है, इसलिए हम तुम्हें दिखाने जा रहे हैं कि इस समस्या को बहुत अधिक सुविधाजनक तरीके से हल करने के लिए **containers** का उपयोग कैसे करें।

एक **container** एक हल्की, स्वतंत्र, executable सॉफ़्टवेयर इकाई है जो एक container **image** से बनाई जाती है और जिसमें एप्लिकेशन चलाने के लिए आवश्यक सब कुछ शामिल होता है जिसमें code, system libraries और settings शामिल हैं।
जैसा कि तुम कल्पना कर सकते हो, यह तुम्हारी pipelines को अधिक reproducible बनाने के लिए बहुत मददगार होने वाला है।

ध्यान दो कि हम इसे [Docker](https://www.docker.com/get-started/) का उपयोग करके सिखा रहे हैं, लेकिन याद रखो कि Nextflow [कई अन्य container technologies](https://nextflow.io/docs/latest/container.html) को भी सपोर्ट करता है।

??? info "इस सेक्शन से कैसे शुरू करें"

    कोर्स का यह सेक्शन मानता है कि तुमने [Hello Nextflow](./index.md) कोर्स के भाग 1-4 पूरे कर लिए हैं और तुम्हारे पास एक पूर्ण working pipeline है।

    अगर तुम कोर्स को इस बिंदु से शुरू कर रहे हो, तो तुम्हें solutions से `modules` डायरेक्टरी को कॉपी करना होगा:

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. वार्मअप: `hello-containers.nf` चलाओ

हम शुरुआती बिंदु के रूप में workflow स्क्रिप्ट `hello-containers.nf` का उपयोग करने जा रहे हैं।
यह इस ट्रेनिंग कोर्स के भाग 4 को पूरा करके बनाई गई स्क्रिप्ट के बराबर है, सिवाय इसके कि हमने आउटपुट destinations बदल दिए हैं:

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

बस यह सुनिश्चित करने के लिए कि सब कुछ काम कर रहा है, कोई भी बदलाव करने से पहले स्क्रिप्ट को एक बार चलाओ:

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

पहले की तरह, तुम्हें आउटपुट फ़ाइलें `output` ब्लॉक में निर्दिष्ट डायरेक्टरी (`results/hello_containers/`) में मिलेंगी।

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

अगर यह तुम्हारे लिए काम कर गया, तो तुम containers का उपयोग करना सीखने के लिए तैयार हो।

---

## 1. Container को 'manually' उपयोग करो

हम अपने workflow में एक स्टेप जोड़ना चाहते हैं जो execution के लिए एक container का उपयोग करेगा।

हालांकि, हम पहले कुछ बेसिक concepts और operations पर जाने वाले हैं ताकि containers क्या हैं इसकी तुम्हारी समझ मजबूत हो, इससे पहले कि हम Nextflow में उनका उपयोग करना शुरू करें।

### 1.1. Container image को pull करो

एक container का उपयोग करने के लिए, तुम आमतौर पर एक container registry से एक container image को डाउनलोड या _pull_ करते हो, और फिर एक container instance बनाने के लिए container image को चलाते हो।

सामान्य syntax इस प्रकार है:

```bash title="Syntax"
docker pull '<container>'
```

`docker pull` भाग container system को एक repository से container image pull करने का निर्देश है।

`'<container>'` भाग container image का URI address है।

एक उदाहरण के रूप में, चलो एक container image pull करें जिसमें [cowpy](https://github.com/jeffbuttars/cowpy) है, एक python implementation जो `cowsay` नामक टूल का है जो मनमाने टेक्स्ट इनपुट को मज़ेदार तरीके से प्रदर्शित करने के लिए ASCII art जेनरेट करता है।

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
हमने [Seqera Containers](https://seqera.io/containers/) सेवा का उपयोग करके `cowpy` Conda पैकेज से यह Docker container image जेनरेट किया: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`।

पूरी pull कमांड चलाओ:

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

अगर तुमने पहले कभी image डाउनलोड नहीं की है, तो इसे पूरा होने में एक मिनट लग सकता है।
एक बार यह हो जाने के बाद, तुम्हारे पास container image की एक local कॉपी है।

### 1.2. `cowpy` को one-off कमांड के रूप में चलाने के लिए container का उपयोग करो

एक बहुत सामान्य तरीका जिससे लोग containers का उपयोग करते हैं वह है उन्हें सीधे चलाना, _यानी_ non-interactively।
यह one-off कमांड चलाने के लिए बहुत अच्छा है।

सामान्य syntax इस प्रकार है:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

`docker run --rm '<container>'` भाग container system को एक container image से एक container instance spin up करने और उसमें एक कमांड execute करने का निर्देश है।
`--rm` flag system को कमांड पूरा होने के बाद container instance को shut down करने के लिए कहता है।

`[tool command]` syntax उस टूल पर निर्भर करता है जिसका तुम उपयोग कर रहे हो और container कैसे सेट अप है।
चलो बस `cowpy` से शुरू करते हैं।

पूरी तरह से assembled, container execution कमांड इस तरह दिखती है; आगे बढ़ो और इसे चलाओ।

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

System ने container को spin up किया, `cowpy` कमांड को उसके parameters के साथ चलाया, आउटपुट को console पर भेजा और अंत में, container instance को shut down कर दिया।

### 1.3. `cowpy` को interactively चलाने के लिए container का उपयोग करो

तुम एक container को interactively भी चला सकते हो, जो तुम्हें container के अंदर एक shell prompt देता है और तुम्हें कमांड के साथ खेलने की अनुमति देता है।

#### 1.3.1. Container को spin up करो

Interactively चलाने के लिए, हम बस `docker run` कमांड में `-it` जोड़ते हैं।
वैकल्पिक रूप से, हम container के अंदर जिस shell का उपयोग करना चाहते हैं उसे निर्दिष्ट कर सकते हैं, _जैसे_ कमांड के अंत में `/bin/bash` जोड़कर।

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

ध्यान दो कि तुम्हारा prompt `(base) root@b645838b3314:/tmp#` जैसा कुछ बदल जाता है, जो दर्शाता है कि तुम अब container के अंदर हो।

तुम filesystem के root से डायरेक्टरी contents को list करने के लिए `ls /` चलाकर इसे verify कर सकते हो:

```bash
ls /
```

??? abstract "कमांड आउटपुट"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

हम यहाँ `tree` के बजाय `ls` का उपयोग करते हैं क्योंकि `tree` utility इस container में उपलब्ध नहीं है।
तुम देख सकते हो कि container के अंदर की filesystem तुम्हारे host system की filesystem से अलग है।

हमने अभी जो किया उसकी एक सीमा यह है कि container डिफ़ॉल्ट रूप से host system से पूरी तरह से isolated है।
इसका मतलब है कि container host system पर किसी भी फ़ाइल को access नहीं कर सकता जब तक कि तुम स्पष्ट रूप से इसे ऐसा करने की अनुमति नहीं देते।

हम तुम्हें एक मिनट में दिखाएंगे कि यह कैसे करना है।

#### 1.3.2. वांछित टूल कमांड(s) चलाओ

अब जब तुम container के अंदर हो, तो तुम `cowpy` कमांड को सीधे चला सकते हो और इसे कुछ parameters दे सकते हो।
उदाहरण के लिए, टूल documentation कहता है कि हम `-c` के साथ character ('cowacter') बदल सकते हैं।

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

अब आउटपुट डिफ़ॉल्ट cow के बजाय Linux penguin, Tux को दिखाता है, क्योंकि हमने `-c tux` parameter निर्दिष्ट किया।

क्योंकि तुम container के अंदर हो, तुम `cowpy` कमांड को जितनी बार चाहो उतनी बार चला सकते हो, इनपुट parameters को vary करते हुए, Docker कमांड्स की परेशानी के बिना।

!!! Tip

    एक अलग character चुनने के लिए '-c' flag का उपयोग करो, जिसमें शामिल हैं:
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

यह अच्छा है। और भी अच्छा होगा अगर हम अपनी `greetings.csv` को इसमें इनपुट के रूप में feed कर सकें।
लेकिन चूंकि हमारे पास filesystem तक access नहीं है, हम नहीं कर सकते।

चलो इसे ठीक करते हैं।

#### 1.3.3. Container से exit करो

Container से exit करने के लिए, तुम prompt पर `exit` टाइप कर सकते हो या ++ctrl+d++ keyboard shortcut का उपयोग कर सकते हो।

```bash
exit
```

तुम्हारा prompt अब वापस वैसा होना चाहिए जैसा container शुरू करने से पहले था।

#### 1.3.4. Container में data को mount करो

जैसा कि पहले बताया गया है, container डिफ़ॉल्ट रूप से host system से isolated है।

Container को host filesystem को access करने की अनुमति देने के लिए, तुम निम्नलिखित syntax का उपयोग करके host system से container में एक **volume** को **mount** कर सकते हो:

```bash title="Syntax"
-v <outside_path>:<inside_path>
```

हमारे मामले में `<outside_path>` वर्तमान working डायरेक्टरी होगी, इसलिए हम बस एक dot (`.`) का उपयोग कर सकते हैं, और `<inside_path>` बस एक alias है जो हम बनाते हैं; चलो इसे `/my_project` कहते हैं (inside path absolute होना चाहिए)।

एक volume mount करने के लिए, हम paths को replace करते हैं और docker run कमांड में volume mounting argument को इस प्रकार जोड़ते हैं:

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

यह वर्तमान working डायरेक्टरी को एक volume के रूप में mount करता है जो container के अंदर `/my_project` के तहत accessible होगी।

तुम `/my_project` की contents को list करके जांच सकते हो कि यह काम करता है:

```bash
ls /my_project
```

??? success "कमांड आउटपुट"

    ```console
    data               hello-config.nf      hello-modules.nf   hello-world.nf  nextflow.config  solutions         work
    hello-channels.nf  hello-containers.nf  hello-workflow.nf  modules         results          test-params.json
    ```

तुम अब container के अंदर से working डायरेक्टरी की contents देख सकते हो, जिसमें `data/` के तहत `greetings.csv` फ़ाइल शामिल है।

इसने effectively container wall के माध्यम से एक tunnel स्थापित की जिसका उपयोग तुम अपने filesystem के उस हिस्से को access करने के लिए कर सकते हो।

#### 1.3.5. Mounted data का उपयोग करो

अब जब हमने working डायरेक्टरी को container में mount कर दिया है, तो हम `greetings.csv` फ़ाइल की contents को प्रदर्शित करने के लिए `cowpy` कमांड का उपयोग कर सकते हैं।

ऐसा करने के लिए, हम CSV फ़ाइल की contents को `cowpy` कमांड में pipe करने के लिए `cat /my_project/data/greetings.csv | ` का उपयोग करेंगे।

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
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
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

यह हमारे उदाहरण greetings को बोलते हुए एक turkey की वांछित ASCII art उत्पन्न करता है!
सिवाय यहाँ turkey पूरी rows को repeat कर रहा है बजाय सिर्फ greetings के।
हम पहले से जानते हैं कि हमारा Nextflow workflow बेहतर काम करेगा!

इस कमांड के साथ खेलने के लिए स्वतंत्र महसूस करो।
जब तुम कर लो, तो पहले की तरह container से exit करो:

```bash
exit
```

तुम अपने सामान्य shell में वापस आ जाओगे।

### सारांश

तुम जानते हो कि एक container को कैसे pull करना है और इसे one-off या interactively कैसे चलाना है। तुम यह भी जानते हो कि अपने data को अपने container के अंदर से कैसे accessible बनाना है, जो तुम्हें अपने system पर कोई सॉफ़्टवेयर इंस्टॉल किए बिना वास्तविक data पर किसी भी टूल को आज़माने देता है जिसमें तुम रुचि रखते हो।

### आगे क्या है?

सीखो कि Nextflow processes के execution के लिए containers का उपयोग कैसे करें।

---

## 2. Nextflow में containers का उपयोग करो

Nextflow में processes को containers के अंदर चलाने के लिए built-in सपोर्ट है ताकि तुम ऐसे टूल्स चला सको जो तुम्हारे compute environment में इंस्टॉल नहीं हैं।
इसका मतलब है कि तुम अपनी processes को चलाने के लिए किसी भी container image का उपयोग कर सकते हो, और Nextflow image को pull करने, data को mount करने, और इसके अंदर process को चलाने का ध्यान रखेगा।

इसे demonstrate करने के लिए, हम `collectGreetings` स्टेप के बाद, हमारे द्वारा विकसित की जा रही pipeline में एक `cowpy` स्टेप जोड़ने जा रहे हैं।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

### 2.1. एक `cowpy` मॉड्यूल लिखो

पहले, चलो `cowpy` process मॉड्यूल बनाते हैं।

#### 2.1.1. नए मॉड्यूल के लिए एक फ़ाइल stub बनाओ

मॉड्यूल के लिए `cowpy.nf` नाम की एक खाली फ़ाइल बनाओ।

```bash
touch modules/cowpy.nf
```

यह हमें process code डालने के लिए एक जगह देता है।

#### 2.1.2. मॉड्यूल फ़ाइल में `cowpy` process code को कॉपी करो

हम अपनी `cowpy` process को उन अन्य processes पर model कर सकते हैं जो हमने पहले लिखी हैं।

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// cowpy के साथ ASCII art जेनरेट करें
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

Process एक `input_file` की अपेक्षा करती है जिसमें greetings हैं साथ ही एक `character` value।

आउटपुट एक नई टेक्स्ट फ़ाइल होगी जिसमें `cowpy` टूल द्वारा जेनरेट की गई ASCII art होगी।

### 2.2. Workflow में cowpy जोड़ो

अब हमें मॉड्यूल को import करना है और process को call करना है।

#### 2.2.1. `hello-containers.nf` में `cowpy` process को import करो

Workflow ब्लॉक के ऊपर import declaration डालो और इसे उचित रूप से भरो।

=== "बाद में"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="5"
    // मॉड्यूल्स को include करें
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

=== "पहले"

    ```groovy title="hello-containers.nf" linenums="3"
    // मॉड्यूल्स को include करें
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    ```

अब `cowpy` मॉड्यूल workflow में उपयोग के लिए उपलब्ध है।

#### 2.2.2. Workflow में `cowpy` process को call जोड़ो

चलो `cowpy()` process को `collectGreetings()` process के आउटपुट से connect करते हैं, जो जैसा कि तुम्हें याद होगा दो आउटपुट उत्पन्न करता है:

- `collectGreetings.out.outfile` में आउटपुट फ़ाइल है <--_जो हम चाहते हैं_
- `collectGreetings.out.report` में प्रति batch greetings की count के साथ report फ़ाइल है

Workflow ब्लॉक में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-containers.nf" linenums="19" hl_lines="12-13"
        main:
        // CSV फ़ाइल से इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
        // सभी अभिवादनों को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // cowpy के साथ अभिवादनों की ASCII art जेनरेट करें
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "पहले"

    ```groovy title="hello-containers.nf" linenums="19"
        main:
        // CSV फ़ाइल से इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
        // सभी अभिवादनों को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

ध्यान दो कि हमने एक नया CLI parameter, `params.character`, declare किया ताकि हम निर्दिष्ट कर सकें कि हम किस character को greetings कहते हुए चाहते हैं।

#### 2.2.3. `params` ब्लॉक में `character` parameter जोड़ो

यह technically वैकल्पिक है लेकिन यह recommended practice है और यह character के लिए एक डिफ़ॉल्ट value सेट करने का एक अवसर है जबकि हम इस पर हैं।

=== "बाद में"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="7"
    /*
    * Pipeline पैरामीटर्स
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
    * Pipeline पैरामीटर्स
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

अब हम आलसी हो सकते हैं और अपनी कमांड lines में character parameter टाइप करना skip कर सकते हैं।

#### 2.2.4. Workflow आउटपुट को अपडेट करो

हमें `cowpy` process के आउटपुट को publish करने के लिए workflow आउटपुट को अपडेट करना होगा।

##### 2.2.4.1. `publish:` सेक्शन को अपडेट करो

`workflow block` में, निम्नलिखित code change करो:

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

`cowpy` process केवल एक आउटपुट उत्पन्न करती है इसलिए हम इसे सामान्य तरीके से `.out` जोड़कर refer कर सकते हैं।

लेकिन अभी के लिए, चलो workflow-level आउटपुट को अपडेट करना समाप्त करते हैं।

##### 2.2.4.2. `output` ब्लॉक को अपडेट करो

हमें अंतिम `cowpy_art` आउटपुट को `output` ब्लॉक में जोड़ना होगा। जबकि हम इस पर हैं, चलो publishing destinations को भी edit करते हैं क्योंकि अब हमारी pipeline पूरी हो गई है और हम जानते हैं कि हमें वास्तव में किन आउटपुट की परवाह है।

`output` ब्लॉक में, निम्नलिखित code changes करो:

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

अब published आउटपुट थोड़े अधिक organized होंगे।

#### 2.2.5. Workflow को चलाओ

बस recap करने के लिए, यह वह है जिसका हम लक्ष्य रख रहे हैं:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

क्या तुम्हें लगता है कि यह काम करने वाला है?

चलो एक clean slate रखने के लिए पिछले published आउटपुट को delete करते हैं, और `-resume` flag के साथ workflow चलाते हैं।

```bash
rm -r hello_containers/
nextflow run hello-containers.nf -resume
```

??? failure "कमांड आउटपुट (स्पष्टता के लिए संपादित)"

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

अरे नहीं, एक error है!
`error exit status (127)` द्वारा दिया गया error code का मतलब है कि हमने जो executable माँगा वह नहीं मिला।

यह समझ में आता है, क्योंकि हम `cowpy` टूल को call कर रहे हैं लेकिन हमने वास्तव में अभी तक एक container निर्दिष्ट नहीं किया है (oops)।

### 2.3. `cowpy` process को चलाने के लिए एक container का उपयोग करो

हमें एक container निर्दिष्ट करना होगा और Nextflow को बताना होगा कि `cowpy()` process के लिए इसका उपयोग करें।

#### 2.3.1. `cowpy` के लिए एक container निर्दिष्ट करो

हम वही image उपयोग कर सकते हैं जो हम इस tutorial के पहले सेक्शन में सीधे उपयोग कर रहे थे।

`cowpy.nf` मॉड्यूल को edit करो ताकि process definition में `container` directive को इस प्रकार जोड़ा जा सके:

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

यह Nextflow को बताता है कि _अगर Docker का उपयोग enabled है_, तो इसे process को execute करने के लिए यहाँ निर्दिष्ट container image का उपयोग करना चाहिए।

#### 2.3.2. `nextflow.config` फ़ाइल के माध्यम से Docker के उपयोग को enable करो

ध्यान दो कि हमने कहा _'अगर Docker का उपयोग enabled है'_। डिफ़ॉल्ट रूप से, यह नहीं है, इसलिए हमें Nextflow को बताना होगा कि इसे Docker का उपयोग करने की अनुमति है।
इस उद्देश्य के लिए, हम इस कोर्स के अगले और अंतिम भाग (भाग 6) के विषय को थोड़ा anticipate करने जा रहे हैं, जो configuration को cover करता है।

Nextflow workflow execution को configure करने के लिए जो मुख्य तरीके प्रदान करता है उनमें से एक `nextflow.config` फ़ाइल का उपयोग करना है।
जब ऐसी फ़ाइल वर्तमान डायरेक्टरी में मौजूद होती है, तो Nextflow स्वचालित रूप से इसे load करेगा और इसमें मौजूद किसी भी configuration को apply करेगा।

हमने एक `nextflow.config` फ़ाइल प्रदान की जिसमें code की एक single line है जो स्पष्ट रूप से Docker को disable करती है: `docker.enabled = false`।

अब, चलो Docker को enable करने के लिए इसे `true` में बदलते हैं:

=== "बाद में"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "पहले"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = false
    ```

!!! tip

    कमांड-लाइन से, प्रति-रन आधार पर, `-with-docker <container>` parameter का उपयोग करके Docker execution को enable करना संभव है।
    हालांकि, यह हमें केवल पूरे workflow के लिए एक container निर्दिष्ट करने की अनुमति देता है, जबकि हमने अभी जो approach दिखाया वह हमें प्रति process एक अलग container निर्दिष्ट करने की अनुमति देता है।
    यह modularity, code maintenance और reproducibility के लिए बेहतर है।

#### 2.3.3. Docker enabled के साथ workflow चलाओ

`-resume` flag के साथ workflow चलाओ:

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

इस बार यह वास्तव में काम करता है!
हमेशा की तरह तुम संबंधित results डायरेक्टरी में workflow आउटपुट पा सकते हो, हालांकि इस बार वे थोड़े अधिक neatly organized हैं, केवल report और final आउटपुट top level पर हैं, और सभी intermediate फ़ाइलें एक subdirectory में रास्ते से बाहर shoved हैं।

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

Final ASCII art आउटपुट `results/hello_containers/` डायरेक्टरी में है, `cowpy-COLLECTED-batch-output.txt` नाम के तहत।

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
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
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

और वहाँ यह है, हमारा सुंदर turkey वांछित रूप से greetings कह रहा है।

#### 2.3.4. Inspect करो कि Nextflow ने containerized task को कैसे launch किया

इस सेक्शन के अंतिम coda के रूप में, चलो `cowpy` process calls में से एक के लिए work subdirectory पर एक नज़र डालते हैं ताकि Nextflow containers के साथ hood के नीचे कैसे काम करता है इस पर थोड़ी अधिक insight मिल सके।

`cowpy` process के लिए work subdirectory का path खोजने के लिए अपने `nextflow run` कमांड से आउटपुट check करो।
ऊपर दिखाए गए रन के लिए हमें जो मिला उसे देखते हुए, `cowpy` process के लिए console log line `[98/656c6c]` से शुरू होती है।
यह निम्नलिखित truncated डायरेक्टरी path से मेल खाती है: `work/98/656c6c`।

उस डायरेक्टरी में, तुम्हें `.command.run` फ़ाइल मिलेगी जिसमें वे सभी कमांड हैं जो Nextflow ने pipeline को execute करने के दौरान तुम्हारी ओर से चलाए।

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
        # इनपुट फ़ाइलों को stage करें
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

अगर तुम इस फ़ाइल में `nxf_launch` search करते हो, तो तुम्हें कुछ ऐसा दिखना चाहिए:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
}
```

जैसा कि तुम देख सकते हो, Nextflow process call को launch करने के लिए `docker run` कमांड का उपयोग कर रहा है।
यह संबंधित work subdirectory को भी container में mount करता है, container के अंदर working डायरेक्टरी को तदनुसार सेट करता है, और `.command.sh` फ़ाइल में हमारी templated bash स्क्रिप्ट को चलाता है।

वह सारा कठिन काम जो हमें पहले सेक्शन में manually करना पड़ा? Nextflow हमारे लिए पर्दे के पीछे करता है!

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

### सारांश

तुम जानते हो कि processes को चलाने के लिए Nextflow में containers का उपयोग कैसे करें।

### आगे क्या है?

एक ब्रेक लो!

जब तुम तैयार हो, तो [**भाग 6: Hello Config**](./06_hello_config.md) पर जाओ ताकि सीख सको कि अपनी pipeline के execution को अपने infrastructure के अनुसार configure कैसे करें साथ ही inputs और parameters के configuration को कैसे मैनेज करें।

यह बिल्कुल अंतिम भाग है, और फिर तुम इस कोर्स के साथ समाप्त हो जाओगे!

---

## क्विज़

<quiz>
एक container क्या है?
- [ ] एक प्रकार की virtual machine
- [ ] एक फ़ाइल compression format
- [x] एक हल्की, स्वतंत्र executable इकाई जिसमें एप्लिकेशन चलाने के लिए आवश्यक सब कुछ शामिल है
- [ ] एक network protocol
</quiz>

<quiz>
एक container image और एक container instance के बीच क्या अंतर है?
- [ ] वे एक ही चीज़ हैं
- [x] एक image एक template है; एक instance उस image से बनाया गया एक running container है
- [ ] एक instance एक template है; एक image एक running container है
- [ ] Images Docker के लिए हैं; instances Singularity के लिए हैं
</quiz>

<quiz>
`docker run` कमांड में `-v` flag क्या करता है?
- [ ] Verbose आउटपुट enable करता है
- [ ] Container को validate करता है
- [x] Host system से container में एक volume mount करता है
- [ ] Container के version को निर्दिष्ट करता है

और जानें: [1.3.4. Container में data को mount करो](#134-container-में-data-को-mount-करो)
</quiz>

<quiz>
Containers का उपयोग करते समय तुम्हें volumes mount करने की आवश्यकता क्यों है?
- [ ] Container performance को बेहतर बनाने के लिए
- [ ] Disk space बचाने के लिए
- [x] क्योंकि containers डिफ़ॉल्ट रूप से host filesystem से isolated हैं
- [ ] Networking enable करने के लिए

और जानें: [1.3.4. Container में data को mount करो](#134-container-में-data-को-mount-करो)
</quiz>

<quiz>
तुम एक Nextflow process के लिए एक container कैसे निर्दिष्ट करते हो?
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

और जानें: [2.3.1. cowpy के लिए एक container निर्दिष्ट करो](#231-cowpy-के-लिए-एक-container-निर्दिष्ट-करो)
</quiz>

<quiz>
कौन सी `nextflow.config` setting तुम्हारे workflow के लिए Docker को enable करती है?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

और जानें: [2.3.2. `nextflow.config` फ़ाइल के माध्यम से Docker के उपयोग को enable करो](#232-nextflowconfig-फ़ाइल-के-माध्यम-से-docker-के-उपयोग-को-enable-करो)
</quiz>

<quiz>
जब एक container में एक process चलाते हैं तो Nextflow स्वचालित रूप से क्या handle करता है? (सभी लागू का चयन करें)
- [x] आवश्यकता होने पर container image को pull करना
- [x] Work डायरेक्टरी को mount करना
- [x] Container के अंदर process स्क्रिप्ट को चलाना
- [x] Execution के बाद container instance को clean up करना

और जानें: [2.3.4. Inspect करो कि Nextflow ने containerized task को कैसे launch किया](#234-inspect-करो-कि-nextflow-ने-containerized-task-को-कैसे-launch-किया)
</quiz>
