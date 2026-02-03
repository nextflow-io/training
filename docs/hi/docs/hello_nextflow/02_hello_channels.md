# भाग 2: Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube channel पर [पूरी playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) देखें।

:green_book: वीडियो transcript [यहाँ](./transcripts/02_hello_channels.md) उपलब्ध है।
///
-->

इस course के Part 1 (Hello World) में, हमने तुम्हें दिखाया कि process call में directly input provide करके process को variable input कैसे provide करें: `sayHello(params.input)`।
यह जानबूझकर simplified approach था।
Practice में, उस approach की major limitations हैं; namely कि यह केवल बहुत simple cases के लिए काम करता है जहाँ हम process को केवल एक बार, single value पर run करना चाहते हैं।
अधिकांश realistic workflow use cases में, हम multiple values (उदाहरण के लिए, multiple samples के लिए experimental data) process करना चाहते हैं, इसलिए हमें inputs handle करने का एक अधिक sophisticated तरीका चाहिए।

इसके लिए Nextflow [**channels**](https://nextflow.io/docs/latest/channel.html) हैं।
Channels ऐसी queues हैं जो inputs को efficiently handle करने और उन्हें multi-step workflows में एक step से दूसरे में shuttle करने के लिए designed हैं, जबकि built-in parallelism और कई additional benefits provide करते हैं।

इस course के इस भाग में, तुम सीखोगे कि विभिन्न स्रोतों से multiple inputs handle करने के लिए channel कैसे use करें।
तुम channel contents को आवश्यकतानुसार transform करने के लिए [**operators**](https://nextflow.io/docs/latest/reference/operator.html) use करना भी सीखोगे।

??? info "इस section से कैसे शुरू करें"

    Course का यह section मानता है कि तुमने [Hello Nextflow](./index.md) course का Part 1 complete कर लिया है, लेकिन यदि तुम उस section में covered basics से comfortable हो, तो तुम बिना कुछ special किए यहाँ से शुरू कर सकते हो।

---

## 0. Warmup: `hello-channels.nf` चलाएं

हम starting point के रूप में workflow script `hello-channels.nf` use करेंगे।
यह इस training course के Part 1 में काम करके produce की गई script के equivalent है, सिवाय इसके कि हमने output destination बदल दी है:

```groovy title="hello-channels.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
```

यह sure करने के लिए कि सब कुछ काम कर रहा है, कोई भी changes करने से पहले script को एक बार run करो:

```bash
nextflow run hello-channels.nf --input 'Hello Channels!'
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [wise_jennings] DSL2 - revision: b24f4902d6

    executor >  local (1)
    [6f/824bc1] process > sayHello [100%] 1 of 1 ✔
    ```

पहले की तरह, तुम `results/hello_channels` directory में `output.txt` नामक output file पाओगे (जैसा कि workflow script के `output` block में specify किया गया है, ऊपर दिखाया गया है)।

??? abstract "Directory contents"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "File contents"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

यदि यह तुम्हारे लिए काम किया, तो तुम channels के बारे में सीखने के लिए ready हो।

---

## 1. Channel के माध्यम से explicitly variable inputs provide करें

हम implicit handling पर rely करने के बजाय `sayHello()` process को variable input pass करने के लिए एक **channel** बनाएंगे, जिसकी certain limitations हैं।

### 1.1. Input channel बनाएं

Channel set up करने के लिए हम विभिन्न प्रकार के [**channel factories**](https://nextflow.io/docs/latest/reference/channel.html) use कर सकते हैं।
अभी के लिए चीजों को simple रखने के लिए, हम सबसे basic channel factory use करेंगे, जिसे [`channel.of`](https://nextflow.io/docs/latest/reference/channel.html#of) कहते हैं, जो single value वाला channel बनाएगा।
Functionally यह पहले जैसे set up के similar होगा, लेकिन Nextflow को implicitly channel बनाने देने के बजाय, अब हम यह explicitly कर रहे हैं।

यह code की वह line है जिसे हम use करेंगे:

```console title="Syntax"
greeting_ch = channel.of('Hello Channels!')
```

यह `channel.of()` channel factory का उपयोग करके `greeting_ch` नामक एक channel बनाता है, जो एक simple queue channel set up करता है, और greeting value के रूप में use करने के लिए string `'Hello Channels!'` load करता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note "नोट"

    हम readability के लिए temporarily CLI parameter use करने के बजाय hardcoded strings पर वापस switch कर रहे हैं। Channel के level पर जो हो रहा है उसे cover करने के बाद हम CLI parameters use करने पर वापस जाएंगे।

Workflow block में, channel factory code add करो:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of('Hello Channels!')
        // एक अभिवादन emit करें
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // एक अभिवादन emit करें
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

यह अभी functional नहीं है क्योंकि हमने अभी तक process call को input switch नहीं किया है।

### 1.2. Process call में input के रूप में channel add करें

अब हमें अपने newly created channel को `sayHello()` process call में plug करना होगा, जिस CLI parameter को हम पहले directly provide कर रहे थे उसे replace करते हुए।

Workflow block में, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of('Hello Channels!')
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of('Hello Channels!')
        // एक अभिवादन emit करें
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

यह Nextflow को बताता है कि `greeting_ch` channel की contents पर `sayHello` process run करे।

अब हमारा workflow properly functional है; यह `sayHello('Hello Channels!')` लिखने का explicit equivalent है।

### 1.3. Workflow चलाएं

चलो इसे चलाते हैं!

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [fabulous_crick] DSL2 - revision: 23e20f76e8

    executor >  local (1)
    [c0/4f1872] process > sayHello (1) [100%] 1 of 1 ✔
    ```

यदि तुमने दोनों edits correctly किए, तो तुम्हें एक successful execution मिलनी चाहिए।
तुम results directory check कर सकते हो यह satisfy करने के लिए कि outcome अभी भी पहले जैसा ही है।

??? abstract "File contents"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

तो हमने same end result achieve करते हुए अपने workflow की flexibility बढ़ा दी है।
यह ऐसा लग सकता है कि हम बिना किसी tangible benefit के more code लिख रहे हैं, लेकिन जैसे ही हम more inputs handle करना शुरू करेंगे value clear हो जाएगी।

उसके preview के रूप में, move on करने से पहले एक और चीज़ देखते हैं: data input manage करने के लिए explicit channel use करने का एक छोटा लेकिन convenient benefit।

### 1.4. Channel contents inspect करने के लिए `view()` use करें

Nextflow channels इस तरह built हैं कि हम operators का उपयोग करके उनकी contents पर operate कर सकते हैं, जिसे हम इस chapter में बाद में detail में cover करेंगे।

अभी के लिए, हम तुम्हें बस दिखाएंगे कि channel की contents inspect करने के लिए [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) नामक एक super simple operator कैसे use करें।
तुम `view()` को एक debugging tool के रूप में सोच सकते हो, जैसे Python में `print()` statement, या अन्य languages में इसके equivalent।

Workflow block में यह tiny line add करो:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of('Hello Channels!')
                             .view()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of('Hello Channels!')
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Exact spaces की amount matter नहीं करती जब तक यह 4 का multiple है; हम बस `.view()` statement की start को channel construction के `.of()` part से align करने का aim कर रहे हैं।

अब workflow फिर से run करो:

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [scruffy_shaw] DSL2 - revision: 2ede41e14a

    executor >  local (1)
    [ef/f7e40a] sayHello (1) [100%] 1 of 1 ✔
    Hello Channels!
    ```

जैसा कि तुम देख सकते हो, यह channel contents को console पर output करता है।
यहाँ हमारे पास केवल एक element है, लेकिन जब हम next section में channel में multiple values load करना शुरू करेंगे, तुम देखोगे कि यह एक element per line output करने के लिए set है।

### सीख

तुम जानते हो कि process को input provide करने के लिए basic channel factory कैसे use करें।

### आगे क्या?

सीखो कि workflow को multiple input values पर iterate करने के लिए channels कैसे use करें।

---

## 2. Multiple input values पर run करने के लिए workflow modify करें

Workflows typically inputs के batches पर run होते हैं जो bulk में process होने के लिए meant हैं, इसलिए हम workflow को upgrade करना चाहते हैं ताकि multiple input values accept कर सके।

### 2.1. Input channel में multiple greetings load करें

Conveniently, `channel.of()` channel factory जो हम use कर रहे हैं वह एक से अधिक value accept करने में काफी खुश है, इसलिए हमें उसे modify करने की need नहीं है।
हम बस channel में multiple values load कर सकते हैं।

चलो उन्हें `'Hello'`, `'Bonjour'` और `'Holà'` बनाते हैं।

#### 2.1.1. More greetings add करें

Workflow block से पहले, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // इनपुट के लिए एक channel बनाएं
    greeting_ch = channel.of('Hello','Bonjour','Holà')
                         .view()
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // इनपुट के लिए एक channel बनाएं
    greeting_ch = channel.of('Hello Channels')
                         .view()
    ```

Documentation हमें बताती है कि यह काम करना चाहिए। क्या यह सच में इतना simple हो सकता है?

#### 2.1.2. Command run करें और log output देखें

चलो try करते हैं।

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console hl_lines="6"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [amazing_crick] DSL2 - revision: 59a9a5888a

    executor >  local (3)
    [f4/c9962c] process > sayHello (1) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

यह certainly ठीक से run हुआ लगता है।
Execution monitor दिखाता है कि `sayHello` process के लिए `3 of 3` calls किए गए, और हम `view()` statement द्वारा enumerate किए गए तीन greetings देखते हैं, जैसा promised एक per line।

हालाँकि, results directory में अभी भी केवल एक output है:

??? abstract "Directory contents"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "File contents"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

तुम्हें वहाँ तीन greetings में से एक दिखनी चाहिए, लेकिन जो तुम्हें मिली वह यहाँ दिखाई गई से different हो सकती है।
क्या तुम सोच सकते हो कि ऐसा क्यों हो सकता है?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_Diagram में, channel को green में represent किया गया है, और elements का order pipe में marbles की तरह represent किया गया है: पहले load किया गया right पर है, फिर दूसरा middle में है, फिर तीसरा left पर है।_

Execution monitor को वापस देखते हुए, इसने हमें केवल एक subdirectory path (`f4/c9962c`) दी।
चलो वहाँ देखते हैं।

??? abstract "Directory contents"

    ```console hl_lines="9"
    work/f4/c9962ce91ef87480babcb86b2b9042/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "File contents"

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

यह तो वह greeting भी नहीं है जो हमें results directory में मिली! क्या हो रहा है?

इस point पर, हमें तुम्हें बताना होगा कि default रूप से, ANSI logging system same process के multiple calls से logging को same line पर लिखता है।
तो sayHello() process के तीनों calls से status same spot पर land कर रही है।

सौभाग्य से, हम process calls की full list देखने के लिए उस behavior को disable कर सकते हैं।

#### 2.1.3. `-ansi-log false` option के साथ command फिर से run करें

Logging को expand करके per process call एक line display करने के लिए, command में `-ansi-log false` add करो।

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "Command output"

    ```console
     N E X T F L O W  ~  version 25.10.2
    Launching `hello-channels.nf` [desperate_monod] DSL2 - revision: 59a9a5888a
    Hello
    Bonjour
    Holà
    [23/871c7e] Submitted process > sayHello (2)
    [7f/21e2c2] Submitted process > sayHello (1)
    [f4/ea10a6] Submitted process > sayHello (3)
    ```

इस बार हम output में listed तीनों process runs और उनके associated work subdirectories देखते हैं।

यह बहुत better है, कम से कम एक simple workflow के लिए।
एक complex workflow, या बड़ी संख्या में inputs के लिए, terminal पर full list output होने से थोड़ा overwhelming हो जाएगा।
इसीलिए `-ansi-log false` default behavior नहीं है।

!!! tip "सुझाव"

    Status report करने का तरीका दो logging modes के बीच थोड़ा different है।
    Condensed mode में, Nextflow report करता है कि calls successfully complete हुईं या नहीं।
    इस expanded mode में, यह केवल report करता है कि वे submitted की गईं।

Anyway, अब जबकि हमारे पास प्रत्येक process call की subdirectories हैं, हम उनके logs और outputs खोज सकते हैं।

??? abstract "Directory contents"

    ```console
    work/23/871c7ec3642a898ecd5e6090d21300/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/7f/21e2c2f3cc8833ef3858b236e5575c/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/f4/ea10a680d5687596d3eaa3fcf69272/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "File contents"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

यह दिखाता है कि तीनों processes successfully run हुईं (yay)।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-task-dirs.svg"
</figure>

उस ने कहा, हमारे पास अभी भी problem है कि results directory में केवल एक output file है।

तुम्हें याद होगा कि हमने `sayHello` process के लिए output file name hardcode किया था, इसलिए तीनों calls ने `output.txt` नामक file produce की।

जब तक output files work subdirectories में रहती हैं, अन्य processes से isolated, तब तक यह okay है।
लेकिन जब वे same results directory में publish होती हैं, जो भी पहले वहाँ copy की गई वह अगली द्वारा overwrite हो जाती है, और इसी तरह।

### 2.2. सुनिश्चित करें कि output file names unique होंगे

हम सभी outputs को same results directory में publish करना जारी रख सकते हैं, लेकिन हमें ensure करना होगा कि उनके unique names होंगे।
Specifically, हमें first process को dynamically file name generate करने के लिए modify करना होगा ताकि final file names unique हों।

तो हम file names को unique कैसे बनाएं?
ऐसा करने का एक common तरीका output file name के भाग के रूप में inputs (input channel से received) से कुछ unique metadata use करना है।
यहाँ, convenience के लिए, हम greeting itself use करेंगे क्योंकि यह बस एक short string है, और इसे base output filename से पहले prepend करेंगे।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi-unique.svg"
</figure>

#### 2.2.1. Dynamic output file name construct करें

Process block में, निम्नलिखित code changes करो:

=== "After"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }
    ```

Output definition और `script:` command block दोनों में `output.txt` replace करना sure करो।

!!! tip "सुझाव"

    Output definition में, तुम्हें output filename expression के आसपास double quotes use करना MUST है (single quotes नहीं), otherwise यह fail होगा।

यह हर बार process call होने पर एक unique output file name produce करेगा, ताकि इसे output directory में same process के अन्य calls के outputs से distinguish किया जा सके।

#### 2.2.2. Workflow चलाएं

चलो इसे run करते हैं। Note करो कि हम default ANSI log settings के साथ वापस running पर हैं।

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sharp_minsky] DSL2 - revision: 16a291febe

    executor >  local (3)
    [e8/33ee64] sayHello (2) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Summary view पर वापस आते हुए, output फिर से एक line पर summarize हो गया है।
यह देखने के लिए `results` directory पर नज़र डालो कि क्या सभी output greetings वहाँ हैं।

??? abstract "Directory contents"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

हाँ! और प्रत्येक में expected contents हैं।

??? abstract "File contents"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

Success! अब हम जितनी चाहें उतनी greetings add कर सकते हैं बिना output files के overwrite होने की चिंता किए।

!!! tip "सुझाव"

    Practice में, input data itself के आधार पर files name करना almost हमेशा impractical है।
    Dynamic filenames generate करने का better तरीका input files के साथ metadata को process में pass करना है।
    Metadata typically 'sample sheet' या equivalents के माध्यम से provide किया जाता है।
    तुम यह बाद में अपने Nextflow training में सीखोगे ([Metadata side quest](../side_quests/metadata.md) देखें)।

### सीख

तुम जानते हो कि channel के माध्यम से multiple input elements कैसे feed करें।

### आगे क्या?

सीखो कि channel की contents transform करने के लिए operator कैसे use करें।

---

## 3. Array के माध्यम से multiple inputs provide करें

हमने अभी तुम्हें दिखाया कि multiple input elements कैसे handle करें जो directly channel factory में hardcoded थे।
क्या होगा यदि हम उन multiple inputs को different way में provide करना चाहते?

उदाहरण के लिए, imagine करो कि हम इस तरह elements का array वाला एक input variable set up करते हैं:

`greetings_array = ['Hello','Bonjour','Holà']`

क्या हम उसे अपने output channel में load कर सकते हैं और expect कर सकते हैं कि यह काम करे?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

चलो पता लगाते हैं।

### 3.1. Channel को input के रूप में values का array provide करें

Common sense suggest करता है कि हमें single value के बजाय simply values का array pass करने में सक्षम होना चाहिए।
चलो try करते हैं; हमें input variable set up करना होगा और इसे channel factory में load करना होगा।

#### 3.1.1. Input variable set up करें

चलो `greetings_array` variable जो हमने अभी imagine किया उसे workflow block में add करके reality बनाते हैं:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // इनपुट अभिवादनों की एक array declare करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

यह अभी functional नहीं है, हमने बस array के लिए declaration add किया है।

#### 3.1.2. Channel factory को input के रूप में greetings का array set करें

अब हम channel factory में currently hardcoded values `'Hello','Bonjour','Holà'` को अभी बनाए गए `greetings_array` से replace करेंगे।

Workflow block में, निम्नलिखित change करो:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // इनपुट अभिवादनों की एक array declare करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of(greetings_array)
                             .view()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // इनपुट अभिवादनों की एक array declare करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

यह अब functional होना चाहिए।

#### 3.1.3. Workflow चलाएं

चलो इसे running try करते हैं:

```bash
nextflow run hello-channels.nf
```

??? failure "Command output"

    ```console hl_lines="7 11 16"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

    executor >  local (1)
    [a8/1f6ead] sayHello (1) | 0 of 1
    [Hello, Bonjour, Holà]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`


    Command executed:

      echo '[Hello, Bonjour, Holà]' > '[Hello, Bonjour, Holà]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/a8/1f6ead5f3fa30a3c508e2e7cf83ffb

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

अरे नहीं! Error है!

`view()` का output और error messages देखो।

ऐसा लगता है Nextflow ने single process call run करने की कोशिश की, `[Hello, Bonjour, Holà]` को string value के रूप में use करते हुए, array में तीन strings को separate values के रूप में use करने के बजाय।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-fail.svg"
</figure>

तो यह 'packaging' है जो problem cause कर रही है।
हम Nextflow को array unpack करवाकर individual strings को channel में load कैसे करवाएं?

### 3.2. Channel contents transform करने के लिए operator use करें

यहीं [**operators**](https://nextflow.io/docs/latest/reference/operator.html) play में आते हैं।
तुम पहले से `.view()` operator use कर चुके हो, जो बस देखता है कि वहाँ क्या है।
अब हम उन operators को देखेंगे जो हमें channel की contents पर act करने की अनुमति देते हैं।

यदि तुम Nextflow documentation में [operators की list](https://www.nextflow.io/docs/latest/reference/operator.html) skim through करते हो, तो तुम [`flatten()`](https://www.nextflow.io/docs/latest/reference/operator.html#flatten) पाओगे, जो exactly वही करता है जो हमें चाहिए: array की contents unpack करना और उन्हें individual items के रूप में emit करना।

#### 3.2.1. `flatten()` operator add करें

हमारे input channel पर `flatten()` operator apply करने के लिए, हम इसे channel factory declaration में append करते हैं।

Workflow block में, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // इनपुट अभिवादनों की एक array declare करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // इनपुट अभिवादनों की एक array declare करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of(greetings_array)
                             .view()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

यहाँ हमने readability के लिए operator को next line पर add किया, लेकिन तुम prefer करो तो operators को channel factory के same line पर add कर सकते हो, इस तरह:
`greeting_ch = channel.of(greetings_array).view().flatten()`

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-success.svg"
</figure>

#### 3.2.2. `view()` statement(s) refine करें

हम इसे test करने के लिए तुरंत run कर सकते हैं, लेकिन जब हम इस पर हैं, हम channel contents को inspect करने के तरीके को refine करेंगे।

हम यह contrast करने में सक्षम होना चाहते हैं कि `flatten()` operator apply होने से पहले और बाद में contents कैसी दिखती हैं, इसलिए हम दूसरा add करेंगे, AND हम उन्हें output में अधिक clearly labeled करने के लिए थोड़ा code add करेंगे।

Workflow block में, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // इनपुट अभिवादनों की एक array declare करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // इनपुट अभिवादनों की एक array declare करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

तुम देखोगे कि हमने दूसरा `.view` statement add किया है, और प्रत्येक के लिए, हमने empty parentheses (`()`) को curly braces से replace किया है जिसमें कुछ code है, जैसे `{ greeting -> "Before flatten: $greeting" }`।

इन्हें _closures_ कहते हैं। इनमें contained code channel में प्रत्येक item के लिए execute होगा।
हम inner value के लिए एक temporary variable define करते हैं, यहाँ `greeting` कहलाता है (लेकिन यह कोई भी arbitrary name हो सकता है), जो केवल उस closure के scope के भीतर use होता है।

इस example में, `$greeting` channel में load किए गए प्रत्येक individual item को represent करता है।
इसका result neatly labeled console output होगा।

!!! info "जानकारी"

    कुछ pipelines में तुम operator closures के अंदर `$it` नामक एक special variable देख सकते हो।
    यह एक _implicit_ variable है जो inner variable तक short-hand access की अनुमति देती है,
    बिना इसे `->` के साथ define करने की need के।

    हम clarity में मदद के लिए explicit होना prefer करते हैं, इसलिए `$it` syntax discouraged है और slowly Nextflow language से phase out हो जाएगी।

#### 3.2.3. Workflow चलाएं

Finally, तुम workflow को फिर से running try कर सकते हो!

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    Before flatten: [Hello, Bonjour, Holà]
    After flatten: Hello
    After flatten: Bonjour
    After flatten: Holà
    ```

इस बार यह काम करता है AND हमें `flatten()` operator run करने से पहले और बाद में channel की contents कैसी दिखती हैं इसकी additional insight देता है।

- एक single `Before flatten:` statement क्योंकि उस point पर channel में एक item है, original array।
- तीन separate `After flatten:` statements, प्रत्येक greeting के लिए एक, जो अब channel में individual items हैं।

महत्वपूर्ण रूप से, इसका मतलब है कि प्रत्येक item अब workflow द्वारा separately process किया जा सकता है।

!!! tip "सुझाव"

    एक different channel factory, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist) use करके technically same results achieve करना संभव है, जिसमें इसके operation में एक implicit mapping step शामिल है।
    यहाँ हमने वह use न करने का choice किया ताकि एक simple use case पर operator के use को demonstrate किया जा सके।

### सीख

तुम जानते हो कि channel की contents transform करने के लिए `flatten()` जैसे operator कैसे use करें, और operator apply करने से पहले और बाद में channel contents inspect करने के लिए `view()` operator कैसे use करें।

### आगे क्या?

सीखो कि workflow को input values के source के रूप में file कैसे लेने दें।

---

## 4. CSV file से input values पढ़ें

Realistically, हम शायद ही कभी values के array से शुरू करेंगे।
Most likely, हमारे पास एक या अधिक files होंगी जिनमें वह data है जिसे process करना है, किसी प्रकार के structured format में।

हमने `greetings.csv` नामक एक CSV file prepare की है जिसमें कई input greetings हैं, उस तरह के columnar data को mimic करते हुए जो तुम real data analysis में process करना चाह सकते हो, `data/` के तहत stored।
(Numbers meaningful नहीं हैं, वे बस illustrative purposes के लिए हैं।)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

हमारा next task अपने workflow को इस file से values पढ़ने के लिए adapt करना है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

चलो देखते हैं कि हम यह कैसे कर सकते हैं।

### 4.1. Greetings के source के रूप में CSV file expect करने के लिए script modify करें

शुरू करने के लिए, हमें script में दो key changes करने होंगे:

- Input parameter को CSV file point करने के लिए switch करें
- Channel factory को file handle करने के लिए designed किसी में switch करें

#### 4.1.1. Input parameter को CSV file point करने के लिए switch करें

Part 1 में हमने जो `params.input` parameter set up किया था याद है?
हम इसे update करेंगे ताकि हमारी greetings वाली CSV file point करे।

Parameter declaration में निम्नलिखित edit करो:

=== "After"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * Pipeline पैरामीटर
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Pipeline parameters
     */
    input: String = 'Holà mundo!'
    ```

यह मानता है कि file workflow code के साथ co-located है।
तुम बाद में अपनी Nextflow journey में अन्य data locations के साथ deal करना सीखोगे।

#### 4.1.2. File handle करने के लिए designed channel factory में switch करें

चूंकि अब हम simple strings के बजाय file को input के रूप में use करना चाहते हैं, हम पहले वाली `channel.of()` channel factory use नहीं कर सकते।
हमें एक new channel factory, [`channel.fromPath()`](https://nextflow.io/docs/latest/reference/channel.html#frompath) use करने में switch करना होगा, जिसमें file paths handle करने के लिए कुछ built-in functionality है।

Workflow block में, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten() को uncomment करें
                             // .view { greeting -> "Flatten के बाद: $greeting" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // इनपुट अभिवादनों की एक array declare करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

तुम notice करोगे कि हमने channel input को `param.input` पर वापस switch किया, और `greetings_array` declaration delete कर दी क्योंकि अब हमें इसकी need नहीं होगी।
हमने `flatten()` और दूसरे `view()` statement को भी comment out कर दिया है।

#### 4.1.3. Workflow चलाएं

चलो new channel factory और input file के साथ workflow running try करते हैं।

```bash
nextflow run hello-channels.nf
```

??? failure "Command output"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    Before flatten: /workspaces/training/hello-nextflow/data/greetings.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f


    Command executed:

      echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.csv-output.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

अरे नहीं, यह काम नहीं करता। Console output और error message की start देखो।
`Command executed:` bit यहाँ विशेष रूप से helpful है।

यह थोड़ा familiar लग सकता है।
ऐसा लगता है Nextflow ने file path itself को string value के रूप में use करके single process call run करने की कोशिश की।
तो इसने file path correctly resolve किया है, लेकिन इसने actually इसकी contents parse नहीं की, जो हम चाहते थे।

हम Nextflow को file open करवाकर इसकी contents को channel में load कैसे करवाएं?

Sounds like हमें एक और [operator](https://www.nextflow.io/docs/latest/reference/operator.html) चाहिए!

### 4.2. File parse करने के लिए `splitCsv()` operator use करें

Operators की list को फिर से देखते हुए, हमें [`splitCsv()`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) मिलता है, जो CSV-formatted text को parse और split करने के लिए designed है।

#### 4.2.1. Channel पर `splitCsv()` apply करें

Operator apply करने के लिए, हम इसे पहले की तरह channel factory line में append करते हैं।

Workflow block में, `flatten()` को `splitcsv()` (uncommented) से replace करने के लिए निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten() को uncomment करें
                             // .view { greeting -> "Flatten के बाद: $greeting" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

जैसा तुम देख सकते हो, हमने before/after `view()` statements भी update किए हैं।
Technically हम same variable name (`greeting`) use कर सकते थे लेकिन हमने इसे कुछ अधिक appropriate (`csv`) में update किया ताकि code दूसरों द्वारा अधिक readable हो।

#### 4.2.2. Workflow फिर से चलाएं

चलो added CSV-parsing logic के साथ workflow running try करते हैं।

```bash
nextflow run hello-channels.nf
```

??? failure "Command output"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    ERROR ~ Error executing process > 'sayHello (2)'

    Caused by:
      Missing output file(s) `[Bonjour, French, 456]-output.txt` expected by process `sayHello (2)`


    Command executed:

      echo '[Bonjour, French, 456]' > '[Bonjour, French, 456]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/24/76da2fcc4876b61632749f99e26a50

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Interestingly, यह भी fail होता है, लेकिन एक different error के साथ।
इस बार Nextflow ने file की contents parse की हैं (yay!) लेकिन इसने प्रत्येक row को एक array के रूप में load किया है, और प्रत्येक array channel में एक element है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-fail.svg"
</figure>

हमें इसे बताना होगा कि प्रत्येक row में केवल first column ले।
तो हम इसे कैसे unpack करें?

हमने पहले channel की contents unpack करने के लिए `flatten()` use किया है, लेकिन यह यहाँ काम नहीं करेगा क्योंकि flatten _everything_ unpack करता है (यदि तुम खुद देखना चाहते हो तो try करो)।

इसके बजाय, हम `map()` नामक एक और operator use करेंगे जो वास्तव में useful है और Nextflow pipelines में बहुत pop up होता है।

### 4.3. Greetings extract करने के लिए `map()` operator use करें

[`map()`](https://www.nextflow.io/docs/latest/reference/operator.html#map) operator एक बहुत handy little tool है जो हमें channel की contents पर सभी प्रकार की mappings करने की अनुमति देता है।

इस case में, हम इसे अपनी data file में प्रत्येक row से उस एक element को extract करने के लिए use करेंगे जो हम चाहते हैं।
Syntax ऐसा दिखता है:

```groovy title="Syntax"
.map { row -> row[0] }
```

इसका मतलब है 'channel में प्रत्येक row के लिए, उसमें contained 0th (first) item लो'।

तो चलो इसे अपने CSV parsing पर apply करते हैं।

#### 4.3.1. Channel पर `map()` apply करें

Workflow block में, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "After map: $csv" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

तुम देखोगे कि हमने confirm करने के लिए एक और `view()` call add किया कि operator वही करता है जो हम expect करते हैं।

#### 4.3.2. Workflow चलाएं

चलो इसे एक बार और run करते हैं:

```bash
nextflow run hello-channels.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [focused_volhard] DSL2 - revision: de435e45be

    executor >  local (3)
    [54/6eebe3] sayHello (3) [100%] 3 of 3 ✔
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    After map: Hello
    After map: Bonjour
    After map: Holà
    ```

इस बार यह बिना error के run होना चाहिए।

`view()` statements के output को देखते हुए, तुम निम्नलिखित देखोगे:

- एक single `Before splitCsv:` statement: उस point पर channel में एक item है, original file path।
- तीन separate `After splitCsv:` statements: प्रत्येक greeting के लिए एक, लेकिन प्रत्येक एक array के भीतर contained है जो file में उस line से correspond करता है।
- तीन separate `After map:` statements: प्रत्येक greeting के लिए एक, जो अब channel में individual elements हैं।

_Note करो कि lines तुम्हारे output में different order में appear हो सकती हैं।_

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-and-map.svg"
</figure>

तुम यह verify करने के लिए output files भी देख सकते हो कि प्रत्येक greeting correctly extract और workflow के माध्यम से process हुई।

हमने पहले जैसा same result achieve किया है, लेकिन अब हमारे पास process करने के लिए greetings के channel में more elements add करने की बहुत अधिक flexibility है एक input file modify करके, बिना कोई code modify किए।
तुम बाद की training में complex inputs handle करने के लिए more sophisticated approaches सीखोगे।

### सीख

तुम जानते हो कि `.fromPath()` channel constructor और operators `splitCsv()` और `map()` का उपयोग करके input values की file पढ़ना और उन्हें appropriately handle करना।

अधिक generally, तुम्हारे पास basic understanding है कि Nextflow कैसे processes को inputs manage करने के लिए **channels** और उनकी contents transform करने के लिए **operators** use करता है।
तुमने यह भी देखा है कि channels कैसे parallel execution को implicitly handle करते हैं।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-parallel.svg"
</figure>

### आगे क्या?

एक big break लो, तुमने इसमें hard work किया!

जब तुम ready हो, तो [**Part 3: Hello Workflow**](./03_hello_workflow.md) पर move करो यह सीखने के लिए कि more steps कैसे add करें और उन्हें एक proper workflow में कैसे connect करें।

---

## Quiz

<quiz>
Nextflow में channel क्या है?
- [ ] एक file path specification
- [ ] एक process definition
- [x] Processes के बीच data pass करने के लिए एक queue-जैसी structure
- [ ] एक configuration setting

और जानें: [1.1. Create an input channel](#11-create-an-input-channel)
</quiz>

<quiz>
यह code क्या output करेगा?

```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

- [ ] `['Hello', 'Bonjour', 'Hola']` (एक single list)
- [x] प्रत्येक element separate line पर: `Hello`, `Bonjour`, `Hola`
- [ ] कुछ नहीं (channels default रूप से print नहीं करते)
- [ ] एक error (invalid syntax)

और जानें: [1.1. Create an input channel](#11-create-an-input-channel)
</quiz>

<quiz>
जब channel में multiple values होती हैं, Nextflow process execution कैसे handle करता है?
- [ ] Process सभी values के साथ एक बार run होता है
- [x] Process channel में प्रत्येक value के लिए एक बार run होता है
- [ ] Process केवल first value के साथ run होता है
- [ ] Process केवल last value के साथ run होता है

और जानें: [2. Modify the workflow to run on multiple input values](#2-modify-the-workflow-to-run-on-multiple-input-values)
</quiz>

<quiz>
`flatten()` operator क्या करता है?
- [ ] Multiple channels को एक में combine करता है
- [ ] Channel elements को sort करता है
- [x] Arrays को individual elements में unpack करता है
- [ ] Duplicate elements remove करता है

और जानें: [3.2.1. Add the `flatten()` operator](#321-add-the-flatten-operator)
</quiz>

<quiz>
`view()` operator का purpose क्या है?
- [ ] Channel contents filter करना
- [ ] Channel elements transform करना
- [x] Channel contents inspect और debug करना
- [ ] Channel contents को file में save करना

और जानें: [1.4. Use `view()` to inspect the channel contents](#14-use-view-to-inspect-the-channel-contents)
</quiz>

<quiz>
`splitCsv()` क्या करता है?
- [ ] Channel contents से CSV file बनाता है
- [ ] String को commas से split करता है
- [x] CSV file को प्रत्येक row represent करने वाले arrays में parse करता है
- [ ] Multiple CSV files merge करता है

और जानें: [4.2. Use the `splitCsv()` operator to parse the file](#42-use-the-splitcsv-operator-to-parse-the-file)
</quiz>

<quiz>
`map()` operator का purpose क्या है?
- [ ] Channel से elements filter करना
- [ ] Multiple channels combine करना
- [x] Channel में प्रत्येक element को transform करना
- [ ] Channel में elements count करना

और जानें: [4.3. Use the `map()` operator to extract the greetings](#43-use-the-map-operator-to-extract-the-greetings)
</quiz>

<quiz>
Multiple inputs process करते समय dynamic output filenames use करना क्यों important है?
- [ ] Performance improve करने के लिए
- [ ] Disk space reduce करने के लिए
- [x] Output files को एक दूसरे को overwrite करने से रोकने के लिए
- [ ] Resume functionality enable करने के लिए

और जानें: [2.2. Ensure the output file names will be unique](#22-ensure-the-output-file-names-will-be-unique)
</quiz>
