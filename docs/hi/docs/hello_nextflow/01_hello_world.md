# भाग 1: Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube channel पर [पूरी playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) देखें।

:green_book: वीडियो transcript [यहाँ](./transcripts/01_hello_world.md) उपलब्ध है।
///
-->

Hello Nextflow training course के इस पहले भाग में, हम एक बहुत ही बुनियादी domain-agnostic Hello World उदाहरण के साथ topic में आसानी से प्रवेश करते हैं, जिसे हम foundational Nextflow logic और components के उपयोग को demonstrate करने के लिए progressively build up करेंगे।

??? info "Hello World उदाहरण क्या है?"

    एक "Hello World!" एक minimalist उदाहरण है जो programming language या software framework की basic syntax और structure को demonstrate करने के लिए है।
    उदाहरण में आमतौर पर output device, जैसे console या terminal, पर "Hello, World!" phrase print करना, या इसे एक file में लिखना शामिल है।

---

## 0. Warmup: Hello World उदाहरण को सीधे चलाएं

चलो इसे एक simple command के साथ demonstrate करते हैं जिसे हम सीधे terminal में चलाते हैं, यह दिखाने के लिए कि यह क्या करता है इससे पहले कि हम इसे Nextflow में wrap करें।

!!! tip "सुझाव"

    याद रखो कि तुम्हें अब `hello-nextflow/` directory के अंदर होना चाहिए जैसा कि [Getting Started](00_orientation.md) page पर describe किया गया है।

### 0.1. Terminal को hello कहलवाएं

अपने terminal में निम्नलिखित command चलाओ।

```bash
echo 'Hello World!'
```

??? success "Command output"

    ```console
    Hello World!
    ```

यह terminal में 'Hello World' text output करता है।

### 0.2. Output को एक file में लिखें

Pipelines चलाने में ज्यादातर files से data पढ़ना और results को अन्य files में लिखना शामिल है, तो चलो command को modify करते हैं ताकि text output को एक file में लिखें ताकि उदाहरण थोड़ा अधिक relevant हो।

```bash
echo 'Hello World!' > output.txt
```

??? success "Command output"

    ```console

    ```

यह terminal में कुछ भी output नहीं करता।

### 0.3. Output खोजें

'Hello World' text अब उस output file में होना चाहिए जो हमने specify की थी, जिसका नाम `output.txt` है।
तुम इसे file explorer में या command line से `cat` utility का उपयोग करके खोल सकते हो, उदाहरण के लिए।

??? abstract "File contents"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

यही वह है जिसे हम अपने पहले Nextflow workflow के साथ replicate करने की कोशिश करेंगे।

### सारांश

तुम अब जानते हो कि terminal में एक simple command कैसे चलाया जाए जो कुछ text output करता है, और optionally, इसे एक file में output कैसे लिखवाया जाए।

### आगे क्या है?

पता लगाओ कि यह Nextflow workflow के रूप में कैसा दिखेगा।

---

## 1. Script को examine करें और चलाएं

हम तुम्हें `hello-world.nf` नामक एक fully functional अगर minimalist workflow script प्रदान करते हैं जो पहले जैसा ही काम करता है ('Hello World!' लिखता है) लेकिन Nextflow के साथ।

शुरू करने के लिए, चलो workflow script खोलते हैं ताकि तुम समझ सको कि यह कैसे structured है।
फिर हम इसे चलाएंगे और इसके outputs खोजेंगे।

### 1.1. Code को examine करें

तुम `hello-world.nf` script अपनी current directory में पाओगे, जो `hello-nextflow` होनी चाहिए। इसे editor pane में खोलो।

??? full-code "पूरी code file"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' को एक फ़ाइल में प्रिंट करने के लिए echo का उपयोग करें
    */
    process sayHello {

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }

    workflow {

        main:
        // एक अभिवादन emit करें
        sayHello()
    }
    ```

एक Nextflow workflow script में आमतौर पर एक या अधिक [**process**](https://nextflow.io/docs/latest/process.html) definitions और [**workflow**](https://nextflow.io/docs/latest/workflow.html) itself शामिल होते हैं, साथ ही कुछ optional blocks (यहाँ मौजूद नहीं) जिन्हें हम बाद में introduce करेंगे।

प्रत्येक **process** describe करता है कि pipeline में corresponding step को क्या operation(s) accomplish करने चाहिए, जबकि **workflow** dataflow logic describe करता है जो विभिन्न steps को connect करता है।

हम पहले **process** block पर एक closer look लेंगे, फिर हम **workflow** block को देखेंगे।

#### 1.1.1. `process` definition

Code का पहला block एक **process** describe करता है।

Process definition keyword `process` से शुरू होती है, उसके बाद process name और अंत में curly braces द्वारा delimited process body।
Process body में एक script block होना चाहिए जो चलाने के लिए command specify करता है, जो कुछ भी हो सकता है जो तुम command line terminal में चला सकते हो।

```groovy title="hello-world.nf" linenums="3"
/*
* 'Hello World!' को एक फ़ाइल में प्रिंट करने के लिए echo का उपयोग करें
*/
process sayHello {

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

यहाँ हमारे पास `sayHello` नामक एक **process** है जो अपना **output** `output.txt` नामक एक file में लिखता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

यह एक बहुत ही minimal process definition है जिसमें बस एक `output` definition और execute करने के लिए `script` है।

`output` definition में `path` qualifier शामिल है, जो Nextflow को बताता है कि इसे एक path के रूप में handle किया जाना चाहिए (जिसमें directory paths और files दोनों शामिल हैं)।
एक अन्य common qualifier `val` है।

महत्वपूर्ण रूप से, output definition यह _determine_ नहीं करती कि क्या output बनाया जाएगा।
यह बस _declare_ करती है कि expected output क्या है, ताकि Nextflow execution complete होने के बाद इसे खोज सके।
यह verify करने के लिए आवश्यक है कि command successfully execute हुई और यदि आवश्यक हो तो output को downstream processes को pass करने के लिए। जो output produce होता है जो output block में declare किए गए से match नहीं करता वह downstream processes को pass नहीं किया जाएगा।

!!! warning "चेतावनी"

    यह उदाहरण brittle है क्योंकि हमने output filename को दो अलग-अलग जगहों (script और output blocks) में hardcode किया है।
    यदि हम एक को बदलते हैं लेकिन दूसरे को नहीं, तो script break हो जाएगी।
    बाद में, तुम इस problem को mitigate करने के लिए variables का उपयोग करने के तरीके सीखोगे।

एक real-world pipeline में, एक process में आमतौर पर additional blocks जैसे directives और inputs होते हैं, जिन्हें हम थोड़ी देर में introduce करेंगे।

#### 1.1.2. `workflow` definition

Code का दूसरा block **workflow** itself describe करता है।
Workflow definition keyword `workflow` से शुरू होती है, उसके बाद एक optional name, फिर curly braces द्वारा delimited workflow body।

यहाँ हमारे पास एक **workflow** है जिसमें एक `main:` block है (जो कहता है 'यह workflow का main body है') जिसमें `sayHello` process का call है।

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // एक अभिवादन emit करें
    sayHello()
}
```

यह एक बहुत ही minimal **workflow** definition है।
एक real-world pipeline में, workflow में आमतौर पर **channels** द्वारा connected **processes** के multiple calls होते हैं, और processes एक या अधिक variable **input(s)** expect करते हैं।

तुम इस training module में बाद में variable inputs जोड़ना सीखोगे; और तुम इस course के Part 3 में more processes जोड़ना और उन्हें channels द्वारा connect करना सीखोगे।

!!! tip "सुझाव"

    तकनीकी रूप से `main:` line इस तरह के simple workflows के लिए आवश्यक नहीं है, इसलिए तुम ऐसे workflows देख सकते हो जिनमें यह नहीं है।
    लेकिन workflow-level outputs का लाभ उठाने के लिए हमें इसकी आवश्यकता होगी, इसलिए हम इसे शुरू से ही शामिल कर सकते हैं।

### 1.2. Workflow चलाएं

Code देखना इसे चलाने जितना fun नहीं है, तो चलो इसे practice में try करते हैं।

#### 1.2.1. Workflow launch करें और execution monitor करें

Terminal में, निम्नलिखित command चलाओ:

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console hl_lines="7"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [65/7be2fa] sayHello | 1 of 1 ✔
    ```

यदि तुम्हारा console output ऐसा कुछ दिखता है, तो बधाई हो, तुमने अभी अपना पहला Nextflow workflow चलाया!

यहाँ सबसे महत्वपूर्ण output अंतिम line है, जो ऊपर के output में highlighted है:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

यह हमें बताता है कि `sayHello` process successfully एक बार execute हुआ (`1 of 1 ✔`)।

महत्वपूर्ण रूप से, यह line तुम्हें यह भी बताती है कि `sayHello` process call का output कहाँ मिलेगा।
चलो अब उसे देखते हैं।

#### 1.2.2. `work` directory में output और logs खोजें

जब तुम किसी दी गई directory में पहली बार Nextflow चलाते हो, तो यह `work` नामक एक directory बनाता है जहाँ यह execution के दौरान generate की गई सभी files (और कोई भी symlinks) लिखेगा।

`work` directory के भीतर, Nextflow outputs और logs को प्रति process call organize करता है।
प्रत्येक process call के लिए, Nextflow एक nested subdirectory बनाता है, जिसे unique बनाने के लिए hash के साथ named किया जाता है, जहाँ यह सभी आवश्यक inputs को stage करेगा (default रूप से symlinks का उपयोग करके), helper files लिखेगा, और logs और process के किसी भी outputs को लिखेगा।

उस subdirectory का path console output में square brackets में truncated form में दिखाया जाता है।
ऊपर दिखाए गए run के लिए हमें जो मिला उसे देखते हुए, sayHello process के लिए console log line `[65/7be2fa]` से शुरू होती है। यह निम्नलिखित directory path से correspond करता है: `work/65/7be2fad5e71e5f49998f795677fd68`

चलो देखते हैं कि वहाँ क्या है।

??? abstract "Directory contents"

    ```console
    work
    └── 65
        └── 7be2fad5e71e5f49998f795677fd68
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "वही चीज़ नहीं दिख रही?"

    तुम्हारे system पर exact subdirectory names अलग होंगे।

    यदि तुम VSCode file explorer में task subdirectory की contents browse करते हो, तो तुम सभी files तुरंत देखोगे।
    हालाँकि, log files terminal में invisible होने के लिए set हैं, इसलिए यदि तुम उन्हें देखने के लिए `ls` या `tree` का उपयोग करना चाहते हो, तो तुम्हें invisible files display करने के लिए relevant option set करना होगा।

    ```bash
    tree -a work
    ```

पहली चीज़ जो तुम देखना चाहते हो वह है workflow का actual output, यानी `sayHello` process द्वारा produce की गई `output.txt` file।
इसे खोलो और तुम `Hello World!` greeting पाओगे, जो हमारे minimalist workflow का point था।

??? abstract "File contents"

    ```console title="output.txt"
    Hello World!
    ```

यह काम किया!

माना, इतने छोटे result के लिए यह बहुत सारा wrapper code लग सकता है, लेकिन उस सारे wrapper code का value अधिक obvious हो जाएगा जब हम input files पढ़ना और multiple steps को string together करना शुरू करेंगे।

उस ने कहा, चलो उस directory में अन्य files को भी देखते हैं। ये task execution के भाग के रूप में Nextflow द्वारा produce की गई helper और log files हैं।

- **`.command.begin`**: Process call के execution की शुरुआत से संबंधित Metadata
- **`.command.err`**: Process call द्वारा emit किए गए Error messages (`stderr`)
- **`.command.log`**: Process call द्वारा emit किया गया Complete log output
- **`.command.out`**: Process call द्वारा Regular output (`stdout`)
- **`.command.run`**: Process call को execute करने के लिए Nextflow द्वारा चलाई गई Full script
- **`.command.sh`**: वह command जो वास्तव में process call द्वारा चलाई गई
- **`.exitcode`**: Command से resulting exit code

`.command.sh` file विशेष रूप से useful है क्योंकि यह तुम्हें बताती है कि Nextflow ने main command क्या execute की, सभी bookkeeping और task/environment setup को शामिल नहीं करते हुए।

??? abstract "File contents"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

यह उससे match करता है जो हमने पहले manually चलाया था।

इस case में यह बहुत straightforward है क्योंकि process command hardcoded थी, लेकिन बाद में course में तुम process commands देखोगे जिनमें variables का कुछ interpolation शामिल है।
यह विशेष रूप से valuable बनाता है कि तुम exactly देख सको कि Nextflow ने code को कैसे interpret किया और जब तुम failed run को troubleshoot कर रहे हो तो क्या command produce हुई।

### 1.3. Workflow फिर से चलाएं

Workflow को कुछ बार फिर से चलाने का try करो, फिर `work/` के तहत task directories देखो।

??? abstract "Directory contents"

    ```console
    work
    ├── 0f
    │   └── 52b7e07b0e274a80843fca48ed21b8
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 65
        └── 7be2fad5e71e5f49998f795677fd68
    │   │   ├── .command.begin
    │   │   ├── .command.err
    │   │   ├── .command.log
    │   │   ├── .command.out
    │   │   ├── .command.run
    │   │   ├── .command.sh
    │   │   ├── .exitcode
    │   │   └── output.txt
    │   └── e029f2e75305874a9ab263d21ebc2c
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 6c
    │   └── d4fd787e0b01b3c82e85696c297500
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── e8
        └── ab99fad46ade52905ec973ff39bb80
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

तुम देखोगे कि प्रत्येक run के लिए output और log files के complete set के साथ एक new subdirectory बनाई गई है।
यह तुम्हें दिखाता है कि same workflow को कई बार चलाने से previous runs के results overwrite नहीं होंगे।

### सारांश

तुम जानते हो कि एक simple Nextflow script को कैसे decipher करना है, इसे कैसे चलाना है और work directory में output और relevant log files कैसे खोजनी हैं।

### आगे क्या है?

सीखो कि workflow outputs को एक अधिक convenient location पर कैसे publish करें।

---

## 2. Outputs publish करें

जैसा कि तुमने अभी सीखा, हमारे pipeline द्वारा produce किया गया output कई layers deep एक working directory में buried है।
यह जानबूझकर किया गया है; Nextflow इस directory का control में है और हमें इसके साथ interact नहीं करना चाहिए।
हालाँकि, यह उन outputs को retrieve करने में असुविधाजनक बनाता है जिनकी हमें care है।

सौभाग्य से, Nextflow [workflow output definitions](https://nextflow.io/docs/latest/workflow.html#workflow-outputs) का उपयोग करके outputs को एक designated directory में publish करने का एक तरीका प्रदान करता है।

### 2.1. Basic usage

इसमें code के दो नए pieces शामिल होंगे:

1. `workflow` body के अंदर एक `publish:` block, जो process outputs declare करता है।
2. Script में एक `output` block जो output options जैसे mode और location specify करता है।

#### 2.1.1. `sayHello` process का output declare करें

हमें workflow body में एक `publish:` block जोड़ना होगा (उसी प्रकार का code element जैसा `main:` block) और `sayHello()` process का output list करना होगा।

Workflow script file `hello-world.nf` में, code की निम्नलिखित lines जोड़ो:

=== "After"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // एक अभिवादन emit करें
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // एक अभिवादन emit करें
        sayHello()
    }
    ```

तुम देखोगे कि हम process के output को simply `sayHello().out` करके refer कर सकते हैं, और इसे एक arbitrary name, `first_output` assign कर सकते हैं।

#### 2.1.2. Script में एक `output:` block जोड़ें

अब हमें बस वह `output:` block जोड़ना है जहाँ output directory path specify किया जाएगा। Note करो कि यह new block script के भीतर `workflow` block के **बाहर** और **नीचे** बैठता है।

Workflow script file `hello-world.nf` में, code की निम्नलिखित lines जोड़ो:

=== "After"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // एक अभिवादन emit करें
        sayHello()

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '.'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // एक अभिवादन emit करें
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

हम इसका उपयोग `workflow` block में declare किए गए किसी भी process outputs को specific paths assign करने के लिए कर सकते हैं।
बाद में, तुम sophisticated output directory structures generate करने के तरीके सीखोगे, लेकिन अभी के लिए, हम simplicity के लिए एक minimal path hardcode कर रहे हैं।

#### 2.1.3. Workflow चलाएं

अब modified workflow script चलाओ:

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

    executor >  local (1)
    [9f/48ef97] sayHello | 1 of 1 ✔
    ```

Terminal output परिचित दिखना चाहिए। बाहरी रूप से, कुछ भी नहीं बदला है।

हालाँकि, अपना file explorer check करो: इस बार, Nextflow ने `results/` नामक एक new directory बनाई है।

??? abstract "Directory contents"

    ```console hl_lines="10-11 22"
    .
    ├── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── results
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/9f/48ef97f110b0dbd83635d7cbe288d2/output.txt
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── work
        ├── 65
        └── 9f
    ```

`results` directory के अंदर, हमें work directory में command द्वारा produce की गई `output.txt` का एक symbolic link मिलता है जो हमने अभी चलाई।

यह हमें work subdirectory में खोदे बिना आसानी से output files retrieve करने की अनुमति देता है।

### 2.2. Custom location सेट करें

Default location होना great है, लेकिन तुम customize करना चाह सकते हो कि results कहाँ save होते हैं और वे कैसे organized होते हैं।

उदाहरण के लिए, तुम अपने outputs को subdirectories में organize करना चाह सकते हो।
ऐसा करने का सबसे simple तरीका प्रति output specific output path assign करना है।

#### 2.2.1. Output path modify करें

एक बार फिर, specific output के लिए publish behavior modify करना वास्तव में straightforward है।
Custom location set करने के लिए, बस `path` को accordingly edit करो:

=== "After"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path '.'
        }
    }
    ```

चूंकि यह individual output के level पर set है, तुम अपनी needs के अनुसार different locations और subdirectories specify कर सकते हो।

#### 2.2.2. Workflow फिर से चलाएं

चलो try करते हैं।

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [8c/79499c] process > sayHello [100%] 1 of 1 ✔
    ```

इस बार result specified subdirectory के तहत लिखा जाता है।

??? abstract "Directory contents"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

तुम देखोगे कि previous execution का result अभी भी वहाँ है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

तुम जितने चाहो उतने levels of nesting का उपयोग कर सकते हो।
Process name या अन्य variables का उपयोग करके results को organize करने के लिए उपयोग की जाने वाली directories को name करना भी संभव है, और top-level output directory का default name बदलना संभव है (जो `-o` CLI flag या config variable `outputDir` द्वारा controlled है)।
हम इन options को बाद में training में cover करेंगे।

### 2.3. Publish mode को copy पर सेट करें

Default रूप से, outputs `work` directory से symbolic links के रूप में publish होते हैं।
इसका मतलब है कि filesystem पर केवल एक single file है।

यह great है जब तुम बहुत large files के साथ deal कर रहे हो, जिनकी तुम multiple copies store नहीं करना चाहते।
हालाँकि, यदि तुम किसी भी समय work directory delete करते हो (हम cleanup operations को shortly cover करेंगे), तो तुम file तक access खो दोगे।
इसलिए तुम्हें किसी भी important files की copies को एक secure place पर save करने की plan होनी चाहिए।

एक easy option उन outputs के लिए publish mode को copy में switch करना है जिनकी तुम care करते हो।

#### 2.3.1. Mode directive जोड़ें

यह bit वास्तव में straightforward है।
बस relevant workflow-level output definition में `mode 'copy'` जोड़ो:

=== "After"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="4"
    output {
        first_output {
            path 'hello_world'
            mode 'copy'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

यह उस specific output के लिए publish mode set करता है।

#### 2.3.2. Workflow फिर से चलाएं

चलो try करते हैं।

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [df/521638] process > sayHello [100%] 1 of 1 ✔
    ```

इस बार, यदि तुम results देखो, तो file एक proper copy है बजाय सिर्फ एक symlink के।

??? abstract "Directory contents"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

चूंकि यह भी individual output के level पर set है, यह तुम्हें publish mode को एक granular way में set करने की अनुमति देता है।
यह विशेष रूप से बाद में handy होगा जब हम multi-step pipelines पर move करेंगे, जहाँ तुम शायद केवल final outputs को copy करना और intermediate outputs को symlinks के रूप में छोड़ना चाहो, उदाहरण के लिए।

जैसा कि पहले noted किया गया, outputs कैसे publish होते हैं इसे control करने के लिए अन्य, अधिक sophisticated options हैं।
हम तुम्हें तुम्हारी Nextflow journey में due time में उनका उपयोग करना दिखाएंगे।

### 2.4. Process-level `publishDir` directives पर note

बहुत हाल तक, outputs publish करने का established तरीका प्रत्येक individual process के level पर `publishDir` directive का उपयोग करके करना था।

जो हमने अभी `sayHello` process के outputs के लिए किया, उसे achieve करने के लिए, हमने इसके बजाय process definition में निम्नलिखित line जोड़ी होती:

```groovy title="hello-world.nf" linenums="6" hl_lines="3"
process sayHello {

    publishDir 'results/hello_world', mode: 'copy'

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

तुम अभी भी इस code pattern को पुराने Nextflow pipelines और process modules में हर जगह पाओगे, इसलिए इसके बारे में aware होना important है।
हालाँकि, हम किसी भी new work में इसका उपयोग करने की recommend नहीं करते क्योंकि यह eventually Nextflow language के future versions में disallow हो जाएगा।

### सारांश

तुम जानते हो कि workflow outputs को एक अधिक convenient location पर कैसे publish करें।

### आगे क्या है?

सीखो कि command-line parameter के माध्यम से variable input कैसे provide करें और default values को effectively कैसे utilize करें।

---

## 3. Command line पर passed variable input का उपयोग करें

अपनी current state में, हमारा workflow process command में hardcoded greeting का उपयोग करता है।
हम एक input variable का उपयोग करके कुछ flexibility add करना चाहते हैं, ताकि हम runtime पर greeting को अधिक आसानी से बदल सकें।

इसके लिए हमें अपनी script में तीन sets of changes करने होंगे:

1. Process को variable input expect करने के लिए change करें
2. User input capture करने के लिए एक command-line parameter set up करें
3. Workflow body में process को input pass करें

चलो ये changes एक-एक करके करते हैं।

### 3.1. `sayHello` process को variable input expect करने के लिए change करें

हमें process definition को edit करना होगा (1) एक input variable accept करने के लिए और (2) उस variable को command line में use करने के लिए।

#### 3.1.1. Process definition में एक input block जोड़ें

पहले, चलो process definition को `greeting` नामक एक input accept करने के लिए adapt करते हैं।

Process block में, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3-4"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
        path 'output.txt'
    ```

`greeting` variable `val` से prefixed है Nextflow को बताने के लिए कि यह एक value है (न कि एक path)।

#### 3.1.2. Input variable use करने के लिए process command edit करें

अब हम original hardcoded value को उस input variable की value से swap करते हैं जो हम receive करने की expect करते हैं।

Process block में, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo '${greeting}' > output.txt
    """
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

`$` symbol और curly braces (`{ }`) Nextflow को बताते हैं कि यह एक variable name है जिसे actual input value से replace किया जाना चाहिए (=interpolated)।

!!! tip "सुझाव"

    Nextflow के previous versions में curly braces (`{ }`) technically optional थे, इसलिए तुम पुराने workflows देख सकते हो जहाँ यह `echo '$greeting' > output.txt` के रूप में लिखा गया है।

अब जबकि `sayHello()` process variable input accept करने के लिए ready है, हमें workflow level पर process call को input value provide करने का एक तरीका चाहिए।

### 3.2. User input capture करने के लिए command-line parameter set up करें

हम simply एक input को directly hardcode कर सकते हैं process call `sayHello('Hello World!')` बनाकर।
हालाँकि, जब हम अपने workflow के साथ real work कर रहे होते हैं, तो हम command line से इसके inputs को control करने में सक्षम होना चाहेंगे, ताकि हम कुछ ऐसा कर सकें:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

सौभाग्य से, Nextflow में [`params`](https://nextflow.io/docs/latest/config.html#params) नामक एक built-in workflow parameter system है, जो CLI parameters को declare और use करना आसान बनाता है।

General syntax है `params.<parameter_name>` declare करना Nextflow को बताने के लिए कि command line पर एक `--<parameter_name>` parameter expect करें।

यहाँ, हम `--input` नामक एक parameter बनाना चाहते हैं, इसलिए हमें workflow में कहीं `params.input` declare करना होगा।
सिद्धांत रूप में हम इसे कहीं भी लिख सकते हैं; लेकिन चूंकि हम इसे `sayHello()` process call को देना चाहते हैं, हम इसे सीधे वहाँ plug कर सकते हैं `sayHello(params.input)` लिखकर।

Workflow block में, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // एक अभिवादन emit करें
    sayHello(params.input)
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // एक अभिवादन emit करें
    sayHello()
    ```

यह Nextflow को बताता है कि `--input` parameter के माध्यम से provide की गई value पर `sayHello` process चलाएं।

Effect में, हमने section की शुरुआत में outlined steps (2) और (3) एक ही बार में accomplish कर लिए हैं।

### 3.3. Workflow command चलाएं

चलो इसे चलाते हैं!

```bash
nextflow run hello-world.nf --input 'Bonjour le monde!'
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elated_lavoisier] DSL2 - revision: 7c031b42ea

    executor >  local (1)
    [4b/654319] sayHello | 1 of 1 ✔
    ```

यदि तुमने ये सभी edits correctly किए, तो तुम्हें एक और successful execution मिलनी चाहिए।

यह check करने के लिए output file खोलना sure करो कि तुम्हारे पास अब greeting का new version है।

??? abstract "File contents"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Et voilà!

Note करो कि new execution ने `results` directory में publish की गई output file को overwrite कर दिया है।
हालाँकि, previous runs के results अभी भी `work` के तहत task directories में preserved हैं।

!!! tip "सुझाव"

    तुम Nextflow-level parameters को pipeline-level parameters से readily distinguish कर सकते हो।

    - Pipeline पर apply होने वाले Parameters हमेशा double hyphen (`--`) लेते हैं।
    - Nextflow setting को modify करने वाले Parameters, _जैसे_ `-resume` feature जो हमने पहले use किया, single hyphen (`-`) लेते हैं।

### 3.4. Command line parameters के लिए default values use करें

ठीक है, यह convenient था, लेकिन कई cases में, दिए गए parameter के लिए default value supply करना sense बनाता है ताकि तुम्हें हर run के लिए इसे specify न करना पड़े।

#### 3.4.1. CLI parameter के लिए default value set करें

चलो workflow definition से पहले declare करके `input` parameter को एक default value देते हैं।

```groovy title="hello-world.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String = 'Holà mundo!'
}
```

जैसा कि तुम देखते हो, हम specify कर सकते हैं कि workflow किस type के input की expect करता है (Nextflow 25.10.2 और बाद में)।
Syntax है `name: Type = default_value`।
Supported types में `String`, `Integer`, `Float`, `Boolean`, और `Path` शामिल हैं।

!!! info "जानकारी"

    पुराने workflows में, तुम देख सकते हो कि पूरा `params` block सिर्फ `input = 'Holà mundo!'` के रूप में लिखा गया है।

जैसे-जैसे तुम अपने pipeline में more parameters add करते हो, तुम्हें उन सभी को इस block में add करना चाहिए, चाहे तुम्हें उन्हें default value देने की need हो या नहीं।
यह एक glance में सभी configurable parameters खोजना आसान बना देगा।

#### 3.4.2. Parameter specify किए बिना workflow फिर से चलाएं

अब जबकि तुम्हारे पास एक default value set है, तुम command line में value specify किए बिना workflow को फिर से चला सकते हो।

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

    executor >  local (1)
    [72/394147] sayHello | 1 of 1 ✔
    ```

Output पहले जैसी same place में होगा, लेकिन contents new text के साथ updated होनी चाहिए।

??? abstract "File contents"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow ने output बनाने के लिए greeting parameter की default value use की।

#### 3.4.3. Default value override करें

यदि तुम command line पर parameter provide करते हो, तो CLI value default value को override करेगी।

Try करो:

```bash
nextflow run hello-world.nf --input 'Konnichiwa!'
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

    executor >  local (1)
    [6f/a12a91] sayHello | 1 of 1 ✔
    ```

एक बार फिर, तुम्हें अपनी results directory में corresponding updated output मिलनी चाहिए।

??? abstract "File contents"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note "नोट"

    Nextflow में, ऐसी कई places हैं जहाँ तुम parameters के लिए values specify कर सकते हो।
    यदि same parameter multiple places में different values पर set है, तो Nextflow [यहाँ](https://www.nextflow.io/docs/latest/config.html) described precedence के order के आधार पर determine करेगा कि कौन सी value use करनी है।

    हम इसे Part 6 (Configuration) में अधिक detail में cover करेंगे।

### सारांश

तुम जानते हो कि command-line parameter के माध्यम से runtime पर provide किए गए simple variable input का उपयोग कैसे करें, साथ ही default values को set up, use और override कैसे करें।

### आगे क्या है?

सीखो कि executions को अधिक conveniently कैसे manage करें।

---

## 4. Workflow executions manage करें

Workflows launch करना और outputs retrieve करना जानना great है, लेकिन तुम जल्दी ही पाओगे कि workflow management के कुछ अन्य aspects हैं जो तुम्हारी life आसान बना देंगे, खासकर यदि तुम अपने खुद के workflows develop कर रहे हो।

यहाँ हम तुम्हें दिखाते हैं कि जब तुम्हें same workflow re-launch करना हो तो [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) feature कैसे use करें, [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) के साथ past executions का log कैसे inspect करें, और [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean) के साथ older work directories कैसे delete करें।

### 4.1. `-resume` के साथ workflow re-launch करें

कभी-कभी, तुम एक pipeline को re-run करना चाहोगे जो तुम पहले launch कर चुके हो बिना उन steps को redo किए जो पहले से successfully complete हो चुके हैं।

Nextflow में [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) नामक एक option है जो तुम्हें ऐसा करने की अनुमति देता है।
Specifically, इस mode में, कोई भी processes जो पहले से exact same code, settings और inputs के साथ run हो चुके हैं, skip हो जाएंगे।
इसका मतलब है Nextflow केवल वे processes run करेगा जो तुमने last run के बाद से add या modify किए हैं, या जिन्हें तुम new settings या inputs provide कर रहे हो।

ऐसा करने के दो key advantages हैं:

- यदि तुम अपने pipeline develop करने के middle में हो, तो तुम अधिक rapidly iterate कर सकते हो क्योंकि तुम्हें अपने changes test करने के लिए केवल वह process(es) run करने होंगे जिन पर तुम actively काम कर रहे हो।
- यदि तुम production में pipeline run कर रहे हो और कुछ wrong हो जाता है, तो कई cases में तुम issue fix कर सकते हो और pipeline relaunch कर सकते हो, और यह failure के point से running resume करेगा, जो तुम्हें बहुत time और compute बचा सकता है।

इसे use करने के लिए, simply अपने command में `-resume` add करो और इसे run करो:

```bash
nextflow run hello-world.nf -resume
```

??? success "Command output"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

Console output परिचित दिखना चाहिए, लेकिन एक चीज़ है जो पहले की तुलना में थोड़ी different है।

Process status line (line 5) में add हुए `cached:` bit को देखो, जिसका मतलब है कि Nextflow ने recognize किया है कि यह पहले से यह work कर चुका है और simply previous successful run के result को re-use किया है।

तुम यह भी देख सकते हो कि work subdirectory hash previous run जैसा ही है।
Nextflow literally तुम्हें previous execution की ओर point कर रहा है और कह रहा है "मैंने पहले से वह वहाँ कर लिया था।"

!!! tip "सुझाव"

    जब तुम `resume` के साथ pipeline re-run करते हो, Nextflow work directory के बाहर publish की गई किसी भी files को overwrite नहीं करता जो पहले successfully run हुई executions द्वारा publish की गई थीं।

### 4.2. Past executions का log inspect करें

चाहे तुम new pipeline develop कर रहे हो या production में pipelines run कर रहे हो, किसी point पर तुम्हें शायद past runs के बारे में information देखनी होगी।
यहाँ बताया गया है कि ऐसा कैसे करें।

जब भी तुम nextflow workflow launch करते हो, current working directory में `.nextflow` नामक एक hidden directory के तहत `history` नामक एक log file में एक line लिखी जाती है।

??? abstract "File contents"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

यह file तुम्हें current working directory के भीतर से launch की गई हर Nextflow run के लिए timestamp, run name, status, revision ID, session ID और full command line देती है।

इस information को access करने का एक अधिक convenient तरीका `nextflow log` command use करना है।

```bash
nextflow log
```

??? success "Command output"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

यह log file की contents को terminal में output करेगा, एक header line के साथ augmented।

तुम notice करोगे कि session ID तब बदलती है जब तुम एक new `nextflow run` command run करते हो, EXCEPT यदि तुम `-resume` option use कर रहे हो।
उस case में, session ID same रहती है।

Nextflow session ID का उपयोग `.nextflow` के तहत स्थित `cache` directory के तहत run caching information को group करने के लिए करता है।

### 4.3. Older work directories delete करें

Development process के दौरान, तुम typically अपनी draft pipeline को बड़ी संख्या में बार run करोगे, जो कई subdirectories में कई files के accumulation का कारण बन सकता है।

सौभाग्य से Nextflow में एक helpful `clean` subcommand शामिल है जो past runs की work subdirectories को automatically delete कर सकता है जिनकी तुम्हें अब care नहीं है।

#### 4.3.1. Deletion criteria determine करें

यह determine करने के लिए कई [options](https://nextflow.io/docs/latest/reference/cli.html#clean) हैं कि क्या delete करना है।

यहाँ हम तुम्हें एक example दिखाते हैं जो given run से पहले के runs की सभी subdirectories delete करता है, इसके run name का उपयोग करके specified।

Most recent successful run देखो जहाँ तुमने `-resume` use नहीं किया; हमारे case में run name `golden_cantor` था।

Run name machine-generated two-part string है जो `Launching (...)` console output line में square brackets में दिखाया जाता है।
तुम Nextflow log का उपयोग करके इसके timestamp और/या command line के आधार पर run look up भी कर सकते हो।

#### 4.3.2. Dry run करें

पहले हम dry run flag `-n` use करते हैं यह check करने के लिए कि command दिए जाने पर क्या delete होगा:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Command output"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

तुम्हारे output में different task directory names होंगे और lines की different number हो सकती है, लेकिन यह example के समान दिखना चाहिए।

यदि तुम कोई lines output नहीं देखते, तो या तो तुमने valid run name provide नहीं किया या delete करने के लिए कोई past runs नहीं हैं। Example command में `golden_cantor` को तुम्हारे log में जो भी corresponding latest run name है उसमें change करना sure करो।

#### 4.3.3. Deletion के साथ proceed करें

यदि output expected दिखता है और तुम deletion के साथ proceed करना चाहते हो, तो `-n` के बजाय `-f` flag के साथ command re-run करो:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Command output"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Output पहले जैसा similar होना चाहिए, लेकिन अब 'Would remove' के बजाय 'Removed' कह रहा है।
Note करो कि यह two-character subdirectories (जैसे ऊपर `a3/`) को remove नहीं करता लेकिन यह उनकी contents को empty कर देता है।

!!! warning "चेतावनी"

    Past runs की work subdirectories delete करना उन्हें Nextflow के cache से remove करता है और उन directories में stored किसी भी outputs को delete करता है।
    इसका मतलब है कि यह corresponding processes को re-run किए बिना execution resume करने की Nextflow की ability break करता है।

    तुम किसी भी outputs को save करने के लिए responsible हो जिनकी तुम care करते हो या जिन पर rely करने की plan है! यही main reason है कि हम `publish` directive के लिए `symlink` mode के बजाय `copy` mode use करना prefer करते हैं।

### सारांश

तुम जानते हो कि outputs को एक specific directory में कैसे publish करें, पहले से identical way में run हो चुके steps को repeat किए बिना pipeline कैसे relaunch करें, और old work directories को clean up करने के लिए `nextflow clean` command कैसे use करें।

अधिक generally, तुम जानते हो कि एक simple Nextflow workflow को कैसे interpret करें, इसके execution को manage करें, और outputs retrieve करें।

### आगे क्या है?

थोड़ा break लो, तुमने इसे earn किया है!

जब तुम ready हो, तो [**Part 2: Hello Channels**](./02_hello_channels.md) पर move करो यह सीखने के लिए कि अपने workflow में inputs feed करने के लिए channels कैसे use करें, जो तुम्हें Nextflow के built-in dataflow parallelism और अन्य powerful features का लाभ उठाने की अनुमति देगा।

---

## Quiz

<quiz>
Nextflow process के minimum required components क्या हैं?
- [ ] केवल Input और output blocks
- [x] Output और script blocks
- [ ] Input, output, और script blocks
- [ ] केवल एक script block

और जानें: [1.1.1. The `process` definition](#111-the-process-definition)
</quiz>

<quiz>
Process में output block का purpose क्या है?
- [ ] Console पर results print करना
- [ ] Files को work directory में save करना
- [x] Process से expected outputs declare करना
- [ ] Environment variables define करना

और जानें: [1.1.1. The `process` definition](#111-the-process-definition)
</quiz>

<quiz>
Nextflow workflow run करने के लिए कौन सा command use होता है?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
Task की work directory देखते हुए, कौन सी file वह actual command contain करती है जो execute हुई थी?

```
work/a3/7be2fa.../
├── .command.begin
├── .command.err
├── .command.log
├── .command.out
├── .command.run
├── .command.sh
├── .exitcode
└── output.txt
```

- [ ] `.command.run`
- [x] `.command.sh`
- [ ] `.command.log`
- [ ] `.command.out`

और जानें: [1.2.2. Find the output and logs in the `work` directory](#122-find-the-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
`-resume` flag क्या करता है?
- [ ] Workflow को beginning से restart करता है
- [ ] Workflow को pause करता है
- [x] उन processes को skip करता है जो पहले से successfully complete हो चुके हैं
- [ ] Workflow का backup बनाता है

और जानें: [4.1. Re-launch a workflow with `-resume`](#41-re-launch-a-workflow-with--resume)
</quiz>

<quiz>
Workflow outputs publish करने के लिए default mode क्या है?
- [ ] Files को output directory में copy करना
- [x] Output directory में symbolic links बनाना
- [ ] Files को output directory में move करना
- [ ] Output directory में files compress करना

और जानें: [2.3. Set the publish mode to copy](#23-set-the-publish-mode-to-copy)
</quiz>

<quiz>
Command line से Nextflow workflow को parameter value कैसे pass करते हो?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

और जानें: [3.2. Set up a command-line parameter to capture user input](#32-set-up-a-command-line-parameter-to-capture-user-input)
</quiz>

<quiz>
Nextflow script block के अंदर variable कैसे reference करते हो?
- [ ] `%variable%` syntax use करें
- [x] `#!groovy ${variable}` syntax use करें
- [ ] `{{variable}}` syntax use करें
- [ ] `[variable]` syntax use करें
</quiz>
