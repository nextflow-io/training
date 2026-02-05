# भाग 3: Hello Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube channel पर [पूरी playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) देखें।

:green_book: वीडियो transcript [यहाँ](./transcripts/03_hello_workflow.md) उपलब्ध है।
///
-->

अधिकांश real-world workflows में एक से अधिक step शामिल होते हैं।
इस training module में, तुम सीखोगे कि multi-step workflow में processes को एक साथ कैसे connect करें।

यह तुम्हें निम्नलिखित achieve करने का Nextflow तरीका सिखाएगा:

1. एक process से दूसरे में data flow करवाना
2. Multiple process calls से outputs को एक single process call में collect करना
3. एक process को additional parameters pass करना
4. एक process से आने वाले multiple outputs handle करना

Demonstrate करने के लिए, हम Parts 1 और 2 के domain-agnostic Hello World example पर build करना जारी रखेंगे।
इस बार, हम अपने workflow में निम्नलिखित changes करेंगे जो बेहतर reflect करते हैं कि लोग actual workflows कैसे बनाते हैं:

1. Greeting को uppercase में convert करने वाला एक second step add करें।
2. सभी transformed greetings collect करने और उन्हें एक single file में लिखने वाला एक third step add करें।
3. Final output file को name करने के लिए एक parameter add करें और उसे collection step को secondary input के रूप में pass करें।
4. Collection step को process किए गए greetings के बारे में एक simple statistic भी report करवाएं।

??? info "इस section से कैसे शुरू करें"

    Course का यह section मानता है कि तुमने [Hello Nextflow](./index.md) course के Parts 1-2 complete कर लिए हैं, लेकिन यदि तुम उन sections में covered basics से comfortable हो, तो तुम बिना कुछ special किए यहाँ से शुरू कर सकते हो।

---

## 0. वार्मअप: `hello-workflow.nf` चलाएं

हम starting point के रूप में workflow script `hello-workflow.nf` use करेंगे।
यह इस training course के Part 2 में काम करके produce की गई script के equivalent है, सिवाय इसके कि हमने `view()` statements remove कर दिए हैं और output destination बदल दी है:

```groovy title="hello-workflow.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
```

यह diagram workflow के current operation को summarize करता है।
यह familiar लगना चाहिए, सिवाय इसके कि अब हम explicitly दिखा रहे हैं कि process के outputs channel में package किए जाते हैं, ठीक वैसे ही जैसे inputs थे।
हम एक मिनट में उस output channel का अच्छा उपयोग करने वाले हैं।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-workflow-channels.svg"
</figure>

यह sure करने के लिए कि सब कुछ काम कर रहा है, कोई भी changes करने से पहले script को एक बार run करो:

```bash
nextflow run hello-workflow.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [admiring_lamarr] DSL2 - revision: 4d4053520d

    executor >  local (3)
    [b1/5826b5] process > sayHello (2) [100%] 3 of 3 ✔
    ```

पहले की तरह, तुम `output` block में specified location पर output files पाओगे।
इस chapter के लिए, यह `results/hello_workflow/` के तहत है।

??? abstract "Directory contents"

    ```console
    results/hello_workflow
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    └── Holà-output.txt
    ```

यदि यह तुम्हारे लिए काम किया, तो तुम multi-step workflow assemble करना सीखने के लिए ready हो।

---

## 1. Workflow में एक second step add करें

हम प्रत्येक greeting को uppercase में convert करने के लिए एक step add करेंगे।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

इसके लिए, हमें तीन चीजें करनी होंगी:

- Uppercase conversion करने के लिए जो command use करेंगे उसे define करें।
- Uppercasing command wrap करने वाला एक new process लिखें।
- Workflow block में new process call करें और इसे `sayHello()` process के output को input के रूप में लेने के लिए set up करें।

### 1.1. Uppercasing command define करें और terminal में test करें

Greetings को uppercase में conversion करने के लिए, हम `tr` नामक एक classic UNIX tool use करेंगे 'text replacement' के लिए, निम्नलिखित syntax के साथ:

```bash title="Syntax"
tr '[a-z]' '[A-Z]'
```

यह एक बहुत naive text replacement one-liner है जो accented letters के लिए account नहीं करती, इसलिए उदाहरण के लिए 'Holà' 'HOLà' बन जाएगा, लेकिन यह Nextflow concepts demonstrate करने के लिए काफी अच्छा job करेगी और यही matter करता है।

इसे test करने के लिए, हम `echo 'Hello World'` command run कर सकते हैं और इसका output `tr` command को pipe कर सकते हैं:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

Output `UPPER-output.txt` नामक एक text file है जिसमें `Hello World` string का uppercase version है।

??? abstract "File contents"

    ```console title="UPPER-output.txt"
    HELLO WORLD
    ```

यही basically हम अपने workflow के साथ करने की कोशिश करेंगे।

### 1.2. Uppercasing step को Nextflow process के रूप में लिखें

हम अपने new process को पहले वाले पर model कर सकते हैं, क्योंकि हम सभी same components use करना चाहते हैं।

Workflow script में, पहले वाले के ठीक नीचे निम्नलिखित process definition add करो:

```groovy title="hello-workflow.nf" linenums="20"
/*
 * अभिवादन को uppercase में बदलने के लिए text replacement tool का उपयोग करें
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

इसमें, हम input filename के आधार पर second output filename compose करते हैं, जैसा हमने originally first process के output के लिए किया था।

### 1.3. Workflow block में new process का call add करें

अब हमें Nextflow को बताना होगा कि actually हमने जो process define किया उसे call करे।

Workflow block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="44" hl_lines="10-11"
    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)
        // अभिवादन को uppercase में बदलें
        convertToUpper()

        publish:
        first_output = sayHello.out
    }
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="44"
    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

यह अभी functional नहीं है क्योंकि हमने specify नहीं किया कि `convertToUpper()` process को क्या input होना चाहिए।

### 1.4. First process का output second process को pass करें

अब हमें `sayHello()` process का output `convertToUpper()` process में flow करवाना होगा।

Conveniently, Nextflow automatically process का output channel में package करता है, जैसा warmup section के diagram में दिखाया गया है।
हम किसी process के output channel को `<process>.out` के रूप में refer कर सकते हैं।

तो `sayHello` process का output `sayHello.out` नामक channel है, जिसे हम सीधे `convertToUpper()` के call में plug कर सकते हैं।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep-connector.svg"
</figure>

Workflow block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // अभिवादन को uppercase में बदलें
        convertToUpper()
    ```

इस तरह के simple case (एक output से एक input) के लिए, दो processes connect करने के लिए बस इतना करना होगा!

### 1.5. Workflow output publishing set up करें

Finally, चलो workflow outputs update करते हैं ताकि second process के results भी publish हों।

#### 1.5.1. `workflow` block का `publish:` section update करें

`workflow` block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="56" hl_lines="3"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
    }
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="56"
        publish:
        first_output = sayHello.out
    }
    ```

Logic पहले जैसा ही है।

#### 1.5.2. `output` block update करें

`output` block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="61" hl_lines="6-9"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="61"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

एक बार फिर, logic पहले जैसा ही है।

यह तुम्हें दिखाता है कि तुम बहुत granular level पर output settings control कर सकते हो, प्रत्येक individual output के लिए।
बेझिझक processes में से एक के लिए paths या publish mode बदलकर देखो कि क्या होता है।

बेशक, इसका मतलब है कि हम यहाँ कुछ information repeat कर रहे हैं, जो inconvenient हो सकती है यदि हम सभी outputs के लिए location को same तरीके से update करना चाहें।
Course में बाद में, तुम सीखोगे कि multiple outputs के लिए इन settings को structured तरीके से कैसे configure करें।

### 1.6. `-resume` के साथ workflow चलाएं

चलो `-resume` flag का उपयोग करके इसे test करते हैं, क्योंकि हम पहले से workflow का first step successfully run कर चुके हैं।

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [high_cantor] DSL2 - revision: d746983511

    executor >  local (3)
    [ab/816321] process > sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [e0/ecf81b] process > convertToUpper (3) [100%] 3 of 3 ✔
    ```

Console output में अब एक extra line है जो हमने अभी add किए new process से correspond करती है।

तुम `results/hello_workflow` directory में outputs पाओगे जैसा `output` block में set है।

??? abstract "Directory contents"

    ```console
    results/hello_workflow/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

यह convenient है! लेकिन यह अभी भी second process की calls में से एक की work directory के अंदर देखना worth है।

??? abstract "Directory contents"

    ```console
    work/e0/ecf81b4cacc648b9b994218d5b29d7/
    ├── Holà-output.txt -> /workspaces/training/hello-nextflow/work/ab/81632178cd37e9e815959278808819/Holà-output.txt
    └── UPPER-Holà-output.txt
    ```

Notice करो कि दो `*-output` files हैं: first process का output साथ ही second का output।

First process का output वहाँ है क्योंकि Nextflow ने execution के लिए जरूरी सब कुछ same subdirectory के अंदर रखने के लिए इसे वहाँ **staged** किया।

हालांकि, यह actually first process call की subdirectory में original file को point करने वाला symbolic link है।
By default, जब एक single machine पर run कर रहे हों जैसा हम यहाँ कर रहे हैं, तो Nextflow input और intermediate files stage करने के लिए copies के बजाय symbolic links use करता है।

अब, आगे बढ़ने से पहले, सोचो कि हमने केवल `sayHello` के output को `convertToUpper` के input से connect किया और दोनों processes series में run हो सकीं।
Nextflow ने individual input और output files handle करने और उन्हें दो commands के बीच pass करने का hard work हमारे लिए किया।

यह एक कारण है कि Nextflow channels इतने powerful हैं: वे workflow steps को एक साथ connect करने में involved busywork का ख्याल रखते हैं।

### सारांश

तुम जानते हो कि एक step का output अगले step को input के रूप में provide करके processes को chain कैसे करें।

### आगे क्या है?

सीखो कि batched process calls से outputs collect करके एक single process में कैसे feed करें।

---

## 2. सभी greetings collect करने के लिए third step add करें

जब हम यहाँ कर रहे हैं उस तरह channel में प्रत्येक element पर transformation apply करने के लिए process use करते हैं, तो कभी-कभी हम उस process के output channel से elements collect करना और उन्हें किसी प्रकार का analysis या summation perform करने वाले दूसरे process में feed करना चाहते हैं।

Demonstrate करने के लिए, हम अपनी pipeline में एक new step add करेंगे जो `convertToUpper` process द्वारा produce किए गए सभी uppercase greetings collect करता है और उन्हें एक single file में लिखता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

Surprise spoil नहीं करना, लेकिन इसमें एक बहुत useful operator शामिल होगा।

### 2.1. Collection command define करें और terminal में test करें

जो collection step हम अपने workflow में add करना चाहते हैं वह `cat` command use करेगा multiple uppercased greetings को एक single file में concatenate करने के लिए।

Terminal में command by itself run करो यह verify करने के लिए कि यह expected तरीके से काम करता है, जैसा हमने पहले किया था।

अपने terminal में निम्नलिखित run करो:

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

Output `COLLECTED-output.txt` नामक एक text file है जिसमें original greetings के uppercase versions हैं।

??? abstract "File contents"

    ```console title="COLLECTED-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

यही result है जो हम अपने workflow से achieve करना चाहते हैं।

### 2.2. Collection step करने के लिए new process बनाएं

चलो एक new process बनाते हैं और इसे `collectGreetings()` call करते हैं।
हम इसे लिखना शुरू कर सकते हैं जो हमने पहले देखा है उसके आधार पर।

#### 2.2.1. Process के 'obvious' parts लिखें

Workflow script में निम्नलिखित process definition add करो:

```groovy title="hello-workflow.nf" linenums="37"
/*
 * Uppercase अभिवादनों को एक single output फ़ाइल में collect करें
 */
process collectGreetings {

    input:
    ???

    output:
    path "COLLECTED-output.txt"

    script:
    """
    cat ??? > 'COLLECTED-output.txt'
    """
}
```

यह वह है जो हम confidence के साथ लिख सकते हैं जो तुमने अब तक सीखा है उसके आधार पर।
लेकिन यह functional नहीं है!
यह input definition(s) और script command का first half छोड़ देता है क्योंकि हमें यह figure out करना होगा कि इसे कैसे लिखें।

#### 2.2.2. `collectGreetings()` के inputs define करें

हमें `convertToUpper()` process की सभी calls से greetings collect करनी होंगी।
हम जानते हैं कि हम workflow के previous step से क्या पा सकते हैं?

`convertToUpper()` द्वारा output किया गया channel उन individual files के paths contain करेगा जिनमें uppercased greetings हैं।
यह एक input slot की राशि है; चलो इसे simplicity के लिए `input_files` call करते हैं।

Process block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          path input_files
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          ???
    ```

Notice करो कि हम `path` prefix use करते हैं भले ही हम expect करते हैं कि इसमें multiple files होंगी।

#### 2.2.3. Concatenation command compose करें

यह वह जगह है जहाँ चीजें थोड़ी tricky हो सकती हैं, क्योंकि हमें arbitrary number की input files handle करने में सक्षम होना होगा।
Specifically, हम command को up front नहीं लिख सकते, इसलिए हमें Nextflow को बताना होगा कि runtime पर इसे कैसे compose करें inputs के आधार पर जो process में flow करते हैं।

दूसरे शब्दों में, यदि हमारे पास element `[file1.txt, file2.txt, file3.txt]` containing करने वाला input channel है, तो हमें Nextflow को उसे `cat file1.txt file2.txt file3.txt` में turn करना होगा।

सौभाग्य से, Nextflow हमारे लिए ऐसा करने में quite happy है यदि हम script command में simply `cat ${input_files}` लिखें।

Process block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="54" hl_lines="3"
        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="54"
        script:
        """
        cat ??? > 'COLLECTED-output.txt'
        """
    ```

Theory में यह किसी भी arbitrary number की input files handle करना चाहिए।

!!! tip

    कुछ command-line tools को प्रत्येक input file के लिए argument (जैसे `-input`) provide करने की आवश्यकता होती है।
    उस case में, हमें command compose करने के लिए थोड़ा extra work करना होगा।
    तुम इसका उदाहरण [Nextflow for Genomics](../../nf4_science/genomics/) training course में देख सकते हो।

### 2.3. Workflow में collection step add करें

अब हमें बस uppercasing step के output पर collection process call करना चाहिए।
वह भी एक channel है, जिसे `convertToUpper.out` call किया जाता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-connector.svg"
</figure>

#### 2.3.1. Process calls connect करें

Workflow block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)

        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out)
    }
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="75"
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
    }
    ```

यह `convertToUpper()` के output को `collectGreetings()` के input से connect करता है।

#### 2.3.2. `-resume` के साथ workflow चलाएं

चलो इसे try करते हैं।

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Command output"

    ```console hl_lines="8"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [mad_gilbert] DSL2 - revision: 6acfd5e28d

    executor >  local (3)
    [79/33b2f0] sayHello (2)         | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (3)   | 3 of 3, cached: 3 ✔
    [47/50fe4a] collectGreetings (1) | 3 of 3 ✔
    ```

यह successfully run होती है, third step सहित।

हालाँकि, last line पर `collectGreetings()` की calls की संख्या देखो।
हम केवल एक expect कर रहे थे, लेकिन तीन हैं।

अब final output file की contents पर नज़र डालो।

??? abstract "File contents"

    ```console title="results/COLLECTED-output.txt"
    Holà
    ```

Oh no. Collection step प्रत्येक greeting पर individually run किया गया था, जो हम नहीं चाहते थे।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-no-operator.svg"
</figure>

हमें Nextflow को explicitly बताने के लिए कुछ करना होगा कि हम चाहते हैं कि third step `convertToUpper()` द्वारा output किए गए channel में सभी elements पर run हो।

### 2.4. Greetings को single input में collect करने के लिए operator use करें

हाँ, एक बार फिर हमारी problem का answer operator है।

Specifically, हम aptly-named [`collect()`](https://nextflow.io/docs/latest/reference/operator.html#collect) operator use करने वाले हैं।

#### 2.4.1. `collect()` operator add करें

इस बार यह थोड़ा different दिखने वाला है क्योंकि हम channel factory के context में operator add नहीं कर रहे हैं; हम इसे output channel में add कर रहे हैं।

हम `convertToUpper.out` लेते हैं और `collect()` operator append करते हैं, जो हमें `convertToUpper.out.collect()` देता है।
हम इसे सीधे `collectGreetings()` process call में plug कर सकते हैं।

Workflow block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out)
    }
    ```

#### 2.4.2. कुछ `view()` statements add करें

चलो channel contents की before और after states visualize करने के लिए कुछ `view()` statements भी include करते हैं।

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out.collect())

        // optional view statements
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    }
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="73"
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out.collect())
    }
    ```

`view()` statements कहीं भी जा सकते हैं जहाँ तुम चाहो; हमने उन्हें readability के लिए call के ठीक बाद रखा।

#### 2.4.3. `-resume` के साथ workflow फिर से चलाएं

चलो इसे try करते हैं:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    Before collect: /workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt
    Before collect: /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt
    Before collect: /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt
    After collect: [/workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt, /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt, /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt]
    ```

यह successfully run होती है, हालाँकि log output थोड़ा messier लग सकता है (हमने readability के लिए इसे clean up किया)।

इस बार third step केवल एक बार call किया गया!
`view()` statements के output को देखते हुए, हम निम्नलिखित देखते हैं:

- तीन `Before collect:` statements, प्रत्येक greeting के लिए एक: उस point पर file paths channel में individual items हैं।
- एक single `After collect:` statement: तीनों file paths अब एक single element में package हैं।

हम इसे निम्नलिखित diagram के साथ summarize कर सकते हैं:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-WITH-operator.svg"
</figure>

Finally, तुम output file की contents पर नज़र डाल सकते हो यह satisfy करने के लिए कि सब कुछ correctly काम किया।

??? abstract "File contents"

    ```console title="results/COLLECTED-output.txt"
    BONJOUR
    HELLO
    HOLà
    ```

इस बार हमारे पास final output file में तीनों greetings हैं। Success!

!!! note

    यदि तुम इसे `-resume` के बिना कई बार run करते हो, तो तुम देखोगे कि greetings का order एक run से दूसरे में बदलता है।
    यह तुम्हें दिखाता है कि जिस order में elements process calls के through flow करते हैं वह consistent होने की guarantee नहीं है।

#### 2.4.4. Readability के लिए `view()` statements remove करें

अगले section में move करने से पहले, हम recommend करते हैं कि तुम console output को clutter करने से बचने के लिए `view()` statements delete कर दो।

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="73"
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out.collect())
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out.collect())

        // optional view statements
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    ```

यह basically point 2.4.2 का reverse operation है।

### सारांश

तुम जानते हो कि batch of process calls से outputs collect करके joint analysis या summation step में कैसे feed करें।

Recap करने के लिए, यह है जो तुमने अब तक build किया है:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

### आगे क्या है?

सीखो कि process को एक से अधिक input कैसे pass करें।

---

## 3. Process को additional parameters pass करें

हम final output file को कुछ specific name देने में सक्षम होना चाहते हैं ताकि greetings के subsequent batches को previous results overwrite किए बिना process किया जा सके।

इसके लिए, हम workflow में निम्नलिखित refinements करने वाले हैं:

- Collector process को output file के लिए user-defined name accept करने के लिए modify करें (`batch_name`)
- Workflow में command-line parameter add करें (`--batch`) और इसे collector process को pass करें

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-batch.svg"
</figure>

### 3.1. Collector process modify करें

हमें additional input declare करना होगा और इसे output file name में integrate करना होगा।

#### 3.1.1. Additional input declare करें

Good news: हम process definition में जितने चाहें उतने input variables declare कर सकते हैं।
चलो इसे `batch_name` call करते हैं।

Process block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="3"
        input:
        path input_files
        val batch_name
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="42"
        input:
        path input_files
    ```

तुम अपने processes को जितने चाहो उतने inputs expect करने के लिए set up कर सकते हो।
Right now, ये सभी required inputs हैं; workflow काम करने के लिए तुम्हें value provide करनी _होगी_।

तुम अपनी Nextflow journey में बाद में सीखोगे कि required vs. optional inputs कैसे manage करें।

#### 3.1.2. Output file name में `batch_name` variable use करें

हम output file name में variable insert कर सकते हैं उसी तरह से जैसा हमने पहले dynamic file names compose किए हैं।

Process block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-${batch_name}-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

यह process को workflow के final output के लिए specific filename generate करने के लिए `batch_name` value use करने के लिए set up करता है।

### 3.2. `batch` command-line parameter add करें

अब हमें `batch_name` के लिए value supply करने और इसे process call में feed करने का तरीका चाहिए।

#### 3.2.1. Parameter set up करने के लिए `params` use करें

तुम already जानते हो कि CLI parameters declare करने के लिए `params` system कैसे use करें।
चलो इसे `batch` parameter declare करने के लिए use करते हैं (एक default value के साथ क्योंकि हम lazy हैं)।

Pipeline parameters section में, निम्नलिखित code changes करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="6"
    /*
     * Pipeline पैरामीटर
     */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="55"
    /*
     * Pipeline पैरामीटर
     */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

बिल्कुल जैसे हमने `--input` के लिए demonstrate किया, तुम command line पर `--batch` के साथ value specify करके उस default value को override कर सकते हो।

#### 3.2.2. Process को `batch` parameter pass करें

Parameter की value process को provide करने के लिए, हमें process call में इसे add करना होगा।

Workflow block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out.collect())
    ```

तुम देखते हो कि process को multiple inputs provide करने के लिए, तुम बस call parentheses में उन्हें list करते हो, commas द्वारा separated।

!!! warning

    तुम्हें process को inputs उसी EXACT ORDER में provide करनी होंगी जिसमें वे process के input definition block में listed हैं।

### 3.3. Workflow चलाएं

चलो command line पर batch name के साथ इसे run करके try करते हैं।

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [confident_rutherford] DSL2 - revision: bc58af409c

    executor >  local (1)
    [79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [b5/f19efe] collectGreetings   | 1 of 1 ✔
    ```

यह successfully run होती है और desired output produce करती है:

??? abstract "File contents"

    ```console title="results/COLLECTED-trio-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

अब, जब तक हम parameter को appropriately specify करते हैं, inputs के other batches पर subsequent runs previous results को clobber नहीं करेंगी।

### सारांश

तुम जानते हो कि process को एक से अधिक input कैसे pass करें।

### आगे क्या है?

सीखो कि multiple outputs कैसे emit करें और उन्हें conveniently कैसे handle करें।

---

## 4. Collector step में output add करें

अब तक हम ऐसे processes use कर रहे थे जो केवल एक-एक output produce करते थे।
हम उनके respective outputs को `<process>.out` syntax use करके बहुत conveniently access कर सके, जिसे हमने output को next process को pass करने के context में use किया (जैसे `convertToUpper(sayHello.out)`) और `publish:` section के context में (जैसे `first_output = sayHello.out`)।

क्या होता है जब process एक से अधिक produce करती है?
हम multiple outputs कैसे handle करते हैं?
क्या हम specific output select और use कर सकते हैं?

सभी excellent questions हैं, और short answer है हाँ हम कर सकते हैं!

Multiple outputs separate channels में package होंगे।
हम या तो उन output channels को names दे सकते हैं, जो बाद में उन्हें individually refer करना आसान बनाता है, या हम उन्हें index द्वारा refer कर सकते हैं।

Demonstration purposes के लिए, चलो कहते हैं कि हम inputs के given batch के लिए collect की जा रही greetings की संख्या count करना और इसे एक file में report करना चाहते हैं।

### 4.1. Process को greetings count और output करने के लिए modify करें

इसके लिए process definition में दो key changes की आवश्यकता होगी: हमें greetings count करने और report file लिखने का तरीका चाहिए, फिर हमें process के `output` block में उस report file को add करना होगा।

#### 4.1.1. Collected greetings की संख्या count करें

Conveniently, Nextflow हमें process definition के `script:` block में arbitrary code add करने देता है, जो इस तरह की चीजें करने के लिए really handy आता है।

इसका मतलब है कि हम `input_files` array में files की संख्या पाने के लिए Nextflow की built-in `size()` function use कर सकते हैं, और result को `echo` command के साथ file में लिख सकते हैं।

`collectGreetings` process block में, निम्नलिखित code changes करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="55"
        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

`count_greetings` variable runtime पर compute होगा।

#### 4.1.2. Report file emit करें और outputs name करें

Principle में हमें बस `output:` block में report file add करनी होगी।

हालाँकि, जबकि हम यह कर रहे हैं, हम अपने output declarations में कुछ `emit:` tags भी add करने वाले हैं। ये हमें indices use करने के बजाय name द्वारा outputs select करने में enable करेंगे।

Process block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 3"
        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="46"
        output:
        path "COLLECTED-${batch_name}-output.txt"
    ```

`emit:` tags optional हैं, और हम केवल outputs में से एक को tag add कर सकते थे।
लेकिन जैसा saying है, why not both?

!!! tip

    यदि तुम `emit:` use करके process के outputs को name नहीं करते, तो तुम अभी भी उन्हें उनके respective (zero-based) index use करके individually access कर सकते हो।
    उदाहरण के लिए, तुम first output पाने के लिए `<process>.out[0]` use करोगे, second output पाने के लिए `<process>.out[1]`, और so on।

    हम outputs name करना prefer करते हैं क्योंकि otherwise, error द्वारा wrong index grab करना बहुत आसान है, especially जब process बहुत सारे outputs produce करती है।

### 4.2. Workflow outputs update करें

अब जब हमारे पास `collectGreetings` process से दो outputs आ रहे हैं, तो `collectGreetings.out` output दो channels contain करता है:

- `collectGreetings.out.outfile` final output file contain करता है
- `collectGreetings.out.report` report file contain करता है

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-report.svg"
</figure>

हमें workflow outputs को accordingly update करना होगा।

#### 4.2.1. `publish:` section update करें

`workflow block` में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4 5"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out
    ```

जैसा तुम देख सकते हो, specific process outputs को refer करना अब trivial है।
जब हम Part 5 (Containers) में अपनी pipeline में एक और step add करने जाएंगे, तो हम आसानी से `collectGreetings.out.outfile` को refer कर सकेंगे और इसे new process को hand कर सकेंगे (spoiler: new process को `cowpy` call किया जाता है)।

लेकिन अभी के लिए, चलो workflow-level outputs update करना finish करते हैं।

#### 4.2.2. `output` block update करें

`output` block में, निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="86" hl_lines="14-17"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
        batch_report {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="80"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

हमें `collected` output definition update करने की जरूरत नहीं है क्योंकि वह name बदला नहीं है।
हमें बस new output add करना होगा।

### 4.3. Workflow चलाएं

चलो greetings के current batch के साथ इसे run करके try करते हैं।

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [ecstatic_wilson] DSL2 - revision: c80285f8c8

    executor >  local (1)
    [c5/4c6ca9] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [0e/6cbc59] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [02/61ead2] collectGreetings   [100%] 1 of 1 ✔
    ```

यदि तुम `results/hello_workflow/` directory में देखो, तो तुम्हें new report file, `trio-report.txt` मिलेगी।
इसे खोलो यह verify करने के लिए कि workflow ने correctly greetings की count report की जो process की गईं।

??? abstract "File contents"

    ```txt title="trio-report.txt"
    There were 3 greetings in this batch.
    ```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-4-way.svg"
</figure>

CSV में और greetings add करके test करने में feel free करो कि क्या होता है।

### सारांश

तुम जानते हो कि process को multiple named outputs emit करवाना और उन्हें workflow level पर appropriately handle करना।

अधिक generally, तुम common ways में processes को एक साथ connect करने में involved key principles समझते हो।

### आगे क्या है?

Extra long break लो, तुमने इसे earn किया है।

जब तुम ready हो, तो [**Part 4: Hello Modules**](./04_hello_modules.md) पर move करो यह सीखने के लिए कि better maintainability और code efficiency के लिए अपने code को कैसे modularize करें।

---

## Quiz

<quiz>
Workflow block में process का output कैसे access करते हो?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

और जानें: [1.4. Pass the output of the first process to the second process](#14-first-process-का-output-second-process-को-pass-करें)
</quiz>

<quiz>
Nextflow में process execution का order क्या determine करता है?
- [ ] Workflow block में processes लिखे जाने का order
- [ ] Process name द्वारा alphabetical order
- [x] Processes के बीच data dependencies
- [ ] Parallel execution के लिए random order

और जानें: [1.4. Pass the output of the first process to the second process](#14-first-process-का-output-second-process-को-pass-करें)
</quiz>

<quiz>
Downstream process के लिए सभी outputs को single list में gather करने के लिए `???` को कौन सा operator replace करना चाहिए?

```groovy hl_lines="4"
workflow {
    greetings_ch = Channel.of('Hello', 'Bonjour', 'Hola')
    SAYHELLO(greetings_ch)
    GATHER_ALL(SAYHELLO.out.???)
}
```

- [ ] `flatten()`
- [x] `collect()`
- [ ] `mix()`
- [ ] `join()`

और जानें: [2.4. Use an operator to collect the greetings into a single input](#24-greetings-को-single-input-में-collect-करने-के-लिए-operator-use-करें)
</quiz>

<quiz>
`collect()` operator कब use करना चाहिए?
- [ ] जब तुम items को parallel में process करना चाहते हो
- [ ] जब तुम्हें channel contents filter करने की जरूरत हो
- [x] जब downstream process को upstream process से सभी items की जरूरत हो
- [ ] जब तुम data को multiple processes में split करना चाहते हो

और जानें: [2.4. Use an operator to collect the greetings into a single input](#24-greetings-को-single-input-में-collect-करने-के-लिए-operator-use-करें)
</quiz>

<quiz>
Process से named output कैसे access करते हो?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

और जानें: [4.1.2. Emit the report file and name outputs](#412-report-file-emit-करें-और-outputs-name-करें)
</quiz>

<quiz>
Process में output name करने के लिए correct syntax क्या है?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

और जानें: [4.1.2. Emit the report file and name outputs](#412-report-file-emit-करें-और-outputs-name-करें)
</quiz>

<quiz>
Process को multiple inputs provide करते समय क्या true होना चाहिए?
- [ ] सभी inputs same type के होने चाहिए
- [ ] Inputs alphabetical order में provide होने चाहिए
- [x] Inputs का order input block में defined order से match होना चाहिए
- [ ] एक बार में केवल दो inputs provide किए जा सकते हैं

और जानें: [3. Pass more than one input to a process](#3-process-को-additional-parameters-pass-करें)
</quiz>
