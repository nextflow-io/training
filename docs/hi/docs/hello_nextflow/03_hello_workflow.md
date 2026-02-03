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
3. एक process को एक से अधिक input pass करना
4. एक process से आने वाले multiple outputs handle करना

Demonstrate करने के लिए, हम Parts 1 और 2 के domain-agnostic Hello World example पर build करना जारी रखेंगे।
इस बार, हम अपने workflow में निम्नलिखित changes करेंगे जो बेहतर reflect करते हैं कि लोग actual workflows कैसे बनाते हैं:

1. Greeting को uppercase में convert करने वाला एक second step add करें।
2. सभी transformed greetings collect करने और उन्हें एक single file में लिखने वाला एक third step add करें।
3. Final output file को name करने के लिए एक parameter add करें और उसे collection step को secondary input के रूप में pass करें।
4. Collection step को process किए गए के बारे में एक simple statistic भी report करवाएं।

??? info "इस section से कैसे शुरू करें"

    Course का यह section मानता है कि तुमने [Hello Nextflow](./index.md) course के Parts 1-2 complete कर लिए हैं, लेकिन यदि तुम उन sections में covered basics से comfortable हो, तो तुम बिना कुछ special किए यहाँ से शुरू कर सकते हो।

---

## 0. Warmup: `hello-workflow.nf` चलाएं

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
    cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

इसमें, हम input filename के आधार पर second output filename compose करते हैं, जैसा हमने originally first process के output के लिए किया था।

### 1.3. Workflow block में new process का call add करें

अब हमें Nextflow को बताना होगा कि actually हमने जो process define किया उसे call करे।

Workflow block में, निम्नलिखित code change करो:

=== "After"

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

=== "Before"

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

Conveniently, Nextflow automatically process का output `<process>.out` नामक channel में package करता है।
तो `sayHello` process का output `sayHello.out` नामक channel है, जिसे हम सीधे `convertToUpper()` के call में plug कर सकते हैं।

Workflow block में, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // अभिवादन को uppercase में बदलें
        convertToUpper()
    ```

इस तरह के simple case (एक output से एक input) के लिए, दो processes connect करने के लिए बस इतना करना होगा!

### 1.5. Workflow output publishing set up करें

Finally, चलो workflow outputs update करते हैं ताकि second process के results भी publish हों।

#### 1.5.1. `workflow` block का `publish:` section update करें

`workflow` block में, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="56" hl_lines="3"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="56"
        publish:
        first_output = sayHello.out
    }
    ```

Logic पहले जैसा ही है।

#### 1.5.2. `output` block update करें

`output` block में, निम्नलिखित code change करो:

=== "After"

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

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="61"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

एक बार फिर, logic पहले जैसा ही है।

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

### सीख

तुम जानते हो कि एक step का output अगले step को input के रूप में provide करके processes को chain कैसे करें।

### आगे क्या?

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

```groovy title="hello-workflow.nf" linenums="37"
/*
 * Uppercase अभिवादनों को एक single output फ़ाइल में collect करें
 */
process collectGreetings {

    input:
    path input_files

    output:
    path "COLLECTED-output.txt"

    script:
    """
    cat ${input_files} > 'COLLECTED-output.txt'
    """
}
```

### 2.3. Workflow में collection step add करें

अब हमें बस uppercasing step के output पर collection process call करना चाहिए।

Workflow block में, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)

        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out)
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="75"
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
    }
    ```

### 2.4. Greetings को single input में collect करने के लिए operator use करें

हमें aptly-named [`collect()`](https://www.nextflow.io/docs/latest/reference/operator.html#collect) operator use करना होगा।

Workflow block में, निम्नलिखित code change करो:

=== "After"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "Before"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out)
    }
    ```

### सीख

तुम जानते हो कि batch of process calls से outputs collect करके joint analysis या summation step में कैसे feed करें।

### आगे क्या?

सीखो कि process को एक से अधिक input कैसे pass करें।

---

## 3. Process को एक से अधिक input pass करें

हम final output file को कुछ specific name देने में सक्षम होना चाहते हैं ताकि greetings के subsequent batches को previous results overwrite किए बिना process किया जा सके।

### 3.1. Collector process modify करें

हमें additional input declare करना होगा और इसे output file name में integrate करना होगा।

```groovy title="hello-workflow.nf" linenums="42" hl_lines="3"
    input:
    path input_files
    val batch_name
```

### 3.2. `batch` command-line parameter add करें

```groovy title="hello-workflow.nf" linenums="55" hl_lines="6"
/*
 * Pipeline पैरामीटर
 */
params {
    input: Path = 'data/greetings.csv'
    batch: String = 'batch'
}
```

### 3.3. Workflow चलाएं

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

### सीख

तुम जानते हो कि process को एक से अधिक input कैसे pass करें।

### आगे क्या?

सीखो कि multiple outputs emit और conveniently handle कैसे करें।

---

## 4. Collector step में output add करें

Multiple outputs separate channels में package होंगे।
हम या तो उन output channels को names दे सकते हैं, जो बाद में उन्हें individually refer करना आसान बनाता है, या हम उन्हें index द्वारा refer कर सकते हैं।

### 4.1. Greetings count और output करने के लिए process modify करें

```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 3"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

### 4.2. Workflow outputs update करें

```groovy title="hello-workflow.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

### 4.3. Workflow चलाएं

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

### सीख

तुम जानते हो कि process को multiple named outputs emit करवाना और उन्हें workflow level पर appropriately handle करना।

अधिक generally, तुम common ways में processes को एक साथ connect करने में involved key principles समझते हो।

### आगे क्या?

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

और जानें: [1.4. Pass the output of the first process to the second process](#14-pass-the-output-of-the-first-process-to-the-second-process)
</quiz>

<quiz>
Nextflow में process execution का order क्या determine करता है?
- [ ] Workflow block में processes लिखे जाने का order
- [ ] Process name द्वारा alphabetical order
- [x] Processes के बीच data dependencies
- [ ] Parallel execution के लिए random order

और जानें: [1.4. Pass the output of the first process to the second process](#14-pass-the-output-of-the-first-process-to-the-second-process)
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

और जानें: [2.4. Use an operator to collect the greetings into a single input](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
Process से named output कैसे access करते हो?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

और जानें: [4.1.2. Emit the report file and name outputs](#412-emit-the-report-file-and-name-outputs)
</quiz>

<quiz>
Process में output name करने के लिए correct syntax क्या है?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

और जानें: [4.1.2. Emit the report file and name outputs](#412-emit-the-report-file-and-name-outputs)
</quiz>

<quiz>
Process को multiple inputs provide करते समय क्या true होना चाहिए?
- [ ] सभी inputs same type के होने चाहिए
- [ ] Inputs alphabetical order में provide होने चाहिए
- [x] Inputs का order input block में defined order से match होना चाहिए
- [ ] एक बार में केवल दो inputs provide किए जा सकते हैं

और जानें: [3. Pass more than one input to a process](#3-pass-more-than-one-input-to-a-process)
</quiz>
