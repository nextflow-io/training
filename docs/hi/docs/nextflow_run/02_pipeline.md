# भाग 2: असली pipelines चलाएं

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

इस कोर्स के भाग 1 (बुनियादी संचालन चलाएं) में, हमने एक उदाहरण workflow से शुरुआत की थी जिसमें code complexity को कम रखने के लिए केवल न्यूनतम features थीं।
उदाहरण के लिए, `1-hello.nf` ने एक समय में एक single value प्रदान करने के लिए एक कमांड-लाइन parameter (`--input`) का उपयोग किया।

हालांकि, अधिकांश real-world pipelines scale पर बड़ी मात्रा में डेटा की कुशल processing enable करने के लिए अधिक sophisticated features का उपयोग करती हैं, और कभी-कभी complex logic द्वारा एक साथ जुड़े कई processing steps apply करती हैं।

प्रशिक्षण के इस भाग में, हम मूल Hello World pipeline के expanded versions आज़माकर real-world pipelines की मुख्य features प्रदर्शित करते हैं।

## 1. एक फ़ाइल से इनपुट डेटा प्रोसेस करना

एक real-world pipeline में, हम आमतौर पर एक या अधिक इनपुट फ़ाइलों में contained कई data points (या data series) को process करना चाहते हैं।
और जहां संभव हो, हम independent data की processing को parallel में चलाना चाहते हैं, ताकि analysis के लिए wait करने में लगने वाले समय को कम किया जा सके।

यह प्रदर्शित करने के लिए कि Nextflow यह कैसे करता है, हमने `greetings.csv` नाम की एक CSV फ़ाइल तैयार की है जिसमें कई इनपुट greetings हैं, जो उस प्रकार के columnar data की नकल करती है जिसे तुम एक real data analysis में process करना चाहोगे।
ध्यान दो कि numbers meaningful नहीं हैं, वे सिर्फ illustrative purposes के लिए हैं।

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

हमने मूल workflow का एक improved version भी लिखा है, जिसे अब `2a-inputs.nf` कहा जाता है, जो CSV फ़ाइल में पढ़ेगा, greetings extract करेगा और उनमें से प्रत्येक को एक अलग फ़ाइल में लिखेगा।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

आइए पहले workflow चलाएं, और बाद में हम relevant Nextflow code पर नज़र डालेंगे।

### 1.1. Workflow चलाओ

अपने terminal में निम्नलिखित कमांड चलाओ।

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

रोमांचक रूप से, यह indicate करता है कि process के लिए '3 of 3' calls किए गए, जो encouraging है, क्योंकि हमने जो CSV इनपुट के रूप में प्रदान किया उसमें तीन rows of data थीं।
यह suggest करता है कि `sayHello()` process को तीन बार call किया गया, प्रत्येक इनपुट row पर एक बार।

### 1.2. `results` डायरेक्टरी में published outputs खोजो

आइए 'results' डायरेक्टरी देखें यह देखने के लिए कि क्या हमारा workflow अभी भी वहां हमारे outputs की एक copy लिख रहा है।

??? abstract "डायरेक्टरी सामग्री"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

हां! हम `2a-inputs` नाम की एक नई डायरेक्टरी देखते हैं जिसमें अलग-अलग नामों वाली तीन output फ़ाइलें हैं, सुविधाजनक रूप से।

तुम उनमें से प्रत्येक को खोल कर satisfy हो सकते हो कि उनमें appropriate greeting string है।

??? abstract "फ़ाइल सामग्री"

    ```console title="results/2a-inputs/Hello-output.txt"
    Hello
    ```

    ```console title="results/2a-inputs/Bonjour-output.txt"
    Bonjour
    ```

    ```console title="results/2a-inputs/Holà-output.txt"
    Holà
    ```

यह confirm करता है कि इनपुट फ़ाइल में प्रत्येक greeting को appropriately process किया गया है।

### 1.3. Original outputs और logs खोजो

तुमने शायद देखा होगा कि ऊपर के console आउटपुट में केवल एक task directory का reference था।
क्या इसका मतलब है कि `sayHello()` के सभी तीन calls उस एक task directory में execute किए गए?

#### 1.3.1. Terminal में दी गई task directory की जांच करो

आइए उस `8e/0eb066` task directory के अंदर देखें।

??? abstract "डायरेक्टरी सामग्री"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

हमें केवल एक greeting के corresponding आउटपुट मिलता है (साथ ही accessory फ़ाइलें यदि हम hidden फ़ाइलों का display enable करें)।

तो यहां क्या हो रहा है?

Default रूप से, ANSI logging system एक ही process के सभी calls के लिए status information एक ही पंक्ति पर लिखती है।
परिणामस्वरूप, इसने हमें console आउटपुट में तीन task directory paths (`8e/0eb066`) में से केवल एक दिखाया।
दो अन्य हैं जो वहां listed नहीं हैं।

#### 1.3.2. Terminal को अधिक details दिखाने दो

हम logging behavior को modify करके process calls की पूरी list देख सकते हैं जैसा कि follows कमांड में `-ansi-log false` जोड़कर:

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "कमांड आउटपुट"

    ```console linenums="1"
    N E X T F L O W  ~  version 25.10.2
    Launching `2a-inputs.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
    [ab/1a8ece] Submitted process > sayHello (1)
    [0d/2cae24] Submitted process > sayHello (2)
    [b5/0df1d6] Submitted process > sayHello (3)
    ```

इस बार हम सभी तीन process runs और उनकी associated work subdirectories आउटपुट में listed देखते हैं।
ANSI logging disable करने से Nextflow को terminal आउटपुट में colours का उपयोग करने से भी रोका गया।

ध्यान दो कि दो logging modes के बीच status report करने का तरीका थोड़ा अलग है।
Condensed mode में, Nextflow report करता है कि calls सफलतापूर्वक पूरी हुईं या नहीं।
इस expanded mode में, यह केवल report करता है कि वे submitted हुईं।

यह confirm करता है कि `sayHello()` process को तीन बार call किया जाता है, और प्रत्येक के लिए एक अलग task directory बनाई जाती है।

यदि हम वहां listed प्रत्येक task directory के अंदर देखें, हम verify कर सकते हैं कि प्रत्येक एक greeting से correspond करती है।

??? abstract "डायरेक्टरी सामग्री"

    ```console title="ab/1a8ece"
    work/ab/1a8ece307e53f03fce689dde904b64/
    └── Hello-output.txt
    ```

    ```console title="0d/2cae24"
    work/0d/2cae2481a53593bc607077c80c9466/
    └── Bonjour-output.txt
    ```

    ```console title="b5/0df1d6"
    work/b5/0df1d642353269909c2ce23fc2a8fa/
    └── Holà-output.txt
    ```

यह confirm करता है कि प्रत्येक process call को अन्य सभी से isolation में execute किया जाता है।
इसके कई फायदे हैं, जिसमें collisions से बचना शामिल है यदि process कोई intermediate files produce करता है जिनके नाम unique नहीं हैं।

!!! tip "सुझाव"

    एक complex workflow, या बड़ी संख्या में inputs के लिए, terminal में पूरी list आउटपुट होना थोड़ा overwhelming हो सकता है, इसलिए लोग सामान्य तौर पर routine usage में `-ansi-log false` का उपयोग नहीं करते।

### 1.4. Workflow code की जांच करो

तो workflow का यह version एक CSV फ़ाइल of inputs को पढ़ने, inputs को separately process करने, और outputs को uniquely नाम देने में सक्षम है।

आइए देखें कि workflow code में यह क्या संभव बनाता है।

??? full-code "पूर्ण code फ़ाइल"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' को एक फ़ाइल में प्रिंट करने के लिए echo का उपयोग करें
    */
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

    /*
    * Pipeline पैरामीटर
    */
    params {
        input: Path
    }

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

    output {
        first_output {
            path '2a-inputs'
            mode 'copy'
        }
    }
    ```

एक बार फिर, तुम्हें code syntax याद करने की जरूरत नहीं है, लेकिन workflow के key components को पहचानना सीखना अच्छा है जो important functionality प्रदान करते हैं।

#### 1.4.1. CSV से इनपुट डेटा load करना

यह सबसे interesting भाग है: हमने कमांड-लाइन से single value लेने से, CSV फ़ाइल लेने, इसे parse करने और इसमें contained individual greetings को process करने में कैसे switch किया?

Nextflow में, हम यह एक **channel** के साथ करते हैं: inputs को कुशलता से handle करने और उन्हें multi-step workflows में एक step से दूसरे में shuttle करने के लिए डिज़ाइन किया गया एक construct, जो built-in parallelism और कई अतिरिक्त benefits प्रदान करता है।

आइए इसे break down करें।

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // एक अभिवादन emit करें
    sayHello(greeting_ch)
```

यह code `greeting_ch` नाम का एक channel बनाता है जो CSV फ़ाइल पढ़ता है, इसे parse करता है, और प्रत्येक row से पहला column extract करता है।
परिणाम एक channel है जिसमें `Hello`, `Bonjour`, और `Holà` हैं।

??? tip "यह कैसे काम करता है?"

    यहां बताया गया है कि उस पंक्ति का plain English में क्या मतलब है:

    - `channel.fromPath` एक **channel factory** है जो file path(s) से एक channel बनाता है
    - `(params.input)` specifies करता है कि filepath कमांड लाइन पर `--input` द्वारा प्रदान किया गया है

    दूसरे शब्दों में, वह पंक्ति Nextflow को बताती है: `--input` के साथ दिए गए filepath को लो और इसकी contents को इनपुट डेटा के रूप में treat करने के लिए तैयार हो जाओ।

    फिर अगली दो पंक्तियां **operators** apply करती हैं जो फ़ाइल की actual parsing और appropriate data structure में data की loading करती हैं:

    - `.splitCsv()` Nextflow को बताता है कि CSV फ़ाइल को rows और columns को represent करने वाले array में parse करें
    - `.map { line -> line[0] }` Nextflow को बताता है कि प्रत्येक row से केवल पहले column में element लें

    तो व्यवहार में, निम्नलिखित CSV फ़ाइल से शुरू करते हुए:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    हमने उसे एक array में transform किया है जो इस तरह दिखता है:

    ```txt title="Array contents"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    और फिर हमने तीन rows में से प्रत्येक से पहला element लिया है और उन्हें एक Nextflow channel में load किया है जिसमें अब `Hello`, `Bonjour`, और `Holà` हैं।

    यदि तुम channels और operators को गहराई से समझना चाहते हो, जिसमें उन्हें खुद कैसे लिखना है, [Hello Nextflow Part 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file) देखो।

#### 1.4.2. प्रत्येक greeting पर process call करो

अगला, workflow के `main:` block की अंतिम पंक्ति में, हम loaded `greeting_ch` channel को `sayHello()` process को इनपुट के रूप में प्रदान करते हैं।

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // एक अभिवादन emit करें
    sayHello(greeting_ch)
```

यह Nextflow को बताता है कि channel में प्रत्येक element पर, _i.e._ प्रत्येक greeting पर individually process run करे।
और क्योंकि Nextflow इस तरह से smart है, यह इन process calls को संभव होने पर parallel में run करेगा, available computing infrastructure पर निर्भर करते हुए।

इस तरह तुम comparatively बहुत कम code के साथ बहुत सारे डेटा (कई samples, या data points, जो भी तुम्हारी research की unit है) की कुशल और scalable processing achieve कर सकते हो।

#### 1.4.3. Outputs को कैसे नाम दिया जाता है

अंत में, यह देखने के लिए process code पर एक quick look लेने लायक है कि हम output फ़ाइलों को uniquely नाम कैसे देते हैं।

```groovy title="2a-inputs.nf" linenums="6" hl_lines="7 11"
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

तुम देखते हो कि, `1-hello.nf` में इस process के version की तुलना में, output declaration और command का relevant bit output फ़ाइल नाम में greeting value शामिल करने के लिए बदल गए हैं।

यह ensure करने का एक तरीका है कि output फ़ाइल नाम common results डायरेक्टरी में publish होने पर collide नहीं होंगे।

और यही एकमात्र change है जो हमें process declaration के अंदर करना पड़ा!

### सीख

तुम basic level पर समझते हो कि channels और operators हमें कई inputs को कुशलता से process करने में कैसे enable करते हैं।

### आगे क्या?

जानो कि multi-step workflows कैसे constructed होते हैं और वे कैसे operate करते हैं।

---

## 2. Multi-step workflows चलाना

अधिकांश real-world workflows में एक से अधिक step होते हैं।
आइए जो हमने अभी channels के बारे में सीखा उस पर build करें, और देखें कि Nextflow multi-step workflow में processes को एक साथ जोड़ने के लिए channels और operators का उपयोग कैसे करता है।

उस end के लिए, हम तुम्हें एक example workflow प्रदान करते हैं जो तीन अलग-अलग steps को एक साथ chain करता है और निम्नलिखित demonstrate करता है:

1. एक process से अगले में data flow कराना
2. कई process calls से outputs को एक single process call में collect करना

Specifically, हमने workflow का एक expanded version बनाया है जिसे `2b-multistep.nf` कहा जाता है जो प्रत्येक इनपुट greeting लेता है, इसे uppercase में convert करता है, फिर सभी uppercased greetings को एक single output फ़ाइल में collect करता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

पहले की तरह, हम पहले workflow चलाएंगे फिर code देखेंगे कि क्या नया है।

### 2.1. Workflow चलाओ

अपने terminal में निम्नलिखित कमांड चलाओ:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "कमांड आउटपुट"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

तुम देखते हो कि जैसा कि promised था, workflow के भाग के रूप में कई steps run किए गए; पहले दो (`sayHello` और `convertToUpper`) presumably प्रत्येक individual greeting पर run किए गए, और तीसरा (`collectGreetings`) सभी तीन `convertToUpper` calls के outputs पर केवल एक बार run किया गया होगा।

### 2.2. Outputs खोजो

आइए `results` डायरेक्टरी पर एक नज़र डालकर verify करें कि वास्तव में यही हुआ।

??? abstract "डायरेक्टरी सामग्री"

    ```console linenums="1" hl_lines="8-16"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── batch-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt

    ```

जैसा कि तुम देख सकते हो, हमारे पास `2b-multistep` नाम की एक नई डायरेक्टरी है, और इसमें पहले की तुलना में काफी अधिक फ़ाइलें हैं।
कुछ फ़ाइलों को `intermediates` नाम की एक subdirectory में group किया गया है, जबकि दो फ़ाइलें top level पर स्थित हैं।

वे दो multi-step workflow के final results हैं।
फ़ाइल नाम देखने और उनकी contents check करने के लिए एक मिनट लो ताकि confirm हो सके कि वे वही हैं जो तुम expect करते हो।

??? abstract "फ़ाइल सामग्री"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

पहली में हमारी तीन greetings हैं, जैसा कि promised था uppercased और एक single फ़ाइल में वापस collect की गई।
दूसरी एक report फ़ाइल है जो run के बारे में कुछ information summarize करती है।

### 2.3. Code की जांच करो

आइए code देखें और multi-step workflows के लिए key patterns identify करें।

??? full-code "पूर्ण code फ़ाइल"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' को एक फ़ाइल में प्रिंट करने के लिए echo का उपयोग करें
    */
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

    /*
    * Use a text replacement tool to convert the greeting to uppercase
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

    /*
    * Collect uppercase greetings into a single output file
    */
    process collectGreetings {

        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report

        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    }

    /*
    * Pipeline पैरामीटर
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        collected {
            path '2b-multistep'
            mode 'copy'
        }
        batch_report {
            path '2b-multistep'
            mode 'copy'
        }
    }
    ```

वहां बहुत कुछ हो रहा है, लेकिन workflow के पिछले version की तुलना में सबसे obvious difference यह है कि अब कई process definitions हैं, और correspondingly, workflow block में कई process calls हैं।

आइए करीब से देखें और देखें कि क्या हम सबसे interesting pieces identify कर सकते हैं।

#### 2.3.1. Workflow structure को visualize करना

यदि तुम Nextflow extension के साथ VSCode का उपयोग कर रहे हो, तो किसी भी Nextflow script में workflow block के ठीक ऊपर displayed छोटे `DAG preview` link पर click करके तुम processes कैसे connected हैं इसका एक helpful diagram प्राप्त कर सकते हो।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

यह तुम्हें एक अच्छा overview देता है कि processes कैसे connected हैं और वे क्या produce करते हैं।

तुम देखते हो कि original `sayHello` process के अलावा, अब हमारे पास `convertToUpper` और `collectGreetings` भी हैं, जो उन processes के नामों से match करते हैं जो हमने console आउटपुट में देखे।
दो नई process definitions `sayHello` process के समान तरीके से structured हैं, सिवाय इसके कि `collectGreetings` एक अतिरिक्त input parameter `batch` लेता है और दो outputs produce करता है।

हम प्रत्येक के लिए code में detail में नहीं जाएंगे, लेकिन यदि तुम curious हो, तो तुम [Hello Nextflow के Part 2](../hello_nextflow/03_hello_workflow.md) में details देख सकते हो।

अभी के लिए, आइए dig करें कि processes एक दूसरे से कैसे connected हैं।

#### 2.3.2. Processes कैसे connected हैं

यहां देखने की वास्तव में interesting बात यह है कि workflow के `main:` block में process calls कैसे एक साथ chained हैं।

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // एक अभिवादन emit करें
    sayHello(greeting_ch)
    // अभिवादन को uppercase में बदलें
    convertToUpper(sayHello.out)
    // सभी अभिवादनों को एक फ़ाइल में collect करें
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

तुम देख सकते हो कि पहला process call, `sayHello(greeting_ch)`, unchanged है।
फिर अगला process call, `convertToUpper` को, `sayHello` के output को `sayHello.out` के रूप में refer करता है।

Pattern सरल है: `processName.out` एक process के output channel को refer करता है, जिसे सीधे अगले process को pass किया जा सकता है।
इस तरह हम Nextflow में एक step से दूसरे में data shuttle करते हैं।

#### 2.3.3. एक process कई inputs ले सकता है

तीसरा process call, `collectGreetings` को, थोड़ा अलग है।

```groovy title="2b-multistep.nf" linenums="77"
    // सभी अभिवादनों को एक फ़ाइल में collect करें
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

तुम देखते हो कि इस call को दो inputs दिए गए हैं, `convertToUpper.out.collect()` और `params.batch`।
अभी के लिए `.collect()` bit को ignore करते हुए, हम इसे `collectGreetings(input1, input2)` के रूप में generalize कर सकते हैं।

यह process module में दो input declarations से match करता है:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

जब Nextflow इसे parse करता है, यह call में पहले input को `path input_files` को assign करेगा, और दूसरे को `val batch_name` को।

तो अब तुम जानते हो कि एक process कई inputs ले सकता है, और workflow block में call कैसी दिखती है।

अब आइए उस पहले input, `convertToUpper.out.collect()` पर करीब से नज़र डालें।

#### 2.3.4. `collectGreetings` call में `collect()` क्या करता है

`sayHello` के output को `convertToUpper` को pass करने के लिए, हमने simply `sayHello` के output channel को `sayHello.out` के रूप में refer किया। लेकिन अगले step के लिए, हम `convertToUpper.out.collect()` का reference देख रहे हैं।

यह `collect()` bit क्या है और यह क्या करता है?

यह एक operator है, बिल्कुल। जैसे `splitCsv` और `map` operators जो हमने पहले encounter किए थे।
इस बार operator को `collect` कहा जाता है, और यह `convertToUpper` द्वारा produced output channel पर apply होता है।

`collect` operator का उपयोग एक ही process के कई calls से outputs को collect करने और उन्हें एक single channel element में package करने के लिए किया जाता है।

इस workflow के context में, यह `convertToUpper.out` channel में तीन uppercased greetings ले रहा है --जो तीन अलग-अलग channel items हैं, और normally अगले process द्वारा अलग-अलग calls में handle की जाएंगी-- और उन्हें एक single item में package कर रहा है।

अधिक practical terms में: यदि हम `collectGreetings()` को feed करने से पहले `convertToUpper()` के output पर `collect()` apply नहीं करते, Nextflow simply प्रत्येक greeting पर independently `collectGreetings()` run करेगा, जो हमारे goal को achieve नहीं करेगा।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

इसके विपरीत, `collect()` का उपयोग करना हमें workflow के दूसरे step द्वारा produced सभी अलग-अलग uppercased greetings लेने और उन्हें सब एक साथ pipeline के तीसरे step में एक single call को feed करने की अनुमति देता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

इस तरह हम सभी greetings को वापस एक ही फ़ाइल में लाते हैं।

Process calls के बीच channels की contents पर transformations apply करने के लिए कई अन्य [operators](https://www.nextflow.io/docs/latest/reference/operator.html#operator-page) उपलब्ध हैं।

यह pipeline developers को उनकी pipeline की flow logic customize करने में बहुत flexibility देता है।
Downside यह है कि यह कभी-कभी pipeline क्या कर रही है यह decipher करना कठिन बना सकता है।

#### 2.3.5. एक input parameter का default value हो सकता है

तुमने शायद देखा होगा कि `collectGreetings` एक दूसरा input लेता है, `params.batch`:

```groovy title="2b-multistep.nf" linenums="77"
    // सभी अभिवादनों को एक फ़ाइल में collect करें
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

यह workflow को `--batch` नाम का एक CLI parameter pass करता है।
हालांकि, जब हमने पहले workflow launch किया, हमने `--batch` parameter specify नहीं किया।

वहां क्या हो रहा है?
`params` block पर एक नज़र डालो:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

Workflow में एक default value configured है, इसलिए हमें इसे provide करने की जरूरत नहीं है।
लेकिन यदि हम कमांड लाइन पर एक provide करते हैं, तो default के बजाय जो value हम specify करते हैं वह use होगी।

इसे try करो:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv --batch test
```

??? success "कमांड आउटपुट"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [a5/cdff26] sayHello (1)       | 3 of 3 ✔
    [c5/78794f] convertToUpper (2) | 3 of 3 ✔
    [d3/b4d86c] collectGreetings   | 1 of 1 ✔
    ```

तुम्हें अपने custom batch name के साथ named नए final outputs देखने चाहिए।

??? abstract "डायरेक्टरी सामग्री"

    ```console linenums="1" hl_lines="10 12"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── COLLECTED-test-output.txt
        ├── batch-report.txt
        ├── test-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

यह input configuration का एक aspect है, जिसे हम Part 3 में अधिक detail में cover करेंगे, लेकिन अभी के लिए important बात यह जानना है कि input parameters को default values दी जा सकती हैं।

#### 2.3.6. एक process कई outputs produce कर सकता है

`collectGreetings` process definition में, हम निम्नलिखित output declarations देखते हैं:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

जिन्हें फिर `publish:` block में `emit:` के साथ दिए गए name से refer किया जाता है:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

यह specific outputs को workflow में अन्य processes को individually pass करना आसान बनाता है, विभिन्न operators के combination में।

#### 2.3.7. Published outputs को organize किया जा सकता है

`output` block में, हमने intermediate results को group करने के लिए custom paths का उपयोग किया है ताकि workflow के सिर्फ final outputs को pick out करना आसान हो।

```groovy title="2b-multistep.nf" linenums="87" hl_lines="3 7 11 15"
output {
    first_output {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    uppercased {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    collected {
        path '2b-multistep'
        mode 'copy'
    }
    batch_report {
        path '2b-multistep'
        mode 'copy'
    }
}
```

Published outputs को organize करने के और भी sophisticated तरीके हैं; हम configuration पर part में कुछ touch करेंगे।

!!! tip "Workflows build करने के बारे में अधिक जानना चाहते हो?"

    Multi-step workflows build करने की detailed coverage के लिए, [Hello Nextflow Part 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md) देखो।

### सीख

तुम basic level पर समझते हो कि channels और operators का उपयोग करके multi-step workflows कैसे constructed होते हैं और वे कैसे operate करते हैं।
तुमने यह भी देखा है कि processes कई inputs ले सकते हैं और कई outputs produce कर सकते हैं, और इन्हें structured तरीके से publish किया जा सकता है।

### आगे क्या?

जानो कि code reuse और maintainability को promote करने के लिए Nextflow pipelines को कैसे modularize किया जा सकता है।

---

## 3. Modularized pipelines चलाना

अब तक, हमने जो सभी workflows देखे हैं उनमें एक single workflow फ़ाइल थी जिसमें सभी relevant code था।

हालांकि, real-world pipelines आमतौर पर _modularized_ होने से benefit करती हैं, मतलब code अलग-अलग फ़ाइलों में split है।
यह उनके development और maintenance को अधिक efficient और sustainable बना सकता है।

यहां हम Nextflow में code modularity के सबसे common form को demonstrate करने जा रहे हैं, जो **modules** का उपयोग है।

Nextflow में, एक **module** एक single process definition है जो एक standalone code फ़ाइल में अपने आप में encapsulated है।
Workflow में module use करने के लिए, तुम बस अपनी workflow code फ़ाइल में एक single-line import statement जोड़ते हो; फिर तुम process को workflow में normally integrate कर सकते हो।
यह workflow code की multiple copies produce किए बिना process definitions को multiple workflows में reuse करना possible बनाता है।

अब तक हम ऐसे workflows चला रहे थे जिनके सभी processes एक monolithic code फ़ाइल में included थे।
अब हम देखने जा रहे हैं कि जब processes individual modules में stored होते हैं तो यह कैसा दिखता है।

हमने बेशक demonstration purposes के लिए एक suitable workflow फिर से तैयार की है, जिसे `2c-modules.nf` कहा जाता है, साथ ही `modules/` डायरेक्टरी में located modules का एक set।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "डायरेक्टरी सामग्री"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

तुम देखते हो कि चार Nextflow फ़ाइलें हैं, प्रत्येक processes में से एक के नाम पर।
तुम अभी के लिए `cowpy.nf` फ़ाइल को ignore कर सकते हो; हम बाद में उस पर आएंगे।

### 3.1. Code की जांच करो

इस बार हम पहले code देखने जा रहे हैं।
`2c-modules.nf` workflow फ़ाइल open करके शुरू करो।

??? full-code "पूर्ण code फ़ाइल"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Modules को include करें

include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline पैरामीटर
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        collected {
            path '2c-modules'
            mode 'copy'
        }
        batch_report {
            path '2c-modules'
            mode 'copy'
        }
    }
    ```

तुम देखते हो कि workflow logic पिछले workflow version के exactly same है।
हालांकि, process code workflow फ़ाइल से गायब है, और इसके बजाय `modules` के अंतर्गत अलग फ़ाइलों की ओर point करने वाले `include` statements हैं।

```groovy title="hello-modules.nf" linenums="3"
// Modules को include करें
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

उन फ़ाइलों में से एक open करो और तुम्हें corresponding process का code मिलेगा।

??? full-code "पूर्ण code फ़ाइल"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' को एक फ़ाइल में प्रिंट करने के लिए echo का उपयोग करें
    */
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

जैसा कि तुम देख सकते हो, process code change नहीं हुआ है; इसे बस main workflow फ़ाइल में होने के बजाय एक individual module फ़ाइल में copy किया गया है।
वही अन्य दो processes पर भी लागू होता है।

तो आइए देखें कि इस नए version को run करना कैसा दिखता है।

### 3.2. Workflow चलाओ

अपने terminal में यह कमांड चलाओ, `-resume` flag के साथ:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [j6/cdfa66] sayHello (1)       | 3 of 3, cached: ✔
    [95/79484f] convertToUpper (2) | 3 of 3, cached: ✔
    [5e/4358gc] collectGreetings   | 1 of 1, cached: ✔
    ```

तुम देखोगे कि सभी process executions successfully cached हुए, मतलब Nextflow ने recognize किया कि यह requested work पहले ही कर चुका है, भले ही code split हो गया है और main workflow फ़ाइल का नाम बदल गया है।

इसमें से कुछ भी Nextflow को matter नहीं करता; जो matter करता है वह job script है जो सभी code को एक साथ pull और evaluate करने के बाद generate होती है।

!!! tip "सुझाव"

    Workflow के एक section को 'subworkflow' के रूप में encapsulate करना भी possible है जिसे एक larger pipeline में import किया जा सकता है, लेकिन यह इस course के scope के बाहर है।

    तुम [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/) पर Side Quest में composable workflows develop करने के बारे में अधिक जान सकते हो।

### सीख

तुम जानते हो कि code reuse को promote करने और maintainability improve करने के लिए processes को standalone modules में कैसे store किया जा सकता है।

### आगे क्या?

Software dependencies manage करने के लिए containers का उपयोग करना सीखो।

---

## 4. Containerized software का उपयोग करना

अब तक जो workflows हम examples के रूप में use कर रहे थे उन्हें बस हमारे environment में उपलब्ध UNIX tools का उपयोग करके बहुत basic text processing operations run करने की जरूरत थी।

हालांकि, real-world pipelines को आमतौर पर specialized tools और packages की आवश्यकता होती है जो अधिकांश environments में default रूप से included नहीं होते।
आमतौर पर, तुम्हें इन tools को install करना होगा, उनकी dependencies manage करनी होंगी, और किसी भी conflicts को resolve करना होगा।

यह सब बहुत tedious और annoying है।
इस समस्या को address करने का एक बेहतर तरीका **containers** का उपयोग करना है।

एक **container** एक lightweight, standalone, executable unit of software है जो container **image** से बनाई जाती है जिसमें code, system libraries और settings सहित application run करने के लिए आवश्यक सब कुछ शामिल होता है।

!!! tip "सुझाव"

    हम यह [Docker](https://www.docker.com/get-started/) technology का उपयोग करके सिखाते हैं, लेकिन Nextflow [कई अन्य container technologies](https://www.nextflow.io/docs/latest/container.html#) को भी support करता है।

### 4.1. Container को directly use करो

पहले, आइए एक container के साथ directly interact करने की कोशिश करें।
यह तुम्हारी understanding को solidify करने में मदद करेगा कि containers क्या हैं इससे पहले कि हम उन्हें Nextflow में use करना शुरू करें।

#### 4.1.1. Container image pull करो

Container use करने के लिए, तुम आमतौर पर container registry से container image download या "pull" करते हो, और फिर container instance create करने के लिए container image run करते हो।

General syntax इस प्रकार है:

```bash title="Syntax"
docker pull '<container>'
```

- `docker pull` container system को container repository से container image pull करने का instruction है।
- `'<container>'` container image का URI address है।

एक example के रूप में, आइए एक container image pull करें जिसमें [cowpy](https://github.com/jeffbuttars/cowpy) है, `cowsay` नाम के tool का एक python implementation जो मज़ेदार तरीके से arbitrary text inputs display करने के लिए ASCII art generate करता है।

विभिन्न repositories हैं जहां तुम published containers पा सकते हो।
हमने इस Docker container image को `cowpy` Conda package से generate करने के लिए [Seqera Containers](https://seqera.io/containers/) service का उपयोग किया: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`।

Complete pull कमांड चलाओ:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "कमांड आउटपुट"

    ```console
    Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
    131d6a1b707a8e65: Pulling from library/cowpy
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
    622dd7f15040: Pull complete
    895fb5d0f4df: Pull complete
    Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

यह system को specified image download करने के लिए कहता है।
Download complete होने के बाद, तुम्हारे पास container image की एक local copy है।

#### 4.1.2. Container spin up करो

Containers को one-off कमांड के रूप में run किया जा सकता है, लेकिन तुम उन्हें interactively भी use कर सकते हो, जो तुम्हें container के अंदर एक shell prompt देता है और तुम्हें command के साथ play करने की अनुमति देता है।

General syntax इस प्रकार है:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'` container system को container image से container instance spin up करने और उसमें एक कमांड execute करने का instruction है।
- `--rm` system को कमांड complete होने के बाद container instance shut down करने के लिए कहता है।

Fully assembled, container execution कमांड इस तरह दिखता है:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

वह कमांड चलाओ, और तुम्हें अपना prompt कुछ इस तरह बदलता देखना चाहिए `(base) root@b645838b3314:/tmp#`, जो indicate करता है कि तुम अब container के अंदर हो।

तुम directory contents list करने के लिए `ls` चलाकर इसे verify कर सकते हो:

```bash
ls /
```

??? success "कमांड आउटपुट"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

तुम देखते हो कि container के अंदर filesystem तुम्हारे host system पर filesystem से अलग है।

!!! tip "सुझाव"

    जब तुम एक container run करते हो, यह default रूप से host system से isolated होता है।
    इसका मतलब है कि container host system पर किसी भी files को access नहीं कर सकता जब तक तुम explicitly इसे ऐसा करने की अनुमति नहीं देते यह specify करके कि तुम `docker run` कमांड के भाग के रूप में following syntax का उपयोग करके एक volume mount करना चाहते हो:

    ```bash title="Syntax"
    -v <outside_path>:<inside_path>
    ```

    यह effectively container wall के through एक tunnel establish करता है जिसका उपयोग तुम अपने filesystem के उस part को access करने के लिए कर सकते हो।

    यह [Hello Nextflow के Part 5](../hello_nextflow/05_hello_containers.md) में अधिक detail में cover किया गया है।

#### 4.1.3. `cowpy` tool चलाओ

Container के अंदर से, तुम `cowpy` कमांड directly run कर सकते हो।

```bash
cowpy "Hello Containers"
```

??? success "कमांड आउटपुट"

    ```console
    ______________________________________________________
    < Hello Containers >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

यह default cow character (या 'cowacter') की ASCII art produce करता है जिसमें हमारे द्वारा specified text वाला speech bubble है।

अब जब तुमने basic usage test कर लिया है, तुम इसे कुछ parameters देने की कोशिश कर सकते हो।
उदाहरण के लिए, tool documentation कहता है कि हम `-c` के साथ character set कर सकते हैं।

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

इस बार ASCII art आउटपुट Linux penguin, Tux, दिखाता है, क्योंकि हमने `-c tux` parameter specify किया।

चूंकि तुम container के अंदर हो, तुम input parameters vary करते हुए cowpy कमांड जितनी बार चाहो run कर सकते हो, अपने system पर किसी भी libraries को install करने की चिंता किए बिना।

??? tip "अन्य उपलब्ध characters"

    एक अलग character pick करने के लिए '-c' flag use करो, जिसमें शामिल हैं:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

इसके साथ free feel करके play करो।
जब तुम done हो जाओ, `exit` कमांड का उपयोग करके container exit करो:

```bash
exit
```

तुम अपने normal shell में वापस पाओगे।

### 4.2. Workflow में container use करो

जब हम एक pipeline run करते हैं, हम Nextflow को बताना चाहते हैं कि प्रत्येक step पर कौन सा container use करना है, और importantly, हम चाहते हैं कि यह वह सारा काम handle करे जो हमने अभी किया: container pull करो, इसे spin up करो, कमांड run करो और जब यह done हो जाए तो container tear down करो।

अच्छी खबर: यही exactly है जो Nextflow हमारे लिए करने जा रहा है।
हमें बस प्रत्येक process के लिए एक container specify करना है।

यह काम कैसे करता है demonstrate करने के लिए, हमने अपनी workflow का एक और version बनाया है जो तीसरे step में produced collected greetings की फ़ाइल पर `cowpy` run करता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-cowpy.svg"
</figure>

इसे speech bubble में तीन greetings के साथ ASCII art containing एक फ़ाइल output करनी चाहिए।

#### 4.2.1. Code की जांच करो

Workflow पिछले वाले के बहुत similar है, plus `cowpy` run करने का extra step।

??? full-code "पूर्ण code फ़ाइल"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // Modules को include करें

include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

    /*
    * Pipeline पैरामीटर
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
        // सभी अभिवादनों को एक फ़ाइल में collect करें
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // cowpy के साथ अभिवादनों का ASCII art जनरेट करें
        cowpy(collectGreetings.out.outfile, params.character)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }

    output {
        first_output {
            path '2d-container/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2d-container/intermediates'
            mode 'copy'
        }
        collected {
            path '2d-container/intermediates'
            mode 'copy'
        }
        batch_report {
            path '2d-container'
            mode 'copy'
        }
        cowpy_art {
            path '2d-container'
            mode 'copy'
        }
    }
    ```

तुम देखते हो कि यह workflow एक module फ़ाइल से `cowpy` process import करती है, और इसे `collectGreetings()` call के output पर call करती है, plus `params.character` नाम का एक input parameter।

```groovy title="2d-container.nf" linenums="25"
// cowpy के साथ ASCII art जनरेट करें
cowpy(collectGreetings.out, params.character)
```

`cowpy` process, जो ASCII art generate करने के लिए cowpy कमांड wrap करता है, `cowpy.nf` module में defined है।

??? full-code "पूर्ण code फ़ाइल"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // cowpy के साथ ASCII art जनरेट करें (https://github.com/jeffbuttars/cowpy)
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

`cowpy` process को दो inputs चाहिए: speech bubble में डालने के लिए text containing input फ़ाइल का path (`input_file`), और character variable के लिए एक value।

Importantly, इसमें `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'` पंक्ति भी शामिल है, जो उस container URI की ओर point करती है जिसका हमने पहले use किया था।

#### 4.2.2. Check करो कि Docker configuration में enabled है

हम इस प्रशिक्षण course के Part 3 को slightly anticipate करने जा रहे हैं `nextflow.config` configuration फ़ाइल introduce करके, जो Nextflow workflow execution configure करने के लिए offer करता है main ways में से एक है।
जब `nextflow.config` नाम की एक फ़ाइल current डायरेक्टरी में present होती है, Nextflow automatically इसे load करेगा और इसमें contained किसी भी configuration को apply करेगा।

उस end के लिए, हमने एक `nextflow.config` फ़ाइल include की है जिसमें Docker enable करने वाली code की एक single पंक्ति है।

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

यह configuration Nextflow को किसी भी process के लिए Docker use करने के लिए कहती है जो एक compatible container specify करता है।

!!! tip "सुझाव"

    `-with-docker <container>` parameter का उपयोग करके कमांड-लाइन से, per-run basis पर Docker execution enable करना technically possible है।
    हालांकि, यह हमें केवल entire workflow के लिए एक container specify करने की अनुमति देता है, जबकि जो approach हमने तुम्हें अभी दिखाया वह हमें प्रति process एक अलग container specify करने की अनुमति देता है।
    बाद वाला modularity, code maintenance और reproducibility के लिए बहुत बेहतर है।

#### 4.2.3. Workflow चलाओ

Recap करने के लिए, यह वह है जो हम run करने वाले हैं:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

क्या तुम्हें लगता है यह काम करेगा?

आइए `-resume` flag के साथ workflow run करें, और specify करें कि हम character को turkey चाहते हैं।

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

पहले तीन steps cached हुए क्योंकि हमने उन्हें पहले run कर लिया था, लेकिन `cowpy` process नया है इसलिए वह actually run होता है।

तुम `results` डायरेक्टरी में `cowpy` step का output पा सकते हो।

??? abstract "फ़ाइल सामग्री"

    ```console title="results/2d-container/cowpy-COLLECTED-batch-output.txt"
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

तुम देखते हो कि character सभी greetings बोल रहा है, क्योंकि यह collected uppercased greetings की फ़ाइल पर run किया गया।

अधिक महत्वपूर्ण बात, हम इसे cowpy और इसकी सभी dependencies का proper installation किए बिना अपनी pipeline के भाग के रूप में run करने में सक्षम थे।
और अब हम pipeline को collaborators के साथ share कर सकते हैं और उन्हें इसे उनके infrastructure पर run करवा सकते हैं बिना उन्हें Docker या इसके alternatives में से किसी एक (जैसे Singularity/Apptainer) के अलावा कुछ भी install करने की जरूरत के, जैसा कि ऊपर mentioned है।

#### 4.2.4. Inspect करो कि Nextflow ने containerized task कैसे launch किया

इस section के एक final coda के रूप में, आइए `cowpy` process calls में से एक के लिए work subdirectory पर एक नज़र डालें ताकि यह थोड़ा और insight मिले कि Nextflow hood के नीचे containers के साथ कैसे काम करता है।

`cowpy` process के लिए work subdirectory का path खोजने के लिए अपने `nextflow run` कमांड का output check करो।
ऊपर दिखाए गए run के लिए जो हमें मिला उसे देखते हुए, `cowpy` process के लिए console log पंक्ति `[7f/caf718]` से शुरू होती है।
यह निम्नलिखित truncated डायरेक्टरी path से correspond करता है: `work/7f/caf718`।

उस डायरेक्टरी में, तुम्हें `.command.run` फ़ाइल मिलेगी जिसमें वे सभी commands हैं जो Nextflow ने pipeline execute करते समय तुम्हारी ओर से run किए।

??? abstract "फ़ाइल सामग्री"

    ```console title="work/7f/caf71890cce1667c094d880f4b6dcc/.command.run"
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
    ...
    ```

तुम देख सकते हो कि इसमें container को pull और mount करने के लिए setup instructions शामिल हैं, साथ ही script अंदर चलाने के instructions।

### सीख

तुम जानते हो कि Nextflow processes को container images specify करके software dependencies provide करना कैसे संभव बनाता है।

### आगे क्या?

यह Part 2 के अंत को mark करता है।
अगले part में, तुम workflow execution को configure करना सीखोगे जैसे inputs और outputs manage करना, software packaging technologies के बीच switch करना, और HPC या cloud जैसे different execution platforms पर run करना।

[Part 3: Run configuration](./03_config.md) पर जारी रखो।
