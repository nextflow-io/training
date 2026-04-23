# भाग 3: Run configuration

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

यह section explore करेगा कि pipeline के व्यवहार को customize करने, इसे different environments में adapt करने, और resource usage को optimize करने के लिए Nextflow pipeline की configuration कैसे manage करें _workflow code की एक भी पंक्ति बदले बिना_।

ऐसा करने के कई तरीके हैं, जिन्हें combination में use किया जा सकता है और [Configuration](https://nextflow.io/docs/latest/config.html) documentation में described precedence के order के अनुसार interpret किया जाता है।

इस course के इस part में, हम तुम्हें सबसे simple और सबसे common configuration फ़ाइल mechanism, `nextflow.config` फ़ाइल, दिखाने जा रहे हैं, जिसे तुमने Part 2 में containers पर section में पहले ही encounter किया था।

हम Nextflow configuration के essential components जैसे process directives, executors, profiles, और parameter files पर जाएंगे।
इन configuration options को effectively utilize करना सीखकर, तुम Nextflow pipelines की flexibility, scalability, और performance का पूरा लाभ उठा सकते हो।

इन configuration elements को exercise करने के लिए, हम इस प्रशिक्षण course के Part 2 के अंत में जो workflow आखिरी बार चलाया था उसकी एक fresh copy चलाने जा रहे हैं, जिसका नाम `3-main.nf` रखा गया है।

यदि तुम Hello pipeline से familiar नहीं हो या तुम्हें reminder की जरूरत हो, [यह info page](../info/hello_pipeline.md) देखो।

---

## 1. Workflow input parameters manage करें

??? example "परिदृश्य"

    तुमने एक pipeline download की है और इसे बार-बार same input files और settings के साथ run करना चाहते हो, लेकिन तुम हर बार सभी parameters type नहीं करना चाहते।
    या शायद तुम pipeline को एक colleague के लिए set up कर रहे हो जो command-line arguments के साथ comfortable नहीं है।

हम configuration के एक aspect से शुरू करने जा रहे हैं जो simply अब तक हम जो कर रहे थे उसका एक extension है: input parameters का management।

Currently, हमारी workflow कई parameter values command-line के माध्यम से accept करने के लिए set up है, जो workflow script में ही एक `params` block में declared हैं।
एक की default value उसकी declaration के भाग के रूप में set है।

हालांकि, तुम सभी के लिए defaults set करना चाह सकते हो, या existing default को override करना चाह सकते हो बिना कमांड लाइन पर parameters specify किए, या original script फ़ाइल modify किए।

ऐसा करने के कई तरीके हैं; हम तुम्हें तीन basic तरीके दिखाने जा रहे हैं जो बहुत commonly use होते हैं।

### 1.1. `nextflow.config` में values set up करें

यह सबसे simple approach है, हालांकि यह possibly least flexible है क्योंकि main `nextflow.config` फ़ाइल कुछ ऐसी नहीं है जिसे तुम हर run के लिए edit करना चाहते हो।
लेकिन इसका advantage है कि यह workflow में parameters _declare_ करने (जो definitely वहां belong करता है) बनाम _default values_ supply करने की concerns को separate करता है, जो configuration फ़ाइल में अधिक home पर हैं।

आइए इसे दो steps में करें।

#### 1.1.1. Configuration फ़ाइल में एक `params` block बनाओ

`nextflow.config` फ़ाइल में निम्नलिखित code changes करो:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

ध्यान दो कि हमने simply workflow से `params` block को configuration फ़ाइल में copy नहीं किया।
`batch` parameter जिसकी default value पहले से declared थी, उसके लिए syntax थोड़ा different है।
Workflow फ़ाइल में, वह एक typed declaration है।
Configuration में, वे value assignments हैं।

Technically, यह workflow फ़ाइल में still specified default values को override करने के लिए sufficient है।
तुम `batch` के लिए default value modify कर सकते हो और workflow run करके satisfy हो सकते हो कि configuration फ़ाइल में set value workflow फ़ाइल में set को override करती है।

लेकिन configuration को पूरी तरह से configuration फ़ाइल में move करने की spirit में, आइए उस default value को workflow फ़ाइल से entirely remove करें।

#### 1.1.2. Workflow फ़ाइल में `batch` के लिए default value remove करो

`3-main.nf` workflow फ़ाइल में निम्नलिखित code change करो:

=== "बाद में"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "पहले"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

अब workflow फ़ाइल itself इन parameters के लिए कोई default values set नहीं करती।

#### 1.1.3. Pipeline चलाओ

आइए test करें कि यह command line में कोई parameters specify किए बिना correctly काम करता है।

```bash
nextflow run 3-main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

यह अभी भी पहले जैसा ही output produce करता है।

Final ASCII art output `results/3-main/` डायरेक्टरी में है, `cowpy-COLLECTED-batch-output.txt` नाम के अंतर्गत, पहले जैसा ही।

??? abstract "फ़ाइल सामग्री"

    ```console title="results/3-main/cowpy-COLLECTED-batch-output.txt"
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

Functionally, इस move ने कुछ नहीं बदला, लेकिन conceptually configuration फ़ाइल में default values set होना थोड़ा cleaner है।

### 1.2. Run-specific configuration फ़ाइल use करें

??? example "परिदृश्य"

    तुम अपनी main configuration फ़ाइल modify किए बिना different settings के साथ experiment करना चाहते हो।

तुम ऐसा एक subdirectory में एक नई `nextflow.config` फ़ाइल बनाकर कर सकते हो जिसे तुम अपने experiments के लिए working directory के रूप में use करोगे।

#### 1.2.1. Blank configuration के साथ working directory बनाओ

आइए एक नई डायरेक्टरी बनाकर और उसमें move करके शुरू करें:

```bash
mkdir -p tux-run
cd tux-run
```

फिर, उस डायरेक्टरी में एक blank configuration फ़ाइल बनाओ:

```bash
touch nextflow.config
```

यह एक empty फ़ाइल produce करती है।

#### 1.2.2. Experimental configuration set up करो

अब नई फ़ाइल open करो और जो parameters तुम customize करना चाहते हो वे add करो:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

ध्यान दो कि input फ़ाइल का path directory structure को reflect करना चाहिए।

#### 1.2.3. Pipeline चलाओ

अब हम अपनी नई working directory के भीतर से pipeline run कर सकते हैं।
Path को accordingly adapt करना सुनिश्चित करो!

```bash
nextflow run ../3-main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

यह `tux-run/` के अंतर्गत directories का एक नया set create करेगा जिसमें `tux-run/work/` और `tux-run/results/` शामिल हैं।

इस run में, Nextflow हमारी current डायरेक्टरी में `nextflow.config` को pipeline की root डायरेक्टरी में `nextflow.config` के साथ combine करता है, और इस तरह default character (turkey) को tux character से override करता है।

Final output फ़ाइल में greetings बोलते हुए tux character होना चाहिए।

??? abstract "फ़ाइल सामग्री"

    ```console title="tux-run/results/3-main/cowpy-COLLECTED-experiment-output.txt"
    _________
    / HELLO   \
    | BONJOUR |
    \ HOLà    /
    ---------
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

बस इतना ही; अब तुम्हारे पास अपनी 'normal' configuration modify किए बिना experimenting के लिए एक space है।

!!! warning "चेतावनी"

    अगले section में जाने से पहले पिछली डायरेक्टरी में वापस change करना सुनिश्चित करो!

    ```bash
    cd ..
    ```

अब आइए parameter values set करने का एक और useful तरीका देखें।

### 1.3. Parameter फ़ाइल use करें

??? example "परिदृश्य"

    तुम्हें exact run parameters किसी collaborator के साथ share करने होंगे, या उन्हें publication के लिए record करना होगा।

Subdirectory approach experimenting के लिए बढ़िया काम करता है, लेकिन इसमें थोड़ा setup involve है और requires कि तुम paths को accordingly adapt करो।
जब तुम अपनी pipeline को values के specific set के साथ run करना चाहते हो, या किसी और को minimal effort के साथ ऐसा करने में enable करना चाहते हो, तो एक simpler approach है।

Nextflow हमें [parameter file](https://nextflow.io/docs/latest/config.html#parameter-file) के माध्यम से YAML या JSON format में parameters specify करने की अनुमति देता है, जो default values के alternative sets manage और distribute करना बहुत convenient बनाता है, साथ ही run-specific parameter values भी।

#### 1.3.1. Example parameter फ़ाइल की जांच करो

इसे demonstrate करने के लिए, हम current डायरेक्टरी में एक example parameter फ़ाइल provide करते हैं, जिसे `test-params.yaml` कहा जाता है:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

इस parameter फ़ाइल में प्रत्येक input के लिए एक key-value pair है जिसे हम specify करना चाहते हैं।
ध्यान दो कि यदि तुम configuration फ़ाइल से syntax compare करो तो equal signs (`=`) के बजाय colons (`:`) का use है।
Config फ़ाइल Groovy में लिखी है, जबकि parameter फ़ाइल YAML में लिखी है।

!!! info "जानकारी"

    हम एक example के रूप में parameter फ़ाइल का JSON version भी provide करते हैं लेकिन हम यहां इसके साथ run नहीं करने जा रहे।
    उसे अपने आप try करने के लिए free feel करो।

#### 1.3.2. Pipeline चलाओ

इस parameter फ़ाइल के साथ workflow run करने के लिए, simply base कमांड में `-params-file <filename>` add करो।

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [2b/9a7d1e] sayHello (2)       | 3 of 3 ✔
    [5c/8f3b2a] convertToUpper (3) | 3 of 3 ✔
    [a3/29d8fb] collectGreetings   | 1 of 1 ✔
    [b7/83ef12] cowpy              | 1 of 1 ✔
    ```

Final output फ़ाइल में greetings बोलते हुए stegosaurus character होना चाहिए।

??? abstract "फ़ाइल सामग्री"

    ```console title="results/3-main/cowpy-COLLECTED-yaml-output.txt"
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
    \                             .       .
    \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
      \                 |    \  \ - ~ ~ - /  /    |
            _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
          ---------   \__/                        \__/
          .'  O    \     /               /       \  "
        (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Parameter फ़ाइल use करना overkill लग सकता है जब तुम्हारे पास specify करने के लिए केवल कुछ parameters हों, लेकिन कुछ pipelines दर्जनों parameters expect करती हैं।
उन cases में, parameter फ़ाइल use करना हमें massive command lines type किए बिना और workflow script modify किए बिना runtime पर parameter values provide करने की अनुमति देगा।

यह collaborators को parameters के sets distribute करना भी आसान बनाता है, या publication के लिए supporting information के रूप में, उदाहरण के लिए।
यह तुम्हारे काम को दूसरों द्वारा अधिक reproducible बनाता है।

### सीख

तुम जानते हो कि workflow inputs manage करने के लिए key configuration options का लाभ कैसे उठाएं।

### आगे क्या?

सीखो कि तुम्हारे workflow outputs कहां और कैसे publish होते हैं इसे कैसे manage करें।

---

## 2. Workflow outputs manage करें

??? example "परिदृश्य"

    तुम्हारी pipeline outputs को एक hardcoded डायरेक्टरी में publish करती है, लेकिन तुम हर बार workflow code edit किए बिना project या experiment name से results organize करना चाहते हो।

हमें जो workflow मिली है वह workflow-level output declarations के लिए paths use करती है, जो terribly flexible नहीं है और इसमें बहुत repetition involve है।

आइए कुछ common तरीके देखें जिनसे तुम इसे अधिक flexible होने के लिए configure कर सकते हो।

### 2.1. `outputDir` directory name customize करो

अब तक हमने जो workflow का प्रत्येक version run किया है उसने अपने outputs को output definitions में hardcoded एक different subdirectory में publish किया है।

हमने Part 1 में `-output-dir` CLI flag use करके उस subdirectory को change किया था, लेकिन वह अभी भी बस एक static string है।
आइए इसे एक config फ़ाइल में configure करें, जहां हम अधिक complex dynamic paths define कर सकते हैं।
हम इसके लिए एक whole new parameter create कर सकते थे, लेकिन आइए `batch` parameter use करें क्योंकि यह right there है।

#### 2.1.1. Configuration फ़ाइल में `outputDir` के लिए एक value set करो

Nextflow outputs publish करने के लिए जो path use करता है वह `outputDir` option द्वारा controlled है।
सभी outputs के लिए path change करने के लिए, तुम `nextflow.config` configuration फ़ाइल में इस option के लिए एक value set कर सकते हो।

`nextflow.config` फ़ाइल में निम्नलिखित code add करो:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

यह built-in default path, `results/`, को `results_config/` plus subdirectory के रूप में `batch` parameter की value से replace करेगा।

याद रहे कि तुम अपने command में `-output-dir` parameter use करके command-line से भी यह option set कर सकते हो (संक्षिप्त रूप में `-o`), लेकिन फिर तुम `batch` parameter value use नहीं कर पाओगे।
CLI flag use करने से config में `outputDir` overwrite हो जाएगा यदि वह set है।

#### 2.1.2. Hardcoded path का repeated part remove करो

हमारे पास output options में अभी भी एक subdirectory hardcoded है, तो आइए अब उससे छुटकारा पाएं।

Workflow फ़ाइल में निम्नलिखित code changes करो:

=== "बाद में"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "पहले"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path '3-main/intermediates'
            mode 'copy'
        }
        uppercased {
            path '3-main/intermediates'
            mode 'copy'
        }
        collected {
            path '3-main/intermediates'
            mode 'copy'
        }
        batch_report {
            path '3-main'
            mode 'copy'
        }
        cowpy_art {
            path '3-main'
            mode 'copy'
        }
    }
    ```

हम `outputDir` default modify करने के बजाय प्रत्येक path में बस `${params.batch}` भी add कर सकते थे, लेकिन यह अधिक concise है।

#### 2.1.3. Pipeline चलाओ

आइए test करें कि यह correctly काम करता है, command line से batch name को `outdir` set करते हुए।

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [amazing_church] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [9c/6a03ea] sayHello (2)       [100%] 3 of 3 ✔
    [11/9e58a6] convertToUpper (3) [100%] 3 of 3 ✔
    [c8/1977e5] collectGreetings   [100%] 1 of 1 ✔
    [38/f01eda] cowpy              [100%] 1 of 1 ✔
    ```

यह अभी भी पहले जैसा ही output produce करता है, सिवाय इसके कि इस बार हम अपने outputs `results_config/outdir/` के अंतर्गत पाते हैं।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results_config/outdir
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── outdir-report.txt
    ```

तुम जो भी directory hierarchy चाहो उसे construct करने के लिए इस approach को custom path definitions के साथ combine कर सकते हो।

### 2.2. Process के अनुसार outputs organize करो

Outputs को और organize करने का एक popular तरीका है इसे process के अनुसार करना, _i.e._ pipeline में run होने वाले प्रत्येक process के लिए subdirectories create करना।

#### 2.2.1. Output paths को process names के reference से replace करो

तुम्हें बस output path declaration में process के name को `<process>.name` के रूप में reference करना है।

Workflow फ़ाइल में निम्नलिखित changes करो:

=== "बाद में"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "पहले"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

यह output path configuration से remaining hardcoded elements remove करता है।

#### 2.2.2. Pipeline चलाओ

आइए test करें कि यह correctly काम करता है, command line से batch name को `pnames` set करते हुए।

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [4a/c2e6b8] sayHello (2)       | 3 of 3 ✔
    [6f/d4a172] convertToUpper (3) | 3 of 3 ✔
    [e8/4f19d7] collectGreetings   | 1 of 1 ✔
    [f2/a85c36] cowpy              | 1 of 1 ✔
    ```

यह अभी भी पहले जैसा ही output produce करता है, सिवाय इसके कि इस बार हम अपने outputs `results_config/pnames/` के अंतर्गत पाते हैं, और वे process के अनुसार grouped हैं।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results_config/pnames/
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

!!! note "नोट"

    ध्यान दो कि यहां हमने `intermediates` बनाम top level पर final outputs के बीच distinction erase कर दिया है।
    तुम बेशक इन approaches को mix और match कर सकते हो और multiple variables भी include कर सकते हो, उदाहरण के लिए पहले output का path `#!groovy "${params.batch}/intermediates/${sayHello.name}"` set करके।

### 2.3. Workflow level पर publish mode set करो

अंत में, repetitive code की मात्रा reduce करने की spirit में, हम per-output `mode` declarations को configuration में एक single पंक्ति से replace कर सकते हैं।

#### 2.3.1. Configuration फ़ाइल में `workflow.output.mode` add करो

`nextflow.config` फ़ाइल में निम्नलिखित code add करो:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

`outputDir` option की तरह ही, configuration फ़ाइल में `workflow.output.mode` को एक value देना workflow फ़ाइल में जो set है उसे override करने के लिए sufficient होगा, लेकिन आइए unnecessary code anyway remove करें।

#### 2.3.2. Workflow फ़ाइल से output mode remove करो

Workflow फ़ाइल में निम्नलिखित changes करो:

=== "बाद में"

    ```groovy title="3-main.nf" linenums="42"
    output {
        first_output {
            path { sayHello.name }
        }
        uppercased {
            path { convertToUpper.name }
        }
        collected {
            path { collectGreetings.name }
        }
        batch_report {
            path { collectGreetings.name }
        }
        cowpy_art {
            path { cowpy.name }
        }
    }
    ```

=== "पहले"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

यह अधिक concise है, है ना?

#### 2.3.3. Pipeline चलाओ

आइए test करें कि यह correctly काम करता है, command line से batch name को `outmode` set करते हुए।

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [5b/d91e3c] sayHello (2)       | 3 of 3 ✔
    [8a/f6c241] convertToUpper (3) | 3 of 3 ✔
    [89/cd3a48] collectGreetings   | 1 of 1 ✔
    [9e/71fb52] cowpy              | 1 of 1 ✔
    ```

यह अभी भी पहले जैसा ही output produce करता है, सिवाय इसके कि इस बार हम अपने outputs `results_config/outmode/` के अंतर्गत पाते हैं।
वे अभी भी सब proper copies हैं, symlinks नहीं।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results_config/outmode/
    ├── collectGreetings
    │   ├── COLLECTED-outmode-output.txt
    │   └── outmode-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-outmode-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Per-output तरीके से mode set करने का मुख्य कारण जो तुम अभी भी चाह सकते हो वह है यदि तुम same workflow के भीतर mix और match करना चाहते हो, _i.e._ कुछ outputs copied हों और कुछ symlinked।

बहुत सारे अन्य options हैं जिन्हें तुम इस तरह customize कर सकते हो, लेकिन hopefully यह तुम्हें options की range और उन्हें अपनी preferences के अनुसार effectively utilize करने का sense देता है।

### सीख

तुम जानते हो कि तुम्हारे outputs कहां publish होते हैं उन directories की naming और structure को कैसे control करें, साथ ही workflow output publishing mode।

### आगे क्या?

सीखो कि अपनी workflow configuration को अपने compute environment में कैसे adapt करें, software packaging technology से शुरू करते हुए।

---

## 3. Software packaging technology select करें

अब तक हम configuration elements देख रहे थे जो control करते हैं कि inputs कैसे जाते हैं और outputs कहां आते हैं। अब अपनी workflow configuration को अपने compute environment में adapt करने पर अधिक specifically focus करने का समय है।

उस path पर पहला step यह specify करना है कि प्रत्येक step में run होने वाले software packages कहां से आने वाले हैं।
क्या वे local compute environment में पहले से installed हैं?
क्या हमें images retrieve करके उन्हें container system के माध्यम से run करने की जरूरत है?
या क्या हमें Conda packages retrieve करके एक local Conda environment build करने की जरूरत है?

इस प्रशिक्षण course के पहले part (Parts 1-4) में हमने बस अपनी workflow में locally installed software use किया।
फिर Part 5 में, हमने Docker containers और `nextflow.config` फ़ाइल introduce की, जिसे हमने Docker containers का use enable करने के लिए use किया।

अब आइए देखें कि हम `nextflow.config` फ़ाइल के माध्यम से एक alternative software packaging option कैसे configure कर सकते हैं।

### 3.1. Config फ़ाइल में Docker disable करो और Conda enable करो

??? example "परिदृश्य"

    तुम अपनी pipeline को एक HPC cluster पर move कर रहे हो जहां security reasons के लिए Docker allowed नहीं है।
    Cluster Singularity और Conda support करता है, इसलिए तुम्हें अपनी configuration accordingly switch करनी होगी।

जैसा कि पहले noted था, Nextflow Singularity (जो HPC पर अधिक widely use होता है) सहित multiple container technologies support करता है, साथ ही software package managers जैसे Conda भी।

हम Docker के बजाय Conda use करने के लिए अपनी configuration फ़ाइल change कर सकते हैं।
ऐसा करने के लिए, आइए `docker.enabled` की value को `false` में switch करें, और Conda का use enable करने वाला एक directive add करें:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

यह Nextflow को उन processes के लिए Conda environments create और utilize करने की अनुमति देगा जिनके पास Conda packages specified हैं।
जिसका मतलब है कि अब हमें अपने `cowpy` process में उनमें से एक add करना होगा!

### 3.2. Process definition में एक Conda package specify करो

हमने `cowpy` tool containing Conda package के लिए URI पहले ही retrieve कर लिया है: `conda-forge::cowpy==1.1.5`

अब हम `conda` directive use करके `cowpy` process definition में URI add करते हैं:

=== "बाद में"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "पहले"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

Clear करने के लिए, हम `docker` directive को _replace_ नहीं कर रहे, हम एक alternative option _add_ कर रहे हैं।

!!! tip "सुझाव"

    किसी given conda package के लिए URI प्राप्त करने के कुछ different तरीके हैं।
    हम [Seqera Containers](https://seqera.io/containers/) search query use करने की recommend करते हैं, जो तुम्हें एक URI देगा जिसे तुम copy और paste कर सकते हो, भले ही तुम इससे container create करने का plan नहीं कर रहे हो।

### 3.3. Verify करने के लिए workflow run करो कि यह Conda use कर सकता है

आइए try करें।

```bash
nextflow run 3-main.nf --batch conda
```

??? success "कमांड आउटपुट"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

यह बिना किसी issue के काम करना चाहिए और `results_config/conda` के अंतर्गत पहले जैसे ही outputs produce करना चाहिए।

Behind the scenes, Nextflow ने Conda packages retrieve किए हैं और environment create किया है, जिसमें normally थोड़ा काम लगता है; तो यह अच्छा है कि हमें वह सब खुद नहीं करना पड़ा!

!!! info "जानकारी"

    यह quickly run होता है क्योंकि `cowpy` package quite small है, लेकिन यदि तुम large packages के साथ काम कर रहे हो, पहली बार इसमें usual से थोड़ा अधिक समय लग सकता है, और तुम complete होने से पहले console output को एक मिनट या उससे अधिक के लिए 'stuck' देख सकते हो।
    यह normal है और पहली बार नया package use करने पर Nextflow जो extra काम करता है उसके कारण है।

हमारे standpoint से, यह exactly same दिखता है जैसे Docker के साथ run करना, भले ही backend पर mechanics थोड़ी different हैं।

इसका मतलब है कि हम जरूरत पड़ने पर Conda environments के साथ run करने के लिए all set हैं।

??? info "Docker और Conda को mix और match करना"

    चूंकि ये directives प्रति process assign की जाती हैं, 'mix और match' करना possible है, _i.e._ अपनी workflow में कुछ processes को Docker के साथ और अन्य को Conda के साथ run करने के लिए configure करना, उदाहरण के लिए, यदि तुम जो compute infrastructure use कर रहे हो वह दोनों support करती है।
    उस case में, तुम अपनी configuration फ़ाइल में Docker और Conda दोनों enable करोगे।
    यदि किसी given process के लिए दोनों available हैं, Nextflow containers को prioritize करेगा।

    और जैसा कि पहले noted था, Nextflow multiple other software packaging और container technologies support करता है, इसलिए तुम just उन दो तक limited नहीं हो।

### सीख

तुम जानते हो कि प्रत्येक process को कौन सा software package use करना चाहिए configure कैसे करें, और technologies के बीच switch कैसे करें।

### आगे क्या?

सीखो कि Nextflow द्वारा actually काम करने के लिए use किया जाने वाला execution platform कैसे change करें।

---

## 4. Execution platform select करें

??? example "परिदृश्य"

    तुम अपने laptop पर अपनी pipeline develop और test कर रहे थे, लेकिन अब तुम्हें इसे हजारों samples पर run करना है।
    तुम्हारी institution के पास एक Slurm scheduler वाला HPC cluster है जिसे तुम इसके बजाय use करना चाहोगे।

अब तक, हम अपनी pipeline local executor के साथ run कर रहे थे।
यह प्रत्येक task को उस machine पर execute करता है जिस पर Nextflow run हो रहा है।
जब Nextflow begin होता है, यह available CPUs और memory देखता है।
यदि run करने के लिए ready tasks के resources available resources से exceed करते हैं, Nextflow last tasks को execution से hold back रखेगा जब तक कि एक या अधिक earlier tasks finish नहीं हो जाते, necessary resources free करते हुए।

Local executor convenient और efficient है, लेकिन यह उस single machine तक limited है। बहुत large workloads के लिए, तुम discover कर सकते हो कि तुम्हारी local machine एक bottleneck है, या तो इसलिए कि तुम्हारे पास एक single task है जिसे तुम्हारे पास available से अधिक resources चाहिए, या इसलिए कि तुम्हारे पास इतने सारे tasks हैं कि single machine के उन्हें run करने का wait करना बहुत long लेगा।

Nextflow [कई different execution backends](https://nextflow.io/docs/latest/executor.html) support करता है, जिसमें HPC schedulers (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor और अन्य) के साथ-साथ cloud execution backends जैसे (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes और अधिक) शामिल हैं।

### 4.1. Different backend target करना

Executor का choice `executor` नाम की एक process directive द्वारा set होता है।
Default रूप से यह `local` पर set है, इसलिए following configuration implied है:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Different backend target करने के लिए executor set करने के लिए, तुम simply वह executor specify करोगे जो तुम चाहते हो resource allocations के लिए ऊपर described similar syntax use करके (सभी options के लिए [Executors](https://nextflow.io/docs/latest/executor.html) देखो)।

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "चेतावनी"

    हम actually training environment में इसे test नहीं कर सकते क्योंकि यह HPC से connect करने के लिए set up नहीं है।

### 4.2. Execution parameters के लिए backend-specific syntax से deal करना

अधिकांश high-performance computing platforms allow करते हैं (और कभी-कभी require करते हैं) कि तुम resource allocation requests और limitations (जैसे number of CPUs और memory) और use करने के लिए job queue का name जैसे certain parameters specify करो।

दुर्भाग्य से, इनमें से प्रत्येक system different technologies, syntaxes और configurations use करती है यह define करने के लिए कि job को कैसे define और relevant scheduler को submit किया जाना चाहिए।

??? abstract "उदाहरण"

    उदाहरण के लिए, same job जिसे 8 CPUs और 4GB RAM की आवश्यकता है "my-science-work" queue पर execute होने के लिए backend के आधार पर following different तरीकों से express करने की जरूरत है।

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Fortunately, Nextflow यह सब simplify करता है।
यह एक standardized syntax provide करता है ताकि तुम `cpus`, `memory` और `queue` (सभी available options के लिए [Process directives](https://nextflow.io/docs/latest/reference/process.html#process-directives) देखो) जैसी relevant properties बस एक बार specify कर सको।
फिर, runtime पर, Nextflow executor setting के आधार पर appropriate backend-specific scripts generate करने के लिए उन settings use करेगा।

हम अगले section में उस standardized syntax को cover करेंगे।

### सीख

तुम अब जानते हो कि different kinds की computing infrastructure use करने के लिए executor कैसे change करें।

### आगे क्या?

सीखो कि Nextflow में resource allocations और limitations कैसे evaluate और express करें।

---

## 5. Compute resource allocations control करें

??? example "परिदृश्य"

    तुम्हारी pipeline cluster पर keep failing हो रही है क्योंकि tasks memory limits exceed करने के लिए kill हो रहे हैं।
    या शायद तुम्हें ऐसे resources के लिए charge किया जा रहा है जिनका तुम use नहीं कर रहे और costs optimize करना चाहते हो।

अधिकांश high-performance computing platforms allow करते हैं (और कभी-कभी require करते हैं) कि तुम certain resource allocation parameters जैसे number of CPUs और memory specify करो।

Default रूप से, Nextflow प्रत्येक process के लिए single CPU और 2GB memory use करेगा।
Corresponding process directives `cpus` और `memory` कहलाती हैं, इसलिए following configuration implied है:

```groovy title="Built-in configuration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

तुम इन values को modify कर सकते हो, या तो सभी processes के लिए या specific named processes के लिए, अपनी configuration फ़ाइल में additional process directives use करके।
Nextflow उन्हें chosen executor के लिए appropriate instructions में translate करेगा।

लेकिन तुम कैसे जानोगे कि कौन सी values use करनी हैं?

### 5.1. Resource utilization report generate करने के लिए workflow run करो

??? example "परिदृश्य"

    तुम नहीं जानते कि तुम्हारे processes को कितनी memory या CPU चाहिए और resources waste करने या jobs kill होने से बचना चाहते हो।

यदि तुम up front नहीं जानते कि तुम्हारे processes को कितनी CPU और memory likely चाहिए होगी, तुम कुछ resource profiling कर सकते हो, meaning तुम workflow को कुछ default allocations के साथ run करते हो, record करते हो कि प्रत्येक process ने कितना use किया, और वहां से, estimate करते हो कि base allocations कैसे adjust करनी हैं।

Conveniently, Nextflow में ऐसा करने के लिए built-in tools शामिल हैं, और request पर तुम्हारे लिए happily एक report generate करेगा।

ऐसा करने के लिए, अपनी command line में `-with-report <filename>.html` add करो।

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

Report एक html फ़ाइल है, जिसे तुम download करके अपने browser में open कर सकते हो। तुम training environment में इसे view करने के लिए left पर file explorer में इस पर right click करके `Show preview` पर click भी कर सकते हो।

Report को देखने में कुछ minutes लो और देखो कि क्या तुम resources adjust करने के लिए कुछ opportunities identify कर सकते हो।
जो tabs utilization results को allocated के percentage के रूप में show करते हैं उन पर click करना सुनिश्चित करो।

सभी available features describe करने वाला [Reports](https://nextflow.io/docs/latest/reports.html) documentation देखो।

### 5.2. सभी processes के लिए resource allocations set करो

Profiling show करती है कि हमारी training workflow में processes बहुत lightweight हैं, इसलिए आइए default memory allocation को प्रति process 1GB तक reduce करें।

अपनी `nextflow.config` फ़ाइल में, pipeline parameters section से पहले, निम्नलिखित add करो:

```groovy title="nextflow.config" linenums="4"
/*
* Process settings
*/
process {
    memory = 1.GB
}
```

यह हमारे द्वारा consume किए जाने वाले compute की मात्रा reduce करने में help करेगा।

### 5.3. Specific process के लिए resource allocations set करो

साथ ही, हम pretend करने जा रहे हैं कि `cowpy` process को others से अधिक resources चाहिए, बस इसलिए कि हम demonstrate कर सकें कि individual process के लिए allocations कैसे adjust करें।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }
    ```

इस configuration के साथ, सभी processes 1GB memory और single CPU (implied default) request करेंगे, सिवाय `cowpy` process के, जो 2GB और 2 CPUs request करेगा।

!!! info "जानकारी"

    यदि तुम्हारे पास कम CPUs वाली machine है और तुम प्रति process high number allocate करते हो, तुम process calls को एक दूसरे के पीछे queue होते देख सकते हो।
    ऐसा इसलिए है क्योंकि Nextflow ensure करता है कि हम available से अधिक CPUs request न करें।

### 5.4. Updated configuration के साथ workflow run करो

आइए try करें, profiling report के लिए एक different filename supply करते हुए ताकि हम configuration changes से पहले और बाद की performance compare कर सकें।

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

तुम शायद कोई real difference notice नहीं करोगे क्योंकि यह इतना small workload है, लेकिन यह वह approach है जिसे तुम real-world workflow की performance और resource requirements analyze करने के लिए use करोगे।

यह बहुत useful है जब तुम्हारे processes की different resource requirements हैं। यह तुम्हें guesswork नहीं, actual data के आधार पर प्रत्येक process के लिए सेट की गई resource allocations को right-size करने का power देता है।

!!! tip "सुझाव"

    यह सिर्फ एक tiny taster है जो तुम resources के अपने use को optimize करने के लिए कर सकते हो।
    Nextflow में itself कुछ really neat [dynamic retry logic](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) built in है जो resource limitations के कारण fail होने वाले jobs को retry करता है।
    Additionally, Seqera Platform आपकी resource allocations को automatically optimize करने के लिए AI-driven tooling भी offer करता है।

### 5.5. Resource limits add करो

तुम जो computing executor और compute infrastructure use कर रहे हो उसके आधार पर, कुछ constraints हो सकती हैं कि तुम क्या (या must) allocate कर सकते हो।
उदाहरण के लिए, तुम्हारे cluster को require हो सकता है कि तुम certain limits के भीतर रहो।

तुम `resourceLimits` directive use करके relevant limitations set कर सकते हो। Syntax इस तरह दिखता है जब यह process block में अकेला है:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow तुम्हारे द्वारा specified executor के आधार पर इन values को appropriate instructions में translate करेगा।

हम इसे run नहीं करने जा रहे, क्योंकि training environment में हमारे पास relevant infrastructure तक access नहीं है।
हालांकि, यदि तुम resource allocations के साथ workflow run करने की try करते जो इन limits से exceed करती हैं, फिर `.command.run` script फ़ाइल में `sbatch` command look up करते, तुम देखोगे कि executor को actually भेजे जाने वाले requests `resourceLimits` द्वारा specified values पर capped हैं।

??? info "Institutional reference configurations"

    nf-core project ने world भर में various institutions द्वारा shared [configuration files का collection](https://nf-co.re/configs/) compile किया है, जो HPC और cloud executors की एक wide range cover करती है।

    वे shared configs दोनों के लिए valuable हैं जो लोग वहां काम करते हैं और इसलिए बस अपनी institution की configuration out of the box utilize कर सकते हैं, और एक model के रूप में उन लोगों के लिए जो अपनी own infrastructure के लिए configuration develop करना चाहते हैं।

### सीख

तुम जानते हो कि resource utilization assess करने के लिए profiling report कैसे generate करें और सभी processes और/या individual processes के लिए resource allocations कैसे modify करें, साथ ही HPC पर run करने के लिए resource limitations set करें।

### आगे क्या?

सीखो कि preset configuration profiles कैसे set up करें और runtime पर उनके बीच switch करें।

---

## 6. Preset configurations के बीच switch करने के लिए profiles use करें

??? example "परिदृश्य"

    तुम regularly अपने laptop पर development के लिए और अपनी institution के HPC पर production runs के लिए pipelines run करने के बीच switch करते हो।
    तुम हर बार environments switch करने पर manually configuration settings change करने से थक गए हो।

हमने तुम्हें कई तरीके दिखाए हैं जिनसे तुम अपनी pipeline configuration को उस project के आधार पर customize कर सकते हो जिस पर तुम काम कर रहे हो या जो compute environment तुम use कर रहे हो।

तुम शायद जो computing infrastructure use कर रहे हो उसके आधार पर alternative settings के बीच switch करना चाहो। उदाहरण के लिए, तुम अपने laptop पर locally develop और small-scale tests run करना चाहोगे, फिर HPC या cloud पर full-scale workloads run करना।

Nextflow तुम्हें any number of [**profiles**](https://nextflow.io/docs/latest/config.html#profiles) set up करने देता है जो different configurations describe करते हैं, जिन्हें तुम फिर configuration फ़ाइल itself modify करने के बजाय एक command-line argument use करके runtime पर select कर सकते हो।

### 6.1. Local development और HPC पर execution के बीच switch करने के लिए profiles create करो

आइए दो alternative profiles set up करें; एक regular computer पर small scale loads run करने के लिए, जहां हम Docker containers use करेंगे, और एक Slurm scheduler वाले university HPC पर run करने के लिए, जहां हम Conda packages use करेंगे।

#### 6.1.1. Profiles set up करो

अपनी `nextflow.config` फ़ाइल में, pipeline parameters section के बाद लेकिन output settings से पहले, निम्नलिखित add करो:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
}
```

तुम देखते हो कि university HPC के लिए, हम resource limitations भी specify कर रहे हैं।

#### 6.1.2. Profile के साथ workflow run करो

अपनी Nextflow command line में profile specify करने के लिए, हम `-profile` argument use करते हैं।

आइए `my_laptop` configuration के साथ workflow run करने की try करें।

```bash
nextflow run 3-main.nf -profile my_laptop
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

जैसा कि तुम देख सकते हो, यह हमें runtime पर configurations के बीच बहुत conveniently toggle करने की अनुमति देता है।

!!! warning "चेतावनी"

    `univ_hpc` profile training environment में properly run नहीं होगी क्योंकि हमारे पास Slurm scheduler तक access नहीं है।

यदि future में हम configuration के अन्य elements पाते हैं जो always इनके साथ co-occurring हैं, हम simply उन्हें corresponding profile(s) में add कर सकते हैं।
हम additional profiles भी create कर सकते हैं यदि configuration के अन्य elements हैं जिन्हें हम एक साथ group करना चाहते हैं।

### 6.2. Test parameters का एक profile create करो

??? example "परिदृश्य"

    तुम चाहते हो कि others अपना own input data gather किए बिना quickly तुम्हारी pipeline try कर सकें।

Profiles सिर्फ infrastructure configuration के लिए नहीं हैं।
हम उन्हें workflow parameters के लिए default values set करने के लिए भी use कर सकते हैं, ताकि others के लिए workflow try करना आसान हो बिना उन्हें appropriate input values खुद gather करने के।
तुम इसे parameter फ़ाइल use करने का एक alternative consider कर सकते हो।

#### 6.2.1. Profile set up करो

इस context में default values express करने का syntax इस तरह दिखता है, एक profile के लिए जिसे हम `test` name देते हैं:

```groovy title="Syntax example"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

यदि हम अपनी workflow के लिए एक test profile add करते हैं, `profiles` block इस तरह बन जाता है:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

Technical configuration profiles की तरह ही, तुम जो भी arbitrary name चाहो उसके अंतर्गत parameters specify करने वाले multiple different profiles set up कर सकते हो।

#### 6.2.2. Test profile के साथ locally workflow run करो

Conveniently, profiles mutually exclusive नहीं हैं, इसलिए हम अपनी command line में following syntax `-profile <profile1>,<profile2>` (any number of profiles के लिए) use करके multiple profiles specify कर सकते हैं।

यदि तुम ऐसे profiles combine करते हो जो configuration के same elements के लिए values set करते हैं और same configuration फ़ाइल में described हैं, Nextflow conflict resolve करेगा जो भी value उसने last में read की (_i.e._ जो भी फ़ाइल में बाद में आती है) उसे use करके।
यदि conflicting settings different configuration sources में set हैं, default [order of precedence](https://www.nextflow.io/docs/latest/config.html) apply होता है।

आइए अपने previous command में test profile add करने की try करें:

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

यह जहां possible हो Docker use करेगा और `results_config/test` के अंतर्गत outputs produce करेगा, और इस बार character comedic duo `dragonandcow` है।

??? abstract "फ़ाइल सामग्री"

    ```console title="results_config/test/"
     _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
                \                    ^    /^
                  \                  / \  // \
                  \   |\___/|      /   \//  .\
                    \  /O  O  \__  /    //  | \ \           *----*
                      /     /  \/_/    //   |  \  \          \   |
                      \@___\@`    \/_   //    |   \   \         \/\ \
                    0/0/|       \/_ //     |    \    \         \ \
                0/0/0/0/|        \///      |     \     \       | |
              0/0/0/0/0/_|_ /   (  //       |      \     _\     |  /
          0/0/0/0/0/0/`/,_ _ _/  ) ; -.    |    _ _\.-~       /   /
                      ,-}        _      *-.|.-~-.           .~    ~
      \     \__/        `/\      /                 ~-. _ .-~      /
      \____(oo)           *.   }            {                   /
      (    (--)          .----~-.\        \-`                 .~
      //__\\  \__ Ack!   ///.----..<        \             _ -~
      //    \\               ///-._ _ _ _ _ _ _{^ - - - - ~
    ```

इसका मतलब है कि जब तक हम workflow code के साथ कोई test data files distribute करते हैं, कोई भी command line या parameter फ़ाइल के माध्यम से अपने own inputs supply किए बिना quickly workflow try कर सकता है।

!!! tip "सुझाव"

    हम externally stored larger files के लिए URLs point कर सकते हैं।
    Nextflow उन्हें automatically download करेगा जब तक open connection है।

    अधिक details के लिए, Side Quest [Working with Files](../side_quests/working_with_files.md) देखो।

### 6.3. Resolved configuration देखने के लिए `nextflow config` use करो

जैसा कि ऊपर noted था, कभी-कभी same parameter को different values पर profiles में set किया जा सकता है जिन्हें तुम combine करना चाहते हो।
और अधिक generally, numerous places हैं जहां configuration के elements stored हो सकते हैं, और कभी-कभी same properties को different places में different values पर set किया जा सकता है।

Nextflow किसी भी conflicts को resolve करने के लिए set [order of precedence](https://www.nextflow.io/docs/latest/config.html#configuration-file) apply करता है, लेकिन यह खुद determine करना tricky हो सकता है।
और भले ही कुछ भी conflicting न हो, सभी possible places look up करना tedious हो सकता है जहां चीजें configured हो सकती हैं।

Fortunately, Nextflow में `config` नाम का एक convenient utility tool शामिल है जो तुम्हारे लिए उस whole process को automate कर सकता है।

`config` tool तुम्हारी current working directory में सभी contents explore करेगा, किसी भी configuration files को hoover up करेगा, और fully resolved configuration produce करेगा जो Nextflow workflow run करने के लिए use करेगा।
यह तुम्हें कुछ भी launch किए बिना find out करने की अनुमति देता है कि कौन सी settings use होंगी।

#### 6.3.1. Default configuration resolve करो

Default रूप से apply होने वाली configuration resolve करने के लिए यह कमांड run करो।

```bash
nextflow config
```

??? success "कमांड आउटपुट"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }

    docker {
      enabled = false
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
    }

    outputDir = 'results_config/batch'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

यह तुम्हें base configuration दिखाता है जो तुम्हें मिलती है यदि तुम command line में कुछ भी extra specify नहीं करते।

#### 6.3.2. Specific settings activated के साथ configuration resolve करो

यदि तुम command-line parameters provide करते हो, जैसे एक या अधिक profiles enable करना या parameter फ़ाइल load करना, command additionally उन्हें account में लेगा।

```bash
nextflow config -profile my_laptop,test
```

??? success "कमांड आउटपुट"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }

    docker {
      enabled = true
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
      executor = 'local'
    }

    outputDir = 'results_config/test'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

यह especially उन complex projects के लिए useful हो जाता है जिनमें configuration की multiple layers involve हैं।

### सीख

तुम जानते हो कि minimal hassle के साथ runtime पर preset configuration select करने के लिए profiles का use कैसे करें।
अधिक generally, तुम जानते हो कि different compute platforms suit करने और अपने analyses की reproducibility enhance करने के लिए अपनी workflow executions को कैसे configure करें।

### आगे क्या?

सीखो कि GitHub जैसे remote repositories से directly pipelines कैसे run करें।

---

## 7. Remote repositories से pipelines run करें

??? example "परिदृश्य"

    तुम एक well-established pipeline जैसे nf-core से run करना चाहते हो बिना code खुद download और manage किए।

अब तक हम current directory में located workflow scripts run कर रहे थे।
Practice में, तुम अक्सर remote repositories में stored pipelines run करना चाहोगे, जैसे GitHub।

Nextflow इसे straightforward बनाता है: तुम किसी भी pipeline को directly एक Git repository URL से run कर सकते हो बिना manually पहले इसे download किए।

### 7.1. GitHub से pipeline run करो

Remote pipeline run करने का basic syntax `nextflow run <repository>` है, जहां `<repository>` एक GitHub repository path हो सकता है जैसे `nextflow-io/hello`, एक full URL, या GitLab, Bitbucket, या other Git hosting services का path।

Official Nextflow "hello" demo pipeline run करने की try करो:

```bash
nextflow run nextflow-io/hello
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] DSL2 - revision: 2ce0b0e294 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

पहली बार जब तुम एक remote pipeline run करते हो, Nextflow इसे download करता है और locally cache करता है।
Subsequent runs cached version use करते हैं जब तक तुम explicitly एक update request नहीं करते।

### 7.2. Reproducibility के लिए version specify करो

Default रूप से, Nextflow default branch से latest version run करता है।
तुम `-r` flag use करके particular version (tag), branch, या commit specify कर सकते हो:

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `https://github.com/nextflow-io/hello` [sick_carson] DSL2 - revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Exact versions specify करना reproducibility के लिए essential है।

### सीख

तुम जानते हो कि GitHub और other remote repositories से directly pipelines कैसे run करें, और reproducibility के लिए versions कैसे specify करें।

### आगे क्या?

अपने आप को एक big pat on the back दो!
तुम सब कुछ जानते हो जो तुम्हें Nextflow pipelines run और manage करना शुरू करने के लिए जानना जरूरी है।

यह इस course को conclude करता है, लेकिन यदि तुम सीखना जारी रखने के लिए eager हो, हमारे पास दो main recommendations हैं:

- यदि तुम अपनी own pipelines develop करने में deeper dig करना चाहते हो, [Hello Nextflow](../hello_nextflow/index.md) देखो, beginners के लिए एक course जो इस वाले जैसी ही general progression cover करता है लेकिन channels और operators के बारे में बहुत अधिक detail में जाता है।
- यदि तुम code में deeper जाए बिना Nextflow pipelines run करना सीखना जारी रखना चाहते हो, [Hello nf-core](../hello_nf-core/index.md) के first part को देखो, जो hugely popular [nf-core](https://nf-co.re/) project से pipelines खोजने और run करने के लिए tooling introduce करता है।

मज़े करो!

---

## Quiz

<quiz>
जब parameter values workflow फ़ाइल और `nextflow.config` दोनों में set हैं, तो कौन precedence लेता है?
- [ ] Workflow फ़ाइल value
- [x] Configuration फ़ाइल value
- [ ] पहली encounter की गई value
- [ ] यह error cause करता है

अधिक जानें: [1.1. `nextflow.config` में values set up करें](#11-nextflowconfig-में-values-set-up-करें)
</quiz>

<quiz>
Workflow फ़ाइल में parameter default set करने बनाम config फ़ाइल में syntax difference क्या है?
- [ ] वे same syntax use करते हैं
- [x] Workflow typed declaration (`#!groovy param: Type = value`) use करती है, config assignment (`#!groovy param = value`) use करती है
- [ ] Config typed declaration use करती है, workflow assignment use करती है
- [ ] केवल config files default values set कर सकती हैं

अधिक जानें: [1.1. `nextflow.config` में values set up करें](#11-nextflowconfig-में-values-set-up-करें)
</quiz>

<quiz>
Workflow run करते समय parameter फ़ाइल कैसे specify करते हो?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

अधिक जानें: [1.3. Parameter फ़ाइल use करें](#13-parameter-फ़ाइल-use-करें)
</quiz>

<quiz>
`outputDir` configuration option क्या control करता है?
- [ ] Work directory का location
- [x] Base path जहां workflow outputs publish होते हैं
- [ ] Log files के लिए directory
- [ ] Module files का location

अधिक जानें: [2.1. `outputDir` directory name customize करो](#21-outputdir-directory-name-customize-करो)
</quiz>

<quiz>
Output path configuration में process name को dynamically कैसे reference करते हो?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<process>.name"`
- [x] `#!groovy path { <process>.name }`
- [ ] `@processName`

अधिक जानें: [2.2. Process के अनुसार outputs organize करो](#22-process-के-अनुसार-outputs-organize-करो)
</quiz>

<quiz>
यदि Docker और Conda दोनों enable हैं और process में दोनों directives हैं, तो कौन prioritized है?
- [x] Docker (containers)
- [ ] Conda
- [ ] Process में पहले defined जो
- [ ] यह error cause करता है

अधिक जानें: [3. Software packaging technology select करें](#3-software-packaging-technology-select-करें)
</quiz>

<quiz>
Nextflow में default executor क्या है?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

अधिक जानें: [4. Execution platform select करें](#4-execution-platform-select-करें)
</quiz>

<quiz>
Resource utilization report generate करने वाली कमांड क्या है?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

अधिक जानें: [5.1. Resource utilization report generate करने के लिए workflow run करो](#51-resource-utilization-report-generate-करने-के-लिए-workflow-run-करो)
</quiz>

<quiz>
Config फ़ाइल में `cowpy` नाम के specific process के लिए resource requirements कैसे set करते हो?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

अधिक जानें: [5.3. Specific process के लिए resource allocations set करो](#53-specific-process-के-लिए-resource-allocations-set-करो)
</quiz>

<quiz>
`resourceLimits` directive क्या करता है?
- [ ] Minimum resource requirements set करता है
- [ ] Processes को resources allocate करता है
- [x] Maximum resources cap करता है जो request की जा सकती हैं
- [ ] Real-time में resource usage monitor करता है

अधिक जानें: [5.5. Resource limits add करो](#55-resource-limits-add-करो)
</quiz>

<quiz>
एक single command में multiple profiles कैसे specify करते हो?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

अधिक जानें: [6. Preset configurations के बीच switch करने के लिए profiles use करें](#6-preset-configurations-के-बीच-switch-करने-के-लिए-profiles-use-करें)
</quiz>

<quiz>
वह command जो fully resolved configuration show करती है जो Nextflow use करेगा?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

अधिक जानें: [6.3. Resolved configuration देखने के लिए `nextflow config` use करो](#63-resolved-configuration-देखने-के-लिए-nextflow-config-use-करो)
</quiz>

<quiz>
Profiles किसके लिए use की जा सकती हैं? (सभी लागू select करो)
- [x] Infrastructure-specific settings define करना (executors, containers)
- [x] Different environments के लिए resource limits set करना
- [x] Easy workflow testing के लिए test parameters provide करना
- [ ] नई processes define करना

अधिक जानें: [6. Preset configurations के बीच switch करने के लिए profiles use करें](#6-preset-configurations-के-बीच-switch-करने-के-लिए-profiles-use-करें)
</quiz>
