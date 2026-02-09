# भाग 2: nf-core के लिए Hello को फिर से लिखना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [और जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core प्रशिक्षण पाठ्यक्रम के इस दूसरे भाग में, हम तुम्हें दिखाएंगे कि [Hello Nextflow](../hello_nextflow/index.md) शुरुआती पाठ्यक्रम द्वारा बनाई गई पाइपलाइन का nf-core संगत संस्करण कैसे बनाया जाए।

तुमने प्रशिक्षण के पहले खंड में देखा होगा कि nf-core पाइपलाइन काफी विस्तृत संरचना का पालन करती हैं जिसमें बहुत सारी सहायक फ़ाइलें होती हैं।
यह सब शुरू से बनाना बहुत थकाऊ होगा, इसलिए nf-core समुदाय ने इसके बजाय एक टेम्पलेट से यह करने के लिए टूलिंग विकसित की है, ताकि प्रक्रिया को बूटस्ट्रैप किया जा सके।

हम तुम्हें दिखाएंगे कि इस टूलिंग का उपयोग करके पाइपलाइन स्कैफोल्ड कैसे बनाया जाए, फिर मौजूदा 'नियमित' पाइपलाइन कोड को nf-core स्कैफोल्ड पर कैसे अनुकूलित किया जाए।

यदि तुम Hello पाइपलाइन से परिचित नहीं हो या तुम्हें याद दिलाने की जरूरत है, तो [यह जानकारी पृष्ठ](../info/hello_pipeline.md) देखो।

---

## 1. एक नई पाइपलाइन प्रोजेक्ट बनाना

सबसे पहले, हम नई पाइपलाइन के लिए स्कैफोल्ड बनाते हैं।

!!! note "नोट"

    सुनिश्चित करो कि तुम अपने टर्मिनल में `hello-nf-core` डायरेक्टरी में हो।

### 1.1. टेम्पलेट-आधारित पाइपलाइन निर्माण टूल चलाना

चलो `nf-core pipelines create` कमांड के साथ एक नई पाइपलाइन बनाकर शुरू करते हैं।
यह nf-core बेस टेम्पलेट का उपयोग करके एक नई पाइपलाइन स्कैफोल्ड बनाएगा, जिसे पाइपलाइन नाम, विवरण और लेखक के साथ अनुकूलित किया जाएगा।

```bash
nf-core pipelines create
```

इस कमांड को चलाने से पाइपलाइन निर्माण के लिए एक Text User Interface (TUI) खुलेगा:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/VwjXNXONHlY?si=d0HkFSISnKn76TeI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

यह TUI तुमसे अपनी पाइपलाइन के बारे में बुनियादी जानकारी प्रदान करने के लिए कहेगा और तुम्हें अपनी पाइपलाइन स्कैफोल्ड में शामिल करने या बाहर करने के लिए फीचर्स की पसंद प्रदान करेगा।

- स्वागत स्क्रीन पर, **Let's go!** पर क्लिक करो।
- `Choose pipeline type` स्क्रीन पर, **Custom** पर क्लिक करो।
- अपनी पाइपलाइन विवरण निम्नानुसार दर्ज करो (`< YOUR NAME >` को अपने नाम से बदलते हुए), फिर **Next** पर क्लिक करो।

```
[ ] GitHub organisation: core
[ ] Workflow name: hello
[ ] A short description of your pipeline: A basic nf-core style version of Hello Nextflow
[ ] Name of the main author(s): < YOUR NAME >
```

- Template features स्क्रीन पर, `Toggle all features` को **off** पर सेट करो, फिर चुनिंदा रूप से निम्नलिखित को **enable** करो। अपने चयन की जांच करो और **Continue** पर क्लिक करो।

```
[ ] Add testing profiles
[ ] Use nf-core components
[ ] Use nf-schema
[ ] Add configuration files
[ ] Add documentation
```

- `Final details` स्क्रीन पर, **Finish** पर क्लिक करो। पाइपलाइन बनने की प्रतीक्षा करो, फिर **Continue** पर क्लिक करो।
- Create GitHub repository स्क्रीन पर, **Finish without creating a repo** पर क्लिक करो। यह बाद में GitHub रिपॉजिटरी बनाने के लिए निर्देश प्रदर्शित करेगा। इन्हें अनदेखा करो और **Close** पर क्लिक करो।

एक बार TUI बंद हो जाने पर, तुम्हें निम्नलिखित कंसोल आउटपुट दिखाई देना चाहिए।

??? success "कमांड आउटपुट"

    ```console
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    INFO     Launching interactive nf-core pipeline creation tool.
    ```

कंसोल आउटपुट में कोई स्पष्ट पुष्टि नहीं है कि पाइपलाइन निर्माण काम कर गया, लेकिन तुम्हें `core-hello` नाम की एक नई डायरेक्टरी दिखाई देनी चाहिए।

नई डायरेक्टरी की सामग्री देखो कि तुमने टेम्पलेट का उपयोग करके कितना काम बचाया।

```bash
tree core-hello
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    core-hello/
    ├── assets
    │   ├── samplesheet.csv
    │   └── schema_input.json
    ├── conf
    │   ├── base.config
    │   ├── modules.config
    │   ├── test.config
    │   └── test_full.config
    ├── docs
    │   ├── output.md
    │   ├── README.md
    │   └── usage.md
    ├── main.nf
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── README.md
    ├── subworkflows
    │   ├── local
    │   │   └── utils_nfcore_hello_pipeline
    │   │       └── main.nf
    │   └── nf-core
    │       ├── utils_nextflow_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       └── nextflow.config
    │       ├── utils_nfcore_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       ├── main.workflow.nf.test.snap
    │       │       └── nextflow.config
    │       └── utils_nfschema_plugin
    │           ├── main.nf
    │           ├── meta.yml
    │           └── tests
    │               ├── main.nf.test
    │               ├── nextflow.config
    │               └── nextflow_schema.json
    └── workflows
        └── hello.nf

    14 directories, 34 files
    ```

यह बहुत सारी फ़ाइलें हैं!

उम्मीद है कि तुम उनमें से बहुत सारी को पहचान लोगे जो हमने `nf-core/demo` पाइपलाइन संरचना की खोज करते समय देखी थीं।
लेकिन चिंता मत करो अगर तुम अभी भी थोड़ा खोया हुआ महसूस कर रहे हो; हम इस प्रशिक्षण के दौरान महत्वपूर्ण भागों के माध्यम से एक साथ चलेंगे।

!!! note "नोट"

    इस प्रशिक्षण के पहले भाग में हमने जिस `nf-core/demo` पाइपलाइन की जांच की थी, उसकी तुलना में एक महत्वपूर्ण अंतर यह है कि कोई `modules` डायरेक्टरी नहीं है।
    ऐसा इसलिए है क्योंकि हमने डिफ़ॉल्ट nf-core मॉड्यूल में से किसी को भी शामिल करने का चुनाव नहीं किया।

### 1.2. जांचना कि स्कैफोल्ड कार्यात्मक है

विश्वास करो या न करो, भले ही तुमने अभी तक वास्तविक काम करने के लिए कोई मॉड्यूल नहीं जोड़ा है, पाइपलाइन स्कैफोल्ड वास्तव में test प्रोफाइल का उपयोग करके चलाया जा सकता है, उसी तरह जैसे हमने `nf-core/demo` पाइपलाइन चलाई थी।

```bash
nextflow run ./core-hello -profile docker,test --outdir core-hello-results
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `./core-hello/main.nf` [scruffy_marconi] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-47-18

    Core Nextflow options
      runName                   : scruffy_marconi
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    -[core/hello] Pipeline completed successfully-
    ```

यह तुम्हें दिखाता है कि सभी बुनियादी वायरिंग मौजूद है।
तो आउटपुट कहां हैं? क्या कोई हैं?

वास्तव में, `core-hello-results` नाम की परिणामों की एक नई डायरेक्टरी बनाई गई जिसमें मानक निष्पादन रिपोर्ट शामिल हैं:

```bash
tree core-hello-results
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-18.json
        └── pipeline_dag_2025-11-21_04-47-18.html

    1 directory, 6 files
    ```

तुम रिपोर्ट पर एक नज़र डाल सकते हो कि क्या चलाया गया था, और जवाब है: बिल्कुल कुछ नहीं!

![खाली निष्पादन टाइमलाइन रिपोर्ट](./img/execution_timeline_empty.png)

चलो देखते हैं कि वास्तव में कोड में क्या है।

### 1.3. प्लेसहोल्डर वर्कफ़्लो की जांच करना

यदि तुम `main.nf` फ़ाइल के अंदर देखो, तो तुम्हें दिखाई देगा कि यह `workflows/hello` से `HELLO` नाम का एक वर्कफ़्लो इम्पोर्ट करता है।

यह भाग 1 में हमारे सामने आए `workflows/demo.nf` वर्कफ़्लो के बराबर है, और हमारे रुचि के वर्कफ़्लो के लिए एक प्लेसहोल्डर वर्कफ़्लो के रूप में कार्य करता है, जिसमें कुछ nf-core कार्यक्षमता पहले से ही मौजूद है।

```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="15 17 19 35"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

[Hello Nextflow](../hello_nextflow/index.md) में विकसित किए गए जैसे बुनियादी Nextflow वर्कफ़्लो की तुलना में, तुम्हें कुछ चीजें दिखाई देंगी जो यहां नई हैं (ऊपर हाइलाइट की गई पंक्तियां):

- वर्कफ़्लो ब्लॉक का एक नाम है
- वर्कफ़्लो इनपुट `take:` कीवर्ड का उपयोग करके घोषित किए जाते हैं और चैनल निर्माण को पैरेंट वर्कफ़्लो में ले जाया जाता है
- वर्कफ़्लो सामग्री को `main:` ब्लॉक के अंदर रखा जाता है
- आउटपुट `emit:` कीवर्ड का उपयोग करके घोषित किए जाते हैं

ये Nextflow की वैकल्पिक विशेषताएं हैं जो वर्कफ़्लो को **composable** बनाती हैं, जिसका अर्थ है कि इसे किसी अन्य वर्कफ़्लो के भीतर से कॉल किया जा सकता है।

!!! note "गहराई से Composable वर्कफ़्लो"

    [Workflows of Workflows](../side_quests/workflows_of_workflows.md) साइड क्वेस्ट वर्कफ़्लो कंपोजिशन को बहुत अधिक गहराई से खोजता है, जिसमें कई वर्कफ़्लो को एक साथ कैसे कंपोज़ किया जाए और उनके बीच जटिल डेटा फ्लो को कैसे प्रबंधित किया जाए शामिल है। हम यहां composability का परिचय दे रहे हैं क्योंकि यह nf-core टेम्पलेट आर्किटेक्चर की एक मौलिक आवश्यकता है, जो पाइपलाइन इनिशियलाइज़ेशन, मुख्य विश्लेषण वर्कफ़्लो, और पूर्णता कार्यों को अलग, पुन: प्रयोज्य घटकों में व्यवस्थित करने के लिए नेस्टेड वर्कफ़्लो का उपयोग करता है।

हमें उस संरचना में अपने रुचि के वर्कफ़्लो से प्रासंगिक लॉजिक को प्लग करने की आवश्यकता होगी।
इसके लिए पहला कदम हमारे मूल वर्कफ़्लो को composable बनाना है।

### सारांश

तुम अब जानते हो कि nf-core टूल का उपयोग करके पाइपलाइन स्कैफोल्ड कैसे बनाया जाए।

### आगे क्या है?

सीखो कि एक सरल वर्कफ़्लो को nf-core संगत बनाने की प्रस्तावना के रूप में composable कैसे बनाया जाए।

---

## 2. मूल Hello Nextflow वर्कफ़्लो को composable बनाना

अब हमारे वर्कफ़्लो को nf-core स्कैफोल्ड में एकीकृत करने का काम शुरू करने का समय है।
याद दिलाने के लिए, हम उस वर्कफ़्लो के साथ काम कर रहे हैं जो हमारे [Hello Nextflow](../hello_nextflow/index.md) प्रशिक्षण पाठ्यक्रम में प्रदर्शित है।

!!! tip "सुझाव"

    यदि तुम उस पाइपलाइन से परिचित नहीं हो या तुम्हें याद दिलाने की जरूरत है, तो [The Hello pipeline](../info/hello_pipeline.md) देखो।

हम तुम्हें `original-hello` डायरेक्टरी में पूर्ण Hello Nextflow वर्कफ़्लो की एक साफ, पूरी तरह से कार्यात्मक प्रति प्रदान करते हैं, इसके मॉड्यूल और डिफ़ॉल्ट CSV फ़ाइल के साथ जिसे यह इनपुट के रूप में उपयोग करने की अपेक्षा करता है।

```bash
tree original-hello/
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    original-hello/
    ├── hello.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config

    1 directory, 6 files
    ```

इसे चलाने के लिए स्वतंत्र महसूस करो कि यह काम करता है:

```bash
nextflow run original-hello/hello.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/hello.nf` [goofy_babbage] DSL2 - revision: e9e72441e9

    executor >  local (8)
    [a4/081cec] sayHello (1)       | 3 of 3 ✔
    [e7/7e9058] convertToUpper (3) | 3 of 3 ✔
    [0c/17263b] collectGreetings   | 1 of 1 ✔
    [94/542280] cowpy              | 1 of 1 ✔
    ```

चलो `hello.nf` वर्कफ़्लो फ़ाइल खोलते हैं और कोड का निरीक्षण करते हैं, जो नीचे पूर्ण रूप से दिखाया गया है (प्रोसेस को छोड़कर, जो मॉड्यूल में हैं):

```groovy title="original-hello/hello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow {

  // create a channel for inputs from a CSV file
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // emit a greeting
  sayHello(greeting_ch)

  // convert the greeting to uppercase
  convertToUpper(sayHello.out)

  // collect all the greetings into one file
  collectGreetings(convertToUpper.out.collect(), params.batch)

  // generate ASCII art of the greetings with cowpy
  cowpy(collectGreetings.out.outfile, params.character)
}
```

जैसा कि तुम देख सकते हो, यह वर्कफ़्लो एक सरल अनाम वर्कफ़्लो के रूप में लिखा गया था जिसे अपने आप चलाया जा सकता है।
इसे पैरेंट वर्कफ़्लो के भीतर से चलाने योग्य बनाने के लिए जैसा कि nf-core टेम्पलेट की आवश्यकता है, हमें इसे **composable** बनाने की आवश्यकता है।

चलो आवश्यक परिवर्तनों के माध्यम से एक-एक करके चलते हैं।

### 2.1. वर्कफ़्लो को नाम देना

सबसे पहले, चलो वर्कफ़्लो को एक नाम देते हैं ताकि हम इसे पैरेंट वर्कफ़्लो से संदर्भित कर सकें।

=== "बाद में"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow HELLO {
    ```

=== "पहले"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow {
    ```

वर्कफ़्लो नामों पर वही सम्मेलन लागू होते हैं जो मॉड्यूल नामों पर लागू होते हैं।

### 2.2. चैनल निर्माण को `take` से बदलना

अब, चैनल निर्माण को एक सरल `take` स्टेटमेंट से बदलो जो अपेक्षित इनपुट घोषित करता है।

=== "बाद में"

    ```groovy title="original-hello/hello.nf" linenums="18"
        take:
        // channel of greetings
        greeting_ch
    ```

=== "पहले"

    ```groovy title="original-hello/hello.nf" linenums="18"
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.greeting)
                            .splitCsv()
                            .map { line -> line[0] }
    ```

यह इनपुट कैसे प्रदान किए जाते हैं इसके विवरण को पैरेंट वर्कफ़्लो पर छोड़ देता है।

जब हम इस पर हैं, हम `params.greeting = 'greetings.csv'` लाइन को भी कमेंट आउट कर सकते हैं

=== "बाद में"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        //params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

=== "पहले"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

!!! note "नोट"

    यदि तुम्हारे पास Nextflow language server एक्सटेंशन इंस्टॉल है, तो सिंटैक्स चेकर तुम्हारे कोड को लाल squiggles के साथ रोशन करेगा।
    ऐसा इसलिए है क्योंकि यदि तुम `take:` स्टेटमेंट डालते हो, तो तुम्हें `main:` भी होना चाहिए।

    हम इसे अगले चरण में जोड़ेंगे।

### 2.3. वर्कफ़्लो ऑपरेशन को `main` स्टेटमेंट से पहले रखना

अगला, वर्कफ़्लो के बॉडी में कॉल किए गए बाकी ऑपरेशन से पहले एक `main` स्टेटमेंट जोड़ो।

=== "बाद में"

    ```groovy title="original-hello/hello.nf" linenums="22" hl_lines="1"
        main:

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "पहले"

    ```groovy title="original-hello/hello.nf" linenums="21"
        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

यह मूल रूप से कहता है 'यह वह है जो यह वर्कफ़्लो _करता है_'।

### 2.4. `emit` स्टेटमेंट जोड़ना

अंत में, एक `emit` स्टेटमेंट जोड़ो जो घोषित करता है कि वर्कफ़्लो के अंतिम आउटपुट क्या हैं।

```groovy title="original-hello/hello.nf" linenums="35"
    emit:
    cowpy_hellos = cowpy.out
```

यह मूल वर्कफ़्लो की तुलना में कोड में एक नया जोड़ है।

### 2.5. पूर्ण परिवर्तनों का सारांश

यदि तुमने सभी परिवर्तन वर्णित के अनुसार किए हैं, तो तुम्हारा वर्कफ़्लो अब इस तरह दिखना चाहिए:

```groovy title="original-hello/hello.nf" linenums="1" hl_lines="16 18-20 22 36-37"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
// params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow HELLO {

    take:
    // channel of greetings
    greeting_ch

    main:

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // generate ASCII art of the greetings with cowpy
    cowpy(collectGreetings.out.outfile, params.character)

    emit:
    cowpy_hellos = cowpy.out
}
```

यह सब कुछ वर्णन करता है जो Nextflow को चाहिए सिवाय इनपुट चैनल में क्या फीड करना है।
वह पैरेंट वर्कफ़्लो में परिभाषित किया जाएगा, जिसे **entrypoint** वर्कफ़्लो भी कहा जाता है।

### 2.6. एक डमी entrypoint वर्कफ़्लो बनाना

हमारे composable वर्कफ़्लो को जटिल nf-core स्कैफोल्ड में एकीकृत करने से पहले, चलो सत्यापित करते हैं कि यह सही तरीके से काम करता है।
हम composable वर्कफ़्लो को अलगाव में परीक्षण करने के लिए एक सरल डमी entrypoint वर्कफ़्लो बना सकते हैं।

उसी `original-hello` डायरेक्टरी में `main.nf` नाम की एक खाली फ़ाइल बनाओ।

```bash
touch original-hello/main.nf
```

निम्नलिखित कोड को `main.nf` फ़ाइल में कॉपी करो।

```groovy title="original-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow

// import the workflow code from the hello.nf file
include { HELLO } from './hello.nf'

// declare input parameter
params.greeting = 'greetings.csv'

workflow {
  // create a channel for inputs from a CSV file
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // call the imported workflow on the channel of greetings
  HELLO(greeting_ch)

  // view the outputs emitted by the workflow
  HELLO.out.view { output -> "Output: $output" }
}
```

यहां दो महत्वपूर्ण अवलोकन करने हैं:

- इम्पोर्ट किए गए वर्कफ़्लो को कॉल करने के लिए सिंटैक्स मूल रूप से मॉड्यूल को कॉल करने के लिए सिंटैक्स के समान है।
- वर्कफ़्लो में इनपुट खींचने से संबंधित सब कुछ (इनपुट पैरामीटर और चैनल निर्माण) अब इस पैरेंट वर्कफ़्लो में घोषित किया गया है।

!!! note "नोट"

    entrypoint वर्कफ़्लो फ़ाइल को `main.nf` नाम देना एक सम्मेलन है, आवश्यकता नहीं।

    यदि तुम इस सम्मेलन का पालन करते हो, तो तुम अपने `nextflow run` कमांड में वर्कफ़्लो फ़ाइल नाम निर्दिष्ट करना छोड़ सकते हो।
    Nextflow स्वचालित रूप से निष्पादन डायरेक्टरी में `main.nf` नाम की फ़ाइल की तलाश करेगा।

    हालांकि, यदि तुम चाहो तो entrypoint वर्कफ़्लो फ़ाइल को कुछ और नाम दे सकते हो।
    उस स्थिति में, अपने `nextflow run` कमांड में वर्कफ़्लो फ़ाइल नाम निर्दिष्ट करना सुनिश्चित करो।

### 2.7. जांचना कि वर्कफ़्लो चलता है

हमारे पास अंततः सभी टुकड़े हैं जो हमें यह सत्यापित करने की आवश्यकता है कि composable वर्कफ़्लो काम करता है।
चलो इसे चलाते हैं!

```bash
nextflow run ./original-hello
```

यहां तुम `main.nf` नामकरण सम्मेलन का उपयोग करने का लाभ देखते हो।
यदि हमने entrypoint वर्कफ़्लो को `something_else.nf` नाम दिया होता, तो हमें `nextflow run original-hello/something_else.nf` करना पड़ता।

यदि तुमने सभी परिवर्तन सही तरीके से किए हैं, तो यह पूर्णता तक चलना चाहिए।

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/main.nf` [friendly_wright] DSL2 - revision: 1ecd2d9c0a

    executor >  local (8)
    [24/c6c0d8] HELLO:sayHello (3)       | 3 of 3 ✔
    [dc/721042] HELLO:convertToUpper (3) | 3 of 3 ✔
    [48/5ab2df] HELLO:collectGreetings   | 1 of 1 ✔
    [e3/693b7e] HELLO:cowpy              | 1 of 1 ✔
    Output: /workspaces/training/hello-nf-core/work/e3/693b7e48dc119d0c54543e0634c2e7/cowpy-COLLECTED-test-batch-output.txt
    ```

इसका मतलब है कि हमने सफलतापूर्वक अपने HELLO वर्कफ़्लो को composable बनाने के लिए अपग्रेड किया है।

### सारांश

तुम जानते हो कि एक वर्कफ़्लो को नाम देकर और `take`, `main` और `emit` स्टेटमेंट जोड़कर composable कैसे बनाया जाए, और इसे entrypoint वर्कफ़्लो से कैसे कॉल किया जाए।

### आगे क्या है?

सीखो कि एक बुनियादी composable वर्कफ़्लो को nf-core स्कैफोल्ड पर कैसे ग्राफ्ट किया जाए।

---

## 3. अपडेट किए गए वर्कफ़्लो लॉजिक को प्लेसहोल्डर वर्कफ़्लो में फिट करना

अब जब हमने सत्यापित कर लिया है कि हमारा composable वर्कफ़्लो सही तरीके से काम करता है, चलो खंड 1 में बनाए गए nf-core पाइपलाइन स्कैफोल्ड पर वापस आते हैं।
हम अभी विकसित किए गए composable वर्कफ़्लो को nf-core टेम्पलेट संरचना में एकीकृत करना चाहते हैं, इसलिए अंतिम परिणाम कुछ इस तरह दिखना चाहिए।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/core-hello.svg"
</figure>

तो हम यह कैसे करें? चलो `core-hello/workflows/hello.nf` (nf-core स्कैफोल्ड) में `HELLO` वर्कफ़्लो की वर्तमान सामग्री पर एक नज़र डालते हैं।

```groovy title="core-hello/workflows/hello.nf" linenums="1"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

कुल मिलाकर यह कोड सॉफ्टवेयर टूल के संस्करण को कैप्चर करने से संबंधित कुछ हाउसकीपिंग के अलावा बहुत कम करता है जो पाइपलाइन में चलते हैं।

हमें खंड 2 में विकसित किए गए मूल वर्कफ़्लो के composable संस्करण से प्रासंगिक कोड जोड़ने की आवश्यकता है।

हम इसे निम्नलिखित चरणों में निपटाने जा रहे हैं:

1. मॉड्यूल को कॉपी करना और मॉड्यूल इम्पोर्ट सेट अप करना
2. `take` घोषणा को जैसा है वैसा छोड़ना
3. `main` ब्लॉक में वर्कफ़्लो लॉजिक जोड़ना
4. `emit` ब्लॉक को अपडेट करना

!!! note "नोट"

    हम इस पहले पास के लिए संस्करण कैप्चर को अनदेखा करने जा रहे हैं और बाद में इस प्रशिक्षण के एक भाग में देखेंगे कि इसे कैसे वायर किया जाए।

### 3.1. मॉड्यूल को कॉपी करना और मॉड्यूल इम्पोर्ट सेट अप करना

हमारे Hello Nextflow वर्कफ़्लो से चार प्रोसेस `original-hello/modules/` में मॉड्यूल के रूप में संग्रहीत हैं।
हमें उन मॉड्यूल को nf-core प्रोजेक्ट संरचना (`core-hello/modules/local/` के तहत) में कॉपी करने और nf-core वर्कफ़्लो फ़ाइल में इम्पोर्ट स्टेटमेंट जोड़ने की आवश्यकता है।

पहले चलो मॉड्यूल फ़ाइलों को `original-hello/` से `core-hello/` में कॉपी करते हैं:

```bash
mkdir -p core-hello/modules/local/
cp original-hello/modules/* core-hello/modules/local/.
```

तुम्हें अब `core-hello/` के तहत सूचीबद्ध मॉड्यूल की डायरेक्टरी दिखाई देनी चाहिए।

```bash
tree core-hello/modules
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    core-hello/modules
    └── local
        ├── collectGreetings.nf
        ├── convertToUpper.nf
        ├── cowpy.nf
        └── sayHello.nf

    1 directory, 4 files
    ```

अब चलो मॉड्यूल इम्पोर्ट स्टेटमेंट सेट अप करते हैं।

ये `original-hello/hello.nf` वर्कफ़्लो में इम्पोर्ट स्टेटमेंट थे:

```groovy title="original-hello/hello.nf" linenums="9"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'
```

`core-hello/workflows/hello.nf` फ़ाइल खोलो और उन इम्पोर्ट स्टेटमेंट को इसमें नीचे दिखाए अनुसार ट्रांसपोज़ करो।

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="8-11"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    ```

यहां दो और दिलचस्प अवलोकन हैं:

- हमने nf-core स्टाइल सम्मेलन का पालन करने के लिए इम्पोर्ट स्टेटमेंट की फॉर्मेटिंग को अनुकूलित किया है।
- हमने मॉड्यूल के सापेक्ष पथ को अपडेट किया है ताकि यह प्रतिबिंबित हो कि वे अब नेस्टिंग के एक अलग स्तर पर संग्रहीत हैं।

### 3.2. `take` घोषणा को जैसा है वैसा छोड़ना

nf-core प्रोजेक्ट में samplesheet की अवधारणा के आसपास बहुत सारी पूर्व-निर्मित कार्यक्षमता है, जो आमतौर पर स्तंभीय डेटा युक्त एक CSV फ़ाइल होती है।
चूंकि यह मूल रूप से हमारी `greetings.csv` फ़ाइल है, हम वर्तमान `take` घोषणा को जैसा है वैसा रखेंगे, और अगले चरण में बस इनपुट चैनल का नाम अपडेट करेंगे।

```groovy title="core-hello/workflows/hello.nf" linenums="21"
    take:
    ch_samplesheet // channel: samplesheet read in from --input
```

इनपुट हैंडलिंग इस वर्कफ़्लो के अपस्ट्रीम (इस कोड फ़ाइल में नहीं) की जाएगी।

### 3.3. `main` ब्लॉक में वर्कफ़्लो लॉजिक जोड़ना

अब जब हमारे मॉड्यूल वर्कफ़्लो के लिए उपलब्ध हैं, हम `main` ब्लॉक में वर्कफ़्लो लॉजिक को प्लग कर सकते हैं।

याद दिलाने के लिए, यह मूल वर्कफ़्लो में प्रासंगिक कोड है, जो जब हमने इसे composable बनाया तो बहुत अधिक नहीं बदला (हमने बस `main:` लाइन जोड़ी):

```groovy title="original-hello/hello.nf" linenums="22"
    main:

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // generate ASCII art of the greetings with cowpy
    cowpy(collectGreetings.out.outfile, params.character)
```

हमें `main:` के बाद आने वाले कोड को वर्कफ़्लो के नए संस्करण में कॉपी करने की आवश्यकता है।

वहां पहले से ही कुछ कोड है जो वर्कफ़्लो द्वारा चलाए जाने वाले टूल के संस्करणों को कैप्चर करने से संबंधित है। हम अभी के लिए इसे अकेला छोड़ने जा रहे हैं (हम बाद में टूल संस्करणों से निपटेंगे)।
हम शीर्ष पर `ch_versions = channel.empty()` इनिशियलाइज़ेशन रखेंगे, फिर अपना वर्कफ़्लो लॉजिक डालेंगे, अंत में संस्करण संकलन कोड रखेंगे।
यह क्रम समझ में आता है क्योंकि एक वास्तविक पाइपलाइन में, प्रोसेस संस्करण जानकारी emit करेंगे जो वर्कफ़्लो चलने के दौरान `ch_versions` चैनल में जोड़ी जाएगी।

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19" hl_lines="10-20"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input

        main:

        ch_versions = Channel.empty()

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input
        main:

        ch_versions = Channel.empty()

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

तुम देखोगे कि हमने कोड को अधिक पठनीय बनाने के लिए `main:` से पहले एक खाली लाइन भी जोड़ी।

यह बहुत अच्छा लगता है, लेकिन हमें अभी भी चैनल का नाम अपडेट करने की आवश्यकता है जिसे हम `sayHello()` प्रोसेस को `greeting_ch` से `ch_samplesheet` में पास कर रहे हैं जैसा कि नीचे दिखाया गया है, `take:` कीवर्ड के तहत लिखे गए से मेल खाने के लिए।

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emit a greeting (updated to use the nf-core convention for samplesheets)
        sayHello(ch_samplesheet)
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emit a greeting
        sayHello(greeting_ch)
    ```

अब वर्कफ़्लो लॉजिक सही तरीके से वायर्ड है।

### 3.4. `emit` ब्लॉक को अपडेट करना

अंत में, हमें वर्कफ़्लो के अंतिम आउटपुट की घोषणा शामिल करने के लिए `emit` ब्लॉक को अपडेट करने की आवश्यकता है।

=== "बाद में"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "पहले"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55"
        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

यह HELLO वर्कफ़्लो में ही करने के लिए आवश्यक संशोधनों को पूरा करता है।
इस बिंदु पर, हमने समग्र कोड संरचना प्राप्त कर ली है जिसे हमने लागू करने के लिए निर्धारित किया था।

### सारांश

तुम जानते हो कि एक composable वर्कफ़्लो के मुख्य टुकड़ों को nf-core प्लेसहोल्डर वर्कफ़्लो में कैसे फिट किया जाए।

### आगे क्या है?

सीखो कि nf-core पाइपलाइन स्कैफोल्ड में इनपुट को कैसे संभाला जाए।

---

## 4. इनपुट हैंडलिंग को अनुकूलित करना

अब जब हमने सफलतापूर्वक अपने वर्कफ़्लो लॉजिक को nf-core स्कैफोल्ड में एकीकृत कर लिया है, हमें एक और महत्वपूर्ण टुकड़े को संबोधित करने की आवश्यकता है: यह सुनिश्चित करना कि हमारा इनपुट डेटा सही तरीके से प्रोसेस किया गया है।
nf-core टेम्पलेट जटिल जीनोमिक्स डेटासेट के लिए डिज़ाइन की गई परिष्कृत इनपुट हैंडलिंग के साथ आता है, इसलिए हमें इसे हमारी सरल `greetings.csv` फ़ाइल के साथ काम करने के लिए अनुकूलित करने की आवश्यकता है।

### 4.1. पहचानना कि इनपुट कहां संभाले जाते हैं

पहला कदम यह पता लगाना है कि इनपुट हैंडलिंग कहां की जाती है।

तुम्हें याद होगा कि जब हमने Hello Nextflow वर्कफ़्लो को composable होने के लिए फिर से लिखा, तो हमने इनपुट पैरामीटर घोषणा को एक स्तर ऊपर, `main.nf` entrypoint वर्कफ़्लो में ले जाया।
तो चलो शीर्ष स्तर `main.nf` entrypoint वर्कफ़्लो पर एक नज़र डालते हैं जो पाइपलाइन स्कैफोल्ड के हिस्से के रूप में बनाया गया था:

```groovy title="core-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core/hello
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/core/hello
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HELLO  } from './workflows/hello'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hello_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hello_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CORE_HELLO {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    HELLO (
        samplesheet
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    CORE_HELLO (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

nf-core प्रोजेक्ट नेस्टेड subworkflows का भारी उपयोग करता है, इसलिए यह बिट पहले दृष्टिकोण पर थोड़ा भ्रमित करने वाला हो सकता है।

यहां जो मायने रखता है वह यह है कि दो वर्कफ़्लो परिभाषित हैं:

- `CORE_HELLO` HELLO वर्कफ़्लो को चलाने के लिए एक पतला रैपर है जिसे हमने अभी `core-hello/workflows/hello.nf` में अनुकूलित करना समाप्त किया।
- एक अनाम वर्कफ़्लो जो `CORE_HELLO` के साथ-साथ दो अन्य subworkflows, `PIPELINE_INITIALISATION` और `PIPELINE_COMPLETION` को कॉल करता है।

यहां एक आरेख है कि वे एक दूसरे से कैसे संबंधित हैं:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/hello-nested-workflows.svg"
</figure>

महत्वपूर्ण रूप से, हम इस स्तर पर इनपुट चैनल का निर्माण करने वाला कोई कोड नहीं पा सकते हैं, केवल `--input` पैरामीटर के माध्यम से प्रदान की गई samplesheet के संदर्भ।

थोड़ा खोजने से पता चलता है कि इनपुट हैंडलिंग `PIPELINE_INITIALISATION` subworkflow द्वारा की जाती है, उचित रूप से पर्याप्त, जिसे `core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf` से इम्पोर्ट किया गया है।

यदि हम उस फ़ाइल को खोलते हैं और नीचे स्क्रॉल करते हैं, तो हम इस कोड के टुकड़े पर आते हैं:

```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76"
    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
```

यह चैनल फैक्ट्री है जो samplesheet को पार्स करती है और इसे एक रूप में पास करती है जो HELLO वर्कफ़्लो द्वारा उपभोग के लिए तैयार है।

!!! note "नोट"

    ऊपर का सिंटैक्स हमने पहले उपयोग किए गए से थोड़ा अलग है, लेकिन मूल रूप से यह:

    ```groovy
    channel.<...>.set { ch_samplesheet }
    ```

    इसके बराबर है:

    ```groovy
    ch_samplesheet = channel.<...>
    ```

इस कोड में कुछ पार्सिंग और सत्यापन चरण शामिल हैं जो nf-core पाइपलाइन टेम्पलेट के साथ शामिल उदाहरण samplesheet के लिए अत्यधिक विशिष्ट हैं, जो लेखन के समय बहुत डोमेन-विशिष्ट है और हमारे सरल पाइपलाइन प्रोजेक्ट के लिए उपयुक्त नहीं है।

### 4.2. टेम्पलेट किए गए इनपुट चैनल कोड को बदलना

अच्छी खबर यह है कि हमारी पाइपलाइन की जरूरतें बहुत सरल हैं, इसलिए हम उस सभी को उस चैनल निर्माण कोड से बदल सकते हैं जिसे हमने मूल Hello Nextflow वर्कफ़्लो में विकसित किया था।

याद दिलाने के लिए, यह चैनल निर्माण कैसा दिखता था (जैसा कि समाधान डायरेक्टरी में देखा गया):

```groovy title="solutions/composable-hello/main.nf" linenums="10" hl_lines="4"
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)
        .splitCsv()
        .map { line -> line[0] }
```

तो हमें बस इसे इनिशियलाइज़ेशन वर्कफ़्लो में प्लग करने की आवश्यकता है, मामूली परिवर्तनों के साथ: हम चैनल नाम को `greeting_ch` से `ch_samplesheet` में अपडेट करते हैं, और पैरामीटर नाम को `params.greeting` से `params.input` में (हाइलाइट की गई लाइन देखें)।

=== "बाद में"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-7"
        //
        // Create channel from input file provided through params.input
        //

        ch_samplesheet = channel.fromPath(params.input)
            .splitCsv()
            .map { line -> line[0] }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

=== "पहले"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-23"
        //
        // Create channel from input file provided through params.input
        //

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

यह इनपुट प्रोसेसिंग को काम करने के लिए आवश्यक परिवर्तनों को पूरा करता है।

इसके वर्तमान रूप में, यह हमें स्कीमा सत्यापन के लिए nf-core की अंतर्निहित क्षमताओं का लाभ उठाने नहीं देगा, लेकिन हम इसे बाद में जोड़ सकते हैं।
अभी के लिए, हम इसे यथासंभव सरल रखने पर ध्यान केंद्रित कर रहे हैं ताकि हम परीक्षण डेटा पर सफलतापूर्वक चला सकें।

### 4.3. टेस्ट प्रोफाइल को अपडेट करना

परीक्षण डेटा और पैरामीटर की बात करते हुए, चलो इस पाइपलाइन के लिए टेस्ट प्रोफाइल को टेम्पलेट में प्रदान की गई उदाहरण samplesheet के बजाय `greetings.csv` मिनी-samplesheet का उपयोग करने के लिए अपडेट करते हैं।

`core-hello/conf` के तहत, हमें दो टेम्पलेट किए गए टेस्ट प्रोफाइल मिलते हैं: `test.config` और `test_full.config`, जो एक छोटे डेटा नमूने और एक पूर्ण आकार के परीक्षण के लिए हैं।
हमारी पाइपलाइन के उद्देश्य को देखते हुए, पूर्ण आकार के टेस्ट प्रोफाइल को सेट अप करने का वास्तव में कोई मतलब नहीं है, इसलिए `test_full.config` को अनदेखा करने या हटाने के लिए स्वतंत्र महसूस करो।
हम कुछ डिफ़ॉल्ट पैरामीटर के साथ हमारी `greetings.csv` फ़ाइल पर चलाने के लिए `test.config` सेट अप करने पर ध्यान केंद्रित करने जा रहे हैं।

#### 4.3.1. `greetings.csv` फ़ाइल को कॉपी करना

सबसे पहले हमें `greetings.csv` फ़ाइल को हमारे पाइपलाइन प्रोजेक्ट में एक उपयुक्त स्थान पर कॉपी करने की आवश्यकता है।
आमतौर पर छोटी परीक्षण फ़ाइलें `assets` डायरेक्टरी में संग्रहीत की जाती हैं, इसलिए चलो फ़ाइल को हमारी कार्य डायरेक्टरी से कॉपी करते हैं।

```bash
cp greetings.csv core-hello/assets/.
```

अब `greetings.csv` फ़ाइल परीक्षण इनपुट के रूप में उपयोग के लिए तैयार है।

#### 4.3.2. `test.config` फ़ाइल को अपडेट करना

अब हम `test.config` फ़ाइल को निम्नानुसार अपडेट कर सकते हैं:

=== "बाद में"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-10"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Input data
        input  = "${projectDir}/assets/greetings.csv"

        // Other parameters
        batch     = 'test'
        character = 'tux'
    }
    ```

=== "पहले"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-8"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Input data
        // TODO nf-core: Specify the paths to your test data on nf-core/test-datasets
        // TODO nf-core: Give any required params for the test so that command line flags are not needed
        input  = params.pipelines_testdata_base_path + 'viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
    }
    ```

मुख्य बिंदु:

- **`${projectDir}` का उपयोग करना**: यह एक Nextflow implicit variable है जो उस डायरेक्टरी की ओर इशारा करता है जहां मुख्य वर्कफ़्लो स्क्रिप्ट स्थित है (पाइपलाइन रूट)। इसका उपयोग सुनिश्चित करता है कि पथ काम करता है चाहे पाइपलाइन कहीं से भी चलाई जाए।
- **Absolute paths**: `${projectDir}` का उपयोग करके, हम एक absolute path बनाते हैं, जो पाइपलाइन के साथ शिप होने वाले परीक्षण डेटा के लिए महत्वपूर्ण है।
- **परीक्षण डेटा स्थान**: nf-core पाइपलाइन आमतौर पर छोटी परीक्षण फ़ाइलों के लिए पाइपलाइन रिपॉजिटरी के भीतर `assets/` डायरेक्टरी में परीक्षण डेटा संग्रहीत करती हैं, या बड़ी फ़ाइलों के लिए बाहरी परीक्षण डेटासेट का संदर्भ देती हैं।

और जब हम इस पर हैं, चलो डिफ़ॉल्ट संसाधन सीमाओं को कसते हैं ताकि यह सुनिश्चित हो सके कि यह बहुत बुनियादी मशीनों (जैसे Github Codespaces में न्यूनतम VM) पर चलेगा:

=== "बाद में"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 2,
            memory: '4.GB',
            time: '1.h'
        ]
    }
    ```

=== "पहले"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 4,
            memory: '15.GB',
            time: '1.h'
        ]
    }
    ```

यह हमें करने के लिए आवश्यक कोड संशोधनों को पूरा करता है।

### 4.4. टेस्ट प्रोफाइल के साथ पाइपलाइन चलाना

यह बहुत कुछ था, लेकिन हम अंततः पाइपलाइन चलाने की कोशिश कर सकते हैं!
ध्यान दें कि हमें कमांड लाइन में `--validate_params false` जोड़ना होगा क्योंकि हमने अभी तक सत्यापन सेट अप नहीं किया है (वह बाद में आएगा)।

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

यदि तुमने सभी संशोधन सही तरीके से किए हैं, तो यह पूर्णता तक चलना चाहिए।

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `core-hello/main.nf` [condescending_allen] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-11-21_07-29-37

    Core Nextflow options
      runName                   : condescending_allen
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (1)
    [ed/727b7e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [45/bb6096] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [81/7e2e34] CORE_HELLO:HELLO:collectGreetings   [100%] 1 of 1 ✔
    [96/9442a1] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

जैसा कि तुम देख सकते हो, इसने इनिशियलाइज़ेशन subworkflow के लिए धन्यवाद शुरुआत में विशिष्ट nf-core सारांश तैयार किया, और प्रत्येक मॉड्यूल के लिए लाइनें अब पूर्ण PIPELINE:WORKFLOW:module नाम दिखाती हैं।

### 4.5. पाइपलाइन आउटपुट खोजना

अब सवाल यह है: पाइपलाइन के आउटपुट कहां हैं?
और जवाब काफी दिलचस्प है: अब परिणामों को देखने के लिए दो अलग-अलग स्थान हैं।

जैसा कि तुम्हें पहले से याद होगा, हमारे नए बनाए गए वर्कफ़्लो के पहले रन ने `core-hello-results/` नाम की एक डायरेक्टरी तैयार की जिसमें विभिन्न निष्पादन रिपोर्ट और मेटाडेटा शामिल थे।

```bash
tree core-hello-results
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_report_2025-11-21_07-29-37.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_07-29-37.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── execution_trace_2025-11-21_07-29-37.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-13.json
        ├── params_2025-11-21_07-29-41.json
        └── pipeline_dag_2025-11-21_04-47-18.html
        └── pipeline_dag_2025-11-21_07-29-37.html

    1 directory, 12 files
    ```

तुम देखते हो कि हमें पहले रन से मिली रिपोर्ट के अलावा निष्पादन रिपोर्ट का एक और सेट मिला, जब वर्कफ़्लो अभी भी सिर्फ एक प्लेसहोल्डर था।
इस बार तुम सभी कार्यों को देखते हो जो अपेक्षित रूप से चलाए गए थे।

![Hello पाइपलाइन के लिए निष्पादन टाइमलाइन रिपोर्ट](./img/execution_timeline_hello.png)

!!! note "नोट"

    एक बार फिर कार्यों को समानांतर में नहीं चलाया गया क्योंकि हम Github Codespaces में एक न्यूनतम मशीन पर चल रहे हैं।
    इन्हें समानांतर में चलते हुए देखने के लिए, अपने codespace के CPU आवंटन और टेस्ट कॉन्फ़िगरेशन में संसाधन सीमाओं को बढ़ाने का प्रयास करो।

यह बहुत अच्छा है, लेकिन हमारे वास्तविक पाइपलाइन परिणाम वहां नहीं हैं!

यहां क्या हुआ: हमने मॉड्यूल में ही कुछ भी नहीं बदला, इसलिए मॉड्यूल-स्तर `publishDir` निर्देशों द्वारा संभाले गए आउटपुट अभी भी मूल पाइपलाइन में निर्दिष्ट `results` डायरेक्टरी में जा रहे हैं।

```bash
tree results
```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results
    ├── Bonjour-output.txt
    ├── COLLECTED-test-batch-output.txt
    ├── COLLECTED-test-output.txt
    ├── cowpy-COLLECTED-test-batch-output.txt
    ├── cowpy-COLLECTED-test-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt

    0 directories, 10 files
    ```

आह, वे वहां हैं, मूल Hello पाइपलाइन के पहले के रन के आउटपुट के साथ मिश्रित।

यदि हम चाहते हैं कि वे डेमो पाइपलाइन के आउटपुट की तरह बड़े करीने से व्यवस्थित हों, तो हमें यह बदलने की आवश्यकता होगी कि हम आउटपुट को प्रकाशित करने के लिए कैसे सेट अप करते हैं।
हम तुम्हें इस प्रशिक्षण पाठ्यक्रम में बाद में दिखाएंगे कि यह कैसे करना है।

<!-- TODO: Update this once we've updated Hello Nextflow to use workflow-level outputs -->

और यह है! यह मूल पाइपलाइन के समान परिणाम प्राप्त करने के लिए बहुत काम की तरह लग सकता है, लेकिन तुम्हें वे सभी प्यारी रिपोर्ट स्वचालित रूप से उत्पन्न होती हैं, और अब तुम्हारे पास इनपुट सत्यापन और कुछ साफ मेटाडेटा हैंडलिंग क्षमताओं सहित nf-core की अतिरिक्त सुविधाओं का लाभ उठाने के लिए एक ठोस नींव है जिसे हम बाद के खंड में कवर करेंगे।

---

### सारांश

तुम जानते हो कि nf-core टेम्पलेट का उपयोग करके एक नियमित Nextflow पाइपलाइन को nf-core स्टाइल पाइपलाइन में कैसे परिवर्तित किया जाए।
इसके हिस्से के रूप में, तुमने सीखा कि एक वर्कफ़्लो को composable कैसे बनाया जाए, और nf-core टेम्पलेट के तत्वों की पहचान कैसे की जाए जिन्हें आमतौर पर कस्टम nf-core स्टाइल पाइपलाइन विकसित करते समय अनुकूलित करने की आवश्यकता होती है।

### आगे क्या है?

एक ब्रेक लो, यह कठिन काम था! जब तुम तैयार हो, तो [भाग 3: एक nf-core मॉड्यूल का उपयोग करें](./03_use_module.md) पर जाओ ताकि nf-core/modules रिपॉजिटरी से समुदाय-रखरखाव वाले मॉड्यूल का लाभ उठाना सीख सको।
