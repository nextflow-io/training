# भाग 1: Nextflow चलाना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

इस भाग में, हम Nextflow पाइपलाइन चलाने की मुख्य अवधारणाओं से परिचय कराते हैं।
हम एक सरल Hello World वर्कफ़्लो से शुरू करते हैं, फिर एक पूर्ण बहु-चरण पाइपलाइन की ओर बढ़ते हैं जो कंटेनर का उपयोग करके समानांतर में कई इनपुट प्रोसेस करती है।

---

## 1. Hello World

वर्कफ़्लो `1-hello.nf` एक कमांड-लाइन आर्गुमेंट के ज़रिए एक अभिवादन लेता है और उसे एक फ़ाइल में लिखता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

### 1.1. वर्कफ़्लो लॉन्च करना

अपने टर्मिनल में निम्नलिखित कमांड चलाओ।

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `1-hello.nf` [infallible_volhard] revision: 82d40dbf94

    executor >  local (1)
    [6d/740edd] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output: 1-hello/Hello World!-output.txt
    ```

आउटपुट में मुख्य लाइन प्रोसेस स्टेटस लाइन है:

```console
[6d/740edd] sayHello | 1 of 1 ✔
```

यह हमें बताती है कि `sayHello` प्रोसेस एक बार सफलतापूर्वक चली।
`[6d/740edd]` प्रीफ़िक्स कार्य की work directory का एक छोटा रूप है — इसके बारे में नीचे और जानकारी है।
इसके बाद आने वाला `Outputs:` ब्लॉक हर उस फ़ाइल को सूचीबद्ध करता है जिसे पाइपलाइन ने प्रकाशित किया, जो नीचे [1.4](#14-optional-code-walkthrough) में बताए गए `output` ब्लॉक के अनुसार लेबल की गई है।

### 1.2. आउटपुट खोजना

यह वर्कफ़्लो अपना आउटपुट एक `results` डायरेक्टरी में प्रकाशित करने के लिए कॉन्फ़िगर किया गया है।
चलाने के बाद, तुम्हें वहाँ आउटपुट मिलना चाहिए:

```console title="results/"
results
└── 1-hello
    └── Hello World!-output.txt
```

फ़ाइल खोलो और पुष्टि करो कि उसमें `Hello World!` है।

### 1.3. `work/` डायरेक्टरी एक्सप्लोर करना

पर्दे के पीछे, Nextflow `work/` नाम की डायरेक्टरी के अंदर हर प्रोसेस कॉल के लिए एक अनोखी कार्य डायरेक्टरी बनाता है।
कंसोल आउटपुट में दिखाया गया हैश (`[6d/740edd]`) उस डायरेक्टरी का पाथ है।

```bash
ls work/6d/740edd*
```

अंदर तुम्हें आउटपुट फ़ाइल के साथ कई छिपी हुई लॉग फ़ाइलें मिलेंगी:

- **`.command.sh`**: वह सटीक कमांड जो Nextflow ने चलाई
- **`.command.out`** / **`.command.err`**: प्रोसेस से stdout और stderr
- **`.command.log`**: संयुक्त लॉग आउटपुट
- **`.exitcode`**: प्रोसेस का exit code

`.command.sh` फ़ाइल डीबगिंग के समय विशेष रूप से उपयोगी है — यह ठीक-ठीक दिखाती है कि क्या execute हुआ।

### 1.4. वैकल्पिक: कोड वॉकथ्रू

अगर तुम सिर्फ पाइपलाइन चलाना चाहते हो तो कोड समझना ज़रूरी नहीं है, लेकिन अगर उत्सुक हो, तो एक नज़र डालना उचित है।

??? optional "इस अभ्यास से जुड़े कोड को एक्सप्लोर करने के लिए क्लिक करो"

    आइए `1-hello.nf` खोलें और इसके मुख्य घटकों को देखें।

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { sayHello } from './modules/sayHello.nf'

    /*
     * पाइपलाइन पैरामीटर
     */
    params {
        input: String
    }

    workflow {

        main:
        // एक अभिवादन emit करें
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    हम निम्नलिखित देखते हैं:

    - एक `include` स्टेटमेंट जो एक `process` मॉड्यूल की ओर इशारा करता है
    - एक `params` ब्लॉक जो पाइपलाइन पैरामीटर परिभाषित करता है
    - एक `workflow` ब्लॉक जो किए जाने वाले काम का वर्णन करता है
    - एक `output` ब्लॉक जो बताता है कि आउटपुट के साथ क्या करना है

    आइए प्रत्येक को बारी-बारी से देखें।

    ### `process` मॉड्यूल

    `include` स्टेटमेंट Nextflow को एक अलग कोड फ़ाइल से `sayHello` नाम की चीज़ लोड करने के लिए कहता है।

    ```groovy title="1-hello.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    ```

    उस फ़ाइल में, हमें `sayHello` नाम के एक प्रोसेस की परिभाषा मिलती है:

    ```groovy title="modules/sayHello.nf" linenums="4"
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

    एक **process** पाइपलाइन में एक एकल चरण परिभाषित करता है।
    यह अपने इनपुट, आउटपुट और execute होने वाली script घोषित करता है।
    `val` qualifier का मतलब है कि इनपुट एक सामान्य मान (string, number, आदि) है।
    `path` qualifier का मतलब है कि आउटपुट एक फ़ाइल पाथ है।

    प्रोसेस की परिभाषा मुख्य वर्कफ़्लो फ़ाइल में लिखना संभव है, लेकिन उन्हें अलग मॉड्यूल फ़ाइलों में रखने से वे पुन: उपयोग योग्य बनती हैं: एक ही मॉड्यूल को कई वर्कफ़्लो स्क्रिप्ट import कर सकती हैं।

    ### `params` ब्लॉक

    `params` ब्लॉक उन कमांड-लाइन पैरामीटर की घोषणा करता है जिन्हें वर्कफ़्लो स्वीकार करता है:

    ```groovy title="1-hello.nf" linenums="8"
    params {
        input: String
    }
    ```

    यहाँ घोषित कोई भी पैरामीटर कमांड लाइन पर double-dash (`--input`) के साथ उपलब्ध हो जाता है।
    समर्थित प्रकारों में `String`, `Integer`, `Float`, `Boolean`, और `Path` शामिल हैं।

    !!! tip "सुझाव"

        वर्कफ़्लो पैरामीटर हमेशा दो dashes (`--input`) का उपयोग करते हैं ताकि उन्हें Nextflow के अपने CLI flags से अलग किया जा सके, जो एक dash का उपयोग करते हैं (जैसे `-resume`)।

    ### `workflow` ब्लॉक

    **workflow** ब्लॉक dataflow logic परिभाषित करता है: कौन से प्रोसेस चलाने हैं और किस क्रम में।

    ```groovy title="1-hello.nf" linenums="12"
    workflow {

        main:
        // एक अभिवादन emit करें
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

    यहाँ केवल एक प्रोसेस को कॉल किया जा रहा है इसलिए यह बहुत सरल है; हम बाद में अधिक वास्तविक उदाहरण देखेंगे।

    `main:` सेक्शन `--input` मान के साथ `sayHello` प्रोसेस को कॉल करता है।
    `publish:` सेक्शन सूचीबद्ध करता है कि कौन से आउटपुट results डायरेक्टरी में कॉपी किए जाने चाहिए।

    ### `output` ब्लॉक

    फ़ाइल के नीचे `output` ब्लॉक गंतव्य पाथ और copy mode निर्दिष्ट करता है।

    ```groovy title="1-hello.nf" linenums="22"
    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    प्रत्येक नामित एंट्री वर्कफ़्लो में एक `publish:` लेबल से मेल खाती है और उसे `results/` के अंतर्गत एक उपडायरेक्टरी से मैप करती है।

### सारांश

तुम जानते हो कि Nextflow पाइपलाइन कैसे चलाएं और उसके आउटपुट कैसे खोजें, और तुम जानते हो कि काम `work/` के अंतर्गत कार्य डायरेक्टरी में execute होता है।

### आगे क्या है?

जानो कि Nextflow कई इनपुट को कुशलतापूर्वक कैसे संभालता है।

---

## 2. कई इनपुट प्रोसेस करना

वास्तविक दुनिया की पाइपलाइन आमतौर पर कई डेटा प्रोसेस करती हैं, न कि सिर्फ एक।
वर्कफ़्लो `2-inputs.nf` एक CSV फ़ाइल से पढ़ता है और प्रत्येक पंक्ति के लिए समानांतर में `sayHello` एक बार चलाता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

आइए पहले वर्कफ़्लो चलाएं, फिर देखें कि Nextflow इन कई इनपुट को संभालने के लिए किस तंत्र का उपयोग करता है।

### 2.1. वर्कफ़्लो चलाना

अपने टर्मिनल में निम्नलिखित कमांड चलाओ।

```bash
nextflow run 2-inputs.nf --input data/greetings.csv
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [nauseous_babbage] revision: b90778224d

    executor >  local (3)
    [66/de7844] sayHello (3) | 3 of 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
    ```

`3 of 3` हमें बताता है कि `sayHello` प्रोसेस तीन बार कॉल हुई, CSV में प्रत्येक पंक्ति के लिए एक बार।

`results` डायरेक्टरी में, तुम्हें अब तीन आउटपुट फ़ाइलें दिखनी चाहिए, प्रत्येक अभिवादन के लिए एक:

```console title="results/2-inputs/"
2-inputs
├── Bonjour-output.txt
├── Hello-output.txt
└── Hola-output.txt
```

किसी भी आउटपुट फ़ाइल को खोलो और पुष्टि करो कि प्रत्येक में एक अभिवादन है।

ऊपर दिखाया गया संक्षिप्त आउटपुट `sayHello` के लिए एक सारांश लाइन दिखाता है, लेकिन Nextflow ने वास्तव में तीन अलग-अलग कार्य executions लॉन्च किए, CSV में प्रत्येक पंक्ति के लिए एक, और जैसे ही तुम्हारी मशीन में संसाधन उपलब्ध हुए उन्हें समानांतर में चलाया।

जैसे [1.3](#13-explore-the-work-directory) में तुमने जो एकल कार्य एक्सप्लोर किया, इन तीनों executions में से प्रत्येक को `work/` के अंतर्गत अपनी खुद की कार्य डायरेक्टरी मिलती है, जो दूसरों से पूरी तरह अलग होती है:

```console title="work/"
work
├── 2d/276c63.../
│   ├── .command.sh
│   └── Hola-output.txt
├── ab/007682.../
│   ├── .command.sh
│   └── Bonjour-output.txt
└── d9/2476082.../
    ├── .command.sh
    └── Hello-output.txt
```

प्रत्येक `.command.sh` में केवल उस एक अभिवादन के लिए कमांड होती है:

```console title="work/2d/276c63.../.command.sh"
#!/bin/bash -ue
echo 'Hola' > 'Hola-output.txt'
```

यह अलगाव ही समानांतर execution को सुरक्षित बनाता है: एक ही समय पर चल रहे तीन कार्य कभी भी एक working directory साझा नहीं करते, इसलिए एक कार्य जो लिखता है वह दूसरे कार्य के लिखे से टकरा या उसे overwrite नहीं कर सकता, भले ही वे एक ही नाम की फ़ाइलें बनाएं।
यही कारण है कि `-resume` (आगे बताया गया) अलग-अलग कार्यों को स्वतंत्र रूप से cache और reuse कर सकता है: प्रत्येक कार्य के इनपुट, आउटपुट और लॉग पूरी तरह उसकी अपनी डायरेक्टरी के अंदर रहते हैं, कार्यों के बीच कुछ भी साझा नहीं होता जो sync से बाहर हो सके।

### 2.2. `-ansi-log false` के साथ वर्कफ़्लो फिर से चलाना

डिफ़ॉल्ट रूप से, Nextflow प्रति प्रोसेस एक सारांश लाइन में आउटपुट संक्षिप्त करता है।
प्रत्येक प्रोसेस कॉल को अलग-अलग सूचीबद्ध देखने के लिए, `-ansi-log false` जोड़ो:

```bash
nextflow run 2-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 26.04.4
    Launching `2-inputs.nf` [extravagant_bardeen] - revision: b90778224d
    [43/0bac1c] Submitted process > sayHello (1)
    [2d/99f604] Submitted process > sayHello (2)
    [6d/7578d7] Submitted process > sayHello (3)

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Hello-output.txt
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hola-output.txt
    ```

यह तीनों प्रोसेस कॉल और प्रत्येक के लिए बनाई गई अनोखी work उपडायरेक्टरी दिखाता है।

### 2.3. पूर्ण किए गए काम को छोड़ने के लिए `-resume` का उपयोग करना

अब विस्तारित इनपुट फ़ाइल पर स्विच करो, जो दो और अभिवादन जोड़ती है, और कमांड लाइन में `-resume` जोड़ो:

```bash
nextflow run 2-inputs.nf --input data/greetings-extended.csv -resume
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [adoring_mayer] revision: b90778224d

    executor >  local (2)
    [84/2f3067] sayHello (5) | 5 of 5, cached: 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Ciao-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
        - 2-inputs/Ola-output.txt
    ```

Nextflow ने केवल दो नए इनपुट चलाए।
पिछले रन में प्रोसेस किए गए तीन अभिवादन cache हो गए और स्वचालित रूप से reuse हुए।

यह बहु-चरण पाइपलाइन में उन चरणों के लिए execution छोड़ने के लिए भी काम करता है जो पहले से सफलतापूर्वक चल चुके हैं।
उदाहरण के लिए, अगर किसी सिस्टम एरर से पाइपलाइन रन बाधित हो गई, या अगर तुमने विकास में पाइपलाइन में नए चरण जोड़े।

`-resume` क्षमता लंबी पाइपलाइन में विशेष रूप से मूल्यवान है जहाँ विफलता से उबरना महत्वपूर्ण समय और संसाधन बचा सकता है।

### 2.4. वैकल्पिक: कोड वॉकथ्रू

अगर तुम सिर्फ पाइपलाइन चलाना चाहते हो तो कोड समझना ज़रूरी नहीं है, लेकिन अगर उत्सुक हो, तो एक नज़र डालना उचित है।

??? optional "इस अभ्यास से जुड़े कोड को एक्सप्लोर करने के लिए क्लिक करो"

    `2-inputs.nf` में मुख्य बदलाव वर्कफ़्लो के `main:` सेक्शन में है:

    ```groovy title="2-inputs.nf" linenums="14" hl_lines="3 4 5"
        main:
        // CSV फ़ाइल से इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)
    ```

    यहाँ जो तुम देख रहे हो उसे **channel** कहते हैं: एक queue construct जो इनपुट डेटा को इस तरह संभालता है जिससे operations को parallelize करना आसान हो जाता है।

    - `channel.fromPath(params.input)` `--input` के साथ दिए गए फ़ाइल पाथ से एक चैनल बनाता है
    - `.splitCsv()` CSV को पंक्तियों में parse करता है
    - `#!groovy .map { line -> line[0] }` प्रत्येक पंक्ति से पहला कॉलम निकालता है

    परिणाम एक चैनल है जिसमें `Hello`, `Bonjour`, और `Hola` हैं।
    जब `sayHello(greeting_ch)` को पास किया जाता है, तो Nextflow स्वचालित रूप से प्रत्येक आइटम के लिए प्रोसेस को एक बार कॉल करता है, और जब संसाधन उपलब्ध हों तो उन्हें समानांतर में चलाता है।

### सारांश

तुम जानते हो कि CSV फ़ाइल से कई इनपुट को समानांतर में कैसे प्रोसेस करें, और पूर्ण किए गए काम को दोहराने से बचने के लिए `-resume` का उपयोग कैसे करें।

### आगे क्या है?

जानो कि एक पूर्ण बहु-चरण पाइपलाइन चैनल का उपयोग करके प्रोसेस को एक साथ कैसे जोड़ती है, और विश्लेषण टूल और उनकी dependencies को प्रबंधित करने के लिए कंटेनर का उपयोग कैसे करें।

---

## 3. बहु-चरण पाइपलाइन चलाना

अब तक तुमने एक एकल प्रोसेस चलाई, फिर उसे इनपुट के एक सेट पर समानांतर में कई बार चलाया।
वास्तविक पाइपलाइन आमतौर पर आगे जाती हैं: वे कई प्रोसेस को एक साथ जोड़ती हैं, एक का आउटपुट अगले में फीड करती हैं, और अक्सर रास्ते में एक से अधिक सॉफ़्टवेयर पर निर्भर करती हैं।
वर्कफ़्लो `main.nf` इन दोनों को एक पूर्ण पाइपलाइन में एक साथ रखता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

प्रत्येक इनपुट अभिवादन सभी चार चरणों से गुज़रता है: `sayHello` इसे एक फ़ाइल में लिखता है, `convertToUpper` टेक्स्ट को uppercase में बदलता है, `collectGreetings` सभी परिणामों को एक फ़ाइल में मर्ज करता है, और `cowpy` एक containerized टूल का उपयोग करके मर्ज किए गए आउटपुट से ASCII art बनाता है।
Nextflow इन चरणों को चैनल से जोड़ता है: एक प्रोसेस का आउटपुट अगले का इनपुट बन जाता है, इसलिए जैसे ही डेटा उपलब्ध होता है पूरी chain स्वचालित रूप से चलती है, बिना तुम्हें प्रत्येक चरण को हाथ से orchestrate किए।

ध्यान दो कि यह वर्कफ़्लो मॉड्यूल का उपयोग करता है: प्रत्येक प्रोसेस `modules/` के अंतर्गत अपनी खुद की फ़ाइल में परिभाषित है, और `main.nf` उन्हें inline परिभाषित करने के बजाय `include` स्टेटमेंट से import करता है।
यह प्रत्येक प्रोसेस को कोड duplicate किए बिना कई वर्कफ़्लो में पुन: उपयोग योग्य बनाता है। अधिक जानने के लिए, नीचे कोड exploration सेक्शन देखो।

### 3.1. वर्कफ़्लो चलाना

अपने टर्मिनल में निम्नलिखित कमांड चलाओ।

```bash
nextflow run main.nf --input data/greetings.csv
```

`character` पैरामीटर `nextflow.config` में डिफ़ॉल्ट रूप से `turkey` है, इसलिए ASCII art एक turkey का उपयोग करती है जब तक तुम इसे override न करो (कोशिश करो `--character tux` जोड़ना)।

??? success "कमांड आउटपुट"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nostalgic_brahmagupta] revision: ce74f81996

    executor >  local (8)
    [56/8499f6] sayHello (3)       | 3 of 3 ✔
    [cc/0ee42a] convertToUpper (3) | 3 of 3 ✔
    [eb/0f2e24] collectGreetings   | 1 of 1 ✔
    [b5/34e07f] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

चार प्रोसेस चलीं, लेकिन समान संख्या में नहीं।
`sayHello` और `convertToUpper` प्रत्येक इनपुट के लिए एक बार चले (3 of 3): प्रत्येक अभिवादन को अपने आप लिखा और uppercase किया जाना ज़रूरी है।
`collectGreetings` और `cowpy` केवल एक बार चले (1 of 1): अभिवादन मर्ज करना और ASCII art बनाना तभी समझ में आता है जब हर एक व्यक्तिगत परिणाम तैयार हो।
यह fan-out-then-fan-in आकार, कई समानांतर कार्य जो कम संख्या में downstream कार्यों में फीड होते हैं, वास्तविक पाइपलाइन में सामान्य है।

Nextflow अगला शुरू करने से पहले पूरे चरण के समाप्त होने का इंतज़ार नहीं करता।
जैसे ही एक `sayHello` आउटपुट तैयार होता है, मेल खाने वाला `convertToUpper` कार्य शुरू हो सकता है, इसलिए विभिन्न प्रोसेस के कार्य सख्त batches के बजाय एक साथ चलते हैं।
`collectGreetings` और `cowpy` को इंतज़ार करना पड़ता है, क्योंकि उनमें से प्रत्येक पहले हर upstream परिणाम के उपलब्ध होने पर निर्भर करता है।

`results` डायरेक्टरी उस fan-in को दर्शाती है, साथ ही पाइपलाइन लेखक ने जो प्रकाशित करना चुना और कहाँ: 1.4 के कोड वॉकथ्रू से `output` ब्लॉक याद करो, जो इस संरचना को परिभाषित करता है।

```console title="results/"
results
└── batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
```

शीर्ष-स्तरीय डायरेक्टरी का नाम `batch` पैरामीटर के नाम पर है, जो डिफ़ॉल्ट रूप से `batch` है; तुम इसे बाद के अभ्यासों में बदलते देखोगे।

ASCII art फ़ाइल के लिए `cowpy-COLLECTED-batch-output.txt` देखो।

??? abstract "फ़ाइल सामग्री"

    ```console title="results/batch/cowpy-COLLECTED-batch-output.txt"
     _________
    / HELLO   \
    | BONJOUR |
    \ HOLA    /
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

जैसे [2.1](#21-run-the-workflow) में, इन 8 कार्य executions में से प्रत्येक, सभी चार प्रोसेस में, `work/` के अंतर्गत अपनी खुद की डायरेक्टरी पाती है, जो दूसरों से पूरी तरह अलग होती है।
`collectGreetings` इस बात का एक अच्छा उदाहरण है कि यह क्यों मायने रखता है: यह तीन अलग-अलग कार्य डायरेक्टरी में रहने वाले सभी तीन `convertToUpper` कार्यों के आउटपुट पर निर्भर करता है, इसलिए Nextflow उन फ़ाइलों के symlinks को `collectGreetings` की अपनी डायरेक्टरी के अंदर stage करता है बजाय इसके कि वह सीधे अपने upstream कार्यों की डायरेक्टरी से पढ़े:

```console title="work/eb/0f2e24.../"
COLLECTED-batch-output.txt
UPPER-Bonjour-output.txt -> ../../69/e6c057.../UPPER-Bonjour-output.txt
UPPER-Hello-output.txt   -> ../../cc/ba7d19.../UPPER-Hello-output.txt
UPPER-Hola-output.txt    -> ../../cc/0ee42a.../UPPER-Hola-output.txt
batch-report.txt
.command.sh
```

प्रत्येक कार्य केवल उन विशिष्ट फ़ाइलों को देखता है जिनकी उसे ज़रूरत है, चाहे वे कहीं से भी आई हों, और कभी भी किसी दूसरे कार्य की डायरेक्टरी की आंतरिक सामग्री नहीं देखता।
पूरी पाइपलाइन में, [2.1](#21-run-the-workflow) में एकल प्रोसेस के साथ तुमने जो वही अलगाव देखा, वही Nextflow को हर प्रोसेस के हर कार्य को एक साथ, सुरक्षित रूप से चलाने देता है।

!!! note "नोट"

    `cowpy` चरण स्थानीय रूप से इंस्टॉल सॉफ़्टवेयर पर निर्भर करने के बजाय एक Docker कंटेनर के अंदर चलता है।
    एक कंटेनर एक application को उसे चलाने के लिए ज़रूरी हर चीज़ के साथ package करता है, इसलिए तुम्हें खुद dependencies इंस्टॉल और प्रबंधित नहीं करनी पड़तीं, और पाइपलाइन किसी भी मशीन पर जो कंटेनर चला सकती है उसी तरह व्यवहार करती है।
    Nextflow कंटेनर के विकल्प के रूप में Conda को भी सपोर्ट करता है; उनके बीच स्विच करने के तरीके के लिए [भाग 2](./02_configure_pipeline.md) देखो।

### 3.2. वैकल्पिक: कोड वॉकथ्रू

अगर तुम सिर्फ पाइपलाइन चलाना चाहते हो तो कोड समझना ज़रूरी नहीं है, लेकिन अगर उत्सुक हो, तो एक नज़र डालना उचित है।

??? optional "इस अभ्यास से जुड़े कोड को एक्सप्लोर करने के लिए क्लिक करो"

    ### एक चरण से दूसरे चरण में डेटा कैसे प्रवाहित होता है

    प्रत्येक प्रोसेस अपना आउटपुट चैनल अगले को पास करती है:

    ```groovy title="main.nf" linenums="19" hl_lines="7 8 9"
        main:
        // CSV फ़ाइल से इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        sayHello(greeting_ch)
        convertToUpper(sayHello.out)
        collectGreetings(convertToUpper.out.collect(), params.batch)
        cowpy(collectGreetings.out.outfile, params.character)
    ```

    पैटर्न `processName.out` एक प्रोसेस के आउटपुट चैनल को संदर्भित करता है।

    `.collect()` ऑपरेटर `convertToUpper` के सभी व्यक्तिगत आउटपुट को `collectGreetings` को पास करने से पहले एक single channel item में इकट्ठा करता है।

    ### प्रोसेस मॉड्यूल का उपयोग करना

    `main.nf` सीधे कोई प्रोसेस कोड परिभाषित नहीं करता।
    इसके बजाय, यह `modules/` के अंतर्गत अपनी खुद की फ़ाइल से प्रत्येक प्रोसेस import करता है:

    ```groovy title="main.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

    प्रत्येक मॉड्यूल फ़ाइल में एक एकल प्रोसेस परिभाषा होती है, जो [1.4](#14-optional-code-walkthrough) में `sayHello` मॉड्यूल की तरह ही संरचित होती है।
    प्रोसेस को अलग फ़ाइलों में रखने से वे कोड duplicate किए बिना कई वर्कफ़्लो में पुन: उपयोग योग्य बनती हैं।

    <figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
    </figure>

    ### Containerized सॉफ़्टवेयर का उपयोग करना

    `cowpy` प्रोसेस अपनी मॉड्यूल फ़ाइल में निर्दिष्ट एक Docker कंटेनर के अंदर चलती है:

    ```groovy title="modules/cowpy.nf" linenums="2" hl_lines="3"
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

    Nextflow स्वचालित रूप से image pull करता है, script को कंटेनर के अंदर चलाता है, और बाद में cleanup करता है।
    Docker इस प्रोजेक्ट के लिए `nextflow.config` में enabled है:

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

    यह एकल लाइन पाइपलाइन में किसी भी प्रोसेस के लिए Docker enable करती है जिसमें कंटेनर निर्दिष्ट है।

### सारांश

तुमने एक पूर्ण बहु-चरण पाइपलाइन चलाई जो एक containerized टूल का उपयोग करके समानांतर में कई इनपुट प्रोसेस करती है।

### आगे क्या है?

[भाग 2](./02_configure_pipeline.md) पर जाओ, जहाँ तुम सीखोगे कि `nextflow.config` का उपयोग करके पाइपलाइन व्यवहार को कैसे कॉन्फ़िगर करें।

---

## सारांश

इस भाग में तुमने सीखा:

- Nextflow वर्कफ़्लो चलाना और उसके आउटपुट खोजना
- `work/` डायरेक्टरी और उसकी लॉग फ़ाइलें एक्सप्लोर करना
- CSV फ़ाइल से कई इनपुट को समानांतर में प्रोसेस करना
- नए इनपुट जोड़ते समय पूर्ण किए गए काम को छोड़ने के लिए `-resume` का उपयोग करना
- एक बहु-चरण पाइपलाइन चलाना जो एक containerized टूल का उपयोग करती है
