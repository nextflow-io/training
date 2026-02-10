# भाग 3: Hello Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [पूरी प्लेलिस्ट](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) Nextflow YouTube चैनल पर देखें।

:green_book: वीडियो ट्रांसक्रिप्ट [यहाँ](./transcripts/03_hello_workflow.md) उपलब्ध है।
///

अधिकांश वास्तविक दुनिया के workflows में एक से अधिक स्टेप शामिल होते हैं।
इस प्रशिक्षण मॉड्यूल में, तुम सीखोगे कि processes को एक multi-step workflow में कैसे जोड़ा जाए।

यह तुम्हें निम्नलिखित को प्राप्त करने का Nextflow तरीका सिखाएगा:

1. एक process से दूसरे में डेटा का प्रवाह बनाना
2. कई process calls से outputs को एकत्र करके एक single process call में भेजना
3. एक process को अतिरिक्त पैरामीटर पास करना
4. एक process से निकलने वाले कई outputs को संभालना

प्रदर्शन के लिए, हम भाग 1 और 2 के domain-agnostic Hello World उदाहरण पर निर्माण करना जारी रखेंगे।
इस बार, हम अपने workflow में निम्नलिखित बदलाव करेंगे ताकि यह बेहतर तरीके से दर्शाए कि लोग वास्तविक workflows कैसे बनाते हैं:

1. एक दूसरा स्टेप जोड़ें जो greeting को uppercase में बदल दे।
2. एक तीसरा स्टेप जोड़ें जो सभी transformed greetings को एकत्र करके एक single फ़ाइल में लिखे।
3. अंतिम output फ़ाइल को नाम देने के लिए एक पैरामीटर जोड़ें और उसे collection स्टेप में secondary input के रूप में पास करें।
4. Collection स्टेप को यह भी रिपोर्ट करवाएं कि क्या प्रोसेस किया गया, इसके बारे में एक सरल आंकड़ा।

??? info "इस खंड से कैसे शुरू करें"

    कोर्स का यह सेक्शन मानता है कि तुमने [Hello Nextflow](./index.md) कोर्स के भाग 1-2 पूरे कर लिए हैं, लेकिन अगर तुम उन सेक्शन में कवर की गई बेसिक्स के साथ सहज हो, तो तुम यहाँ से बिना कुछ विशेष किए शुरू कर सकते हो।

---

## 0. वार्मअप: `hello-workflow.nf` चलाएं

हम शुरुआती बिंदु के रूप में workflow स्क्रिप्ट `hello-workflow.nf` का उपयोग करने जा रहे हैं।
यह इस प्रशिक्षण कोर्स के भाग 2 को पूरा करके बनाई गई स्क्रिप्ट के बराबर है, सिवाय इसके कि हमने `view()` statements हटा दिए हैं और output destination बदल दिया है:

```groovy title="hello-workflow.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
```

यह डायग्राम workflow के वर्तमान संचालन को सारांशित करता है।
यह परिचित दिखना चाहिए, सिवाय इसके कि अब हम स्पष्ट रूप से दिखा रहे हैं कि process के outputs को एक channel में पैकेज किया जाता है, ठीक वैसे ही जैसे inputs थे।
हम एक मिनट में उस output channel का अच्छा उपयोग करने जा रहे हैं।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-workflow-channels.svg"
</figure>

बस यह सुनिश्चित करने के लिए कि सब कुछ काम कर रहा है, कोई भी बदलाव करने से पहले स्क्रिप्ट को एक बार चलाएं:

```bash
nextflow run hello-workflow.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [admiring_lamarr] DSL2 - revision: 4d4053520d

    executor >  local (3)
    [b1/5826b5] process > sayHello (2) [100%] 3 of 3 ✔
    ```

पहले की तरह, तुम्हें output फ़ाइलें `output` ब्लॉक में निर्दिष्ट स्थान पर मिलेंगी।
इस अध्याय के लिए, यह `results/hello_workflow/` के अंदर है।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results/hello_workflow
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    └── Holà-output.txt
    ```

अगर यह तुम्हारे लिए काम कर गया, तो तुम multi-step workflow को असेंबल करना सीखने के लिए तैयार हो।

---

## 1. Workflow में दूसरा स्टेप जोड़ें

हम प्रत्येक greeting को uppercase में बदलने के लिए एक स्टेप जोड़ने जा रहे हैं।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

इसके लिए, हमें तीन काम करने होंगे:

- वह कमांड परिभाषित करें जिसका उपयोग हम uppercase conversion के लिए करने जा रहे हैं।
- एक नया process लिखें जो uppercasing कमांड को wrap करे।
- Workflow ब्लॉक में नए process को कॉल करें और इसे `sayHello()` process के output को input के रूप में लेने के लिए सेट करें।

### 1.1. Uppercasing कमांड परिभाषित करें और टर्मिनल में टेस्ट करें

Greetings को uppercase में बदलने के लिए, हम `tr` नामक एक क्लासिक UNIX टूल का उपयोग करने जा रहे हैं जो 'text replacement' के लिए है, निम्नलिखित सिंटैक्स के साथ:

```bash title="Syntax"
tr '[a-z]' '[A-Z]'
```

यह एक बहुत ही सरल text replacement one-liner है जो accented letters को ध्यान में नहीं रखता, इसलिए उदाहरण के लिए 'Holà' 'HOLà' बन जाएगा, लेकिन यह Nextflow अवधारणाओं को प्रदर्शित करने के लिए काफी अच्छा काम करेगा और यही मायने रखता है।

इसे टेस्ट करने के लिए, हम `echo 'Hello World'` कमांड चला सकते हैं और इसके output को `tr` कमांड में pipe कर सकते हैं:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

Output एक टेक्स्ट फ़ाइल है जिसे `UPPER-output.txt` कहा जाता है जिसमें `Hello World` स्ट्रिंग का uppercase संस्करण होता है।

??? abstract "फ़ाइल सामग्री"

    ```console title="UPPER-output.txt"
    HELLO WORLD
    ```

मूल रूप से यही हम अपने workflow के साथ करने की कोशिश करने जा रहे हैं।

### 1.2. Uppercasing स्टेप को Nextflow process के रूप में लिखें

हम अपने नए process को पहले वाले पर मॉडल कर सकते हैं, क्योंकि हम सभी समान components का उपयोग करना चाहते हैं।

Workflow स्क्रिप्ट में निम्नलिखित process परिभाषा जोड़ें, पहले वाले के ठीक नीचे:

```groovy title="hello-workflow.nf" linenums="20"
/*
 * Greeting को uppercase में बदलने के लिए एक text replacement टूल का उपयोग करें
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

इसमें, हम दूसरी output फ़ाइलनाम को input फ़ाइलनाम के आधार पर compose करते हैं, वैसे ही जैसे हमने मूल रूप से पहले process के output के लिए किया था।

### 1.3. Workflow ब्लॉक में नए process को कॉल करें

अब हमें Nextflow को बताना होगा कि वास्तव में उस process को कॉल करे जिसे हमने अभी परिभाषित किया है।

Workflow ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="44" hl_lines="10-11"
    workflow {

        main:
        // CSV फ़ाइल से inputs के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक greeting emit करें
        sayHello(greeting_ch)
        // greeting को uppercase में बदलें
        convertToUpper()

        publish:
        first_output = sayHello.out
    }
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="44"
    workflow {

        main:
        // CSV फ़ाइल से inputs के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक greeting emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

यह अभी तक functional नहीं है क्योंकि हमने यह निर्दिष्ट नहीं किया है कि `convertToUpper()` process को क्या input होना चाहिए।

### 1.4. पहले process के output को दूसरे process में पास करें

अब हमें `sayHello()` process के output को `convertToUpper()` process में flow करना होगा।

सुविधाजनक रूप से, Nextflow स्वचालित रूप से एक process के output को एक channel में पैकेज करता है, जैसा कि वार्मअप सेक्शन में डायग्राम में दिखाया गया है।
हम एक process के output channel को `<process>.out` के रूप में refer कर सकते हैं।

तो `sayHello` process का output एक channel है जिसे `sayHello.out` कहा जाता है, जिसे हम सीधे `convertToUpper()` की कॉल में plug कर सकते हैं।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep-connector.svg"
</figure>

Workflow ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // greeting को uppercase में बदलें
        convertToUpper(sayHello.out)
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // greeting को uppercase में बदलें
        convertToUpper()
    ```

इस तरह के एक सरल मामले के लिए (एक output से एक input), दो processes को जोड़ने के लिए हमें बस इतना ही करना होगा!

### 1.5. Workflow output publishing सेट करें

अंत में, आइए workflow outputs को अपडेट करें ताकि दूसरे process के परिणाम भी publish हों।

#### 1.5.1. `workflow` ब्लॉक के `publish:` सेक्शन को अपडेट करें

`workflow` ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

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

तर्क पहले जैसा ही है।

#### 1.5.2. `output` ब्लॉक को अपडेट करें

`output` ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

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

एक बार फिर, तर्क पहले जैसा ही है।

यह दिखाता है कि तुम बहुत granular स्तर पर output सेटिंग्स को नियंत्रित कर सकते हो, हर व्यक्तिगत output के लिए।
यह देखने के लिए स्वतंत्र महसूस करो कि क्या होता है अगर तुम processes में से एक के लिए paths या publish mode बदलते हो।

बेशक, इसका मतलब है कि हम यहाँ कुछ जानकारी दोहरा रहे हैं, जो असुविधाजनक हो सकती है अगर हम सभी outputs के लिए location को एक ही तरीके से अपडेट करना चाहते हैं।
कोर्स में बाद में, तुम सीखोगे कि इन सेटिंग्स को कई outputs के लिए structured तरीके से कैसे configure किया जाए।

### 1.6. `-resume` के साथ workflow चलाएं

आइए `-resume` फ्लैग का उपयोग करके इसे टेस्ट करें, क्योंकि हमने पहले ही workflow के पहले स्टेप को सफलतापूर्वक चलाया है।

```bash
nextflow run hello-workflow.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [high_cantor] DSL2 - revision: d746983511

    executor >  local (3)
    [ab/816321] process > sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [e0/ecf81b] process > convertToUpper (3) [100%] 3 of 3 ✔
    ```

Console output में अब एक अतिरिक्त लाइन है जो हमने अभी जोड़े गए नए process से मेल खाती है।

तुम्हें outputs `results/hello_workflow` डायरेक्टरी में मिलेंगे जैसा कि `output` ब्लॉक में सेट किया गया है।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results/hello_workflow/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

यह सुविधाजनक है! लेकिन दूसरे process की calls में से एक की work डायरेक्टरी के अंदर देखना अभी भी उचित है।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    work/e0/ecf81b4cacc648b9b994218d5b29d7/
    ├── Holà-output.txt -> /workspaces/training/hello-nextflow/work/ab/81632178cd37e9e815959278808819/Holà-output.txt
    └── UPPER-Holà-output.txt
    ```

ध्यान दें कि दो `*-output` फ़ाइलें हैं: पहले process का output साथ ही दूसरे का output।

पहले process का output वहाँ है क्योंकि Nextflow ने इसे वहाँ **stage** किया ताकि execution के लिए आवश्यक सब कुछ एक ही subdirectory में हो।

हालाँकि, यह वास्तव में एक symbolic link है जो पहले process call की subdirectory में मूल फ़ाइल की ओर इशारा करता है।
डिफ़ॉल्ट रूप से, जब एक single मशीन पर चल रहे हों जैसा कि हम यहाँ कर रहे हैं, Nextflow input और intermediate फ़ाइलों को stage करने के लिए copies के बजाय symbolic links का उपयोग करता है।

अब, आगे बढ़ने से पहले, सोचो कि हमने बस `sayHello` के output को `convertToUpper` के input से कैसे जोड़ा और दोनों processes को series में चलाया जा सकता था।
Nextflow ने व्यक्तिगत input और output फ़ाइलों को संभालने और उन्हें दो कमांड्स के बीच पास करने का कठिन काम हमारे लिए किया।

यह एक कारण है कि Nextflow channels इतने शक्तिशाली हैं: वे workflow स्टेप्स को एक साथ जोड़ने में शामिल busywork का ध्यान रखते हैं।

### सारांश

तुम जानते हो कि एक स्टेप के output को अगले स्टेप के input के रूप में प्रदान करके processes को कैसे chain किया जाए।

### आगे क्या है?

सीखो कि batched process calls से outputs को कैसे एकत्र किया जाए और उन्हें एक single process में कैसे feed किया जाए।

---

## 2. सभी greetings को एकत्र करने के लिए तीसरा स्टेप जोड़ें

जब हम एक channel में elements में से प्रत्येक पर transformation लागू करने के लिए एक process का उपयोग करते हैं, जैसा कि हम यहाँ कई greetings के लिए कर रहे हैं, हम कभी-कभी उस process के output channel से elements को एकत्र करना चाहते हैं, और उन्हें एक अन्य process में feed करना चाहते हैं जो किसी प्रकार का विश्लेषण या summation करता है।

प्रदर्शन के लिए, हम अपनी pipeline में एक नया स्टेप जोड़ेंगे जो `convertToUpper` process द्वारा उत्पादित सभी uppercase greetings को एकत्र करता है और उन्हें एक single फ़ाइल में लिखता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

आश्चर्य को खराब न करने के लिए, लेकिन इसमें एक बहुत उपयोगी operator शामिल होने वाला है।

### 2.1. Collection कमांड परिभाषित करें और टर्मिनल में टेस्ट करें

हमारे workflow में जोड़ने वाला collection स्टेप कई uppercased greetings को एक single फ़ाइल में concatenate करने के लिए `cat` कमांड का उपयोग करेगा।

आइए टर्मिनल में कमांड को अकेले चलाएं यह सत्यापित करने के लिए कि यह अपेक्षित रूप से काम करता है, ठीक वैसे ही जैसे हमने पहले किया है।

अपने टर्मिनल में निम्नलिखित चलाएं:

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

Output एक टेक्स्ट फ़ाइल है जिसे `COLLECTED-output.txt` कहा जाता है जिसमें मूल greetings के uppercase संस्करण होते हैं।

??? abstract "फ़ाइल सामग्री"

    ```console title="COLLECTED-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

यही वह परिणाम है जिसे हम अपने workflow के साथ प्राप्त करना चाहते हैं।

### 2.2. Collection स्टेप करने के लिए एक नया process बनाएं

आइए एक नया process बनाएं और इसे `collectGreetings()` कहें।
हम इसे लिखना शुरू कर सकते हैं जो हमने पहले देखा है उसके आधार पर।

#### 2.2.1. Process के 'स्पष्ट' भागों को लिखें

Workflow स्क्रिप्ट में निम्नलिखित process परिभाषा जोड़ें:

```groovy title="hello-workflow.nf" linenums="37"
/*
 * Uppercase greetings को एक single output फ़ाइल में एकत्र करें
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

यह वह है जो हम अब तक सीखे गए के आधार पर आत्मविश्वास के साथ लिख सकते हैं।
लेकिन यह functional नहीं है!
यह input परिभाषा(ओं) और script कमांड के पहले आधे हिस्से को छोड़ देता है क्योंकि हमें यह पता लगाना होगा कि इसे कैसे लिखा जाए।

#### 2.2.2. `collectGreetings()` के inputs परिभाषित करें

हमें `convertToUpper()` process की सभी calls से greetings को एकत्र करना होगा।
हम क्या जानते हैं कि हम workflow के पिछले स्टेप से क्या प्राप्त कर सकते हैं?

`convertToUpper()` द्वारा output किया गया channel उन व्यक्तिगत फ़ाइलों के paths को contain करेगा जिनमें uppercased greetings हैं।
यह एक input slot के बराबर है; आइए इसे सरलता के लिए `input_files` कहें।

Process ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

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

ध्यान दें कि हम `path` prefix का उपयोग करते हैं भले ही हम इसमें कई फ़ाइलें होने की उम्मीद करते हैं।

#### 2.2.3. Concatenation कमांड compose करें

यहाँ चीजें थोड़ी मुश्किल हो सकती हैं, क्योंकि हमें input फ़ाइलों की एक arbitrary संख्या को संभालने में सक्षम होना चाहिए।
विशेष रूप से, हम कमांड को पहले से नहीं लिख सकते, इसलिए हमें Nextflow को बताना होगा कि runtime पर इसे कैसे compose किया जाए जो inputs process में flow होते हैं उसके आधार पर।

दूसरे शब्दों में, अगर हमारे पास एक input channel है जिसमें element `[file1.txt, file2.txt, file3.txt]` है, तो हमें Nextflow को इसे `cat file1.txt file2.txt file3.txt` में बदलने की आवश्यकता है।

सौभाग्य से, Nextflow हमारे लिए ऐसा करने में काफी खुश है अगर हम बस script कमांड में `cat ${input_files}` लिखते हैं।

Process ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

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

सिद्धांत रूप में यह input फ़ाइलों की किसी भी arbitrary संख्या को संभाल सकता है।

!!! tip "सुझाव"

    कुछ command-line टूल्स को प्रत्येक input फ़ाइल के लिए एक argument (जैसे `-input`) प्रदान करने की आवश्यकता होती है।
    उस स्थिति में, हमें कमांड को compose करने के लिए थोड़ा अतिरिक्त काम करना होगा।
    तुम इसका एक उदाहरण [Nextflow for Genomics](../../nf4_science/genomics/) प्रशिक्षण कोर्स में देख सकते हो।

### 2.3. Workflow में collection स्टेप जोड़ें

अब हमें बस uppercasing स्टेप के output पर collection process को कॉल करना चाहिए।
वह भी एक channel है, जिसे `convertToUpper.out` कहा जाता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-connector.svg"
</figure>

#### 2.3.1. Process calls को connect करें

Workflow ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // greeting को uppercase में बदलें
        convertToUpper(sayHello.out)

        // सभी greetings को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out)
    }
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="75"
        // greeting को uppercase में बदलें
        convertToUpper(sayHello.out)
    }
    ```

यह `convertToUpper()` के output को `collectGreetings()` के input से जोड़ता है।

#### 2.3.2. `-resume` के साथ workflow चलाएं

आइए इसे आज़माएं।

```bash
nextflow run hello-workflow.nf -resume
```

??? success "कमांड आउटपुट"

    ```console hl_lines="8"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [mad_gilbert] DSL2 - revision: 6acfd5e28d

    executor >  local (3)
    [79/33b2f0] sayHello (2)         | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (3)   | 3 of 3, cached: 3 ✔
    [47/50fe4a] collectGreetings (1) | 3 of 3 ✔
    ```

यह सफलतापूर्वक चलता है, तीसरे स्टेप सहित।

हालाँकि, आखिरी लाइन पर `collectGreetings()` के लिए calls की संख्या देखो।
हम केवल एक की उम्मीद कर रहे थे, लेकिन तीन हैं।

अब अंतिम output फ़ाइल की सामग्री पर एक नज़र डालो।

??? abstract "फ़ाइल सामग्री"

    ```console title="results/COLLECTED-output.txt"
    Holà
    ```

ओह नहीं। Collection स्टेप प्रत्येक greeting पर व्यक्तिगत रूप से चलाया गया, जो हम नहीं चाहते थे।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-no-operator.svg"
</figure>

हमें Nextflow को स्पष्ट रूप से बताने के लिए कुछ करना होगा कि हम चाहते हैं कि तीसरा स्टेप `convertToUpper()` द्वारा output किए गए channel में सभी elements पर चले।

### 2.4. Greetings को एक single input में एकत्र करने के लिए एक operator का उपयोग करें

हाँ, एक बार फिर हमारी समस्या का उत्तर एक operator है।

विशेष रूप से, हम उपयुक्त नाम वाले [`collect()`](https://nextflow.io/docs/latest/reference/operator.html#collect) operator का उपयोग करने जा रहे हैं।

#### 2.4.1. `collect()` operator जोड़ें

इस बार यह थोड़ा अलग दिखने वाला है क्योंकि हम channel factory के संदर्भ में operator नहीं जोड़ रहे हैं; हम इसे एक output channel में जोड़ रहे हैं।

हम `convertToUpper.out` लेते हैं और `collect()` operator जोड़ते हैं, जो हमें `convertToUpper.out.collect()` देता है।
हम इसे सीधे `collectGreetings()` process call में plug कर सकते हैं।

Workflow ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // सभी greetings को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // सभी greetings को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out)
    }
    ```

#### 2.4.2. कुछ `view()` statements जोड़ें

आइए channel सामग्री की before और after states को visualize करने के लिए कुछ `view()` statements भी शामिल करें।

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // सभी greetings को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out.collect())

        // वैकल्पिक view statements
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    }
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="73"
        // सभी greetings को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out.collect())
    }
    ```

`view()` statements कहीं भी जा सकते हैं जहाँ तुम चाहते हो; हमने उन्हें readability के लिए call के ठीक बाद रखा।

#### 2.4.3. `-resume` के साथ फिर से workflow चलाएं

आइए इसे आज़माएं:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "कमांड आउटपुट"

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

यह सफलतापूर्वक चलता है, हालाँकि log output इससे थोड़ा अधिक messy दिख सकता है (हमने readability के लिए इसे साफ किया)।

इस बार तीसरा स्टेप केवल एक बार कॉल किया गया!
`view()` statements के output को देखते हुए, हम निम्नलिखित देखते हैं:

- तीन `Before collect:` statements, प्रत्येक greeting के लिए एक: उस बिंदु पर फ़ाइल paths channel में व्यक्तिगत items हैं।
- एक single `After collect:` statement: तीन फ़ाइल paths अब एक single element में पैकेज किए गए हैं।

हम इसे निम्नलिखित डायग्राम के साथ सारांशित कर सकते हैं:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-WITH-operator.svg"
</figure>

अंत में, तुम output फ़ाइल की सामग्री पर एक नज़र डाल सकते हो यह संतुष्ट करने के लिए कि सब कुछ सही तरीके से काम किया।

??? abstract "फ़ाइल सामग्री"

    ```console title="results/COLLECTED-output.txt"
    BONJOUR
    HELLO
    HOLà
    ```

इस बार हमारे पास अंतिम output फ़ाइल में सभी तीन greetings हैं। सफलता!

!!! note "नोट"

    अगर तुम इसे `-resume` के बिना कई बार चलाते हो, तो तुम देखोगे कि greetings का क्रम एक run से दूसरे में बदलता है।
    यह तुम्हें दिखाता है कि जिस क्रम में elements process calls के माध्यम से flow होते हैं वह consistent होने की गारंटी नहीं है।

#### 2.4.4. Readability के लिए `view()` statements हटाएं

अगले सेक्शन पर जाने से पहले, हम अनुशंसा करते हैं कि तुम console output को cluttering से बचाने के लिए `view()` statements को delete कर दो।

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="73"
        // सभी greetings को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out.collect())
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // सभी greetings को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out.collect())

        // वैकल्पिक view statements
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    ```

यह मूल रूप से बिंदु 2.4.2 से उल्टा ऑपरेशन है।

### सारांश

तुम जानते हो कि process calls के एक batch से outputs को कैसे एकत्र किया जाए और उन्हें एक joint विश्लेषण या summation स्टेप में कैसे feed किया जाए।

संक्षेप में, यह वह है जो तुमने अब तक बनाया है:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

### आगे क्या है?

सीखो कि एक process को एक से अधिक input कैसे पास किया जाए।

---

## 3. एक process को अतिरिक्त पैरामीटर पास करें

हम अंतिम output फ़ाइल को कुछ विशिष्ट नाम देने में सक्षम होना चाहते हैं ताकि अंतिम परिणामों को overwrite किए बिना greetings के बाद के batches को process कर सकें।

इसके लिए, हम workflow में निम्नलिखित refinements करने जा रहे हैं:

- Collector process को output फ़ाइल के लिए एक user-defined नाम स्वीकार करने के लिए संशोधित करें (`batch_name`)
- Workflow में एक command-line पैरामीटर जोड़ें (`--batch`) और इसे collector process में पास करें

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-batch.svg"
</figure>

### 3.1. Collector process को संशोधित करें

हमें अतिरिक्त input को declare करना होगा और इसे output फ़ाइल नाम में integrate करना होगा।

#### 3.1.1. अतिरिक्त input declare करें

अच्छी खबर: हम process परिभाषा में जितने चाहें उतने input variables declare कर सकते हैं।
आइए इसे `batch_name` कहें।

Process ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

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

तुम अपने processes को जितने चाहो उतने inputs की उम्मीद करने के लिए सेट कर सकते हो।
अभी, ये सभी required inputs होने के लिए सेट हैं; workflow के काम करने के लिए तुम्हें एक value प्रदान _करनी होगी_।

तुम अपनी Nextflow यात्रा में बाद में सीखोगे कि required बनाम optional inputs को कैसे manage किया जाए।

#### 3.1.2. Output फ़ाइल नाम में `batch_name` variable का उपयोग करें

हम output फ़ाइल नाम में variable को उसी तरह insert कर सकते हैं जैसे हमने पहले dynamic फ़ाइल नाम compose किए हैं।

Process ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

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

यह process को workflow के अंतिम output के लिए एक विशिष्ट फ़ाइलनाम generate करने के लिए `batch_name` value का उपयोग करने के लिए सेट करता है।

### 3.2. एक `batch` command-line पैरामीटर जोड़ें

अब हमें `batch_name` के लिए value supply करने और इसे process call में feed करने का एक तरीका चाहिए।

#### 3.2.1. पैरामीटर सेट करने के लिए `params` का उपयोग करें

तुम पहले से जानते हो कि CLI पैरामीटर declare करने के लिए `params` सिस्टम का उपयोग कैसे किया जाए।
आइए इसका उपयोग एक `batch` पैरामीटर declare करने के लिए करें (एक डिफ़ॉल्ट value के साथ क्योंकि हम आलसी हैं)।

Pipeline पैरामीटर सेक्शन में, निम्नलिखित कोड परिवर्तन करें:

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

ठीक वैसे ही जैसे हमने `--input` के लिए प्रदर्शित किया, तुम command line पर `--batch` के साथ एक value निर्दिष्ट करके उस डिफ़ॉल्ट value को override कर सकते हो।

#### 3.2.2. Process को `batch` पैरामीटर पास करें

Process को पैरामीटर की value प्रदान करने के लिए, हमें इसे process call में जोड़ना होगा।

Workflow ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // सभी greetings को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "पहले"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // सभी greetings को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out.collect())
    ```

तुम देखते हो कि एक process को कई inputs प्रदान करने के लिए, तुम बस उन्हें call parentheses में list करते हो, commas से अलग।

!!! warning "चेतावनी"

    तुम्हें process को inputs को बिल्कुल उसी क्रम में प्रदान करना होगा जैसे वे process के input definition ब्लॉक में listed हैं।

### 3.3. Workflow चलाएं

आइए इसे command line पर एक batch नाम के साथ चलाने की कोशिश करें।

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [confident_rutherford] DSL2 - revision: bc58af409c

    executor >  local (1)
    [79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [b5/f19efe] collectGreetings   | 1 of 1 ✔
    ```

यह सफलतापूर्वक चलता है और वांछित output उत्पन्न करता है:

??? abstract "फ़ाइल सामग्री"

    ```console title="results/COLLECTED-trio-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

अब, जब तक हम पैरामीटर को उचित रूप से निर्दिष्ट करते हैं, inputs के अन्य batches पर बाद के runs पिछले परिणामों को clobber नहीं करेंगे।

### सारांश

तुम जानते हो कि एक process को एक से अधिक input कैसे पास किया जाए।

### आगे क्या है?

सीखो कि कई outputs को कैसे emit किया जाए और उन्हें सुविधाजनक तरीके से कैसे संभाला जाए।

---

## 4. Collector स्टेप में एक output जोड़ें

अब तक हम ऐसे processes का उपयोग कर रहे हैं जो प्रत्येक केवल एक output उत्पन्न करते थे।
हम उनके respective outputs को बहुत सुविधाजनक तरीके से `<process>.out` सिंटैक्स का उपयोग करके access कर सकते थे, जिसका उपयोग हमने अगले process को output पास करने के संदर्भ में (जैसे `convertToUpper(sayHello.out)`) और `publish:` सेक्शन के संदर्भ में (जैसे `first_output = sayHello.out`) दोनों में किया।

जब एक process एक से अधिक उत्पन्न करता है तो क्या होता है?
हम कई outputs को कैसे संभालते हैं?
क्या हम एक विशिष्ट output को select और उपयोग कर सकते हैं?

सभी उत्कृष्ट प्रश्न, और संक्षिप्त उत्तर है हाँ हम कर सकते हैं!

कई outputs को अलग channels में पैकेज किया जाएगा।
हम या तो उन output channels को नाम देना चुन सकते हैं, जो बाद में उन्हें व्यक्तिगत रूप से refer करना आसान बनाता है, या हम उन्हें index द्वारा refer कर सकते हैं।

प्रदर्शन उद्देश्यों के लिए, मान लें कि हम greetings की संख्या को गिनना चाहते हैं जो inputs के एक दिए गए batch के लिए एकत्र किए जा रहे हैं और इसे एक फ़ाइल में रिपोर्ट करना चाहते हैं।

### 4.1. Greetings की संख्या को गिनने और output करने के लिए process को संशोधित करें

इसके लिए process परिभाषा में दो प्रमुख परिवर्तनों की आवश्यकता होगी: हमें greetings को गिनने और एक report फ़ाइल लिखने का एक तरीका चाहिए, फिर हमें उस report फ़ाइल को process के `output` ब्लॉक में जोड़ना होगा।

#### 4.1.1. एकत्रित greetings की संख्या गिनें

सुविधाजनक रूप से, Nextflow हमें process परिभाषा के `script:` ब्लॉक में arbitrary कोड जोड़ने देता है, जो इस तरह की चीजें करने के लिए वास्तव में काम आता है।

इसका मतलब है कि हम `input_files` array में फ़ाइलों की संख्या प्राप्त करने के लिए Nextflow के built-in `size()` function का उपयोग कर सकते हैं, और परिणाम को एक `echo` कमांड के साथ फ़ाइल में लिख सकते हैं।

`collectGreetings` process ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

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

`count_greetings` variable runtime पर compute किया जाएगा।

#### 4.1.2. Report फ़ाइल emit करें और outputs को नाम दें

सिद्धांत रूप में हमें बस report फ़ाइल को `output:` ब्लॉक में जोड़ना होगा।

हालाँकि, जब हम इस पर हैं, हम अपने output declarations में कुछ `emit:` tags भी जोड़ने जा रहे हैं। ये हमें positional indices का उपयोग करने के बजाय नाम से outputs को select करने में सक्षम बनाएंगे।

Process ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

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

`emit:` tags वैकल्पिक हैं, और हम केवल एक output में tag जोड़ सकते थे।
लेकिन जैसा कि कहावत है, दोनों क्यों नहीं?

!!! tip "सुझाव"

    अगर तुम `emit:` का उपयोग करके एक process के outputs को नाम नहीं देते हो, तो तुम अभी भी उन्हें उनके respective (zero-based) index का उपयोग करके व्यक्तिगत रूप से access कर सकते हो।
    उदाहरण के लिए, तुम पहले output को प्राप्त करने के लिए `<process>.out[0]` का उपयोग करोगे, दूसरे output को प्राप्त करने के लिए `<process>.out[1]`, और इसी तरह।

    हम outputs को नाम देना पसंद करते हैं क्योंकि अन्यथा, गलती से गलत index grab करना बहुत आसान है, खासकर जब process बहुत सारे outputs उत्पन्न करता है।

### 4.2. Workflow outputs को अपडेट करें

अब जब हमारे पास `collectGreetings` process से दो outputs आ रहे हैं, तो `collectGreetings.out` output में दो channels हैं:

- `collectGreetings.out.outfile` में अंतिम output फ़ाइल है
- `collectGreetings.out.report` में report फ़ाइल है

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-report.svg"
</figure>

हमें तदनुसार workflow outputs को अपडेट करना होगा।

#### 4.2.1. `publish:` सेक्शन को अपडेट करें

`workflow ब्लॉक` में, निम्नलिखित कोड परिवर्तन करें:

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

जैसा कि तुम देख सकते हो, विशिष्ट process outputs को refer करना अब तुच्छ है।
जब हम भाग 5 (Containers) में अपनी pipeline में एक और स्टेप जोड़ने जाते हैं, तो हम आसानी से `collectGreetings.out.outfile` को refer कर सकेंगे और इसे नए process को सौंप सकेंगे (स्पॉइलर: नया process `cowpy` कहलाता है)।

लेकिन अभी के लिए, आइए workflow-level outputs को अपडेट करना समाप्त करें।

#### 4.2.2. `output` ब्लॉक को अपडेट करें

`output` ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

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

हमें `collected` output परिभाषा को अपडेट करने की आवश्यकता नहीं है क्योंकि वह नाम नहीं बदला है।
हमें बस नया output जोड़ना होगा।

### 4.3. Workflow चलाएं

आइए इसे greetings के वर्तमान batch के साथ चलाने की कोशिश करें।

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [ecstatic_wilson] DSL2 - revision: c80285f8c8

    executor >  local (1)
    [c5/4c6ca9] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [0e/6cbc59] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [02/61ead2] collectGreetings   [100%] 1 of 1 ✔
    ```

अगर तुम `results/hello_workflow/` डायरेक्टरी में देखते हो, तो तुम्हें नई report फ़ाइल मिलेगी, `trio-report.txt`।
यह सत्यापित करने के लिए इसे खोलो कि workflow ने सही तरीके से greetings की गिनती की रिपोर्ट की जो process की गई थी।

??? abstract "फ़ाइल सामग्री"

    ```txt title="trio-report.txt"
    There were 3 greetings in this batch.
    ```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-4-way.svg"
</figure>

CSV में और greetings जोड़ने और यह टेस्ट करने के लिए स्वतंत्र महसूस करो कि क्या होता है।

### सारांश

तुम जानते हो कि एक process को कई named outputs कैसे emit करवाया जाए और workflow स्तर पर उन्हें उचित रूप से कैसे संभाला जाए।

अधिक सामान्य रूप से, तुम processes को सामान्य तरीकों से एक साथ जोड़ने में शामिल प्रमुख सिद्धांतों को समझते हो।

### आगे क्या है?

एक अतिरिक्त लंबा ब्रेक लो, तुमने इसे अर्जित किया है।

जब तुम तैयार हो, तो बेहतर maintainability और कोड efficiency के लिए अपने कोड को modularize करना सीखने के लिए [**भाग 4: Hello Modules**](./04_hello_modules.md) पर जाओ।

---

## क्विज़

<quiz>
Workflow ब्लॉक में एक process के output को कैसे access करते हैं?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

और जानें: [1.4. पहले process के output को दूसरे process में पास करें](#14-pass-the-output-of-the-first-process-to-the-second-process)
</quiz>

<quiz>
Nextflow में process execution का क्रम क्या निर्धारित करता है?
- [ ] Workflow ब्लॉक में processes लिखे जाने का क्रम
- [ ] Process नाम से alphabetical क्रम
- [x] Processes के बीच डेटा निर्भरताएं
- [ ] Parallel execution के लिए random क्रम

और जानें: [1.4. पहले process के output को दूसरे process में पास करें](#14-pass-the-output-of-the-first-process-to-the-second-process)
</quiz>

<quiz>
Downstream process के लिए सभी outputs को एक single list में इकट्ठा करने के लिए `???` को किस operator से बदलना चाहिए?

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

और जानें: [2.4. Greetings को एक single input में एकत्र करने के लिए एक operator का उपयोग करें](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
`collect()` operator का उपयोग कब करना चाहिए?
- [ ] जब तुम items को parallel में process करना चाहते हो
- [ ] जब तुम्हें channel सामग्री को filter करने की आवश्यकता हो
- [x] जब एक downstream process को एक upstream process से सभी items की आवश्यकता हो
- [ ] जब तुम कई processes में डेटा को split करना चाहते हो

और जानें: [2.4. Greetings को एक single input में एकत्र करने के लिए एक operator का उपयोग करें](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
एक process से एक named output को कैसे access करते हैं?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

और जानें: [4.1.2. Report फ़ाइल emit करें और outputs को नाम दें](#412-emit-the-report-file-and-name-outputs)
</quiz>

<quiz>
एक process में एक output को नाम देने के लिए सही सिंटैक्स क्या है?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

और जानें: [4.1.2. Report फ़ाइल emit करें और outputs को नाम दें](#412-emit-the-report-file-and-name-outputs)
</quiz>

<quiz>
एक process को कई inputs प्रदान करते समय, क्या सत्य होना चाहिए?
- [ ] सभी inputs एक ही प्रकार के होने चाहिए
- [ ] Inputs को alphabetical क्रम में प्रदान किया जाना चाहिए
- [x] Inputs का क्रम input ब्लॉक में परिभाषित क्रम से मेल खाना चाहिए
- [ ] एक समय में केवल दो inputs प्रदान किए जा सकते हैं

और जानें: [3. एक process को एक से अधिक input पास करें](#3-pass-more-than-one-input-to-a-process)
</quiz>
