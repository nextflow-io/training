# वर्कफ़्लो के वर्कफ़्लो

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

जब आप एक pipeline विकसित कर रहे होते हैं, तो आप अक्सर खुद को विभिन्न डेटा प्रकारों या विश्लेषण चरणों के लिए समान प्रक्रियाओं के अनुक्रम बनाते हुए पाते हैं। आप इन process अनुक्रमों को कॉपी और पेस्ट करते हुए समाप्त हो सकते हैं, जिससे डुप्लिकेट कोड बनता है जिसे बनाए रखना कठिन होता है; या आप एक विशाल workflow बना सकते हैं जो समझना और संशोधित करना मुश्किल है।

Nextflow की सबसे शक्तिशाली विशेषताओं में से एक छोटे, पुन: प्रयोज्य workflow मॉड्यूल से जटिल pipelines बनाने की इसकी क्षमता है। यह modular दृष्टिकोण pipelines को विकसित करना, परीक्षण करना और बनाए रखना आसान बनाता है।

### सीखने के लक्ष्य

इस side quest में, हम ऐसे workflow मॉड्यूल विकसित करने का पता लगाएंगे जिन्हें अलग से परीक्षण और उपयोग किया जा सकता है, उन मॉड्यूल को एक बड़े pipeline में संयोजित करेंगे, और मॉड्यूल के बीच डेटा प्रवाह का प्रबंधन करेंगे।

इस side quest के अंत तक, आप निम्न करने में सक्षम होंगे:

- जटिल pipelines को तार्किक, पुन: प्रयोज्य इकाइयों में विभाजित करना
- प्रत्येक workflow मॉड्यूल का स्वतंत्र रूप से परीक्षण करना
- नए pipelines बनाने के लिए workflows को मिलाना और मिलान करना
- विभिन्न pipelines में सामान्य workflow मॉड्यूल साझा करना
- अपने कोड को अधिक maintainable और समझने में आसान बनाना

ये कौशल आपको स्वच्छ, maintainable कोड संरचना बनाए रखते हुए जटिल pipelines बनाने में मदद करेंगे।

### पूर्वापेक्षाएँ

इस side quest को शुरू करने से पहले आपको:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष शुरुआती पाठ्यक्रम पूरा किया होना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (processes, channels, operators, modules) का उपयोग करने में सहज होना चाहिए

---

## 0. शुरुआत करें

#### प्रशिक्षण codespace खोलें

यदि आपने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित अनुसार प्रशिक्षण वातावरण खोलना सुनिश्चित करें।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाएं

चलिए उस डायरेक्टरी में चलते हैं जहां इस ट्यूटोरियल के लिए फ़ाइलें स्थित हैं।

```bash
cd side-quests/workflows_of_workflows
```

आप VSCode को इस डायरेक्टरी पर ध्यान केंद्रित करने के लिए सेट कर सकते हैं:

```bash
code .
```

#### सामग्री की समीक्षा करें

आपको एक `modules` डायरेक्टरी मिलेगी जिसमें कई process परिभाषाएं हैं जो 'Hello Nextflow' में आपने जो सीखा था उस पर आधारित हैं:

```console title="डायरेक्टरी सामग्री"
modules/
├── say_hello.nf             # एक अभिवादन बनाता है (Hello Nextflow से)
├── say_hello_upper.nf       # uppercase में परिवर्तित करता है (Hello Nextflow से)
├── timestamp_greeting.nf    # अभिवादन में timestamps जोड़ता है
├── validate_name.nf         # इनपुट नामों को validate करता है
└── reverse_text.nf          # टेक्स्ट सामग्री को उलट देता है
```

#### असाइनमेंट की समीक्षा करें

आपकी चुनौती इन मॉड्यूल को दो अलग workflows में असेंबल करना है जिन्हें हम फिर एक मुख्य workflow में संयोजित करेंगे:

- एक `GREETING_WORKFLOW` जो नामों को validate करता है, अभिवादन बनाता है, और timestamps जोड़ता है
- एक `TRANSFORM_WORKFLOW` जो टेक्स्ट को uppercase में परिवर्तित करता है और इसे उलट देता है

#### तैयारी चेकलिस्ट

लगता है कि आप शुरू करने के लिए तैयार हैं?

- [ ] मैं इस पाठ्यक्रम के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता/समझती हूं
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी working डायरेक्टरी उचित रूप से सेट की है
- [ ] मैं असाइनमेंट को समझता/समझती हूं

यदि आप सभी बॉक्स को चेक कर सकते हैं, तो आप जाने के लिए तैयार हैं।

---

## 1. Greeting Workflow बनाएं

चलिए एक workflow बनाकर शुरू करते हैं जो नामों को validate करता है और timestamped अभिवादन उत्पन्न करता है।

### 1.1. workflow संरचना बनाएं

```bash title="workflow डायरेक्टरी और फ़ाइल बनाएं"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. पहले (sub)workflow कोड जोड़ें

यह कोड `workflows/greeting.nf` में जोड़ें:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Processes की chain: validate -> greeting बनाएं -> timestamp जोड़ें
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

यह एक पूर्ण workflow है, जिसकी संरचना 'Hello Nextflow' ट्यूटोरियल में आपने देखे गए workflows के समान है, जिसे हम स्वतंत्र रूप से परीक्षण कर सकते हैं। चलिए अब इसे आज़माते हैं:

```bash
nextflow run workflows/greeting.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
    executor >  local (9)
    [51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
    [2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
    [8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
    ```

यह अपेक्षा के अनुसार काम करता है, लेकिन इसे composable बनाने के लिए हमें कुछ चीजें बदलने की आवश्यकता है।

### 1.3. workflow को composable बनाएं

Composable workflows में 'Hello Nextflow' ट्यूटोरियल में आपने देखे गए workflows से कुछ अंतर हैं:

- workflow block को नाम दिए जाने की आवश्यकता है
- इनपुट `take:` keyword का उपयोग करके घोषित किए जाते हैं
- Workflow सामग्री `main:` block के अंदर रखी जाती है
- आउटपुट `emit:` keyword का उपयोग करके घोषित किए जाते हैं

चलिए greeting workflow को इस संरचना से मेल खाने के लिए अपडेट करते हैं। कोड को निम्नलिखित में बदलें:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // इनपुट channel नामों के साथ

    main:
        // Processes की chain: validate -> greeting बनाएं -> timestamp जोड़ें
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // मूल अभिवादन
        timestamped = timestamped_ch  // Timestamped अभिवादन
}
```

आप देख सकते हैं कि workflow को अब नाम दिया गया है और इसमें `take:` और `emit:` block है, और ये वे कनेक्शन हैं जिनका उपयोग हम उच्च स्तर के workflow को संयोजित करने के लिए करेंगे।
Workflow सामग्री भी `main:` block के अंदर रखी गई है। यह भी नोट करें कि हमने `names_ch` इनपुट channel घोषणा को हटा दिया है, क्योंकि इसे अब workflow को एक argument के रूप में पास किया जाता है।

चलिए यह देखने के लिए workflow को फिर से परीक्षण करते हैं कि क्या यह अपेक्षा के अनुसार काम करता है:

```bash
nextflow run workflows/greeting.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

यह आपको एक अन्य नई अवधारणा के बारे में बताता है, एक 'entry workflow'। Entry workflow वह workflow है जो तब कॉल किया जाता है जब आप एक Nextflow script चलाते हैं। डिफ़ॉल्ट रूप से, Nextflow एक unnamed workflow को entry workflow के रूप में उपयोग करेगा, जब मौजूद हो, और यह वही है जो आप अब तक कर रहे थे, workflow blocks इस तरह शुरू करते हुए:

```groovy title="hello.nf" linenums="1"
workflow {
```

लेकिन हमारे greeting workflow में un-named workflow नहीं है, बल्कि हमारे पास एक named workflow है:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

इसलिए Nextflow ने एक error फेंकी और वह नहीं किया जो हम चाहते थे।

हमने `take:`/`emit:` सिंटैक्स इसलिए नहीं जोड़ा कि हम workflow को सीधे कॉल कर सकें - हमने इसे इसलिए किया कि हम इसे अन्य workflows के साथ संयोजित कर सकें। समाधान एक मुख्य script बनाना है जिसमें एक unnamed entry workflow हो जो हमारे named workflow को import और कॉल करे।

### 1.4. मुख्य workflow बनाएं और परीक्षण करें

अब हम एक मुख्य workflow बनाएंगे जो `greeting` workflow को import और उपयोग करता है।

`main.nf` बनाएं:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

नोट करें कि इस फ़ाइल में हमारी workflow entry un-named है, और ऐसा इसलिए है क्योंकि हम इसे entry workflow के रूप में उपयोग करने जा रहे हैं।

इसे चलाएं और आउटपुट देखें:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
    executor >  local (9)
    [05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
    [b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
    [ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
    ```

यह काम करता है! हमने named greeting workflow को एक मुख्य workflow में लपेटा है जिसमें un-named entry `workflow` block है। मुख्य workflow `GREETING_WORKFLOW` workflow का उपयोग लगभग (बिल्कुल नहीं) एक process की तरह कर रहा है, और `names` channel को एक argument के रूप में पास कर रहा है।

### मुख्य बात

इस खंड में, आपने कई महत्वपूर्ण अवधारणाएं सीखी हैं:

- **Named Workflows**: एक named workflow (`GREETING_WORKFLOW`) बनाना जिसे import और पुन: उपयोग किया जा सकता है
- **Workflow Interfaces**: एक composable workflow बनाने के लिए `take:` के साथ स्पष्ट इनपुट और `emit:` के साथ आउटपुट परिभाषित करना
- **Entry Points**: यह समझना कि Nextflow को एक script चलाने के लिए एक unnamed entry workflow की आवश्यकता होती है
- **Workflow Composition**: एक अन्य workflow के भीतर named workflow को import और उपयोग करना
- **Workflow Namespaces**: `.out` namespace का उपयोग करके workflow आउटपुट तक पहुंचना (`GREETING_WORKFLOW.out.greetings`)

अब आपके पास एक काम करने वाला greeting workflow है जो:

- इनपुट के रूप में नामों के channel लेता है
- प्रत्येक नाम को validate करता है
- प्रत्येक valid नाम के लिए एक अभिवादन बनाता है
- अभिवादन में timestamps जोड़ता है
- मूल और timestamped दोनों अभिवादन को आउटपुट के रूप में expose करता है

यह modular दृष्टिकोण आपको greeting workflow को स्वतंत्र रूप से परीक्षण करने या बड़े pipelines में एक component के रूप में उपयोग करने की अनुमति देता है।

---

## 2. Transform Workflow जोड़ें

अब चलिए एक workflow बनाते हैं जो अभिवादन पर text transformations लागू करता है।

### 2.1. workflow फ़ाइल बनाएं

```bash
touch workflows/transform.nf
```

### 2.2. workflow कोड जोड़ें

यह कोड `workflows/transform.nf` में जोड़ें:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // संदेशों के साथ इनपुट channel

    main:
        // क्रम में transformations लागू करें
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Uppercase अभिवादन
        reversed = reversed_ch  // Reversed uppercase अभिवादन
}
```

हम यहां composable सिंटैक्स की व्याख्या नहीं दोहराएंगे, लेकिन ध्यान दें कि named workflow फिर से `take:` और `emit:` block के साथ घोषित किया गया है, और workflow सामग्री `main:` block के अंदर रखी गई है।

### 2.3. मुख्य workflow को अपडेट करें

दोनों workflows का उपयोग करने के लिए `main.nf` को अपडेट करें:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // greeting workflow चलाएं
    GREETING_WORKFLOW(names)

    // transform workflow चलाएं
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // परिणाम देखें
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

पूर्ण pipeline चलाएं:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (13)
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/a0/d4f5df4d6344604498fa47a6084a11/UPPER-timestamped_Bob-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/69/b5e37f6c79c2fd38adb75d0eca8f87/UPPER-timestamped_Charlie-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/cd/e6a7e0b17e7d5a2f71bb8123cd53a7/UPPER-timestamped_Alice-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/7a/7a222f7957b35d1d121338566a24ac/REVERSED-UPPER-timestamped_Bob-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/46/8d19af62e33a5a6417c773496e0f90/REVERSED-UPPER-timestamped_Charlie-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

यदि आप उन reversed फ़ाइलों में से एक पर नज़र डालते हैं, तो आप देखेंगे कि यह अभिवादन का uppercase संस्करण उलटा है:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed फ़ाइल सामग्री"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### मुख्य बात

अब आपके पास एक पूर्ण pipeline होना चाहिए जो:

- greeting workflow के माध्यम से नामों को प्रोसेस करता है
- timestamped अभिवादन को transform workflow में भेजता है
- अभिवादन के uppercase और reversed दोनों संस्करण उत्पन्न करता है

---

## सारांश

इस side quest में, हमने Nextflow में workflow composition की शक्तिशाली अवधारणा का पता लगाया है, जो हमें छोटे, पुन: प्रयोज्य components से जटिल pipelines बनाने की अनुमति देती है।

यह modular दृष्टिकोण monolithic pipelines की तुलना में कई लाभ प्रदान करता है:

- प्रत्येक workflow को स्वतंत्र रूप से विकसित, परीक्षण और debug किया जा सकता है
- Workflows को विभिन्न pipelines में पुन: उपयोग किया जा सकता है
- समग्र pipeline संरचना अधिक पठनीय और maintainable बन जाती है
- एक workflow में परिवर्तन अन्य को जरूरी नहीं प्रभावित करते हैं यदि interfaces सुसंगत रहें
- Entry points को आपके pipeline के विभिन्न भागों को आवश्यकतानुसार चलाने के लिए कॉन्फ़िगर किया जा सकता है

_यह ध्यान रखना महत्वपूर्ण है कि workflows को कॉल करना थोड़ा processes को कॉल करने जैसा है, लेकिन यह वास्तव में एक ही चीज नहीं है। उदाहरण के लिए, आप N आकार के channel के साथ इसे कॉल करके एक workflow को N बार नहीं चला सकते - आपको workflow को N आकार के channel को पास करने और आंतरिक रूप से iterate करने की आवश्यकता होगी।_

अपने स्वयं के काम में इन तकनीकों को लागू करने से आप अधिक परिष्कृत Nextflow pipelines बनाने में सक्षम होंगे जो maintainable और scalable रहते हुए जटिल bioinformatics कार्यों को संभाल सकते हैं।

### मुख्य पैटर्न

1.  **Workflow संरचना**: हमने `take:` और `emit:` सिंटैक्स का उपयोग करके प्रत्येक workflow के लिए स्पष्ट इनपुट और आउटपुट परिभाषित किए, components के बीच अच्छी तरह से परिभाषित interfaces बनाए, और `main:` block के भीतर workflow तर्क को लपेटा।

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // इनपुट channels यहां घोषित किए जाते हैं
            input_ch

        main:
            // Workflow logic यहां जाता है
            // यह वह जगह है जहां processes को कॉल किया जाता है और channels को manipulate किया जाता है
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // आउटपुट channels यहां घोषित किए जाते हैं
            output_ch = result_ch
    }
    ```

2.  **Workflow imports:** हमने दो स्वतंत्र workflow मॉड्यूल बनाए और उन्हें include statements के साथ एक मुख्य pipeline में import किया।

    - एकल workflow include करें

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - एकाधिक workflows include करें

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - नाम संघर्षों से बचने के लिए alias के साथ include करें

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Entry points**: Nextflow को यह जानने के लिए एक unnamed entry workflow की आवश्यकता होती है कि निष्पादन कहां से शुरू करना है। यह entry workflow आपके named workflows को कॉल करता है।

    - Unnamed workflow (entry point)

    ```groovy
    workflow {
        // यह entry point है जब script को चलाया जाता है
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Named workflow (entry workflow से कॉल किया जाता है)

    ```groovy
    workflow NAMED_WORKFLOW {
        // entry workflow से कॉल किया जाना चाहिए
    }
    ```

4.  **डेटा प्रवाह का प्रबंधन:** हमने सीखा कि namespace notation (`WORKFLOW_NAME.out.channel_name`) का उपयोग करके workflow आउटपुट तक कैसे पहुंचें और उन्हें अन्य workflows को पास करें।

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### अतिरिक्त संसाधन

- [Nextflow Workflow Documentation](https://www.nextflow.io/docs/latest/workflow.html)
- [Channel Operators Reference](https://www.nextflow.io/docs/latest/operator.html)
- [Error Strategy Documentation](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## आगे क्या है?

[Side Quests के मेनू](./index.md) पर वापस जाएं या सूची में अगले विषय पर जाने के लिए पृष्ठ के निचले दाहिने कोने में बटन पर क्लिक करें।
