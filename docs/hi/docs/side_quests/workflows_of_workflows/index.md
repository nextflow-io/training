# वर्कफ़्लो के वर्कफ़्लो

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

जब तुम एक पाइपलाइन विकसित कर रहे होते हो, तो अक्सर तुम खुद को अलग-अलग डेटा प्रकारों या विश्लेषण चरणों के लिए प्रोसेस के समान अनुक्रम बनाते हुए पाते हो। तुम इन प्रोसेस अनुक्रमों को कॉपी-पेस्ट करते रह सकते हो, जिससे डुप्लीकेट कोड बन जाता है जिसे बनाए रखना मुश्किल होता है; या तुम एक विशाल वर्कफ़्लो बना सकते हो जिसे समझना और बदलना कठिन हो।

Nextflow की सबसे शक्तिशाली विशेषताओं में से एक है छोटे, पुन: उपयोग योग्य वर्कफ़्लो मॉड्यूल से जटिल पाइपलाइन बनाने की क्षमता। यह मॉड्यूलर दृष्टिकोण पाइपलाइन को विकसित करना, परीक्षण करना और बनाए रखना आसान बनाता है।

### सीखने के लक्ष्य

इस साइड क्वेस्ट में, हम यह जानेंगे कि वर्कफ़्लो मॉड्यूल कैसे विकसित करें जिन्हें अलग से परीक्षण और उपयोग किया जा सके, उन मॉड्यूल को एक बड़ी पाइपलाइन में कैसे जोड़ें, और मॉड्यूल के बीच डेटा प्रवाह कैसे प्रबंधित करें।

इस साइड क्वेस्ट के अंत तक, तुम यह कर पाओगे:

- जटिल पाइपलाइन को तार्किक, पुन: उपयोग योग्य इकाइयों में विभाजित करना
- प्रत्येक वर्कफ़्लो मॉड्यूल को स्वतंत्र रूप से परीक्षण करना
- नई पाइपलाइन बनाने के लिए वर्कफ़्लो को मिलाना और मिलान करना
- विभिन्न पाइपलाइन में सामान्य वर्कफ़्लो मॉड्यूल साझा करना
- अपने कोड को अधिक रखरखाव योग्य और समझने में आसान बनाना

ये कौशल तुम्हें स्वच्छ, रखरखाव योग्य कोड संरचना बनाए रखते हुए जटिल पाइपलाइन बनाने में मदद करेंगे।

### पूर्वापेक्षाएँ

इस साइड क्वेस्ट को शुरू करने से पहले तुम्हें:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष शुरुआती कोर्स पूरा करना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (प्रोसेस, चैनल, ऑपरेटर, मॉड्यूल) का उपयोग करने में सहज होना चाहिए।

---

## 0. शुरू करना

#### ट्रेनिंग codespace खोलो

अगर तुमने अभी तक नहीं किया है, तो [Environment Setup](../envsetup/index.md) में बताए अनुसार ट्रेनिंग वातावरण खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाओ

चलो उस डायरेक्टरी में जाते हैं जहाँ इस ट्यूटोरियल की फ़ाइलें हैं।

```bash
cd side-quests/workflows_of_workflows
```

तुम VSCode को इस डायरेक्टरी पर फ़ोकस करने के लिए सेट कर सकते हो:

```bash
code .
```

#### सामग्री की समीक्षा करो

तुम्हें एक `modules` डायरेक्टरी मिलेगी जिसमें कई प्रोसेस परिभाषाएँ हैं जो 'Hello Nextflow' में सीखी गई बातों पर आधारित हैं:

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

#### असाइनमेंट की समीक्षा करो

तुम्हारी चुनौती इन मॉड्यूल को दो अलग वर्कफ़्लो में जोड़ना है जिन्हें हम फिर एक मुख्य वर्कफ़्लो में जोड़ेंगे:

- एक `GREETING_WORKFLOW` जो नामों को सत्यापित करता है, अभिवादन बनाता है, और टाइमस्टैम्प जोड़ता है
- एक `TRANSFORM_WORKFLOW` जो टेक्स्ट को अपरकेस में बदलता है और उसे उलटता है

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### तैयारी की जाँच सूची

क्या तुम शुरू करने के लिए तैयार हो?

- [ ] मैं इस कोर्स का लक्ष्य और इसकी पूर्वापेक्षाएँ समझता/समझती हूँ
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी वर्किंग डायरेक्टरी उचित रूप से सेट की है
- [ ] मैं असाइनमेंट समझता/समझती हूँ

अगर तुम सभी बॉक्स चेक कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. Greeting Workflow बनाओ

चलो एक ऐसा वर्कफ़्लो बनाकर शुरू करते हैं जो नामों को सत्यापित करता है और टाइमस्टैम्प वाले अभिवादन उत्पन्न करता है।

### 1.1. वर्कफ़्लो संरचना बनाओ

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. पहला (सब)वर्कफ़्लो कोड जोड़ो

`workflows/greeting.nf` में यह कोड जोड़ो:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // प्रोसेस को जंजीर में जोड़ो: सत्यापित करो -> अभिवादन बनाओ -> टाइमस्टैम्प जोड़ो
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

यह एक पूर्ण वर्कफ़्लो है, जिसकी संरचना 'Hello Nextflow' ट्यूटोरियल में देखी गई संरचनाओं के समान है, जिसे हम स्वतंत्र रूप से परीक्षण कर सकते हैं। चलो अभी वह करते हैं:

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

यह अपेक्षा के अनुसार काम करता है, लेकिन इसे composable बनाने के लिए हमें कुछ चीज़ें बदलनी होंगी।

### 1.3. वर्कफ़्लो को composable बनाओ

Composable वर्कफ़्लो में 'Hello Nextflow' ट्यूटोरियल में देखे गए वर्कफ़्लो से कुछ अंतर होते हैं:

- वर्कफ़्लो ब्लॉक का नाम होना चाहिए
- इनपुट `take:` कीवर्ड का उपयोग करके घोषित किए जाते हैं
- वर्कफ़्लो की सामग्री `main:` ब्लॉक के अंदर रखी जाती है
- आउटपुट `emit:` कीवर्ड का उपयोग करके घोषित किए जाते हैं

चलो greeting वर्कफ़्लो को इस संरचना से मेल खाने के लिए अपडेट करते हैं। कोड को निम्नलिखित में बदलो:

<!-- TODO: switch to before/after tabs -->

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // नामों के साथ इनपुट चैनल

    main:
        // प्रोसेस को जंजीर में जोड़ो: सत्यापित करो -> अभिवादन बनाओ -> टाइमस्टैम्प जोड़ो
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // मूल अभिवादन
        timestamped = timestamped_ch  // टाइमस्टैम्प वाले अभिवादन
}
```

तुम देख सकते हो कि वर्कफ़्लो का अब नाम है और इसमें `take:` और `emit:` ब्लॉक हैं, और ये वे कनेक्शन हैं जिनका उपयोग हम एक उच्च स्तरीय वर्कफ़्लो बनाने के लिए करेंगे।
वर्कफ़्लो की सामग्री भी `main:` ब्लॉक के अंदर रखी गई है। ध्यान दो कि हमने `names_ch` इनपुट चैनल घोषणा हटा दी है, क्योंकि अब इसे वर्कफ़्लो में एक आर्गुमेंट के रूप में पास किया जाता है।

चलो वर्कफ़्लो को फिर से परीक्षण करते हैं यह देखने के लिए कि यह अपेक्षा के अनुसार काम करता है या नहीं:

```bash
nextflow run workflows/greeting.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

यह तुम्हें एक नई अवधारणा के बारे में बताता है, एक 'entry workflow'। Entry workflow वह वर्कफ़्लो है जिसे Nextflow स्क्रिप्ट चलाने पर कॉल किया जाता है। डिफ़ॉल्ट रूप से, Nextflow एक बिना नाम वाले वर्कफ़्लो को entry workflow के रूप में उपयोग करेगा, जब वह मौजूद हो, और यही तुम अब तक करते आए हो, वर्कफ़्लो ब्लॉक इस तरह शुरू होते हैं:

```groovy title="hello.nf" linenums="1"
workflow {
```

लेकिन हमारे greeting वर्कफ़्लो में कोई बिना नाम वाला वर्कफ़्लो नहीं है, बल्कि हमारे पास एक नामित वर्कफ़्लो है:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

इसीलिए Nextflow ने एक त्रुटि दी और वह नहीं किया जो हम चाहते थे।

हमने `take:`/`emit:` सिंटैक्स इसलिए नहीं जोड़ा ताकि हम वर्कफ़्लो को सीधे कॉल कर सकें - हमने इसे इसलिए किया ताकि हम इसे अन्य वर्कफ़्लो के साथ जोड़ सकें। समाधान यह है कि एक मुख्य स्क्रिप्ट बनाई जाए जिसमें एक बिना नाम वाला entry workflow हो जो हमारे नामित वर्कफ़्लो को import और कॉल करे।

### 1.4. मुख्य वर्कफ़्लो बनाओ और परीक्षण करो

अब हम एक मुख्य वर्कफ़्लो बनाएंगे जो `greeting` वर्कफ़्लो को import और उपयोग करता है।

`main.nf` बनाओ:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

ध्यान दो कि इस फ़ाइल में हमारा वर्कफ़्लो entry बिना नाम का है, और ऐसा इसलिए है क्योंकि हम इसे entry workflow के रूप में उपयोग करने वाले हैं।

इसे चलाओ और आउटपुट देखो:

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

यह काम करता है! हमने नामित greeting वर्कफ़्लो को एक बिना नाम वाले entry `workflow` ब्लॉक वाले मुख्य वर्कफ़्लो में लपेटा है। मुख्य वर्कफ़्लो `GREETING_WORKFLOW` वर्कफ़्लो का उपयोग लगभग (पूरी तरह नहीं) एक प्रोसेस की तरह कर रहा है, और `names` चैनल को एक आर्गुमेंट के रूप में पास कर रहा है।

### सारांश

इस अनुभाग में, तुमने कई महत्वपूर्ण अवधारणाएँ सीखी हैं:

- **Named Workflows**: एक नामित वर्कफ़्लो (`GREETING_WORKFLOW`) बनाना जिसे import और पुन: उपयोग किया जा सके
- **Workflow Interfaces**: एक composable वर्कफ़्लो बनाने के लिए `take:` के साथ स्पष्ट इनपुट और `emit:` के साथ आउटपुट परिभाषित करना
- **Entry Points**: यह समझना कि Nextflow को स्क्रिप्ट चलाने के लिए एक बिना नाम वाले entry workflow की आवश्यकता है
- **Workflow Composition**: किसी अन्य वर्कफ़्लो के भीतर एक नामित वर्कफ़्लो को import और उपयोग करना
- **Workflow Namespaces**: `.out` namespace (`GREETING_WORKFLOW.out.greetings`) का उपयोग करके वर्कफ़्लो आउटपुट तक पहुँचना

अब तुम्हारे पास एक काम करने वाला greeting वर्कफ़्लो है जो:

- इनपुट के रूप में नामों का एक चैनल लेता है
- प्रत्येक नाम को सत्यापित करता है
- प्रत्येक वैध नाम के लिए एक अभिवादन बनाता है
- अभिवादन में टाइमस्टैम्प जोड़ता है
- मूल और टाइमस्टैम्प वाले दोनों अभिवादन को आउटपुट के रूप में उजागर करता है

यह मॉड्यूलर दृष्टिकोण तुम्हें greeting वर्कफ़्लो को स्वतंत्र रूप से परीक्षण करने या बड़ी पाइपलाइन में एक घटक के रूप में उपयोग करने की अनुमति देता है।

---

## 2. Transform Workflow जोड़ो

अब चलो एक ऐसा वर्कफ़्लो बनाते हैं जो अभिवादन पर टेक्स्ट परिवर्तन लागू करता है।

### 2.1. वर्कफ़्लो फ़ाइल बनाओ

```bash
touch workflows/transform.nf
```

### 2.2. वर्कफ़्लो कोड जोड़ो

`workflows/transform.nf` में यह कोड जोड़ो:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // संदेशों के साथ इनपुट चैनल

    main:
        // क्रम में परिवर्तन लागू करो
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // अपरकेस अभिवादन
        reversed = reversed_ch  // उलटे अपरकेस अभिवादन
}
```

हम यहाँ composable सिंटैक्स की व्याख्या नहीं दोहराएंगे, लेकिन ध्यान दो कि नामित वर्कफ़्लो फिर से `take:` और `emit:` ब्लॉक के साथ घोषित किया गया है, और वर्कफ़्लो की सामग्री `main:` ब्लॉक के अंदर रखी गई है।

### 2.3. मुख्य वर्कफ़्लो अपडेट करो

दोनों वर्कफ़्लो का उपयोग करने के लिए `main.nf` अपडेट करो:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // greeting वर्कफ़्लो चलाओ
    GREETING_WORKFLOW(names)

    // transform वर्कफ़्लो चलाओ
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // परिणाम देखो
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

पूरी पाइपलाइन चलाओ:

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

अगर तुम उन उलटी फ़ाइलों में से एक देखो, तो तुम पाओगे कि यह अभिवादन का अपरकेस संस्करण उलटा है:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### सारांश

अब तुम्हारे पास एक पूर्ण पाइपलाइन होनी चाहिए जो:

- greeting वर्कफ़्लो के माध्यम से नामों को प्रोसेस करती है
- टाइमस्टैम्प वाले अभिवादन को transform वर्कफ़्लो में भेजती है
- अभिवादन के अपरकेस और उलटे दोनों संस्करण उत्पन्न करती है

---

## सारांश

इस साइड क्वेस्ट में, हमने Nextflow में वर्कफ़्लो composition की शक्तिशाली अवधारणा का पता लगाया है, जो हमें छोटे, पुन: उपयोग योग्य घटकों से जटिल पाइपलाइन बनाने की अनुमति देती है।

यह मॉड्यूलर दृष्टिकोण monolithic पाइपलाइन की तुलना में कई फायदे प्रदान करता है:

- प्रत्येक वर्कफ़्लो को स्वतंत्र रूप से विकसित, परीक्षण और डीबग किया जा सकता है
- वर्कफ़्लो को विभिन्न पाइपलाइन में पुन: उपयोग किया जा सकता है
- समग्र पाइपलाइन संरचना अधिक पठनीय और रखरखाव योग्य बन जाती है
- एक वर्कफ़्लो में बदलाव जरूरी नहीं कि दूसरों को प्रभावित करें यदि interfaces सुसंगत रहें
- Entry points को तुम्हारी पाइपलाइन के विभिन्न हिस्सों को चलाने के लिए कॉन्फ़िगर किया जा सकता है

_यह ध्यान रखना महत्वपूर्ण है कि वर्कफ़्लो को कॉल करना प्रोसेस को कॉल करने जैसा थोड़ा लगता है, लेकिन यह वास्तव में एक जैसा नहीं है। उदाहरण के लिए, तुम N आकार के चैनल के साथ कॉल करके किसी वर्कफ़्लो को N बार नहीं चला सकते - तुम्हें N आकार का चैनल वर्कफ़्लो को पास करना होगा और आंतरिक रूप से iterate करना होगा।_

अपने काम में इन तकनीकों को लागू करने से तुम अधिक परिष्कृत Nextflow पाइपलाइन बना पाओगे जो जटिल bioinformatics कार्यों को संभाल सकती हैं और साथ ही रखरखाव योग्य और scalable भी रहती हैं।

### मुख्य पैटर्न

1.  **Workflow structure**: हमने `take:` और `emit:` सिंटैक्स का उपयोग करके प्रत्येक वर्कफ़्लो के लिए स्पष्ट इनपुट और आउटपुट परिभाषित किए, घटकों के बीच अच्छी तरह से परिभाषित interfaces बनाए, और वर्कफ़्लो लॉजिक को `main:` ब्लॉक के अंदर लपेटा।

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // इनपुट चैनल यहाँ घोषित किए जाते हैं
            input_ch

        main:
            // वर्कफ़्लो लॉजिक यहाँ जाता है
            // यहाँ प्रोसेस कॉल किए जाते हैं और चैनल manipulate किए जाते हैं
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // आउटपुट चैनल यहाँ घोषित किए जाते हैं
            output_ch = result_ch
    }
    ```

2.  **Workflow imports:** हमने दो स्वतंत्र वर्कफ़्लो मॉड्यूल बनाए और उन्हें include statements के साथ एक मुख्य पाइपलाइन में import किया।

    - एक वर्कफ़्लो include करो

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - कई वर्कफ़्लो include करो

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - नाम संघर्ष से बचने के लिए alias के साथ include करो

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Entry points**: Nextflow को यह जानने के लिए एक बिना नाम वाले entry workflow की आवश्यकता है कि execution कहाँ से शुरू करनी है। यह entry workflow तुम्हारे नामित वर्कफ़्लो को कॉल करता है।

    - बिना नाम वाला वर्कफ़्लो (entry point)

    ```groovy
    workflow {
        // स्क्रिप्ट चलाने पर यह entry point है
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - नामित वर्कफ़्लो (entry workflow से कॉल किया जाता है)

    ```groovy
    workflow NAMED_WORKFLOW {
        // entry workflow से कॉल किया जाना चाहिए
    }
    ```

4.  **Managing data flow:** हमने सीखा कि namespace notation (`WORKFLOW_NAME.out.channel_name`) का उपयोग करके वर्कफ़्लो आउटपुट तक कैसे पहुँचें और उन्हें अन्य वर्कफ़्लो को कैसे पास करें।

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

[Side Quests के मेनू](../) पर वापस जाओ या अगले विषय पर जाने के लिए पृष्ठ के नीचे दाईं ओर बटन पर क्लिक करो।
