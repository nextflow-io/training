# वर्कफ़्लो के वर्कफ़्लो

जब तुम एक पाइपलाइन विकसित कर रहे होते हो, तो अक्सर तुम विभिन्न डेटा प्रकारों या विश्लेषण चरणों के लिए समान प्रोसेस अनुक्रम बनाते हुए पाते हो। तुम इन प्रोसेस अनुक्रमों को कॉपी और पेस्ट करते हुए समाप्त हो सकते हो, जिससे डुप्लिकेट कोड बनता है जिसे बनाए रखना मुश्किल होता है; या तुम एक विशाल वर्कफ़्लो बना सकते हो जिसे समझना और संशोधित करना कठिन होता है।

Nextflow की सबसे शक्तिशाली विशेषताओं में से एक छोटे, पुन: प्रयोज्य वर्कफ़्लो मॉड्यूल से जटिल पाइपलाइन बनाने की क्षमता है। यह मॉड्यूलर दृष्टिकोण पाइपलाइन को विकसित करना, परीक्षण करना और बनाए रखना आसान बनाता है।

### सीखने के लक्ष्य

इस साइड क्वेस्ट में, हम वर्कफ़्लो मॉड्यूल विकसित करने का पता लगाएंगे जिन्हें अलग से परीक्षण और उपयोग किया जा सकता है, उन मॉड्यूल को एक बड़ी पाइपलाइन में संयोजित करेंगे, और मॉड्यूल के बीच डेटा प्रवाह का प्रबंधन करेंगे।

इस साइड क्वेस्ट के अंत तक, तुम सक्षम होगे:

- जटिल पाइपलाइन को तार्किक, पुन: प्रयोज्य इकाइयों में विभाजित करना
- प्रत्येक वर्कफ़्लो मॉड्यूल का स्वतंत्र रूप से परीक्षण करना
- नई पाइपलाइन बनाने के लिए वर्कफ़्लो को मिक्स और मैच करना
- विभिन्न पाइपलाइन में सामान्य वर्कफ़्लो मॉड्यूल साझा करना
- अपने कोड को अधिक रखरखाव योग्य और समझने में आसान बनाना

ये कौशल तुम्हें स्वच्छ, रखरखाव योग्य कोड संरचना बनाए रखते हुए जटिल पाइपलाइन बनाने में मदद करेंगे।

### पूर्वापेक्षाएँ

इस साइड क्वेस्ट को शुरू करने से पहले तुम्हें:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष शुरुआती पाठ्यक्रम पूरा कर लेना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (प्रोसेस, चैनल, ऑपरेटर, मॉड्यूल) का उपयोग करने में सहज होना चाहिए

---

## 0. शुरू करना

#### ट्रेनिंग कोडस्पेस खोलें

यदि तुमने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित अनुसार ट्रेनिंग वातावरण खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाएं

चलो उस डायरेक्टरी में चलते हैं जहाँ इस ट्यूटोरियल के लिए फ़ाइलें स्थित हैं।

```bash
cd side-quests/workflows_of_workflows
```

तुम VSCode को इस डायरेक्टरी पर फ़ोकस करने के लिए सेट कर सकते हो:

```bash
code .
```

#### सामग्री की समीक्षा करें

तुम्हें एक `modules` डायरेक्टरी मिलेगी जिसमें कई प्रोसेस परिभाषाएँ हैं जो 'Hello Nextflow' में सीखी गई बातों पर आधारित हैं:

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

#### असाइनमेंट की समीक्षा करें

तुम्हारी चुनौती इन मॉड्यूल को दो अलग वर्कफ़्लो में असेंबल करना है जिन्हें हम फिर एक मुख्य वर्कफ़्लो में संयोजित करेंगे:

- एक `GREETING_WORKFLOW` जो नामों को मान्य करता है, अभिवादन बनाता है, और टाइमस्टैम्प जोड़ता है
- एक `TRANSFORM_WORKFLOW` जो टेक्स्ट को अपरकेस में बदलता है और इसे उलट देता है

#### तैयारी चेकलिस्ट

क्या तुम गोता लगाने के लिए तैयार हो?

- [ ] मैं इस पाठ्यक्रम के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता हूँ
- [ ] मेरा कोडस्पेस चल रहा है
- [ ] मैंने अपनी वर्किंग डायरेक्टरी उचित रूप से सेट कर ली है
- [ ] मैं असाइनमेंट को समझता हूँ

यदि तुम सभी बॉक्स चेक कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. Greeting Workflow बनाएं

चलो एक वर्कफ़्लो बनाकर शुरू करते हैं जो नामों को मान्य करता है और टाइमस्टैम्प वाले अभिवादन उत्पन्न करता है।

### 1.1. वर्कफ़्लो संरचना बनाएं

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. पहला (सब)वर्कफ़्लो कोड जोड़ें

यह कोड `workflows/greeting.nf` में जोड़ें:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

यह एक पूर्ण वर्कफ़्लो है, जिसकी संरचना 'Hello Nextflow' ट्यूटोरियल में देखे गए वर्कफ़्लो के समान है, जिसे हम स्वतंत्र रूप से परीक्षण कर सकते हैं। चलो अब इसे आज़माते हैं:

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

यह अपेक्षा के अनुसार काम करता है, लेकिन इसे संयोज्य बनाने के लिए कुछ चीजें हैं जिन्हें हमें बदलने की आवश्यकता है।

### 1.3. वर्कफ़्लो को संयोज्य बनाएं

संयोज्य वर्कफ़्लो में 'Hello Nextflow' ट्यूटोरियल में देखे गए वर्कफ़्लो से कुछ अंतर होते हैं:

- वर्कफ़्लो ब्लॉक को नाम दिया जाना चाहिए
- इनपुट `take:` कीवर्ड का उपयोग करके घोषित किए जाते हैं
- वर्कफ़्लो सामग्री `main:` ब्लॉक के अंदर रखी जाती है
- आउटपुट `emit:` कीवर्ड का उपयोग करके घोषित किए जाते हैं

चलो greeting वर्कफ़्लो को इस संरचना से मेल खाने के लिए अपडेट करते हैं। कोड को निम्नलिखित में बदलें:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Input channel with names

    main:
        // Chain processes: validate -> create greeting -> add timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Original greetings
        timestamped = timestamped_ch  // Timestamped greetings
}
```

तुम देख सकते हो कि वर्कफ़्लो अब नामित है और इसमें `take:` और `emit:` ब्लॉक हैं, और ये वे कनेक्शन हैं जिनका उपयोग हम उच्च स्तरीय वर्कफ़्लो बनाने के लिए करेंगे।
वर्कफ़्लो सामग्री भी `main:` ब्लॉक के अंदर रखी गई है। यह भी ध्यान दें कि हमने `names_ch` इनपुट चैनल घोषणा को हटा दिया है, क्योंकि यह अब वर्कफ़्लो को एक तर्क के रूप में पास किया जाता है।

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

यह तुम्हें एक और नई अवधारणा के बारे में बताता है, एक 'entry workflow'। entry workflow वह वर्कफ़्लो है जो तब कॉल किया जाता है जब तुम एक Nextflow स्क्रिप्ट चलाते हो। डिफ़ॉल्ट रूप से, Nextflow एक अनाम वर्कफ़्लो को entry workflow के रूप में उपयोग करेगा, जब मौजूद हो, और यही तुम अब तक कर रहे हो, वर्कफ़्लो ब्लॉक इस तरह शुरू होते हैं:

```groovy title="hello.nf" linenums="1"
workflow {
```

लेकिन हमारे greeting वर्कफ़्लो में एक अनाम वर्कफ़्लो नहीं है, बल्कि हमारे पास एक नामित वर्कफ़्लो है:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

इसलिए Nextflow ने एक त्रुटि फेंकी और वह नहीं किया जो हम चाहते थे।

हमने `take:`/`emit:` सिंटैक्स इसलिए नहीं जोड़ा ताकि हम वर्कफ़्लो को सीधे कॉल कर सकें - हमने इसे इसलिए किया ताकि हम इसे अन्य वर्कफ़्लो के साथ संयोजित कर सकें। समाधान एक मुख्य स्क्रिप्ट बनाना है जिसमें एक अनाम entry workflow हो जो हमारे नामित वर्कफ़्लो को आयात और कॉल करता है।

### 1.4. मुख्य वर्कफ़्लो बनाएं और परीक्षण करें

अब हम एक मुख्य वर्कफ़्लो बनाएंगे जो `greeting` वर्कफ़्लो को आयात और उपयोग करता है।

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

ध्यान दें कि इस फ़ाइल में हमारी वर्कफ़्लो एंट्री अनाम है, और ऐसा इसलिए है क्योंकि हम इसे entry workflow के रूप में उपयोग करने जा रहे हैं।

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

यह काम करता है! हमने नामित greeting वर्कफ़्लो को एक मुख्य वर्कफ़्लो में लपेट दिया है जिसमें एक अनाम entry `workflow` ब्लॉक है। मुख्य वर्कफ़्लो `GREETING_WORKFLOW` वर्कफ़्लो का उपयोग लगभग (बिल्कुल नहीं) एक प्रोसेस की तरह कर रहा है, और `names` चैनल को एक तर्क के रूप में पास कर रहा है।

### सारांश

इस खंड में, तुमने कई महत्वपूर्ण अवधारणाएँ सीखी हैं:

- **नामित वर्कफ़्लो**: एक नामित वर्कफ़्लो (`GREETING_WORKFLOW`) बनाना जिसे आयात और पुन: उपयोग किया जा सकता है
- **वर्कफ़्लो इंटरफ़ेस**: एक संयोज्य वर्कफ़्लो बनाने के लिए `take:` के साथ स्पष्ट इनपुट और `emit:` के साथ आउटपुट परिभाषित करना
- **एंट्री पॉइंट**: यह समझना कि Nextflow को एक स्क्रिप्ट चलाने के लिए एक अनाम entry workflow की आवश्यकता होती है
- **वर्कफ़्लो संयोजन**: किसी अन्य वर्कफ़्लो के भीतर एक नामित वर्कफ़्लो को आयात और उपयोग करना
- **वर्कफ़्लो नेमस्पेस**: `.out` नेमस्पेस का उपयोग करके वर्कफ़्लो आउटपुट तक पहुँचना (`GREETING_WORKFLOW.out.greetings`)

अब तुम्हारे पास एक कार्यशील greeting वर्कफ़्लो है जो:

- नामों के एक चैनल को इनपुट के रूप में लेता है
- प्रत्येक नाम को मान्य करता है
- प्रत्येक मान्य नाम के लिए एक अभिवादन बनाता है
- अभिवादन में टाइमस्टैम्प जोड़ता है
- मूल और टाइमस्टैम्प वाले दोनों अभिवादन को आउटपुट के रूप में उजागर करता है

यह मॉड्यूलर दृष्टिकोण तुम्हें greeting वर्कफ़्लो का स्वतंत्र रूप से परीक्षण करने या इसे बड़ी पाइपलाइन में एक घटक के रूप में उपयोग करने की अनुमति देता है।

---

## 2. Transform Workflow जोड़ें

अब चलो एक वर्कफ़्लो बनाते हैं जो अभिवादन पर टेक्स्ट ट्रांसफ़ॉर्मेशन लागू करता है।

### 2.1. वर्कफ़्लो फ़ाइल बनाएं

```bash
touch workflows/transform.nf
```

### 2.2. वर्कफ़्लो कोड जोड़ें

यह कोड `workflows/transform.nf` में जोड़ें:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Input channel with messages

    main:
        // Apply transformations in sequence
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Uppercase greetings
        reversed = reversed_ch  // Reversed uppercase greetings
}
```

हम यहाँ संयोज्य सिंटैक्स की व्याख्या दोहराएंगे नहीं, लेकिन ध्यान दें कि नामित वर्कफ़्लो फिर से `take:` और `emit:` ब्लॉक के साथ घोषित किया गया है, और वर्कफ़्लो सामग्री `main:` ब्लॉक के अंदर रखी गई है।

### 2.3. मुख्य वर्कफ़्लो को अपडेट करें

दोनों वर्कफ़्लो का उपयोग करने के लिए `main.nf` को अपडेट करें:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Run the greeting workflow
    GREETING_WORKFLOW(names)

    // Run the transform workflow
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // View results
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

पूर्ण पाइपलाइन चलाएं:

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

यदि तुम उन उलटी फ़ाइलों में से एक पर नज़र डालते हो, तो तुम देखोगे कि यह अभिवादन का अपरकेस संस्करण उलटा है:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### सारांश

अब तुम्हारे पास एक पूर्ण पाइपलाइन होनी चाहिए जो:

- greeting वर्कफ़्लो के माध्यम से नामों को प्रोसेस करती है
- टाइमस्टैम्प वाले अभिवादन को transform वर्कफ़्लो में फीड करती है
- अभिवादन के अपरकेस और उलटे दोनों संस्करण उत्पन्न करती है

---

## सारांश

इस साइड क्वेस्ट में, हमने Nextflow में वर्कफ़्लो संयोजन की शक्तिशाली अवधारणा का पता लगाया है, जो हमें छोटे, पुन: प्रयोज्य घटकों से जटिल पाइपलाइन बनाने की अनुमति देती है।

यह मॉड्यूलर दृष्टिकोण मोनोलिथिक पाइपलाइन की तुलना में कई फायदे प्रदान करता है:

- प्रत्येक वर्कफ़्लो को स्वतंत्र रूप से विकसित, परीक्षण और डीबग किया जा सकता है
- वर्कफ़्लो को विभिन्न पाइपलाइन में पुन: उपयोग किया जा सकता है
- समग्र पाइपलाइन संरचना अधिक पठनीय और रखरखाव योग्य हो जाती है
- एक वर्कफ़्लो में परिवर्तन जरूरी नहीं कि दूसरों को प्रभावित करें यदि इंटरफ़ेस सुसंगत रहते हैं
- एंट्री पॉइंट को आवश्यकतानुसार तुम्हारी पाइपलाइन के विभिन्न भागों को चलाने के लिए कॉन्फ़िगर किया जा सकता है

_यह ध्यान रखना महत्वपूर्ण है कि जबकि वर्कफ़्लो को कॉल करना प्रोसेस को कॉल करने जैसा है, यह वास्तव में एक ही चीज़ नहीं है। उदाहरण के लिए, तुम एक वर्कफ़्लो को N बार नहीं चला सकते हैं इसे N आकार के चैनल के साथ कॉल करके - तुम्हें वर्कफ़्लो को N आकार का चैनल पास करना होगा और आंतरिक रूप से पुनरावृत्ति करनी होगी।_

अपने स्वयं के काम में इन तकनीकों को लागू करने से तुम्हें अधिक परिष्कृत Nextflow पाइपलाइन बनाने में सक्षम होगे जो जटिल बायोइन्फॉर्मेटिक्स कार्यों को संभाल सकती हैं जबकि रखरखाव योग्य और स्केलेबल रहती हैं।

### मुख्य पैटर्न

1.  **वर्कफ़्लो संरचना**: हमने `take:` और `emit:` सिंटैक्स का उपयोग करके प्रत्येक वर्कफ़्लो के लिए स्पष्ट इनपुट और आउटपुट परिभाषित किए, घटकों के बीच अच्छी तरह से परिभाषित इंटरफ़ेस बनाए, और `main:` ब्लॉक के भीतर वर्कफ़्लो लॉजिक को लपेटा।

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Input channels are declared here
            input_ch

        main:
            // Workflow logic goes here
            // This is where processes are called and channels are manipulated
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Output channels are declared here
            output_ch = result_ch
    }
    ```

2.  **वर्कफ़्लो आयात:** हमने दो स्वतंत्र वर्कफ़्लो मॉड्यूल बनाए और उन्हें include स्टेटमेंट के साथ एक मुख्य पाइपलाइन में आयात किया।

    - एकल वर्कफ़्लो शामिल करें

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - एकाधिक वर्कफ़्लो शामिल करें

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - नाम संघर्ष से बचने के लिए उपनाम के साथ शामिल करें

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **एंट्री पॉइंट**: Nextflow को यह जानने के लिए एक अनाम entry workflow की आवश्यकता होती है कि निष्पादन कहाँ से शुरू करना है। यह entry workflow तुम्हारे नामित वर्कफ़्लो को कॉल करता है।

    - अनाम वर्कफ़्लो (entry point)

    ```groovy
    workflow {
        // This is the entry point when the script is run
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - नामित वर्कफ़्लो (entry workflow से कॉल किया गया)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Must be called from the entry workflow
    }
    ```

4.  **डेटा प्रवाह का प्रबंधन:** हमने सीखा कि नेमस्पेस नोटेशन (`WORKFLOW_NAME.out.channel_name`) का उपयोग करके वर्कफ़्लो आउटपुट तक कैसे पहुँचें और उन्हें अन्य वर्कफ़्लो में कैसे पास करें।

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

[साइड क्वेस्ट के मेनू](./index.md) पर वापस जाएं या सूची में अगले विषय पर जाने के लिए पृष्ठ के निचले दाएं कोने में बटन पर क्लिक करें।
