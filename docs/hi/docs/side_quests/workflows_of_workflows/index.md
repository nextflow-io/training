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

- [Hello Nextflow](../../hello_nextflow/index.md) ट्यूटोरियल या समकक्ष शुरुआती कोर्स पूरा करना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (प्रोसेस, चैनल, ऑपरेटर, मॉड्यूल) का उपयोग करने में सहज होना चाहिए।

---

## 0. शुरू करना

#### ट्रेनिंग codespace खोलो

अगर तुमने अभी तक नहीं किया है, तो [Environment Setup](../../envsetup/index.md) में बताए अनुसार ट्रेनिंग वातावरण खोलना सुनिश्चित करो।

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

एडिटर प्रोजेक्ट डायरेक्टरी पर फ़ोकस के साथ खुलता है।

#### सामग्री की समीक्षा करो

तुम्हें एक `modules` डायरेक्टरी मिलेगी जिसमें प्रोसेस परिभाषाएँ हैं, एक `workflows` डायरेक्टरी जिसमें दो पहले से लिखे गए वर्कफ़्लो स्क्रिप्ट हैं, और एक `main.nf` फ़ाइल जिसे तुम धीरे-धीरे अपडेट करोगे:

```console title="Directory contents"
├── main.nf
├── workflows/
│   ├── greeting.nf              # Standalone greeting workflow (to be made composable)
│   └── transform.nf             # Standalone transform workflow (to be made composable)
└── modules/
    ├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
    ├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
    ├── timestamp_greeting.nf    # Adds timestamps to greetings
    ├── validate_name.nf         # Validates input names
    └── reverse_text.nf          # Reverses text content
```

`modules/` डायरेक्टरी में अलग-अलग प्रोसेस परिभाषाएँ हैं, और `workflows/` डायरेक्टरी में दो पहले से लिखे गए वर्कफ़्लो स्क्रिप्ट हैं जिन पर तुम इस साइड क्वेस्ट में काम करोगे।

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

## 1. Greeting Workflow को पाइपलाइन में जोड़ो

Greeting वर्कफ़्लो नामों को सत्यापित करता है और टाइमस्टैम्प वाले अभिवादन उत्पन्न करता है।

### 1.1. Greeting Workflow की समीक्षा करो और चलाओ

`workflows/greeting.nf` खोलो और कोड देखो:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {
    main:
    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    publish:
    greetings = greetings_ch
    timestamped = timestamped_ch
}

output {
    greetings {
    }
    timestamped {
    }
}
```

यह एक पूर्ण, स्व-निहित वर्कफ़्लो है जिसकी संरचना 'Hello Nextflow' ट्यूटोरियल में देखी गई संरचनाओं के समान है।
यह इनपुट नामों को hardcode करता है, तीन प्रोसेस को जंजीर में जोड़ता है, और दो आउटपुट publish करता है।

इसे चलाकर सत्यापित करो कि सब कुछ काम करता है:

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

इसे अन्य वर्कफ़्लो के साथ composable बनाने के लिए, कुछ चीज़ें बदलनी होंगी।

### 1.2. वर्कफ़्लो को composable बनाओ

एक वर्कफ़्लो को composable बनाने के लिए, चार चीज़ें बदलनी होती हैं:
वर्कफ़्लो को एक नाम मिलता है, इनपुट `take:` ब्लॉक में जाते हैं, आउटपुट `emit:` ब्लॉक में जाते हैं,
और standalone `publish:`/`output {}` ब्लॉक हटा दिए जाते हैं (वे entry workflow में होने चाहिए)।

चलो इन बदलावों को एक-एक करके देखते हैं।

#### 1.2.1. वर्कफ़्लो को नाम दो

वर्कफ़्लो को एक नाम दो ताकि इसे parent workflow से import किया जा सके।

=== "बाद में"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow GREETING_WORKFLOW {
    ```

=== "पहले"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow {
    ```

नाम के साथ, वर्कफ़्लो को अन्य स्क्रिप्ट में import किया जा सकता है।

#### 1.2.2. `take:` के साथ इनपुट घोषित करो

Hardcoded चैनल घोषणा को एक `take:` ब्लॉक से बदलो जो घोषित करे कि वर्कफ़्लो किन इनपुट की अपेक्षा करता है।
`take:` ब्लॉक `main:` से पहले जाता है, और `names_ch = channel.of(...)` लाइन हटा दी जाती है।

=== "बाद में"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="2 3 5"
    workflow GREETING_WORKFLOW {
        take:
        names_ch // नामों के साथ इनपुट चैनल

        main:
        // प्रोसेस को जंजीर में जोड़ो: सत्यापित करो -> अभिवादन बनाओ -> टाइमस्टैम्प जोड़ो
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

=== "पहले"

    ```groovy title="workflows/greeting.nf" linenums="5"
    workflow GREETING_WORKFLOW {
        main:
        names_ch = channel.of('Alice', 'Bob', 'Charlie')

        // प्रोसेस को जंजीर में जोड़ो: सत्यापित करो -> अभिवादन बनाओ -> टाइमस्टैम्प जोड़ो
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

`take:` ब्लॉक केवल नाम से चैनल घोषित करता है — इसमें क्या जाएगा यह parent workflow द्वारा परिभाषित किया जाएगा।

#### 1.2.3. `emit:` के साथ आउटपुट घोषित करो

`publish:` सेक्शन को हटाओ और `output {}` ब्लॉक को भी हटाओ, उन्हें एक `emit:` ब्लॉक से बदलो जो आउटपुट को नाम दे।

=== "बाद में"

    ```groovy title="workflows/greeting.nf" linenums="14" hl_lines="2 3 4"

        emit:
        greetings = greetings_ch // मूल अभिवादन
        timestamped = timestamped_ch // टाइमस्टैम्प वाले अभिवादन
    }
    ```

=== "पहले"

    ```groovy title="workflows/greeting.nf" linenums="14"

        publish:
        greetings = greetings_ch
        timestamped = timestamped_ch
    }

    output {
        greetings {
        }
        timestamped {
        }
    }
    ```

`emit:` ब्लॉक नामित आउटपुट उजागर करता है जिन्हें parent workflow `GREETING_WORKFLOW.out.greetings` और `GREETING_WORKFLOW.out.timestamped` के माध्यम से एक्सेस कर सकते हैं।

#### 1.2.4. परिणाम सत्यापित करो और परीक्षण करो

तीनों बदलावों के बाद, पूरी फ़ाइल इस तरह दिखनी चाहिए:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="5 6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
    names_ch // नामों के साथ इनपुट चैनल

    main:
    // प्रोसेस को जंजीर में जोड़ो: सत्यापित करो -> अभिवादन बनाओ -> टाइमस्टैम्प जोड़ो
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings = greetings_ch // मूल अभिवादन
    timestamped = timestamped_ch // टाइमस्टैम्प वाले अभिवादन
}
```

अब इसे सीधे चलाने की कोशिश करो:

```bash
nextflow run workflows/greeting.nf
```

??? failure "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

यह एक महत्वपूर्ण अवधारणा का परिचय देता है: **entry workflow**।
Nextflow एक बिना नाम वाले `workflow {}` ब्लॉक को entry point के रूप में उपयोग करता है जब तुम किसी स्क्रिप्ट को सीधे चलाते हो।
`GREETING_WORKFLOW` का नाम है, इसलिए Nextflow नहीं जानता कि इसे अपने आप कैसे चलाएं।

यह जानबूझकर है — composable वर्कफ़्लो को entry workflow से कॉल किए जाने के लिए डिज़ाइन किया गया है, सीधे चलाने के लिए नहीं।
समाधान यह है कि `main.nf` में एक entry workflow बनाई जाए जो `GREETING_WORKFLOW` को import और कॉल करे।

### 1.3. मुख्य वर्कफ़्लो अपडेट करो और परीक्षण करो

अब चलो मुख्य वर्कफ़्लो को greeting वर्कफ़्लो कॉल करने के लिए अपडेट करते हैं।

#### 1.3.1. Greeting Workflow को include करो और कॉल करो

`include` स्टेटमेंट जोड़ो, वर्कफ़्लो बॉडी को `GREETING_WORKFLOW` कॉल करने के लिए अपडेट करो और `publish:` में `channel.empty()` placeholder को बदलो:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="1 7 8 11"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // greeting वर्कफ़्लो चलाओ
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        publish:
        greetings = channel.empty()
    }
    ```

Entry workflow बिना नाम का रहता है ताकि Nextflow इसे पाइपलाइन entry point के रूप में उपयोग करे।

#### 1.3.2. Output ब्लॉक अपडेट करो

Published greetings को `greetings/` सबडायरेक्टरी में रूट करने के लिए `path` निर्देश जोड़ो:

=== "बाद में"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="14" hl_lines="2 3"
    output {
        greetings {
        }
    }
    ```

#### 1.3.3. वर्कफ़्लो चलाओ

वर्कफ़्लो चलाओ और परीक्षण करो कि यह काम करता है:

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
    ```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results/
    └── greetings
        ├── Alice-output.txt
        ├── Bob-output.txt
        └── Charlie-output.txt
    ```

??? abstract "फ़ाइल सामग्री"

    ```console title="results/greetings/Alice-output.txt"
    Hello, Alice!
    ```

Greeting फ़ाइलें `results/greetings/` में publish होती हैं।
मुख्य वर्कफ़्लो `GREETING_WORKFLOW` को कॉल करता है और इसके आउटपुट को सीधे `publish:` सेक्शन से जोड़ता है।

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

## 2. Transform Workflow को पाइपलाइन में जोड़ो

Transform वर्कफ़्लो टाइमस्टैम्प वाले अभिवादन पर टेक्स्ट परिवर्तन लागू करता है।

### 2.1. वर्कफ़्लो की समीक्षा करो और चलाओ

`workflows/transform.nf` खोलो और कोड देखो:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped_*.txt')

    // क्रम में परिवर्तन लागू करो
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
```

यह standalone वर्कफ़्लो `greeting.nf` द्वारा उत्पन्न `results/` डायरेक्टरी से टाइमस्टैम्प वाली greeting फ़ाइलें पढ़ता है, उन्हें अपरकेस में बदलता है, फिर टेक्स्ट को उलटता है।

इसे चलाकर सत्यापित करो कि यह सेक्शन 1.1 के greeting परिणामों के साथ काम करता है:

```bash
nextflow run workflows/transform.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/transform.nf` [blissful_curie] DSL2 - revision: 4e7b1c9f02
    executor >  local (6)
    [3e/a14c29] process > SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [c8/51b9e3] process > REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

इसे `GREETING_WORKFLOW` के साथ composable बनाने के लिए, सेक्शन 1.2 के समान तीन बदलाव लागू होते हैं।

### 2.2. इसे composable बनाओ

सेक्शन 1.2 के समान तीन बदलाव लागू करो: वर्कफ़्लो को नाम दो, hardcoded इनपुट को `take:` से बदलो, और `publish:`/`output {}` को `emit:` से बदलो।

पूरी फ़ाइल इस तरह दिखनी चाहिए:

```groovy title="workflows/transform.nf" linenums="1" hl_lines="4 5 6 8 13 14 15"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
    input_ch // संदेशों के साथ इनपुट चैनल

    main:
    // क्रम में परिवर्तन लागू करो
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper = upper_ch // अपरकेस अभिवादन
    reversed = reversed_ch // उलटे अपरकेस अभिवादन
}
```

Transform वर्कफ़्लो अब composable है और मुख्य वर्कफ़्लो में import होने के लिए तैयार है।

### 2.3. मुख्य वर्कफ़्लो अपडेट करो और परीक्षण करो

अब चलो मुख्य वर्कफ़्लो को transformation वर्कफ़्लो कॉल करने के लिए अपडेट करते हैं।

#### 2.3.1. Transformation Workflow को include करो और कॉल करो

Include स्टेटमेंट जोड़ो, टाइमस्टैम्प वाले अभिवादन पर chained `TRANSFORM_WORKFLOW` कॉल जोड़ो, और दो नए `publish:` एंट्री जोड़ो:

=== "बाद में"

    ```groovy title="main.nf" linenums="1" hl_lines="2 11 12 16 17"
    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // greeting वर्कफ़्लो चलाओ
        GREETING_WORKFLOW(names)

        // transform वर्कफ़्लो चलाओ
        TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
        upper = TRANSFORM_WORKFLOW.out.upper
        reversed = TRANSFORM_WORKFLOW.out.reversed
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="1"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // greeting वर्कफ़्लो चलाओ
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

यह टाइमस्टैम्प वाले अभिवादन पर transformation वर्कफ़्लो चलाएगा।

#### 2.3.2. Output ब्लॉक अपडेट करो

`output {}` ब्लॉक में `upper` और `reversed` एंट्री जोड़ो, प्रत्येक के लिए अपनी सबडायरेक्टरी के साथ `path` निर्देश:

=== "बाद में"

    ```groovy title="main.nf" linenums="20" hl_lines="5 6 7 8 9 10"
    output {
        greetings {
            path 'greetings'
        }
        upper {
            path 'upper'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

=== "पहले"

    ```groovy title="main.nf" linenums="20" hl_lines="2 3 4 5"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

यह अंतिम आउटपुट को उचित डायरेक्टरी में publish करेगा।

#### 2.3.3. पूरी पाइपलाइन चलाओ

पाइपलाइन चलाओ और परीक्षण करो कि सब कुछ काम करता है:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results/
    ├── greetings
    │   ├── Alice-output.txt
    │   ├── Bob-output.txt
    │   └── Charlie-output.txt
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── upper
        ├── UPPER-timestamped_Alice-output.txt
        ├── UPPER-timestamped_Bob-output.txt
        └── UPPER-timestamped_Charlie-output.txt
    ```

??? abstract "फ़ाइल सामग्री"

    ```console title="results/reversed/REVERSED-UPPER-timestamped_Alice-output.txt"
    !ECILA ,OLLEH ]04:50:71 60-30-5202[
    ```

पाइपलाइन end-to-end काम कर रही है: अभिवादन को अपरकेस किया गया है और उलटा किया गया है।

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

अपने काम में इन तकनीकों को लागू करने से तुम अधिक परिष्कृत Nextflow पाइपलाइन बना पाओगे जो जटिल डेटा प्रोसेसिंग कार्यों को संभाल सकती हैं और साथ ही रखरखाव योग्य और scalable भी रहती हैं।

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

[Side Quests के मेनू](../index.md) पर वापस जाओ या अगले विषय पर जाने के लिए पृष्ठ के नीचे दाईं ओर बटन पर क्लिक करो।
