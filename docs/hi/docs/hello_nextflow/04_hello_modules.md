# भाग 4: Hello Modules

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [पूरी प्लेलिस्ट](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) Nextflow YouTube चैनल पर देखें।

:green_book: वीडियो ट्रांसक्रिप्ट [यहाँ](./transcripts/04_hello_modules.md) उपलब्ध है।
///

यह सेक्शन बताता है कि अपने वर्कफ़्लो कोड को कैसे व्यवस्थित करें ताकि तुम्हारी पाइपलाइन का विकास और रखरखाव अधिक कुशल और टिकाऊ हो सके।
विशेष रूप से, हम यह प्रदर्शित करने जा रहे हैं कि [**modules**](https://nextflow.io/docs/latest/module.html) का उपयोग कैसे करें।

Nextflow में, एक **module** एक स्वतंत्र कोड फ़ाइल है, जो अक्सर एक single process परिभाषा को समाहित करती है।
किसी वर्कफ़्लो में module का उपयोग करने के लिए, तुम बस अपनी वर्कफ़्लो कोड फ़ाइल में एक single-line `include` स्टेटमेंट जोड़ते हो; फिर तुम process को वर्कफ़्लो में उसी तरह integrate कर सकते हो जैसे तुम सामान्य रूप से करते हो।
यह कोड की कई प्रतियाँ बनाए बिना कई वर्कफ़्लो में process परिभाषाओं को पुनः उपयोग करना संभव बनाता है।

जब हमने अपना वर्कफ़्लो विकसित करना शुरू किया, तो हमने सब कुछ एक ही कोड फ़ाइल में लिखा।
अब हम processes को अलग-अलग modules में ले जाने वाले हैं।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

इससे हमारा कोड अधिक शेयर करने योग्य, लचीला और रखरखाव योग्य हो जाएगा।

??? info "इस खंड से कैसे शुरू करें"

    कोर्स का यह सेक्शन मानता है कि तुमने [Hello Nextflow](./index.md) कोर्स के भाग 1-3 पूरे कर लिए हैं, लेकिन अगर तुम उन सेक्शन में शामिल बुनियादी बातों से सहज हो, तो तुम बिना कुछ विशेष किए यहाँ से शुरू कर सकते हो।

---

## 0. वार्मअप: `hello-modules.nf` चलाएँ

हम वर्कफ़्लो स्क्रिप्ट `hello-modules.nf` को शुरुआती बिंदु के रूप में उपयोग करने जा रहे हैं।
यह इस प्रशिक्षण कोर्स के भाग 3 के माध्यम से काम करके बनाई गई स्क्रिप्ट के बराबर है, सिवाय इसके कि हमने आउटपुट गंतव्यों को बदल दिया है:

```groovy title="hello-modules.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_modules'
        mode 'copy'
    }
    uppercased {
        path 'hello_modules'
        mode 'copy'
    }
    collected {
        path 'hello_modules'
        mode 'copy'
    }
    batch_report {
        path 'hello_modules'
        mode 'copy'
    }
}
```

बस यह सुनिश्चित करने के लिए कि सब कुछ काम कर रहा है, कोई भी बदलाव करने से पहले स्क्रिप्ट को एक बार चलाएँ:

```bash
nextflow run hello-modules.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [hopeful_avogadro] DSL2 - revision: b09af1237d

    executor >  local (7)
    [0f/8795c9] sayHello (3)       [100%] 3 of 3 ✔
    [6a/eb2510] convertToUpper (3) [100%] 3 of 3 ✔
    [af/479117] collectGreetings   [100%] 1 of 1 ✔
    ```

पहले की तरह, तुम्हें आउटपुट फ़ाइलें `output` ब्लॉक में निर्दिष्ट डायरेक्टरी में मिलेंगी (यहाँ, `results/hello_modules/`)।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results/hello_modules/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

अगर यह तुम्हारे लिए काम कर गया, तो तुम अपने वर्कफ़्लो कोड को modularize करना सीखने के लिए तैयार हो।

---

## 1. modules स्टोर करने के लिए एक डायरेक्टरी बनाएँ

अपने modules को एक विशिष्ट डायरेक्टरी में स्टोर करना सबसे अच्छा अभ्यास है।
तुम उस डायरेक्टरी को जो चाहो नाम दे सकते हो, लेकिन परंपरा इसे `modules/` कहने की है।

```bash
mkdir modules
```

---

## 2. `sayHello()` के लिए एक module बनाएँ

अपने सबसे सरल रूप में, किसी मौजूदा process को module में बदलना copy-paste ऑपरेशन से थोड़ा अधिक है।
हम module के लिए एक फ़ाइल stub बनाने जा रहे हैं, प्रासंगिक कोड को कॉपी करेंगे फिर इसे मुख्य वर्कफ़्लो फ़ाइल से हटा देंगे।

फिर हमें बस एक `include` स्टेटमेंट जोड़ना होगा ताकि Nextflow को पता चले कि runtime पर प्रासंगिक कोड को pull करना है।

### 2.1. नए module के लिए एक फ़ाइल stub बनाएँ

चलो `sayHello.nf` नाम की module के लिए एक खाली फ़ाइल बनाते हैं।

```bash
touch modules/sayHello.nf
```

यह हमें process कोड रखने के लिए एक जगह देता है।

### 2.2. `sayHello` process कोड को module फ़ाइल में ले जाएँ

पूरी process परिभाषा को वर्कफ़्लो फ़ाइल से module फ़ाइल में कॉपी करें।

```groovy title="modules/sayHello.nf" linenums="1"
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

एक बार यह हो जाने के बाद, वर्कफ़्लो फ़ाइल से process परिभाषा को हटा दें।

### 2.3. workflow ब्लॉक से पहले एक include declaration जोड़ें

किसी module से process को include करने का syntax काफी सीधा है:

```groovy title="Syntax: include declaration"
include { <PROCESS_NAME> } from '<path_to_module>'
```

चलो इसे `params` ब्लॉक के ऊपर insert करें और इसे उचित रूप से भरें।

=== "बाद में"

    ```groovy title="hello-modules.nf" linenums="44" hl_lines="1 2"
    // Modules को include करें
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Pipeline पैरामीटर
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "पहले"

    ```groovy title="hello-modules.nf" linenums="44"
    /*
    * Pipeline पैरामीटर
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

तुम देख सकते हो कि हमने process का नाम, `sayHello`, और module कोड वाली फ़ाइल का path, `./modules/sayHello.nf`, भर दिया है।

### 2.4. वर्कफ़्लो चलाएँ

हम अनिवार्य रूप से पहले जैसे ही कोड और inputs के साथ वर्कफ़्लो चला रहे हैं, तो चलो `-resume` flag के साथ चलाते हैं और देखते हैं कि क्या होता है।

```bash
nextflow run hello-modules.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

यह बहुत जल्दी चलना चाहिए क्योंकि सब कुछ cached है।
प्रकाशित outputs को चेक करने के लिए स्वतंत्र महसूस करो।

Nextflow ने पहचान लिया कि यह अभी भी वही काम है जो किया जाना है, भले ही कोड कई फ़ाइलों में विभाजित हो।

### सारांश

तुम जानते हो कि किसी process को local module में कैसे extract करना है और तुम जानते हो कि ऐसा करने से वर्कफ़्लो की resumability नहीं टूटती।

### आगे क्या है?

और modules बनाने का अभ्यास करो।
एक बार जब तुमने एक कर लिया, तो तुम एक मिलियन और कर सकते हो...
लेकिन चलो अभी के लिए बस दो और करते हैं।

---

## 3. `convertToUpper()` process को modularize करें

### 3.1. नए module के लिए एक फ़ाइल stub बनाएँ

`convertToUpper.nf` नाम की module के लिए एक खाली फ़ाइल बनाएँ।

```bash
touch modules/convertToUpper.nf
```

### 3.2. `convertToUpper` process कोड को module फ़ाइल में ले जाएँ

पूरी process परिभाषा को वर्कफ़्लो फ़ाइल से module फ़ाइल में कॉपी करें।

```groovy title="modules/convertToUpper.nf" linenums="1"
/*
 * greeting को uppercase में बदलने के लिए text replacement tool का उपयोग करें
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

एक बार यह हो जाने के बाद, वर्कफ़्लो फ़ाइल से process परिभाषा को हटा दें।

### 3.3. `params` ब्लॉक से पहले एक include declaration जोड़ें

`params` ब्लॉक के ऊपर include declaration insert करें और इसे उचित रूप से भरें।

=== "बाद में"

    ```groovy title="hello-modules.nf" linenums="23" hl_lines="3"
    // Modules को include करें
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Pipeline पैरामीटर
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "पहले"

    ```groovy title="hello-modules.nf" linenums="23"
    // Modules को include करें
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Pipeline पैरामीटर
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

यह बहुत परिचित दिखना शुरू हो जाना चाहिए।

### 3.4. वर्कफ़्लो फिर से चलाएँ

इसे `-resume` flag के साथ चलाएँ।

```bash
nextflow run hello-modules.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [nauseous_heisenberg] DSL2 - revision: a04a9f2da0

    [c9/763d42] sayHello (3)       | 3 of 3, cached: 3 ✔
    [60/bc6831] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

यह अभी भी पहले जैसा ही आउटपुट उत्पन्न करना चाहिए।

दो हो गए, एक और बाकी है!

---

## 4. `collectGreetings()` process को modularize करें

### 4.1. नए module के लिए एक फ़ाइल stub बनाएँ

`collectGreetings.nf` नाम की module के लिए एक खाली फ़ाइल बनाएँ।

```bash
touch modules/collectGreetings.nf
```

### 4.2. `collectGreetings` process कोड को module फ़ाइल में ले जाएँ

पूरी process परिभाषा को वर्कफ़्लो फ़ाइल से module फ़ाइल में कॉपी करें।

```groovy title="modules/collectGreetings.nf" linenums="1"
/*
 * uppercase greetings को एक single output फ़ाइल में collect करें
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
```

एक बार यह हो जाने के बाद, वर्कफ़्लो फ़ाइल से process परिभाषा को हटा दें।

### 4.3. `params` ब्लॉक से पहले एक include declaration जोड़ें

`params` ब्लॉक के ऊपर include declaration insert करें और इसे उचित रूप से भरें।

=== "बाद में"

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
    // Modules को include करें
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline पैरामीटर
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "पहले"

    ```groovy title="hello-modules.nf" linenums="3"
    // Modules को include करें
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Pipeline पैरामीटर
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

आखिरी वाला!

### 4.4. वर्कफ़्लो चलाएँ

इसे `-resume` flag के साथ चलाएँ।

```bash
nextflow run hello-modules.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [friendly_coulomb] DSL2 - revision: 7aa2b9bc0f

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

यह अभी भी पहले जैसा ही आउटपुट उत्पन्न करना चाहिए।

### सारांश

तुम जानते हो कि किसी वर्कफ़्लो में कई processes को कैसे modularize करना है।

बधाई हो, तुमने यह सारा काम किया और पाइपलाइन के काम करने के तरीके में बिल्कुल कुछ भी नहीं बदला!

मज़ाक को एक तरफ रखते हुए, अब तुम्हारा कोड अधिक modular है, और अगर तुम कोई अन्य पाइपलाइन लिखने का निर्णय लेते हो जो उन processes में से किसी एक को call करती है, तो तुम्हें बस प्रासंगिक module का उपयोग करने के लिए एक छोटा `include` स्टेटमेंट टाइप करना होगा।
यह कोड को copy-paste करने से बेहतर है, क्योंकि अगर बाद में तुम module को सुधारने का निर्णय लेते हो, तो तुम्हारी सभी पाइपलाइनें सुधारों को inherit करेंगी।

### आगे क्या है?

अगर तुम्हें ऐसा लगता है तो थोड़ा ब्रेक लो।

जब तुम तैयार हो, तो [**भाग 5: Hello Containers**](./05_hello_containers.md) पर जाओ ताकि सीख सको कि software dependencies को अधिक सुविधाजनक और reproducibly तरीके से प्रबंधित करने के लिए containers का उपयोग कैसे करें।

---

## क्विज़

<quiz>
Nextflow में module क्या है?
- [ ] एक configuration फ़ाइल
- [x] एक स्वतंत्र फ़ाइल जिसमें process परिभाषाएँ हो सकती हैं
- [ ] एक workflow परिभाषा
- [ ] एक channel operator

और जानें: [2. `sayHello()` के लिए एक module बनाएँ](#2-create-a-module-for-sayhello)
</quiz>

<quiz>
module फ़ाइलों को स्टोर करने के लिए आमतौर पर किस परंपरा का उपयोग किया जाता है?
- [ ] वर्कफ़्लो के समान डायरेक्टरी में
- [ ] एक `bin/` डायरेक्टरी में
- [x] एक `modules/` डायरेक्टरी में
- [ ] एक `lib/` डायरेक्टरी में

और जानें: [1. modules स्टोर करने के लिए एक डायरेक्टरी बनाएँ](#1-create-a-directory-to-store-modules)
</quiz>

<quiz>
module का उपयोग करने के लिए सही syntax क्या है?

- [ ] `#!groovy import { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy require { SAYHELLO } from './modules/sayhello.nf'`
- [x] `#!groovy include { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy load { SAYHELLO } from './modules/sayhello.nf'`

और जानें: [2.3. एक include declaration जोड़ें](#23-add-an-include-declaration-before-the-workflow-block)
</quiz>

<quiz>
modules का उपयोग करते समय `-resume` functionality का क्या होता है?
- [ ] यह अब काम नहीं करती
- [ ] इसे अतिरिक्त configuration की आवश्यकता होती है
- [x] यह पहले की तरह ही काम करती है
- [ ] यह केवल local modules के लिए काम करती है
</quiz>

<quiz>
modules का उपयोग करने के क्या लाभ हैं? (सभी लागू विकल्प चुनें)
- [x] वर्कफ़्लो में कोड की पुनः उपयोगिता
- [x] आसान रखरखाव
- [x] वर्कफ़्लो कोड का बेहतर संगठन
- [ ] तेज़ execution गति
</quiz>
