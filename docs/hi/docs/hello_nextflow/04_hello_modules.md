# भाग 4: Hello Modules

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube channel पर [पूरी playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) देखें।

:green_book: वीडियो transcript [यहाँ](./transcripts/04_hello_modules.md) उपलब्ध है।
///
-->

यह section cover करता है कि अपने pipeline की development और maintenance को अधिक efficient और sustainable बनाने के लिए अपने workflow code को कैसे organize करें।
Specifically, हम demonstrate करेंगे कि **modules** कैसे use करें।

Nextflow में, एक **module** एक single process definition है जो एक standalone code file में खुद से encapsulated है।
Workflow में module use करने के लिए, तुम बस अपनी workflow code file में एक single-line import statement add करते हो; फिर तुम process को उसी तरह workflow में integrate कर सकते हो जैसे normally करते।
यह multiple workflows में process definitions को reuse करना possible बनाता है बिना code की multiple copies produce किए।

जब हमने अपना workflow develop करना शुरू किया, हमने सब कुछ एक single code file में लिखा।
अब हम processes को individual modules में move करेंगे।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

यह हमारे code को अधिक shareable, flexible और maintainable बनाएगा।

??? info "इस section से कैसे शुरू करें"

    Course का यह section मानता है कि तुमने [Hello Nextflow](./index.md) course के Parts 1-3 complete कर लिए हैं, लेकिन यदि तुम उन sections में covered basics से comfortable हो, तो तुम बिना कुछ special किए यहाँ से शुरू कर सकते हो।

---

## 0. Warmup: `hello-modules.nf` चलाएं

हम starting point के रूप में workflow script `hello-modules.nf` use करेंगे।
यह इस training course के Part 3 में काम करके produce की गई script के equivalent है, सिवाय इसके कि हमने output destinations बदल दी हैं।

यह sure करने के लिए कि सब कुछ काम कर रहा है, कोई भी changes करने से पहले script को एक बार run करो:

```bash
nextflow run hello-modules.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [hopeful_avogadro] DSL2 - revision: b09af1237d

    executor >  local (7)
    [0f/8795c9] sayHello (3)       [100%] 3 of 3 ✔
    [6a/eb2510] convertToUpper (3) [100%] 3 of 3 ✔
    [af/479117] collectGreetings   [100%] 1 of 1 ✔
    ```

यदि यह तुम्हारे लिए काम किया, तो तुम अपने workflow code को modularize करना सीखने के लिए ready हो।

---

## 1. Modules store करने के लिए directory बनाएं

अपने modules को specific directory में store करना best practice है।
तुम उस directory को कुछ भी कह सकते हो, लेकिन convention इसे `modules/` कहना है।

```bash
mkdir modules
```

!!! tip "सुझाव"

    यहाँ हम तुम्हें दिखा रहे हैं कि **local modules** कैसे use करें, meaning modules जो workflow code के बाकी हिस्से के same repository में locally stored हैं, remote modules के contrast में, जो अन्य (remote) repositories में stored हैं।
    **Remote modules** के बारे में अधिक जानकारी के लिए, [documentation](https://www.nextflow.io/docs/latest/module.html) देखें।

---

## 2. `sayHello()` के लिए module बनाएं

इसके simplest form में, existing process को module में turn करना little more than copy-paste operation है।
हम module के लिए एक file stub बनाएंगे, relevant code copy करेंगे फिर इसे main workflow file से delete करेंगे।

फिर हमें बस एक import statement add करना होगा ताकि Nextflow जान सके कि runtime पर relevant code pull करना है।

### 2.1. New module के लिए file stub बनाएं

चलो `sayHello.nf` नामक module के लिए एक empty file बनाते हैं।

```bash
touch modules/sayHello.nf
```

### 2.2. `sayHello` process code को module file में move करें

Workflow file से पूरी process definition को module file में copy करो, `#!/usr/bin/env nextflow` shebang भी copy करना sure करो।

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

एक बार यह हो जाए, workflow file से process definition delete करो, लेकिन shebang को जगह पर छोड़ना sure करो।

### 2.3. Workflow block से पहले import declaration add करें

Local module import करने के लिए syntax काफी straightforward है:

```groovy title="Syntax: Import declaration"
include { <MODULE_NAME> } from '<path_to_module>'
```

चलो इसे `params` block के ऊपर insert करते हैं और इसे appropriately fill out करते हैं।

=== "After"

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

=== "Before"

    ```groovy title="hello-modules.nf" linenums="44"
    /*
    * Pipeline पैरामीटर
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

### 2.4. Workflow चलाएं

हम पहले जैसे essentially same code और inputs के साथ workflow run कर रहे हैं, तो चलो `-resume` flag के साथ run करते हैं और देखते हैं क्या होता है।

```bash
nextflow run hello-modules.nf -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

यह बहुत जल्दी run होना चाहिए क्योंकि सब कुछ cached है।

Nextflow ने recognize किया कि यह अभी भी same work है, भले ही code multiple files में split हो।

### सीख

तुम जानते हो कि process को local module में extract कैसे करें और तुम जानते हो कि ऐसा करना workflow की resumability break नहीं करता।

### आगे क्या?

More modules बनाने की practice करो।

---

## 3. `convertToUpper()` process को modularize करें

### 3.1. New module के लिए file stub बनाएं

```bash
touch modules/convertToUpper.nf
```

### 3.2. `convertToUpper` process code को module file में move करें

```groovy title="modules/convertToUpper.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * अभिवादन को uppercase में बदलने के लिए text replacement tool का उपयोग करें
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

### 3.3. `params` block से पहले import declaration add करें

```groovy title="hello-modules.nf" linenums="23" hl_lines="3"
// Modules को include करें
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
```

### 3.4. Workflow फिर से चलाएं

```bash
nextflow run hello-modules.nf -resume
```

---

## 4. `collectGreetings()` process को modularize करें

### 4.1. New module के लिए file stub बनाएं

```bash
touch modules/collectGreetings.nf
```

### 4.2. `collectGreetings` process code को module file में move करें

```groovy title="modules/collectGreetings.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Uppercase अभिवादनों को एक single output फ़ाइल में collect करें
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

### 4.3. `params` block से पहले import declaration add करें

```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
// Modules को include करें
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

### 4.4. Workflow चलाएं

```bash
nextflow run hello-modules.nf -resume
```

### सीख

तुम जानते हो कि workflow में multiple processes को modularize कैसे करें।

बधाई हो, तुमने यह सारा काम किया और pipeline कैसे काम करती है उसमें absolutely कुछ भी नहीं बदला!

Jokes aside, अब तुम्हारा code अधिक modular है, और यदि तुम एक और pipeline लिखने का decide करते हो जो उन processes में से किसी को call करती है, तुम्हें relevant module use करने के लिए बस एक short import statement type करना होगा।

### आगे क्या?

जब तुम ready हो, तो [**Part 5: Hello Containers**](./05_hello_containers.md) पर move करो यह सीखने के लिए कि software dependencies को अधिक conveniently और reproducibly manage करने के लिए containers कैसे use करें।

---

## Quiz

<quiz>
Nextflow में module क्या है?
- [ ] एक configuration file
- [x] Single process definition वाली एक standalone file
- [ ] एक workflow definition
- [ ] एक channel operator

और जानें: [2. Create a module for `sayHello()`](#2-create-a-module-for-sayhello)
</quiz>

<quiz>
Module files के लिए recommended naming convention क्या है?
- [ ] `module_processName.nf`
- [ ] `processName_module.nf`
- [x] `processName.nf`
- [ ] `mod_processName.nf`
</quiz>

<quiz>
Module files कहाँ store होनी चाहिए?
- [ ] Workflow के same directory में
- [ ] `bin/` directory में
- [x] `modules/` directory में
- [ ] `lib/` directory में

और जानें: [1. Create a directory to store modules](#1-create-a-directory-to-store-modules)
</quiz>

<quiz>
Module import करने के लिए correct syntax क्या है?

- [ ] `#!groovy import { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy require { SAYHELLO } from './modules/sayhello.nf'`
- [x] `#!groovy include { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy load { SAYHELLO } from './modules/sayhello.nf'`

और जानें: [2.3. Add an import declaration](#23-add-an-import-declaration-before-the-workflow-block)
</quiz>

<quiz>
Modules use करने पर `-resume` functionality का क्या होता है?
- [ ] यह अब काम नहीं करती
- [ ] इसे additional configuration की need होती है
- [x] यह पहले जैसी ही काम करती है
- [ ] यह केवल local modules के लिए काम करती है
</quiz>

<quiz>
Modules use करने के benefits क्या हैं? (सभी लागू select करें)
- [x] Workflows में code reusability
- [x] आसान maintenance
- [x] Workflow code का बेहतर organization
- [ ] तेज़ execution speed
</quiz>
