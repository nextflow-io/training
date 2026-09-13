# भाग 2: कंप्यूट संसाधनों और विफलताओं को प्रबंधित करना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[भाग 1](./01_packaging_and_execution.md) में, तुमने यह अनुकूलित किया कि पाइपलाइन के कार्य कहाँ और कैसे चलते हैं।
यहाँ तुम यह अनुकूलित करोगे कि प्रत्येक कार्य को कितना कंप्यूट मिलता है, और जब कोई कार्य तुम्हारे सर्वोत्तम अनुमान के बावजूद विफल हो जाता है तो क्या होता है।

---

## 1. कंप्यूट संसाधन आवंटन को नियंत्रित करना

डिफ़ॉल्ट रूप से, Nextflow `cpus` निर्देश के माध्यम से प्रत्येक प्रोसेस को एक CPU आवंटित करता है, और जब तक तुम एक सेट नहीं करते, मेमोरी सीमा नहीं लगाता:

```groovy title="Built-in configuration"
process {
    cpus = 1
}
```

तुम [Nextflow Run](../nextflow_run/index.md) से पहले से जानते हो कि इस पाइपलाइन का कॉन्फ़िगरेशन सभी प्रोसेस के लिए `memory` को 1 GB सेट करता है।
लेकिन तुम्हें कैसे पता चलेगा कि अपनी खुद की पाइपलाइन के लिए वास्तव में कौन से मान उपयोग करने हैं?

### 1.1. संसाधन उपयोग रिपोर्ट तैयार करना

तुमने [Nextflow Run](../nextflow_run/02_configure_pipeline.md) में `-with-report` के साथ पहले से एक execution रिपोर्ट तैयार की है।
यही रिपोर्ट तुम्हें बताती है कि तुम्हारे प्रोसेस को वास्तव में कितने CPU और मेमोरी की ज़रूरत है: वर्कफ़्लो को कुछ डिफ़ॉल्ट आवंटन के साथ चलाओ, वास्तविक उपयोग रिकॉर्ड करो, फिर वहाँ से समायोजित करो।

```bash
nextflow run main.nf -with-report report-config-1.html
```

रिपोर्ट एक HTML फ़ाइल है जिसे तुम ब्राउज़र में खोल सकते हो।
यह प्रत्येक प्रोसेस के लिए रनटाइम और संसाधन उपयोग को विभाजित करती है, जिसमें यह भी शामिल है कि आवंटित संसाधनों का वास्तव में कितना प्रतिशत उपयोग किया गया।
वर्तमान डिफ़ॉल्ट (1 CPU, 1 GB मेमोरी) के साथ `cowpy` के लिए यह क्या दिखाती है:

| मेट्रिक          | मान    |
| ---------------- | ------ |
| CPU उपयोग        | 116%   |
| पीक मेमोरी उपयोग | 6.4 MB |
| आवंटित मेमोरी    | 1 GB   |

`cowpy` अपने 1 GB आवंटन का 1% से काफी कम उपयोग करता है; 100% से ऊपर `%cpu` का मतलब है कि यह कंटेनर के अंदर संक्षिप्त विस्फोटों में एक CPU के प्रोसेसिंग से अधिक का उपयोग करता है।

उपलब्ध सुविधाओं की पूरी सूची के लिए [Reports](https://nextflow.io/docs/latest/reports.html) देखो।

### 1.2. किसी विशिष्ट प्रोसेस के लिए संसाधन आवंटन सेट करना

ऊपर की रिपोर्ट दिखाती है कि `cowpy` अपने वर्तमान आवंटन के भीतर आराम से है, लेकिन मान लो तुम इसे और अधिक हेडरूम देना चाहते हो, उदाहरण के लिए क्योंकि तुम प्रोडक्शन में बड़े इनपुट की उम्मीद करते हो।
तुम `withName` के साथ एकल प्रोसेस के लिए डिफ़ॉल्ट को ओवरराइड कर सकते हो।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

इसके साथ, प्रत्येक प्रोसेस 1 GB मेमोरी और एक CPU का अनुरोध करता है, सिवाय `cowpy` के, जो 2 GB और 2 CPUs का अनुरोध करता है ([भाग 1](./01_packaging_and_execution.md) से `conda` सेटिंग के अलावा)।

!!! info "जानकारी"

    अगर तुम्हारी मशीन में कम CPU हैं और तुम प्रति प्रोसेस अधिक संख्या आवंटित करते हो, तो कार्य कॉल एक-दूसरे के पीछे कतार में लग सकते हैं, क्योंकि Nextflow उपलब्ध से अधिक CPU का अनुरोध नहीं करेगा।

इसे एक अलग रिपोर्ट फ़ाइलनाम के साथ फिर से चलाओ, ताकि तुम पहले और बाद की तुलना कर सको।

```bash
nextflow run main.nf -with-report report-config-2.html
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [voluminous_venter] revision: c3c85dec78

    executor >  local (8)
    [a1/0e96d4] sayHello (1)       | 3 of 3 ✔
    [a3/7173a3] convertToUpper (2) | 3 of 3 ✔
    [4f/a8ae3d] collectGreetings   | 1 of 1 ✔
    [91/3724f8] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

`cowpy` के लिए दोनों रिपोर्टों की तुलना:

| मेट्रिक          | पहले (1 CPU, 1 GB) | बाद में (2 CPUs, 2 GB) |
| ---------------- | ------------------ | ---------------------- |
| पीक मेमोरी उपयोग | 6.4 MB             | 6.4 MB                 |
| CPU उपयोग        | 116%               | 118%                   |

आवंटन को दोगुना करने से वास्तविक उपयोग बिल्कुल नहीं बदला, जो बताता है कि मूल 1 GB / 1 CPU इस खिलौना वर्कलोड के लिए पहले से ही उदार था।
एक वास्तविक पाइपलाइन पर गैर-तुच्छ डेटा प्रोसेस करते समय, तुम उम्मीद करोगे कि संख्याएँ स्वयं प्रोसेस के बीच सार्थक रूप से भिन्न होंगी, यही कारण है कि तुम अनुमान लगाने के बजाय आवंटित करने का निर्णय लेने से पहले प्रोफ़ाइल करते हो।

### 1.3. संसाधन सीमाएँ जोड़ना

तुम्हारे कंप्यूट इन्फ्रास्ट्रक्चर के आधार पर, तुम जो अनुरोध कर सकते हो उस पर कठोर बाधाएँ हो सकती हैं, उदाहरण के लिए क्लस्टर-व्यापी सीमा।
`resourceLimits` निर्देश तुम्हें वे सीमाएँ सेट करने देता है:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow इन्हें जो भी लक्ष्य executor अपेक्षा करता है उसमें अनुवाद करता है।
अगर कोई प्रोसेस सीमा से अधिक का अनुरोध करता है, तो अनुरोध को अस्वीकार करने के बजाय सीमित कर दिया जाता है।

!!! warning "चेतावनी"

    यह कुछ ऐसा नहीं है जिसे तुम प्रशिक्षण वातावरण में चला सकते हो, क्योंकि इसका प्रभाव होने के लिए HPC इन्फ्रास्ट्रक्चर की आवश्यकता है।

??? info "संस्थागत संदर्भ कॉन्फ़िगरेशन"

    nf-core प्रोजेक्ट दुनिया भर के संस्थानों द्वारा साझा किए गए [कॉन्फ़िगरेशन फ़ाइलों का संग्रह](https://nf-co.re/configs/) बनाए रखता है, जो HPC और क्लाउड executors की एक विस्तृत श्रृंखला को कवर करता है।
    चाहे तुम्हारा अपना संस्थान उनमें से हो या नहीं, ये एक उपयोगी शुरुआती बिंदु हैं।

### सारांश

तुम जानते हो कि संसाधन उपयोग का आकलन करने के लिए प्रोफ़ाइलिंग रिपोर्ट कैसे तैयार करें, किसी विशिष्ट प्रोसेस के लिए संसाधन आवंटन को कैसे ओवरराइड करें, और `resourceLimits` के साथ आवंटन को कैसे सीमित करें।

### आगे क्या है?

जानो कि जब कोई कार्य विफल होता है तो पाइपलाइन को स्वचालित रूप से कैसे पुनर्प्राप्त किया जाए, चाहे तुम्हारा संसाधन आवंटन अनुमान सही हो या नहीं।

---

## 2. रिट्राई के साथ कार्य विफलताओं को संभालना

प्रोफ़ाइलिंग तुम्हें बताती है कि एक प्रोसेस को अधिकांश समय क्या चाहिए, लेकिन वास्तविक वर्कलोड भिन्न होते हैं: एक आवंटन जो अधिकांश इनपुट के लिए आरामदायक है, असामान्य रूप से बड़े इनपुट के लिए अभी भी बहुत तंग हो सकता है, और अनुमान बस गलत हो सकते हैं।
एक विफल कार्य को पूरे रन को बर्बाद करने देने के बजाय, Nextflow एक विफल कार्य को स्वचालित रूप से रिट्राई कर सकता है, वैकल्पिक रूप से प्रत्येक प्रयास पर इसे अधिक संसाधन देता है।

### 2.1. विफल कार्य को स्वचालित रूप से रिट्राई करना

इसे क्रिया में देखने के लिए, जानबूझकर `cowpy` का मेमोरी आवंटन उससे कम सेट करो जो इसे वास्तव में चाहिए: [1.1](#11-generate-a-resource-utilization-report) से याद करो कि यह लगभग 6.4 MB पर पीक करता है, इसलिए 6 MB पर्याप्त से थोड़ा कम होना चाहिए।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

`errorStrategy` Nextflow को बताता है कि जब कोई कार्य विफल होता है तो क्या करना है: `'retry'` पूरी पाइपलाइन को रोकने के बजाय कार्य को फिर से सबमिट करता है।
`maxRetries` यह सीमित करता है कि Nextflow के हार मानने से पहले इसे कितने अतिरिक्त प्रयास मिलते हैं।

```bash
nextflow run main.nf
```

??? failure "कमांड आउटपुट (संक्षिप्त)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [desperate_brazil] revision: c3c85dec78

    executor >  local (10)
    [67/fe1f49] sayHello (1)       | 3 of 3 ✔
    [8a/f13335] convertToUpper (1) | 3 of 3 ✔
    [39/1b24ed] collectGreetings   | 1 of 1 ✔
    [7a/d5eb6f] cowpy              | 0 of 1, retries: 2 ✘
    [9d/b79eb3] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)
    [6d/1d9d84] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (2)
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (137)

    Command executed:
      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      /usr/local/bin/_activate_current_env.sh: line 35:    14 Killed                  micromamba activate "${ENV_NAME:-base}"

    Work dir:
      /workspaces/training/execution-config/work/7a/d5eb6feeac0eed18d95d3da7a7aeb4

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Exit code 137 आउट-ऑफ-मेमोरी किल के लिए मानक संकेत है: कंटेनर में `cowpy` को बिल्कुल भी चलाने के लिए पर्याप्त मेमोरी नहीं थी।
Nextflow ने कार्य को दो बार रिट्राई किया, कुल तीन प्रयास, `maxRetries = 2` से मेल खाते हुए।
चूँकि प्रयासों के बीच मेमोरी आवंटन कभी नहीं बदला, हर प्रयास उसी दीवार से टकराया; एक बार रिट्राई समाप्त हो जाने पर, Nextflow पूरी तरह से विफलता की रिपोर्ट करता है और पाइपलाइन को रोक देता है, गैर-शून्य स्थिति के साथ बाहर निकलता है।

अपने आप रिट्राई करना कुछ भी ठीक नहीं करता अगर प्रयासों के बीच अंतर्निहित कारण नहीं बदलता।

### 2.2. प्रत्येक रिट्राई पर संसाधन बढ़ाना

एक प्रोसेस निर्देश के अंदर, `task.attempt` वर्तमान प्रयास संख्या रखता है, 1 से शुरू होकर।
तुम इसे closure में उपयोग कर सकते हो ताकि प्रत्येक रिट्राई के साथ संसाधन आवंटन को बढ़ाया जा सके।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5 7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = { 6.MB * task.attempt }
            errorStrategy = 'retry'
            maxRetries = 3
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

वर्कफ़्लो को फिर से चलाओ:

```bash
nextflow run main.nf
```

??? success "कमांड आउटपुट (संक्षिप्त)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [grave_joliot] revision: c3c85dec78

    executor >  local (9)
    [0f/211b8a] sayHello (2)       | 3 of 3 ✔
    [26/301ea2] convertToUpper (3) | 3 of 3 ✔
    [22/5a895b] collectGreetings   | 1 of 1 ✔
    [e1/beee86] cowpy              | 1 of 1, retries: 1 ✔
    [b6/7aed6a] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hola-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

पहला प्रयास अभी भी 6 MB पर विफल होता है, लेकिन रिट्राई 12 MB (`6.MB * 2`) के साथ चलता है और सफल होता है, और पाइपलाइन सभी आउटपुट प्रकाशित होने के साथ पूरी होती है।

!!! warning "चेतावनी"

    कंसोल आउटपुट में अभी भी एक `NOTE:` लाइन शामिल है जो विफल पहले प्रयास की रिपोर्ट करती है, भले ही पाइपलाइन समग्र रूप से सफल रही: Nextflow प्रत्येक रिट्राई को व्यक्तिगत रूप से लॉग करता है, लेकिन एक रिट्राई विफलता समग्र परिणाम को प्रभावित नहीं करती।
    यह पुष्टि करने के लिए `Outputs:` सारांश, या कमांड की exit status जाँचो कि रन वास्तव में सफल हुआ या नहीं।

अधिक उन्नत रिट्राई पैटर्न के लिए, जिसमें किस विशिष्ट त्रुटि के आधार पर स्केलिंग शामिल है, Nextflow दस्तावेज़ीकरण में [Dynamic computing resources](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) देखो।

### सारांश

तुम जानते हो कि पाइपलाइन को विफल कार्यों को स्वचालित रूप से रिट्राई कैसे करना है, और `task.attempt` का उपयोग करके प्रत्येक रिट्राई के साथ संसाधन आवंटन को कैसे बढ़ाना है।

### आगे क्या है?

[भाग 3](./03_profiles.md) पर जाओ, जहाँ तुम सीखोगे कि इस तरह के कॉन्फ़िगरेशन को स्विच करने योग्य प्रोफ़ाइल में कैसे बंडल किया जाए।

---

## सारांश

इस भाग में तुमने सीखा:

- संसाधन प्रोफ़ाइलिंग रिपोर्ट तैयार करना और प्रति-प्रोसेस संसाधन आवंटन सेट करना
- `resourceLimits` के साथ संसाधन अनुरोधों को सीमित करना
- `errorStrategy` और `maxRetries` के साथ विफल कार्य को स्वचालित रूप से रिट्राई करना
- `task.attempt` का उपयोग करके प्रत्येक रिट्राई के साथ संसाधन आवंटन को बढ़ाना
