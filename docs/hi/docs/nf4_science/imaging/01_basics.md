# भाग 1: बुनियादी ऑपरेशन चलाएं

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow for Bioimaging प्रशिक्षण पाठ्यक्रम के इस पहले भाग में, हम एक बहुत ही बुनियादी domain-agnostic Hello World उदाहरण का उपयोग करेंगे ताकि आवश्यक ऑपरेशन प्रदर्शित कर सकें और संबंधित Nextflow कोड घटकों को इंगित कर सकें।

## 1. workflow चलाएं

हम आपको `hello-world.nf` नाम की एक workflow स्क्रिप्ट प्रदान करते हैं जो `--greeting` नाम के command-line argument के माध्यम से इनपुट लेती है और उस greeting वाली एक text फ़ाइल बनाती है।
हम अभी कोड को नहीं देखेंगे; पहले देखते हैं कि इसे चलाना कैसा लगता है।

### 1.1. workflow लॉन्च करें और execution की निगरानी करें

टर्मिनल में, निम्नलिखित कमांड चलाएं:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

आपका console आउटपुट कुछ इस तरह दिखना चाहिए:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

बधाई हो, आपने अभी-अभी अपना पहला Nextflow workflow चलाया!

यहाँ सबसे महत्वपूर्ण आउटपुट आखिरी लाइन है (line 6):

```console title="Output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

यह हमें बताता है कि `sayHello` process एक बार सफलतापूर्वक executed हुआ (`1 of 1 ✔`)।

यह बढ़िया है, लेकिन आप सोच रहे होंगे: आउटपुट कहाँ है?

### 1.2. `results` डायरेक्टरी में आउटपुट फ़ाइल खोजें

यह workflow अपने आउटपुट को `results` नाम की डायरेक्टरी में publish करने के लिए configured है।
यदि आप अपनी वर्तमान डायरेक्टरी को देखें, तो आप देखेंगे कि जब आपने workflow चलाई, तो Nextflow ने `results` नाम की एक नई डायरेक्टरी बनाई, जिसमें `output.txt` नाम की एक फ़ाइल है।

```console title="results/" linenums="1"
results
└── output.txt
```

फ़ाइल खोलें; सामग्री आपके द्वारा command line पर दी गई greeting से मेल खानी चाहिए।

<details>
  <summary>फ़ाइल सामग्री</summary>

```console title="results/output.txt" linenums="1"
Hello World!
```

</details>

यह बढ़िया है, हमारी workflow ने वह किया जो उसे करना था!

हालांकि, ध्यान रखें कि 'published' परिणाम वास्तविक आउटपुट की एक copy (या कुछ मामलों में symlink) है जो Nextflow ने तब उत्पन्न किया जब उसने workflow को execute किया।

तो अब, हम यह देखने के लिए अंदर झांकेंगे कि Nextflow ने वास्तव में कार्य कहाँ execute किया।

!!! warning "चेतावनी"

    सभी workflows अपने आउटपुट को results डायरेक्टरी में publish करने के लिए set up नहीं होंगे, और/या डायरेक्टरी का नाम अलग हो सकता है।
    इस section में थोड़ा आगे, हम आपको दिखाएंगे कि यह व्यवहार कहाँ निर्दिष्ट है।

### 1.3. `work/` डायरेक्टरी में मूल आउटपुट और logs खोजें

जब आप एक workflow चलाते हैं, तो Nextflow workflow में प्रत्येक process के हर एक invocation के लिए एक अलग 'task directory' बनाता है (=pipeline में प्रत्येक चरण)।
प्रत्येक के लिए, यह आवश्यक inputs को stage करेगा, संबंधित instruction(s) को execute करेगा और उस एक डायरेक्टरी के भीतर आउटपुट और log फ़ाइलें लिखेगा, जिसका नाम automatically hash का उपयोग करके दिया जाता है ताकि इसे unique बनाया जा सके।

ये सभी task directories आपकी वर्तमान डायरेक्टरी (जहाँ आप कमांड चला रहे हैं) के अंदर `work` नाम की डायरेक्टरी के अंतर्गत रहेंगी।

यह भ्रमित करने वाला लग सकता है, तो देखते हैं कि यह व्यवहार में कैसा दिखता है।

पहले चलाई गई workflow के console आउटपुट पर वापस जाते हुए, हमारे पास यह लाइन थी:

```console title="Excerpt of command output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

देखें कि लाइन `[a3/7be2fa]` से कैसे शुरू होती है?
यह उस एक process call के लिए task directory path का एक संक्षिप्त रूप है, और आपको बताता है कि `work/` डायरेक्टरी path के भीतर `sayHello` process call का आउटपुट कहाँ मिलेगा।

आप निम्नलिखित कमांड टाइप करके (अपने terminal में दिखाई देने वाले `a3/7be2fa` के साथ बदलते हुए) और path को autocomplete करने के लिए tab key दबाकर या asterisk जोड़कर पूरा path पा सकते हैं:

```bash
tree work/a3/7be2fa*
```

यह पूरा path directory path देना चाहिए: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

देखते हैं कि वहाँ क्या है।

!!! Tip "सुझाव"

    यदि आप VSCode file explorer में task subdirectory की सामग्री ब्राउज़ करते हैं, तो आप सभी फ़ाइलें तुरंत देखेंगे।
    हालांकि, log फ़ाइलें terminal में invisible होने के लिए set हैं, इसलिए यदि आप उन्हें देखने के लिए `ls` या `tree` का उपयोग करना चाहते हैं, तो आपको invisible फ़ाइलें प्रदर्शित करने के लिए संबंधित विकल्प set करना होगा।

    ```bash
    tree -a work
    ```

आपके सिस्टम पर सटीक subdirectory नाम अलग होंगे।

<details>
  <summary>डायरेक्टरी सामग्री</summary>

```console title="work/"
work
└── a3
    └── 7be2fad5e71e5f49998f795677fd68
        ├── .command.begin
        ├── .command.err
        ├── .command.log
        ├── .command.out
        ├── .command.run
        ├── .command.sh
        ├── .exitcode
        └── output.txt
```

</details>

आपको तुरंत `output.txt` फ़ाइल पहचाननी चाहिए, जो वास्तव में `sayHello` process का मूल आउटपुट है जो `results` डायरेक्टरी में published हुआ।
यदि आप इसे खोलते हैं, तो आपको फिर से `Hello World!` greeting मिलेगी।

<details>
  <summary>output.txt की फ़ाइल सामग्री</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

तो उन अन्य सभी फ़ाइलों का क्या?

ये helper और log फ़ाइलें हैं जो Nextflow ने task execution के हिस्से के रूप में लिखीं:

- **`.command.begin`**: Sentinel फ़ाइल जो task launch होते ही बन जाती है।
- **`.command.err`**: Process call द्वारा उत्सर्जित error संदेश (`stderr`)
- **`.command.log`**: Process call द्वारा उत्सर्जित पूर्ण log आउटपुट
- **`.command.out`**: Process call द्वारा नियमित आउटपुट (`stdout`)
- **`.command.run`**: Nextflow द्वारा process call को execute करने के लिए चलाई गई पूर्ण स्क्रिप्ट
- **`.command.sh`**: वह कमांड जो वास्तव में process call द्वारा चलाई गई थी
- **`.exitcode`**: कमांड से परिणामी exit code

`.command.sh` फ़ाइल विशेष रूप से उपयोगी है क्योंकि यह आपको मुख्य कमांड दिखाती है जो Nextflow ने execute किया, सभी bookkeeping और task/environment setup को शामिल नहीं करते हुए।

<details>
  <summary>फ़ाइल सामग्री</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip "सुझाव"

    जब कुछ गलत हो जाता है और आपको troubleshoot करने की आवश्यकता होती है कि क्या हुआ, तो यह देखने के लिए `command.sh` स्क्रिप्ट को देखना उपयोगी हो सकता है कि Nextflow ने workflow instructions, variable interpolation आदि के आधार पर बिल्कुल कौन सी कमांड बनाई।

### 1.4. वैकल्पिक अभ्यास: अलग-अलग greetings के साथ फिर से चलाएं

`--greeting` argument के लिए अलग-अलग मानों के साथ workflow को कुछ बार फिर से चलाने की कोशिश करें, फिर `results/` डायरेक्टरी और task directories दोनों की सामग्री देखें।

देखें कि कैसे isolated task directories के आउटपुट और logs संरक्षित हैं, जबकि `results` डायरेक्टरी की सामग्री बाद के executions के आउटपुट द्वारा overwrite हो जाती है।

### निष्कर्ष

आप जानते हैं कि एक साधारण Nextflow स्क्रिप्ट कैसे चलाएं, इसके execution की निगरानी कैसे करें और इसके आउटपुट कैसे खोजें।

### आगे क्या है?

सीखें कि बुनियादी Nextflow स्क्रिप्ट कैसे पढ़ें और पहचानें कि इसके घटक इसकी functionality से कैसे संबंधित हैं।

---

## 2. Hello World workflow starter स्क्रिप्ट की जांच करें

हमने वहाँ मूल रूप से workflow स्क्रिप्ट को एक black box की तरह treat किया।
अब जब हमने देख लिया है कि यह क्या करती है, तो चलिए box खोलते हैं और अंदर देखते हैं।

_यहाँ लक्ष्य Nextflow कोड की syntax को याद रखना नहीं है, बल्कि यह समझना है कि मुख्य घटक क्या हैं और वे कैसे organized हैं।_

### 2.1. समग्र कोड संरचना की जांच करें

चलिए editor pane में `hello-world.nf` स्क्रिप्ट खोलते हैं।

<details>
  <summary>कोड</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

workflow {

    // एक अभिवादन emit करें
    sayHello(params.greeting)
}
```

</details>

एक Nextflow स्क्रिप्ट में दो मुख्य प्रकार के core घटक शामिल होते हैं: एक या अधिक **processes**, और **workflow** स्वयं।
प्रत्येक **process** वर्णन करता है कि pipeline में संबंधित step को क्या operation(s) पूरा करना चाहिए, जबकि **workflow** dataflow logic का वर्णन करता है जो विभिन्न steps को जोड़ता है।

आइए पहले **process** block को करीब से देखें, फिर हम **workflow** block को देखेंगे।

### 2.2. `process` परिभाषा

कोड का पहला block एक **process** का वर्णन करता है।
Process परिभाषा keyword `process` से शुरू होती है, इसके बाद process name और अंत में curly braces द्वारा सीमांकित process body आती है।
Process body में एक script block होना चाहिए जो चलाने के लिए कमांड निर्दिष्ट करता है, जो कुछ भी हो सकता है जिसे आप command line terminal में चला सकते हैं।

यहाँ हमारे पास `sayHello` नाम का एक **process** है जो `greeting` नामक एक **input** variable लेता है और अपने **output** को `output.txt` नामक फ़ाइल में लिखता है।

<details>
  <summary>कोड</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}
```

</details>

यह एक बहुत ही minimal process परिभाषा है जिसमें केवल एक `input` परिभाषा, एक `output` परिभाषा और execute करने के लिए `script` है।

`input` परिभाषा में `val` qualifier शामिल है, जो Nextflow को किसी प्रकार का मान (string, number, जो भी हो) अपेक्षित करने के लिए कहता है।

`output` परिभाषा में `path` qualifier शामिल है, जो Nextflow को बताता है कि इसे path के रूप में संभाला जाना चाहिए (directory paths और files दोनों शामिल हैं)।

!!! Tip "सुझाव"

    Output परिभाषा _निर्धारित नहीं करती_ कि क्या आउटपुट बनाया जाएगा।
    यह केवल _घोषित करती है_ कि अपेक्षित आउटपुट फ़ाइल(ओं) को कहाँ खोजना है, ताकि execution पूर्ण होने के बाद Nextflow इसे देख सके।

    यह सफलतापूर्वक execute की गई कमांड को verify करने और यदि आवश्यक हो तो downstream processes को आउटपुट पास करने के लिए आवश्यक है।
    उत्पादित आउटपुट जो output block में घोषित के साथ मेल नहीं खाता, downstream processes को पास नहीं किया जाएगा।

वास्तविक दुनिया की pipeline में, एक process में आमतौर पर अतिरिक्त जानकारी होती है जैसे process directives, जिन्हें हम थोड़ी देर में पेश करेंगे।

### 2.3. `workflow` परिभाषा

कोड का दूसरा block **workflow** स्वयं का वर्णन करता है।
Workflow परिभाषा keyword `workflow` से शुरू होती है, इसके बाद एक वैकल्पिक नाम, फिर curly braces द्वारा सीमांकित workflow body आती है।

यहाँ हमारे पास एक **workflow** है जो `sayHello` process की एक call से बना है, जो एक इनपुट लेता है, `params.greeting`, जो उस मान को रखता है जो हमने `--greeting` पैरामीटर को दिया।

```groovy title="hello-world.nf" linenums="22"
workflow {

    // एक अभिवादन emit करें
    sayHello(params.greeting)
}
```

यह एक बहुत ही minimal **workflow** परिभाषा है।
वास्तविक दुनिया की pipeline में, workflow में आमतौर पर **channels** द्वारा जुड़े **processes** की कई calls होती हैं, और variable inputs के लिए default मान set up हो सकते हैं।

हम इसे action में देखेंगे जब हम पाठ्यक्रम के भाग 2 में nf-core/molkart चलाएंगे।

### 2.4. command-line parameters की `params` प्रणाली

`params.greeting` जो हम `sayHello()` process call को प्रदान करते हैं, Nextflow कोड का एक सुव्यवस्थित हिस्सा है और इस पर एक अतिरिक्त मिनट बिताने योग्य है।

जैसा कि ऊपर उल्लेख किया गया है, यही तरीका है जिससे हम `--greeting` command-line पैरामीटर के मान को `sayHello()` process call को पास करते हैं।
वास्तव में, केवल `params.someParameterName` घोषित करने से हम workflow को command-line से `--someParameterName` नामक पैरामीटर दे सकेंगे।

!!! Tip "सुझाव"

    `params` प्रणाली का उपयोग करके घोषित ये workflow parameters हमेशा दो dashes (`--`) लेते हैं।
    यह उन्हें Nextflow-level parameters से अलग करता है, जो केवल एक dash (`-`) लेते हैं।

### निष्कर्ष

अब आप जानते हैं कि एक साधारण Nextflow workflow कैसे structured है, और बुनियादी घटक इसकी functionality से कैसे संबंधित हैं।

### आगे क्या है?

अपने workflow executions को सुविधाजनक रूप से प्रबंधित करना सीखें।

---

## 3. workflow executions प्रबंधित करें

Workflows लॉन्च करना और आउटपुट प्राप्त करना जानना बढ़िया है, लेकिन आप जल्दी पाएंगे कि workflow management के कुछ अन्य पहलू हैं जो आपके जीवन को आसान बना देंगे।

यहाँ हम आपको दिखाते हैं कि जब आपको उसी workflow को फिर से लॉन्च करना हो तो `resume` feature का लाभ कैसे उठाएं, `nextflow log` के साथ execution logs की जांच कैसे करें, और `nextflow clean` के साथ पुरानी work directories को कैसे delete करें।

### 3.1. `-resume` के साथ workflow को फिर से लॉन्च करें

कभी-कभी, आप एक pipeline को फिर से चलाना चाहेंगे जिसे आपने पहले लॉन्च किया था, बिना उस काम को फिर से किए जो पहले ही सफलतापूर्वक पूरा हो चुका है।

Nextflow में `-resume` नामक एक विकल्प है जो आपको ऐसा करने की अनुमति देता है।
विशेष रूप से, इस mode में, कोई भी processes जो पहले से ही बिल्कुल वही कोड, settings और inputs के साथ चलाई जा चुकी हैं, skip की जाएंगी।
इसका मतलब है कि Nextflow केवल उन processes को चलाएगा जिन्हें आपने पिछली बार से जोड़ा या संशोधित किया है, या जिन्हें आप नई settings या inputs प्रदान कर रहे हैं।

ऐसा करने के दो मुख्य फायदे हैं:

- यदि आप एक pipeline विकसित करने के बीच में हैं, तो आप अधिक तेजी से iterate कर सकते हैं क्योंकि आपको अपने परिवर्तनों का परीक्षण करने के लिए केवल उस process(es) को चलाना होगा जिस पर आप सक्रिय रूप से काम कर रहे हैं।
- यदि आप production में एक pipeline चला रहे हैं और कुछ गलत हो जाता है, तो कई मामलों में आप समस्या को ठीक कर सकते हैं और pipeline को फिर से लॉन्च कर सकते हैं, और यह failure के बिंदु से चलना resume करेगा, जो आपका बहुत समय और compute बचा सकता है।

इसका उपयोग करने के लिए, बस अपनी कमांड में `-resume` जोड़ें और इसे चलाएं:

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

Process status line (line 5) में जोड़े गए `cached:` bit को देखें, जिसका अर्थ है कि Nextflow ने पहचान लिया है कि उसने यह काम पहले ही कर लिया है और बस पिछले सफल run के परिणाम का पुन: उपयोग किया।

आप यह भी देख सकते हैं कि work subdirectory hash पिछले run के समान है।
Nextflow सचमुच आपको पिछले execution की ओर इशारा कर रहा है और कह रहा है "मैंने वह पहले ही वहाँ कर दिया।"

!!! Tip "सुझाव"

    जब आप `resume` के साथ pipeline को फिर से चलाते हैं, तो Nextflow किसी भी process call द्वारा `publishDir` डायरेक्टरी में लिखी गई किसी भी फ़ाइल को overwrite नहीं करता है जो पहले सफलतापूर्वक चलाई गई थी।

### 3.2. पिछले executions के log की जांच करें

जब भी आप nextflow workflow लॉन्च करते हैं, तो वर्तमान working डायरेक्टरी में `.nextflow` नामक एक hidden डायरेक्टरी के अंतर्गत `history` नामक log फ़ाइल में एक line लिखी जाती है।

इस जानकारी तक पहुँचने का एक अधिक सुविधाजनक तरीका `nextflow log` कमांड का उपयोग करना है।

```bash
nextflow log
```

यह log फ़ाइल की सामग्री को terminal में आउटपुट करेगा, आपको timestamp, run name, status, और वर्तमान working डायरेक्टरी के भीतर से लॉन्च किए गए हर Nextflow run के लिए पूरी command line दिखाएगा।

### 3.3. पुरानी work directories delete करें

विकास प्रक्रिया के दौरान, आप आमतौर पर अपनी draft pipelines को बड़ी संख्या में चलाएंगे, जो कई subdirectories में बहुत सारी फ़ाइलों के संचय का कारण बन सकता है।
चूंकि subdirectories का नाम randomly दिया जाता है, इसलिए उनके नामों से यह बताना मुश्किल है कि पुरानी बनाम अधिक हाल की runs क्या हैं।

Nextflow में एक सुविधाजनक `clean` subcommand शामिल है जो automatically पिछली runs के लिए work subdirectories को delete कर सकता है जिनकी आपको अब परवाह नहीं है, कई [options](https://www.nextflow.io/docs/latest/reference/cli.html#clean) के साथ यह नियंत्रित करने के लिए कि क्या deleted किया जाएगा।

आप Nextflow log का उपयोग करके timestamp और/या command line के आधार पर एक run को देख सकते हैं, फिर पहले की runs से work directories को delete करने के लिए `nextflow clean -before <run_name> -f` का उपयोग कर सकते हैं।

!!! Warning "चेतावनी"

    पिछली runs से work subdirectories को delete करना उन्हें Nextflow के cache से हटा देता है और उन directories में संग्रहीत किसी भी आउटपुट को delete कर देता है।
    इसका मतलब है कि यह संबंधित processes को फिर से चलाए बिना execution को resume करने की Nextflow की क्षमता को तोड़ देता है।

    किसी भी आउटपुट को सहेजना आपकी जिम्मेदारी है जिसकी आपको परवाह है या जिस पर आप निर्भर रहने की योजना बनाते हैं! यदि आप इस उद्देश्य के लिए `publishDir` directive का उपयोग कर रहे हैं, तो सुनिश्चित करें कि `copy` mode का उपयोग करें, `symlink` mode का नहीं।

### निष्कर्ष

आप जानते हैं कि pipeline को उन steps को दोहराए बिना कैसे फिर से लॉन्च करें जो पहले ही identical तरीके से चलाए गए थे, execution log की जांच करें, और पुरानी work directories को साफ करने के लिए `nextflow clean` कमांड का उपयोग करें।

### आगे क्या है?

अब जब आप बुनियादी Nextflow operations समझते हैं, तो आप nf-core/molkart के साथ एक वास्तविक bioimaging pipeline चलाने के लिए तैयार हैं।
