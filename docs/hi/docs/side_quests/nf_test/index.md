# nf-test के साथ परीक्षण

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

यह व्यवस्थित रूप से जाँचने में सक्षम होना कि तुम्हारे वर्कफ़्लो का हर हिस्सा वही कर रहा है जो उसे करना चाहिए, पुनरुत्पादनीयता और दीर्घकालिक रखरखाव के लिए बेहद ज़रूरी है, और विकास प्रक्रिया के दौरान यह बहुत मददगार हो सकता है।

आइए एक मिनट के लिए बात करें कि परीक्षण इतना महत्वपूर्ण क्यों है। अगर तुम एक वर्कफ़्लो विकसित कर रहे हो, तो सबसे पहले तुम कुछ टेस्ट डेटा लोगे जो तुम्हें पता है कि वैध है और एक परिणाम देना चाहिए। तुम पाइपलाइन में पहला प्रोसेस जोड़ते हो और इसे अपने इनपुट से जोड़ते हो ताकि यह काम करे। फिर, यह जाँचने के लिए कि सब कुछ काम कर रहा है, तुम इसे टेस्ट डेटा पर चलाते हो। यह मान लेते हुए कि यह काम करता है, तुम अगले प्रोसेस पर जाते हो और टेस्ट डेटा फिर से चलाते हो। तुम यह प्रक्रिया तब तक दोहराते हो जब तक तुम्हारे पास एक ऐसी पाइपलाइन न हो जिससे तुम खुश हो।

फिर, शायद तुम एक सरल true या false पैरामीटर जैसे `--skip_process` जोड़ते हो। अब तुम्हें पाइपलाइन दो बार चलानी होगी, एक बार प्रत्येक पैरामीटर के साथ यह सुनिश्चित करने के लिए कि यह अपेक्षित रूप से काम करती है। लेकिन रुको, हम कैसे जाँचें कि `--skip_process` वास्तव में प्रोसेस को छोड़ता है? हमें आउटपुट में खोजना होगा या लॉग फ़ाइलें जाँचनी होंगी! यह परेशानी भरा है और त्रुटि का खतरा है।

जैसे-जैसे तुम अपनी पाइपलाइन विकसित करते हो, यह जल्दी ही इतनी जटिल हो जाएगी कि हर बदलाव को मैन्युअल रूप से परखना धीमा और त्रुटि-प्रवण हो जाएगा। इसके अलावा, अगर तुम्हें कोई त्रुटि मिलती है तो यह पता लगाना बहुत मुश्किल होगा कि पाइपलाइन में त्रुटि कहाँ से आ रही है। यहीं पर परीक्षण काम आता है।

परीक्षण तुम्हें व्यवस्थित रूप से जाँचने की अनुमति देता है कि तुम्हारी पाइपलाइन का हर हिस्सा अपेक्षित रूप से काम कर रहा है। अच्छी तरह से लिखे गए परीक्षणों के एक डेवलपर को बहुत फायदे हैं:

- **आत्मविश्वास**: क्योंकि परीक्षण पूरी पाइपलाइन को कवर करते हैं, तुम आश्वस्त हो सकते हो कि कुछ बदलने से कुछ और प्रभावित नहीं होता
- **विश्वास**: जब कई डेवलपर पाइपलाइन पर काम करते हैं, तो वे जानते हैं कि दूसरे डेवलपर ने पाइपलाइन और हर घटक को नहीं तोड़ा है
- **पारदर्शिता**: परीक्षण दिखाते हैं कि पाइपलाइन कहाँ विफल हो रही है और समस्या को ट्रैक करना आसान बनाते हैं। वे दस्तावेज़ीकरण के एक रूप के रूप में भी काम करते हैं, यह दिखाते हुए कि किसी प्रोसेस या वर्कफ़्लो को कैसे चलाया जाए
- **गति**: क्योंकि परीक्षण स्वचालित हैं, उन्हें बहुत जल्दी और बार-बार चलाया जा सकता है। तुम नए बग्स पेश करने के कम डर के साथ तेज़ी से काम कर सकते हो

हम कई अलग-अलग प्रकार के परीक्षण लिख सकते हैं:

1. **मॉड्यूल-स्तरीय परीक्षण**: अलग-अलग प्रोसेस के लिए
2. **वर्कफ़्लो-स्तरीय परीक्षण**: एकल वर्कफ़्लो के लिए
3. **पाइपलाइन-स्तरीय परीक्षण**: पूरी पाइपलाइन के लिए
4. **प्रदर्शन परीक्षण**: पाइपलाइन की गति और दक्षता के लिए
5. **स्ट्रेस परीक्षण**: इसकी सीमाएं निर्धारित करने के लिए अत्यधिक परिस्थितियों में पाइपलाइन के प्रदर्शन का आकलन

अलग-अलग प्रोसेस का परीक्षण करना अन्य भाषाओं में यूनिट परीक्षणों के समान है। वर्कफ़्लो या पूरी पाइपलाइन का परीक्षण करना अन्य भाषाओं में इंटीग्रेशन परीक्षण के समान है, जहाँ हम घटकों की परस्पर क्रियाओं का परीक्षण करते हैं।

[**nf-test**](https://www.nf-test.com/) एक ऐसा टूल है जो तुम्हें मॉड्यूल, वर्कफ़्लो और पाइपलाइन स्तर के परीक्षण लिखने की अनुमति देता है। संक्षेप में, यह तुम्हें व्यवस्थित रूप से जाँचने की अनुमति देता है कि पाइपलाइन का हर अलग हिस्सा अपेक्षित रूप से काम कर रहा है, _अलगाव में_।

### सीखने के लक्ष्य

इस साइड क्वेस्ट में, तुम पाइपलाइन के लिए वर्कफ़्लो-स्तरीय परीक्षण और इसके द्वारा बुलाए जाने वाले तीन प्रोसेस के लिए मॉड्यूल-स्तरीय परीक्षण लिखने के लिए nf-test का उपयोग करना सीखोगे।

इस साइड क्वेस्ट के अंत तक, तुम निम्नलिखित तकनीकों का प्रभावी ढंग से उपयोग करने में सक्षम होगे:

- अपने प्रोजेक्ट में nf-test को इनिशियलाइज़ करना
- मॉड्यूल-स्तरीय और वर्कफ़्लो-स्तरीय परीक्षण जनरेट करना
- सामान्य प्रकार के assertions जोड़ना
- यह समझना कि snapshots बनाम content assertions का उपयोग कब करना है
- पूरे प्रोजेक्ट के लिए परीक्षण चलाना

ये कौशल तुम्हें अपनी पाइपलाइन परियोजनाओं में एक व्यापक परीक्षण रणनीति लागू करने में मदद करेंगे, यह सुनिश्चित करते हुए कि वे अधिक मज़बूत और रखरखाव योग्य हैं।

### पूर्वापेक्षाएँ

इस साइड क्वेस्ट को शुरू करने से पहले, तुम्हें:

- [Hello Nextflow](../hello_nextflow/README.md) ट्यूटोरियल या समकक्ष शुरुआती कोर्स पूरा करना चाहिए।
- बुनियादी Nextflow अवधारणाओं और तंत्रों (प्रोसेस, चैनल, ऑपरेटर, फ़ाइलों के साथ काम करना, मेटा डेटा) का उपयोग करने में सहज होना चाहिए।

---

## 0. शुरू करना

#### ट्रेनिंग codespace खोलो

अगर तुमने अभी तक ऐसा नहीं किया है, तो सुनिश्चित करो कि [पर्यावरण सेटअप](../envsetup/index.md) में वर्णित अनुसार ट्रेनिंग वातावरण खोलो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### प्रोजेक्ट डायरेक्टरी में जाओ

आइए उस डायरेक्टरी में जाएं जहाँ इस ट्यूटोरियल की फ़ाइलें स्थित हैं।

```bash
cd side-quests/nf-test
```

तुम VSCode को इस डायरेक्टरी पर फ़ोकस करने के लिए सेट कर सकते हो:

```bash
code .
```

#### सामग्री की समीक्षा करो

तुम्हें एक main वर्कफ़्लो फ़ाइल और `greetings.csv` नामक एक CSV फ़ाइल मिलेगी जिसमें पाइपलाइन का इनपुट है।

```console title="Directory contents"
.
├── greetings.csv
└── main.nf
```

फ़ाइलों के विस्तृत विवरण के लिए, [Hello Nextflow का वार्मअप](../hello_nextflow/00_orientation.md) देखो।

जिस वर्कफ़्लो का हम परीक्षण करेंगे वह [Hello Workflow](../hello_nextflow/03_hello_workflow.md) में बनाए गए Hello वर्कफ़्लो का एक उपसमूह है।

??? example "Hello Nextflow वर्कफ़्लो क्या करता है?"

    अगर तुमने [Hello Nextflow](../hello_nextflow/index.md) प्रशिक्षण नहीं किया है, तो यहाँ एक त्वरित अवलोकन है कि यह सरल वर्कफ़्लो क्या करता है।

    वर्कफ़्लो एक CSV फ़ाइल लेता है जिसमें अभिवादन हैं, उन पर चार लगातार परिवर्तन चरण चलाता है, और एक एकल टेक्स्ट फ़ाइल आउटपुट करता है जिसमें एक मज़ेदार पात्र का ASCII चित्र है जो अभिवादन कह रहा है।

    चार चरण Nextflow प्रोसेस (`sayHello`, `convertToUpper`, `collectGreetings`, और `cowpy`) के रूप में लागू किए गए हैं जो अलग मॉड्यूल फ़ाइलों में संग्रहीत हैं।

    1. **`sayHello`:** प्रत्येक अभिवादन को अपनी आउटपुट फ़ाइल में लिखता है (जैसे, "Hello-output.txt")
    2. **`convertToUpper`:** प्रत्येक अभिवादन को अपरकेस में बदलता है (जैसे, "HELLO")
    3. **`collectGreetings`:** सभी अपरकेस अभिवादनों को एक एकल बैच फ़ाइल में एकत्र करता है
    4. **`cowpy`:** `cowpy` टूल का उपयोग करके ASCII आर्ट जनरेट करता है

    परिणाम `results/` नामक डायरेक्टरी में प्रकाशित किए जाते हैं, और पाइपलाइन का अंतिम आउटपुट (डिफ़ॉल्ट पैरामीटर के साथ चलाने पर) एक सादा टेक्स्ट फ़ाइल है जिसमें अपरकेस अभिवादन कहने वाले एक पात्र की ASCII आर्ट है।

    इस साइड क्वेस्ट में, हम Hello वर्कफ़्लो के एक मध्यवर्ती रूप का उपयोग करते हैं जिसमें केवल पहले दो प्रोसेस हैं। <!-- TODO: change this to use the full finished workflow as suggested in https://github.com/nextflow-io/training/issues/735 -->

जिस उपसमूह पर हम काम करेंगे वह दो प्रोसेस से बना है: `sayHello` और `convertToUpper`।
तुम नीचे पूरा वर्कफ़्लो कोड देख सकते हो।

??? example "वर्कफ़्लो कोड"

    ```groovy title="main.nf"
    /*
    * पाइपलाइन पैरामीटर
    */
    params.input_file = "greetings.csv"

    /*
    * 'Hello World!' को standard out पर प्रिंट करने के लिए echo का उपयोग करें
    */
    process sayHello {

        publishDir 'results', mode: 'copy'

        input:
            val greeting

        output:
            path "${greeting}-output.txt"

        script:
        """
        echo '$greeting' > '$greeting-output.txt'
        """
    }

    /*
    * अभिवादन को अपरकेस में बदलने के लिए text replace utility का उपयोग करें
    */
    process convertToUpper {

        publishDir 'results', mode: 'copy'

        input:
            path input_file

        output:
            path "UPPER-${input_file}"

        script:
        """
        cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
        """
    }

    workflow {

        // CSV फ़ाइल से इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        // अभिवादन को अपरकेस में बदलें
        convertToUpper(sayHello.out)
    }
    ```

#### वर्कफ़्लो चलाओ

आइए वर्कफ़्लो चलाएं यह सुनिश्चित करने के लिए कि यह अपेक्षित रूप से काम कर रहा है।

```bash
nextflow run main.nf
```

```console title="Result of running the workflow"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

बधाई हो! तुमने अभी एक परीक्षण चलाया!

"रुको, क्या? मैंने बस वर्कफ़्लो चलाया और यह काम किया! यह परीक्षण कैसे है?"

अच्छा सवाल!

आइए समझें कि अभी क्या हुआ।

तुमने डिफ़ॉल्ट पैरामीटर के साथ वर्कफ़्लो चलाया, तुमने पुष्टि की कि यह काम किया और तुम परिणामों से खुश हो। यही परीक्षण का सार है। अगर तुमने Hello Nextflow प्रशिक्षण कोर्स के माध्यम से काम किया है, तो तुमने देखा होगा कि हमने हमेशा हर अनुभाग की शुरुआत उस वर्कफ़्लो को चलाकर की जिसे हम शुरुआती बिंदु के रूप में उपयोग कर रहे थे, यह पुष्टि करने के लिए कि सब कुछ सही तरीके से सेट है।

परीक्षण सॉफ़्टवेयर अनिवार्य रूप से हमारे लिए यह प्रक्रिया करता है।

#### असाइनमेंट की समीक्षा करो

तुम्हारी चुनौती nf-test का उपयोग करके इस वर्कफ़्लो में मानकीकृत परीक्षण जोड़ना है, ताकि यह सत्यापित करना आसान हो सके कि अगर कोई और बदलाव किए जाते हैं तो हर हिस्सा अपेक्षित रूप से काम करता रहे।

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### तैयारी की जाँच सूची

क्या तुम शुरू करने के लिए तैयार हो?

- [ ] मैं इस कोर्स के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता/समझती हूँ
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी working directory उचित रूप से सेट की है
- [ ] मैंने वर्कफ़्लो सफलतापूर्वक चलाया है
- [ ] मैं असाइनमेंट समझता/समझती हूँ

अगर तुम सभी बॉक्स चेक कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. `nf-test` को इनिशियलाइज़ करो

`nf-test` पैकेज एक initialization कमांड प्रदान करता है जो हमारे प्रोजेक्ट के लिए परीक्षण विकसित करना शुरू करने के लिए कुछ चीज़ें सेट करता है।

```bash
nf-test init
```

इससे निम्नलिखित आउटपुट मिलना चाहिए:

```bash
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

यह एक configuration फ़ाइल stub वाली `tests` डायरेक्टरी भी बनाता है।

### 1.1. एक nf-test जनरेट करो

`nf-test` में nf-test फ़ाइलें बनाने के लिए टूल का एक सेट आता है, जो हमारे अधिकांश काम को बचाता है। ये `generate` सबकमांड के अंतर्गत आते हैं। आइए पाइपलाइन के लिए एक परीक्षण जनरेट करें:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

यह `tests` डायरेक्टरी के भीतर एक `main.nf.test` फ़ाइल बनाएगा। यह हमारी पाइपलाइन स्तर की परीक्षण फ़ाइल है। अगर तुम `tree tests/` चलाते हो तो तुम्हें कुछ ऐसा दिखना चाहिए:

```console title="Test directory contents"
tests/
├── main.nf.test
└── nextflow.config
```

`main.nf.test` फ़ाइल हमारी पाइपलाइन स्तर की परीक्षण फ़ाइल है। आइए इसे खोलें और सामग्री देखें।

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

    test("Should run without failures") {

        when {
            params {
                // यहाँ पैरामीटर परिभाषित करें। उदाहरण:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

आइए परीक्षण फ़ाइल की संरचना को समझने के लिए एक पल लें।

`nextflow_pipeline` ब्लॉक सभी पाइपलाइन स्तर के परीक्षणों के लिए प्रवेश बिंदु है। इसमें निम्नलिखित शामिल हैं:

- `name`: परीक्षण का नाम।
- `script`: पाइपलाइन स्क्रिप्ट का पथ।

`test` ब्लॉक वास्तविक परीक्षण है। इसमें निम्नलिखित शामिल हैं:

- `when`: वे शर्तें जिनके तहत परीक्षण चलाया जाना चाहिए। इसमें वे पैरामीटर शामिल हैं जो पाइपलाइन चलाने के लिए उपयोग किए जाएंगे।
- `then`: वे assertions जो की जानी चाहिए। इसमें पाइपलाइन के अपेक्षित परिणाम शामिल हैं।

सरल भाषा में, परीक्षण का तर्क इस प्रकार पढ़ा जाता है:
"**जब** ये _पैरामीटर_ इस _पाइपलाइन_ को प्रदान किए जाते हैं, **तब** हम इन परिणामों की अपेक्षा करते हैं।"

यह एक कार्यात्मक परीक्षण नहीं है, हम अगले अनुभाग में दिखाएंगे कि इसे कैसे बनाया जाए।

### परीक्षण नामों पर एक नोट

ऊपर के उदाहरण में, हमने डिफ़ॉल्ट नाम "Should run without failures" का उपयोग किया जो एक बुनियादी परीक्षण के लिए उपयुक्त है जो केवल जाँचता है कि पाइपलाइन सफलतापूर्वक चलती है। हालाँकि, जैसे-जैसे हम अधिक विशिष्ट परीक्षण मामले जोड़ते हैं, हमें अधिक वर्णनात्मक नामों का उपयोग करना चाहिए जो यह दर्शाते हों कि हम वास्तव में क्या परीक्षण कर रहे हैं। उदाहरण के लिए:

- "Should convert input to uppercase" - जब विशिष्ट कार्यक्षमता का परीक्षण करना हो
- "Should handle empty input gracefully" - जब edge cases का परीक्षण करना हो
- "Should respect max memory parameter" - जब resource constraints का परीक्षण करना हो
- "Should create expected output files" - जब फ़ाइल जनरेशन का परीक्षण करना हो

अच्छे परीक्षण नामों को:

1. "Should" से शुरू होना चाहिए ताकि यह स्पष्ट हो कि अपेक्षित व्यवहार क्या है
2. परीक्षण की जा रही विशिष्ट कार्यक्षमता या परिदृश्य का वर्णन करना चाहिए
3. इतना स्पष्ट होना चाहिए कि अगर परीक्षण विफल हो जाए, तो तुम जान सको कि कौन सी कार्यक्षमता टूटी है

जैसे-जैसे हम बाद में अधिक assertions और विशिष्ट परीक्षण मामले जोड़ते हैं, हम इन अधिक वर्णनात्मक नामों का उपयोग करेंगे ताकि यह स्पष्ट हो कि प्रत्येक परीक्षण क्या सत्यापित कर रहा है।

### 1.2. परीक्षण चलाओ

आइए परीक्षण चलाएं और देखें क्या होता है।

```bash
nf-test test tests/main.nf.test
```

```console title="nf-test pipeline fail"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures' FAILED (4.652s)

  Assertion failed:

  assert workflow.success
         |        |
         workflow false

  Nextflow stdout:

  ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv

   -- Check '/workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/meta/nextflow.log' file for details
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.679s (1 failed)
```

परीक्षण विफल हो गया! क्या हुआ?

1. nf-test ने `when` ब्लॉक में सेटिंग्स का उपयोग करके पाइपलाइन को चलाने की कोशिश की:

```groovy title="tests/main.nf.test"
when {
    params {
        // यहाँ पैरामीटर परिभाषित करें। उदाहरण:
        // outdir = "tests/results"
    }
}
```

2. nf-test ने पाइपलाइन की स्थिति जाँची और इसे `when` ब्लॉक से तुलना की:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

ध्यान दो कि nf-test ने बताया कि पाइपलाइन विफल हुई और Nextflow से त्रुटि संदेश प्रदान किया:

```console title="Error"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

तो समस्या क्या थी? याद करो कि पाइपलाइन में प्रोजेक्ट डायरेक्टरी में एक greetings.csv फ़ाइल है। जब nf-test पाइपलाइन चलाता है, तो यह इस फ़ाइल को खोजेगा, लेकिन यह इसे नहीं ढूंढ सकता। फ़ाइल वहाँ है, तो क्या हो रहा है? खैर, अगर हम पथ को देखें तो हम देख सकते हैं कि परीक्षण `./nf-test/tests/longHashString/` पथ में हो रहा है। Nextflow की तरह, nf-test हर परीक्षण के लिए एक नई डायरेक्टरी बनाता है ताकि सब कुछ अलग रहे। डेटा फ़ाइल वहाँ नहीं है इसलिए हमें मूल परीक्षण में फ़ाइल के पथ को सही करना होगा।

आइए परीक्षण फ़ाइल पर वापस जाएं और `when` ब्लॉक में फ़ाइल के पथ को बदलें।

तुम सोच रहे होगे कि हम परीक्षण में पाइपलाइन की root को कैसे इंगित करेंगे। चूँकि यह एक सामान्य स्थिति है, nf-test में global variables की एक श्रृंखला है जिनका उपयोग हम अपना काम आसान बनाने के लिए कर सकते हैं। तुम पूरी सूची [यहाँ](https://www.nf-test.com/docs/testcases/global_variables/) पा सकते हो, लेकिन अभी के लिए हम `projectDir` variable का उपयोग करेंगे, जिसका अर्थ है पाइपलाइन प्रोजेक्ट की root।

=== "बाद में"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
    when {
        params {
            input_file = "${projectDir}/greetings.csv"
        }
    }
    ```

=== "पहले"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
    when {
        params {
            // यहाँ पैरामीटर परिभाषित करें। उदाहरण:
            // outdir = "tests/results"
        }
    }
    ```

आइए परीक्षण फिर से चलाएं और देखें कि यह काम करता है या नहीं।

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

सफलता! पाइपलाइन सफलतापूर्वक चलती है और परीक्षण पास होता है। इसे जितनी बार चाहो चलाओ और तुम्हें हमेशा वही परिणाम मिलेगा!

डिफ़ॉल्ट रूप से, Nextflow आउटपुट छिपा होता है, लेकिन खुद को यह विश्वास दिलाने के लिए कि nf-test निश्चित रूप से वर्कफ़्लो चला रहा है, तुम `--verbose` flag का उपयोग कर सकते हो:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Pipeline runs all processes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures'
    > Nextflow 24.10.4 is available - Please consider updating your version to it
    > N E X T F L O W  ~  version 24.10.0
    > Launching `/workspaces/training/side-quests/nf-test/main.nf` [zen_ampere] DSL2 - revision: bbf79d5c31
    > [2b/61e453] Submitted process > sayHello (2)
    > [31/4e1606] Submitted process > sayHello (1)
    > [bb/5209ee] Submitted process > sayHello (3)
    > [83/83db6f] Submitted process > convertToUpper (2)
    > [9b/3428b1] Submitted process > convertToUpper (1)
    > [ca/0ba51b] Submitted process > convertToUpper (3)
    PASSED (5.206s)


SUCCESS: Executed 1 tests in 5.239s
```

### 1.3. Assertions जोड़ो

एक सरल जाँच यह सुनिश्चित करना है कि हमारी पाइपलाइन उन सभी प्रोसेस को चला रही है जिनकी हम अपेक्षा करते हैं और चुपचाप किसी को छोड़ नहीं रही है। याद करो कि हमारी पाइपलाइन 6 प्रोसेस चलाती है, एक `sayHello` और एक `convertToUpper` प्रत्येक 3 अभिवादनों के लिए।

आइए अपने परीक्षण में एक assertion जोड़ें यह जाँचने के लिए कि पाइपलाइन अपेक्षित संख्या में प्रोसेस चलाती है। हम परीक्षण का नाम भी अपडेट करेंगे ताकि यह बेहतर ढंग से दर्शाए कि हम क्या परीक्षण कर रहे हैं।

=== "बाद में"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="1 11"
        test("Should run successfully with correct number of processes") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert workflow.success
                assert workflow.trace.tasks().size() == 6
            }

        }
    ```

=== "पहले"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="1"
        test("Should run without failures") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert workflow.success
            }

        }
    ```

परीक्षण का नाम अब बेहतर ढंग से दर्शाता है कि हम वास्तव में क्या सत्यापित कर रहे हैं - न केवल यह कि पाइपलाइन बिना विफल हुए चलती है, बल्कि यह कि यह अपेक्षित संख्या में प्रोसेस चलाती है।

आइए परीक्षण फिर से चलाएं और देखें कि यह काम करता है या नहीं।

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with assertions"
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

सफलता! पाइपलाइन सफलतापूर्वक चलती है और परीक्षण पास होता है। अब हमने पाइपलाइन के विवरण के साथ-साथ समग्र स्थिति का परीक्षण करना शुरू कर दिया है।

### 1.4. आउटपुट का परीक्षण करो

आइए अपने परीक्षण में एक assertion जोड़ें यह जाँचने के लिए कि आउटपुट फ़ाइल बनाई गई थी। हम इसे एक अलग परीक्षण के रूप में जोड़ेंगे, एक सूचनाप्रद नाम के साथ, ताकि परिणामों की व्याख्या करना आसान हो।

=== "बाद में"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="14-33"
        test("Should run successfully with correct number of processes") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert workflow.success
                assert workflow.trace.tasks().size() == 6
            }

        }

        test("Should produce correct output files") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert file("$launchDir/results/Bonjour-output.txt").exists()
                assert file("$launchDir/results/Hello-output.txt").exists()
                assert file("$launchDir/results/Holà-output.txt").exists()
                assert file("$launchDir/results/UPPER-Bonjour-output.txt").exists()
                assert file("$launchDir/results/UPPER-Hello-output.txt").exists()
                assert file("$launchDir/results/UPPER-Holà-output.txt").exists()
            }

        }
    ```

=== "पहले"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="14"
        test("Should run successfully with correct number of processes") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert workflow.success
                assert workflow.trace.tasks().size() == 6
            }

        }
    ```

परीक्षण फिर से चलाओ और देखो कि यह काम करता है या नहीं।

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with file assertions"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

सफलता! परीक्षण पास होते हैं क्योंकि पाइपलाइन सफलतापूर्वक पूरी हुई, सही संख्या में प्रोसेस चले और आउटपुट फ़ाइलें बनाई गईं। यह तुम्हें यह भी दिखाना चाहिए कि अपने परीक्षणों के लिए वे सूचनाप्रद नाम प्रदान करना कितना उपयोगी है।

यह तो बस शुरुआत है, हम पाइपलाइन के विवरण की जाँच करने के लिए assertions लिखते रह सकते हैं, लेकिन अभी के लिए आइए पाइपलाइन के आंतरिक भागों के परीक्षण पर आगे बढ़ें।

### सारांश

तुम जानते हो कि पाइपलाइन के लिए nf-test कैसे लिखें।

### आगे क्या है?

Nextflow प्रोसेस का परीक्षण करना सीखो।

---

## 2. Nextflow प्रोसेस का परीक्षण करो

हमें पाइपलाइन के हर हिस्से के लिए परीक्षण लिखने की ज़रूरत नहीं है, लेकिन जितने अधिक परीक्षण हमारे पास होंगे, उतना ही अधिक हम पाइपलाइन के बारे में व्यापक हो सकते हैं और उतना ही अधिक आश्वस्त हो सकते हैं कि यह अपेक्षित रूप से काम कर रहा है। इस अनुभाग में हम पाइपलाइन में दोनों प्रोसेस को अलग-अलग इकाइयों के रूप में परीक्षण करने जा रहे हैं।

### 2.1. `sayHello` प्रोसेस का परीक्षण करो

आइए `sayHello` प्रोसेस से शुरू करें।

आइए प्रोसेस के लिए परीक्षण जनरेट करने के लिए `nf-test generate` कमांड का फिर से उपयोग करें।

```bash
nf-test generate process main.nf
```

```console title="Output"
> nf-test generate process main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.sayhello.nf.test
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

आइए अभी के लिए `main.sayhello.nf.test` फ़ाइल में `sayhello` प्रोसेस पर ध्यान केंद्रित करें।

आइए फ़ाइल खोलें और सामग्री देखें।

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

    test("Should run without failures") {

        when {
            params {
                // यहाँ पैरामीटर परिभाषित करें। उदाहरण:
                // outdir = "tests/results"
            }
            process {
                """
                // यहाँ प्रोसेस के इनपुट परिभाषित करें। उदाहरण:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

पहले की तरह, हम परीक्षण विवरण से शुरू करते हैं, उसके बाद `when` और `then` ब्लॉक। हालाँकि, हमारे पास एक अतिरिक्त `process` ब्लॉक भी है जो हमें प्रोसेस के इनपुट परिभाषित करने की अनुमति देता है।

आइए परीक्षण चलाएं और देखें कि यह काम करता है या नहीं।

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [1eaad118] 'Should run without failures' FAILED (4.876s)

  Assertion failed:

  assert process.success
         |       |
         |       false
         sayHello

  Nextflow stdout:

  Process `sayHello` declares 1 input but was called with 0 arguments
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.884s (1 failed)
```

परीक्षण विफल होता है क्योंकि `sayHello` प्रोसेस 1 इनपुट घोषित करता है लेकिन 0 arguments के साथ बुलाया गया था। आइए प्रोसेस में एक इनपुट जोड़कर इसे ठीक करें। [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (और ऊपर वार्मअप अनुभाग) से याद करो कि हमारा `sayHello` प्रोसेस एक single value इनपुट लेता है, जिसे हमें प्रदान करना होगा। हमें परीक्षण का नाम भी ठीक करना चाहिए ताकि यह बेहतर ढंग से दर्शाए कि हम क्या परीक्षण कर रहे हैं।

=== "बाद में"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10"
        test("Should run without failures and produce correct output") {

            when {
                params {
                    // यहाँ पैरामीटर परिभाषित करें। उदाहरण:
                    // outdir = "tests/results"
                }
                process {
                    """
                    input[0] = "hello"
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

=== "पहले"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
        test("Should run without failures") {

            when {
                params {
                    // यहाँ पैरामीटर परिभाषित करें। उदाहरण:
                    // outdir = "tests/results"
                }
                process {
                    """
                    // यहाँ प्रोसेस के इनपुट परिभाषित करें। उदाहरण:
                    // input[0] = file("test-file.txt")
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

आइए परीक्षण फिर से चलाएं और देखें कि यह काम करता है या नहीं।

```console title="nf-test pipeline pass"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.604s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.611s
```

सफलता! परीक्षण पास होता है क्योंकि `sayHello` प्रोसेस सफलतापूर्वक चला और आउटपुट बनाया गया।

### 2.2. परीक्षण द्वारा बनाए गए snapshot को देखो

अगर हम `tests/main.sayhello.nf.test` फ़ाइल को देखें, तो हम देख सकते हैं कि यह assertion ब्लॉक में `snapshot()` method का उपयोग करती है:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

यह nf-test को `sayHello` प्रोसेस के आउटपुट का snapshot बनाने के लिए कह रहा है। आइए snapshot फ़ाइल की सामग्री देखें।

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

हम इसे यहाँ प्रिंट नहीं करेंगे, लेकिन तुम्हें प्रोसेस और प्रोसेस आउटपुट के विवरण वाली एक JSON फ़ाइल दिखनी चाहिए। विशेष रूप से, हम एक लाइन देख सकते हैं जो इस तरह दिखती है:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

यह `sayHello` प्रोसेस द्वारा बनाए गए आउटपुट का प्रतिनिधित्व करता है, जिसका हम स्पष्ट रूप से परीक्षण कर रहे हैं। अगर हम परीक्षण फिर से चलाते हैं, तो प्रोग्राम जाँच करेगा कि नया आउटपुट मूल रूप से रिकॉर्ड किए गए आउटपुट से मेल खाता है। यह परीक्षण करने का एक त्वरित, सरल तरीका है कि प्रोसेस आउटपुट नहीं बदलते, यही कारण है कि nf-test इसे डिफ़ॉल्ट के रूप में प्रदान करता है।

!!!warning "चेतावनी"

    इसका मतलब है कि हमें यह सुनिश्चित करना होगा कि मूल रन में हम जो आउटपुट रिकॉर्ड करते हैं वह सही है!

अगर, भविष्य के विकास के दौरान, कोड में कुछ बदलता है जिससे आउटपुट अलग हो जाता है, तो परीक्षण विफल हो जाएगा और हमें यह निर्धारित करना होगा कि परिवर्तन अपेक्षित है या नहीं।

- अगर यह पता चलता है कि कोड में कुछ टूट गया, तो हमें इसे ठीक करना होगा, इस अपेक्षा के साथ कि ठीक किया गया कोड परीक्षण पास करेगा।
- अगर यह एक अपेक्षित परिवर्तन है (जैसे, टूल में सुधार हुआ है और परिणाम बेहतर हैं) तो हमें नए आउटपुट को मिलान के लिए संदर्भ के रूप में स्वीकार करने के लिए snapshot को अपडेट करना होगा। nf-test में इस उद्देश्य के लिए `--update-snapshot` पैरामीटर है।

हम परीक्षण फिर से चला सकते हैं और देख सकते हैं कि परीक्षण पास होना चाहिए:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

सफलता! परीक्षण पास होता है क्योंकि `sayHello` प्रोसेस सफलतापूर्वक चला और आउटपुट snapshot से मेल खाया।

### 2.3. Snapshots का विकल्प: Direct Content Assertions

जबकि snapshots आउटपुट में किसी भी बदलाव को पकड़ने के लिए बहुत अच्छे हैं, कभी-कभी तुम पूरी फ़ाइल के मिलान के बारे में इतने सख्त हुए बिना विशिष्ट सामग्री को सत्यापित करना चाहते हो। उदाहरण के लिए:

- जब आउटपुट के कुछ हिस्से बदल सकते हैं (timestamps, random IDs, आदि) लेकिन कुछ मुख्य सामग्री मौजूद होनी चाहिए
- जब तुम आउटपुट में विशिष्ट patterns या values की जाँच करना चाहते हो
- जब तुम परीक्षण को इस बारे में अधिक स्पष्ट बनाना चाहते हो कि सफलता क्या है

यहाँ बताया गया है कि हम विशिष्ट सामग्री की जाँच करने के लिए अपने परीक्षण को कैसे संशोधित कर सकते हैं:

=== "बाद में"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 16 17"
        test("Should run without failures and contain expected greeting") {

            when {
                params {
                    // यहाँ पैरामीटर परिभाषित करें
                }
                process {
                    """
                    input[0] = "hello"
                    """
                }
            }

            then {
                assert process.success
                assert path(process.out[0][0]).readLines().contains('hello')
                assert !path(process.out[0][0]).readLines().contains('HELLO')
            }

        }
    ```

=== "पहले"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
        test("Should run without failures and produce correct output") {

            when {
                params {
                    // यहाँ पैरामीटर परिभाषित करें। उदाहरण:
                    // outdir = "tests/results"
                }
                process {
                    """
                    input[0] = "hello"
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

ध्यान दो कि nf-test प्रोसेस आउटपुट को lists की list के रूप में देखता है, इसलिए `process.out[0][0]` इस प्रोसेस से पहले चैनल आइटम (या 'emission') का पहला भाग प्राप्त कर रहा है।

यह दृष्टिकोण:

- यह स्पष्ट करता है कि हम आउटपुट में क्या अपेक्षा करते हैं
- आउटपुट में अप्रासंगिक परिवर्तनों के प्रति अधिक लचीला है
- परीक्षण विफल होने पर बेहतर त्रुटि संदेश प्रदान करता है
- अधिक जटिल validations (regex patterns, numerical comparisons, आदि) की अनुमति देता है

आइए परीक्षण चलाएं और देखें कि यह काम करता है या नहीं।

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (7.196s)


SUCCESS: Executed 1 tests in 7.208s
```

### 2.4. `convertToUpper` प्रोसेस का परीक्षण करो

आइए `tests/main.converttoupper.nf.test` फ़ाइल खोलें और सामग्री देखें:

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

    test("Should run without failures") {

        when {
            params {
                // यहाँ पैरामीटर परिभाषित करें। उदाहरण:
                // outdir = "tests/results"
            }
            process {
                """
                // यहाँ प्रोसेस के इनपुट परिभाषित करें। उदाहरण:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

यह `sayHello` प्रोसेस के समान परीक्षण है, लेकिन यह `convertToUpper` प्रोसेस का परीक्षण कर रहा है। हम जानते हैं कि यह विफल होगा क्योंकि `sayHello` की तरह, `convertToUpper` प्रोसेस एक single path इनपुट लेता है, लेकिन हमने एक निर्दिष्ट नहीं किया है।

अब हमें convertToUpper प्रोसेस को एक single इनपुट फ़ाइल प्रदान करनी होगी, जिसमें कुछ टेक्स्ट है जिसे हम अपरकेस में बदलना चाहते हैं। हम इसे कई तरीकों से कर सकते हैं:

- हम परीक्षण के लिए एक समर्पित फ़ाइल बना सकते हैं
- हम मौजूदा data/greetings.csv फ़ाइल का पुनः उपयोग कर सकते हैं
- हम इसे परीक्षण के भीतर ही बना सकते हैं

अभी के लिए, आइए पाइपलाइन स्तर के परीक्षण के साथ उपयोग किए गए उदाहरण का उपयोग करके मौजूदा data/greetings.csv फ़ाइल का पुनः उपयोग करें। पहले की तरह, हम परीक्षण का नाम बेहतर ढंग से दर्शाने के लिए रख सकते हैं कि हम क्या परीक्षण कर रहे हैं, लेकिन इस बार आइए इसे सामग्री को 'snapshot' करने दें बजाय विशिष्ट strings की जाँच करने के (जैसा हमने दूसरे प्रोसेस में किया था)।

=== "बाद में"

    ```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10"
        test("Should run without failures and produce correct output") {

            when {
                params {
                    // यहाँ पैरामीटर परिभाषित करें। उदाहरण:
                    // outdir = "tests/results"
                }
                process {
                    """
                    input[0] = "${projectDir}/greetings.csv"
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

=== "पहले"

    ```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
        test("Should run without failures") {

            when {
                params {
                    // यहाँ पैरामीटर परिभाषित करें। उदाहरण:
                    // outdir = "tests/results"
                }
                process {
                    """
                    // यहाँ प्रोसेस के इनपुट परिभाषित करें। उदाहरण:
                    // input[0] = file("test-file.txt")
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

और परीक्षण चलाओ!

```bash title="nf-test pipeline pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.755s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.764s
```

ध्यान दो, हमने `convertToUpper` प्रोसेस के लिए `tests/main.converttoupper.nf.test.snap` पर एक snapshot फ़ाइल बनाई है। अगर हम परीक्षण फिर से चलाते हैं, तो हमें nf-test फिर से पास होते देखना चाहिए।

```bash title="nf-test process convertToUpper pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s
```

### सारांश

तुम जानते हो कि Nextflow प्रोसेस के लिए परीक्षण कैसे लिखें और उन्हें कैसे चलाएं।

### आगे क्या है?

एक साथ सब कुछ के लिए परीक्षण चलाना सीखो!

## 3. पूरे repository के लिए परीक्षण चलाओ

प्रत्येक घटक पर nf-test चलाना ठीक है, लेकिन यह श्रमसाध्य और त्रुटि-प्रवण है। क्या हम एक साथ सब कुछ परीक्षण नहीं कर सकते?

हाँ, हम कर सकते हैं!

आइए पूरे repo पर nf-test चलाएं।

### 3.1. पूरे repo पर nf-test चलाओ

हम `nf-test test` कमांड चलाकर पूरे repo पर nf-test चला सकते हैं।

```bash
nf-test test .
```

ध्यान दो, हम अपनी वर्तमान डायरेक्टरी से सब कुछ चलाने के लिए बस `.` का उपयोग कर रहे हैं। इसमें हर परीक्षण शामिल होगा!

```console title="nf-test repo pass"
> nf-test test .

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [3d26d9af] 'Should run without failures and produce correct output' PASSED (4.155s)

Test Workflow main.nf

  Test [f183df37] 'Should run successfully with correct number of processes' PASSED (3.33s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (3.102s)

Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (2.614s)


SUCCESS: Executed 4 tests in 13.481s
```

देखो! हमने एक ही कमांड से 4 परीक्षण चलाए, प्रत्येक प्रोसेस के लिए 1 और पूरी पाइपलाइन के लिए 2। कल्पना करो कि एक बड़े codebase पर यह कितना शक्तिशाली है!

---

## सारांश

इस साइड क्वेस्ट में, तुमने अलग-अलग प्रोसेस के लिए परीक्षण बनाने और चलाने के साथ-साथ पूरी पाइपलाइन के end-to-end परीक्षण के लिए nf-test की सुविधाओं का लाभ उठाना सीखा है।
तुम अब आउटपुट validation के दो मुख्य दृष्टिकोणों, snapshots और direct content assertions, और कब किसका उपयोग करना है, से परिचित हो।
तुम यह भी जानते हो कि परीक्षण एक-एक करके या पूरे प्रोजेक्ट के लिए कैसे चलाएं।

अपने काम में इन तकनीकों को लागू करने से तुम यह सुनिश्चित कर सकोगे कि:

- तुम्हारा कोड अपेक्षित रूप से काम करता है
- बदलाव मौजूदा कार्यक्षमता को नहीं तोड़ते
- अन्य डेवलपर आत्मविश्वास के साथ योगदान कर सकते हैं
- समस्याओं को जल्दी पहचाना और ठीक किया जा सकता है
- आउटपुट सामग्री अपेक्षाओं से मेल खाती है

### मुख्य patterns

<!-- TODO: Can we add snippets of code below to illustrate? -->

1. पाइपलाइन-स्तरीय परीक्षण:
   - बुनियादी सफलता परीक्षण
   - प्रोसेस गणना सत्यापन
   - आउटपुट फ़ाइल अस्तित्व जाँच
2. प्रोसेस-स्तरीय परीक्षण
3. आउटपुट validation के दो दृष्टिकोण:
   - पूर्ण आउटपुट सत्यापन के लिए snapshots का उपयोग
   - विशिष्ट सामग्री जाँच के लिए direct content assertions का उपयोग
4. एक ही कमांड से repository में सभी परीक्षण चलाना

### अतिरिक्त संसाधन

अधिक उन्नत परीक्षण सुविधाओं और सर्वोत्तम प्रथाओं के लिए [nf-test दस्तावेज़ीकरण](https://www.nf-test.com/) देखो। तुम शायद चाहोगे:

- अपने परीक्षणों में अधिक व्यापक assertions जोड़ना
- edge cases और त्रुटि स्थितियों के लिए परीक्षण लिखना
- परीक्षणों को स्वचालित रूप से चलाने के लिए continuous integration सेट करना
- वर्कफ़्लो और मॉड्यूल परीक्षणों जैसे अन्य प्रकार के परीक्षणों के बारे में जानना
- अधिक उन्नत content validation तकनीकों का पता लगाना

**याद रखो:** परीक्षण इस बात का जीवंत दस्तावेज़ीकरण है कि तुम्हारा कोड कैसे व्यवहार करना चाहिए। जितने अधिक परीक्षण तुम लिखते हो, और जितने अधिक विशिष्ट तुम्हारे assertions होते हैं, उतना ही अधिक तुम अपनी पाइपलाइन की विश्वसनीयता पर भरोसा कर सकते हो।

---

## आगे क्या है?

[साइड क्वेस्ट के मेनू](../) पर वापस जाओ या सूची में अगले विषय पर जाने के लिए पृष्ठ के नीचे दाईं ओर बटन पर क्लिक करो।
