# nf-test के साथ टेस्टिंग

यह सुनिश्चित करने के लिए कि आपके वर्कफ़्लो का हर हिस्सा वही कर रहा है जो उसे करना चाहिए, व्यवस्थित रूप से टेस्ट करना reproducibility और लंबे समय तक maintenance के लिए महत्वपूर्ण है, और development प्रक्रिया के दौरान बहुत मददगार हो सकता है।

आइए एक मिनट के लिए बात करें कि टेस्टिंग इतनी महत्वपूर्ण क्यों है। अगर तुम एक वर्कफ़्लो develop कर रहे हो, तो सबसे पहले तुम कुछ टेस्ट डेटा लोगे जो valid है और जिसका परिणाम मिलना चाहिए। तुम पाइपलाइन में पहला प्रोसेस जोड़ते हो और इसे अपने इनपुट से जोड़कर काम करवाते हो। फिर, यह चेक करने के लिए कि सब काम कर रहा है, तुम इसे टेस्ट डेटा पर चलाते हो। मान लो यह काम करता है, तो तुम अगले प्रोसेस पर जाते हो और टेस्ट डेटा फिर से चलाते हो। तुम यह प्रक्रिया तब तक दोहराते हो जब तक तुम्हें एक पाइपलाइन नहीं मिल जाती जिससे तुम खुश हो।

फिर, शायद तुम एक सरल true या false पैरामीटर जैसे `--skip_process` जोड़ते हो। अब तुम्हें पाइपलाइन दो बार चलानी होगी, हर पैरामीटर के साथ एक बार यह सुनिश्चित करने के लिए कि यह expected तरीके से काम करती है। लेकिन रुको, हम कैसे चेक करें कि `--skip_process` वास्तव में प्रोसेस को skip करता है? हमें आउटपुट में खोदना होगा या लॉग फ़ाइलें चेक करनी होंगी! यह परेशानी भरा है और गलती होने की संभावना है।

जैसे-जैसे तुम अपनी पाइपलाइन develop करते हो, यह जल्दी ही इतनी जटिल हो जाएगी कि हर iteration को manually टेस्ट करना धीमा और गलती-प्रवण होगा। इसके अलावा, अगर तुम्हें कोई error मिलती है तो यह पता लगाना बहुत मुश्किल होगा कि तुम्हारी पाइपलाइन में error कहां से आ रही है। यहीं पर टेस्टिंग काम आती है।

टेस्टिंग तुम्हें व्यवस्थित रूप से चेक करने देती है कि तुम्हारी पाइपलाइन का हर हिस्सा expected तरीके से काम कर रहा है। अच्छी तरह से लिखे गए टेस्ट के developer के लिए फायदे बहुत बड़े हैं:

- **विश्वास**: क्योंकि टेस्ट पूरी पाइपलाइन को कवर करते हैं, तुम confident हो सकते हो कि कुछ बदलने से कुछ और प्रभावित नहीं होता
- **भरोसा**: जब कई developers पाइपलाइन पर काम करते हैं, तो वे जानते हैं कि दूसरे developers ने पाइपलाइन और हर component को नहीं तोड़ा है
- **पारदर्शिता**: टेस्ट दिखाते हैं कि पाइपलाइन कहां fail हो रही है और समस्या को track करना आसान बनाते हैं। वे documentation के एक रूप के रूप में भी काम करते हैं, दिखाते हैं कि प्रोसेस या वर्कफ़्लो कैसे चलाएं
- **गति**: क्योंकि टेस्ट automated हैं, वे बहुत जल्दी और बार-बार चलाए जा सकते हैं। तुम नए bugs introduce करने के डर के बिना जल्दी iterate कर सकते हो

हम कई अलग-अलग प्रकार के टेस्ट लिख सकते हैं:

1. **Module-level tests**: व्यक्तिगत processes के लिए
2. **Workflow-level tests**: एक single वर्कफ़्लो के लिए
3. **Pipeline-level tests**: पूरी पाइपलाइन के लिए
4. **Performance tests**: पाइपलाइन की गति और efficiency के लिए
5. **Stress tests**: extreme conditions के तहत पाइपलाइन के performance का आकलन करना इसकी सीमाएं निर्धारित करने के लिए

व्यक्तिगत processes को टेस्ट करना अन्य भाषाओं में unit tests के समान है। वर्कफ़्लो या पूरी पाइपलाइन को टेस्ट करना अन्य भाषाओं में integration tests कहलाता है, जहां हम components के interactions को टेस्ट करते हैं।

[**nf-test**](https://www.nf-test.com/) एक tool है जो तुम्हें module, workflow और pipeline level test लिखने देता है। संक्षेप में, यह तुम्हें व्यवस्थित रूप से चेक करने देता है कि पाइपलाइन का हर व्यक्तिगत हिस्सा expected तरीके से काम कर रहा है, _isolation में_।

### सीखने के लक्ष्य

इस side quest में, तुम nf-test का उपयोग करना सीखोगे पाइपलाइन के लिए workflow-level test लिखने के साथ-साथ तीन processes के लिए module-level tests लिखने के लिए जो यह call करती है।

इस side quest के अंत तक, तुम निम्नलिखित तकनीकों का प्रभावी ढंग से उपयोग कर सकोगे:

- अपने project में nf-test initialize करना
- Module-level और workflow-level tests generate करना
- सामान्य प्रकार के assertions जोड़ना
- समझना कि snapshots बनाम content assertions का उपयोग कब करें
- पूरे project के लिए tests चलाना

ये skills तुम्हें अपने pipeline projects में एक comprehensive testing strategy implement करने में मदद करेंगी, यह सुनिश्चित करते हुए कि वे अधिक robust और maintainable हैं।

### पूर्वापेक्षाएँ

इस side quest को लेने से पहले, तुम्हें:

- [Hello Nextflow](../hello_nextflow/README.md) tutorial या equivalent beginner's course पूरा करना चाहिए।
- बुनियादी Nextflow concepts और mechanisms (processes, channels, operators, फ़ाइलों के साथ काम करना, meta data) का उपयोग करने में सहज होना चाहिए

---

## 0. शुरू करो

#### Training codespace खोलो

अगर तुमने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित training environment खोलना सुनिश्चित करो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Project directory में जाओ

आइए उस directory में चलें जहां इस tutorial के लिए फ़ाइलें स्थित हैं।

```bash
cd side-quests/nf-test
```

तुम VSCode को इस directory पर focus करने के लिए सेट कर सकते हो:

```bash
code .
```

#### Materials की समीक्षा करो

तुम्हें एक main workflow फ़ाइल और एक CSV फ़ाइल मिलेगी जिसे `greetings.csv` कहा जाता है जिसमें पाइपलाइन का इनपुट है।

```console title="Directory contents"
.
├── greetings.csv
└── main.nf
```

फ़ाइलों के विस्तृत विवरण के लिए, [Hello Nextflow से warmup](../hello_nextflow/00_orientation.md) देखो।

जिस वर्कफ़्लो को हम टेस्ट करेंगे वह [Hello Workflow](../hello_nextflow/03_hello_workflow.md) में बनाए गए Hello वर्कफ़्लो का एक subset है।

??? example "Hello Nextflow वर्कफ़्लो क्या करता है?"

    अगर तुमने [Hello Nextflow](../hello_nextflow/index.md) training नहीं की है, तो यहां एक quick overview है कि यह simple वर्कफ़्लो क्या करता है।

    वर्कफ़्लो एक CSV फ़ाइल लेता है जिसमें greetings हैं, उन पर चार consecutive transformation steps चलाता है, और एक single text फ़ाइल output करता है जिसमें एक fun character की ASCII picture है जो greetings कह रहा है।

    चार steps को Nextflow processes (`sayHello`, `convertToUpper`, `collectGreetings`, और `cowpy`) के रूप में implement किया गया है जो अलग module फ़ाइलों में stored हैं।

    1. **`sayHello`:** हर greeting को अपनी output फ़ाइल में लिखता है (जैसे, "Hello-output.txt")
    2. **`convertToUpper`:** हर greeting को uppercase में convert करता है (जैसे, "HELLO")
    3. **`collectGreetings`:** सभी uppercase greetings को एक single batch फ़ाइल में collect करता है
    4. **`cowpy`:** `cowpy` tool का उपयोग करके ASCII art generate करता है

    परिणाम `results/` नामक directory में publish किए जाते हैं, और पाइपलाइन का final output (जब default parameters के साथ चलाया जाता है) एक plain text फ़ाइल है जिसमें एक character की ASCII art है जो uppercased greetings कह रहा है।

    इस side quest में, हम Hello वर्कफ़्लो के एक intermediate form का उपयोग करते हैं जिसमें केवल पहले दो processes हैं। <!-- TODO: change this to use the full finished workflow as suggested in https://github.com/nextflow-io/training/issues/735 -->

जिस subset के साथ हम काम करेंगे वह दो processes से बना है: `sayHello` और `convertToUpper`।
तुम नीचे पूरा वर्कफ़्लो code देख सकते हो।

??? example "Workflow code"

    ```groovy title="main.nf"
    /*
    * Pipeline parameters
    */
    params.input_file = "greetings.csv"

    /*
    * Use echo to print 'Hello World!' to standard out
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
    * Use a text replace utility to convert the greeting to uppercase
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

        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
    }
    ```

#### वर्कफ़्लो चलाओ

आइए वर्कफ़्लो चलाएं यह सुनिश्चित करने के लिए कि यह expected तरीके से काम कर रहा है।

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

बधाई हो! तुमने अभी एक टेस्ट चलाया!

"रुको, क्या? मैंने अभी वर्कफ़्लो चलाया और यह काम किया! यह टेस्ट कैसे है?"

अच्छा सवाल!

आइए break down करें कि अभी क्या हुआ।

तुमने default parameters के साथ वर्कफ़्लो चलाया, तुमने confirm किया कि यह काम किया और तुम परिणामों से खुश हो। यह टेस्टिंग का सार है। अगर तुमने Hello Nextflow training course किया है, तो तुमने देखा होगा कि हमने हमेशा हर section की शुरुआत उस वर्कफ़्लो को चलाकर की जिसे हम starting point के रूप में उपयोग कर रहे थे, यह confirm करने के लिए कि सब कुछ सही तरीके से सेट है।

Software को टेस्ट करना essentially हमारे लिए यह प्रक्रिया करता है।

#### Assignment की समीक्षा करो

तुम्हारी चुनौती है इस वर्कफ़्लो में nf-test का उपयोग करके standardized tests जोड़ना, ताकि यह verify करना आसान हो जाए कि हर हिस्सा expected तरीके से काम करना जारी रखता है अगर कोई और बदलाव किए जाते हैं।

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### तैयारी checklist

लगता है तुम dive in करने के लिए तैयार हो?

- [ ] मैं इस course के लक्ष्य और इसकी पूर्वापेक्षाओं को समझता हूं
- [ ] मेरा codespace चल रहा है
- [ ] मैंने अपनी working directory उचित रूप से सेट की है
- [ ] मैंने वर्कफ़्लो सफलतापूर्वक चलाया है
- [ ] मैं assignment को समझता हूं

अगर तुम सभी boxes को check कर सकते हो, तो तुम जाने के लिए तैयार हो।

---

## 1. `nf-test` initialize करो

`nf-test` package एक initialization command प्रदान करता है जो हमारे project के लिए tests develop करना शुरू करने के लिए कुछ चीजें सेट करता है।

```bash
nf-test init
```

यह निम्नलिखित output produce करना चाहिए:

```bash
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

यह एक `tests` directory भी बनाता है जिसमें एक configuration फ़ाइल stub है।

### 1.1. एक nf-test generate करो

`nf-test` nf-test फ़ाइलें बनाने के लिए tools के एक set के साथ आता है, जो हमें अधिकांश काम बचाता है। ये subcommand `generate` के तहत आते हैं। आइए पाइपलाइन के लिए एक टेस्ट generate करें:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

यह `tests` directory के भीतर एक `main.nf.test` फ़ाइल बनाएगा। यह हमारी pipeline level test फ़ाइल है। अगर तुम `tree tests/` चलाते हो तो तुम्हें कुछ ऐसा दिखना चाहिए:

```console title="Test directory contents"
tests/
├── main.nf.test
└── nextflow.config
```

`main.nf.test` फ़ाइल हमारी pipeline level test फ़ाइल है। आइए इसे खोलें और contents को देखें।

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

हम test फ़ाइल की structure को समझने के लिए एक second लेंगे।

`nextflow_pipeline` block सभी pipeline level tests के लिए entry point है। इसमें निम्नलिखित शामिल हैं:

- `name`: टेस्ट का नाम।
- `script`: पाइपलाइन script का path।

`test` block वास्तविक टेस्ट है। इसमें निम्नलिखित शामिल हैं:

- `when`: वे conditions जिनके तहत टेस्ट चलाया जाना चाहिए। इसमें वे parameters शामिल हैं जो पाइपलाइन चलाने के लिए उपयोग किए जाएंगे।
- `then`: वे assertions जो किए जाने चाहिए। इसमें पाइपलाइन के expected outcomes शामिल हैं।

सरल हिंदी में, टेस्ट का logic इस प्रकार पढ़ता है:
"**जब** ये _parameters_ इस _pipeline_ को प्रदान किए जाते हैं, **तो** हम इन परिणामों को देखने की उम्मीद करते हैं।"

यह एक functional test नहीं है, हम अगले section में demonstrate करेंगे कि इसे एक में कैसे बदलें।

### Test Names पर एक नोट

ऊपर के उदाहरण में, हमने default नाम "Should run without failures" का उपयोग किया जो एक basic test के लिए उपयुक्त है जो सिर्फ चेक करता है कि पाइपलाइन सफलतापूर्वक चलती है या नहीं। हालांकि, जैसे-जैसे हम अधिक specific test cases जोड़ते हैं, हमें अधिक descriptive नाम उपयोग करने चाहिए जो indicate करें कि हम वास्तव में क्या टेस्ट कर रहे हैं। उदाहरण के लिए:

- "Should convert input to uppercase" - जब specific functionality को टेस्ट कर रहे हों
- "Should handle empty input gracefully" - जब edge cases को टेस्ट कर रहे हों
- "Should respect max memory parameter" - जब resource constraints को टेस्ट कर रहे हों
- "Should create expected output files" - जब फ़ाइल generation को टेस्ट कर रहे हों

अच्छे test names को चाहिए:

1. "Should" से शुरू होना ताकि यह clear हो कि expected behavior क्या है
2. Specific functionality या scenario का वर्णन करना जो टेस्ट किया जा रहा है
3. इतना clear होना कि अगर टेस्ट fail होता है, तो तुम जानते हो कि कौन सी functionality broken है

जैसे-जैसे हम बाद में अधिक assertions और specific test cases जोड़ते हैं, हम इन अधिक descriptive नामों का उपयोग करेंगे यह clear करने के लिए कि हर टेस्ट क्या verify कर रहा है।

### 1.2. टेस्ट चलाओ

आइए टेस्ट चलाएं यह देखने के लिए कि क्या होता है।

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

टेस्ट fail होता है! क्या हुआ?

1. nf-test ने पाइपलाइन को जैसा है वैसे चलाने की कोशिश की, `when` block में settings का उपयोग करके:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test ने पाइपलाइन की status चेक की और इसे `when` block से compare किया:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

ध्यान दो कि nf-test ने कैसे report किया है कि पाइपलाइन fail हुई और Nextflow से error message प्रदान किया:

```console title="Error"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

तो issue क्या था? याद रखो पाइपलाइन में project directory में एक greetings.csv फ़ाइल है। जब nf-test पाइपलाइन चलाता है, तो यह इस फ़ाइल को ढूंढेगा, लेकिन यह इसे नहीं ढूंढ पाता। फ़ाइल वहां है, क्या हो रहा है? खैर, अगर हम path को देखें तो हम देख सकते हैं कि टेस्ट path `./nf-test/tests/longHashString/` में हो रहा है। Nextflow की तरह, nf-test हर टेस्ट के लिए एक नई directory बनाता है सब कुछ isolated रखने के लिए। डेटा फ़ाइल वहां स्थित नहीं है इसलिए हमें original test में फ़ाइल के path को सही करना होगा।

आइए test फ़ाइल पर वापस जाएं और `when` block में फ़ाइल के path को बदलें।

तुम सोच रहे होगे कि हम टेस्ट में पाइपलाइन के root को कैसे point करने जा रहे हैं। चूंकि यह एक common situation है, nf-test के पास global variables की एक range है जिसका हम अपने जीवन को आसान बनाने के लिए उपयोग कर सकते हैं। तुम पूरी list [यहां](https://www.nf-test.com/docs/testcases/global_variables/) पा सकते हो लेकिन इस बीच हम `projectDir` variable का उपयोग करेंगे, जिसका मतलब है पाइपलाइन project का root।

_पहले:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

_बाद में:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
when {
    params {
        input_file = "${projectDir}/greetings.csv"
    }
}
```

आइए टेस्ट फिर से चलाएं यह देखने के लिए कि यह काम करता है या नहीं।

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

सफलता! पाइपलाइन सफलतापूर्वक चलती है और टेस्ट pass होता है। इसे जितनी बार चाहो चलाओ और तुम्हें हमेशा वही परिणाम मिलेगा!

Default रूप से, Nextflow output छिपा होता है, लेकिन खुद को convince करने के लिए कि nf-test निश्चित रूप से वर्कफ़्लो चला रहा है, तुम `--verbose` flag का उपयोग कर सकते हो:

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

एक simple check यह सुनिश्चित करना है कि हमारी पाइपलाइन सभी processes चला रही है जिनकी हम उम्मीद करते हैं और किसी को चुपचाप skip नहीं कर रही है। याद रखो हमारी पाइपलाइन 6 processes चलाती है, एक `sayHello` कहलाती है और एक `convertToUpper` कहलाती है 3 greetings में से प्रत्येक के लिए।

आइए अपने टेस्ट में एक assertion जोड़ें यह चेक करने के लिए कि पाइपलाइन expected संख्या में processes चलाती है। हम अपने test name को भी update करेंगे ताकि यह बेहतर reflect करे कि हम क्या टेस्ट कर रहे हैं।

**पहले:**

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

**बाद में:**

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

Test name अब बेहतर reflect करता है कि हम वास्तव में क्या verify कर रहे हैं - सिर्फ यह नहीं कि पाइपलाइन fail हुए बिना चलती है, बल्कि यह कि यह expected संख्या में processes चलाती है।

आइए टेस्ट फिर से चलाएं यह देखने के लिए कि यह काम करता है या नहीं।

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

सफलता! पाइपलाइन सफलतापूर्वक चलती है और टेस्ट pass होता है। अब हमने पाइपलाइन के details को टेस्ट करना शुरू कर दिया है, साथ ही overall status को भी।

### 1.4. Output को टेस्ट करो

आइए अपने टेस्ट में एक assertion जोड़ें यह चेक करने के लिए कि output फ़ाइल बनाई गई थी। हम इसे एक अलग test के रूप में जोड़ेंगे, एक informative name के साथ, ताकि परिणामों को interpret करना आसान हो।

**पहले:**

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

**बाद में:**

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

टेस्ट फिर से चलाओ यह देखने के लिए कि यह काम करता है या नहीं।

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

सफलता! Tests pass होते हैं क्योंकि पाइपलाइन सफलतापूर्वक पूरी हुई, सही संख्या में processes चलीं और output फ़ाइलें बनाई गईं। यह तुम्हें यह भी दिखाना चाहिए कि अपने tests के लिए वे informative names प्रदान करना कितना उपयोगी है।

यह सिर्फ surface है, हम पाइपलाइन के details को चेक करने के लिए assertions लिखते रह सकते हैं, लेकिन अभी के लिए आइए पाइपलाइन के internals को टेस्ट करने की ओर बढ़ें।

### सारांश

तुम जानते हो कि पाइपलाइन के लिए nf-test कैसे लिखें।

### आगे क्या है?

Nextflow process को टेस्ट करना सीखो।

---

## 2. Nextflow process को टेस्ट करो

हमें पाइपलाइन के हर हिस्से के लिए tests लिखने की जरूरत नहीं है, लेकिन जितने अधिक tests हमारे पास होंगे उतना ही हम पाइपलाइन के बारे में comprehensive हो सकते हैं और उतना ही confident हो सकते हैं कि यह expected तरीके से काम कर रही है। इस section में हम पाइपलाइन में दोनों processes को individual units के रूप में टेस्ट करने जा रहे हैं।

### 2.1. `sayHello` process को टेस्ट करो

आइए `sayHello` process से शुरू करें।

आइए process के लिए tests generate करने के लिए फिर से `nf-test generate` command का उपयोग करें।

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

आइए अभी के लिए `main.sayhello.nf.test` फ़ाइल में `sayhello` process पर focus करें।

आइए फ़ाइल खोलें और contents को देखें।

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
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

पहले की तरह, हम test details से शुरू करते हैं, इसके बाद `when` और `then` blocks आते हैं। हालांकि, हमारे पास एक अतिरिक्त `process` block भी है जो हमें process के inputs को define करने देता है।

आइए टेस्ट चलाएं यह देखने के लिए कि यह काम करता है या नहीं।

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

टेस्ट fail होता है क्योंकि `sayHello` process 1 input declare करता है लेकिन 0 arguments के साथ call किया गया था। आइए process में एक input जोड़कर इसे ठीक करें। [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (और ऊपर warmup section) से याद रखो कि हमारा `sayHello` process एक single value input लेता है, जो हमें प्रदान करना होगा। हमें test name को भी ठीक करना चाहिए ताकि यह बेहतर reflect करे कि हम क्या टेस्ट कर रहे हैं।

**पहले:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
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

**बाद में:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
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

आइए टेस्ट फिर से चलाएं यह देखने के लिए कि यह काम करता है या नहीं।

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

सफलता! टेस्ट pass होता है क्योंकि `sayHello` process सफलतापूर्वक चला और output बनाया गया।

### 2.2. Test द्वारा बनाए गए snapshot को देखो

अगर हम `tests/main.sayhello.nf.test` फ़ाइल को देखें, तो हम देख सकते हैं कि यह assertion block में एक method `snapshot()` का उपयोग करता है:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

यह nf-test को `sayHello` process के output का एक snapshot बनाने के लिए कह रहा है। आइए snapshot फ़ाइल के contents को देखें।

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

हम इसे यहां print नहीं करेंगे, लेकिन तुम्हें एक JSON फ़ाइल दिखनी चाहिए जिसमें process और process outputs के details हैं। विशेष रूप से, हम एक line देख सकते हैं जो इस तरह दिखती है:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

यह `sayHello` process द्वारा बनाए गए outputs को represent करता है, जिसे हम explicitly टेस्ट कर रहे हैं। अगर हम टेस्ट फिर से चलाते हैं, तो program चेक करेगा कि नया output उस output से match करता है जो originally record किया गया था। यह process outputs को चेक करने का एक quick, simple तरीका है, यही कारण है कि nf-test इसे default के रूप में प्रदान करता है।

!!!warning

    इसका मतलब है कि हमें यह सुनिश्चित करना होगा कि जो output हम original run में record करते हैं वह सही है!

अगर, भविष्य के development के दौरान, code में कुछ बदलता है जो output को अलग बनाता है, तो टेस्ट fail होगा और हमें यह निर्धारित करना होगा कि परिवर्तन expected है या नहीं।

- अगर यह पता चलता है कि code में कुछ टूट गया, तो हमें इसे ठीक करना होगा, इस expectation के साथ कि fixed code टेस्ट pass करेगा।
- अगर यह एक expected change है (जैसे, tool को improve किया गया है और परिणाम बेहतर हैं) तो हमें snapshot को update करना होगा ताकि नए output को reference के रूप में स्वीकार किया जा सके। nf-test के पास इस उद्देश्य के लिए एक parameter `--update-snapshot` है।

हम टेस्ट फिर से चला सकते हैं और देख सकते हैं कि टेस्ट pass होना चाहिए:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

सफलता! टेस्ट pass होता है क्योंकि `sayHello` process सफलतापूर्वक चला और output snapshot से match किया।

### 2.3. Snapshots का विकल्प: Direct Content Assertions

जबकि snapshots output में किसी भी बदलाव को पकड़ने के लिए बहुत अच्छे हैं, कभी-कभी तुम पूरी फ़ाइल के matching के बारे में इतने strict हुए बिना specific content को verify करना चाहते हो। उदाहरण के लिए:

- जब output के कुछ हिस्से बदल सकते हैं (timestamps, random IDs, आदि) लेकिन कुछ key content मौजूद होना चाहिए
- जब तुम output में specific patterns या values को चेक करना चाहते हो
- जब तुम टेस्ट को इस बारे में अधिक explicit बनाना चाहते हो कि सफलता क्या है

यहां बताया गया है कि हम specific content को चेक करने के लिए अपने टेस्ट को कैसे modify कर सकते हैं:

**पहले:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
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

**बाद में:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 16 17"
    test("Should run without failures and contain expected greeting") {

        when {
            params {
                // define parameters here
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

ध्यान दो कि nf-test process outputs को lists की lists के रूप में देखता है, इसलिए `process.out[0][0]` इस process से पहले channel item (या 'emission') के पहले हिस्से को fetch कर रहा है।

यह approach:

- यह clear करता है कि हम output में exactly क्या expect करते हैं
- Output में irrelevant changes के लिए अधिक resilient है
- जब tests fail होते हैं तो बेहतर error messages प्रदान करता है
- अधिक complex validations की अनुमति देता है (regex patterns, numerical comparisons, आदि)

आइए टेस्ट चलाएं यह देखने के लिए कि यह काम करता है या नहीं।

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

### 2.4. `convertToUpper` process को टेस्ट करो

आइए `tests/main.converttoupper.nf.test` फ़ाइल खोलें और contents को देखें:

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
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

यह `sayHello` process के समान एक टेस्ट है, लेकिन यह `convertToUpper` process को टेस्ट कर रहा है। हम जानते हैं कि यह fail होगा क्योंकि `sayHello` की तरह, `convertToUpper` process एक single path input लेता है, लेकिन हमने एक specify नहीं किया है।

अब हमें convertToUpper process को एक single input फ़ाइल supply करनी होगी, जिसमें कुछ text शामिल है जिसे हम uppercase में convert करना चाहते हैं। हम यह कई तरीकों से कर सकते हैं:

- हम टेस्ट करने के लिए एक dedicated फ़ाइल बना सकते हैं
- हम existing data/greetings.csv फ़ाइल को फिर से उपयोग कर सकते हैं
- हम इसे टेस्ट के भीतर on the fly बना सकते हैं

अभी के लिए, आइए pipeline level test के साथ उपयोग किए गए उदाहरण का उपयोग करके existing data/greetings.csv फ़ाइल को फिर से उपयोग करें। पहले की तरह, हम test को बेहतर reflect करने के लिए name दे सकते हैं कि हम क्या टेस्ट कर रहे हैं, लेकिन इस बार आइए इसे content को 'snapshot' करने दें बजाय specific strings को चेक करने के (जैसा कि हमने दूसरे process में किया)।

**पहले:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
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

**बाद में:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
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

और टेस्ट चलाओ!

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

ध्यान दो, हमने `tests/main.converttoupper.nf.test.snap` पर `convertToUpper` process के लिए एक snapshot फ़ाइल बनाई है। अगर हम टेस्ट फिर से चलाते हैं, तो हमें nf-test फिर से pass होता दिखना चाहिए।

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

तुम जानते हो कि Nextflow process के लिए tests कैसे लिखें और उन्हें कैसे चलाएं।

### आगे क्या है?

एक साथ सब कुछ के लिए tests चलाना सीखो!

## 3. पूरे repository के लिए tests चलाओ

हर component पर nf-test चलाना ठीक है, लेकिन laborious और error prone है। क्या हम एक साथ सब कुछ टेस्ट नहीं कर सकते?

हां हम कर सकते हैं!

आइए पूरे repo पर nf-test चलाएं।

### 3.1. पूरे repo पर nf-test चलाओ

हम `nf-test test` command चलाकर पूरे repo पर nf-test चला सकते हैं।

```bash
nf-test test .
```

ध्यान दो, हम अपनी current directory से सब कुछ चलाने के लिए सिर्फ `.` का उपयोग कर रहे हैं। इसमें हर टेस्ट शामिल होगा!

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

इसे देखो! हमने 4 tests चलाए, हर process के लिए 1 और पूरी पाइपलाइन के लिए 2 एक single command के साथ। कल्पना करो कि यह एक बड़े codebase पर कितना powerful है!

---

## सारांश

इस side quest में, तुमने nf-test की features का leverage करना सीखा है individual processes के लिए tests बनाने और चलाने के साथ-साथ पूरी पाइपलाइन के लिए end-to-end tests बनाने के लिए।
तुम अब output validation के मुख्य दो approaches, snapshots और direct content assertions, और किसी एक का उपयोग कब करना है, के बारे में aware हो।
तुम यह भी जानते हो कि tests को एक-एक करके या पूरे project के लिए कैसे चलाएं।

अपने काम में इन तकनीकों को apply करने से तुम यह सुनिश्चित कर सकोगे कि:

- तुम्हारा code expected तरीके से काम करता है
- Changes existing functionality को नहीं तोड़ते
- दूसरे developers confidence के साथ contribute कर सकते हैं
- समस्याओं को जल्दी identify और fix किया जा सकता है
- Output content expectations से match करता है

### मुख्य patterns

<!-- TODO: Can we add snippets of code below to illustrate? -->

1. Pipeline-level tests:
   - Basic success testing
   - Process count verification
   - Output फ़ाइल existence checks
2. Process-level tests
3. Output validation के लिए दो approaches:
   - Complete output verification के लिए snapshots का उपयोग
   - Specific content checks के लिए direct content assertions का उपयोग
4. एक single command के साथ repository में सभी tests चलाना

### अतिरिक्त संसाधन

अधिक advanced testing features और best practices के लिए [nf-test documentation](https://www.nf-test.com/) देखो। तुम शायद चाहोगे:

- अपने tests में अधिक comprehensive assertions जोड़ना
- Edge cases और error conditions के लिए tests लिखना
- Tests को automatically चलाने के लिए continuous integration सेट करना
- Workflow और module tests जैसे अन्य प्रकार के tests के बारे में सीखना
- अधिक advanced content validation techniques explore करना

**याद रखो:** Tests तुम्हारे code के व्यवहार का living documentation हैं। जितने अधिक tests तुम लिखते हो, और तुम्हारे assertions जितने अधिक specific होते हैं, तुम अपनी पाइपलाइन की reliability में उतने ही confident हो सकते हो।

---

## आगे क्या है?

[Side Quests के menu](./index.md) पर वापस जाओ या list में अगले topic पर जाने के लिए page के नीचे दाईं ओर button पर click करो।
