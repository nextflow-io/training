# nf-test के साथ परीक्षण

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

यह सुनिश्चित करने के लिए कि आपके workflow का हर भाग वही कर रहा है जो उसे करना चाहिए, व्यवस्थित रूप से परीक्षण करना reproducibility और दीर्घकालिक रखरखाव के लिए महत्वपूर्ण है, और development process के दौरान एक बड़ी मदद हो सकती है।

आइए एक मिनट के लिए बात करें कि परीक्षण इतना महत्वपूर्ण क्यों है। यदि आप एक workflow develop कर रहे हैं, तो आप सबसे पहले कुछ test data लेंगे जो आपको पता है कि valid है और result produce करना चाहिए। आप pipeline में पहली process जोड़ते हैं और इसे काम करने के लिए अपने inputs के साथ wire करते हैं। फिर, यह जांचने के लिए कि सब कुछ काम कर रहा है, आप इसे test data पर run करते हैं। यह मानते हुए कि यह काम करता है, आप अगली process पर जाते हैं और फिर से test data run करते हैं। आप इस process को तब तक दोहराते हैं जब तक आपके पास एक pipeline न हो जिससे आप संतुष्ट हों।

फिर, शायद आप एक सरल true या false parameter जैसे `--skip_process` जोड़ते हैं। अब आपको pipeline को दो बार run करना होगा, एक बार प्रत्येक parameter के साथ यह सुनिश्चित करने के लिए कि यह अपेक्षा के अनुसार काम करता है। लेकिन रुकिए, हम यह कैसे जांचें कि `--skip_process` वास्तव में process को skip करता है? हमें outputs में जाना होगा या log files check करनी होंगी! यह एक परेशानी है और error-prone है।

जैसे-जैसे आप अपनी pipeline develop करते हैं, यह जल्दी ही इतनी complex हो जाएगी कि manually हर iteration को test करना slow और error-prone होगा। इसके अलावा, यदि आपको error मिलती है तो यह पता लगाना बहुत मुश्किल होगा कि आपकी pipeline में error कहां से आ रही है। यहीं पर testing काम आती है।

Testing आपको systematically check करने की अनुमति देती है कि आपकी pipeline का हर भाग अपेक्षा के अनुसार काम कर रहा है। developer के लिए अच्छी तरह से लिखे गए tests के फायदे बहुत बड़े हैं:

- **Confidence**: क्योंकि tests पूरी pipeline को cover करते हैं, आप confident रह सकते हैं कि कुछ बदलने से कुछ और प्रभावित नहीं होता है
- **Trust**: जब कई developers pipeline पर काम करते हैं, तो वे जानते हैं कि दूसरे developers ने pipeline और हर component को break नहीं किया है
- **Transparency**: Tests दिखाते हैं कि pipeline कहां fail हो रही है और problem को track करना आसान बनाते हैं। वे documentation के एक रूप के रूप में भी काम करते हैं, दिखाते हैं कि process या workflow कैसे run करें
- **Speed**: क्योंकि tests automated हैं, वे बहुत जल्दी और बार-बार run किए जा सकते हैं। आप नए bugs introduce करने के कम डर के साथ जल्दी iterate कर सकते हैं

हम कई अलग-अलग प्रकार के tests लिख सकते हैं:

1. **Module-level tests**: व्यक्तिगत processes के लिए
2. **Workflow-level tests**: एकल workflow के लिए
3. **Pipeline-level tests**: पूरी pipeline के लिए
4. **Performance tests**: Pipeline की speed और efficiency के लिए
5. **Stress tests**: extreme conditions में pipeline के performance का आकलन उसकी सीमाओं को निर्धारित करने के लिए

व्यक्तिगत processes का परीक्षण अन्य languages में unit tests के समान है। Workflow या पूरी pipeline का परीक्षण अन्य languages में integration tests कहलाता है, जहां हम components के interactions का परीक्षण करते हैं।

[**nf-test**](https://www.nf-test.com/) एक tool है जो आपको module, workflow और pipeline level tests लिखने की अनुमति देता है। संक्षेप में, यह आपको systematically check करने की अनुमति देता है कि pipeline का हर व्यक्तिगत भाग अपेक्षा के अनुसार काम कर रहा है, _isolation में_।

### सीखने के लक्ष्य

इस side quest में, आप pipeline के लिए workflow-level test के साथ-साथ तीन processes के लिए module-level tests लिखने के लिए nf-test का उपयोग करना सीखेंगे।

इस side quest के अंत तक, आप निम्नलिखित techniques को effectively उपयोग करने में सक्षम होंगे:

- अपने project में nf-test initialize करना
- Module-level और workflow-level tests generate करना
- Assertions के सामान्य प्रकार जोड़ना
- समझना कि snapshots vs. content assertions कब उपयोग करें
- पूरे project के लिए tests run करना

ये skills आपको अपने pipeline projects में एक comprehensive testing strategy implement करने में मदद करेंगी, यह सुनिश्चित करते हुए कि वे अधिक robust और maintainable हैं।

### पूर्वापेक्षाएँ

इस side quest को लेने से पहले, आपको:

- [Hello Nextflow](../hello_nextflow/README.md) tutorial या equivalent beginner's course पूरा किया हो
- Basic Nextflow concepts और mechanisms (processes, channels, operators, files के साथ काम करना, meta data) के साथ comfortable होना चाहिए

---

## 0. शुरू करें

#### Training codespace खोलें

यदि आपने अभी तक ऐसा नहीं किया है, तो [Environment Setup](../envsetup/index.md) में वर्णित अनुसार training environment खोलना सुनिश्चित करें।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Project directory में जाएं

आइए उस directory में चलें जहां इस tutorial के लिए files स्थित हैं।

```bash
cd side-quests/nf-test
```

आप VSCode को इस directory पर focus करने के लिए set कर सकते हैं:

```bash
code .
```

#### Materials की समीक्षा करें

आपको एक main workflow file और `greetings.csv` नामक एक CSV file मिलेगी जिसमें pipeline का input होगा।

```console title="Directory contents"
.
├── greetings.csv
└── main.nf
```

Files के विस्तृत description के लिए, [Hello Nextflow से warmup](../hello_nextflow/00_orientation.md) देखें।

जिस workflow का हम परीक्षण करेंगे वह [Hello Workflow](../hello_nextflow/03_hello_workflow.md) में बनाए गए Hello workflow का एक subset है।

??? example "Hello Nextflow workflow क्या करता है?"

    यदि आपने [Hello Nextflow](../hello_nextflow/index.md) training नहीं की है, तो यहां एक quick overview है कि यह simple workflow क्या करता है।

    Workflow एक CSV file लेता है जिसमें greetings होते हैं, उन पर चार consecutive transformation steps run करता है, और एक single text file output करता है जिसमें एक fun character की ASCII picture होती है जो greetings कह रहा है।

    चार steps को Nextflow processes (`sayHello`, `convertToUpper`, `collectGreetings`, और `cowpy`) के रूप में implement किया गया है जो अलग-अलग module files में stored हैं।

    1. **`sayHello`:** प्रत्येक greeting को अपनी output file में लिखता है (जैसे, "Hello-output.txt")
    2. **`convertToUpper`:** प्रत्येक greeting को uppercase में convert करता है (जैसे, "HELLO")
    3. **`collectGreetings`:** सभी uppercase greetings को एक single batch file में collect करता है
    4. **`cowpy`:** `cowpy` tool का उपयोग करके ASCII art generate करता है

    Results `results/` नामक directory में publish किए जाते हैं, और pipeline का final output (जब default parameters के साथ run किया जाता है) एक plain text file होती है जिसमें एक character का ASCII art होता है जो uppercased greetings कह रहा है।

    इस side quest में, हम Hello workflow के एक intermediate form का उपयोग करते हैं जिसमें केवल पहली दो processes होती हैं।

जिस subset के साथ हम काम करेंगे वह दो processes से बना है: `sayHello` और `convertToUpper`।
आप नीचे पूरा workflow code देख सकते हैं।

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

        // CSV फ़ाइल से इनपुट के लिए एक channel बनाएं
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
    }
    ```

#### Workflow चलाएं

आइए workflow को run करें यह सुनिश्चित करने के लिए कि यह अपेक्षा के अनुसार काम कर रहा है।

```bash
nextflow run main.nf
```

```console title="Workflow run करने का परिणाम"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

बधाई हो! आपने अभी-अभी एक test run किया!

"रुकिए, क्या? मैंने अभी workflow run किया और यह काम किया! यह test कैसे है?"

अच्छा सवाल!

आइए समझें कि अभी क्या हुआ।

आपने default parameters के साथ workflow run किया, आपने confirm किया कि यह काम किया और आप results से खुश हैं। यह testing का सार है। यदि आपने Hello Nextflow training course में काम किया है, तो आपने देखा होगा कि हमने हमेशा हर section की शुरुआत उस workflow को run करके की जिसे हम शुरुआती बिंदु के रूप में उपयोग कर रहे थे, यह confirm करने के लिए कि सब कुछ सही तरीके से set up है।

Software का परीक्षण अनिवार्य रूप से हमारे लिए यह process करता है।

#### Assignment की समीक्षा करें

आपकी चुनौती इस workflow में nf-test का उपयोग करके standardized tests जोड़ना है, ताकि यह verify करना आसान हो जाए कि हर भाग अपेक्षा के अनुसार काम करना जारी रखता है यदि कोई और changes किए जाते हैं।

#### Readiness checklist

क्या आपको लगता है कि आप dive in करने के लिए तैयार हैं?

- [ ] मैं इस course के goal और इसकी prerequisites को समझता हूं
- [ ] मेरा codespace up और running है
- [ ] मैंने अपनी working directory को appropriately set कर लिया है
- [ ] मैंने workflow को successfully run किया है
- [ ] मैं assignment को समझता हूं

यदि आप सभी boxes check कर सकते हैं, तो आप जाने के लिए तैयार हैं।

---

## 1. `nf-test` Initialize करें

`nf-test` package एक initialization command provide करता है जो हमारे लिए अपने project के लिए tests develop करना शुरू करने के लिए कुछ चीजें set up करता है।

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

यह एक configuration file stub युक्त `tests` directory भी बनाता है।

### 1.1. nf-test generate करें

`nf-test` nf-test files बनाने के लिए tools के एक set के साथ आता है, जो हमें अधिकांश काम बचाता है। ये `generate` subcommand के अंतर्गत आते हैं। आइए pipeline के लिए एक test generate करें:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

यह `tests` directory के भीतर एक `main.nf.test` file बनाएगा। यह हमारी pipeline level test file है। यदि आप `tree tests/` run करते हैं तो आपको कुछ ऐसा दिखना चाहिए:

```console title="Test directory contents"
tests/
├── main.nf.test
└── nextflow.config
```

`main.nf.test` file हमारी pipeline level test file है। आइए इसे खोलें और contents को देखें।

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

    test("Should run without failures") {

        when {
            params {
                // यहाँ parameters define करें। उदाहरण:
                // outdir = "tests/results" (output directory)
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

हम test file की structure को समझने के लिए एक second लेंगे।

`nextflow_pipeline` block सभी pipeline level tests के लिए entry point है। इसमें निम्नलिखित होते हैं:

- `name`: Test का नाम
- `script`: Pipeline script का path

`test` block actual test है। इसमें निम्नलिखित होते हैं:

- `when`: वे conditions जिनके तहत test run होना चाहिए। इसमें वे parameters शामिल हैं जो pipeline को run करने के लिए उपयोग किए जाएंगे
- `then`: वे assertions जो बनाए जाने चाहिए। इसमें pipeline के expected outcomes शामिल हैं

Plain English में, test का logic निम्नानुसार पढ़ा जाता है:
"**जब** ये _parameters_ इस _pipeline_ को provide किए जाते हैं, **तो** हम इन results को देखने की उम्मीद करते हैं।"

यह एक functional test नहीं है, हम अगले section में demonstrate करेंगे कि इसे कैसे बनाया जाए।

### Test Names पर एक नोट

ऊपर दिए गए उदाहरण में, हमने default name "Should run without failures" का उपयोग किया जो एक basic test के लिए उपयुक्त है जो सिर्फ यह check करता है कि pipeline successfully runs होता है। हालांकि, जैसे ही हम अधिक specific test cases जोड़ते हैं, हमें अधिक descriptive names का उपयोग करना चाहिए जो indicate करते हैं कि हम वास्तव में क्या test कर रहे हैं। उदाहरण के लिए:

- "Should convert input to uppercase" - जब specific functionality test कर रहे हों
- "Should handle empty input gracefully" - जब edge cases test कर रहे हों
- "Should respect max memory parameter" - जब resource constraints test कर रहे हों
- "Should create expected output files" - जब file generation test कर रहे हों

अच्छे test names को:

1. "Should" से शुरू होना चाहिए यह clear करने के लिए कि expected behavior क्या है
2. उस specific functionality या scenario का describe करना चाहिए जो test किया जा रहा है
3. इतना clear होना चाहिए कि यदि test fail होता है, तो आपको पता चल जाए कि कौन सी functionality broken है

जैसे-जैसे हम बाद में अधिक assertions और specific test cases जोड़ते हैं, हम इन अधिक descriptive names का उपयोग करेंगे यह clear करने के लिए कि प्रत्येक test क्या verify कर रहा है।

### 1.2. Test run करें

आइए test run करें यह देखने के लिए कि क्या होता है।

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

Test fail हो जाता है! क्या हुआ?

1. nf-test ने pipeline को as is run करने की कोशिश की, `when` block में settings का उपयोग करते हुए:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test ने pipeline की status check की और इसे `when` block से compare किया:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

ध्यान दें कि कैसे nf-test ने report किया कि pipeline fail हो गई और Nextflow से error message provide किया:

```console title="Error"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

तो issue क्या था? याद रखें कि pipeline में project directory में एक greetings.csv file है। जब nf-test pipeline run करता है, तो यह इस file को look करेगा, लेकिन यह इसे नहीं ढूंढ सकता। File वहां है, क्या हो रहा है? खैर, यदि हम path को देखें तो हम देख सकते हैं कि test path `./nf-test/tests/longHashString/` में हो रहा है। Nextflow की तरह, nf-test हर test के लिए एक नई directory बनाता है सब कुछ isolated रखने के लिए। Data file वहां स्थित नहीं है इसलिए हमें original test में file का path correct करना होगा।

आइए test file पर वापस जाएं और `when` block में file का path change करें।

आप सोच सकते हैं कि हम test में pipeline के root को कैसे point करने जा रहे हैं। चूंकि यह एक common situation है, nf-test के पास global variables की एक range है जिसे हम अपने जीवन को आसान बनाने के लिए उपयोग कर सकते हैं। आप पूरी list [यहां](https://www.nf-test.com/docs/testcases/global_variables/) पा सकते हैं लेकिन इस बीच हम `projectDir` variable का उपयोग करेंगे, जिसका अर्थ है pipeline project का root।

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

आइए test को फिर से run करें यह देखने के लिए कि यह काम करता है या नहीं।

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

सफलता! Pipeline successfully runs होती है और test pass हो जाता है। इसे जितनी बार चाहें run करें और आपको हमेशा same result मिलेगा!

Default रूप से, Nextflow output hidden है, लेकिन खुद को convince करने के लिए कि nf-test definitely workflow run कर रहा है, आप `--verbose` flag का उपयोग कर सकते हैं:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Pipeline सभी processes चलाती है"
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

### 1.3. Assertions जोड़ें

एक simple check यह सुनिश्चित करना है कि हमारी pipeline उन सभी processes को run कर रही है जिनकी हम उम्मीद करते हैं और कोई भी silently skip नहीं कर रही है। याद रखें हमारी pipeline 6 processes run करती है, एक `sayHello` और एक `convertToUpper` 3 greetings में से प्रत्येक के लिए।

आइए अपने test में एक assertion जोड़ें यह check करने के लिए कि pipeline expected number of processes run करती है। हम अपने test name को भी update करेंगे ताकि यह better reflect करे कि हम क्या test कर रहे हैं।

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

Test name अब better reflect करता है कि हम वास्तव में क्या verify कर रहे हैं - न केवल यह कि pipeline failing के बिना runs होती है, बल्कि यह कि यह expected number of processes run करती है।

आइए test को फिर से run करें यह देखने के लिए कि यह काम करता है या नहीं।

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline assertions के साथ passes होती है"
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

सफलता! Pipeline successfully runs होती है और test pass हो जाता है। अब हमने pipeline के details का परीक्षण करना शुरू कर दिया है, साथ ही overall status का भी।

### 1.4. Output का परीक्षण करें

आइए अपने test में एक assertion जोड़ें यह check करने के लिए कि output file बनाई गई थी। हम इसे एक अलग test के रूप में जोड़ेंगे, एक informative name के साथ, ताकि results को interpret करना आसान हो जाए।

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

Test को फिर से run करें यह देखने के लिए कि यह काम करता है या नहीं।

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline file assertions के साथ passes होती है"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

सफलता! Tests pass होते हैं क्योंकि pipeline successfully complete हुई, correct number of processes run हुईं और output files बनाई गईं। यह आपको यह भी दिखाना चाहिए कि अपने tests के लिए वे informative names provide करना कितना useful है।

यह सिर्फ surface है, हम pipeline के details को check करने के लिए assertions लिखना जारी रख सकते हैं, लेकिन अभी के लिए आइए pipeline के internals को test करने की ओर बढ़ें।

### निष्कर्ष

आप जानते हैं कि pipeline के लिए nf-test कैसे लिखें।

### आगे क्या है?

Nextflow process को test करना सीखें।

---

## 2. Nextflow process का परीक्षण करें

हमें pipeline के हर भाग के लिए tests लिखने की जरूरत नहीं है, लेकिन हमारे पास जितने अधिक tests हैं हम pipeline के बारे में उतने ही comprehensive हो सकते हैं और उतने ही confident हो सकते हैं कि यह अपेक्षा के अनुसार काम कर रहा है। इस section में हम pipeline में दोनों processes को individual units के रूप में test करने जा रहे हैं।

### 2.1. `sayHello` process का परीक्षण करें

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

आइए अभी के लिए `main.sayhello.nf.test` file में `sayhello` process पर focus करें।

आइए file खोलें और contents को देखें।

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

    test("Should run without failures") {

        when {
            params {
                // यहाँ parameters define करें। उदाहरण:
                // outdir = "tests/results" (output directory)
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

पहले की तरह, हम test details से शुरू करते हैं, इसके बाद `when` और `then` blocks। हालांकि, हमारे पास एक अतिरिक्त `process` block भी है जो हमें process के inputs define करने की अनुमति देता है।

आइए test run करें यह देखने के लिए कि यह काम करता है या नहीं।

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

Test fail होता है क्योंकि `sayHello` process 1 input declare करती है लेकिन 0 arguments के साथ call की गई थी। आइए process में एक input जोड़कर इसे ठीक करें। [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (और ऊपर warmup section) से याद रखें कि हमारी `sayHello` process एक single value input लेती है, जिसे हमें provide करने की आवश्यकता होगी। हमें test name को भी ठीक करना चाहिए ताकि यह better reflect करे कि हम क्या test कर रहे हैं।

**पहले:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // यहाँ parameters define करें। उदाहरण:
                // outdir = "tests/results" (output directory)
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
                // यहाँ parameters define करें। उदाहरण:
                // outdir = "tests/results" (output directory)
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

आइए test को फिर से run करें यह देखने के लिए कि यह काम करता है या नहीं।

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

सफलता! Test pass होता है क्योंकि `sayHello` process successfully run हुई और output बनाया गया।

### 2.2. Test द्वारा बनाए गए snapshot को देखें

यदि हम `tests/main.sayhello.nf.test` file को देखें, तो हम देख सकते हैं कि यह assertion block में एक method `snapshot()` का उपयोग करता है:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

यह nf-test को बता रहा है कि `sayHello` process के output का एक snapshot बनाए। आइए snapshot file के contents को देखें।

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

हम इसे यहां print नहीं करेंगे, लेकिन आपको process और process outputs के details युक्त एक JSON file दिखनी चाहिए। विशेष रूप से, हम एक line देख सकते हैं जो इस तरह दिखती है:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

यह `sayHello` process द्वारा बनाए गए outputs को represent करता है, जिसे हम explicitly test कर रहे हैं। यदि हम test को फिर से run करते हैं, तो program check करेगा कि नया output उस output से match करता है जो originally record किया गया था। यह test करने का एक quick, simple तरीका है कि process outputs change नहीं होते हैं, यही कारण है कि nf-test इसे default के रूप में provide करता है।

!!!warning

    इसका मतलब है कि हमें यह सुनिश्चित करना होगा कि जो output हम original run में record करते हैं वह सही है!

यदि, future development के दौरान, code में कुछ बदलता है जो output को different बनाता है, तो test fail हो जाएगा और हमें यह निर्धारित करना होगा कि change expected है या नहीं।

- यदि यह पता चलता है कि code में कुछ टूट गया है, तो हमें इसे ठीक करना होगा, इस expectation के साथ कि fixed code test pass करेगा
- यदि यह एक expected change है (जैसे, tool को improve किया गया है और results बेहतर हैं) तो हमें new output को reference के रूप में match करने के लिए snapshot को update करने की आवश्यकता होगी। nf-test के पास इस purpose के लिए एक parameter `--update-snapshot` है

हम test को फिर से run कर सकते हैं और देख सकते हैं कि test pass होना चाहिए:

```console title="nf-test process snapshot के साथ pass"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

सफलता! Test pass होता है क्योंकि `sayHello` process successfully run हुई और output snapshot से match हो गया।

### 2.3. Snapshots का विकल्प: Direct Content Assertions

जबकि snapshots output में किसी भी changes को catch करने के लिए बहुत अच्छे हैं, कभी-कभी आप specific content को verify करना चाहते हैं बिना entire file के matching के बारे में इतने strict होने के। उदाहरण के लिए:

- जब output के parts change हो सकते हैं (timestamps, random IDs, आदि) लेकिन कुछ key content present होना चाहिए
- जब आप output में specific patterns या values check करना चाहते हैं
- जब आप test को इस बारे में अधिक explicit बनाना चाहते हैं कि success क्या है

यहां बताया गया है कि हम specific content को check करने के लिए अपने test को कैसे modify कर सकते हैं:

**पहले:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // यहाँ parameters define करें। उदाहरण:
                // outdir = "tests/results" (output directory)
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

ध्यान दें कि nf-test process outputs को lists की lists के रूप में देखता है, इसलिए `process.out[0][0]` इस process से पहले channel item (या 'emission') के पहले भाग को fetch कर रहा है।

यह approach:

- यह clear करता है कि हम output में exactly क्या expect करते हैं
- Output में irrelevant changes के लिए अधिक resilient है
- जब tests fail होते हैं तो बेहतर error messages provide करता है
- अधिक complex validations (regex patterns, numerical comparisons, आदि) की अनुमति देता है

आइए test run करें यह देखने के लिए कि यह काम करता है या नहीं।

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

### 2.4. `convertToUpper` process का परीक्षण करें

आइए `tests/main.converttoupper.nf.test` file खोलें और contents को देखें:

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

    test("Should run without failures") {

        when {
            params {
                // यहाँ parameters define करें। उदाहरण:
                // outdir = "tests/results" (output directory)
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

यह `sayHello` process के समान test है, लेकिन यह `convertToUpper` process को test कर रहा है। हम जानते हैं कि यह fail हो जाएगा क्योंकि `sayHello` की तरह, `convertToUpper` process एक single path input लेती है, लेकिन हमने एक specify नहीं की है।

अब हमें convertToUpper process को एक single input file supply करने की आवश्यकता है, जिसमें कुछ text शामिल है जिसे हम uppercase में convert करना चाहते हैं। हम ऐसा करने के कई तरीके हैं:

- हम test करने के लिए एक dedicated file बना सकते हैं
- हम existing data/greetings.csv file को फिर से उपयोग कर सकते हैं
- हम इसे test के भीतर on the fly बना सकते हैं

अभी के लिए, आइए pipeline level test के साथ उपयोग किए गए example का उपयोग करते हुए existing data/greetings.csv file को फिर से उपयोग करें। पहले की तरह, हम test का नाम बेहतर reflect करने के लिए name दे सकते हैं कि हम क्या test कर रहे हैं, लेकिन इस बार आइए इसे content को 'snapshot' करने दें बजाय specific strings check करने के (जैसा कि हमने दूसरी process में किया था)।

**पहले:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // यहाँ parameters define करें। उदाहरण:
                // outdir = "tests/results" (output directory)
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
                // यहाँ parameters define करें। उदाहरण:
                // outdir = "tests/results" (output directory)
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

और test run करें!

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

ध्यान दें, हमने `tests/main.converttoupper.nf.test.snap` पर `convertToUpper` process के लिए एक snapshot file बनाई है। यदि हम test को फिर से run करते हैं, तो हमें nf-test फिर से pass होता देखना चाहिए।

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

### निष्कर्ष

आप जानते हैं कि Nextflow process के लिए tests कैसे लिखें और उन्हें कैसे run करें।

### आगे क्या है?

सभी चीजों के लिए एक साथ tests run करना सीखें!

## 3. पूरे repository के लिए tests run करें

प्रत्येक component पर nf-test run करना ठीक है, लेकिन laborious और error prone है। क्या हम एक साथ सब कुछ test नहीं कर सकते?

हां हम कर सकते हैं!

आइए पूरे repo पर nf-test run करें।

### 3.1. पूरे repo पर nf-test run करें

हम `nf-test test` command run करके पूरे repo पर nf-test run कर सकते हैं।

```bash
nf-test test .
```

ध्यान दें, हम अपनी current directory से सब कुछ run करने के लिए केवल `.` का उपयोग कर रहे हैं। इसमें हर test शामिल होगा!

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

यह देखिए! हमने 4 tests run किए, प्रत्येक process के लिए 1 और एक single command के साथ पूरी pipeline के लिए 2। कल्पना करें कि एक बड़े codebase पर यह कितना powerful है!

---

## सारांश

इस side quest में, आपने व्यक्तिगत processes के साथ-साथ पूरी pipeline के लिए end-to-end tests बनाने और run करने के लिए nf-test की features का leverage करना सीखा है।
अब आप output validation के मुख्य दो approaches, snapshots और direct content assertions से अवगत हैं, और कब किसका उपयोग करना है।
आप यह भी जानते हैं कि tests को एक-एक करके या पूरे project के लिए कैसे run करें।

अपने स्वयं के काम में इन techniques को apply करना आपको यह सुनिश्चित करने में सक्षम बनाएगा कि:

- आपका code अपेक्षा के अनुसार काम करता है
- Changes existing functionality को break नहीं करते हैं
- अन्य developers confidence के साथ contribute कर सकते हैं
- Problems को जल्दी identify और fix किया जा सकता है
- Output content expectations से match करता है

### प्रमुख patterns

1. Pipeline-level tests:
   - Basic success testing
   - Process count verification
   - Output file existence checks
2. Process-level tests
3. Output validation के दो approaches:
   - Complete output verification के लिए snapshots का उपयोग करना
   - Specific content checks के लिए direct content assertions का उपयोग करना
4. एक single command के साथ repository में सभी tests run करना

### अतिरिक्त संसाधन

अधिक advanced testing features और best practices के लिए [nf-test documentation](https://www.nf-test.com/) देखें। आप चाह सकते हैं:

- अपने tests में अधिक comprehensive assertions जोड़ें
- Edge cases और error conditions के लिए tests लिखें
- Tests को automatically run करने के लिए continuous integration set up करें
- Workflow और module tests जैसे अन्य प्रकार के tests के बारे में जानें
- अधिक advanced content validation techniques का पता लगाएं

**याद रखें:** Tests आपके code के व्यवहार का living documentation हैं। आप जितने अधिक tests लिखते हैं, और आपके assertions जितने अधिक specific हैं, आप अपनी pipeline की reliability में उतने ही confident हो सकते हैं।

---

## आगे क्या है?

[Side Quests के menu](./index.md) पर वापस लौटें या list में अगले topic पर जाने के लिए page के निचले दाईं ओर button पर click करें।
