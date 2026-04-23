# भाग 1: Plugin की मूल बातें

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

इस अनुभाग में, तुम सीखोगे कि plugins Nextflow को कैसे विस्तारित करते हैं, फिर तीन अलग-अलग plugins को काम करते हुए देखने के लिए आज़माओगे।

---

## 1. Plugins कैसे काम करते हैं

Plugins कई प्रकार के एक्सटेंशन के माध्यम से Nextflow को विस्तारित करते हैं:

| एक्सटेंशन का प्रकार | यह क्या करता है                                  | उदाहरण                       |
| ------------------- | ------------------------------------------------ | ---------------------------- |
| Functions           | वर्कफ़्लो से callable कस्टम functions जोड़ता है  | `samplesheetToList()`        |
| Workflow monitors   | कार्य पूर्णता जैसी events पर प्रतिक्रिया देता है | Custom logging, Slack alerts |
| Executors           | कार्य execution backends जोड़ता है               | AWS Batch, Kubernetes        |
| Filesystems         | स्टोरेज backends जोड़ता है                       | S3, Azure Blob               |

Functions और workflow monitors (Nextflow API में "trace observers" कहलाते हैं) plugin लेखकों के लिए सबसे सामान्य प्रकार हैं।
Executors और filesystems आमतौर पर platform vendors द्वारा बनाए जाते हैं।

अगले अभ्यास तुम्हें function plugins और एक observer plugin दिखाएंगे, ताकि तुम दोनों प्रकारों को काम करते हुए देख सको।

---

## 2. Function plugins का उपयोग करना

Function plugins callable functions जोड़ते हैं जिन्हें तुम अपने वर्कफ़्लो में import करते हो।
तुम दो आज़माओगे: nf-hello (एक सरल उदाहरण) और nf-schema (एक व्यापक रूप से उपयोग किया जाने वाला real-world plugin)।
दोनों अभ्यास एक ही `hello.nf` पाइपलाइन को संशोधित करते हैं, ताकि तुम देख सको कि plugins एक मौजूदा वर्कफ़्लो को कैसे बेहतर बनाते हैं।

### 2.1. nf-hello: हाथ से लिखे कोड को बदलना

[nf-hello](https://github.com/nextflow-io/nf-hello) plugin एक `randomString` function प्रदान करता है जो random strings उत्पन्न करता है।
पाइपलाइन पहले से ही इस function का अपना inline संस्करण परिभाषित करती है, जिसे तुम plugin के संस्करण से बदलोगे।

#### 2.1.1. शुरुआती बिंदु देखना

पाइपलाइन देखो:

```bash
cat hello.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

params.input = 'greetings.csv'

/**
 * एक random alphanumeric string उत्पन्न करें
 */
def randomString(int length) {
    def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
}

process SAY_HELLO {
    input:
        val greeting
    output:
        stdout
    script:
    """
    echo '$greeting'
    """
}

workflow {
    greeting_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> "${row.greeting}_${randomString(8)}" }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
```

पाइपलाइन अपना `randomString` function inline परिभाषित करती है, फिर इसका उपयोग प्रत्येक greeting में एक random ID जोड़ने के लिए करती है।

इसे चलाओ:

```bash
nextflow run hello.nf
```

```console title="Output"
Output: Hello_aBcDeFgH
Output: Bonjour_xYzWvUtS
Output: Holà_qRsPdMnK
Output: Ciao_jLhGfEcB
Output: Hallo_tNwOiAuR
```

तुम्हारे आउटपुट का क्रम और random strings अलग होंगे, और यदि तुम script को फिर से चलाओगे तो तुम्हें greetings का एक अलग सेट मिलेगा।

#### 2.1.2. Plugin को configure करना

Inline function को plugin के function से बदलो। अपने `nextflow.config` में यह जोड़ो:

```groovy title="nextflow.config"
// Plugin development अभ्यासों के लिए configuration
plugins {
    id 'nf-hello@0.5.0'
}
```

Plugins को `nextflow.config` में `plugins {}` block का उपयोग करके घोषित किया जाता है।
Nextflow उन्हें [Nextflow Plugin Registry](https://registry.nextflow.io/) से स्वचालित रूप से डाउनलोड करता है, जो community और official plugins का एक केंद्रीय भंडार है।

#### 2.1.3. Plugin function का उपयोग करना

Inline `randomString` function को plugin संस्करण से बदलो:

=== "बाद में"

    ```groovy title="hello.nf" hl_lines="3"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "पहले"

    ```groovy title="hello.nf" hl_lines="5-12"
    #!/usr/bin/env nextflow

    params.input = 'greetings.csv'

    /**
     * एक random alphanumeric string उत्पन्न करें
     */
    def randomString(int length) {
        def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
        def random = new Random()
        return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
    }

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

`include` statement `randomString` को एक ऐसी library से import करता है जो सिद्ध, परीक्षित और contributors के एक व्यापक समूह द्वारा maintained है जो bugs को पकड़ और ठीक कर सकते हैं।
प्रत्येक पाइपलाइन अपनी खुद की function की copy रखने के बजाय, plugin का उपयोग करने वाली हर पाइपलाइन को एक ही परीक्षित implementation मिलती है।
इससे duplicate कोड और उससे जुड़ा maintenance का बोझ कम होता है।
`#!groovy include { function } from 'plugin/plugin-id'` syntax वही `include` है जो Nextflow modules के लिए उपयोग किया जाता है, बस `plugin/` prefix के साथ।
तुम GitHub पर nf-hello repository में [`randomString` का source code](https://github.com/nextflow-io/nf-hello/blob/e67bddebfa589c7ae51f41bf780c92068dc09e93/plugins/nf-hello/src/main/nextflow/hello/HelloExtension.groovy#L110) देख सकते हो।

#### 2.1.4. इसे चलाओ

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_yqvtclcc
Output: Bonjour_vwwpyzcs
Output: Holà_wrghmgab
Output: Ciao_noniajuy
Output: Hallo_tvrtuxtp
Pipeline complete! 👋
```

(तुम्हारी random strings अलग होंगी।)

आउटपुट में अभी भी random suffixes हैं, लेकिन अब `randomString` inline कोड के बजाय nf-hello plugin से आता है।
"Pipeline is starting!" और "Pipeline complete!" संदेश नए हैं।
ये plugin के observer component से आते हैं, जिसे तुम भाग 5 में explore करोगे।

Nextflow plugins को पहली बार उपयोग होने पर स्वचालित रूप से डाउनलोड करता है, इसलिए `nf-hello@0.5.0` घोषित करने वाली कोई भी पाइपलाइन projects के बीच कोड copy किए बिना वही परीक्षित `randomString` function प्राप्त करती है।

तुमने अब function plugin उपयोग करने के तीन चरण देखे हैं: इसे `nextflow.config` में घोषित करो, `include` के साथ function import करो, और इसे अपने वर्कफ़्लो में call करो।
अगला अभ्यास इन्हीं चरणों को एक real-world plugin पर लागू करता है।

### 2.2. nf-schema: validated CSV parsing

[nf-schema](https://github.com/nextflow-io/nf-schema) plugin सबसे व्यापक रूप से उपयोग किए जाने वाले Nextflow plugins में से एक है।
यह `samplesheetToList` प्रदान करता है, एक function जो CSV/TSV फ़ाइलों को JSON schema का उपयोग करके parse करता है जो अपेक्षित columns और types को परिभाषित करता है।

पाइपलाइन वर्तमान में `splitCsv` और एक manual `map` का उपयोग करके `greetings.csv` पढ़ती है, लेकिन nf-schema इसे validated, schema-driven parsing से बदल सकता है।
एक JSON schema फ़ाइल (`greetings_schema.json`) अभ्यास डायरेक्टरी में पहले से प्रदान की गई है।

??? info "Schema क्या है?"

    Schema valid data कैसी दिखती है इसका एक औपचारिक विवरण है।
    यह परिभाषित करता है कि कौन से columns अपेक्षित हैं, प्रत्येक मान का प्रकार क्या होना चाहिए (string, number, आदि), और कौन से fields आवश्यक हैं।

    इसे एक अनुबंध के रूप में सोचो: यदि इनपुट डेटा schema से मेल नहीं खाता, तो tool पाइपलाइन में बाद में भ्रामक errors का कारण बनने देने के बजाय समस्या को जल्दी पकड़ सकता है।

#### 2.2.1. Schema देखो

```bash
cat greetings_schema.json
```

```json title="Output"
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "greeting": {
        "type": "string",
        "description": "The greeting text"
      },
      "language": {
        "type": "string",
        "description": "The language of the greeting"
      }
    },
    "required": ["greeting"]
  }
}
```

Schema दो columns (`greeting` और `language`) परिभाषित करता है और `greeting` को आवश्यक चिह्नित करता है।
यदि कोई `greeting` column के बिना CSV पास करता है, तो nf-schema पाइपलाइन चलने से पहले error पकड़ लेता है।

#### 2.2.2. Config में nf-schema जोड़ो

दोनों plugins शामिल करने के लिए `nextflow.config` अपडेट करो:

=== "बाद में"

    ```groovy title="nextflow.config" hl_lines="3"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
    }
    ```

#### 2.2.3. samplesheetToList उपयोग करने के लिए hello.nf अपडेट करो

`splitCsv` इनपुट को `samplesheetToList` से बदलो:

=== "बाद में"

    ```groovy title="hello.nf" hl_lines="4 20 21 22"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'
    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        def samplesheet_list = samplesheetToList(params.input, 'greetings_schema.json')
        greeting_ch = Channel.fromList(samplesheet_list)
            .map { row -> "${row[0]}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "पहले"

    ```groovy title="hello.nf" hl_lines="19 20 21"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

कस्टम `splitCsv` और `map` parsing कोड को `samplesheetToList` से बदला गया है, एक सिद्ध, परीक्षित function जो पाइपलाइन चलने से पहले samplesheet को schema के विरुद्ध validate भी करता है।
यह हाथ से लिखे parsing logic के maintenance के बोझ को कम करता है और साथ ही पाइपलाइन उपयोगकर्ताओं के अनुभव को बेहतर बनाता है, जिन्हें स्पष्ट error messages मिलते हैं जब उनका इनपुट अपेक्षित format से मेल नहीं खाता।
प्रत्येक row column क्रम में values की एक list बन जाती है, इसलिए `row[0]` greeting है और `row[1]` language है।

#### 2.2.4. इसे चलाओ

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_diozjdwm
Output: Bonjour_speathmm
Output: Holà_dllxnzap
Output: Ciao_wzueddzc
Output: Hallo_hsxwrjbh
Pipeline complete! 👋
```

(तुम्हारी random strings अलग होंगी।)

आउटपुट वही है, लेकिन अब schema पाइपलाइन चलने से पहले CSV structure को validate करता है।
जटिल sample sheets और कई columns वाली real pipelines में, इस प्रकार की validation उन errors को रोकती है जो manual `splitCsv` + `map` से छूट जाती हैं।

#### 2.2.5. Validation को काम करते देखो

Schema validation क्या पकड़ता है यह देखने के लिए, `greetings.csv` में errors डालने की कोशिश करो।

आवश्यक `greeting` column का नाम बदलकर `message` करो:

```csv title="greetings.csv" hl_lines="1"
message,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

पाइपलाइन चलाओ:

```bash
nextflow run hello.nf
```

```console title="Output"
ERROR ~ Validation of samplesheet failed!

The following errors have been detected in greetings.csv:

-> Entry 1: Missing required field(s): greeting
-> Entry 2: Missing required field(s): greeting
-> Entry 3: Missing required field(s): greeting
-> Entry 4: Missing required field(s): greeting
-> Entry 5: Missing required field(s): greeting
```

पाइपलाइन चलने से मना कर देती है क्योंकि schema को `greeting` column चाहिए और वह नहीं मिल रही।

अब आवश्यक column को वापस करो लेकिन optional `language` column का नाम बदलकर `lang` करो:

```csv title="greetings.csv" hl_lines="1"
greeting,lang
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```bash
nextflow run hello.nf
```

इस बार पाइपलाइन चलती है, लेकिन एक warning प्रिंट करती है:

```console title="Output (partial)"
WARN: Found the following unidentified headers in greetings.csv:
	- lang
```

आवश्यक columns hard errors का कारण बनते हैं; optional columns warnings का कारण बनते हैं।
यह वह प्रकार की early feedback है जो दर्जनों columns वाली real pipelines में debugging का समय बचाती है।

#### 2.2.6. Validation behavior configure करो

`lang` के बारे में warning उपयोगी है, लेकिन तुम configuration के माध्यम से इसकी severity को नियंत्रित कर सकते हो।
Plugins अपने स्वयं के configuration scope(s) शामिल कर सकते हैं जो उनके behavior को नियंत्रित करते हैं।
nf-schema plugin में `validation` configuration scope शामिल है; यहाँ settings को संशोधित करके तुम nf-schema के behavior को बदल सकते हो।

अपरिचित headers को warning के बजाय error का कारण बनाने के लिए `nextflow.config` में एक `validation` block जोड़ो:

=== "बाद में"

    ```groovy title="nextflow.config" hl_lines="6-10"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }

    validation {
        logging {
            unrecognisedHeaders = "error"
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

उसी `lang` column के साथ पाइपलाइन फिर से चलाओ:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
Found the following unidentified headers in greetings.csv:
	- lang
 -- Check script 'hello.nf' at line: 20 or see '.nextflow.log' file for more details
```

पाइपलाइन अब warning देने के बजाय fail हो जाती है।
पाइपलाइन कोड नहीं बदला; केवल configuration बदली।

आगे बढ़ने से पहले `greetings.csv` को उसकी मूल स्थिति में वापस करो और `validation` block हटाओ:

```csv title="greetings.csv"
greeting,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```groovy title="nextflow.config"
plugins {
    id 'nf-hello@0.5.0'
    id 'nf-schema@2.6.1'
}
```

nf-hello और nf-schema दोनों function plugins हैं: वे functions प्रदान करते हैं जिन्हें तुम `include` के साथ import करते हो और अपने वर्कफ़्लो कोड में call करते हो।
अगला अभ्यास एक अलग प्रकार का plugin दिखाता है जो बिना किसी `include` statement के काम करता है।

---

## 3. Observer plugin का उपयोग करना: nf-co2footprint

सभी plugins import करने के लिए functions प्रदान नहीं करते।
[nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint) plugin एक **trace observer** का उपयोग करके तुम्हारी पाइपलाइन के resource उपयोग की निगरानी करता है और उसके carbon footprint का अनुमान लगाता है।
तुम्हें कोई भी पाइपलाइन कोड बदलने की ज़रूरत नहीं है; बस इसे config में जोड़ो।

### 3.1. Config में nf-co2footprint जोड़ो

`nextflow.config` अपडेट करो:

=== "बाद में"

    ```groovy title="nextflow.config" hl_lines="4"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

### 3.2. पाइपलाइन चलाओ

```bash
nextflow run hello.nf
```

Plugin execution के दौरान कई INFO और WARN messages उत्पन्न करता है।
ये local machine पर चलने वाले एक छोटे उदाहरण के लिए सामान्य हैं:

```console title="Output (partial)"
nf-co2footprint plugin  ~  version 1.2.0
WARN - [nf-co2footprint] Target zone null not found. Attempting to retrieve carbon intensity for fallback zone GLOBAL.
INFO - [nf-co2footprint] Using fallback carbon intensity from GLOBAL from CI table: 480.0 gCO₂eq/kWh.
WARN - [nf-co2footprint] Executor 'null' not mapped.
WARN - [nf-co2footprint] Fallback to: `machineType = null`, `pue = 1.0`. ...
...
WARN - [nf-co2footprint] No CPU model detected. Using default CPU power draw value (11.41 W).
WARN - [nf-co2footprint] 🔁 Requested memory is null for task 2. Using maximum consumed memory/`peak_rss` (0 GB) for CO₂e footprint computation.
```

Zone, executor, CPU model, और memory के बारे में warnings इसलिए आती हैं क्योंकि plugin local training environment के पूर्ण hardware विवरण का पता नहीं लगा सकता।
एक production environment में (जैसे, HPC cluster या cloud), ये values उपलब्ध होंगी और अनुमान अधिक सटीक होंगे।

अंत में, इस तरह की एक line देखो:

```console title="Output (partial)"
🌱 The workflow run used 126.76 uWh of electricity, resulting in the release of 60.84 ug of CO₂ equivalents into the atmosphere.
```

(तुम्हारे numbers अलग होंगे।)

### 3.3. Report देखो

Plugin तुम्हारी working directory में आउटपुट फ़ाइलें उत्पन्न करता है:

```bash
ls co2footprint_*
```

```console title="Output"
co2footprint_report_<timestamp>.html
co2footprint_summary_<timestamp>.txt
co2footprint_trace_<timestamp>.txt
```

Summary देखो:

```bash
cat co2footprint_summary_*.txt
```

```console title="Output"
Total CO₂e footprint measures of this workflow run (including cached tasks):
  CO₂e emissions: 60.84 ug
  Energy consumption: 126.76 uWh
  CO₂e emissions (market): -

Which equals:
  - 3.48E-7 km travelled by car
  - It takes one tree 0.17s to sequester the equivalent amount of CO₂ from the atmosphere
  - 1.22E-7 % of a flight from Paris to London
```

(तुम्हारे numbers अलग होंगे।)

पहला अनुभाग raw energy और emissions figures दिखाता है।
"Which equals" अनुभाग उन numbers को परिचित equivalents में बदलकर परिप्रेक्ष्य में रखता है।
Summary में plugin के configuration options की सूची और [Green Algorithms](https://doi.org/10.1002/advs.202100707) research paper का citation भी शामिल है जिस पर calculation method आधारित है।

### 3.4. Plugin configure करो

अनुभाग 3.2 से "Target zone null" warning इसलिए आई क्योंकि plugin में कोई location configure नहीं थी।
nf-co2footprint plugin एक `co2footprint` configuration scope परिभाषित करता है जहाँ तुम अपना geographic location सेट कर सकते हो।

`nextflow.config` में एक `co2footprint` block जोड़ो:

=== "बाद में"

    ```groovy title="nextflow.config" hl_lines="7-9"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }

    co2footprint {
        location = 'GB'
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

!!! tip "सुझाव"

    यदि तुम चाहो तो अपने देश का code उपयोग करो (जैसे, `'US'`, `'DE'`, `'FR'`)।

पाइपलाइन चलाओ:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
INFO - [nf-co2footprint] Using fallback carbon intensity from GB from CI table: 163.92 gCO₂eq/kWh.
```

Zone warning चली गई।
Plugin अब global fallback (480.0 gCO₂eq/kWh) के बजाय GB-specific carbon intensity (163.92 gCO₂eq/kWh) का उपयोग करता है।

!!! note "नोट"

    तुम्हें `WARN: Unrecognized config option 'co2footprint.location'` message भी दिख सकता है।
    यह केवल दिखावटी है और इसे सुरक्षित रूप से नज़रअंदाज़ किया जा सकता है; plugin फिर भी value को सही तरीके से पढ़ता है।

भाग 6 में, तुम अपने plugin के लिए एक configuration scope बनाओगे।

यह plugin पूरी तरह से observer mechanism के माध्यम से काम करता है, resource metrics एकत्र करने और पाइपलाइन पूरी होने पर अपनी report उत्पन्न करने के लिए workflow lifecycle events में hook करता है।

तुमने अब function plugins (जो `include` के साथ import किए जाते हैं) और एक observer plugin (जो केवल config के माध्यम से सक्रिय होता है) दोनों आज़माए हैं।
ये दो सबसे सामान्य extension types हैं, लेकिन जैसा कि अनुभाग 1 की तालिका दिखाती है, plugins executors और filesystems भी जोड़ सकते हैं।

---

## 4. Plugins खोजना

[Nextflow Plugin Registry](https://registry.nextflow.io/) उपलब्ध plugins खोजने का केंद्रीय hub है।

![registry.nextflow.io पर nf-hello plugin page](img/plugin-registry-nf-hello.png)

प्रत्येक plugin page उसका विवरण, उपलब्ध versions, installation instructions, और documentation के links दिखाता है।

---

## 5. Plugin development के लिए तैयारी करना

निम्नलिखित अनुभाग (भाग 2-6) एक अलग पाइपलाइन फ़ाइल, `greet.nf`, का उपयोग करते हैं, जो nf-schema पर निर्भर करती है लेकिन nf-hello या nf-co2footprint पर नहीं।

केवल nf-schema रखने के लिए `nextflow.config` अपडेट करो:

```groovy title="nextflow.config"
// Plugin development अभ्यासों के लिए configuration
plugins {
    id 'nf-schema@2.6.1'
}
```

co2footprint आउटपुट फ़ाइलें हटाओ:

```bash
rm -f co2footprint_*
```

`hello.nf` फ़ाइल संदर्भ के लिए तुम्हारे भाग 1 के काम को बनाए रखती है; आगे से, तुम `greet.nf` के साथ काम करोगे।

---

## सारांश

तुमने तीन अलग-अलग plugins का उपयोग किया:

- **nf-hello**: एक function plugin जो `randomString` प्रदान करता है, `include` के साथ import किया गया
- **nf-schema**: एक function plugin जो schema-validated CSV parsing के लिए `samplesheetToList` प्रदान करता है
- **nf-co2footprint**: एक observer plugin जो resource उपयोग की स्वचालित रूप से निगरानी करता है, बिना किसी `include` की ज़रूरत के

मुख्य patterns:

- Plugins को `nextflow.config` में `#!groovy plugins { id 'plugin-name@version' }` के साथ घोषित किया जाता है
- Function plugins को `#!groovy include { function } from 'plugin/plugin-id'` की आवश्यकता होती है
- Observer plugins config में घोषित होने के बाद स्वचालित रूप से काम करते हैं
- Plugins configuration scopes परिभाषित कर सकते हैं (जैसे, `#!groovy validation {}`, `#!groovy co2footprint {}`) behavior को customize करने के लिए
- [Nextflow Plugin Registry](https://registry.nextflow.io/) उपलब्ध plugins की सूची देता है

---

## आगे क्या है?

निम्नलिखित अनुभाग तुम्हें दिखाते हैं कि अपना खुद का plugin कैसे बनाएं।
यदि तुम plugin development में रुचि नहीं रखते, तो तुम यहाँ रुक सकते हो या [सारांश](summary.md) पर आगे जा सकते हो।

[भाग 2 पर जारी रखो :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
