# सारांश

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

तुमने Plugin Development प्रशिक्षण पूरा कर लिया है।
यह पेज हर भाग में तुमने जो बनाया उसका पुनर्कथन करता है, वितरण को कवर करता है, और आगे कहाँ जाना है इस पर मार्गदर्शन प्रदान करता है।

---

## तुमने क्या सीखा

### भाग 1: Plugins का उपयोग करना

तुमने जाना कि Nextflow plugins एक उपयोगकर्ता के नज़रिए से कैसे काम करते हैं।
तुमने nf-schema और nf-co2footprint इंस्टॉल किए, उन्हें `nextflow.config` के ज़रिए कॉन्फ़िगर किया, और देखा कि plugins किस तरह इनपुट को validate कर सकते हैं, फ़ंक्शन जोड़ सकते हैं, और पाइपलाइन lifecycle events में hook कर सकते हैं।

### भाग 2: सेटअप करना

तुमने Java 21+ के साथ एक plugin development वातावरण सेट किया, `nextflow plugin create` कमांड का उपयोग करके एक नया plugin प्रोजेक्ट बनाया, और वह प्रोजेक्ट संरचना सीखी जो Nextflow अपेक्षित करता है: सोर्स फ़ाइलें, build कॉन्फ़िगरेशन, और Makefile वर्कफ़्लो।

### भाग 3: कस्टम फ़ंक्शन

तुमने एक `PluginExtensionPoint` क्लास में `@Function`-annotated मेथड बनाकर अपना पहला extension point लागू किया।
तुमने `reverseGreeting` और `decorateGreeting` बनाए, फिर उन्हें एक पाइपलाइन स्क्रिप्ट से import करके call किया।

### भाग 4: Testing

तुमने Groovy testing framework का उपयोग करके अपने कस्टम फ़ंक्शन के लिए unit tests लिखे।
तुमने सीखा कि `make test` से tests कैसे चलाएं और यह verify करें कि plugin इंस्टॉल करने से पहले सही तरह से काम कर रहा है।

### भाग 5: Observers

तुमने पाइपलाइन lifecycle events में hook करने के लिए `TraceObserver` इंटरफ़ेस लागू किया।
तुमने `GreetingObserver` (पाइपलाइन शुरू होने और पूरी होने पर प्रतिक्रिया देना) और `TaskCounterObserver` (पूर्ण हुए कार्यों की गिनती करना) बनाए, फिर उन्हें एक `TraceObserverFactory` के ज़रिए रजिस्टर किया।

### भाग 6: कॉन्फ़िगरेशन

तुमने runtime पर values पढ़ने के लिए `session.config.navigate()` का उपयोग करके अपने plugin को `nextflow.config` के ज़रिए कॉन्फ़िगर करने योग्य बनाया।
तुमने अपने plugin के विकल्पों को औपचारिक रूप से घोषित करने के लिए एक `@ConfigScope` क्लास जोड़ी, जिससे "Unrecognized config option" चेतावनियाँ समाप्त हुईं और IDE सपोर्ट सक्षम हुआ।

---

## वितरण

एक बार जब तुम्हारा plugin लोकल पर काम करने लगे, तो तुम इसे Nextflow plugin registry के ज़रिए दूसरों के साथ साझा कर सकते हो।

### Versioning

अपने releases के लिए [semantic versioning](https://semver.org/) का पालन करो:

| Version बदलाव             | कब उपयोग करें                     | उदाहरण                                         |
| ------------------------- | --------------------------------- | ---------------------------------------------- |
| **MAJOR** (1.0.0 → 2.0.0) | Breaking changes                  | कोई फ़ंक्शन हटाना, return types बदलना          |
| **MINOR** (1.0.0 → 1.1.0) | नई सुविधाएं, backward compatible  | एक नया फ़ंक्शन जोड़ना                           |
| **PATCH** (1.0.0 → 1.0.1) | Bug fixes, backward compatible    | मौजूदा फ़ंक्शन में bug ठीक करना                |

हर release से पहले `build.gradle` में version अपडेट करो:

```groovy title="build.gradle"
version = '1.0.0'  // Semantic versioning का उपयोग करें: MAJOR.MINOR.PATCH
```

### Registry पर Publishing

[Nextflow plugin registry](https://registry.nextflow.io/) community के साथ plugins साझा करने का आधिकारिक तरीका है।

Publishing वर्कफ़्लो:

1. [registry](https://registry.nextflow.io/) पर अपना plugin नाम claim करो (अपने GitHub account से साइन इन करो)
2. `~/.gradle/gradle.properties` में अपनी API credentials कॉन्फ़िगर करो
3. सब कुछ सही काम कर रहा है यह verify करने के लिए tests चलाओ: `make test`
4. `make release` से publish करो

चरण-दर-चरण निर्देशों के लिए, [आधिकारिक publishing दस्तावेज़ीकरण](https://www.nextflow.io/docs/latest/guides/gradle-plugin.html#publishing-a-plugin) देखो।

एक बार publish होने के बाद, उपयोगकर्ता बिना किसी लोकल सेटअप के तुम्हारा plugin इंस्टॉल कर सकते हैं:

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting@1.0.0'
}
```

Nextflow पहले उपयोग पर registry से plugin अपने आप डाउनलोड कर लेता है।

---

## Plugin development चेकलिस्ट

- [ ] Java 21+ इंस्टॉल है
- [ ] `nextflow plugin create <name> <org>` से प्रोजेक्ट बनाओ
- [ ] `@Function` मेथड के साथ extension क्लास लागू करो
- [ ] Unit tests लिखो और `make test` से चलाओ
- [ ] `make install` से build और install करो
- [ ] वैकल्पिक रूप से workflow events के लिए `TraceObserver` implementations जोड़ो
- [ ] वैकल्पिक रूप से plugin कॉन्फ़िगरेशन के लिए `ConfigScope` जोड़ो
- [ ] `nextflow.config` में `plugins { id 'plugin-id' }` से सक्षम करो
- [ ] `include { fn } from 'plugin/plugin-id'` से फ़ंक्शन import करो
- [ ] Version करो और registry पर publish करो

---

## मुख्य code patterns

**फ़ंक्शन परिभाषा:**

```groovy
@Function
String myFunction(String input, String optional = 'default') {
    return input.transform()
}
```

**Plugin कॉन्फ़िगरेशन:**

```groovy
nextflowPlugin {
    provider = 'my-org'
    className = 'my.org.MyPlugin'
    extensionPoints = ['my.org.MyExtension']
}
```

**Workflows में उपयोग:**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { item -> myFunction(item) }
        .view()
}
```

---

## Extension point सारांश

| प्रकार              | Class/Annotation | उद्देश्य                                              |
| ------------------- | ---------------- | ----------------------------------------------------- |
| Function            | `@Function`      | Workflows से call किया जा सकता है                     |
| Trace Observer      | `TraceObserver`  | Workflow lifecycle events में hook करना               |
| Configuration Scope | `@ScopeName`     | nextflow.config में plugin कॉन्फ़िगरेशन परिभाषित करना |

---

## आगे क्या है?

यहाँ तुम्हारी plugin development यात्रा जारी रखने के लिए कुछ व्यावहारिक अगले कदम हैं।

**कुछ असली बनाओ।**
अपने काम से एक use case चुनो: एक कस्टम फ़ंक्शन जिसे तुम्हारी टीम बार-बार उपयोग करती है, एक observer जो पाइपलाइन पूरी होने पर Slack notifications भेजता है, या एक config scope जो तुम्हारे संगठन की पाइपलाइनों में विकल्पों को मानकीकृत करता है।
किसी असली समस्या से शुरू करना तुम्हारी समझ को गहरा करने का सबसे तेज़ तरीका है।

**nf-hello को संदर्भ के रूप में उपयोग करो।**
[nf-hello](https://github.com/nextflow-io/nf-hello) रिपॉज़िटरी आधिकारिक minimal plugin उदाहरण है।
यह नए प्रोजेक्ट के लिए एक अच्छा शुरुआती बिंदु है और जब तुम्हें यह जाँचना हो कि कुछ कैसे संरचित है तो एक उपयोगी संदर्भ है।

**आधिकारिक दस्तावेज़ीकरण पढ़ो।**
Nextflow docs इस प्रशिक्षण से परे विषयों को कवर करते हैं, जिनमें channel factories, operator overloading, और advanced observer patterns शामिल हैं।
[developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html) गाइड सबसे व्यापक संदर्भ है।

**मौजूदा plugins का अध्ययन करो।**
[Nextflow plugins रिपॉज़िटरी](https://github.com/nextflow-io/plugins) में nf-schema, nf-wave, और nf-tower जैसे आधिकारिक plugins का सोर्स कोड है।
Production plugin code पढ़ना उन patterns और conventions को सीखने के सर्वोत्तम तरीकों में से एक है जो परिचयात्मक उदाहरणों से परे जाते हैं।

---

## अतिरिक्त संसाधन

**आधिकारिक दस्तावेज़ीकरण:**

- [Plugins का उपयोग करना](https://www.nextflow.io/docs/latest/plugins/plugins.html): plugins इंस्टॉल करने और कॉन्फ़िगर करने के लिए व्यापक गाइड
- [Plugins विकसित करना](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html): विस्तृत plugin development संदर्भ
- [Config scopes](https://nextflow.io/docs/latest/developer/config-scopes.html): plugins के लिए configuration scopes बनाना

**Plugin खोज:**

- [Nextflow Plugin Registry](https://registry.nextflow.io/): उपलब्ध plugins ब्राउज़ करो और खोजो
- [Plugin registry docs](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html): registry दस्तावेज़ीकरण

**उदाहरण और संदर्भ:**

- [nf-hello](https://github.com/nextflow-io/nf-hello): सरल उदाहरण plugin (शुरुआत के लिए बढ़िया)
- [Nextflow plugins रिपॉज़िटरी](https://github.com/nextflow-io/plugins): संदर्भ के लिए आधिकारिक plugins का संग्रह
