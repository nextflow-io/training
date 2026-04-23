---
title: Plugin Development
hide:
  - toc
---

# Plugin Development

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow का plugin सिस्टम तुम्हें भाषा को custom फ़ंक्शन, monitoring hooks, execution backends, और बहुत कुछ के साथ extend करने की सुविधा देता है।
Plugins समुदाय को Nextflow के core को modify किए बिना उसमें features जोड़ने में सक्षम बनाते हैं, जिससे वे पाइपलाइनों में reusable functionality साझा करने के लिए आदर्श बन जाते हैं।

इस प्रशिक्षण के दौरान, तुम सीखोगे कि मौजूदा plugins का उपयोग कैसे करें और वैकल्पिक रूप से अपना खुद का plugin कैसे बनाएं।

## Audience & prerequisites

Part 1 में मौजूदा plugins का उपयोग शामिल है और यह सभी Nextflow उपयोगकर्ताओं के लिए प्रासंगिक है।
Parts 2-6 में अपना खुद का plugin बनाना शामिल है और इसमें Groovy कोड और build tools का उपयोग होता है।
पहले से Java या Groovy का अनुभव होना ज़रूरी नहीं है।

**Prerequisites**

- एक GitHub account या [यहाँ](../../envsetup/02_local) बताए अनुसार लोकल इंस्टॉलेशन।
- [Hello Nextflow](../../hello_nextflow/index.md) कोर्स पूरा किया हो या समकक्ष अनुभव हो।
- Java 21 या उससे नया संस्करण (प्रशिक्षण वातावरण में शामिल है; केवल Parts 2-6 के लिए आवश्यक है)।

**Working directory:** `side-quests/plugin_development`

## Learning objectives

इस प्रशिक्षण के अंत तक, तुम यह करने में सक्षम होगे:

**Plugins का उपयोग करना (Part 1):**

- अपने वर्कफ़्लो में मौजूदा plugins को install और configure करना
- Plugin फ़ंक्शन को import और उपयोग करना

**Plugins develop करना (Parts 2-6):**

- Nextflow के built-in project generator का उपयोग करके एक नया plugin प्रोजेक्ट बनाना
- वर्कफ़्लो से callable custom फ़ंक्शन implement करना
- अपने plugin को locally build, test, और install करना
- वर्कफ़्लो events (जैसे, कार्य पूर्णता, पाइपलाइन start/end) को custom logging या notifications के लिए monitor करना
- Plugins को customizable बनाने के लिए configuration options जोड़ना
- अपना plugin distribute करना

## Lesson plan

#### Part 1: Plugin basics

एक Nextflow वर्कफ़्लो में मौजूदा plugins का उपयोग करें और उनके व्यवहार को configure करें।

#### Part 2: Create a plugin project

एक नया plugin प्रोजेक्ट generate करें और उसकी संरचना की जांच करें।

#### Part 3: Custom functions

Custom फ़ंक्शन implement करें, अपना plugin build करें, और इसे एक वर्कफ़्लो में चलाएं।

#### Part 4: Testing

Spock framework का उपयोग करके unit tests लिखें और चलाएं।

#### Part 5: Workflow monitoring

एक task counter बनाने के लिए कार्य पूर्णता जैसे events पर प्रतिक्रिया दें।

#### Part 6: Configuration & Distribution

अपने plugin को customizable बनाने के लिए `nextflow.config` से settings पढ़ें, फिर इसे साझा करना सीखें।

कोर्स शुरू करने के लिए तैयार हो?

[सीखना शुरू करें :material-arrow-right:](01_plugin_basics.md){ .md-button .md-button--primary }
