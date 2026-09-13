---
title: nf-core का उपयोग करें
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - nf-core कम्युनिटी पाइपलाइन खोजना, प्राप्त करना और चलाना
    - पैरामीटर और कॉन्फ़िगरेशन फ़ाइलों का उपयोग करके पाइपलाइन execution को कॉन्फ़िगर करना
    - यह समझना कि nf-core पाइपलाइन पैरामीटर और इनपुट डेटा को कैसे validate करती हैं
    - एक production-scale पाइपलाइन (nf-core/rnaseq) चलाना और उसके डिफ़ॉल्ट resource allocations को override करना
  audience_prerequisites:
    - "**दर्शक:** यह कोर्स उन लोगों के लिए है जो पहले से local Nextflow पाइपलाइन चलाना जानते हैं और nf-core में नए हैं, और मौजूदा कम्युनिटी पाइपलाइन चलाना चाहते हैं।"
    - "**कौशल:** कमांड लाइन, बेसिक scripting अवधारणाओं और सामान्य फ़ाइल फ़ॉर्मेट से कुछ परिचय होना ज़रूरी है।"
    - "**कोर्स:** [Nextflow Run](../nextflow_run/index.md) पूरा किया होना चाहिए या `nextflow run` के साथ local पाइपलाइन चलाने में सहज होना चाहिए।"
    - "**डोमेन:** अभ्यासों में bioinformatics पाइपलाइन का उपयोग किया गया है, लेकिन किसी वैज्ञानिक डोमेन ज्ञान की आवश्यकता नहीं है।"
---

# nf-core का उपयोग करें

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Use nf-core, nf-core कम्युनिटी पाइपलाइन खोजने, चलाने और कॉन्फ़िगर करने का एक व्यावहारिक परिचय है।**

व्यावहारिक उदाहरणों और guided अभ्यासों के माध्यम से, तुम nf-core पाइपलाइन खोजना और प्राप्त करना, उन्हें उनके built-in test profiles का उपयोग करके चलाना, और पैरामीटर व कॉन्फ़िगरेशन फ़ाइलों के ज़रिए उनके execution को customize करना सीखोगे।

तुम अपने खुद के analyses के लिए nf-core पाइपलाइन चलाना शुरू करने के लिए ज़रूरी कौशल और आत्मविश्वास लेकर जाओगे।

<!-- additional_information -->

## कोर्स का अवलोकन

यह कोर्स व्यावहारिक है, जिसमें goal-oriented अभ्यास हैं जो जानकारी को धीरे-धीरे प्रस्तुत करते हैं।

तुम `nf-core/demo` से शुरू करोगे, जो nf-core प्रोजेक्ट द्वारा प्रशिक्षण उद्देश्यों के लिए maintained एक minimal पाइपलाइन है, फिर जो सीखा है उसे `nf-core/rnaseq` पर लागू करोगे — जो bulk RNA sequencing analysis के लिए एक widely-used production पाइपलाइन है।

यह कोर्स पाइपलाइन चलाने पर केंद्रित है।
अगर तुम nf-core-compatible पाइपलाइन develop करने का परिचय ढूंढ रहे हो, तो [Build with nf-core](../nfcore_build/index.md) देखो।

### पाठ योजना

| कोर्स अध्याय                                                         | सारांश                                                                                                    | अनुमानित समय |
| -------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- | ------------ |
| [भाग 1: एक demo पाइपलाइन चलाएं](./01_run_demo.md)                   | एक nf-core पाइपलाइन खोजें और प्राप्त करें, और उसे test profile का उपयोग करके चलाएं                      | 20 मिनट      |
| [भाग 2: पाइपलाइन execution को कॉन्फ़िगर करें](./02_configure_execution.md) | पैरामीटर सेट करें, validation समझें, और resource allocation व tool arguments को customize करें      | 20 मिनट      |
| [भाग 3: एक production पाइपलाइन चलाएं](./03_run_production_pipeline.md) | nf-core/rnaseq को pull और run करें, और उसके डिफ़ॉल्ट resource allocations को override करें           | 20 मिनट      |

इस कोर्स के अंत तक, तुम nf-core प्रोजेक्ट द्वारा प्रदान की जाने वाली कम्युनिटी पाइपलाइन की विशाल संपदा का लाभ उठाने में सक्षम हो जाओगे।

क्या तुम कोर्स शुरू करने के लिए तैयार हो?

[सीखना शुरू करें :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
