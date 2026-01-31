---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Nextflow workflows को लॉन्च और मैनेज करना
    - आउटपुट (परिणाम) और लॉग फ़ाइलों को खोजना और समझना
    - एक साधारण मल्टी-स्टेप workflow में मुख्य Nextflow कॉम्पोनेंट्स को पहचानना
    - HPC और cloud सहित सामान्य कंप्यूटिंग प्लेटफॉर्म पर pipeline execution को कॉन्फ़िगर करना
    - reproducibility, portability और कोड पुन: उपयोग के लिए सर्वोत्तम प्रथाओं को समझना जो pipelines को FAIR बनाती हैं, जिसमें कोड modularity और software containers शामिल हैं
  audience_prerequisites:
    - "**दर्शक:** यह कोर्स उन शिक्षार्थियों के लिए डिज़ाइन किया गया है जो Nextflow में बिल्कुल नए हैं और मौजूदा pipelines चलाना चाहते हैं।"
    - "**कौशल:** कमांड लाइन, बुनियादी scripting अवधारणाओं और सामान्य फ़ाइल फॉर्मेट्स से कुछ परिचितता अपेक्षित है।"
    - "**डोमेन:** सभी अभ्यास डोमेन-अज्ञेयवादी हैं, इसलिए किसी पूर्व वैज्ञानिक ज्ञान की आवश्यकता नहीं है।"
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run reproducible और scalable डेटा विश्लेषण workflows चलाने का एक हैंड्स-ऑन परिचय है।**

व्यावहारिक उदाहरणों और निर्देशित अभ्यासों के माध्यम से, तुम Nextflow का उपयोग करने के मूल सिद्धांतों को सीखोगे, जिसमें pipelines कैसे execute करें, फ़ाइलों और software dependencies को कैसे मैनेज करें, execution को सहजता से parallelize कैसे करें, और विभिन्न कंप्यूटिंग वातावरणों में workflows कैसे चलाएं।

तुम Nextflow के साथ workflows चलाने के लिए कौशल और आत्मविश्वास प्राप्त करोगे।

<!-- additional_information -->

## कोर्स अवलोकन

### तुम क्या करोगे

यह कोर्स हैंड्स-ऑन है, जिसमें लक्ष्य-उन्मुख अभ्यास हैं जो धीरे-धीरे जानकारी प्रस्तुत करने के लिए संरचित हैं।

तुम एक Nextflow pipeline के कई संस्करण execute करोगे जो टेक्स्ट इनपुट को प्रोसेस करती है।
तुम एक साधारण संस्करण से शुरू करोगे जिसमें एक ही स्टेप है, और अंततः एक मल्टी-स्टेप संस्करण तक पहुंचोगे जो टेबुलर टेक्स्ट इनपुट की एक CSV फ़ाइल लेता है, कुछ transformation स्टेप्स चलाता है, और एक टेक्स्ट फ़ाइल आउटपुट करता है जिसमें transformed टेक्स्ट बोलते हुए एक character की ASCII तस्वीर होती है।

यह कोर्स pipelines चलाने पर केंद्रित है (कोर `nextflow run` कमांड के नाम पर)।
यदि तुम Nextflow pipelines विकसित करने का परिचय खोज रहे हो, तो [Hello Nextflow](../hello_nextflow/index.md) देखो।

### पाठ योजना

हमने इसे तीन भागों में विभाजित किया है जो प्रत्येक Nextflow में लिखी गई pipelines को चलाने और मैनेज करने के विशिष्ट पहलुओं पर ध्यान केंद्रित करेंगे।

| कोर्स अध्याय                                    | सारांश                                                                                                                   | अनुमानित अवधि |
| ----------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ | ------------- |
| [भाग 1: बुनियादी संचालन चलाएं](./01_basics.md)  | एक साधारण workflow को लॉन्च और मैनेज करना                                                                                | 30 मिनट       |
| [भाग 2: असली pipelines चलाएं](./02_pipeline.md) | जटिल इनपुट प्रोसेस करना, मल्टी-स्टेप workflows चलाना, containers का उपयोग करना और execution को सहजता से parallelize करना | 60 मिनट       |
| [भाग 3: Run configuration](./03_config.md)      | pipeline व्यवहार को कस्टमाइज़ करना और विभिन्न कम्प्यूटेशनल वातावरणों में उपयोग को optimize करना                          | 60 मिनट       |

इस कोर्स के अंत तक, तुम अपनी वैज्ञानिक कंप्यूटिंग आवश्यकताओं के लिए reproducible workflows चलाने की अपनी यात्रा में अगले कदम उठाने के लिए अच्छी तरह तैयार होगे।

कोर्स लेने के लिए तैयार हो?

[सीखना शुरू करें :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
