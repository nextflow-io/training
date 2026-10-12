---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - कमांड लाइन से Nextflow pipelines को लॉन्च और मैनेज करना
    - यह समझना कि चैनल और ऑपरेटर किस तरह कुशल मल्टी-इनपुट, मल्टी-स्टेप workflows को सक्षम बनाते हैं
    - software dependencies मैनेज करने और reproducibility सुनिश्चित करने के लिए containers का उपयोग करना
    - pipeline execution और आउटपुट को कॉन्फ़िगर करना
    - execution रिपोर्ट जनरेट करना, पिछले runs का इतिहास देखना, और पुरानी work directories को साफ़ करना
    - GitHub जैसे remote repositories से सीधे pipelines चलाना
  audience_prerequisites:
    - "**दर्शक:** यह कोर्स उन शिक्षार्थियों के लिए डिज़ाइन किया गया है जो Nextflow में बिल्कुल नए हैं और मौजूदा pipelines चलाना चाहते हैं।"
    - "**कौशल:** कमांड लाइन, बुनियादी scripting अवधारणाओं और सामान्य फ़ाइल फॉर्मेट्स से कुछ परिचितता अपेक्षित है।"
    - "**डोमेन:** सभी अभ्यास डोमेन-अज्ञेयवादी हैं, इसलिए किसी पूर्व वैज्ञानिक ज्ञान की आवश्यकता नहीं है।"
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run reproducible और scalable डेटा विश्लेषण workflows चलाने का एक हैंड्स-ऑन परिचय है।**

लक्ष्य-उन्मुख अभ्यासों की एक श्रृंखला के माध्यम से, तुम Nextflow pipelines को लॉन्च और मैनेज करने की मूल बातें सीखोगे, यह समझोगे कि चैनल और ऑपरेटर किस तरह कई इनपुट की parallel processing को सक्षम बनाते हैं, और software dependencies मैनेज करने के लिए containers का उपयोग करोगे।

तुम Nextflow के साथ workflows चलाने के लिए कौशल और आत्मविश्वास प्राप्त करोगे।

<!-- additional_information -->

## कोर्स अवलोकन

यह कोर्स हैंड्स-ऑन है, जिसमें लक्ष्य-उन्मुख अभ्यास हैं जो धीरे-धीरे जानकारी प्रस्तुत करने के लिए संरचित हैं।

तुम एक Nextflow pipeline के कई संस्करण execute करोगे जो टेक्स्ट इनपुट को प्रोसेस करती है, एक साधारण सिंगल-स्टेप संस्करण से शुरू करके एक मल्टी-स्टेप संस्करण तक पहुंचोगे जो इनपुट की एक CSV फ़ाइल लेता है, कुछ transformation स्टेप्स चलाता है, और एक containerized टूल द्वारा जनरेट की गई ASCII art वाली एक टेक्स्ट फ़ाइल आउटपुट करता है।

यह कोर्स pipelines चलाने पर केंद्रित है (कोर `nextflow run` कमांड के नाम पर)।
यदि तुम Nextflow pipelines विकसित करने का परिचय खोज रहे हो, तो [Hello Nextflow](../hello_nextflow/index.md) देखो।

!!! note "नोट"

    इस कोर्स का पिछला संस्करण खोज रहे हो? इसे इस पेज पर मौजूद संस्करण ने supersede कर दिया है, लेकिन यह training साइट के [3.6.1 release](https://training.nextflow.io/3.6.1/nextflow_run/) में अभी भी देखा जा सकता है।

### पाठ योजना

| कोर्स अध्याय                                                       | सारांश                                                                                            | अनुमानित अवधि |
| ------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------- | ------------- |
| [भाग 1: Nextflow चलाएं](./01_run_nextflow.md)                      | Nextflow pipelines को लॉन्च और मैनेज करना, और workflow mechanics की मूल बातें समझना               | 25 मिनट       |
| [भाग 2: Pipeline को कॉन्फ़िगर करें](./02_configure_pipeline.md)    | `nextflow.config` का उपयोग करके pipeline execution और आउटपुट को कॉन्फ़िगर करना                    | 20 मिनट       |
| [भाग 3: Workflow executions मैनेज करें](./03_manage_executions.md) | execution रिपोर्ट जनरेट करना, पिछले runs का इतिहास देखना, और पुरानी work directories को साफ़ करना | 10 मिनट       |
| [भाग 4: Remote pipelines चलाएं](./04_remote_repositories.md)       | GitHub से सीधे एक pipeline चलाना और उसे किसी specific revision पर pin करना                        | 10 मिनट       |

इस कोर्स के अंत तक, तुम अपनी वैज्ञानिक कंप्यूटिंग आवश्यकताओं के लिए reproducible workflows चलाने की अपनी यात्रा में अगले कदम उठाने के लिए अच्छी तरह तैयार होगे।

कोर्स लेने के लिए तैयार हो?

[सीखना शुरू करें :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
