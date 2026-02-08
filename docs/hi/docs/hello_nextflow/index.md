---
title: Hello Nextflow
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Nextflow workflows को लॉन्च और मैनेज करना
    - Nextflow द्वारा जनरेट किए गए आउटपुट (परिणाम) और लॉग फ़ाइलों को खोजना और समझना
    - बुनियादी समस्याओं का समाधान करना
    - मुख्य Nextflow कंपोनेंट्स से एक सरल मल्टी-स्टेप workflow बनाना
    - आवश्यक प्रकार के channel factories और operators को पहचानना और उन्हें एक सरल workflow में प्रभावी ढंग से उपयोग करना
    - HPC और cloud सहित सामान्य कंप्यूटिंग प्लेटफॉर्म पर चलाने के लिए pipeline execution को कॉन्फ़िगर करना
    - Reproducibility, portability और code re-use के लिए best practices लागू करना जो pipelines को FAIR बनाते हैं, जिसमें code modularity और software containers शामिल हैं
  audience_prerequisites:
    - "**दर्शक:** यह कोर्स उन learners के लिए डिज़ाइन किया गया है जो Nextflow में बिल्कुल नए हैं और अपने खुद के pipelines विकसित करना चाहते हैं।"
    - "**कौशल:** command line, बुनियादी scripting concepts और सामान्य file formats से कुछ परिचितता मानी गई है।"
    - "**डोमेन:** सभी exercises डोमेन-अज्ञेयवादी हैं, इसलिए कोई पूर्व वैज्ञानिक ज्ञान आवश्यक नहीं है।"
  videos_playlist: https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n
---

# Hello Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello Nextflow reproducible और scalable data analysis workflows बनाने का एक hands-on परिचय है।**

व्यावहारिक उदाहरणों और guided exercises के माध्यम से काम करते हुए, तुम Nextflow के साथ pipelines विकसित करने की मूल बातें सीखोगे, जिसमें processes को परिभाषित करना, उन्हें pipelines में जोड़ना, files और software dependencies को मैनेज करना, execution को आसानी से parallelize करना, और विभिन्न computing environments में workflows चलाना शामिल है।

तुम Nextflow के साथ अपने खुद के workflows विकसित करने और चलाने के कौशल और आत्मविश्वास लेकर जाओगे।

<!-- additional_information -->

## कोर्स अवलोकन

यह कोर्स hands-on होने के लिए डिज़ाइन किया गया है, जिसमें goal-oriented exercises को धीरे-धीरे जानकारी पेश करने के लिए संरचित किया गया है।

तुम एक सरल Nextflow pipeline विकसित करोगे जो कुछ text inputs लेता है, कुछ transformation steps चलाता है, और एक single text file आउटपुट करता है जिसमें एक character का ASCII picture है जो transformed text बोल रहा है।

### पाठ योजना

तुम्हें concepts और code से अभिभूत होने से बचाने के लिए, हमने इसे छह भागों में विभाजित किया है जो प्रत्येक Nextflow के साथ pipelines विकसित करने के विशिष्ट पहलुओं पर ध्यान केंद्रित करेंगे।

| कोर्स अध्याय                                        | सारांश                                                                                                    | अनुमानित अवधि |
| --------------------------------------------------- | --------------------------------------------------------------------------------------------------------- | ------------- |
| [भाग 1: Hello World](./01_hello_world.md)           | Nextflow workflow को assemble और run करने में शामिल बुनियादी components और principles                     | 30 मिनट       |
| [भाग 2: Hello Channels](./02_hello_channels.md)     | Inputs को process करने और execution को आसानी से parallelize करने के लिए channels और operators का उपयोग    | 45 मिनट       |
| [भाग 3: Hello Workflow](./03_hello_workflow.md)     | Multiple steps को एक साथ chain करने और steps के बीच data transfer को handle करने के लिए channels का उपयोग | 60 मिनट       |
| [भाग 4: Hello Modules](./04_hello_modules.md)       | Reusability बढ़ाने और maintenance burden कम करने के लिए code modularity principles लागू करना              | 20 मिनट       |
| [भाग 5: Hello Containers](./05_hello_containers.md) | Software dependencies को मैनेज करने और reproducibility बढ़ाने के लिए containers का उपयोग करना             | 60 मिनट       |
| [भाग 6: Hello Config](./06_hello_config.md)         | Pipeline behavior को customize करना और विभिन्न computational environments में usage को optimize करना      | 60 मिनट       |

इस कोर्स के अंत तक, तुम अपनी scientific computing आवश्यकताओं के लिए reproducible workflows विकसित करने की अपनी यात्रा में अगले कदमों से निपटने के लिए अच्छी तरह से तैयार होगे।

कोर्स लेने के लिए तैयार हो?

[शुरू करें :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
