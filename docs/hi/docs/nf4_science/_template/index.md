---
title: Nextflow for {DOMAIN}
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Write a linear workflow to apply {METHOD} to a single sample
    - Handle accessory files such as {ACCESSORY_FILES} appropriately
    - Leverage Nextflow's dataflow paradigm to parallelize per-sample processing
    - Implement multi-sample aggregation using relevant channel operators
  audience_prerequisites:
    - "**Audience:** This course is designed for researchers in {DOMAIN} and related fields who want to develop or customize data analysis pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts, and common {DOMAIN} file formats is assumed."
    - "**Prerequisites:** Foundational Nextflow concepts and tooling covered in [Hello Nextflow](../../hello_nextflow/)."
---

---
title: {DOMAIN} के लिए Nextflow
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - एक सिंगल सैंपल पर {METHOD} लागू करने के लिए एक लीनियर वर्कफ़्लो लिखें
    - {ACCESSORY_FILES} जैसी सहायक फ़ाइलों को उचित तरीके से हैंडल करें
    - प्रति-सैंपल प्रोसेसिंग को समानांतर बनाने के लिए Nextflow के डेटाफ़्लो पैराडाइम का लाभ उठाएं
    - प्रासंगिक चैनल ऑपरेटरों का उपयोग करके मल्टी-सैंपल एग्रीगेशन लागू करें
  audience_prerequisites:
    - "**दर्शक:** यह कोर्स {DOMAIN} और संबंधित क्षेत्रों के शोधकर्ताओं के लिए डिज़ाइन किया गया है जो डेटा विश्लेषण पाइपलाइनों को विकसित या कस्टमाइज़ करना चाहते हैं।"
    - "**कौशल:** कमांड लाइन, बुनियादी स्क्रिप्टिंग अवधारणाओं और सामान्य {DOMAIN} फ़ाइल फ़ॉर्मेट्स से कुछ परिचितता मानी गई है।"
    - "**पूर्वापेक्षाएं:** [Hello Nextflow](../../hello_nextflow/) में शामिल बुनियादी Nextflow अवधारणाएं और टूलिंग।"
---

# {DOMAIN} के लिए Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**एक व्यावहारिक कोर्स जो Nextflow को वास्तविक दुनिया के {DOMAIN} उपयोग के मामले पर लागू करता है: {METHOD_SHORT_DESCRIPTION}।**

यह कोर्स [Hello Nextflow](../../hello_nextflow/) शुरुआती प्रशिक्षण पर आधारित है और दिखाता है कि {DOMAIN} डोमेन के विशिष्ट संदर्भ में Nextflow का उपयोग कैसे करें।
तुम [{TOOL_A}]({TOOL_A_URL}) और [{TOOL_B}]({TOOL_B_URL}) के साथ एक {METHOD} पाइपलाइन लागू करोगे।

<!-- additional_information -->

## कोर्स अवलोकन

यह कोर्स व्यावहारिक है, जिसमें लक्ष्य-उन्मुख अभ्यास हैं जो जानकारी को धीरे-धीरे पेश करने के लिए संरचित हैं।

तुम पद्धति को समझने के लिए टर्मिनल में विश्लेषण टूल्स को मैन्युअली चलाकर शुरू करोगे, फिर धीरे-धीरे एक Nextflow पाइपलाइन बनाओगे जो विश्लेषण को स्वचालित और स्केल करती है।

### पाठ योजना

हमने इसे तीन भागों में विभाजित किया है जो प्रत्येक {DOMAIN} उपयोग के मामले में Nextflow को लागू करने के विशिष्ट पहलुओं पर ध्यान केंद्रित करते हैं।

| कोर्स अध्याय                                              | सारांश                                                                                                      | अनुमानित अवधि |
| --------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- | ------------- |
| [भाग 1: विधि अवलोकन](./01_method.md)                     | {METHOD} पद्धति को समझना और टूल्स को मैन्युअली चलाना                                                       | 30 मिनट       |
| [भाग 2: सिंगल-सैंपल प्रोसेसिंग](./02_single_sample.md)   | एक पाइपलाइन बनाना जो {PART2_SUMMARY}, फिर कई सैंपलों तक स्केल करना                                         | 60 मिनट       |
| [भाग 3: मल्टी-सैंपल एग्रीगेशन](./03_multi_sample.md)     | प्रति-सैंपल आउटपुट को एग्रीगेट करने के लिए चैनल ऑपरेटरों का उपयोग करके मल्टी-सैंपल {AGGREGATION_SUMMARY} जोड़ना | 45 मिनट       |

इस कोर्स के अंत तक, तुम एक विशिष्ट {DOMAIN} उपयोग के मामले में बुनियादी Nextflow अवधारणाओं और टूलिंग को लागू करने में सक्षम होगे।

कोर्स शुरू करने के लिए तैयार हो?

[शुरू करें :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
