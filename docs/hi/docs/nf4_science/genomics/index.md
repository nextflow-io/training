---
title: जीनोमिक्स के लिए Nextflow
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - एक single sample पर variant calling लागू करने के लिए एक linear workflow लिखना
    - Index files और reference genome resources जैसी accessory files को उचित तरीके से handle करना
    - Per-sample variant calling को parallelize करने के लिए Nextflow के dataflow paradigm का लाभ उठाना
    - Relevant channel operators का उपयोग करके multi-sample joint calling को implement करना
  audience_prerequisites:
    - "**दर्शक:** यह पाठ्यक्रम जीनोमिक्स और संबंधित क्षेत्रों के शोधकर्ताओं के लिए डिज़ाइन किया गया है जो data analysis pipelines विकसित या customize करना चाहते हैं।"
    - "**कौशल:** Command line, बुनियादी scripting अवधारणाओं, और सामान्य genomics file formats से कुछ परिचित होना मान लिया गया है।"
    - "**पूर्वापेक्षाएँ:** [Hello Nextflow](../../hello_nextflow/) में शामिल बुनियादी Nextflow अवधारणाएँ और टूलिंग।"
---

# जीनोमिक्स के लिए Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**एक व्यावहारिक पाठ्यक्रम जो Nextflow को वास्तविक दुनिया के जीनोमिक्स उपयोग के मामले में लागू करता है: GATK के साथ variant calling।**

यह पाठ्यक्रम [Hello Nextflow](../../hello_nextflow/) शुरुआती प्रशिक्षण पर आधारित है और यह प्रदर्शित करता है कि जीनोमिक्स डोमेन के विशिष्ट संदर्भ में Nextflow का उपयोग कैसे करें।
तुम [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit) के साथ एक variant calling pipeline लागू करोगे, जो high-throughput sequencing डेटा के विश्लेषण के लिए व्यापक रूप से उपयोग किया जाने वाला software package है।

<!-- additional_information -->

## पाठ्यक्रम का अवलोकन

यह पाठ्यक्रम व्यावहारिक है, जिसमें लक्ष्य-उन्मुख अभ्यास हैं जो जानकारी को क्रमिक रूप से प्रस्तुत करने के लिए संरचित हैं।

तुम पहले methodology को समझने के लिए terminal में manually variant calling tools चलाओगे, फिर धीरे-धीरे एक Nextflow pipeline बनाओगे जो विश्लेषण को स्वचालित और स्केल करती है।

### पाठ योजना

हमने इसे तीन भागों में विभाजित किया है जो प्रत्येक जीनोमिक्स उपयोग के मामले में Nextflow को लागू करने के विशिष्ट पहलुओं पर ध्यान केंद्रित करते हैं।

| पाठ्यक्रम अध्याय                                                         | सारांश                                                                                                                | अनुमानित अवधि |
| ------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------- | ------------- |
| [Part 1: Method overview](./01_method.md)                                | Variant calling methodology को समझना और tools को manually चलाना                                                       | 30 मिनट       |
| [Part 2: Per-sample variant calling](./02_per_sample_variant_calling.md) | एक pipeline बनाना जो BAM फ़ाइलों को index करती है और variants को call करती है, फिर कई नमूनों तक स्केल करना            | 60 मिनट       |
| [Part 3: Joint calling on a cohort](./03_joint_calling.md)               | Channel operators का उपयोग करके multi-sample joint genotyping जोड़ना ताकि per-sample outputs को aggregate किया जा सके | 45 मिनट       |

इस पाठ्यक्रम के अंत तक, तुम एक विशिष्ट जीनोमिक्स उपयोग के मामले में बुनियादी Nextflow अवधारणाओं और टूलिंग को लागू कर सकोगे।

पाठ्यक्रम लेने के लिए तैयार हो?

[शुरू करें :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
