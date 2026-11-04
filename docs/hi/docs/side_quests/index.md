---
title: साइड क्वेस्ट
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - एक उत्पादक Nextflow डेवलपमेंट वातावरण सेट अप और कॉन्फ़िगर करें
    - जटिल डेटा ट्रांसफ़ॉर्मेशन के लिए उन्नत स्क्रिप्टिंग पैटर्न लागू करें
    - मल्टी-स्टेप वर्कफ़्लो में मेटाडेटा को हैंडल और प्रोपेगेट करें
    - समानांतर और क्रमिक प्रोसेसिंग के लिए डेटा चैनल को विभाजित और समूहित करें
    - nf-test का उपयोग करके Nextflow वर्कफ़्लो का परीक्षण करें
    - { "पुन": "उपयोग योग्य नामित वर्कफ़्लो मॉड्यूल से जटिल पाइपलाइन बनाएं" }
    - Nextflow फ़ाइल ऑपरेशन का उपयोग करके फ़ाइलों के साथ कुशलतापूर्वक काम करें
    - सामान्य वर्कफ़्लो समस्याओं को व्यवस्थित रूप से डीबग करें
    - Nextflow plugins का उपयोग और निर्माण करें
  audience_prerequisites:
    - "**दर्शक:** यह संग्रह उन शिक्षार्थियों के लिए बनाया गया है जिन्होंने Hello Nextflow शुरुआती कोर्स पूरा कर लिया है और विशिष्ट विषयों में गहराई से जाना चाहते हैं।"
    - "**कौशल:** कमांड लाइन का अनुभव और बुनियादी Nextflow अवधारणाओं और टूलिंग से परिचितता मानी जाती है।"
    - "**कोर्स:** [Hello Nextflow](../hello_nextflow/index.md) या समकक्ष पूरा किया होना चाहिए।"
---

# साइड क्वेस्ट

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**साइड क्वेस्ट स्वतंत्र प्रशिक्षण मिनी-कोर्स हैं जो विशिष्ट Nextflow विषयों में गहराई से जाते हैं।**

प्रत्येक साइड क्वेस्ट को स्वतंत्र रूप से, किसी भी क्रम में किया जा सकता है, ताकि तुम उन विशिष्ट क्षेत्रों में अपने कौशल को बढ़ा सको जो तुम्हारे अपने प्रोजेक्ट के लिए सबसे ज़रूरी हैं।
मिलकर, ये उन कौशलों को कवर करते हैं जो तुम्हें सरल पाइपलाइन से प्रोडक्शन-रेडी वर्कफ़्लो तक जाने के लिए चाहिए।

<!-- additional_information -->

### विषय के अनुसार ब्राउज़ करें

यदि यह तुम्हारा हमारे किसी प्रशिक्षण में पहली बार है, तो प्रशिक्षण वातावरण और सामग्री के अवलोकन के लिए [शुरू करना](./orientation.md) पेज से शुरू करो।
अन्यथा, सीधे उस क्वेस्ट में जाने के लिए स्वतंत्र महसूस करो जो तुम्हें बुला रही है।

<table>
<thead>
<tr><th>साइड क्वेस्ट</th><th>सारांश</th><th>अनुमानित समय</th></tr>
</thead>
<tbody>
<tr><th colspan="3">डेवलपर टूल्स &amp; ट्रिक्स</th></tr>
<tr><td><a href="./dev_environment/">Development Environment</a></td><td>एक उत्पादक लोकल Nextflow डेवलपमेंट वातावरण सेट अप और कॉन्फ़िगर करें</td><td>45 मिनट</td></tr>
<tr><td><a href="./debugging/">Troubleshooting Workflows</a></td><td>सामान्य वर्कफ़्लो त्रुटियों की पहचान और सुधार</td><td>1 घंटा</td></tr>
<tr><td><a href="./essential_scripting_patterns/">Essential Scripting Patterns</a></td><td>सामान्य वर्कफ़्लो चुनौतियों के लिए उन्नत स्क्रिप्टिंग तकनीकें</td><td>90 मिनट</td></tr>
<tr><th colspan="3">डेटाफ़्लो में गहरी डुबकी</th></tr>
<tr><td><a href="./working_with_files/">File Input Processing</a></td><td>फ़ाइल हैंडलिंग, पाथ ऑपरेशन और आउटपुट व्यवस्थित करना</td><td>45 मिनट</td></tr>
<tr><td><a href="./metadata/">Metadata and Meta Maps</a></td><td>नमूना जानकारी ट्रैक और प्रोपेगेट करने के लिए मेटाडेटा मैप का उपयोग</td><td>45 मिनट</td></tr>
<tr><td><a href="./splitting_and_grouping/">Splitting and Grouping</a></td><td>डेटा चैनल को विभाजित और पुनः समूहित करने की तकनीकें</td><td>45 मिनट</td></tr>
<tr><th colspan="3">मॉड्यूलर आर्किटेक्चर व्यवहार में</th></tr>
<tr><td><a href="./workflows_of_workflows/">Workflows of Workflows</a></td><td>पुन: उपयोग योग्य नामित वर्कफ़्लो मॉड्यूल से जटिल पाइपलाइन बनाना</td><td>30 मिनट</td></tr>
<tr><th colspan="3">Nextflow एक्सटेंडेड यूनिवर्स</th></tr>
<tr><td><a href="./nf_test/">Testing with nf-test</a></td><td>Nextflow वर्कफ़्लो के लिए परीक्षण लिखना और चलाना</td><td>1 घंटा</td></tr>
<tr><td><a href="./plugin_development/">Plugin Development</a></td><td>Nextflow plugins का उपयोग और निर्माण</td><td>3 घंटे</td></tr>
</tbody>
</table>

[शुरू करें :material-arrow-right:](orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
