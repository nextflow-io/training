---
title: RNAseq के लिए Nextflow
hide:
  - toc
---

# RNAseq के लिए Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

यह प्रशिक्षण कोर्स ट्रांसक्रिप्टोमिक्स और संबंधित क्षेत्रों के शोधकर्ताओं के लिए है जो डेटा विश्लेषण पाइपलाइन विकसित करने या कस्टमाइज़ करने में रुचि रखते हैं।
यह [Hello Nextflow](../../hello_nextflow/) शुरुआती प्रशिक्षण पर आधारित है और दर्शाता है कि बल्क RNAseq विश्लेषण के विशिष्ट संदर्भ में Nextflow का उपयोग कैसे करें।

विशेष रूप से, यह कोर्स दिखाता है कि एडॉप्टर अनुक्रमों को ट्रिम करने, रीड्स को जीनोम संदर्भ के साथ संरेखित करने और कई चरणों में गुणवत्ता नियंत्रण (QC) करने के लिए एक सरल बल्क RNAseq प्रोसेसिंग पाइपलाइन कैसे लागू करें।

चलिए शुरू करते हैं! प्रशिक्षण वातावरण लॉन्च करने के लिए नीचे दिए गए "Open in GitHub Codespaces" बटन पर क्लिक करें (अधिमानतः एक अलग टैब में), फिर लोड होने के दौरान आगे पढ़ें।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## सीखने के उद्देश्य

इस कोर्स के माध्यम से काम करके, आप सीखेंगे कि एक विशिष्ट RNAseq उपयोग के मामले में बुनियादी Nextflow अवधारणाओं और टूलिंग को कैसे लागू करें।

इस वर्कशॉप के अंत तक आप सक्षम होंगे:

- बुनियादी RNAseq प्रोसेसिंग और QC विधियों को लागू करने के लिए एक रैखिक वर्कफ़्लो लिखना
- FASTQ जैसी डोमेन-विशिष्ट फ़ाइलों और जीनोम संदर्भ संसाधनों को उचित रूप से संभालना
- सिंगल-एंड और पेयर्ड-एंड सीक्वेंसिंग डेटा को संभालना
- प्रति-नमूना RNAseq प्रोसेसिंग को समानांतरित करने के लिए Nextflow के डेटाफ़्लो प्रतिमान का लाभ उठाना
- प्रासंगिक channel ऑपरेटरों का उपयोग करके कई चरणों और नमूनों में QC रिपोर्ट एकत्रित करना

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## पूर्वापेक्षाएँ

कोर्स निम्नलिखित के साथ न्यूनतम परिचितता मानता है:

- इस वैज्ञानिक डोमेन में आमतौर पर उपयोग किए जाने वाले टूल और फ़ाइल फॉर्मेट
- कमांड लाइन के साथ अनुभव
- [Hello Nextflow](../../hello_nextflow/) शुरुआती प्रशिक्षण में कवर की गई बुनियादी Nextflow अवधारणाएँ और टूलिंग

तकनीकी आवश्यकताओं और वातावरण सेटअप के लिए, [Environment Setup](../../envsetup/) मिनी-कोर्स देखें।
