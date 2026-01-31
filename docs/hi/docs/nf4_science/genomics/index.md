# जीनोमिक्स के लिए Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

यह प्रशिक्षण पाठ्यक्रम जीनोमिक्स और संबंधित क्षेत्रों के शोधकर्ताओं के लिए है जो डेटा विश्लेषण pipelines विकसित करने या उन्हें अनुकूलित करने में रुचि रखते हैं।
यह [Hello Nextflow](../../hello_nextflow/) शुरुआती प्रशिक्षण पर आधारित है और यह प्रदर्शित करता है कि जीनोमिक्स डोमेन के विशिष्ट संदर्भ में Nextflow का उपयोग कैसे करें।

विशेष रूप से, यह पाठ्यक्रम [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit) के साथ एक साधारण variant calling pipeline को लागू करने का तरीका प्रदर्शित करता है, जो high-throughput sequencing डेटा के विश्लेषण के लिए व्यापक रूप से उपयोग किया जाने वाला software package है।

आइए शुरू करें! प्रशिक्षण वातावरण लॉन्च करने के लिए नीचे "Open in GitHub Codespaces" बटन पर क्लिक करें (अधिमानतः एक अलग टैब में), फिर जब वह लोड हो रहा हो तब आगे पढ़ें।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## सीखने के उद्देश्य

इस पाठ्यक्रम के माध्यम से काम करके, आप सीखेंगे कि एक विशिष्ट जीनोमिक्स उपयोग के मामले में बुनियादी Nextflow अवधारणाओं और टूलिंग को कैसे लागू करें।

इस workshop के अंत तक आप सक्षम होंगे:

- एक single नमूने पर variant calling लागू करने के लिए एक linear workflow लिखना
- index फ़ाइलें और reference genome संसाधनों जैसी सहायक फ़ाइलों को उचित रूप से संभालना
- per-sample variant calling को समानांतर करने के लिए Nextflow के dataflow paradigm का लाभ उठाना
- प्रासंगिक channel operators का उपयोग करके multi-sample variant calling को लागू करना
- per-step और end-to-end pipeline tests को लागू करना जो जीनोमिक्स-विशिष्ट विशेषताओं को उचित रूप से संभालते हैं

<!-- TODO for future expansion: add metadata/samplesheet handling -->

## पूर्वापेक्षाएँ

पाठ्यक्रम निम्नलिखित के साथ न्यूनतम परिचय मानता है:

- इस वैज्ञानिक डोमेन में सामान्यतः उपयोग किए जाने वाले tools और file formats
- command line के साथ अनुभव
- [Hello Nextflow](../../hello_nextflow/) शुरुआती प्रशिक्षण में शामिल बुनियादी Nextflow अवधारणाओं और टूलिंग

तकनीकी आवश्यकताओं और वातावरण सेटअप के लिए, [Environment Setup](../../envsetup/) mini-course देखें।
