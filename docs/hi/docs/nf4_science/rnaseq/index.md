# Nextflow for RNAseq

यह प्रशिक्षण कोर्स ट्रांसक्रिप्टोमिक्स और संबंधित क्षेत्रों के शोधकर्ताओं के लिए है जो डेटा विश्लेषण पाइपलाइन विकसित करने या कस्टमाइज़ करने में रुचि रखते हैं।
यह [Hello Nextflow](../../hello_nextflow/) शुरुआती प्रशिक्षण पर आधारित है और दिखाता है कि bulk RNAseq विश्लेषण के विशिष्ट संदर्भ में Nextflow का उपयोग कैसे करें।

विशेष रूप से, यह कोर्स दिखाता है कि एक सरल bulk RNAseq प्रोसेसिंग पाइपलाइन को कैसे लागू करें जो adapter sequences को trim करती है, reads को genome reference के साथ align करती है और कई चरणों में quality control (QC) करती है।

चलो शुरू करते हैं! प्रशिक्षण वातावरण लॉन्च करने के लिए नीचे "Open in GitHub Codespaces" बटन पर क्लिक करो (अधिमानतः एक अलग टैब में), फिर जब यह लोड हो रहा हो तब पढ़ते रहो।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## सीखने के उद्देश्य

इस कोर्स के माध्यम से काम करके, तुम सीखोगे कि एक विशिष्ट RNAseq उपयोग के मामले में बुनियादी Nextflow अवधारणाओं और टूलिंग को कैसे लागू करें।

इस वर्कशॉप के अंत तक तुम सक्षम होगे:

- बुनियादी RNAseq प्रोसेसिंग और QC विधियों को लागू करने के लिए एक linear workflow लिखना
- FASTQ और reference genome resources जैसी domain-specific फ़ाइलों को उचित रूप से संभालना
- Single-end और paired-end sequencing डेटा को संभालना
- Per-sample RNAseq प्रोसेसिंग को समानांतर करने के लिए Nextflow के dataflow paradigm का लाभ उठाना
- प्रासंगिक channel operators का उपयोग करके कई चरणों और नमूनों में QC रिपोर्ट को एकत्रित करना

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## पूर्वापेक्षाएँ

कोर्स निम्नलिखित के साथ कुछ न्यूनतम परिचितता मानता है:

- इस वैज्ञानिक क्षेत्र में आमतौर पर उपयोग किए जाने वाले टूल और फ़ाइल फ़ॉर्मेट
- कमांड लाइन के साथ अनुभव
- [Hello Nextflow](../../hello_nextflow/) शुरुआती प्रशिक्षण में शामिल बुनियादी Nextflow अवधारणाएँ और टूलिंग

तकनीकी आवश्यकताओं और वातावरण सेटअप के लिए, [Environment Setup](../../envsetup/) मिनी-कोर्स देखो।
