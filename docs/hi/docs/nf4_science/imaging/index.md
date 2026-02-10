---
title: इमेजिंग के लिए Nextflow run
hide:
  - toc
---

# इमेजिंग के लिए Nextflow run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

यह प्रशिक्षण पाठ्यक्रम इमेजिंग और स्थानिक जीवविज्ञान के शोधकर्ताओं के लिए है जो डेटा विश्लेषण pipelines को चलाने और अनुकूलित करने में रुचि रखते हैं।
यह [nf-core/molkart](https://nf-co.re/molkart) का उपयोग करके workflows को चलाने, व्यवस्थित करने और कॉन्फ़िगर करने से संबंधित बुनियादी Nextflow अवधारणाएँ सिखाता है, जो Molecular Cartography स्थानिक transcriptomics डेटा को प्रोसेस करने के लिए एक pipeline है।
यहाँ आप जो कौशल सीखते हैं वे किसी भी Nextflow या nf-core pipeline में लागू होते हैं।

चलिए शुरू करते हैं! प्रशिक्षण वातावरण लॉन्च करने के लिए नीचे दिए गए "Open in GitHub Codespaces" बटन पर क्लिक करें (अधिमानतः एक अलग टैब में), फिर जब यह लोड हो रहा हो तो आगे पढ़ें।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## सीखने के उद्देश्य

इस पाठ्यक्रम के माध्यम से काम करके, आप सीखेंगे कि इमेजिंग विश्लेषण pipelines को चलाने के लिए बुनियादी Nextflow अवधारणाओं और टूलिंग को कैसे लागू किया जाए।

इस कार्यशाला के अंत तक आप सक्षम होंगे:

- Nextflow workflow को स्थानीय रूप से लॉन्च करना और निष्पादन की निगरानी करना
- Nextflow द्वारा उत्पन्न outputs (परिणाम) और log फ़ाइलों को ढूंढना और समझना
- test डेटा और custom inputs के साथ nf-core pipeline चलाना
- profiles और parameter फ़ाइलों का उपयोग करके pipeline निष्पादन को कॉन्फ़िगर करना
- samplesheets और command-line parameters का उपयोग करके inputs को प्रबंधित करना

## दर्शक और पूर्वापेक्षाएँ

यह पाठ्यक्रम निम्नलिखित से कुछ न्यूनतम परिचितता मानता है:

- command line के साथ अनुभव
- इमेजिंग फ़ाइल formats (TIFF images, tabular data) से बुनियादी परिचितता

तकनीकी आवश्यकताओं और वातावरण सेटअप के लिए, [Environment Setup](../../envsetup/) मिनी-पाठ्यक्रम देखें।
