---
title: वातावरण विकल्प
description: Nextflow प्रशिक्षणों के लिए अपना वातावरण सेट करने के विकल्प
hide:
  - toc
  - footer
---

# वातावरण विकल्प

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

हमारा लक्ष्य एक सुसंगत और पूर्णतः परीक्षित वातावरण प्रदान करना है जो शिक्षार्थियों को software प्रबंधन पर समय और प्रयास खर्च किए बिना Nextflow सीखने पर ध्यान केंद्रित करने की अनुमति देता है।
इसके लिए, हमने एक containerized वातावरण विकसित किया है जिसमें हमारे सभी courses में काम करने के लिए सभी आवश्यक software, code files और उदाहरण data शामिल हैं।

यह containerized वातावरण Github Codespaces पर या VS Code में Devcontainers extension के साथ स्थानीय रूप से बॉक्स से बाहर चलाया जा सकता है।

<div class="grid cards" markdown>

- :material-cloud-outline:{ .lg .middle } **Github Codespaces**

  ***

  GitHub Codespaces एक web-आधारित सेवा है जो हमें प्रशिक्षण के लिए पूर्व-निर्मित वातावरण प्रदान करने की अनुमति देती है, जिसमें सभी tools और data शामिल हैं, cloud में virtual machines द्वारा समर्थित। यह Github account वाले किसी भी व्यक्ति के लिए मुफ्त में उपलब्ध है।

  [Github Codespaces का उपयोग करें:material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

- :material-laptop:{ .lg .middle } **Local Devcontainers**

  ***

  Devcontainers के साथ VS Code स्थानीय रूप से चलने वाला containerized development वातावरण प्रदान करता है जिसमें सभी प्रशिक्षण tools पूर्व-configured हैं। यह Codespaces के समान पूर्व-निर्मित वातावरण प्रदान करता है लेकिन पूरी तरह से तुम्हारे स्थानीय hardware पर चलता है।

  [स्थानीय रूप से Devcontainers का उपयोग करें :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## मैन्युअल इंस्टॉलेशन के लिए निर्देश

अगर ऊपर दिए गए विकल्पों में से कोई भी तुम्हारी ज़रूरतों के अनुकूल नहीं है, तो तुम software dependencies को मैन्युअल रूप से install करके और training repository को clone करके अपने स्थानीय system पर इस वातावरण की प्रतिकृति बना सकते हो।

[मैन्युअल इंस्टॉलेशन :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Gitpod का बहिष्करण"

    Nextflow Training फरवरी 2025 तक [Gitpod](https://gitpod.io) का उपयोग करता था।
    हालाँकि, Gitpod के निर्माताओं ने [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex) system के पक्ष में मुफ्त कार्यक्षमता को retire करने का निर्णय लिया।
    इसी कारण से, हमने GitHub Codespaces का उपयोग करना शुरू किया, जो बिना किसी पूर्व setup के one-click developer वातावरण भी प्रदान करता है।

    तुम्हारे Gitpod में साइन अप करने के समय और उनके सेवा को retire करने के समय के आधार पर, तुम अभी भी उनके पुराने cloud IDE में प्रशिक्षण शुरू करने में सक्षम हो सकते हो, हालाँकि हम आगे विश्वसनीय पहुँच की गारंटी नहीं दे सकते:
    [Gitpod में खोलें](https://gitpod.io/#https://github.com/nextflow-io/training)।
