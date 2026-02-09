# पर्यावरण विकल्प

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

हमारा लक्ष्य एक सुसंगत और पूरी तरह से परीक्षित पर्यावरण प्रदान करना है जो सीखने वालों को सॉफ़्टवेयर प्रबंधन में समय और प्रयास खर्च किए बिना Nextflow सीखने पर ध्यान केंद्रित करने की अनुमति देता है।
इस उद्देश्य के लिए, हमने एक containerized पर्यावरण विकसित किया है जिसमें हमारे सभी कोर्सों को पूरा करने के लिए आवश्यक सभी सॉफ़्टवेयर, कोड फ़ाइलें और उदाहरण डेटा शामिल हैं।

यह containerized पर्यावरण Github Codespaces पर या VS Code में Devcontainers एक्सटेंशन के साथ लोकल रूप से चलाया जा सकता है।

<div class="grid cards" markdown>

-   :material-cloud-outline:{ .lg .middle } __Github Codespaces__

    ---

    GitHub Codespaces एक वेब-आधारित सेवा है जो हमें प्रशिक्षण के लिए एक पूर्व-निर्मित पर्यावरण प्रदान करने की अनुमति देती है, जिसमें सभी टूल्स और डेटा शामिल हैं, जो क्लाउड में वर्चुअल मशीनों द्वारा समर्थित है। यह Github अकाउंट वाले किसी भी व्यक्ति के लिए मुफ़्त में उपलब्ध है।

    [Github Codespaces का उपयोग करें:material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

-   :material-laptop:{ .lg .middle } __लोकल Devcontainers__

    ---

    VS Code के साथ Devcontainers एक लोकल रूप से चलने वाला containerized विकास पर्यावरण प्रदान करता है जिसमें सभी प्रशिक्षण टूल्स पूर्व-कॉन्फ़िगर किए गए हैं। यह Codespaces जैसा ही पूर्व-निर्मित पर्यावरण प्रदान करता है लेकिन पूरी तरह से तुम्हारे लोकल हार्डवेयर पर चलता है।

    [Devcontainers को लोकल रूप से उपयोग करें :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## मैनुअल इंस्टॉलेशन के लिए निर्देश

यदि उपरोक्त विकल्पों में से कोई भी तुम्हारी आवश्यकताओं के अनुरूप नहीं है, तो तुम सॉफ़्टवेयर निर्भरताओं को मैनुअल रूप से इंस्टॉल करके और प्रशिक्षण रिपॉजिटरी को क्लोन करके अपने लोकल सिस्टम पर इस पर्यावरण को दोहरा सकते हो।

[मैनुअल इंस्टॉलेशन :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Gitpod का बंद होना"

    Nextflow Training फरवरी 2025 तक [Gitpod](https://gitpod.io) का उपयोग करता था।
    हालांकि, Gitpod के निर्माताओं ने [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex) सिस्टम के पक्ष में मुफ़्त कार्यक्षमता को बंद करने का निर्णय लिया।
    इस कारण से, हमने GitHub Codespaces का उपयोग करना शुरू कर दिया, जो बिना किसी पूर्व सेटअप के एक-क्लिक डेवलपर पर्यावरण भी प्रदान करता है।

    तुमने Gitpod में कब साइन अप किया और वे सेवा को कब बंद करते हैं, इसके आधार पर तुम अभी भी उनके पुराने क्लाउड IDE में प्रशिक्षण लॉन्च कर सकते हो, हालांकि हम आगे विश्वसनीय पहुंच की गारंटी नहीं दे सकते:
    [Gitpod में खोलें](https://gitpod.io/#https://github.com/nextflow-io/training)।
