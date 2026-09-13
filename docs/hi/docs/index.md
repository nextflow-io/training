---
title: होम
description: Nextflow कम्युनिटी ट्रेनिंग पोर्टल में आपका स्वागत है!
hide:
  - toc
  - footer
---

# Nextflow प्रशिक्षण

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __स्व-सेवा कोर्स__

    ---

    **Nextflow कम्युनिटी ट्रेनिंग पोर्टल में आपका स्वागत है!**

    नीचे दिए गए कोर्स अपनी गति से पूरे करो, हमारे वेब-आधारित वातावरण में या अपने खुद के वातावरण में।
    हर कोर्स व्यावहारिक है, जिसमें लक्ष्य-उन्मुख अभ्यास हैं जिन्हें तुम स्वतंत्र रूप से पूरा कर सकते हो।

    [कोर्स देखो :material-arrow-down:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-account-group-outline:{ .lg .middle } __प्रशिक्षण इवेंट__

    ---

    **स्व-सेवा से परे कुछ ढूंढ रहे हो?**

    संरचित प्रशिक्षण इवेंट, अपना खुद का प्रशिक्षण चलाने के लिए मार्गदर्शन, और हमारी ओपन-सोर्स लाइसेंस और योगदान नीति देखो।

    [प्रशिक्षण इवेंट देखो :material-arrow-right:](training_events.md){ .md-button .md-button--secondary .mt-1 }

</div>

!!! note "AI-सहायता प्राप्त अनुवाद"

    यह अनुवाद कृत्रिम बुद्धिमत्ता का उपयोग करके बनाया गया था और मानव अनुवादकों द्वारा समीक्षित किया गया था।
    हम आपकी प्रतिक्रिया और सुधार के सुझावों का स्वागत करते हैं।
    अधिक जानकारी के लिए हमारी [अनुवाद मार्गदर्शिका](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) देखें।

## Nextflow प्रशिक्षण कोर्स की सूची

<div class="grid cards" markdown>

-   :material-account:{ .lg .middle } __उपयोगकर्ताओं के लिए__

    ---

    ### :material-play-circle:{.nextflow-primary} पाइपलाइन चलाओ {.mt-1}

    बिना कोई कोड लिखे मौजूदा पाइपलाइन चलाना सीखो।

    ??? courses "**Nextflow Run:** Nextflow के साथ पाइपलाइन चलाओ"

        Nextflow पाइपलाइन चलाने का एक तेज़ परिचय जिसके लिए कोड समझने की ज़रूरत नहीं है। इसमें पाइपलाइन लॉन्च करना, आउटपुट प्राप्त करना, कंटेनर का उपयोग करना, और बुनियादी स्तर पर execution कॉन्फ़िगर करना शामिल है।

        [प्रशिक्षण देखो :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ??? courses "**Use nf-core:** कम्युनिटी-क्यूरेटेड पाइपलाइन ढूंढो और चलाओ"

        nf-core कम्युनिटी प्रोजेक्ट से पाइपलाइन ढूंढने, चलाने और कॉन्फ़िगर करने का एक तेज़ परिचय, एक न्यूनतम डेमो पाइपलाइन से शुरू करके प्रोडक्शन-स्केल विश्लेषण पाइपलाइन तक।

        [प्रशिक्षण देखो :material-arrow-right:](nfcore_use/index.md){ .md-button .md-button--secondary }

    ??? courses "**Scale with Seqera:** बड़े पैमाने पर पाइपलाइन लॉन्च और मॉनिटर करो"

        Seqera Platform के साथ Nextflow पाइपलाइन लॉन्च और मॉनिटर करने का एक व्यावहारिक परिचय, वेब इंटरफ़ेस और कमांड लाइन दोनों से।

        [प्रशिक्षण देखो :material-arrow-right:](seqera_scale/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-tune:{.nextflow-primary} execution प्रबंधित करो {.mt-1}

    पाइपलाइन execution को प्रभावी ढंग से प्रबंधित करना सीखो।

    ??? courses "**Execution Config:** पाइपलाइन को प्रो की तरह कॉन्फ़िगर करो"

        Nextflow पाइपलाइन execution कॉन्फ़िगर करने का एक व्यावहारिक परिचय: अलग-अलग कंप्यूट वातावरण के अनुसार ढलना, संसाधन आवंटन और retries नियंत्रित करना, और पूर्व-निर्धारित कॉन्फ़िगरेशन प्रोफ़ाइल के बीच स्विच करना।

        [प्रशिक्षण देखो :material-arrow-right:](execution_config/index.md){ .md-button .md-button--secondary }

    !!! info compact "और विषय आने वाले हैं"

        परफॉर्मेंस ट्यूनिंग, HPC/क्लाउड execution, और अन्य विषय इस सेक्शन के लिए योजनाबद्ध हैं।
        हमारे [छोटे इंटरेस्ट पोल](https://seqera.typeform.com/to/JCs91e8v) में वोट करो कि आगे क्या कवर किया जाए।

-   :material-code-tags:{ .lg .middle } __डेवलपर्स के लिए__

    ---

    ### :material-wrench:{.nextflow-primary} पाइपलाइन लिखो {.mt-1}

    अपनी खुद की Nextflow पाइपलाइन विकसित करना सीखो।

    ??? courses "**Hello Nextflow:** शुरू से अपनी खुद की पाइपलाइन विकसित करो"

        यह कोर्स Nextflow भाषा के मुख्य घटकों को इतने विस्तार से कवर करता है कि सरल लेकिन पूरी तरह कार्यात्मक पाइपलाइन विकसित की जा सकें, साथ ही पाइपलाइन डिज़ाइन, विकास और कॉन्फ़िगरेशन प्रथाओं के प्रमुख तत्व भी।

        [प्रशिक्षण देखो :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Build with nf-core:** nf-core टूल और नियमों का उपयोग करो"

        उन Nextflow डेवलपर्स के लिए जो [nf-core](https://nf-co.re/) अनुपालक पाइपलाइन विकसित करना सीखना चाहते हैं।
        यह कोर्स nf-core पाइपलाइन की संरचना को इतने विस्तार से कवर करता है कि सरल लेकिन पूरी तरह कार्यात्मक पाइपलाइन विकसित की जा सकें जो nf-core टेम्पलेट और विकास की सर्वोत्तम प्रथाओं का लाभ उठाएं, साथ ही मौजूदा nf-core मॉड्यूल का उपयोग भी करें।

        [प्रशिक्षण देखो :material-arrow-right:](nfcore_build/index.md){ .md-button .md-button--secondary }

    ??? catalog "**Side Quests:** Nextflow के उन्नत विषयों में गहराई से जाओ"

        स्टैंडअलोन मिनी-कोर्स का एक संग्रह जो उन Nextflow डेवलपर्स के लिए है जो अपनी क्षमताओं का विस्तार करना और/या विशेष विषयों पर अपने कौशल को गहरा करना चाहते हैं।
        इन्हें क्रमिक रूप से प्रस्तुत किया गया है लेकिन किसी भी क्रम में लिया जा सकता है (प्रत्येक मिनी-कोर्स के अवलोकन में निर्भरताएं देखो)।

        [Side Quests देखो :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} विज्ञान के लिए Nextflow {.mt-1}

    विशिष्ट वैज्ञानिक अनुप्रयोगों के लिए Nextflow पाइपलाइन विकसित करना सीखो।

    ??? courses "**Genomics:** वेरिएंट कॉलिंग पाइपलाइन विकसित करो"

        उन शोधकर्ताओं के लिए एक कोर्स जो अपनी खुद की जीनोमिक्स पाइपलाइन विकसित करना सीखना चाहते हैं, आवश्यक Nextflow विकास पैटर्न प्रदर्शित करने के लिए वेरिएंट कॉलिंग उपयोग के मामले का उपयोग करते हुए।

        [प्रशिक्षण देखो :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**RNAseq:** बल्क RNAseq प्रोसेसिंग पाइपलाइन विकसित करो"

        उन शोधकर्ताओं के लिए एक कोर्स जो अपनी खुद की RNAseq पाइपलाइन विकसित करना सीखना चाहते हैं, आवश्यक Nextflow विकास पैटर्न प्रदर्शित करने के लिए बल्क RNAseq प्रोसेसिंग उपयोग के मामले का उपयोग करते हुए।

        [प्रशिक्षण देखो :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Bioimaging:** इमेजिंग पाइपलाइन चलाओ और कॉन्फ़िगर करो"

        उन शोधकर्ताओं के लिए एक कोर्स जो बायोइमेजिंग पाइपलाइन चलाना और कॉन्फ़िगर करना सीखना चाहते हैं, आवश्यक Nextflow उपयोग पैटर्न प्रदर्शित करने के लिए nf-core/molkart का उपयोग करते हुए।

        [प्रशिक्षण देखो :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

</div>

## सेटअप और सहायता

<div class="grid cards mb-4" markdown>

-   :material-cog-outline:{ .lg .middle } __प्रशिक्षण वातावरण__

    ---

    Nextflow प्रशिक्षण के लिए अपना वातावरण सेट करने के विकल्प।

    [प्रशिक्षण वातावरण देखो :material-arrow-right:](envsetup/index.md){ .md-button .md-button--secondary }

-   :material-tag-outline:{ .lg .middle } __Nextflow संस्करण__

    ---

    Nextflow के सिंटैक्स संस्करणों के विकास को समझना और प्रबंधित करना।

    [संस्करण आवश्यकताएं जांचो :material-arrow-right:](info/nxf_versions.md){ .md-button .md-button--secondary }

-   :material-file-code-outline:{ .lg .middle } __Hello पाइपलाइन__

    ---

    Hello पाइपलाइन क्या करती है और यह कैसे संरचित है, इसका सारांश।

    [सारांश पढ़ो :material-arrow-right:](info/hello_pipeline.md){ .md-button .md-button--secondary }

-   :material-lifebuoy:{ .lg .middle } __सहायता प्राप्त करना__

    ---

    Nextflow प्रशिक्षण में समस्या होने पर उपयोगी संसाधन।

    [सहायता ढूंढो :material-arrow-right:](help.md){ .md-button .md-button--secondary }

</div>
