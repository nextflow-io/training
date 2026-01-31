# मैन्युअल इंस्टॉलेशन

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

प्रशिक्षण चलाने के लिए आवश्यक सब कुछ अपने स्थानीय वातावरण में मैन्युअल रूप से install करना संभव है।

यहाँ हमने standard POSIX-compatible systems (व्यक्तिगत machine जैसे laptop मानते हुए) पर यह कैसे करना है इसका documentation किया है।
ध्यान रखो कि कुछ विवरण तुम्हारे विशिष्ट system के आधार पर अलग हो सकते हैं।

!!! tip "सुझाव"

    आगे बढ़ने से पहले, क्या तुमने [Devcontainers approach](03_devcontainer.md) पर विचार किया है?
    यह मैन्युअल installation की आवश्यकता के बिना सभी आवश्यक tools और dependencies प्रदान करता है।

## सामान्य software आवश्यकताएँ

Nextflow किसी भी POSIX-compatible system (Linux, macOS, Windows Subsystem for Linux, आदि) पर Java installed के साथ उपयोग किया जा सकता है।
हमारे प्रशिक्षण courses में कुछ अतिरिक्त आवश्यकताएँ हैं।

कुल मिलाकर, तुम्हें निम्नलिखित software installed होना चाहिए:

- Bash या समकक्ष shell
- [Java 11 (या बाद का, 21 तक)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (या बाद का)
- [VSCode](https://code.visualstudio.com) [Nextflow extension](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow) के साथ

VSCode application तकनीकी रूप से वैकल्पिक है लेकिन हम दृढ़ता से अनुशंसा करते हैं कि तुम इसे courses में काम करने के साथ-साथ सामान्य रूप से अपने Nextflow development कार्य के लिए उपयोग करो।

Nextflow documentation manual [Environment setup](https://www.nextflow.io/docs/latest/developer-env.html) के तहत इन dependencies को install करने के निर्देश प्रदान करता है।

## Nextflow और nf-core tools

तुम्हें Nextflow खुद install करना होगा, साथ ही nf-core tools, जैसा कि नीचे linked articles में विस्तृत है:

- [Nextflow installation](https://www.nextflow.io/docs/latest/install.html)
- [nf-core tools](https://nf-co.re/docs/nf-core-tools/installation)

हम Nextflow के लिए self-install विकल्प और nf-core tools के लिए PyPI विकल्प का उपयोग करने की अनुशंसा करते हैं।

!!! warning "Version compatibility"

    <!-- इस content में कोई भी update home page पर copy करना होगा -->
    **जनवरी 2026 के अनुसार, हमारे सभी Nextflow प्रशिक्षण courses के लिए Nextflow version 25.10.2 या बाद का आवश्यक है, strict v2 syntax सक्रिय के साथ, जब तक अन्यथा उल्लेख न हो।**

    Version आवश्यकताओं और strict v2 syntax के बारे में अधिक जानकारी के लिए, कृपया [Nextflow versions](../info/nxf_versions.md) guide देखो।

    पूर्व syntax के अनुरूप प्रशिक्षण सामग्री के पुराने versions इस webpage के menu bar में version selector के माध्यम से उपलब्ध हैं।

## प्रशिक्षण सामग्री

प्रशिक्षण सामग्री download करने का सबसे आसान तरीका इस command का उपयोग करके पूरी repository को clone करना है:

```bash
git clone https://github.com/nextflow-io/training.git
```

प्रत्येक course की अपनी डायरेक्टरी है।
किसी course में काम करने के लिए, एक terminal window खोलो (आदर्श रूप से, VSCode application के अंदर से) और संबंधित डायरेक्टरी में `cd` करो।

फिर तुम website पर दिए गए course instructions का पालन कर सकते हो।
