---
title: Nextflow संस्करण
description: Nextflow के syntax संस्करणों के विकास को समझना और प्रबंधित करना
hide:
  - toc
  - footer
---

## वर्तमान समर्थित Nextflow syntax संस्करण और आवश्यकताएँ

प्रशिक्षण पोर्टल के संस्करण 3.0 के अनुसार, हमारे सभी प्रशिक्षण पाठ्यक्रम Nextflow के 25.10.2 संस्करण रिलीज़ पर आधारित हैं, जब तक कि पाठ्यक्रम इंडेक्स पेज पर अन्यथा निर्दिष्ट न हो (पुराने या संग्रहीत सामग्री को छोड़कर जिसमें संस्करण सूचना शामिल नहीं हो सकती)।

क्योंकि पाठ्यक्रम अब workflow स्तर पर typed inputs के साथ-साथ workflow-level output directives का उपयोग करते हैं, इसलिए V2 syntax parser का उपयोग आवश्यक है।
यदि तुम [Github Codespaces](../envsetup/01_setup.md) या [local devcontainers](../envsetup/03_devcontainer.md) के माध्यम से हमारे द्वारा प्रदान किए गए वातावरण का उपयोग करने की योजना बना रहे हो, तो तुम्हें कुछ करने की आवश्यकता नहीं है जब तक कि पाठ्यक्रम निर्देशों में विशेष रूप से नोट न किया गया हो।
हालाँकि, यदि तुम अपने स्वयं के वातावरण में प्रशिक्षण के माध्यम से काम करने की योजना बना रहे हो ([Manual install](../envsetup/02_local.md)), तो तुम्हें v2 syntax parser सक्षम करके Nextflow संस्करण 25.10.2 या बाद के संस्करण का उपयोग करना सुनिश्चित करना होगा।

## प्रशिक्षण सामग्री के पुराने संस्करण

हमारी प्रशिक्षण सामग्री फ़रवरी 2025 से संस्करणित की गई है।

तुम प्रत्येक पृष्ठ के शीर्ष पर ड्रॉपडाउन मेनू आइटम के माध्यम से **25.10.2 से पहले** के Nextflow संस्करणों के साथ काम करने वाली प्रशिक्षण सामग्री के पुराने संस्करण तक पहुँच सकते हो जो प्रशिक्षण सामग्री का क्रमांकित संस्करण दिखाता है।
जब तुम प्रशिक्षण सामग्री का पुराना संस्करण चुनते हो, तो प्रशिक्षण वातावरण के लिंक स्वचालित रूप से वातावरण के संबंधित संस्करण को निर्दिष्ट करेंगे।

## Nextflow syntax संस्करणों के बारे में अन्य जानकारी

Nextflow में दो अलग-अलग versioning अवधारणाएँ हैं जो कभी-कभी भ्रमित होती हैं: **DSL संस्करण** और **syntax parser संस्करण**।

**DSL1 बनाम DSL2** Nextflow pipelines लिखने के मूल रूप से भिन्न तरीकों को संदर्भित करता है।
DSL1 मूल syntax था जहाँ processes channels के माध्यम से अंतर्निहित रूप से जुड़े थे।
DSL2, जो Nextflow 20.07 में पेश किया गया था, ने modularity सुविधाएँ जोड़ीं: अन्य फ़ाइलों से processes और workflows आयात करने की क्षमता, स्पष्ट `workflow` ब्लॉक, और named process outputs।
DSL1 को Nextflow 22.03 में deprecated किया गया और 22.12 में हटा दिया गया।
सभी आधुनिक Nextflow कोड DSL2 का उपयोग करते हैं।

**Syntax parser v1 बनाम v2** अलग-अलग parsers को संदर्भित करता है जो दोनों DSL2 कोड के साथ काम करते हैं।
v1 parser मूल, अधिक अनुमोदक parser है।
v2 parser सख्त है और नई भाषा सुविधाओं को सक्षम करता है जैसे static typing (typed inputs और outputs) और workflow-level output directives।
v2 parser बेहतर त्रुटि संदेश भी प्रदान करता है और runtime के बजाय parse समय पर अधिक त्रुटियों को पकड़ता है।
v2 parser Nextflow 26.04 में डिफ़ॉल्ट बन जाएगा।

संक्षेप में: DSL2 वह भाषा है जो तुम लिखते हो; syntax parser संस्करण यह निर्धारित करता है कि उस भाषा की कितनी सख्ती से व्याख्या की जाती है और कौन सी उन्नत सुविधाएँ उपलब्ध हैं।

### Nextflow संस्करण की जाँच और सेटिंग

तुम `nextflow --version` कमांड का उपयोग करके जाँच सकते हो कि तुम्हारे सिस्टम पर Nextflow का कौन सा संस्करण स्थापित है।

Nextflow के संस्करण को अपडेट करने के बारे में अधिक जानकारी के लिए, कृपया [Updating Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html) पर संदर्भ दस्तावेज़ देखें।

### v2 syntax parser सक्षम करना

अपने वर्तमान सत्र के लिए v2 syntax parser को **सक्षम** करने के लिए, अपने टर्मिनल में निम्न कमांड चलाओ:

```bash
export NXF_SYNTAX_PARSER=v2
```

इसे स्थायी बनाने के लिए (Nextflow 26.04 में v2 के डिफ़ॉल्ट बनने तक), export कमांड को अपने shell profile (`~/.bashrc`, `~/.zshrc`, आदि) में जोड़ो:

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

ध्यान दो कि `NXF_SYNTAX_PARSER=v2` environment variable एक अस्थायी आवश्यकता है।
Nextflow 26.04 से आगे, v2 parser डिफ़ॉल्ट बन जाएगा और इस सेटिंग की आवश्यकता नहीं होगी।

### v2 syntax parser अक्षम करना

अपने वर्तमान सत्र के लिए v2 syntax parser को **अक्षम** करने के लिए, अपने टर्मिनल में निम्न कमांड चलाओ:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### मौजूदा कोड का migration

Nextflow के हाल के संस्करणों के अनुरूप मौजूदा कोड के migration के लिए मार्गदर्शन हेतु, कृपया संदर्भ दस्तावेज़ में [Migration Notes](https://www.nextflow.io/docs/latest/migrations/index.html) देखें।

सबसे हाल के रिलीज़ में migrate करने के लिए ये दो लेख विशेष रूप से सहायक हैं:

- [Migrating to workflow outputs](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migrating to static types](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

इन दोनों सुविधाओं को प्रशिक्षण सामग्री के संस्करण 3.0 से शुरू होने वाले beginner training के हिस्से के रूप में कवर किया गया है।

तुम जिस Nextflow कोड को migrate करना चाहते हो उसकी पीढ़ी के आधार पर, तुम `nextflow lint -format` कमांड का उपयोग करके Nextflow linter द्वारा इसका अधिकांश भाग पूरा करने में सक्षम हो सकते हो।
अधिक विवरण के लिए [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) के लिए CLI reference देखें।

हम आशा करते हैं कि यह सहायक होगा।
यदि तुम्हें सहायता की आवश्यकता है, तो Slack या forum पर संपर्क करो।

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
