# शुरुआत करना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube channel पर [पूरी playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) देखें।

:green_book: वीडियो transcript [यहाँ](./transcripts/00_orientation.md) उपलब्ध है।
///

!!! tip "सुझाव"

    YouTube वीडियो में कुछ super powers हैं!

    - :fontawesome-solid-closed-captioning: उच्च गुणवत्ता (manually curated) captions / subtitles। इन्हें :material-subtitles: icon से चालू करो
    - :material-bookmark: Timeline में video chapters जो page headings से मेल खाते हैं।

-->

## Training environment शुरू करें

GitHub Codespaces पर हमारे द्वारा प्रदान किए गए pre-built environment का उपयोग करने के लिए, नीचे "Open in GitHub Codespaces" button पर क्लिक करो। अन्य विकल्पों के लिए, [Environment options](../envsetup/index.md) देखें।

हम training environment को एक नए browser tab या window में खोलने की सलाह देते हैं (अपने equipment के आधार पर right-click, ctrl-click या cmd-click का उपयोग करो) ताकि तुम environment लोड होने के दौरान पढ़ना जारी रख सको।
कोर्स के माध्यम से काम करने के लिए तुम्हें इन instructions को parallel में खुला रखना होगा।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Environment की मूल बातें

इस training environment में training course के माध्यम से काम करने के लिए आवश्यक सभी software, code और data शामिल हैं, इसलिए तुम्हें खुद कुछ भी install करने की आवश्यकता नहीं है।

Codespace एक VSCode interface के साथ सेट किया गया है, जिसमें एक filesystem explorer, एक code editor और एक terminal shell शामिल है।
कोर्स के दौरान दिए गए सभी instructions (जैसे 'file खोलो', 'code edit करो' या 'यह command चलाओ') VSCode interface के उन तीन भागों को संदर्भित करते हैं जब तक अन्यथा निर्दिष्ट न हो।

यदि तुम इस कोर्स को स्वयं कर रहे हो, तो कृपया अधिक विवरण के लिए [environment basics](../envsetup/01_setup.md) से परिचित हो जाओ।

### Version requirements

यह training Nextflow 25.10.2 या बाद के version के लिए डिज़ाइन की गई है **v2 syntax parser ENABLED के साथ**।
यदि तुम local या custom environment का उपयोग कर रहे हो, तो कृपया सुनिश्चित करो कि तुम [यहाँ](../info/nxf_versions.md) documented सही settings का उपयोग कर रहे हो।

## काम करने के लिए तैयार हो जाओ

एक बार जब तुम्हारा codespace चल रहा हो, तो training में dive करने से पहले तुम्हें दो चीजें करनी होंगी: इस specific course के लिए अपनी working directory सेट करो, और प्रदान की गई materials पर एक नज़र डालो।

### Working directory सेट करें

डिफ़ॉल्ट रूप से, codespace सभी training courses के root पर work directory सेट के साथ खुलता है, लेकिन इस course के लिए, हम `hello-nextflow/` directory में काम करेंगे।

Terminal में यह command चलाकर अभी directory बदलो:

```bash
cd hello-nextflow/
```

तुम VSCode को इस directory पर focus करने के लिए सेट कर सकते हो, ताकि file explorer sidebar में केवल relevant files दिखें:

```bash
code .
```

!!! tip "सुझाव"

    यदि किसी भी कारण से तुम इस directory से बाहर चले जाते हो (जैसे तुम्हारा codespace sleep हो जाए), तो तुम हमेशा full path का उपयोग करके इसमें वापस आ सकते हो, यह मानते हुए कि तुम इसे Github Codespaces training environment में चला रहे हो:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

अब चलो contents पर एक नज़र डालते हैं।

### प्रदान की गई materials को explore करें

तुम training workspace के बाईं ओर file explorer का उपयोग करके इस directory की contents को explore कर सकते हो।
वैकल्पिक रूप से, तुम `tree` command का उपयोग कर सकते हो।

पूरे course में, हम directory structure और contents को readable form में represent करने के लिए `tree` के output का उपयोग करते हैं, कभी-कभी clarity के लिए minor modifications के साथ।

यहाँ हम दूसरे level तक table of contents generate करते हैं:

```bash
tree . -L 2
```

??? abstract "Directory contents"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── test-params.yaml
    ```

Section को expand करने और इसकी contents देखने के लिए colored box पर क्लिक करो।
हम expected command output को concise तरीके से शामिल करने के लिए इस तरह के collapsible sections का उपयोग करते हैं।

- **`.nf` files** workflow scripts हैं जिनके नाम इस आधार पर हैं कि वे course के किस भाग में उपयोग किए जाते हैं।

- **`nextflow.config` file** एक configuration file है जो minimal environment properties सेट करती है।
  तुम इसे अभी के लिए ignore कर सकते हो।

- **`data/` के तहत `greetings.csv` file** में input data है जिसे हम course के अधिकांश भाग में उपयोग करेंगे। इसे Part 2 (Channels) में describe किया गया है, जब हम इसे पहली बार introduce करते हैं।

- **`test-params.*` files** configuration files हैं जिनका उपयोग हम Part 6 (Configuration) में करेंगे। तुम इन्हें अभी के लिए ignore कर सकते हो।

- **`solutions` directory** में completed workflow scripts हैं जो course के प्रत्येक step से result होती हैं।
  ये तुम्हारे काम की जाँच करने और किसी भी issue को troubleshoot करने के लिए reference के रूप में उपयोग की जानी हैं।

## तैयारी checklist

सोचते हो कि तुम dive करने के लिए तैयार हो?

- [ ] मैं इस course का goal और इसकी prerequisites समझता/समझती हूँ
- [ ] मेरा environment up और running है
- [ ] मैंने अपनी working directory appropriately सेट की है

यदि तुम सभी boxes check कर सकते हो, तो तुम जाने के लिए तैयार हो।

**[Part 1: Hello World](./01_hello_world.md) पर जारी रखने के लिए, इस page के नीचे दाएँ कोने में arrow पर क्लिक करो।**
