# GitHub Codespaces

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

GitHub Codespaces एक web-आधारित platform है जो हमें प्रशिक्षण के लिए पूर्व-configured वातावरण प्रदान करने की अनुमति देता है, जो cloud में virtual machines द्वारा समर्थित है।
यह platform Github (जो Microsoft के स्वामित्व में है) द्वारा संचालित है, और Github account वाले किसी भी व्यक्ति के लिए मुफ्त में (उपयोग quotas के साथ) उपलब्ध है।

!!! warning "चेतावनी"

    संगठनों से जुड़े accounts कुछ अतिरिक्त प्रतिबंधों के अधीन हो सकते हैं।
    अगर तुम्हारी स्थिति ऐसी है, तो तुम्हें एक स्वतंत्र व्यक्तिगत account का उपयोग करना पड़ सकता है, या इसके बजाय स्थानीय installation का उपयोग करना पड़ सकता है।

## GitHub account बनाना

तुम [GitHub home page](https://github.com/) से एक मुफ्त GitHub account बना सकते हो।

## अपना GitHub Codespace शुरू करना

एक बार जब तुम GitHub में logged in हो, Nextflow प्रशिक्षण वातावरण खोलने के लिए इस link को अपने browser में खोलो: <https://codespaces.new/nextflow-io/training?quickstart=1&ref=master>

वैकल्पिक रूप से, तुम नीचे दिखाए गए button पर click कर सकते हो, जो प्रत्येक प्रशिक्षण course में दोहराया गया है (आमतौर पर Orientation page पर)।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

तुम्हें एक page दिखाई देना चाहिए जहाँ तुम एक नया GitHub Codespace बना सकते हो:

![Create a GitHub Codespace](img/codespaces_create.png)

### Configuration

सामान्य उपयोग के लिए, तुम्हें कुछ भी configure करने की ज़रूरत नहीं होनी चाहिए।
जब तक course में अन्यथा निर्दिष्ट न हो जो तुम शुरू कर रहे हो, तुम बस मुख्य button पर click करके आगे बढ़ सकते हो।

हालाँकि, "Change options" button पर click करके वातावरण को customize करना संभव है।

??? info "Configuration विकल्प"

    अगर तुम "Change options" button पर click करते हो, तो तुम्हें निम्नलिखित को customize करने का विकल्प दिया जाएगा:

    #### Branch

    यह तुम्हें प्रशिक्षण सामग्री का एक अलग version चुनने की अनुमति देता है।
    `master` branch में आम तौर पर bug fixes और सामग्री होती है जो हाल ही में विकसित और approved हुई है लेकिन अभी तक website पर release नहीं हुई है।
    अन्य branches में work in progress होता है जो पूर्णतः functional नहीं हो सकता।

    #### Machine type

    यह तुम्हें प्रशिक्षण में काम करने के लिए उपयोग की जाने वाली virtual machine को customize करने की अनुमति देता है।

    अधिक cores वाली machine का उपयोग करने से तुम Nextflow की workflow execution को parallelize करने की क्षमता का अधिक लाभ उठा सकते हो।
    हालाँकि, यह तुम्हारे मुफ्त quota allocation को तेज़ी से consume करेगा, इसलिए हम इस setting को बदलने की अनुशंसा नहीं करते जब तक कि course के instructions में ऐसा करने की सलाह न दी गई हो जो तुम लेने की योजना बना रहे हो।

    Quotas के बारे में अधिक जानकारी के लिए नीचे 'GitHub Codespaces quotas' देखो।

### Startup समय

पहली बार एक नया GitHub Codespaces वातावरण खोलने में कई मिनट लग सकते हैं, क्योंकि system को तुम्हारी virtual machine सेट करनी होती है, इसलिए अगर कुछ देर लगे तो चिंता मत करो।
हालाँकि, इसमें पाँच मिनट से अधिक नहीं लगना चाहिए।

## प्रशिक्षण interface में navigate करना

एक बार तुम्हारा GitHub Codespaces load हो जाए, तुम्हें कुछ इस तरह दिखना चाहिए (जो तुम्हारी account preferences के आधार पर light mode में खुल सकता है):

![GitHub Codespaces welcome](img/codespaces_welcome.png)

यह VSCode IDE का interface है, एक लोकप्रिय code development application जिसे हम Nextflow development के लिए उपयोग करने की अनुशंसा करते हैं।

- **Main editor** वह जगह है जहाँ Nextflow code और अन्य text files खुलेंगी। यहाँ तुम code edit करोगे। जब तुम codespace खोलते हो, यह तुम्हें `README.md` file का preview दिखाएगा।
- Main editor के नीचे **terminal** तुम्हें commands चलाने की अनुमति देता है। यहाँ तुम course instructions में दी गई सभी command lines चलाओगे।
- **Sidebar** तुम्हें अपना वातावरण customize करने और बुनियादी कार्य करने की अनुमति देता है (copy, paste, files खोलना, search, git, आदि)। डिफ़ॉल्ट रूप से यह file explorer पर खुला होता है, जो तुम्हें repository की सामग्री browse करने की अनुमति देता है। Explorer में किसी file पर click करने से वह main editor window में खुलेगी।

तुम window panes के सापेक्ष अनुपात को अपनी पसंद के अनुसार adjust कर सकते हो।

<!-- TODO (future) Link to development best practices side quest? -->

## GitHub Codespaces के उपयोग के बारे में अन्य नोट्स

### Session फिर से शुरू करना

एक बार जब तुमने वातावरण बना लिया, तुम आसानी से इसे फिर से शुरू या restart कर सकते हो और जहाँ छोड़ा था वहाँ से जारी रख सकते हो।
तुम्हारा वातावरण 30 मिनट की निष्क्रियता के बाद timeout हो जाएगा और तुम्हारे बदलावों को 2 सप्ताह तक save करेगा।

तुम <https://github.com/codespaces/> से वातावरण फिर से खोल सकते हो।
पिछले वातावरण सूचीबद्ध होंगे।
इसे फिर से शुरू करने के लिए किसी session पर click करो।

![List GitHub Codespace sessions](img/codespaces_list.png)

अगर तुमने अपने पिछले GitHub Codespaces वातावरण का URL save किया है, तो तुम बस इसे अपने browser में खोल सकते हो।
वैकल्पिक रूप से, उसी button पर click करो जिसका उपयोग तुमने पहली बार इसे बनाने के लिए किया था:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

तुम्हें पिछला session दिखना चाहिए, डिफ़ॉल्ट विकल्प इसे फिर से शुरू करना है:

![Resume a GitHub Codespace](img/codespaces_resume.png)

### Files को अपनी स्थानीय machine पर save करना

Explorer panel से कोई भी file save करने के लिए, file पर right-click करो और `Download` चुनो।

### GitHub Codespaces quotas का प्रबंधन

GitHub Codespaces तुम्हें प्रति माह 15 GB-month storage और 120 core-hours प्रति माह देता है।
यह standard workspace (2 cores, 8 GB RAM, और 32 GB storage) का उपयोग करके लगभग 60 घंटे के डिफ़ॉल्ट वातावरण runtime के बराबर है।

तुम उन्हें अधिक resources के साथ बना सकते हो (ऊपर explanation देखो), लेकिन यह तुम्हारा मुफ्त उपयोग तेज़ी से consume करेगा और तुम्हारे पास इस space तक पहुँच के कम घंटे होंगे।
उदाहरण के लिए, अगर तुम 2-core डिफ़ॉल्ट के बजाय 4-core machine चुनते हो, तो तुम्हारा quota आधे समय में समाप्त हो जाएगा।

वैकल्पिक रूप से, तुम अधिक resources तक पहुँच खरीद सकते हो।

अधिक जानकारी के लिए, GitHub documentation देखो:
[About billing for GitHub Codespaces](https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)
