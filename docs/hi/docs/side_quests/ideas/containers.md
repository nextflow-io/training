# भाग 1: अधिक कंटेनर

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. कंटेनर इमेज कैसे खोजें या बनाएं

कुछ सॉफ़्टवेयर डेवलपर्स अपने सॉफ़्टवेयर के लिए कंटेनर इमेज प्रदान करते हैं जो Docker Hub जैसे कंटेनर रजिस्ट्री पर उपलब्ध होते हैं, लेकिन कई नहीं करते।
इस वैकल्पिक अनुभाग में, हम आपको उन टूल्स के लिए कंटेनर इमेज प्राप्त करने के दो तरीके दिखाएंगे जिन्हें आप अपने Nextflow pipelines में उपयोग करना चाहते हैं: Seqera Containers का उपयोग करके और कंटेनर इमेज को स्वयं बनाकर।

आप `quote` pip पैकेज के लिए कंटेनर इमेज प्राप्त/निर्माण करेंगे, जिसका उपयोग इस अनुभाग के अंत में अभ्यास में किया जाएगा।

### 1.1. Seqera Containers से कंटेनर इमेज प्राप्त करें

Seqera Containers एक मुफ्त सेवा है जो pip और conda (bioconda सहित) इंस्टॉल करने योग्य टूल्स के लिए कंटेनर इमेज बनाती है।
[Seqera Containers](https://www.seqera.io/containers/) पर जाएं और `quote` pip पैकेज को खोजें।

![Seqera Containers](img/seqera-containers-1.png)

`quote` pip पैकेज के लिए कंटेनर इमेज का अनुरोध करने के लिए "+Add" और फिर "Get Container" पर क्लिक करें।

![Seqera Containers](img/seqera-containers-2.png)

यदि यह पहली बार इस पैकेज के इस संस्करण के लिए community कंटेनर बनाया जा रहा है, तो इसे पूरा होने में कुछ मिनट लग सकते हैं।
आपके लिए बनाई गई कंटेनर इमेज के URI (जैसे `community.wave.seqera.io/library/pip_quote:ae07804021465ee9`) को कॉपी करने के लिए क्लिक करें।

अब आप Grace Hopper से एक यादृच्छिक कथन प्राप्त करने के लिए `quote` कमांड चलाने हेतु कंटेनर इमेज का उपयोग कर सकते हैं।

```bash
docker run --rm community.wave.seqera.io/library/pip_quote:ae07804021465ee9 quote "Grace Hopper"
```

आउटपुट:

```console title="Output"
Humans are allergic to change. They love to say, 'We've always done it
this way.' I try to fight that. That's why I have a clock on my wall
that runs counter-clockwise.
```

### 1.2. कंटेनर इमेज को स्वयं बनाएं

चलिए `quote` pip पैकेज के लिए कंटेनर इमेज को स्वयं बनाने के लिए Seqera Containers वेबसाइट से कुछ बिल्ड विवरण का उपयोग करें।
Seqera Containers वेबसाइट पर वापस जाएं और "Build Details" बटन पर क्लिक करें।

पहला आइटम जिसे हम देखेंगे वह `Dockerfile` है, एक प्रकार की स्क्रिप्ट फ़ाइल जिसमें कंटेनर इमेज बनाने के लिए आवश्यक सभी कमांड होते हैं।
हमने नीचे दिए गए Dockerfile में कुछ व्याख्यात्मक टिप्पणियां जोड़ी हैं ताकि आप समझ सकें कि प्रत्येक भाग क्या करता है।

```Dockerfile title="Dockerfile"
# micromamba बेस docker इमेज से शुरू करें
FROM mambaorg/micromamba:1.5.10-noble
# conda.yml फ़ाइल को कंटेनर में कॉपी करें
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Nextflow के उपयोग के लिए विभिन्न यूटिलिटीज और conda.yml फ़ाइल में पैकेज इंस्टॉल करें
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
# कंटेनर को root उपयोगकर्ता के रूप में चलाएं
USER root
# PATH वातावरण वेरिएबल को micromamba इंस्टॉलेशन डायरेक्टरी को शामिल करने के लिए सेट करें
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
```

दूसरा आइटम जिसे हम देखेंगे वह `conda.yml` फ़ाइल है, जिसमें कंटेनर इमेज में इंस्टॉल किए जाने वाले पैकेजों की सूची होती है।

```conda.yml title="conda.yml"
channels:
- conda-forge
- bioconda
dependencies:
- pip
- pip:
  - quote==3.0.0 #
```

इन फ़ाइलों की सामग्री को `containers/build` डायरेक्टरी में स्थित स्टब्स में कॉपी करें, फिर कंटेनर इमेज को स्वयं बनाने के लिए निम्नलिखित कमांड चलाएं।

!!! note "नोट"

    हम कंटेनर इमेज को `quote` नाम और `latest` टैग के साथ टैग करने के लिए `-t quote:latest` फ्लैग का उपयोग करते हैं।
    हम इस सिस्टम पर इसे चलाते समय कंटेनर इमेज को संदर्भित करने के लिए इस टैग का उपयोग कर सकेंगे।

```bash
docker build -t quote:latest containers/build
```

बनने के बाद, आप अभी बनाई गई कंटेनर इमेज को चला सकते हैं।

```bash
docker run --rm quote:latest quote "Margaret Oakley Dayhoff"
```

### निष्कर्ष

आपने अपने Nextflow pipelines में उपयोग करने के लिए किसी टूल के लिए कंटेनर इमेज प्राप्त करने के दो अलग-अलग तरीके सीखे हैं: Seqera Containers का उपयोग करके और कंटेनर इमेज को स्वयं बनाकर।

### आगे क्या है?

इस प्रशिक्षण श्रृंखला के [अगले अध्याय](./04_hello_genomics.md) में जारी रखने के लिए आपके पास सब कुछ है।
आप `quote` कंटेनर का उपयोग करके कंप्यूटर/जीवविज्ञान के अग्रदूतों पर उद्धरण प्राप्त करने और `cowsay` कंटेनर का उपयोग करके उन्हें आउटपुट करने के लिए एक वैकल्पिक अभ्यास के साथ भी जारी रख सकते हैं।

---

## 2. गाय को प्रसिद्ध वैज्ञानिकों का उद्धरण बोलवाएं

इस अनुभाग में कुछ stretch अभ्यास हैं, जो आपने अब तक सीखा है उसका अभ्यास करने के लिए।
इन अभ्यासों को करना प्रशिक्षण के बाद के भागों को समझने के लिए _आवश्यक नहीं_ है, लेकिन यह पता लगाकर कि गाय को प्रसिद्ध वैज्ञानिकों का उद्धरण कैसे बोलवाया जाए, आपकी सीख को मजबूत करने का एक मजेदार तरीका प्रदान करते हैं।

```console title="cowsay-output-Grace-Hopper.txt"
  _________________________________________________
 /                                                 \
| Humans are allergic to change. They love to       |
| say, 'We've always done it this way.' I try to fi |
| ght that. That's why I have a clock on my wall th |
| at runs counter-clockwise.                        |
| -Grace Hopper                                     |
 \                                                 /
  =================================================
                                                 \
                                                  \
                                                    ^__^
                                                    (oo)\_______
                                                    (__)\       )\/\
                                                        ||----w |
                                                        ||     ||
```

### 2.1. getQuote process का उपयोग करने के लिए `hello-containers.nf` स्क्रिप्ट को संशोधित करें

हमारे पास `containers/data/pioneers.csv` फ़ाइल में कंप्यूटर और जीवविज्ञान के अग्रदूतों की एक सूची है।
उच्च स्तर पर, इस अभ्यास को पूरा करने के लिए आपको यह करना होगा:

- डिफ़ॉल्ट `params.input_file` को `pioneers.csv` फ़ाइल की ओर इंगित करने के लिए संशोधित करें।
- एक `getQuote` process बनाएं जो प्रत्येक इनपुट के लिए उद्धरण प्राप्त करने के लिए `quote` कंटेनर का उपयोग करता है।
- उद्धरण प्रदर्शित करने के लिए `getQuote` process के आउटपुट को `cowsay` process से कनेक्ट करें।

`quote` कंटेनर इमेज के लिए, आप या तो पिछले stretch अभ्यास में स्वयं बनाई गई इमेज का उपयोग कर सकते हैं या Seqera Containers से प्राप्त की गई इमेज का उपयोग कर सकते हैं।

!!! tip "सुझाव"

    आपके getQuote process के `script` ब्लॉक के लिए एक अच्छा विकल्प हो सकता है:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

आप इस अभ्यास का समाधान `containers/solutions/hello-containers-4.1.nf` में पा सकते हैं।

### 2.2. अपने Nextflow pipeline को `quote` और `sayHello` मोड में निष्पादित होने की अनुमति देने के लिए संशोधित करें।

अपने pipeline में कुछ branching तर्क जोड़ें ताकि यह `quote` और `sayHello` दोनों के लिए इच्छित इनपुट स्वीकार कर सके।
यहां Nextflow workflow में `if` स्टेटमेंट का उपयोग कैसे करें इसका एक उदाहरण है:

```groovy title="hello-containers.nf"
workflow {
    if (params.quote) {
        ...
    }
    else {
        ...
    }
    cowSay(text_ch)
}
```

!!! tip "सुझाव"

    आप किसी process के आउटपुट channel को एक नाम असाइन करने के लिए `new_ch = processName.out` का उपयोग कर सकते हैं।

आप इस अभ्यास का समाधान `containers/solutions/hello-containers-4.2.nf` में पा सकते हैं।

### निष्कर्ष

आप जानते हैं कि processes चलाने के लिए Nextflow में कंटेनर का उपयोग कैसे करें, और अपने pipelines में कुछ branching तर्क कैसे बनाएं!

### आगे क्या है?

जश्न मनाएं, एक stretch ब्रेक लें और कुछ पानी पिएं!

जब आप तैयार हों, तो इस प्रशिक्षण श्रृंखला के भाग 3 पर जाएं ताकि आप सीख सकें कि अब तक आपने जो सीखा है उसे अधिक यथार्थवादी डेटा विश्लेषण उपयोग मामले पर कैसे लागू करें।
