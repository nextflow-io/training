# भाग 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. गाय को प्रसिद्ध वैज्ञानिकों के उद्धरण बुलवाएं

इस खंड में कुछ अतिरिक्त अभ्यास हैं, जो आपने अब तक जो सीखा है उसका अभ्यास करने के लिए।
इन अभ्यासों को करना प्रशिक्षण के बाद के भागों को समझने के लिए _आवश्यक नहीं_ है, लेकिन यह गाय को प्रसिद्ध वैज्ञानिकों के उद्धरण बुलवाने का तरीका सीखकर अपनी शिक्षा को मजबूत करने का एक मजेदार तरीका प्रदान करता है।

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

### 1.1. `hello-containers.nf` स्क्रिप्ट को getQuote प्रोसेस का उपयोग करने के लिए संशोधित करें

हमारे पास `containers/data/pioneers.csv` फ़ाइल में कंप्यूटर और जीव विज्ञान के अग्रदूतों की एक सूची है।
उच्च स्तर पर, इस अभ्यास को पूरा करने के लिए आपको निम्नलिखित करना होगा:

- डिफ़ॉल्ट `params.input_file` को `pioneers.csv` फ़ाइल की ओर इंगित करने के लिए संशोधित करें।
- एक `getQuote` process बनाएं जो प्रत्येक इनपुट के लिए उद्धरण लाने के लिए `quote` कंटेनर का उपयोग करे।
- `getQuote` process के आउटपुट को `cowsay` process से कनेक्ट करें ताकि उद्धरण प्रदर्शित हो सके।

`quote` कंटेनर इमेज के लिए, आप या तो पिछले अतिरिक्त अभ्यास में स्वयं बनाई गई इमेज का उपयोग कर सकते हैं या Seqera Containers से प्राप्त की गई इमेज का उपयोग कर सकते हैं।

!!! Hint

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

### 1.2. अपनी Nextflow पाइपलाइन को `quote` और `sayHello` मोड में निष्पादित करने की अनुमति देने के लिए संशोधित करें।

अपनी पाइपलाइन में शाखा तर्क जोड़ें ताकि यह `quote` और `sayHello` दोनों के लिए इच्छित इनपुट स्वीकार कर सके।
यहाँ Nextflow workflow में `if` स्टेटमेंट का उपयोग करने का एक उदाहरण है:

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

!!! Hint

    आप किसी process के आउटपुट channel को नाम देने के लिए `new_ch = processName.out` का उपयोग कर सकते हैं।

आप इस अभ्यास का समाधान `containers/solutions/hello-containers-4.2.nf` में पा सकते हैं।

### निष्कर्ष

आप जानते हैं कि प्रोसेस चलाने के लिए Nextflow में कंटेनर का उपयोग कैसे करें, और अपनी पाइपलाइन में कुछ शाखा तर्क कैसे बनाएं!

### आगे क्या?

जश्न मनाएं, थोड़ा ब्रेक लें और पानी पिएं!

जब आप तैयार हों, तो इस प्रशिक्षण श्रृंखला के भाग 3 पर जाएं और सीखें कि आपने अब तक जो सीखा है उसे अधिक यथार्थवादी डेटा विश्लेषण उपयोग के मामले में कैसे लागू करें।
