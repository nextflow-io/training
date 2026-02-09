# कोर्स सारांश

Hello Nextflow प्रशिक्षण कोर्स पूरा करने पर बधाई! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [Nextflow YouTube चैनल पर पूरी प्लेलिस्ट](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) देखें।

:green_book: तुम वीडियो के साथ [वीडियो ट्रांसक्रिप्ट](./transcripts/07_next_steps.md) पढ़ सकते हो।
///

## तुम्हारी यात्रा

तुमने एक बहुत ही बेसिक वर्कफ़्लो से शुरुआत की थी जो एक हार्डकोडेड कमांड चलाता था।
छह भागों के दौरान, तुमने उस बेसिक वर्कफ़्लो को एक मॉड्यूलर मल्टी-स्टेप पाइपलाइन में बदल दिया जो Nextflow की प्रमुख विशेषताओं का उपयोग करती है जिसमें चैनल, ऑपरेटर, कंटेनर के लिए बिल्ट-इन सपोर्ट, और कॉन्फ़िगरेशन विकल्प शामिल हैं।

### तुमने क्या बनाया

- Hello वर्कफ़्लो का अंतिम रूप इनपुट के रूप में एक CSV फ़ाइल लेता है जिसमें टेक्स्ट अभिवादन होते हैं।
- चार स्टेप Nextflow प्रोसेस (`sayHello`, `convertToUpper`, `collectGreetings`, और `cowpy`) के रूप में लागू किए गए हैं जो अलग-अलग मॉड्यूल फ़ाइलों में स्टोर हैं।
- परिणाम `results/` नामक डायरेक्टरी में प्रकाशित किए जाते हैं।
- पाइपलाइन का अंतिम आउटपुट एक प्लेन टेक्स्ट फ़ाइल है जिसमें एक कैरेक्टर की ASCII आर्ट है जो अपरकेस अभिवादन बोल रहा है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** प्रत्येक अभिवादन को अपनी आउटपुट फ़ाइल में लिखता है (_जैसे_ "Hello-output.txt")
2. **`convertToUpper`:** प्रत्येक अभिवादन को अपरकेस में बदलता है (_जैसे_ "HELLO")
3. **`collectGreetings`:** सभी अपरकेस अभिवादनों को एक सिंगल बैच फ़ाइल में इकट्ठा करता है
4. **`cowpy`:** `cowpy` टूल का उपयोग करके ASCII आर्ट जेनरेट करता है

वर्कफ़्लो कॉन्फ़िगरेशन लचीले, पुनरुत्पादनीय तरीके से इनपुट और पैरामीटर प्रदान करने का समर्थन करता है।

### प्राप्त कौशल

इस हैंड्स-ऑन कोर्स के माध्यम से, तुमने सीखा है कि कैसे:

- एक सरल मल्टी-स्टेप वर्कफ़्लो बनाने के लिए पर्याप्त कोर Nextflow कंपोनेंट का वर्णन और उपयोग करें
- ऑपरेटर और चैनल फ़ैक्टरी जैसी अगले-स्तर की अवधारणाओं का वर्णन करें
- लोकल रूप से Nextflow वर्कफ़्लो लॉन्च करें
- Nextflow द्वारा जेनरेट किए गए आउटपुट (परिणाम) और लॉग फ़ाइलों को खोजें और व्याख्या करें
- बेसिक समस्याओं का निवारण करें

अब तुम Nextflow में अपनी खुद की पाइपलाइन विकसित करना शुरू करने के लिए बुनियादी ज्ञान से लैस हो।

## अपने कौशल को बढ़ाने के लिए अगले कदम

यहाँ हमारे शीर्ष 3 सुझाव हैं कि आगे क्या करना है:

- [Nextflow for Science](../nf4_science/index.md) के साथ एक वैज्ञानिक विश्लेषण उपयोग केस में Nextflow लागू करो
- [Hello nf-core](../hello_nf-core/index.md) के साथ nf-core से शुरुआत करो
- [Side Quests](../side_quests/index.md) के साथ अधिक उन्नत Nextflow विशेषताओं का अन्वेषण करो

अंत में, हम अनुशंसा करते हैं कि तुम [**Seqera Platform**](https://seqera.io/) पर एक नज़र डालो, Nextflow के निर्माताओं द्वारा विकसित एक क्लाउड-आधारित प्लेटफ़ॉर्म जो तुम्हारे वर्कफ़्लो को लॉन्च और प्रबंधित करना और भी आसान बनाता है, साथ ही तुम्हारे डेटा को प्रबंधित करता है और किसी भी वातावरण में इंटरैक्टिव रूप से विश्लेषण चलाता है।

## फ़ीडबैक सर्वे

आगे बढ़ने से पहले, कृपया कोर्स सर्वे पूरा करने में एक मिनट का समय दो! तुम्हारा फ़ीडबैक हमें सभी के लिए अपनी प्रशिक्षण सामग्री को बेहतर बनाने में मदद करता है।

[सर्वे लें :material-arrow-right:](survey.md){ .md-button .md-button--primary }
