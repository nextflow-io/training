# पाठ्यक्रम सारांश

Hello nf-core प्रशिक्षण पाठ्यक्रम पूरा करने पर बधाई! 🎉

<!-- placeholder for video -->

## आपकी यात्रा

आपने एक डेमो pipeline को retrieve और run करना सीखकर शुरुआत की, फिर एक सरल Nextflow workflow को nf-core pipeline में बदलने का कार्य किया।
आपने सीखा कि कैसे template का उपयोग करके pipeline scaffold बनाया जाए और मौजूदा pipeline को उस scaffold पर लगाया जाए।
फिर आपने धीरे-धीरे pipeline को परिष्कृत किया - एक local मॉड्यूल को nf-core मॉड्यूल से बदला, एक और local मॉड्यूल को nf-core मानकों के अनुरूप बनाया, और इनपुट validation जोड़ा।

### आपने क्या बनाया

आपकी अंतिम `core-hello` pipeline में अब है:

- **मानकीकृत संरचना** nf-core template का उपयोग करते हुए workflows, subworkflows, modules, और configuration के लिए व्यवस्थित डायरेक्टरी के साथ
- **कम्युनिटी मॉड्यूल** nf-core repository से (`cat/cat`) आपके custom modules के साथ
- **व्यापक validation** जो pipeline के चलने से पहले पैरामीटर और इनपुट डेटा दोनों की जांच करता है
- **व्यावसायिक configuration** विभिन्न execution वातावरण के लिए profiles के साथ
- **पूर्ण documentation** और metadata nf-core conventions का पालन करते हुए

### प्रमुख कौशल अर्जित किए

इस hands-on पाठ्यक्रम के माध्यम से, आपने सीखा है:

1. मौजूदा pipeline का पता लगाकर nf-core pipeline संरचना को **समझना और नेविगेट करना**
2. Workflows को nf-core template के भीतर **पुनर्संरचित करना** और composable बनाना
3. कम्युनिटी repository से पहले से बने modules को **खोजना और integrate करना**
4. Naming, संरचना और metadata के लिए nf-core मानकों का पालन करते हुए **custom modules बनाना**
5. स्पष्ट feedback के साथ त्रुटियों को जल्दी पकड़ने के लिए nf-schema का उपयोग करके **validation लागू करना**

अब आप कम्युनिटी best practices का पालन करने वाली production-ready nf-core pipelines बनाने के लिए आवश्यक ज्ञान से सुसज्जित हैं।

## अपने कौशल बढ़ाने के अगले कदम

यहाँ हमारे शीर्ष 3 सुझाव हैं कि आगे क्या करें:

- [Nextflow for Science](../nf4_science/index.md) के साथ एक वैज्ञानिक विश्लेषण use case में Nextflow को लागू करें
- [Side Quests](../side_quests/index.md) के साथ अधिक उन्नत Nextflow features देखें
- [nf-core कम्युनिटी में शामिल होकर](https://nf-co.re/join) सक्रिय भागीदार बनें।

अंत में, हम अनुशंसा करते हैं कि आप [**Seqera Platform**](https://seqera.io/) पर एक नज़र डालें, एक cloud-based platform जो Nextflow के निर्माताओं द्वारा विकसित किया गया है जो आपके workflows को launch और manage करना और भी आसान बनाता है, साथ ही आपके डेटा को manage करना और किसी भी वातावरण में interactively विश्लेषण चलाना संभव बनाता है।

## प्रतिक्रिया सर्वेक्षण

आगे बढ़ने से पहले, कृपया पाठ्यक्रम सर्वेक्षण पूरा करने में एक मिनट का समय दें! आपकी प्रतिक्रिया हमें सभी के लिए अपनी प्रशिक्षण सामग्री को बेहतर बनाने में मदद करती है।

[सर्वेक्षण लें :material-arrow-right:](survey.md){ .md-button .md-button--primary }
