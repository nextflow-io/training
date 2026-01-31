# भाग 4: Hello Modules - ट्रांसक्रिप्ट

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "महत्वपूर्ण नोट"

    यह पेज केवल ट्रांसक्रिप्ट दिखाता है। पूर्ण चरण-दर-चरण निर्देशों के लिए, [पाठ्यक्रम सामग्री](../04_hello_modules.md) पर वापस जाएं।

    ट्रांसक्रिप्ट में दिखाए गए अनुभाग क्रमांक केवल संकेत के उद्देश्य से प्रदान किए गए हैं और सामग्री में सभी अनुभाग क्रमांक शामिल नहीं हो सकते हैं।

## स्वागत

नमस्ते, Hello Nextflow प्रशिक्षण पाठ्यक्रम के भाग चार में आपका स्वागत है।

इस अध्याय को Hello Modules कहा जाता है, और हम Nextflow कोड को modularize करने के बारे में बात करेंगे। हम जो करने जा रहे हैं वह है अपनी एक workflow स्क्रिप्ट को लेकर उसे अलग-अलग फ़ाइलों में विभाजित करना।

यह कोड को navigate और maintain करना आसान बनाता है जैसे-जैसे आपका workflow बड़ा होता जाता है, और pipelines के बीच modules को साझा करना भी संभव बनाता है ताकि यदि आपके पास एक ही टूल का उपयोग करने वाले कई pipelines हैं, तो आपको उस process को केवल एक बार लिखने की आवश्यकता हो।

इसका एक क्लासिक उदाहरण nf-core modules रिपॉजिटरी है, जिसमें हजारों विभिन्न टूल्स ready to go modules में हैं, जिन्हें आप अपने workflow में install और उपयोग कर सकते हैं।

Nextflow sub workflows के साथ भी काम कर सकता है, जो modules की तरह हैं, लेकिन उनमें कई processes होते हैं। यह इस प्रशिक्षण के दायरे से बाहर है, लेकिन यह मूल रूप से उसी तरह काम करता है।

ठीक है। चलिए देखते हैं।

हमेशा की तरह, training.nextflow.io पर जाकर शुरू करें।

साइडबार में "Hello Nextflow" पर जाएं, और हम भाग चार कर रहे हैं: "Hello Modules"।

मैं अब अपने GitHub Code Spaces वातावरण में जाऊंगा और "hello-modules" फ़ाइल को देखूंगा।

पहले की तरह, हम पिछले अध्याय के अंतिम बिंदु से शुरू कर रहे हैं, इसलिए यह स्क्रिप्ट परिचित होनी चाहिए। हमारे पास तीन processes हैं, say hello, convert to upper और collect greetings, और एक सरल workflow में, जो इन तीन commands को चलाता है और अंत में एक message emit करता है। हमारे पास greeting और batch नामक दो parameters हैं, जो नाम निर्दिष्ट करता है, जिसका उपयोग अंत में collected आउटपुट फ़ाइल के लिए किया जाता है।

## 0. Warmup: hello-modules.nf रन करें

हम यह सत्यापित कर सकते हैं कि यह workflow अभी भी हमारी अपेक्षा के अनुसार काम करता है nextflow run hello, modules करके।

बढ़िया। इसने इस process के साथ तीन कार्य रन किए, एक collected कार्य, और यह हमें बताया कि इस batch में तीन greetings हैं। अगर हम results में जाते हैं, तो हमारे पास यहां हमारी विभिन्न आउटपुट फ़ाइलें हैं, जिसमें collected test, batch आउटपुट भी शामिल है।

## 1. modules संग्रहीत करने के लिए एक डायरेक्टरी बनाएं

ठीक है। चलिए कुछ modularization करते हैं।

आम तौर पर अपनी pipeline रिपॉजिटरी में modules को एक subfolder में रखना एक अच्छा विचार है, सिर्फ चीजों को व्यवस्थित रखने के लिए। आप इसे जो चाहें कह सकते हैं, लेकिन परंपरागत रूप से हम इसे modules कहते हैं।

तो चलिए आगे बढ़ते हैं, एक टर्मिनल में जाएं और make the modules करें। आप इसे VS Code में साइडबार में पॉप अप होते देख सकते हैं।

## 2. sayHello() के लिए एक module बनाएं

मैं फिर अपने पहले module के लिए एक नई फ़ाइल बनाने जा रहा हूं। आप "touch" या "code" कर सकते हैं या आप इसे साइडबार में कर सकते हैं, वास्तव में कोई फर्क नहीं पड़ता। तो मैं code modules करने जा रहा हूं और मैं इसे process के नाम पर रखूंगा। तो sayHello.nf। NF, Nextflow फ़ाइलों के लिए एक पारंपरिक फ़ाइल एक्सटेंशन है।

यहां save दबाने जा रहा हूं और हम अपनी नई module फ़ाइल को दिखाई देते देख सकते हैं।

## 2.2. sayHello process कोड को module फ़ाइल में ले जाएं

ठीक है, अगला मैं workflow से module कोड लेने जा रहा हूं। मैं यहां hash bang भी लेने जा रहा हूं और पहले उसे कॉपी करूंगा ताकि यह स्पष्ट रूप से एक Nextflow फ़ाइल हो। और फिर मैं इस process को लेने जा रहा हूं और मैं cut करूंगा। तो मैं इसे अपनी मुख्य workflow स्क्रिप्ट से हटा दूंगा और मैं इसे इस नए module में paste करूंगा।

बस इतनी ही सामग्री इस module फ़ाइल में होगी। सिर्फ एक single process, कोई workflow नहीं, कोई logic नहीं, सिर्फ अकेली process।

मैं अब इस फ़ाइल को बंद कर सकता हूं।

## 2.3. workflow ब्लॉक से पहले एक import घोषणा जोड़ें

अब मेरे workflow में वह पहली process गायब है, इसलिए हमें इसे import करके वापस लाने की आवश्यकता है। इसके लिए syntax अन्य प्रोग्रामिंग भाषाओं के समान है, इसलिए यह परिचित लग सकता है। हम include घुंघराले ब्रैकेट्स करते हैं, process का नाम, say hello, और फिर फ़ाइल पथ modules से, say hello, nf। शानदार।

यहां कुछ tricks हैं। VS Code एक्सटेंशन इस बारे में चतुर है। यह इस फ़ाइल पथ को पहचानता है और आप इस पर hover कर सकते हैं और follow link कर सकते हैं। या मैं Mac पर हूं, मैं option click कर सकता हूं और यह इस फ़ाइल को खोलता है। इसलिए हम जल्दी से इस पर जा सकते हैं।

यह process नाम अब यहां नीचे workflow द्वारा उपयोग किया जा रहा है, और हम यहां भी वही काम कर सकते हैं। यह हमें उस process के बारे में थोड़ी जानकारी दिखाता है, और फिर से, मैं option hold कर सकता हूं, इस पर click कर सकता हूं, और यह इसे editor में खोल देगा।

तो यह VS Code में आपके विभिन्न processes के लिए जब आपके पास बहुत सारी फ़ाइलें हों तो अपने code base में जल्दी से navigate करने का वास्तव में तेज तरीका है।

ठीक है। इस अध्याय के लिए मूल रूप से बस इतना ही। अब हम बाकी processes के लिए फिर से वही काम करते हैं।

## 3. convertToUpper() process को modularize करें

तो चलिए यहां एक नई फ़ाइल बनाते हैं। इसे Convert to upper nf कहें। फिर से, hash bang कॉपी करें। और फिर process को cut करें।

वहां process का नाम कॉपी करें, नए process नाम के साथ एक नया include statement शामिल करें।

## 4. collectGreetings() process को modularize करें

और फिर तीसरी process के लिए भी वही करें। नई फ़ाइल, connect। Greetings,

hash bang करें। Process को cut करें, process को paste करें, और एक नया include statement करें।

अब आप यहां देख सकते हैं कि मुझे यहां एक error underline मिली है जो invalid include source कह रही है। और यह वास्तव में एक वास्तविक error है जो मैंने की क्योंकि मैं बहुत तेजी से आगे बढ़ रहा था। अगर आप ध्यान से देखें, तो आप देख सकते हैं कि मैंने T miss किया और convert to upper में

तो VS Code ने बहुत उपयोगी रूप से मुझे बताया है कि मैंने वहां गलती की है। अगर मैं उस फ़ाइल का नाम ठीक करता हूं, तो error दूर हो जाती है। यह एक अच्छा उदाहरण है कि Nextflow कोड लिखने के लिए VS Code के भीतर error checking इतना उपयोगी क्यों है। अन्यथा मैंने इसे नहीं पकड़ा होता और मुझे बहुत बाद में पता चलता जब मैं workflow को चलाने की कोशिश करता।

हमारी मुख्य pipeline स्क्रिप्ट अब बहुत सरल दिख रही है। इसमें कोई processes नहीं हैं, हमारे पास सिर्फ तीन include statements और हमारा workflow है। हमने workflow के logic में कोई बदलाव नहीं किया है। हमने process कोड में कोई बदलाव नहीं किया है, इसलिए उम्मीद है कि यह बिल्कुल उसी तरह काम करना चाहिए।

## 4.4. सत्यापित करने के लिए workflow चलाएं कि यह पहले की तरह ही काम करता है

चलिए जांच करते हैं। मैं एक टर्मिनल खोलने जा रहा हूं और मैं पहले की तरह बिल्कुल वही कमांड चलाने जा रहा हूं।

निश्चित रूप से, इसने हमारी processes चलाई है, say hello, convert to upper collect greetings, और हमें फिर से तीन greetings दी हैं।

तो हमने अपने कोड को इधर-उधर किया है, लेकिन हमने workflow के execute होने के तरीके के बारे में कुछ भी नहीं बदला है और यह पूरी तरह से अपरिवर्तित है। केवल अंतर यह है कि अब हमारे पास cleaner कोड है, maintain करना आसान है, और दूसरों के साथ साझा करना आसान है।

और बस इतना ही। यह एक छोटा अध्याय था। यह एक सरल अवधारणा है, लेकिन यह बहुत शक्तिशाली है और अधिक जटिल Nextflow workflows लिखने के तरीके की कुंजी है। इसलिए यह महत्वपूर्ण है कि आप इसे समझें और इसका उपयोग करने की आदत डालें।

अगले अध्याय में, हम गति में थोड़ा बदलाव करने जा रहे हैं और Nextflow कोड लिखने के syntax के बारे में इतना सोचना बंद करेंगे, और processes में खुद software का उपयोग कैसे करें इसके बारे में थोड़ा सोचेंगे। भाग पांच के लिए हमसे जुड़ें Hello Containers में।

[अगला वीडियो ट्रांसक्रिप्ट :octicons-arrow-right-24:](05_hello_containers.md)
