# भाग 3: Hello Workflow - वीडियो ट्रांसक्रिप्ट

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "महत्वपूर्ण नोट"

    यह पेज केवल ट्रांसक्रिप्ट दिखाता है। पूर्ण चरण-दर-चरण निर्देशों के लिए, [कोर्स सामग्री](../03_hello_workflow.md) पर वापस जाओ।

    ट्रांसक्रिप्ट में दिखाए गए सेक्शन नंबर केवल संकेतात्मक उद्देश्यों के लिए प्रदान किए गए हैं और सामग्री में सभी सेक्शन नंबर शामिल नहीं हो सकते हैं।

## स्वागत और रीकैप

नमस्ते, और Hello Nextflow के भाग तीन में वापस स्वागत है। इस भाग को Hello Workflow कहा जाता है, और कोर्स के इस भाग में हम वास्तव में pipeline या workflow नाम को सही ठहराना शुरू करते हैं।

हम अपनी सरल pipeline script को, जिसमें अभी तक एक process है, लेंगे और अतिरिक्त processes जोड़ना शुरू करेंगे और देखेंगे कि Nextflow इस orchestration और pipeline के माध्यम से data flow को कैसे हैंडल करता है।

चलो वापस अपने code spaces पर जाते हैं। तुम देखोगे कि मैंने अपनी सभी .nextflow\* डायरेक्टरीज और work डायरेक्टरीज और सब कुछ साफ रखने के लिए डिलीट कर दी है। चिंता मत करो अगर तुम्हारे पास अभी भी कोर्स के पिछले भागों से वे फ़ाइलें पड़ी हैं।

हम hello-workflow.nf नामक फ़ाइल से काम करने वाले हैं। पहले की तरह, यह मूल रूप से उस script को दर्शाता है जिसे हमने अब तक बनाया है, और हमें एक साफ शुरुआती बिंदु देता है। और फिर से, output में नीचे हम देख सकते हैं कि path अब hello_workflow है। तो published फ़ाइलें तुम्हारे results फ़ोल्डर में एक अलग subdirectory में जानी चाहिए।

रीकैप करने के लिए कि हम अब तक कहाँ हैं, हमारे पास यहाँ एक single process है, एक input greeting के साथ, एक output greeting फ़ाइल। और फिर simple Bash script, जो बस एक फ़ाइल में echo कमांड करता है।

हमारे पास एक single workflow input है, यहाँ params block, जहाँ हम कहते हैं कि यह एक path की उम्मीद कर रहा है, और डिफ़ॉल्ट data/greetings.csv है, जो यह फ़ाइल यहाँ ऊपर है।

फिर workflow में ही, हमारे पास एक main block है। हम एक चैनल बना रहे हैं। हम CSV को rows में parse कर रहे हैं और फिर प्रत्येक array के पहले element को ले रहे हैं, और हम उस चैनल को उस process में पास कर रहे हैं, जो फिर तीन tasks generate कर रहा है, और हम workflow से, उस process के outputs publish कर रहे हैं।

और फिर अंत में, output block में, हम Nextflow को बता रहे हैं कि इस चैनल से इन फ़ाइलों को hello_workflow नामक डायरेक्टरी में publish करे। और उन फ़ाइलों को soft link करने के बजाय copy करे।

## 1. workflow में दूसरा स्टेप जोड़ो

ठीक है, इस भाग में हम अपने workflow में एक दूसरा process जोड़ने वाले हैं। हम sayHello process के outputs लेंगे, और उन्हें एक दूसरे स्टेप में process करेंगे, जो उन फ़ाइलों के भीतर सभी अक्षरों को convertToUppercase करेगा।

यह बस एक silly उदाहरण है, यह फिर से कुछ सरल string processing है, लेकिन यह तुम्हें दिखाता है कि हम workflow के भीतर logic को कैसे ले सकते हैं।

हम इसके लिए "tr" नामक एक bash कमांड करने वाले हैं, जो translate के लिए short है। यह एक Unix कमांड है जो हमेशा से है। अगर तुम इससे परिचित नहीं हो, तो मैं तुम्हें दोष नहीं देता। मुझे नहीं लगता कि मैंने training से पहले कभी इसका उपयोग किया था, लेकिन तुम इसे terminal पर बहुत जल्दी try कर सकते हो। अगर मैं "echo 'hello world'" करूँ और फिर 'tr' में pipe करूँ और फिर quotes में तुम character range कहो, तो A से Z, lowercase, और फिर तुम A से Z uppercase करना चाहते हो। और यह बस कहता है, इन अक्षरों को इन अक्षरों में translate करो।

और जब मैं enter दबाता हूँ, तुम देख सकते हो कि यह अब सब कुछ capitalize कर चुका है। बहुत अच्छा अगर तुम्हें लोगों पर चिल्लाना पसंद है।

तो यह bash कमांड की एक बहुत सरल शैली है जिसे हम अपने दूसरे process में उपयोग करने वाले हैं।

## 1.2. uppercasing स्टेप को Nextflow process के रूप में लिखो

तो अगर मैं अपनी script पर वापस जाता हूँ, मैं थोड़ा चीट करने जा रहा हूँ और training के docs से कोड कॉपी करूँगा। लेकिन तुम बिल्कुल देख सकते हो कि क्या हो रहा है।

हमारे पास यहाँ एक नया process है। इसे हमने convertToUpper कहा है, लेकिन हम इसे जो चाहें कह सकते हैं।

हमारे पास एक single input path है, जैसा कि हमने पहले किया था। यह एक value चैनल नहीं है, यह एक path चैनल है। और फिर एक single output।

script block में हम input फ़ाइल पर "cat" करते हैं। और हम इसे squiggly brackets में रख सकते हैं अगर हम चाहें। और जो उस variable को लेता है। और हम उसी bash कमांड को pipe में चलाते हैं और हम परिणाम को इस फ़ाइल नाम के साथ एक फ़ाइल में लिखते हैं, और वह output path द्वारा pick up होता है।

अब हमें इस नए process के साथ कुछ करने की आवश्यकता है। तो हम workflow में जाने वाले हैं जहाँ हम एक workflow के विभिन्न logic को build करते हैं, और उस पहले process के बाद, हम अपना दूसरा process चलाने जा रहे हैं। तो convertToUpper यहाँ process का नाम है।

यह एक input लेता है इसलिए हम इसे अकेले नहीं बुला सकते। हम पहले process के output को process करना चाहते हैं। तो जैसे हमने इसके साथ किया, sayHello out जहाँ हम उन results को publish कर रहे हैं। हम उन्हीं results को यहाँ input के रूप में उपयोग करना चाहते हैं, इसलिए हम उन्हें copy कर सकते हैं और उन्हें वहाँ डाल सकते हैं।

हम sayHello process ".out" चाहते हैं, और Nextflow जानता है कि इसका मतलब यहाँ एक simple single output record है, जो यह फ़ाइल है। तो वह फिर एक दूसरे process को input के रूप में pass किया जाएगा।

## 1.5. workflow output publishing सेट करो

ठीक है। और अंत में, ताकि हम वास्तव में इस दूसरे process के results को save करें, हमें उन्हें workflow से भी publish करने की आवश्यकता है, और फिर उन्हें output block में define करना है, पहले जैसा ही syntax। तो हम इसे copy कर सकते हैं और second outputs कह सकते हैं, या जो भी तुम इसे कहना चाहते हो।

उस process का नाम लो जिसमें हम रुचि रखते हैं, convertToUpper out, और फिर यहाँ नीचे output block में। इसे add करो और हम यहाँ same attributes कर सकते थे। तो हम भी चाहते हैं कि ये फ़ाइलें Hello Workflow subdirectory में हों, और हम भी उन्हें copy करना चाहते हैं।

बढ़िया। चलो इसे चलाने की कोशिश करते हैं। तो अगर मैं terminal लाता हूँ और मैं "nextflow run hello-workflow.nf" करता हूँ, और हम देखेंगे कि यह क्या करता है। देखते हैं कि क्या यह पिछले भागों से अलग दिखता है।

तो यह Nextflow launch करता है। docs में, यह "-resume" के साथ करने के लिए कहता है, लेकिन मैंने अपनी सभी work डायरेक्टरी delete कर दी थी, इसलिए यहाँ कोई फर्क नहीं पड़ता। लेकिन अगर तुमने किया होता, तो वह भी काम करेगा।

और यह लगभग बिल्कुल same दिखता है। लेकिन तुम अब देख सकते हो कि यहाँ output की एक दूसरी line है, जहाँ तुम उस दूसरे process का नाम देख सकते हो जिसे हमने अभी जोड़ा है। और वास्तव में, तुम देख सकते हो कि यह तीन बार successfully चला।

शानदार। अगर मेरी पिछली work डायरेक्टरीज आसपास होतीं और मैंने यह "-resume" के साथ किया होता, तो ये pipeline में बस पहले स्टेप को cached किया गया होता। क्योंकि वे outputs बिल्कुल same थे, इसलिए Nextflow जानता होता कि उन्हें फिर से reuse करना है।

और इसलिए तुम देख सकते हो कि तुम -resume का उपयोग करके अपने workflow को iteratively कैसे build कर सकते हो, step by step, अगर तुम्हें ज़रूरत हो।

ठीक है, चलो यहाँ ऊपर results डायरेक्टरी में देखते हैं और देखते हैं कि क्या यह काम किया। हम देख सकते हैं कि हमारे पास यहाँ ऊपर कुछ और फ़ाइलें हैं। हमें हमारी पहले जैसी original फ़ाइलें मिली हैं पहले process से। और वास्तव में, हमारे पास हमारी upper फ़ाइलें हैं और अक्षर सभी uppercase हैं, इसलिए यह काम कर गया। यह देखना वास्तव में अच्छा है।

इन work डायरेक्टरीज के अंदर check करना भी दिलचस्प है। पहले की तरह, यहाँ hash work डायरेक्टरीज से मेल खाता है। तो अगर मैं "ls work" में देखता हूँ, और फिर उसे expand करता हूँ, तो हम यहाँ विभिन्न फ़ाइलें देखेंगे।

हम पहले process से output फ़ाइल देखते हैं, जिसे यहाँ input के रूप में pull किया गया है। और हम नई output फ़ाइल देख सकते हैं जो generate की गई थी।

अब अगर मैं "-la" के साथ list करने और सभी फ़ाइलें दिखाने के लिए यह करता हूँ, तो हम कुछ और चीजें देखेंगे। सबसे पहले, तुम देखोगे कि यह फ़ाइल वास्तव में पहले process के लिए एक soft link है। यह मूल रूप से हमेशा एक soft link है अगर यह हो सकता है, फ़ाइल स्पेस बचाने के लिए। हम यहाँ फ़ाइलों को publish नहीं कर रहे हैं और यह बस उस फ़ाइल को पहले task से दूसरे task में reference करता है ताकि सब कुछ उस एक working डायरेक्टरी के भीतर encapsulate हो, और बाकी सब चीजों से safe और isolated हो।

और यह वहाँ होना चाहिए क्योंकि अगर हम .command.sh फ़ाइल को देखें, तो अगर मैं "cat work/b8/56\*" करता हूँ, तो तुम देख सकते हो कि, यहाँ फ़ाइल parts relative हैं, तो यह उस input फ़ाइल को cat कर रहा है, जिसे उसी working डायरेक्टरी में soft linked किया गया है।

तो हर work डायरेक्टरी ऐसी ही दिखेगी। जब तुम Nextflow में देखते हो, तो तुम्हारे पास वहाँ सभी input फ़ाइलें उस work डायरेक्टरी में staged होंगी। और फिर तुम्हारे पास कोई भी output फ़ाइलें भी होंगी जो बनाई गई थीं। तो यह बढ़िया है। वह ऐसा दिखता है जैसा हम उम्मीद करते हैं।

## 2.1. collection कमांड define करो और terminal में test करो

ठीक है, चलो अपने workflow पर वापस जाते हैं। अगला स्टेप क्या है जो हम करना चाहते हैं?

हमारे पास अब दो processes हैं और वे इस एक CSV फ़ाइल को ले रहे हैं, उसे parse कर रहे हैं और split कर रहे हैं। और फिर हम इन प्रत्येक processes के लिए तीन tasks रखते हैं और Nextflow इन सभी की parallelization को handle करता है, इसलिए यह सब side by side चलता है जहाँ संभव हो।

काम को parallel में चलाने के लिए split करने का यह तरीका बहुत सामान्य है। और उसका inverse फिर सब कुछ वापस gather करना है। तो हम अपने workflow में अपने अंतिम process के साथ यही करने वाले हैं, हमारे पास यहाँ एक तीसरा होगा, जो इन तीन अलग-अलग outputs को लेता है और उन सभी को एक single फ़ाइल में combine करता है।

हम यह terminal में काफी simply कर सकते हैं, बस यह महसूस करने के लिए कि यह कैसा दिखेगा।

अगर मैं results फ़ोल्डर में जाता हूँ। तो, "cd results/hello_workflow/", और हमारे पास यहाँ सभी UPPER फ़ाइलें हैं। मैं बस "cat" का उपयोग कर सकता हूँ, जिसका उपयोग हम उस फ़ाइल की contents को print करने के लिए करते हैं, और तुम "cat" को multiple फ़ाइलें दे सकते हो और यह एक के बाद एक पढ़ेगा।

तो मैं "UPPER-\*" कह सकता हूँ, जो मुझे Bash expansion के साथ तीन फ़ाइल नामों की same list देता है। और मैं combined.txt कह सकता हूँ। मुझे लगता है कि docs में, यह exact फ़ाइल नाम list करता है, लेकिन यह same चीज़ कर रहा है।

अब, अगर मैं "cat combined.txt" का उपयोग करता हूँ, तो हम देख सकते हैं कि हमारे पास उन तीनों फ़ाइलों की फ़ाइल contents हैं।

तो मूल रूप से यह process यही करने वाला है कि हम उसे एक previous process से सभी अलग-अलग output फ़ाइलें एक single process task में देने की कोशिश करने जा रहे हैं, और फिर हम उन्हें "cat" करके एक साथ जोड़ने जा रहे हैं और output फ़ाइल को save करेंगे।

## 2.2. collection स्टेप करने के लिए एक नया process बनाओ

ठीक है, तो चलो अपना नया process add करते हैं। मैं इसे training सामग्री से paste करने वाला हूँ, और तुम देख सकते हो कि इसने हमें इन question marks के साथ reader के लिए थोड़ा exercise छोड़ा है। लेकिन तुम process की सामान्य outline देख सकते हो मूल रूप से वही है जो हमने अभी terminal में किया था, जहाँ हम input फ़ाइलों के एक bunch का "cat" कर रहे हैं और इसे collected नामक एक output फ़ाइल में लिख रहे हैं, और फिर output फिर से उस single path की उम्मीद करता है।

तो हमें यहाँ किसी प्रकार के input की आवश्यकता है और वे paths का एक set होने वाले हैं। तो फिर से, हम एक input path चैनल define करते हैं और चलो इसे input_files कहते हैं। अब, इसने हमें पहले यहाँ एक single path दिया है, लेकिन एक path में यहाँ multiple फ़ाइलें भी हो सकती हैं, भले ही यह अभी भी एक single declaration है।

मैं इसे यहाँ नीचे copy करने वाला हूँ क्योंकि हम इन फ़ाइलों को "cat" करना चाहते हैं। और तुम सोच सकते हो कि हमें यहाँ array print करने या इस तरह की चीजों के साथ कुछ issues हैं, लेकिन Nextflow आम तौर पर इस मामले में काफी sensible है। और अगर इसे इस तरह multiple फ़ाइलों के साथ एक चैनल दिया जाता है, तो यह, उन सभी को space separators के साथ एक साथ रखेगा। तो यह हमें सही syntax देगा।

यह बढ़िया है। तो अब चलो अपने नए process को wire up करते हैं। मैं workflow में जाता हूँ। मैं outputs को combine करने जा रहा हूँ, नया process नाम, और बस पहले की तरह। मैं इस previous process को लूँगा, convertToUpper और ".out" करूँगा।

बढ़िया। चलो इसे try करते हैं और देखते हैं कि क्या यह terminal में काम करता है। अगर मैं बस कुछ डायरेक्टरीज ऊपर जाता हूँ और फिर Nextflow कमांड rerun करता हूँ, और हम देखेंगे कि क्या होता है।

तो workflow launch हो गया है और अब तुम देख सकते हो कि हमारे पास तीन अलग-अलग process नाम हैं, जो बढ़िया है। पहले दोनों पहले जैसे ही दिखते हैं, और तीसरा नया चलता है, जो अच्छा है।

हालाँकि, यहाँ कुछ अजीब है। हम उन output फ़ाइलों को एक single फ़ाइल में combine करना चाहते थे, और फिर भी इस process को हम देख सकते हैं कि तीन बार चला है, एक बार नहीं।

वास्तव में, अगर हम इन work डायरेक्टरीज में से एक में जाते हैं। और "cat work/" "collected" करते हैं, तो हम देखेंगे। यहाँ केवल एक single शब्द है, तीन नहीं।

और इसलिए जो हुआ है वह यह है कि Nextflow ने उस parallelization को जारी रखा है जैसे उसने पिछले steps में किया था। और इस process ने हमें तीन elements के साथ एक चैनल दिया, और वे तीन चैनल elements हमारे downstream process को pass किए गए, जिसने तीन process tasks generate किए।

यह मूल रूप से तीन अलग बार collect करने की कोशिश की और हर बार उसके पास बस एक single फ़ाइल थी, इसलिए यह बस एक output में single फ़ाइल को cat किया, और वास्तव में, हम .command.sh फ़ाइल में भी देख सकते हैं।

अगर मैं .command.sh करता हूँ, तो हम देख सकते हैं कि यहाँ बस एक single फ़ाइल नाम है और केवल एक single फ़ाइल को उस working डायरेक्टरी में staged किया गया था।

## 2.3. workflow में collection स्टेप add करो

तो किसी तरह हमें Nextflow को बताने की ज़रूरत है कि एक previous process से उन सभी outputs को एक साथ लाए और उन्हें इस downstream process को एक single चैनल element के रूप में दे, तीन के बजाय।

हम यह _collect_ नामक एक चैनल ऑपरेटर के साथ करते हैं।

यह एक बहुत useful ऑपरेटर है, जिसे तुम Nextflow pipelines में हर समय देखोगे। यह यहाँ एक चैनल है, यह output चैनल, बिल्कुल वही जो हमने ऊपर बनाया था। और इसलिए हम चैनल ऑपरेटरों को इसमें append कर सकते हैं जैसे हमने पहले किया था। हम बस dot कर सकते हैं, और फिर इस मामले में, collect, brackets।

और बस इतना ही। यह फिर इस process में pass होने से पहले इस चैनल को manipulate करेगा।

अगर तुम देखना चाहते हो कि इसका क्या हो रहा है, तो हम इसे यहाँ view भी कर सकते हैं। तो यहाँ, यह इस process को चलाने से बिल्कुल भी related नहीं है, इसलिए मैं इसे उस process को चलाने के बाद किसी भी बिंदु पर रख सकता हूँ। लेकिन हम same, output चैनल लेते हैं, और हम इसे .view के साथ देख रहे हैं, और फिर हम इसे .collect.view के साथ फिर से देख रहे हैं।

और जब हम यह चलाते हैं, तो यह हमें उस चैनल की दो अलग-अलग structures दिखाएगा, collect से पहले और बाद में। तो चलो अब इसे try करते हैं। ठीक है, मैंने थोड़ा zoom out किया है क्योंकि कुछ outputs काफी लंबे हैं, लेकिन अगर मैं pipeline चलाता हूँ, तो हम देखेंगे कि क्या यह काम करता है।

मैं आशा कर रहा हूँ कि एक तीसरा process सिर्फ एक बार चलेगा, क्योंकि यह outputs को collect कर रहा है और वास्तव में, तुम देख सकते हो कि collectGreetings तीन में से एक है। तो यह बस एक task चला।

और फिर अगर हम view statements को देखें, तो हमारे पास पहले के तीन elements के लिए तीन view statements हैं, प्रत्येक में एक फ़ाइल path के साथ।

और फिर उस collect statement के बाद, यह सिर्फ एक बार triggered हुआ क्योंकि उस चैनल में एक single element है। और अब हमारे पास तीन अलग-अलग फ़ाइल paths की यह, list है।

यह बिल्कुल वही है जो हम चाहते थे। और तुम देख सकते हो उम्मीद है, यह मूल रूप से उस "map" ऑपरेटर का inverse है जो हमने CSV arrays से अलग चैनल elements में जाने के लिए किया था। अब हम अलग चैनल elements ले रहे हैं और उन्हें एक single array में वापस रख रहे हैं।

बढ़िया, हम इन view statements को clear कर सकते हैं। हमें इनकी अब ज़रूरत नहीं है। हम अगले स्टेप पर जा सकते हैं।

इससे आगे बढ़ने से पहले, और भूलने से पहले, मैं यहाँ एक नया publish statement add करने वाला हूँ। Third output। तुम इसे अपने workflow में कुछ अधिक semantic और descriptive कह सकते हो। और फिर मैं इसे output block में फिर से add करूँगा और कहूँगा path 'hello_workflow' mode 'copy'। बस ताकि इस process द्वारा generate की गई output फ़ाइल यहाँ ऊपर हमारे results फ़ोल्डर में save हो जाए।

बस जल्दी से double check करने के लिए कि यह काम करता है। अब थोड़ा cleaner होना चाहिए क्योंकि हमारे पास वे view statements नहीं हैं। और, हम देखेंगे कि क्या हमें यहाँ ऊपर हमारी नई output फ़ाइल मिलती है। एक, एक task चला, collected नामक एक नई फ़ाइल मिली, और अब हमारे पास वे तीनों शब्द हैं। शानदार। आगे क्या है?

## 3. process को अतिरिक्त पैरामीटर pass करो

ठीक है। अब हम एक single process में multiple inputs को handle करना देखने वाले हैं। अब तक तुम देख सकते हो कि हमारे सभी processes बस एक चीज़ को input के रूप में ले रहे हैं। उन सभी के input के नीचे एक single line है।

हम यह demonstrate करने जा रहे हैं कि Nextflow को एक अलग batch identifier specify करने की अनुमति दें ताकि शायद तुम यह, workflow कई बार चला सको और तुम हर बार एक अलग batch ID दे सको।

मैं यहाँ collectGreetings के input में बस एक दूसरी line add करने वाला हूँ। और मैं इसे "val" कहने वाला हूँ, क्योंकि यह एक string है। अब यह एक value है, path नहीं, और मैं इसे "batch_name" कहने वाला हूँ।

फिर मैं इस variable का उपयोग करने के लिए यहाँ नीचे script को edit करने वाला हूँ, और मैं इसे training material में उसी स्थान पर डालने की कोशिश करूँगा। तो मैं इसे इस फ़ाइल path COLLECTED-$\{batch_name\}-output के बीच में डालता हूँ।

अभी पूरा नहीं हुआ। याद रखो कि हमें Nextflow को बताना होगा कि output फ़ाइल नाम क्या होने जा रहे हैं। तो हमें यहाँ ऊपर भी same चीज़ करनी होगी: COLLECTED-$\{batch_name\}-output.txt"।

शानदार। Nextflow अब एक दूसरा variable input प्राप्त कर रहा है और यह उसे script और output में interpolate कर रहा है।

एक आखिरी चीज़, अब हमें पता लगाना होगा कि इसे कहाँ call किया जा रहा है, और हमें process को दूसरा input pass करना होगा। यह किसी भी अन्य भाषा में किसी भी अन्य function में input की तरह है।

जैसे हमने training में पहले किया था, मैं यहाँ special "params" का उपयोग करने जा रहा हूँ, और हम इसे "params.batch" कहने वाले हैं ताकि हमारे पास -- batch CLI option हो सके। और अब तुम देख सकते हो कि यहाँ हमारे process में दो अलग inputs हैं बस comma separated, जो pass किए जा रहे हैं।

order को सही प्राप्त करना वास्तव में महत्वपूर्ण है, इसलिए यहाँ चैनल के लिए arguments का order और फिर param match करना चाहिए। यहाँ, चैनल और batch name। यह सिर्फ positional matching है।

ठीक है। मैं अब सीधे --batch के साथ इस pipeline को चला सकता हूँ, लेकिन चलो पहले सही चीज़ करते हैं और इसे यहाँ Params में input में define करते हैं। तो मैं इसे batch में add करने वाला हूँ और फिर हम कहेंगे कि यह एक string है और चलो इसे एक डिफ़ॉल्ट दें। तो चलो इसे बस batch कहते हैं। ठीक है? अब workflow चलाने की कोशिश करते हैं।

--batch Trio. मुझे लगता है कि यह training material में कहता है, लेकिन हम वहाँ कोई भी string उपयोग कर सकते हैं। और उम्मीद है कि हम उस results output फ़ाइल को यहाँ ऊपर आते हुए देखेंगे।

और वास्तव में, COLLECTED-trio-output - यह ठीक से काम कर गया है। इसने हमारी फ़ाइल को rename कर दिया है। और तुम अब imagine कर सकते हो कि यह उपयोगी है क्योंकि अगर मैं इसे फिर से एक अलग batch नाम के साथ चलाता हूँ, जैसे replicate_two, तो यह हमें यहाँ ऊपर एक अलग batch नाम देने वाला है।

और और यह फिर इस मामले में output फ़ाइलों को clobber नहीं करेगा। तो यह अच्छा है।

## 4. collector स्टेप में एक output add करो

ठीक है, तो अब हमारे पास यहाँ हमारे process में multiple inputs हैं। लेकिन क्या होता है अगर हम multiple outputs बनाना चाहते हैं? फिर हमारा उदाहरण यह है कि हम इस process के लिए एक report बनाने जा रहे हैं, बस यह कहते हुए कि कितनी फ़ाइलें collected की गईं।

और हम यह यहाँ एक echo कमांड के साथ करेंगे। तो हम कह सकते हैं echo। There were, मैं इसे training material से copy करने वाला हूँ, ताकि तुम्हें मुझे इसे type करते हुए नहीं देखना पड़े।

There were $\{count_greetings\} greetings in this batch, और इसे अब $\{batch_name\} नामक एक नई फ़ाइल में save करो, तो same variable, हम इसे जितनी बार चाहें reuse कर सकते हैं, report.txt।

## 4.1.1. collected greetings की संख्या count करो

हमें वास्तव में उसे किसी तरह calculate करने की आवश्यकता है। हम यह logic Bash script में कर सकते थे अगर हम चाहते, Bash logic का उपयोग करके। हालाँकि, हम Nextflow code के भीतर सीधे scripting भी कर सकते हैं, जब तक कि यह process में script block के भीतर है और quoted सेक्शन के ऊपर है।

यहाँ कुछ भी final rendered script में शामिल नहीं होगा, और यह बस Nextflow द्वारा execute किया जाएगा जब यह एक task render करता है।

तो यहाँ हम बस कुछ logic कर रहे हैं। हम count_greetings नामक एक नया variable बना रहे हैं। हम यहाँ input फ़ाइलें चैनल लेते हैं, और हम उस पर .size() call कर रहे हैं।

ठीक है, वह function मुझे यहाँ इस variable में एक नंबर देने वाला है, और अब हमारी warning चली गई है क्योंकि यह variable define किया जा रहा है।

ठीक है, तो हम work डायरेक्टरी में वह दूसरी फ़ाइल बना रहे हैं, लेकिन हमें Nextflow को बताने की ज़रूरत है कि इसे इस process के एक published output के रूप में उम्मीद करे। तो हम यह बिल्कुल उसी syntax से करते हैं जैसे हमने पहली फ़ाइल के लिए किया था।

हम कहते हैं path क्योंकि यह, फिर से, हम यहाँ एक variable publish कर सकते थे अगर हम चाहते "val" के साथ, लेकिन हम कहने वाले हैं "path"। और फिर expected, फ़ाइल नाम। ध्यान दो यह यहाँ highlighted नहीं है। यह इसलिए है क्योंकि मैंने single quotes का उपयोग किया। मुझे double quotes का उपयोग करना होगा।

## 4.1.2. report फ़ाइल emit करो और outputs को नाम दो

ठीक है, यह बढ़िया है। और अब हम इन outputs को यहाँ नीचे access करना शुरू कर सकते हैं जैसे मैंने यहाँ किया था। लेकिन अब यह अलग-अलग objects का एक array है, इसलिए मैं पहला प्राप्त करने के लिए collectGreetings.out[0] कर सकता हूँ, या दूसरा प्राप्त करने के लिए एक, जो हमारी नई report है।

लेकिन मुझे यह करना वास्तव में पसंद नहीं है क्योंकि index counting को mess up करना काफी आसान है। और तुम lines count करते हुए वहाँ बहुत बैठते हो और तुम एक नया output add करते हो और अचानक सब कुछ break हो जाता है। तो

इसके बजाय name से सब कुछ reference करना बहुत अच्छा है। और हम यहाँ "emit" नामक एक special key के साथ यह कर सकते हैं।

तो हम इसे जो चाहें कह सकते हैं। चलो इसे emit outfile कहते हैं, और emit reports। अगर तुम इन्हें define करते हो और तुम इसे एक या कई पर कर सकते हो, यह तुम्हारे ऊपर है। अब मैं यहाँ नीचे जा सकता हूँ और इसके बजाय मैं dot out dot reports जा सकता हूँ और बस इसे name से call कर सकता हूँ, जो तुम्हारे code को पढ़ते समय समझना बहुत आसान है, और code में changes के लिए safer है।

मैंने, यहाँ .out.report add किया है, लेकिन वास्तव में मुझे दो अलग-अलग outputs published होने की ज़रूरत है। तो मैं इसे कुछ अधिक दिलचस्प के रूप में rename करने वाला हूँ जैसे collected और report और क्या मैंने इसे वही कहा? मैंने इसे out फ़ाइल कहा, sorry। तो वह emit नाम यहाँ outfile और report। क्योंकि हम दो अलग-अलग output चैनलों को publish कर रहे हैं और इसलिए हमें publish block में दोनों को reference करने की ज़रूरत है।

फिर हमें output block में भी इन्हें define करने की ज़रूरत है। तो मैंने उसे collected rename किया, और फिर से, reports के लिए, यहाँ थोड़ा verbose, लेकिन यह वास्तव में उपयोगी है जब तुम एक नया workflow पढ़ने आते हो, तो यहाँ सभी अलग-अलग outputs को देखने के लिए, सभी अलग-अलग चैनलों को side by side list किया गया है, और इसे कम verbose बनाने के तरीके हैं, जिन्हें हम बाद में touch करेंगे।

ठीक है, चलो इसे try करते हैं और अपना workflow चलाते हैं और देखते हैं कि क्या होता है।

उम्मीद है कि अब यह मूल रूप से पहले की तरह ही चलना चाहिए। और हम यहाँ ऊपर replicate_two नामक एक नई output फ़ाइल प्राप्त करने जा रहे हैं, report। और यह है। यह खुल गया है और यह कहता है कि batch में तीन greetings हैं, जो हम उम्मीद करते थे, इसलिए यह perfect है।

अगर मैं यहाँ work डायरेक्टरी में जाता हूँ बस तुम्हें prove करने के लिए कि Nextflow, code में execute किया गया था न कि bash script में, मैं cat work/ command.sh में जा सकता हूँ, और तुम यहाँ देखोगे कि यह बस इस string को सीधे echo कर रहा है। There were three greetings in this batch, और इसलिए वह variable Nextflow द्वारा interpolate किया गया था। यह .command.sh फ़ाइल लिखने से पहले script block में calculate किया गया था। तो resulting variable calculation मूल रूप से इसमें hard coded है इससे पहले कि यह इस मामले में तुम्हारे compute environment पर execute किया जाए।

और इसलिए तुम script के बीच उस separation को देख सकते हो। Block यहाँ और इसके ऊपर कुछ भी। मुझे \_hope है कि यह समझ में आता है।

## सारांश और क्विज़

ठीक है, यह Hello Nextflow के इस भाग का अंत है। तो पहले की तरह, जाओ और quiz check करो। इसे webpage पर या CLI में करो, कुछ सवालों से गुज़रो और बस check करो कि तुमने हमने जो सामग्री cover की है उसे समझ लिया है। देखो कि क्या वहाँ कुछ है जो किसी चीज़ को highlight करता है जो तुम शायद नहीं समझे होगे। बहुत सारे सवाल नहीं। करने में अच्छा और आसान। या तुम यहाँ नीचे webpage पर भी कर सकते हो।

और थोड़ा break लो, थोड़ा घूमो और वापस आओ और Hello, Nextflow के भाग चार में हमसे जुड़ो, जहाँ हम modules के बारे में बात करेंगे। बहुत-बहुत धन्यवाद।
