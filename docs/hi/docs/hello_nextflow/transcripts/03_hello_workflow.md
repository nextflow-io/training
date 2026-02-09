# भाग 3: Hello Workflow - वीडियो ट्रांसक्रिप्ट

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "महत्वपूर्ण नोट"

    यह पेज केवल ट्रांसक्रिप्ट दिखाता है। पूर्ण चरण-दर-चरण निर्देशों के लिए, [कोर्स सामग्री](../03_hello_workflow.md) पर वापस जाओ।

    ट्रांसक्रिप्ट में दिखाए गए सेक्शन नंबर केवल संकेत के उद्देश्य से दिए गए हैं और हो सकता है कि सामग्री में सभी सेक्शन नंबर शामिल न हों।

## स्वागत और पुनरावलोकन

नमस्ते, और Hello Nextflow के भाग तीन में वापस स्वागत है। इस भाग को Hello Workflow कहा जाता है, और यह कोर्स का वह हिस्सा है जहां हम वास्तव में pipeline या workflow नाम को सही ठहराना शुरू करते हैं।

हम अपनी सरल pipeline स्क्रिप्ट को, जिसमें अब तक एक process है, लेंगे और इसमें अतिरिक्त processes जोड़ना शुरू करेंगे और देखेंगे कि Nextflow इस orchestration और pipeline के माध्यम से data flow को कैसे संभालता है।

चलो अपने code spaces पर वापस चलते हैं। तुम देखोगे कि मैंने अपनी सभी .nextflow\* डायरेक्टरीज़ और work डायरेक्टरीज़ और सब कुछ डिलीट कर दिया है ताकि इसे साफ रखा जा सके। चिंता मत करो अगर तुम्हारे पास अभी भी कोर्स के पिछले भागों से वे फ़ाइलें मौजूद हैं।

हम hello-workflow.nf नाम की फ़ाइल से काम करने वाले हैं। पहले की तरह, यह मूल रूप से उस स्क्रिप्ट को दर्शाता है जिसे हमने अब तक बनाया है, और हमें एक साफ शुरुआती जगह देता है। और फिर से, output में नीचे हम देख सकते हैं कि path अब hello_workflow है। तो published फ़ाइलें तुम्हारे results फ़ोल्डर में एक अलग subdirectory में जानी चाहिए।

पुनरावलोकन करने के लिए कि हम अब तक कहां हैं, हमारे पास यहां एक single process है, जिसमें एक input greeting, एक output greeting फ़ाइल है। और फिर simple Bash स्क्रिप्ट, जो बस एक फ़ाइल में echo command करती है।

हमारे पास एक single workflow input है, यहां params block, जहां हम कहते हैं कि यह एक path की उम्मीद कर रहा है, और default data/greetings.csv है, जो यहां ऊपर यह फ़ाइल है।

फिर workflow में ही, हमारे पास एक main block है। हम एक channel बना रहे हैं। हम CSV को rows में parse कर रहे हैं और फिर प्रत्येक array का पहला element ले रहे हैं, और हम उस channel को उस process में pass कर रहे हैं, जो फिर तीन tasks generate कर रहा है, और हम workflow से, उस process के outputs को publish कर रहे हैं।

और फिर अंत में, output block में, हम Nextflow को बता रहे हैं कि इन फ़ाइलों को इस channel से hello_workflow नाम की directory में publish करे। और इन फ़ाइलों को soft link करने के बजाय copy करे।

## 1. workflow में दूसरा चरण जोड़ें

ठीक है, इस भाग में हम अपने workflow में एक दूसरा process जोड़ने वाले हैं। हम sayHello process के outputs लेंगे, और उन्हें दूसरे चरण में process करेंगे, जो उन फ़ाइलों के भीतर सभी अक्षरों को convertToUppercase करने वाला है।

यह बस एक मूर्खतापूर्ण उदाहरण है, यह फिर से कुछ सरल string processing है, लेकिन यह तुम्हें दिखाता है कि हम workflow के भीतर logic को कैसे ले सकते हैं।

हम इसके लिए "tr" नाम की एक bash command करने वाले हैं, जो translate के लिए short है। यह एक Unix command है जो हमेशा से रही है। अगर तुम इससे परिचित नहीं हो, तो मैं तुम्हें दोष नहीं देता। मुझे नहीं लगता कि मैंने training से पहले कभी इसका उपयोग किया था, लेकिन तुम इसे terminal पर बहुत जल्दी आज़मा सकते हो। अगर मैं "echo 'hello world'" करता हूं और फिर 'tr' में pipe करता हूं और फिर quotes में तुम character range कहते हो, तो A से Z, lowercase, और फिर तुम A से Z uppercase करना चाहते हो। और यह बस कहता है, इन अक्षरों को इन अक्षरों में translate करो।

और जब मैं enter दबाता हूं, तुम देख सकते हो कि अब यह सब कुछ capitalize कर दिया है। बहुत अच्छा अगर तुम्हें लोगों पर चिल्लाना पसंद है।

तो यह bash command की एक बहुत ही सरल शैली है जिसे हम अपने दूसरे process में उपयोग करने वाले हैं।

## 1.2. uppercasing चरण को Nextflow process के रूप में लिखें

तो अगर मैं अपनी स्क्रिप्ट पर वापस जाता हूं, तो मैं थोड़ा धोखा देने वाला हूं और बस training के docs से code copy करूंगा। लेकिन तुम बिल्कुल देख सकते हो कि क्या हो रहा है।

हमारे पास यहां एक नया process है। इसे हमने convertToUpper कहा है, लेकिन हम इसे जो चाहें कह सकते हैं।

हमारे पास एक single input path है, जैसा कि हमने पहले किया था। यह value channel नहीं है, यह path channel है। और फिर एक single output।

script block में हम input फ़ाइल पर "cat" करते हैं। और हम इसे squiggly brackets में रख सकते हैं अगर हम चाहें। और जो उस variable को लेता है। और हम pipe में वही bash command चलाते हैं और हम परिणामों को इस फ़ाइल नाम के साथ एक फ़ाइल में लिखते हैं, और वह output path द्वारा pick up किया जाता है।

अब हमें इस नए process के साथ कुछ करने की ज़रूरत है। तो हम workflow में नीचे जाने वाले हैं जहां हम workflow के विभिन्न logic को build करते हैं, और उस पहले process के बाद, हम अपना दूसरा process चलाने वाले हैं। तो convertToUpper यहां process का नाम है।

यह एक input लेता है इसलिए हम इसे अपने आप नहीं बुला सकते। हम पहले process के output को process करना चाहते हैं। तो जैसे हमने इसके साथ किया, sayHello out जहां हम उन परिणामों को publish कर रहे हैं। हम यहां input के रूप में उन्हीं परिणामों का उपयोग करना चाहते हैं, इसलिए हम उन्हें copy कर सकते हैं और उन्हें वहां रख सकते हैं।

हम sayHello process ".out" चाहते हैं, और Nextflow जानता है कि इसका मतलब यहां एक simple single output record है, जो यह फ़ाइल है। तो वह फिर दूसरे process में input के रूप में pass किया जाएगा।

## 1.5. workflow output publishing सेट करें

ठीक है। और अंत में, ताकि हम वास्तव में इस दूसरे process के परिणामों को save करें, हमें उन्हें workflow से भी publish करने की ज़रूरत है, और फिर उन्हें output block में define करें, पहले जैसा ही syntax। तो हम इसे copy कर सकते हैं और second outputs कह सकते हैं, या जो भी तुम इसे कहना चाहते हो।

process का नाम लो जिसमें हम रुचि रखते हैं, convertToUpper out, और फिर यहां नीचे output block में। इसे add करो और हम यहां वही attributes कर सकते हैं। तो हम यह भी चाहते हैं कि ये फ़ाइलें Hello Workflow subdirectory में हों, और हम इन्हें भी copy करना चाहते हैं।

बढ़िया। चलो इसे चलाने की कोशिश करते हैं। तो अगर मैं terminal लाता हूं और मैं "nextflow run hello-workflow.nf" करता हूं, और हम देखेंगे कि यह क्या करता है। देखते हैं कि क्या यह पिछले भागों से अलग दिखता है।

तो यह Nextflow launch करता है। docs में, यह "-resume" के साथ करने के लिए कहता है, लेकिन मैंने अपनी सभी work directory delete कर दी, इसलिए यहां कोई फर्क नहीं पड़ता। लेकिन अगर तुमने किया, तो वह भी काम करेगा।

और यह लगभग बिल्कुल वैसा ही दिखता है। लेकिन तुम देख सकते हो कि अब यहां output की एक दूसरी line है, जहां तुम उस दूसरे process का नाम देख सकते हो जिसे हमने अभी add किया है। और निश्चित रूप से, तुम देख सकते हो कि यह तीन बार सफलतापूर्वक चला।

शानदार। अगर मेरे पास मेरी पिछली work directories होतीं और मैंने इसे "-resume" के साथ किया होता, तो ये pipeline में सिर्फ पहले चरण को cached किया होता। क्योंकि वे outputs बिल्कुल वैसे ही थे, इसलिए Nextflow जानता होता कि उन्हें फिर से reuse करना है।

और तो तुम देख सकते हो कि तुम अपने workflow को चरण-दर-चरण iteratively build करने के लिए -resume का उपयोग कैसे कर सकते हो, अगर तुम्हें ज़रूरत हो।

ठीक है, चलो results directory में यहां ऊपर देखते हैं और देखते हैं कि क्या यह काम किया है। हम देख सकते हैं कि हमारे पास यहां ऊपर कुछ और फ़ाइलें हैं। हमारे पास पहले process से हमारी original फ़ाइलें हैं जैसी पहले थीं। और निश्चित रूप से, हमारे पास हमारी upper फ़ाइलें हैं और अक्षर सभी uppercase हैं, इसलिए यह काम किया है। यह देखना वास्तव में अच्छा है।

इन work directories के अंदर जांचना भी दिलचस्प है। पहले की तरह hash यहां work directories से correspond करता है। तो अगर मैं "ls work" में देखता हूं, और फिर उसे expand करता हूं, तो हम यहां विभिन्न फ़ाइलें देखेंगे।

हम पहले process से output फ़ाइल देखते हैं, जिसे यहां input के रूप में pull किया गया है। और हम नई output फ़ाइल देख सकते हैं जो generate की गई थी।

अब अगर मैं इसे "-la" के साथ करता हूं ताकि list करूं और सभी फ़ाइलें दिखाऊं, तो हम कुछ और चीजें देखेंगे। सबसे पहले, तुम देखोगे कि यह फ़ाइल वास्तव में पहले process के लिए एक soft link है। यह मूल रूप से हमेशा एक soft link होता है अगर यह हो सकता है, फ़ाइल space बचाने के लिए। हम यहां फ़ाइलों को publish नहीं कर रहे हैं और यह बस उस फ़ाइल को पहले task से दूसरे task में reference करता है ताकि सब कुछ उस एक working directory के भीतर encapsulated हो, और बाकी सब चीज़ों से safe और isolated हो।

और वह वहां होना ज़रूरी है क्योंकि अगर हम .command.sh फ़ाइल को देखें, तो अगर मैं "cat work/b8/56\*" करता हूं, तो तुम देख सकते हो कि यहां फ़ाइल parts relative हैं, इसलिए यह उस input फ़ाइल को cat कर रहा है, जिसे उसी working directory में soft link किया गया है।

तो हर work directory ऐसी ही दिखेगी। जब तुम इसे Nextflow में देखते हो, तो तुम्हारे पास सभी input फ़ाइलें होंगी जो उस work directory में staged हैं। और फिर तुम्हारे पास कोई भी output फ़ाइलें होंगी जो create की गई थीं। तो यह बढ़िया है। यह वैसा ही दिखता है जैसा हम उम्मीद करते हैं।

## 2.1. collection command को define करें और terminal में test करें

ठीक है, चलो अपने workflow पर वापस चलते हैं। अगला चरण क्या है जो हम करना चाहते हैं?

हमारे पास अब दो processes हैं और वे इस एक CSV फ़ाइल को ले रहे हैं, इसे parse कर रहे हैं और split कर रहे हैं। और फिर हमारे पास इन processes में से प्रत्येक के लिए तीन tasks हैं और Nextflow इस सब के parallelization को संभालता है, इसलिए यह सब side by side चलता है जहां संभव हो।

काम को parallel में चलाने के लिए split करने का वह तरीका बहुत आम है। और इसका inverse फिर सब कुछ वापस gather करना है। तो हम अपने workflow में अपने अंतिम process के साथ यही करने वाले हैं, हमारे पास यहां एक तीसरा होगा, जो इन तीन अलग-अलग outputs को लेता है और उन सभी को एक single फ़ाइल में combine करता है।

हम इसे terminal में काफी सरलता से कर सकते हैं, बस यह महसूस करने के लिए कि यह कैसा दिखेगा।

अगर मैं results फ़ोल्डर में जाता हूं। तो, "cd results/hello_workflow/", और हमारे पास यहां सभी UPPER फ़ाइलें हैं। मैं बस "cat" का उपयोग कर सकता हूं, जिसका उपयोग हम उस फ़ाइल की contents को print करने के लिए करते हैं, और तुम "cat" को कई फ़ाइलें दे सकते हो और यह एक के बाद एक read करेगा।

तो मैं "UPPER-\*" कह सकता हूं, जो मुझे Bash expansion के साथ तीन फ़ाइल नामों की वही list देता है। और मैं combined.txt कह सकता हूं। मुझे लगता है कि docs में, यह exact फ़ाइल नाम list करता है, लेकिन यह वही काम कर रहा है।

अब, अगर मैं "cat combined.txt" का उपयोग करता हूं, तो हम देख सकते हैं कि हमारे पास उन तीनों फ़ाइलों की फ़ाइल contents हैं।

तो मूल रूप से यह process बस यही करने वाला है कि हम इसे एक single process task में पिछले process से सभी अलग-अलग output फ़ाइलें देने की कोशिश करने वाले हैं, और फिर हम उन्हें "cat" करने वाले हैं और output फ़ाइल save करने वाले हैं।

## 2.2. collection चरण करने के लिए एक नया process बनाएं

ठीक है, तो चलो अपना नया process add करते हैं। मैं इसे training materials से paste करने वाला हूं, और तुम देख सकते हो कि इसने हमें यहां इन question marks के साथ reader के लिए थोड़ा exercise छोड़ा है। लेकिन तुम process की general outline देख सकते हो जो मूल रूप से वही है जो हमने अभी terminal में किया था, जहां हम बहुत सारी input फ़ाइलों का "cat" कर रहे हैं और इसे यहां collected नाम की output फ़ाइल में लिख रहे हैं, और फिर output फिर से उस single path की उम्मीद करता है।

तो हमें यहां किसी प्रकार के input की ज़रूरत है और वे paths का एक set होने वाले हैं। तो फिर से, हम एक input path channel define करते हैं और चलो इसे input_files कहते हैं। अब, इसने पहले हमें यहां एक single path दिया है, लेकिन एक path में यहां कई फ़ाइलें भी हो सकती हैं, भले ही यह अभी भी एक single declaration हो।

मैं इसे यहां नीचे copy करने वाला हूं क्योंकि हम इन फ़ाइलों को "cat" करना चाहते हैं। और तुम सोच सकते हो कि हमें यहां array print करने या इस तरह की चीज़ों के साथ कुछ समस्याएं हैं, लेकिन Nextflow आम तौर पर इस मामले में काफी sensible है। और अगर इसे इस तरह कई फ़ाइलों के साथ एक channel दिया जाता है, तो यह उन सभी को space separators के साथ एक साथ रखेगा। तो यह हमें सही syntax देगा।

यह बढ़िया है। तो अब चलो अपने नए process को wire up करते हैं। मैं workflow में नीचे जाता हूं। मैं outputs को combine करने वाला हूं, नया process नाम, और पहले की तरह ही। मैं इस पिछले process को लेने वाला हूं, convertToUpper और ".out" करूंगा।

बढ़िया। चलो इसे आज़माते हैं और देखते हैं कि क्या यह terminal में काम करता है। अगर मैं बस कुछ directories ऊपर वापस जाता हूं और फिर Nextflow command फिर से चलाता हूं, और हम देखेंगे कि क्या होता है।

तो workflow launch हो गया है और अब तुम देख सकते हो कि हमारे पास तीन अलग-अलग process नाम हैं, जो बढ़िया है। पहले दोनों पहले जैसे ही दिखते हैं, और तीसरा नया चलता है, जो अच्छा है।

हालांकि, यहां कुछ अजीब है। हम उन output फ़ाइलों को एक single फ़ाइल में combine करना चाहते थे, और फिर भी इस process को हम देख सकते हैं कि तीन बार चला है, एक बार नहीं।

निश्चित रूप से, अगर हम इन work directories में से एक में जाते हैं। और "cat work/" "collected" करते हैं, तो हम देखेंगे। यहां केवल एक single word है, तीन नहीं।

और तो जो हुआ है वह यह है कि Nextflow ने उस parallelization को जारी रखा है जैसा कि उसने पिछले चरणों में किया था। और इस process ने हमें तीन elements के साथ एक channel दिया, और वे तीन channel elements हमारे downstream process में pass किए गए, जिसने तीन process tasks generate किए।

इसने मूल रूप से तीन अलग-अलग बार collect करने की कोशिश की और हर बार इसके पास बस एक single फ़ाइल थी, इसलिए इसने बस cat single फ़ाइल को एक output में किया, और वास्तव में, हम इसे .command.sh फ़ाइल में भी देख सकते हैं।

अगर मैं .command.sh करता हूं, तो हम देख सकते हैं कि इसमें यहां बस एक single फ़ाइल नाम है और केवल एक single फ़ाइल उस working directory में staged की गई थी।

## 2.3. workflow में collection चरण add करें

तो किसी तरह हमें Nextflow को बताने की ज़रूरत है कि उन सभी outputs को पिछले process से एक साथ लाए और उन्हें इस downstream process को एक single channel element के रूप में दे, तीन के बजाय।

हम यह _collect_ नाम के एक channel operator के साथ करते हैं।

यह एक बेहद उपयोगी operator है, जिसे तुम Nextflow pipelines में हर समय देखोगे। यह यहां एक channel है, यह output channel, बिल्कुल वैसा ही जैसा हमने ऊपर बनाया था। और इसलिए हम इसमें channel operators append कर सकते हैं जैसा कि हमने पहले किया था। हम बस dot कर सकते हैं, और फिर इस मामले में, collect, brackets।

और बस इतना ही हमें चाहिए। यह फिर इस channel को manipulate करने वाला है इससे पहले कि यह इस process में pass किया जाए।

अगर तुम देखना चाहते हो कि इसके साथ क्या हो रहा है, तो हम इसे यहां भी view कर सकते हैं। तो यहां, यह इस process को चलाने से बिल्कुल भी संबंधित नहीं है, इसलिए मैं इसे उस process को चलाने के बाद किसी भी बिंदु पर रख सकता था। लेकिन हम वही output channel लेते हैं, और हम इसे .view के साथ देख रहे हैं, और फिर हम इसे फिर से .collect.view के साथ देख रहे हैं।

और जब हम इसे चलाते हैं, तो यह हमें उस channel की दो अलग-अलग structures दिखाएगा, collect से पहले और बाद में। तो चलो अब इसे आज़माते हैं। ठीक है, मैंने थोड़ा zoom out किया है क्योंकि कुछ outputs काफी लंबे हैं, लेकिन अगर मैं pipeline चलाता हूं, तो हम देखेंगे कि क्या यह काम करता है।

मुझे उम्मीद है कि तीसरा process सिर्फ एक बार चलेगा, क्योंकि यह outputs को collect कर रहा है और निश्चित रूप से, तुम देख सकते हो कि collectGreetings एक में से एक के रूप में। तो वह सिर्फ एक task चला।

और फिर अगर हम view statements को देखें, तो हमारे पास before के तीन elements के लिए तीन view statements हैं, प्रत्येक में एक फ़ाइल path के साथ।

और फिर उस collect statement के बाद, वह सिर्फ एक बार trigger हुआ क्योंकि उस channel में एक single element है। और अब हमारे पास तीन अलग-अलग फ़ाइल paths की यह list है।

यह बिल्कुल वही है जिसकी हमें उम्मीद थी। और तुम उम्मीद कर सकते हो, यह मूल रूप से उस "map" operator का inverse है जो हमने CSV arrays से अलग channel elements में जाने के लिए किया था। अब हम अलग channel elements ले रहे हैं और उन्हें वापस एक single array में डाल रहे हैं।

बढ़िया, हम इन view statements को clear कर सकते हैं। हमें अब इनकी ज़रूरत नहीं है। हम अगले चरण पर जा सकते हैं।

इससे पहले कि मैं आगे बढ़ूं, और इससे पहले कि मैं भूल जाऊं, मैं यहां एक नया publish statement add करने वाला हूं। Third output। तुम इसे अपने workflow में कुछ अधिक semantic और descriptive कह सकते हो। और फिर मैं इसे फिर से output block में add करने वाला हूं और कहूंगा path 'hello_workflow' mode 'copy'। बस ताकि इस process द्वारा generate की गई output फ़ाइल यहां ऊपर हमारे results फ़ोल्डर में save हो जाए।

बस जल्दी से double check करने के लिए कि यह काम करता है। अब थोड़ा cleaner होना चाहिए क्योंकि हमारे पास वे view statements नहीं हैं। और, हम देखेंगे कि क्या हमें यहां ऊपर हमारी नई output फ़ाइल मिलती है। एक में से, एक task चला, collected नाम की एक नई फ़ाइल मिली, और अब हमारे पास वे तीनों words हैं। शानदार। आगे क्या है?

## 3. process में अतिरिक्त parameters pass करें

ठीक है। अगला हम एक single process में कई inputs को handle करने को देखने वाले हैं। अब तक तुम देख सकते हो कि हमारे सभी processes बस एक चीज़ को input के रूप में ले रहे हैं। उन सभी के input के तहत एक single line है।

हम इसे Nextflow को एक अलग batch identifier specify करने की अनुमति देकर demonstrate करने वाले हैं ताकि शायद तुम इस workflow को कई बार चलाओ और तुम इसे हर बार एक अलग batch ID दे सकते हो।

मैं बस collectGreetings के लिए input में यहां एक दूसरी line add करने वाला हूं। और मैं इसे "val" कहने वाला हूं, क्योंकि यह एक string है। अब यह एक value है, path नहीं, और मैं इसे "batch_name" कहने वाला हूं।

फिर मैं इस variable का उपयोग करने के लिए यहां नीचे script को edit करने वाला हूं, और मैं इसे उसी जगह पर रखने की कोशिश करने वाला हूं जहां training material में है। तो मैं इसे इस फ़ाइल path COLLECTED-$\{batch_name\}-output के बीच में रखता हूं।

अभी पूरा नहीं हुआ। याद रखो कि हमें Nextflow को बताना होगा कि output फ़ाइल नाम क्या होने वाले हैं। तो हमें यहां ऊपर भी वही काम करना होगा: COLLECTED-$\{batch_name\}-output.txt"।

शानदार। Nextflow अब एक दूसरा variable input प्राप्त कर रहा है और यह उसे script और output में interpolate कर रहा है।

एक आखिरी चीज़, अब हमें यह पता लगाना होगा कि इसे कहां call किया जा रहा है, और हमें process में दूसरा input pass करना होगा। यह किसी भी अन्य भाषा में किसी भी अन्य function में input की तरह ही है।

जैसा कि हमने training में पहले किया था, मैं यहां special "params" का उपयोग करने वाला हूं, और हम इसे "params.batch" कहने वाले हैं ताकि हमारे पास -- batch CLI option हो सके। और अब तुम देख सकते हो कि हमारे process में यहां दो अलग-अलग inputs हैं बस comma separated, जो pass किए जा रहे हैं।

order को सही करना वास्तव में महत्वपूर्ण है, इसलिए channel के लिए यहां arguments का order और फिर param match होना चाहिए। channel और फिर वहां batch name। यह बस positional matching है।

ठीक है। मैं अब --batch के साथ सीधे इस pipeline को चला सकता हूं, लेकिन चलो पहले सही काम करते हैं और इसे यहां Params में input में define करते हैं। तो मैं इसे batch में add करने वाला हूं और फिर हम कहेंगे कि यह एक string है और चलो इसे एक default देते हैं। तो चलो इसे बस batch कहते हैं। ठीक है? अब चलो workflow चलाने की कोशिश करते हैं।

--batch Trio। मुझे लगता है कि यह training material में कहता है, लेकिन हम वहां कोई भी string उपयोग कर सकते हैं जो हम चाहें। और उम्मीद है कि हम उस results output फ़ाइल को यहां ऊपर आते हुए देखेंगे।

और निश्चित रूप से, COLLECTED-trio-output - वह ठीक से काम किया है। इसने हमारी फ़ाइल का नाम बदल दिया है। और तुम अब कल्पना कर सकते हो कि यह उपयोगी है क्योंकि अगर मैं इसे फिर से एक अलग batch name के साथ चलाता हूं, जैसे replicate_two, तो यह हमें यहां ऊपर एक अलग batch name देने वाला है।

और और यह फिर इस मामले में output फ़ाइलों को clobber नहीं करेगा। तो यह अच्छा है।

## 4. collector चरण में एक output add करें

ठीक है, तो अब हमारे पास यहां हमारे process में कई inputs हैं। लेकिन क्या होता है अगर हम कई outputs create करना चाहते हैं? यहां हमारा उदाहरण फिर यह है कि हम इस process के लिए एक report create करने वाले हैं, बस यह कहते हुए कि कितनी फ़ाइलें collect की गईं।

और हम यह यहां एक echo command के साथ करेंगे। तो हम echo कह सकते हैं। There were, मैं इसे training material से copy करने वाला हूं, ताकि तुम्हें मुझे इसे type करते हुए देखना न पड़े।

There were $\{count_greetings\} greetings in this batch, और इसे अब $\{batch_name\} नाम की एक नई फ़ाइल में save करें, तो वही variable, हम इसे जितनी बार चाहें उतनी बार reuse कर सकते हैं, report.txt।

## 4.1.1. collect की गई greetings की संख्या count करें

हमें वास्तव में किसी तरह से इसकी गणना करने की ज़रूरत है। हम Bash logic का उपयोग करके Bash script में वह logic कर सकते हैं अगर हम चाहें। हालांकि, हम सीधे Nextflow code के भीतर भी scripting कर सकते हैं, जब तक कि यह process में script block के भीतर हो और quoted section के ऊपर हो।

यहां कुछ भी final rendered script में शामिल नहीं होगा, और यह बस Nextflow द्वारा execute किया जाएगा जब यह एक task render करता है।

तो यहां हम बस कुछ logic कर रहे हैं। हम count_greetings नाम का एक नया variable बना रहे हैं। हम यहां input files channel लेते हैं, और हम इस पर .size() call कर रहे हैं।

ठीक है, वह function मुझे इस variable में एक number देने वाला है, और अब हमारी warning चली गई है क्योंकि यह variable define किया जा रहा है।

ठीक है, तो हम work directory में वह दूसरी फ़ाइल create कर रहे हैं, लेकिन हमें Nextflow को बताने की ज़रूरत है कि इसे इस process के published output के रूप में expect करे। तो हम यह बिल्कुल उसी syntax से करते हैं जैसा हमने पहली फ़ाइल के लिए किया था।

हम path कहते हैं क्योंकि यह, फिर से, हम यहां एक variable publish कर सकते हैं अगर हम "val" के साथ चाहें, लेकिन हम "path" कहने वाले हैं। और फिर expected फ़ाइल नाम। ध्यान दो कि यह यहां highlighted नहीं है। ऐसा इसलिए है क्योंकि मैंने single quotes का उपयोग किया। मुझे double quotes का उपयोग करना होगा।

## 4.1.2. report फ़ाइल को emit करें और outputs को नाम दें

ठीक है, यह बढ़िया है। और अब हम यहां नीचे इन outputs को access कर सकते हैं जैसा मैंने यहां किया। लेकिन अब यह अलग-अलग objects का एक array है, इसलिए मैं पहला पाने के लिए collectGreetings.out[0] कर सकता था, या दूसरा पाने के लिए one, जो हमारी नई report है।

लेकिन मुझे वास्तव में ऐसा करना पसंद नहीं है क्योंकि index counting को mess up करना काफी आसान है। और तुम वहां बहुत सारी lines count करते हुए बैठते हो और तुम एक नया output add करते हो और अचानक सब कुछ break हो जाता है। तो

इसके बजाय नाम से सब कुछ reference करना बहुत अच्छा है। और हम यहां "emit" नाम की एक special key के साथ ऐसा कर सकते हैं।

तो हम इसे जो चाहें कह सकते हैं। चलो इसे emit outfile कहते हैं, और emit reports। अगर तुम इन्हें define करते हो और तुम इसे एक या कई पर कर सकते हो, यह तुम पर निर्भर है। अब मैं यहां नीचे जा सकता हूं और इसके बजाय मैं dot out dot reports जा सकता हूं और बस इसे नाम से call कर सकता हूं, जो तुम्हारे code को पढ़ते समय समझना बहुत आसान है, और यह code में changes के लिए safer है।

मैंने यहां .out.report add किया है, लेकिन वास्तव में मुझे दो अलग-अलग outputs publish करने की ज़रूरत है। तो मैं इसे कुछ अधिक interesting के रूप में rename करने वाला हूं जैसे collected और report और क्या मैंने इसे यही कहा था? मैंने इसे out file कहा, sorry। तो वह emit नाम यहां outfile और report। क्योंकि हम दो अलग-अलग output channels publish कर रहे हैं और इसलिए हमें publish block में उन दोनों को reference करने की ज़रूरत है।

फिर हमें इन्हें output block में भी define करने की ज़रूरत है। तो मैंने उसे collected rename किया, और फिर से, reports के लिए, यहां थोड़ा verbose, लेकिन यह वास्तव में उपयोगी है जब तुम एक नए workflow को पढ़ने के लिए आते हो, तो यहां सभी अलग-अलग outputs, सभी अलग-अलग channels को side by side listed देखने के लिए, और इसे कम verbose बनाने के तरीके हैं, जिन पर हम बाद में touch करेंगे।

ठीक है, चलो इसे आज़माते हैं और अपना workflow चलाते हैं और देखते हैं कि क्या होता है।

उम्मीद है कि अब यह मूल रूप से पहले की तरह ही चलना चाहिए। और हमें यहां ऊपर replicate_two नाम की एक नई output फ़ाइल मिलने वाली है, report। और वहां यह है। यह खुल गया है और यह कहता है कि batch में तीन greetings हैं, जो हमने उम्मीद की थी, इसलिए यह perfect है।

अगर मैं यहां work directory में जाता हूं बस तुम्हें साबित करने के लिए कि Nextflow code में execute किया गया था न कि bash script में, मैं cat work/ command.sh में जा सकता हूं, और तुम यहां देखोगे कि यह बस इस string को सीधे echo कर रहा है। There were three greetings in this batch, और इसलिए वह variable Nextflow द्वारा interpolate किया गया था। यह .command.sh फ़ाइल लिखने से पहले script block में calculate किया गया था। तो resulting variable calculation मूल रूप से इसमें hard coded है इससे पहले कि यह इस मामले में तुम्हारे compute environment पर execute किया जाए।

और इसलिए तुम script के बीच उस separation को देख सकते हो। यहां Block और इसके ऊपर कुछ भी। मुझे उम्मीद है कि यह समझ में आता है।

## सारांश और quiz

ठीक है, यह Hello Nextflow के इस भाग का अंत है। तो पहले की तरह, जाओ और quiz check करो। इसे webpage पर करो या CLI में, कुछ सवालों से गुज़रो और बस check करो कि तुमने हमने जो सामग्री cover की है उसमें से कुछ को समझा है। देखो कि क्या वहां कुछ ऐसा है जो कुछ ऐसा highlight करता है जो शायद तुमने नहीं समझा। बहुत सारे सवाल नहीं। करना अच्छा और आसान। या तुम इसे यहां नीचे webpage पर भी कर सकते हो।

और थोड़ा break लो, थोड़ा घूमो और वापस आओ और भाग चार में हमसे जुड़ो Hello, Nextflow का, जहां हम modules के बारे में बात करेंगे। बहुत-बहुत धन्यवाद।
