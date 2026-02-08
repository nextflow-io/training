# भाग 2: Hello Channels - वीडियो ट्रांसक्रिप्ट

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "महत्वपूर्ण नोट"

    यह पेज केवल ट्रांसक्रिप्ट दिखाता है। पूर्ण चरण-दर-चरण निर्देशों के लिए, [कोर्स मटेरियल](../02_hello_channels.md) पर वापस जाओ।

    ट्रांसक्रिप्ट में दिखाए गए सेक्शन नंबर केवल संकेतात्मक उद्देश्यों के लिए प्रदान किए गए हैं और इसमें मटेरियल के सभी सेक्शन नंबर शामिल नहीं हो सकते हैं।

## स्वागत

हैलो और Hello Nextflow के भाग 2 में वापस स्वागत है। इस चैप्टर को Hello Channels कहा जाता है।

Channels तुम्हारी Nextflow पाइपलाइन में गोंद की तरह हैं। ये वो बिट्स हैं जो सभी अलग-अलग processes को एक साथ रखते हैं, जिनका Nextflow सभी जानकारी पास करने और तुम्हारी workflow को orchestrate करने के लिए उपयोग करता है।

Channels का एक और हिस्सा है operators। ये मूल रूप से functions हैं जिनका हम channels पर contents को modify करने के लिए उपयोग कर सकते हैं। चलो VS code में dive करें और देखें कि हम कहां हैं।

मैं इस VS code पर बहुत zoomed in हूँ, इसलिए चीजों को साफ-सुथरा रखने के लिए, मैंने सभी _.nextflow\*_ फ़ाइलें और _work/_ डायरेक्टरी और results/ और Chapter One की हर चीज़ हटा दी है। और मैं यहां fresh शुरू कर रहा हूँ। लेकिन इसके बारे में ज्यादा चिंता मत करो। अगर तुम नहीं चाहते हो, तो उन फ़ाइलों को छोड़ सकते हो। वे कोई समस्या नहीं पैदा करेंगी।

हम इस चैप्टर के लिए _hello-channels.nf_ पर काम करना शुरू करने वाले हैं, और अगर मैं इसे खोलता हूँ, तो यह उस फ़ाइल के बहुत समान दिखनी चाहिए जिस पर हम पहले काम कर रहे थे। हो सकता है कि अलग-अलग हिस्से स्क्रिप्ट के अलग-अलग भागों में हों, लेकिन सब कुछ मूल रूप से समान होना चाहिए।

एक चीज़ जो अलग है वह यह है कि यहां output block में path अब _hello_channels_ है इस भाग के लिए, जिसका मतलब है कि result फ़ाइलें तुम्हारे results में एक अलग subdirectory में store होंगी अगर तुम्हारे पास अभी भी वो है। तो outputs के बारे में confused हुए बिना शुरू करने के लिए यह एक अच्छी और साफ जगह होनी चाहिए।

ठीक है, तो चलो जल्दी से याद करें कि जब हम इस workflow को रन करते हैं तो यह स्क्रिप्ट क्या करती है। हम _"nextflow run hello-channels.nf"_ करते हैं। हम _"--input myinput"_ कर सकते हैं, और जब हम इसे रन करते हैं, तो यह इस parameter, params.input का उपयोग करने वाला है, जो ऊपर sayHello process के लिए variable के रूप में pass किया गया था, जो greeting में जाता है और output.txt में save हो जाता है। और हम इसे results फ़ाइल में देख सकते हैं। बढ़िया।

## 1. चैनल के माध्यम से स्पष्ट रूप से variable inputs प्रदान करना

यह अच्छा है। लेकिन यह काफी सरल है। इस parameter में हमारा एक variable है, जो एक process में जाता है जो एक बार रन होता है, और वास्तव में scale नहीं होता। और हम इसे यहां बनाने के लिए बहुत सारी अलग-अलग फ़ाइलें नहीं दे सकते। हम इसे बहुत सारे अलग-अलग greetings नहीं दे सकते। हमारे पास बस एक है।

वास्तविकता में, Nextflow तुम्हारे analysis को scale up करने के बारे में है। तो शायद तुम चाहोगे कि यह एक से अधिक चीज़ें करे। और हम _channels_ के साथ ऐसा करते हैं।

Channels बहुत से लोगों के लिए एक अनोखी अवधारणा है जो Nextflow को pick up कर रहे हैं। यह functional programming की अवधारणाओं से आती है, और इसे समझने में थोड़ा समय लग सकता है, लेकिन एक बार जब तुम click कर जाते हो, तो वे वास्तव में Nextflow की शक्ति को unlock करते हैं और यह तुम्हारी workflows लिखने के तरीके की कुंजी है।

## 1.1. एक input channel बनाना

चलो इस स्क्रिप्ट को लेकर शुरू करें और इसे केवल एक _param_ के बजाय एक _channel_ का उपयोग करने दें।

हम workflow पर जाते हैं, जो हमारी सभी workflow logic के बारे में है चीजों को एक साथ जोड़ने के बारे में। और मैं यहां अंदर जाने वाला हूँ और मैं एक नया channel बनाने वाला हूँ।

एक नया channel बनाओ।

और मैं इसे "_greeting_ch"_ कहने वाला हूँ। "_\_ch"_ इस तरह करना convention है, बस ताकि तुम याद रख सको कि यह variable एक channel है। लेकिन तुम इसे जो चाहो कह सकते हो।

और फिर मैं equals कहने वाला हूँ, और मैं _"channel.of"_ करने वाला हूँ।

Channel channels से संबंधित सभी चीज़ों के लिए name space की तरह है। छोटा "c" अगर तुम पहले Nextflow का उपयोग कर रहे हो। और _".of"_ कुछ है जिसे Channel factory कहा जाता है, जो मूल रूप से एक channel बनाने का एक तरीका है।

बहुत सारे अलग-अलग channel factories हैं। अगर मैं यहां सिर्फ "." करता हूँ, तो तुम देख सकते हो कि VS Code उनमें से बहुत सारे suggest कर रहा है, लेकिन _".of"_ सबसे सरल है और बस एक input लेता है।

तो मैं कुछ brackets कर सकता हूँ और मैं _"Hello Channels!"_ कहने वाला हूँ।

बढ़िया। मेरे पास एक channel है। शानदार। मैं save hit कर सकता हूँ, मैं इसे फिर से रन कर सकता हूँ, लेकिन कुछ भी interesting होने वाला नहीं है। VS Code ने मुझे यहां एक orange warning line दी है और मुझे बताया है कि यह set up है: तुमने इसे बनाया है, लेकिन तुमने वास्तव में इसका उपयोग किसी भी चीज़ के लिए नहीं किया है। इस channel को consume नहीं किया जा रहा है।

ठीक है, तो हम इसका उपयोग कैसे करें? बहुत आसान। मैं इसे copy करने वाला हूँ, और मैं _params.input_ को delete करने वाला हूँ और मैं इसके बजाय यहां _"greeting_ch"_ डालने वाला हूँ। तो हम इस channel को sayHello के input के रूप में pass करने वाले हैं।

ध्यान दो कि मैंने अभी के लिए इस string को hard code किया है। पिछले चैप्टर के अंत में हमने जो अच्छा param उपयोग किया था, उसके बाद यह थोड़ा backward step है, लेकिन यह सिर्फ चीजों को शुरू करने के लिए सरल रखता है ताकि तुम logic देख सको।

ठीक है, मैं अपने terminal में जाने वाला हूँ और मैं workflow को फिर से रन करने वाला हूँ। इस बार बिना किसी _"--input"_ के, और यह रन होने वाला है और यह उस channel का उपयोग करने वाला है जो हमने बनाया है और उम्मीद है कि हमारे पास यहां _results/hello_channels/_ में एक फ़ाइल होनी चाहिए और अब यह कहती है _"Hello Channels!"_। शानदार। तो यह वही है जो हम उम्मीद कर रहे थे, हमारे channel से। बढ़िया।

## 1.4. channel contents को inspect करने के लिए view() का उपयोग करना

यहां जोड़ने के लिए एक और चीज़, बस channels पर हम जो एक और function उपयोग कर सकते हैं उसका एक त्वरित परिचय जिसे "_.view"_ कहा जाता है।

यह Python या अन्य भाषाओं में _print_ command के समान है जिनका तुम उपयोग कर सकते हो, और यह बस इस channel की contents को terminal पर dump कर देता है जब हम इसे रन करते हैं।

तो "_.view"_ करो, और फिर अगर मैं workflow को फिर से rerun करता हूँ, तो यह terminal पर print करना चाहिए कि उस channel की contents क्या है, उस समय जब हमने इसे बनाया था।

निश्चित रूप से, तुम देख सकते हो कि यह यहां terminal पर print हुआ है। _"Hello Channels!"_।

ध्यान दो कि तुम चाहो तो इन चीजों को lines में break कर सकते हो, और वास्तव में, Nextflow का automatic formatter तुम्हारे लिए ऐसा करने की कोशिश करेगा। White space यहां वास्तव में महत्वपूर्ण नहीं है, तो तुम इन चीजों को एक के बाद एक chain कर सकते हो।

## 2. कई input values पर रन करने के लिए workflow को modify करना

ठीक है, तो हमारे channel में एक चीज़ है जो अच्छी है, लेकिन यह मूल रूप से वैसी ही है जैसी पहले थी। तो चलो इसे थोड़ा और complicated बनाते हैं। चलो अपने channel में कुछ और चीजें add करें।

"_.of()"_ channel factory कई items ले सकता है, तो चलो कुछ और लिखें। हम _Hello, Bonjour, Hej_ करेंगे। और फिर हम इस workflow को फिर से रन कर सकते हैं और देखेंगे कि क्या होता है।

फिर से रन होना चाहिए। और हमने अब print किया है। हमारे view statement के साथ _"Hello", "Bonjour"_ और _"Hej"_ terminal पर। शानदार।

## 2.1.2. command रन करना और log output देखना

तुम सोच सकते हो कि हम इस बिंदु पर हो गए हैं। लेकिन वास्तव में यहां थोड़ा gotcha है, जो हमें trip up करने वाला है। अगर हम अपनी output फ़ाइल यहां देखें। तुम देख सकते हो कि इसमें _"Hello"_ है, लेकिन इसमें कोई अन्य outputs नहीं है। वास्तव में, यह सिर्फ यह एक है।

अगर हम इस workflow को कई बार रन करते हैं, तो हम देख सकते हैं कि कभी-कभी इसमें _"Bonjour"_ है, कभी-कभी इसमें _"Hej"_ है। यह थोड़ा random है।

अगर हम terminal को देखें, तो हम देख सकते हैं कि यह तीन बार रन हुआ और हम अलग-अलग view outputs देख सकते हैं। लेकिन अगर मैं work directory में जाता हूँ, तो मैं _"cat work"_ कर सकता हूँ। इस hash को डालो और उसे expand करो और _output.txt_। तुम देख सकते हो कि work directory में यह फ़ाइल results directory से अलग है, और यह _"Hej"_ है। तो यहां कुछ ठीक से काम नहीं कर रहा है।

और कुंजी यह है कि, हमारे पास तीन tasks रन हुए। Nextflow output processing के दौरान इसे summarize करने की कोशिश करता है, ताकि यह तुम्हारे पूरे terminal को completely take over न करे, और वह ANSI Logging ANSI escape codes का उपयोग करता है, मूल रूप से अन्य tasks को overwrite कर दिया है। तो यह तुम्हें बस अंतिम वाला दिखाता है जो updated होना हुआ।

## 2.1.3. -ansi-log false option के साथ command को फिर से रन करना

कुछ चीजें हैं जो हम वास्तव में इसे थोड़ा बेहतर समझने के लिए कर सकते हैं। हम work directory में ही देख सकते हैं और तुम वहां सभी अलग-अलग work dirs देख सकते हो, लेकिन यह थोड़ा confusing है क्योंकि यह अलग-अलग Nextflow execution runs के साथ mixed up होगा।

या हम Nextflow को ANSI escape codes का उपयोग न करने के लिए कह सकते हैं।

तो अगर मैं command को फिर से रन करता हूँ, लेकिन इस बार मैं इसे बंद करने के लिए _"-ansi-log false"_ कहता हूँ, तो मैं environment variables _$NO_COLOR_ या _"$NXF_ANSI_LOG=false"_ का भी उपयोग कर सकता था। फिर यह बिना किसी clever updates के इन escape codes के बिना Nextflow logging की अधिक पुरानी शैली का उपयोग करता है। यह सिर्फ सीधे terminal पर print करता है।

और अब हम इन सभी तीन processes को देख सकते हैं जो रन हुए। और उनमें से प्रत्येक का अपना task hash। और अगर हम इन work directories में जाते हैं, तो हम तीन अलग-अलग greetings देखेंगे जो हमने specify किए थे।

तो यह अब थोड़ा अधिक समझ में आता है। उम्मीद है कि तुम समझ गए कि Nextflow यह कर रहा था, यह बस उन work directories के साथ तुम्हें terminal में जो दिखाया गया था उसके साथ थोड़ा clever हो रहा था।

हालांकि, इसने work directories के साथ एक समस्या को fix कर दिया है, लेकिन इसने output फ़ाइल के साथ समस्या को fix नहीं किया है। हमारे पास अभी भी केवल एक output फ़ाइल है जो _"Hello"_ कहती है।

## 2.2. सुनिश्चित करना कि output फ़ाइल names unique होंगे

अब इसे समझने के लिए, हमें अपनी workflow script पर वापस जाने की जरूरत है। हम यहां अपना channel generate कर रहे हैं, हम इसे अपनी process में pass कर रहे हैं, और अगर हम process को देखें, तो हम greeting को _"output.txt"_ नाम की एक फ़ाइल में लिख रहे हैं और उस output फ़ाइल को यहां नीचे output block में pass कर रहे हैं, इसे publish कर रहे हैं।

हालांकि, हर बार जब यह process इन तीन अलग-अलग tasks को रन करता है। वे सभी _"output.txt"_ नाम की एक फ़ाइल generate करते हैं, वे सभी output फ़ाइलें results directory में publish की जाती हैं, और वे सभी एक दूसरे को overwrite कर देती हैं। तो जो भी result फ़ाइल तुम्हें वहां मिलती है वह बस अंतिम वाली है जो generate हुई थी, लेकिन बाकी सभी को clobbered कर दिया। यह वास्तव में वह नहीं है जो हम चाहते हैं।

## 2.2.1. एक dynamic output फ़ाइल name construct करना

इसे handle करने के अलग-अलग तरीके हैं, लेकिन अभी के लिए सबसे आसान बस अलग-अलग unique फ़ाइल names बनाना है। तो हर बार जब task एक अलग greeting के साथ रन होता है, तो यह एक अलग output फ़ाइल generate करेगा, जो publish होने पर अब clash नहीं करेगी। और फिर हमें तीन unique output फ़ाइलें मिलेंगी।

हम इसे बिल्कुल उसी तरह करते हैं। हम इस variable का उपयोग script block के भीतर कहीं भी कर सकते हैं और हम इसे कई बार उपयोग कर सकते हैं।

तो मैं इसे यहां paste कर सकता हूँ, _"$\{greeting\}\_output.txt"_, और फिर मुझे इसे यहां भी paste करने की जरूरत है क्योंकि हम अब _output.txt_ नाम की एक फ़ाइल नहीं बना रहे हैं। तो अगर मैं इसे update नहीं करता हूँ, तो Nextflow crash हो जाएगा एक error के साथ यह कहते हुए कि इसने एक फ़ाइल की उम्मीद की थी, जो कभी generate नहीं हुई।

तो मुझे वहां भी वैसा ही करने की जरूरत है और मुझे double quotes का उपयोग करने की जरूरत है, single quotes का नहीं, ताकि इस variable को समझा जाए।

ठीक है, चलो इसे try करें और देखें कि क्या यह काम करता है। हम workflow को फिर से रन करने वाले हैं। उम्मीद है कि यह हमें तीन अलग-अलग work directories के भीतर तीन अलग-अलग tasks दिखाएगा। और निश्चित रूप से, तुम यहां बाईं ओर results folder में देख सकते हो। अब हमारे पास तीन अलग-अलग फ़ाइलें हैं तीन अलग-अलग फ़ाइल names के साथ और प्रत्येक में अलग-अलग contents हैं जो हम उम्मीद करते हैं। तो फ़ाइलें अब एक दूसरे को clobber नहीं कर रही हैं, और सब कुछ वैसा ही है जैसा हम उम्मीद करते हैं।

यह थोड़ा सा trivial setup है जिससे हम यहां गुजरे हैं, लेकिन यह कुछ प्रमुख अवधारणाओं को रेखांकित करता है जिन्हें तुम्हें समझने की जरूरत है कि फ़ाइल publishing कैसे काम करती है, और कुछ चीजें जो तुम traps के रूप में fall कर सकते हो। तो उम्मीद है कि तुम अपनी खुद की workflows में इससे बच सकते हो।

यह भी ध्यान देने योग्य है कि हमने यहां जो किया है वह real life situations में थोड़ा impractical है। हमने कुछ input data लिया है और हम उस data का उपयोग कर रहे हैं, लेकिन हम उस data के नाम पर फ़ाइल का नाम भी रख रहे हैं, जो तुम आमतौर पर नहीं कर सकते।

तो real अधिक mature Nextflow pipelines में, तुम अक्सर एक दिए गए sample से जुड़े सभी metadata के साथ एक meta object pass करोगे। फिर तुम उसके आधार पर dynamic फ़ाइल names बना सकते हो, जो बहुत अधिक practical है।

अगर तुम best practices के साथ यह करने में interested हो, तो _training.nextflow.io_ पर एक side quest है, जो विशेष रूप से metadata और meta maps के बारे में है, तो तुम अधिक detail के लिए वहां dig कर सकते हो।

## 3. एक array के माध्यम से कई inputs प्रदान करना

ठीक है। अब हम थोड़ा explore करने वाले हैं कि channels कैसे structured हैं और वे coding language में अन्य प्रकार के data structures से कैसे अलग हैं। और मैं थोड़ा सोचने वाला हूँ कि मैं potentially एक array का उपयोग कैसे कर सकता हूँ, जो एक familiar अवधारणा हो सकती है अगर तुम अन्य भाषाओं से आए हो।

क्या मैं एक channel में एक array का उपयोग कर सकता हूँ? चलो try करें। मैं एक array बनाने वाला हूँ, और मैंने इसे docs से copy किया है, _"greetings_array"_ और _"Hello", "Bonjour"_ और _"Holà"_। और फिर मैं इसे अपने hardcoded strings के बजाय यहां डालने वाला हूँ। तो मैं "channel.of" _"greetings_array"_ कहने वाला हूँ, इस array को एक channel में pass कर रहा हूँ। चलो try करें।

Terminal खोलो, और pipeline रन करो।

ठीक है। तुम देख सकते हो कि view statement ने यहां हमारे array को expected रूप से print किया, लेकिन फिर यह सारा red text, या यह red नहीं होगा अगर तुम्हारे पास अभी भी _"-ansi-log"_ off है, लेकिन यह सारा red text हमें बता रहा है कि कुछ गलत हो गया।

हमारे पास अब यहां एक अच्छा green tick नहीं है। हमारे पास एक red cross है, और अगर मैं इसे थोड़ा चौड़ा करता हूँ ताकि यह पढ़ना आसान हो, तो Nextflow हमें बता रहा है कि क्या गलत हुआ।

तो चलो इसे section दर section break down करें। यह कहता है कि error के कारण, और फिर error का reason, जो missing output फ़ाइलें हैं। तो मूल रूप से उस output block ने कहा कि यह फ़ाइल बनाई जानी चाहिए और वह नहीं थी। फिर यह कहता है कि यह command executed किया गया था। तो यह मूल रूप से उस _.command.sh_ फ़ाइल की contents है। यह इस तरह दिखती थी जब उन सभी variables को डाला गया था।

और तुम यहां देख सकते हो कि हमारा echo command वास्तव में केवल एक बार रन किया गया है और इसने पूरे array का उपयोग किया है, लेकिन एक string representation में, जो वास्तव में वह नहीं था जो हम चाहते थे।

और फिर command उस तरह exit हो गया, और वह work directory थी जहां हम जाकर फ़ाइलें देख सकते हैं और थोड़ा और समझ सकते हैं।

ठीक है। तो तब क्या हुआ था। Nextflow ने बस इस पूरे array को एक single channel element के रूप में process में pass किया, जिसका मतलब था कि process केवल एक बार रन हुआ। इसमें एक task था और इसने data को उस structure में उपयोग नहीं किया जिसकी हम उम्मीद करते थे।

## 3.2. channel contents को transform करने के लिए एक operator का उपयोग करना

तो हमें इस channel के साथ पहले कुछ करने की जरूरत है, इससे पहले कि इसका उपयोग किया जा सके। और यह operators का उपयोग करने के लिए एक stage set कर रहा है, जो विशेष functions हैं जिनका हम channels पर channel contents को manipulate करने के लिए उपयोग कर सकते हैं।

इस मामले में, हम कुछ उपयोग करने वाले हैं जिसे _flatten_ कहा जाता है। जिसे हम यहां channel के अंत में pass करते हैं। तो हम channel बनाते हैं और फिर हम _flatten_ रन करते हैं। और फिर से, अगर हम इस पर hover करते हैं, तो यह हमें इस command के लिए documentation सीधे VS Code में दिखाता है, जो बहुत helpful है। तुम Nextflow website पर, documentation पर भी ये सभी docs पा सकते हो।

मैं अभी इस code को रन कर सकता हूँ और देख सकता हूँ कि क्या यह काम करता है, लेकिन यह operators और Nextflow code के भीतर dynamic code कैसे करना है यह introduce करने का भी एक अच्छा अवसर है, जिन्हें closures कहा जाता है।

तो मैं _flatten_ रन करने से पहले यहां एक view command add back करने वाला हूँ। और यहां इस एक में ये squiggly brackets हैं, जो dynamic closure है। और इसके भीतर बस कुछ arbitrary code है जो execute किया जाएगा, एक view operator के context के भीतर।

यहां, यह कह रहा है greeting लो, जो view operator का input है, और वह यहां है। मैं इसे जो चाहूं कह सकता था, मैं इसे _"foo"_ कह सकता था और मुझे बस इसे बाद में _"foo"_ के रूप में refer करने की जरूरत है। और फिर मैं इसके साथ कहता हूँ, यह return करो।

और फिर एक string return कर रहा हूँ जो एक variable के लिए flatten से पहले कहता है। बहुत सरल।

मैं अब एक और exactly same add करने वाला हूँ, लेकिन मैं _flatten_ के बाद कहने वाला हूँ।

तो यह क्या करता है, क्योंकि यह sequence में रन होता है, तुम देखने वाले हो कि _flatten_ रन करने से पहले channel कैसा दिखता है, और फिर हम _flatten_ रन करने के बाद फिर से।

और फिर यह greeting channel अभी भी बनाया गया है, तो यह अभी भी process में pass होने वाला है। और उम्मीद है कि अब workflow रन होगा। चलो इसे try करें।

बढ़िया। तो सबसे पहले यह है कि pipeline इस बार crash नहीं हुआ। हमारे पास तीन processes थे जो properly रन हुए और हमें एक छोटा tick mark मिला है। और फिर हम देख सकते हैं कि हमारे view statements ने काम किया।

हमारे पास _flatten_ से पहले है, जो वह array है जो हमने पहले failure से देखा था, और फिर हमारे पास तीन बार _flatten_ के बाद कहा गया था जहां हमारे पास _"Hello", "Bonjour",_ और वे सभी अन्य तीन अलग-अलग elements array में हैं, जो अब जैसा हमने उम्मीद की थी, channel में तीन अलग-अलग elements हैं।

और तुम देख सकते हो कि _view_ operator तीन बार रन किया गया था। और ऐसा इसलिए है क्योंकि _flatten_ के बाद इस channel में अब तीन elements हैं। और इसलिए operator तीन बार called हो जाता है।

बहुत जल्दी, मैं सिर्फ mention करूंगा कि जब मैं पहले channel factories बना रहा था, मैंने _"."_ किया, और फिर हमने देखा कि channels बनाने के बहुत सारे अलग-अलग तरीके थे, और उनमें से एक को "_fromList"_ कहा जाता है। और वह वास्तव में विशेष रूप से इसी operation को करने के लिए designed है। तो हम बस from list greetings away कर सकते थे, और वह काम करेगा। यह थोड़ा clean और अच्छा syntax है। लेकिन इस demonstration के उद्देश्यों के लिए, हम इसे थोड़ा अधिक step-by-step बनाना चाहते थे ताकि तुम देख सको कि channel को कैसे manipulate किया जा रहा है और कैसे अलग-अलग operators एक channel के content में जो है उसे बदल सकते हैं।

## 4. एक CSV फ़ाइल से input values read करना

ठीक है, हम इसे थोड़ा और realistic कैसे बना सकते हैं? तुम शायद अपनी Nextflow pipeline में hard coded arrays के साथ बहुत सारा code बनाना नहीं चाहोगे। तुम शायद data को बाहर से लेना चाहोगे जब तुम launch करते हो, और वह data लगभग निश्चित रूप से फ़ाइलों में होने वाला है।

तो अगली चीज़ जो हम करने वाले हैं वह यह है कि हम इसे replicate करने वाले हैं, लेकिन एक single CLI parameter या एक hardcoded string या array से data लेने के बजाय, हम इसे एक फ़ाइल से लेने वाले हैं।

तो चलो अपने greetings away से छुटकारा पाएं। और अब हम इस channel factory को फिर से change करने वाले हैं। मैंने अभी कहा कि चुनने के लिए एक bunch थे और एक है _".fromPath"_। और मैं इसे बताने वाला हूँ कि, इस मामले में, _params.input_ लो, जो हमारे input पर वापस जा रहा है जिसका हम पहले उपयोग कर रहे थे।

अब वह parameter वास्तव में अभी तक उपयोग किए जाने के लिए तैयार नहीं है। हम अभी भी कह रहे हैं कि यह एक string है और यहां एक default के साथ hard coded है, लेकिन हम उस string को overwrite कर सकते हैं। हम अब चाहते हैं कि यह इसके बजाय एक फ़ाइल हो। तो type अलग है। यह अब एक _String_ नहीं है। यह एक _Path_ है।

और फिर हम default set कर सकते हैं अगर हम चाहें, फिर से एक Path पर। और अगर मैं बाईं ओर explore में देखता हूँ, तो तुम देख सकते हो इस repository में, इस working directory में, मेरे पास data नाम की एक directory है। मेरे पास वहां _"greetings.csv"_ नाम की एक फ़ाइल है।

तो मैं बस यहां default को _"data/greetings.csv"_ पर set कर सकता हूँ। अब, जब मैं इस pipeline को फिर से बिना किसी command line options के रन करता हूँ, तो यह इस default value का उपयोग करेगा। यह जानता है कि यह एक path है, इसलिए यह जानता है कि इसे एक path के रूप में handle करना चाहिए न कि एक string के रूप में।

और फिर यह इस _params.input_ से एक channel factory में pass होने वाला है और हमारा channel बनाएगा, जो फिर इस process में उपयोग किया जाने वाला है जिसे _sayHello_ कहा जाता है। चलो इसे try करें।

ठीक है। Failed। चिंता मत करो। यह expected था। और अगर तुम training material follow कर रहे हो, तो तुम देखोगे कि यह वहां भी expected था। चलो देखें कि क्या हो रहा है।

इसने pipeline रन करने की कोशिश की है। इसने process को execute करने की कोशिश की है, और इसे पहले हमने देखे गए एक के समान error मिली है।

यहां यह कहता है: हमने \_echo रन करने की कोशिश की, लेकिन इस CSV फ़ाइल की contents को echo करने के बजाय, इसने बस path को echo किया। और तुम देख सकते हो कि यह इस CSV फ़ाइल का पूरा absolute path यहां है।

और फिर निश्चित रूप से, क्योंकि इसने इस वास्तव में complicated path पर लिखने की कोशिश की, यह वास्तव में नहीं जानता था कि क्या करना है। और यह process work directory के scope के बाहर था।

मैंने शुरुआत में mention किया था कि Nextflow हर executed task को एक विशेष work directory के भीतर encapsulate करता है। और अगर तुम data लिखने की कोशिश करते हो, जो उस work directory के बाहर है, तो Nextflow तुम्हें एक safety precaution के रूप में रोकेगा। और यही यहां हुआ है। हमने एक absolute path पर लिखने की कोशिश की और Nextflow fail हो गया और हमें रोक दिया।

## 4.2. फ़ाइल को parse करने के लिए splitCsv() operator का उपयोग करना

ठीक है, चलो इस channel को देखें और देखें कि यह कैसा दिखता है। हम _".view"_ कर सकते हैं, और मैंने इसे website से copy किया है। तो _.view_, और हमारे पास यहां एक dynamic closure है और हम input के रूप में एक variable name "_csv"_ कहते हैं। तो वह channel contents है, और हम splitCsv से पहले कहते हैं, और यह ऐसा दिखता है।

अगर मैं इसे फिर से रन करता हूँ, तो यह अभी भी fail होगा, लेकिन यह हमें दिखाएगा कि इस channel के अंदर क्या है। यह विशेष रूप से exciting नहीं है। यह वह _path_ variable है। तो तुम देख सकते हो कि यह यहां सिर्फ एक string है क्योंकि इसे एक terminal पर print किया जा रहा है, लेकिन यह एक _path_ object है, जिसमें इस फ़ाइल के बारे में information और metadata है।

हम फ़ाइल के metadata को input में pass नहीं करना चाहते। हम उस फ़ाइल की contents pass करना चाहते हैं। अगर हम _greetings.csv_ फ़ाइल को देखें, तो तुम यहां देख सकते हो कि इसमें ये अलग-अलग variables हैं। _Hello, Bonjour, Holà_ फिर से। और ये वास्तव में वो चीजें हैं जो हम अपनी process में pass करना चाहते हैं, न कि सिर्फ फ़ाइल को खुद एक single object के रूप में।

तो हमें इस CSV फ़ाइल को parse करने की जरूरत है। हमें इसे unpack करने की जरूरत है, CSV फ़ाइल की contents पर पहुंचने की जरूरत है, और फिर channel के भीतर contents को process में pass करने की जरूरत है।

जैसा कि तुम शायद log message से बता सकते हो, हम _splitCsv_ का उपयोग करना चाहते हैं, जो एक और operator है, एक और channel operator। तो अगर मैं "_dot" "s"_ करता हूँ, और फिर तुम देख सकते हो कि यह auto suggested है। ओह, _splitCsv_ और कुछ brackets।

और फिर _splitCsv_ के बाद, मैं एक और _view_ statement डालने वाला हूँ बस ताकि हम देख सकें कि बाद में यह कैसा दिखता है। चलो pipeline रन करें और देखें कि हमें क्या मिला।

ठीक है। यह अभी भी fail हुआ, लेकिन एक नए और exciting तरीके से, जो progress है।

इस बार फिर से, हमें अपनी script के साथ कुछ समस्या है, जो render हुई है। अब। हमें final path नहीं मिला है, लेकिन हमें variables का एक array मिला है, जो बहुत कुछ उस error की तरह दिखता है जो हमें पहले मिला था जब हम एक array को एक fixed input के रूप में pass कर रहे थे।

View operator से हमारी logging के साथ, हम देख सकते हैं कि _splitCsv_ से पहले path था। और निश्चित रूप से, _splitCsv_ के बाद, हमारे पास तीन अलग-अलग outputs हैं और उन outputs में से प्रत्येक बहुत कुछ _greetings.csv_ फ़ाइल की प्रत्येक row की तरह दिखता है, जो समझ में आता है।

तो यहां क्या हुआ है कि Nextflow ने इस CSV फ़ाइल को parse किया, हमें तीन objects दिए, CSV फ़ाइल की प्रत्येक line के लिए एक array। तो फिर तीन बार हमने variables के एक array को channel में एक single string value के बजाय pass किया।

ठीक है, तो पिछली बार जब हमारी यह समस्या थी, तो हमने _flatten_ का उपयोग किया था। चलो बस बहुत जल्दी। Flatten try करें और देखें कि क्या होता है।

मैं इन variables को कुछ भी कह सकता हूँ। तो मैं इसे _myarray_ कहने वाला हूँ क्योंकि यह अब वास्तव में एक CSV नहीं है। चलो इसे फिर से रन करने की कोशिश करें और देखें कि _flatten_ के साथ क्या होता है।

तो इस बार हम रन करने वाले हैं, हमने CSV को तीन array objects में parse किया, और फिर हमने इसे flatten किया। और इस बार यह pass हुआ। और Nextflow pipeline रन हुआ। हालांकि तुम देख सकते हो कि _flatten_ वास्तव में town पर जाता है और सब कुछ flatten कर देता है। और इसलिए हमें प्रत्येक row के लिए तीन independent array entries मिलती हैं। और इसलिए इसने process को हर CSV की row पर तीन बार रन किया। और अब हमारे पास results फ़ाइलों का एक पूरा bunch है, और 123, 456, और सभी प्रकार की चीजें, न कि सिर्फ CSV के पहले column की, जो वास्तव में हम चाहते थे।

## 4.3. greetings extract करने के लिए map() operator का उपयोग करना

तो हम केवल पहले column पर कैसे पहुंचें? अगर flatten यहां बहुत simplistic है, तो हमें एक अधिक complex operator की जरूरत है जहां हम वास्तव में customize कर सकें और बता सकें कि हम CSV से क्या चाहते हैं।

ऐसा करने के लिए, हम _map_ का उपयोग करने जा रहे हैं। मूल रूप से _map_ बस कहता है, कुछ code, कुछ function रन करो हर element पर जो मुझे दिया जाता है और उस पर किसी प्रकार का transformation करो। और क्योंकि यह इतना flexible है, तुम इसे Nextflow code में हर समय आते हुए देखोगे।

अपने आप में, यह कुछ भी नहीं करता। तो हम regular brackets नहीं चाहते, हम यहां एक closure चाहते हैं और हमें इसे बताने की जरूरत है कि क्या करना है। तो मैं _"row"_ कहने वाला हूँ, क्योंकि यह CSV से rows दी जा रही है, तो यह एक logical variable name है। Input है। और मैं बस उस array के पहले element को return करना चाहता हूँ।

Nextflow में arrays zero based हैं, तो हम बस पहला element कहने वाले हैं, जो row zero है। अगर हम दूसरे column को चाहते हैं, तो मैं एक हो सकता हूँ या तीसरे column को दो, और इसी तरह। हम यहां जो चाहें return कर सकते हैं, लेकिन मैं बस पहली value return करने वाला हूँ।

और अब, हम pipeline को फिर से रन कर सकते हैं और देख सकते हैं कि क्या यह वह करता है जो हम उम्मीद करते हैं।

निश्चित रूप से, _splitCsv_ के बाद हमारे पास हमारे arrays हैं, और फिर _map_ के बाद, हमारे पास हमारी अच्छी clean strings हैं, बस _"Hello", "Bonjour"_ और _"Holà"_। और pipeline अब वह कर रहा है जो हम चाहते हैं। शानदार।

तो हम अब इन सभी view commands से छुटकारा पा सकते हैं। हमें अब उनकी जरूरत नहीं है।

## Recap

हमने अपनी तरह की debugging समाप्त कर दी है और यह वह code है जिसके साथ हम end करते हैं। अपने CLI parameter को _input_ कहा जाता है लेते हुए, जिसे _Path_ के रूप में classed किया गया है। Nextflow path ढूंढता है, इसे load करता है, और CSV फ़ाइल को समझता है। सभी अलग-अलग rows return करता है। और फिर हम channel में उस row के केवल पहले element को map करते हैं जो हमें channel contents देता है, जो process में pass किया जाता है।

और process channel में प्रत्येक element पर रन होता है, जो तीन है। और यह process को तीन बार रन करता है, इसे तीन tasks देता है। और फिर उन results को workflow से publish किया जाता है, process output द्वारा pick किया जाता है। Workflow से published किया जाता है और output block में _"hello_channels"_ नाम की एक subdirectory में save किया जाता है।

बहुत अच्छा। हम अब कुछ ऐसे पर पहुंच रहे हैं जो एक real life Nextflow pipeline के अधिक करीब है जिसे तुम कुछ real analysis के लिए रन कर सकते हो।

## सारांश

ठीक है। उम्मीद है कि अब तुम्हें Nextflow channels और operators क्या हैं और operators channels पर कैसे काम करते हैं और तुम उन्हें कैसे बना सकते हो, इसका एहसास हो रहा है।

Channels, जैसा मैंने इस video की शुरुआत में कहा, Nextflow का गोंद हैं। और तुम यहां देख सकते हो कि हम अलग-अलग inputs ले सकते हैं और उन्हें manipulate कर सकते हैं और उस data को ले सकते हैं और फिर उन्हें downstream workflow logic में pass कर सकते हैं।

और यह workflow block यहां वास्तव में वह जगह है जहां तुम वह सब parallelization और सभी clever logic build करते हो, और Nextflow को explain करते हो कि अपना workflow DAG कैसे build करना है, और अपनी pipeline को कैसे orchestrate करना है।

Channels समझने के लिए सबसे आसान अवधारणा नहीं हैं। तो एक break लो, इसके बारे में थोड़ा सोचो, शायद material को फिर से पढ़ो, और वास्तव में सुनिश्चित करो कि तुमने इन अवधारणाओं को समझ लिया है क्योंकि यह Nextflow की तुम्हारी समझ की कुंजी है और जितना बेहतर तुम channels और अलग-अलग channel operators और अलग-अलग channel factories को समझते हो। Nextflow लिखने में उतना ही मज़ा आएगा और तुम्हारी pipelines उतनी ही अधिक शक्तिशाली होंगी।

यह Python या अन्य भाषाओं में regular programming के समान नहीं है। हम यहां _if_ statements का उपयोग नहीं कर रहे हैं, यह channels और operators का उपयोग करते हुए functional flow programming है। तो यह थोड़ा अलग है, लेकिन यह भी super powerful है।

यह इस चैप्टर का अंत है। जाओ और एक quick break लो और मैं तुम्हें भाग तीन में अगले video में देखूंगा जहां हम Hello Workflow से गुजरने वाले हैं, और workflows के बारे में थोड़ा और बात करने वाले हैं।

बिल्कुल पिछले चैप्टर की तरह, यहां webpage के नीचे कुछ quiz questions हैं, तो तुम इनके माध्यम से एक quick run ले सकते हो और सुनिश्चित कर सकते हो कि तुम material के सभी अलग-अलग भागों को समझते हो जो हमने अभी किया है। और इसके अलावा, मैं तुम्हें अगले video में देखूंगा। बहुत-बहुत धन्यवाद।

ठीक है।
