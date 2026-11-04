# भाग 6: Hello Config - वीडियो ट्रांसक्रिप्ट

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता से अनुवादित - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "महत्वपूर्ण नोट"

    यह पेज केवल ट्रांसक्रिप्ट दिखाता है। पूर्ण चरण-दर-चरण निर्देशों के लिए, [कोर्स मटेरियल](../06_hello_config.md) पर वापस जाओ।

    ट्रांसक्रिप्ट में दिखाए गए सेक्शन नंबर केवल संकेत के लिए प्रदान किए गए हैं और हो सकता है कि मटेरियल में सभी सेक्शन नंबर शामिल न हों।

## स्वागत

नमस्ते, और Hello Nextflow के भाग छह में वापस स्वागत है। यह सेक्शन पूरी तरह से configs के बारे में है, और यह इस कोर्स का आखिरी भाग है।

Nextflow विशेष रूप से दो चीजों में अच्छा है, reproducibility और portability। Configs वह जगह है जहां हम इनमें से दूसरी चीज को वास्तव में चमकते हुए देखते हैं। एक Nextflow पाइपलाइन को विभिन्न तरीकों से चलाने और विभिन्न सिस्टमों पर काम करने के लिए configure करने की क्षमता, बिना अंतर्निहित पाइपलाइन कोड को edit किए।

यह superpower Nextflow पाइपलाइनों को अन्य लोगों द्वारा विभिन्न स्थानों पर reuse करने की अनुमति देती है, या विभिन्न infrastructures पर जिन तक तुम्हारी खुद की पहुंच हो सकती है।

इसका मतलब है कि तुम अपने लैपटॉप पर पाइपलाइन कोड develop कर सकते हो, इसे cloud पर push कर सकते हो, अपने HPC पर चला सकते हो, और यह वही पाइपलाइन कोड है और यह हर जगह चलता है।

इस सेक्शन में, हम कुछ topics से गुजरने वाले हैं। हम यह शुरू करेंगे कि Nextflow config फ़ाइलों को कैसे handle करता है, यह उन्हें कहां से लोड करता है, और तुम उन्हें कैसे लिखते हो और कैसे structure करते हो, और पाइपलाइन और config फ़ाइल में क्या जाना चाहिए के बीच वह separation।

फिर हम कुछ सामान्य use cases पर आएंगे जैसे कि आउटपुट फ़ाइलें कहां store की जाती हैं इसे बदलना, और साथ ही पाइपलाइन को विभिन्न infrastructures पर काम करने के लिए कैसे प्राप्त करें, दोनों सॉफ़्टवेयर packaging के विभिन्न प्रकार का उपयोग करके या विभिन्न infrastructures पर jobs submit करके।

## Config फ़ाइल hierarchies

ठीक है, चलो शुरू करते हैं। जब config फ़ाइलों को लोड करने की बात आती है, तो Nextflow कई अलग-अलग जगहों से pull कर सकता है, जो एक अच्छी बात है और थोड़ा जोखिम भरा भी हो सकता है क्योंकि कभी-कभी यह थोड़ा मुश्किल हो सकता है यह जानना कि यह कहां से config फ़ाइल प्राप्त कर रहा है और किस क्रम में चीजें लोड करता है।

तो मैं वास्तव में recommend करता हूं कि तुम इस लिंक पर क्लिक करो, जो हमें Nextflow docs में ले जाता है। और इस configuration पेज पर, यह मुख्य स्थानों को सूचीबद्ध करता है जहां से config लोड किया जाता है, और महत्वपूर्ण रूप से, वह क्रम जिसमें ये चीजें लोड की जाती हैं।

तो तुम देख सकते हो, तुम अपनी Nextflow home directory में एक config फ़ाइल डाल सकते हो, जो आमतौर पर तुम्हारे home dir में ".nextflow" होती है। और वह फ़ाइल हमेशा तुम्हारे सिस्टम पर हर Nextflow run द्वारा लोड की जाएगी।

देखने के लिए अगली जगह तुम्हारी पाइपलाइन, repository या directory की root में "nextflow.config" नामक फ़ाइल है।

उसके बाद, एक और फ़ाइल जिसे "nextflow.config" कहा जाता है, लेकिन इस बार उस directory में जहां से तुम Nextflow launch कर रहे हो: launch directory।

अंत में, तुम command line पर "-c" argument के साथ config फ़ाइल paths प्रदान कर सकते हो, और तुम इसे कई बार कर सकते हो। और वे उस क्रम में applied किए जाते हैं जो तुम specify करते हो।

तुम चाहो तो इन सभी locations में config फ़ाइलें प्रदान कर सकते हो, और वे iteratively लोड की जाएंगी, प्रत्येक केवल उन config scopes में पिछले को overwrite करते हुए जहां वे clash करते हैं।

यह एक वास्तव में powerful सिस्टम है क्योंकि इसका मतलब है कि तुम sensible defaults सेट कर सकते हो और फिर धीरे-धीरे अधिक से अधिक specific हो सकते हो जैसे-जैसे तुम उस config पर narrow down करते हो।

## 0. Warmup: hello-config.nf चलाओ

ठीक है, चलो इसे बंद करते हैं और अपने Codespaces में jump करते हैं और शुरू करते हैं। पहले की तरह मैंने यहां साफ़ कर दिया है, मैंने अपनी पिछली results directories, मेरे Nextflow, और मेरी work directories वगैरह को हटा दिया है। चिंता न करो अगर तुम्हारे पास अभी भी वे फ़ाइलें पड़ी हैं। यह सिर्फ इसलिए है क्योंकि मैं बहुत zoomed in हूं और इसलिए चीजें बहुत जल्दी messy हो जाती हैं अन्यथा।

हम hello-config.nf के साथ काम करने वाले हैं, हमारी directory में आखिरी फ़ाइल, और यह पिछले सेक्शन में जहां हमने छोड़ा था वहां से follow करना चाहिए।

तो हमारे पास हमारे चार अलग-अलग processes हैं, जो module फ़ाइलों से include किए जा रहे हैं। हमारे पास हमारे पाइपलाइन parameters हैं, हमारा workflow block जहां हम विभिन्न processes को call कर रहे हैं और channels को एक साथ stitch कर रहे हैं, आउटपुट channels को publish कर रहे हैं, और फिर नीचे output block जहां हम define करते हैं कि वे फ़ाइलें कहां store की जानी चाहिए और उन्हें कैसे copy किया जाना चाहिए।

हमारे पास पिछले chapter से पहले से ही एक "nextflow.config" फ़ाइल भी है, जहां हम Docker enable करते हैं, और हम आज इस फ़ाइल में build करने वाले हैं।

पहले की तरह, हमने इस main script में आउटपुट path को hello config में बदल दिया है, बस इसलिए कि यह पिछले results से clash न करे जो तुमने generate किए हैं।

ठीक है, चलो जल्दी से check करते हैं कि सब कुछ अभी भी हमारी अपेक्षा के अनुसार काम कर रहा है। एक terminal ऊपर लाओ और हम nextflow run hello-config.nf करते हैं। Nextflow लोड होता है। हमारे चार अलग-अलग processes को चलाना चाहिए। cowpy का उपयोग करके कुछ अच्छे ascii artwork generate करें और फिर उस directory में हमारे results को हमारी results फ़ाइलों में save करें।

मैं यहां जल्दी से देख सकता हूं बस यह सुनिश्चित करने के लिए कि ये फ़ाइलें जैसी हम उम्मीद करते हैं वैसी दिखती हैं, और निश्चित रूप से, वहां हमारा विशाल Turkey है। बढ़िया।

## 1.1. डिफ़ॉल्ट values को nextflow.config में मूव करें

अब पहली चीज जो हम करने जा रहे हैं वह है कुछ चीजों को हमारी script से हमारी config फ़ाइल में move करना शुरू करना।

और जिस चीज की हमें परवाह है वह ज्यादातर इस stage पर parameters हैं। हम डिफ़ॉल्ट values को config फ़ाइल में लेना चाहते हैं, ताकि यह स्पष्ट हो कि defaults क्या हैं और लोगों के लिए उन्हें overwrite करना आसान हो।

मैं इस params block को यहां से script से लेने जा रहा हूं और इसे config फ़ाइल में डालने जा रहा हूं। और हमें यहां थोड़ा सावधान रहने की जरूरत है, क्योंकि अभी config और scripts के बीच syntax थोड़ा अलग है। Config फ़ाइल type declarations नहीं ले सकती क्योंकि हम वास्तव में इन params को define नहीं कर रहे हैं, हम बस उन्हें reference कर रहे हैं। तो मैं उनसे छुटकारा पाने वाला हूं।

लेकिन अन्यथा यह बहुत समान है। हमारे पास एक params block है और फिर हमारे पास हमारे अलग-अलग input parameters हैं, batch parameter, character parameter।

मैं अब अपनी script पर वापस जा सकता हूं और मुझे अब इन defaults को define करने की जरूरत नहीं है क्योंकि ये values अब मेरी Nextflow config फ़ाइल में हैं।

हालांकि, मैं parameter names और उनके types को छोड़ देता हूं, ताकि Nextflow को वह जानकारी पता हो और अभी भी सभी type safety और सब कुछ कर सके।

ठीक है। हम उन फ़ाइलों को save करते हैं और जल्दी से check करते हैं कि सब कुछ अभी भी पहले जैसा काम करता है। यहां कोई बदलाव नहीं होना चाहिए। हमने values को same रखा है। हमने बस बदल दिया है कि वे कहां define किए गए हैं।

बढ़िया।

## 1.2. एक run-specific configuration फ़ाइल का उपयोग करें

अब, अब तक हम उसी directory से Nextflow launch कर रहे हैं जहां हमारी पाइपलाइन script है। तो हमारी launch dir और हमारी pipeline dir एक ही चीज हैं।

यह दिखाने के लिए कि हम विभिन्न launch directories के साथ विभिन्न config फ़ाइलें कैसे रख सकते हैं, हम अब एक नई subdirectory बनाने जा रहे हैं।

तो मैं mkdir कहने जा रहा हूं, और हम इसे tux-run कहने वाले हैं।

और फिर मैं cd करने जा रहा हूं, tux-run में directory बदलने जा रहा हूं। और ध्यान दो कि हम अब अपनी working directory में अब पाइपलाइन scripts वाली directory में नहीं हैं।

ठीक है, चलो एक नई "nextflow.config" फ़ाइल बनाते हैं। तो touch Nextflow config, और चलो इसे बस VS Code में खोलते हैं। तुम sidebar में यहां भी देख सकते हो कि हम अब इस subdirectory में हैं।

अब हम वही params block ले सकते हैं जो हमारे पास top level nextflow.config में था, इसे copy करें और अब हम इन values को बदल सकते हैं।

सबसे पहले, data अब एक अलग relative path है क्योंकि हम एक subdirectory में हैं, तो हमें इसे update करना होगा। और फिर हम batch को experiment में बदलने वाले हैं, और हम character को Turkey से tux में बदलने वाले हैं।

अब वहां save पर क्लिक करें, और चलो इसे try करते हैं। बस data की तरह, मुझे अब script तक पहुंचने के लिए ../ कहना होगा। तो यह Hello config है। और मैं enter दबाता हूं।

पाइपलाइन कोड बिल्कुल भी नहीं बदला है, लेकिन अब हमारे पास दो sets of config लोड होने वाले हैं, और launch dir config फ़ाइल defaults को overwrite करना चाहिए, जो pipeline nextflow.config में set किए गए थे, और हमें results के विभिन्न sets मिलने चाहिए।

निश्चित रूप से, हमारी directory के भीतर यहां, tux-run के भीतर, तुम देख सकते हो कि हमारे पास एक dot Nextflow directory और एक work directory है और ऐसा इसलिए है क्योंकि ये हमेशा तुम्हारी launch directory में बनाई जाती हैं। तो ये उन work और results वालों से अलग हैं जो हमारे पास पहले के runs से थे।

अब, अगर मैं results में देखता हूं, तो हम अपने collected को देख सकते हैं और वहां हमारा छोटा tux character है। तो तुम देख सकते हो कि वे parameters ठीक से interpret किए गए थे।

## 1.3. एक parameter फ़ाइल का उपयोग करें

ठीक है। पहले जब मैं विभिन्न config फ़ाइलों के बारे में बात कर रहा था जो लोड की जा सकती हैं, तो मैंने एक और जगह छोड़ दी थी जहां से हम config प्राप्त कर सकते हैं।

तुम इसे command line से प्राप्त कर सकते हो जैसा कि हमने dash dash parameter names के साथ देखा है, लेकिन हम एक YAML या एक JSON फ़ाइल भी supply कर सकते हैं, बस params की।

Config फ़ाइल में सभी विभिन्न प्रकार के scopes हो सकते हैं, लेकिन ये फ़ाइलें बस parameters हैं, और यह एक बार में कई parameters supply करने का एक अच्छा user-friendly तरीका है, और शायद थोड़ा अधिक reproducible तरीका क्योंकि तुम उन्हें फ़ाइल में लिखते हो, इसलिए बाद में उन्हें प्राप्त करना आसान है।

तो चलो अपने terminal पर वापस जाते हैं और इससे पहले कि हम भूल जाएं, सुनिश्चित करें कि हम एक directory ऊपर वापस जाएं, तो मैं अब subdirectory में नहीं हूं, और मैं उस YAML फ़ाइल को देखने जा रहा हूं जो हमारे पास यहां test-params.yaml नामक है।

तो अगर मैं बस code test-params.yaml करता हूं, तो तुम देख सकते हो कि यह बस एक regular YAML फ़ाइल है। इसके बारे में कुछ खास नहीं है। keys हमारे parameter names होने के साथ, YAML formatting के साथ तो यहां एक colon, और फिर एक value।

ध्यान दो कि यह Nextflow code नहीं है, इसलिए हम यहां variables जैसी चीजें नहीं डाल सकते। ये बस static values हैं।

इसके अलावा क्योंकि JSON वास्तव में YAML के रूप में parse करता है, हमारे पास एक test-params.json फ़ाइल भी हो सकती है, जो बहुत समान दिखती है। यह बस अलग data format है।

तो हमारे पास यहां दो अलग-अलग test फ़ाइलें हैं और हमारे पास थोड़े अलग variables हैं।

ठीक है, तो हम इन्हें Nextflow को कैसे देते हैं? यह बहुत आसान है। हम Nextflow run hello config करते हैं, जैसा पहले था। और config फ़ाइल के लिए "-c" के बजाय, या उन डिफ़ॉल्ट फ़ाइल names को लोड करने के बजाय, हम -params-file करते हैं। Single hyphen क्योंकि यह एक core Nextflow option है।

और फिर उस फ़ाइल के लिए path pass करें। तो मैं "-params-file test-params.yaml" करने वाला हूं, और हम देखेंगे कि वे ठीक से लोड हुए हैं या नहीं।

ठीक है। यह चला। चलो बस खुद को याद दिलाते हैं कि इस YAML फ़ाइल में क्या था। तो batch को YAML पर set किया गया था, तो इसे यही कहा जाना चाहिए, और इसमें एक stegosaurus होना चाहिए। तो चलो ऊपर जाते हैं और results में देखते हैं। और हमारे पास COLLECTED-yaml है। तो देखते हैं कि क्या हमारे पास एक Stegosaurus है। शानदार, एक टोपी पहने हुए Stegosaurus। यही हम पसंद करते हैं।

तो यह वास्तव में अच्छी तरह से काम कर गया, और JSON फ़ाइल के साथ बिल्कुल वैसा ही है। हम बस यहां फ़ाइल extension को switch out करते हैं और Nextflow जानता है कि इसे कैसे read करना है।

और इस मामले में, हमारे पास JSON नामक एक batch होना चाहिए और हमारे पास एक turtle होना चाहिए। चलो एक नजर डालते हैं। अद्भुत। मेरे पसंदीदा CLI tools में से एक।

## 2.1. -output-dir के साथ आउटपुट directory को customize करें

ठीक है, तो यह ज्यादातर पाइपलाइन के inputs के बारे में और parameters बदलने के बारे में सोच रहा है। outputs के बारे में क्या?

अब, हालांकि हम params का उपयोग करके sub directories को बदल रहे हैं, तुमने देखा होगा कि हमारी सभी फ़ाइलें अभी भी results में जा रही हैं।

हम उस base directory को बदल सकते हैं जिसमें सभी फ़ाइलें -output-dir नामक एक command line flag के साथ publish की जाती हैं। तो अगर मैं Nextflow run hello config करता हूं, और फिर मैं -output-dir करता हूं, और हम इसे "custom-outdir-cli" कहने वाले हैं। टाइप नहीं कर सकता। बस इसलिए कि हम याद रखें कि ये फ़ाइलें कहां से आईं।

यह एक core Nextflow option है और यह एक बहुत नया है। यह हाल ही में जोड़ा गया था, और यह उन चीजों में से एक है जो हम new language parser और सब कुछ के साथ कर सकते हैं।

टाइप करने के लिए यह थोड़ा लंबा है। तुम चाहो तो इसे बस "-o" भी call कर सकते हो। तो अगर मैं बस वापस जाता हूं। बस इसे "-o" में shorten कर सकते हैं, जो थोड़ा आसान है।

ठीक है। हम वह चलाते हैं। हमने इस बिंदु पर अपनी पाइपलाइन में या यहां तक कि अपने config में कुछ भी नहीं बदला है, और यह उम्मीद से हमारे सभी results को एक अलग top level directory में save करना चाहिए। और तुम कल्पना कर सकते हो कि तुम इसे मूल रूप से किसी भी path पर set कर सकते हो जो तुम चाहते हो।

यह अभी top पर पहुंच गया है। हमारे पास एक custom-outdir-cli है, और सभी फ़ाइलें वहां बिल्कुल उसी तरह organize की गई हैं, उनकी same sub directories और फ़ाइल names के साथ। तो यह सिर्फ यह बदलने का एक वास्तव में आसान तरीका है कि पाइपलाइन अपने results को कहां publish करती है, इस बारे में बहुत अधिक सोचे बिना कि वे results कैसे organize किए गए हैं।

## 2.1.2. output block से hardcoded paths को हटाएं

अगर मैं इस directory में देखता हूं, तो हम देख सकते हैं कि हमारे पास अभी भी Hello Config नामक एक subdirectory है, जो अब थोड़ी redundant लगती है।

तो चलो बस अपनी script को फिर से लोड करते हैं और हम अब उस subdirectory को नीचे output block से हटा सकते हैं। क्योंकि हमें अब वास्तव में इसकी जरूरत नहीं है। तो हम अभी बस यह कर सकते हैं, इसे यहां से delete कर सकते हैं। और फिर अगर यह बस यही है, तो तुम या तो इसे पूरी तरह से delete कर सकते हो या इसे एक empty string के रूप में छोड़ सकते हो। मैं इसे अभी के लिए एक empty string के रूप में छोड़ने वाला हूं, क्योंकि हम वापस आने वाले हैं और भविष्य में इसकी जगह कुछ अलग चीजें डालने वाले हैं। लेकिन अगर तुम्हें sub directories की परवाह नहीं है, तो path declaration को पूरी तरह से हटाना cleanest है।

ठीक है, चलो save करें। बस जल्दी से इसे फिर से try करें। मैं वास्तव में अपनी "custom-outdir-cli" directory को हटाने वाला हूं ताकि हम वहां मौजूद किसी भी existing फ़ाइलों से confused न हों। क्योंकि याद रखें, जब तुम चीजें publish करते हो, तो यह उन फ़ाइलों को नहीं हटाता है जो पहले से वहां थीं। यह बस नई add करता है। चलो उस command को फिर से चलाते हैं, custom-outdir-cli।

और अब अगर तुम "ls custom-outdir-cli" करते हो, तो अब वहां Hello Config नामक कोई और directory नहीं है।

## 2.2.1. configuration फ़ाइल में outputDir सेट करें

ठीक है, यहां command line flag, "-o" या "-output-dir" अच्छा है। लेकिन config में इसके लिए defaults सेट करने के बारे में क्या? हम यह कैसे करते हैं?

मैं "nextflow.config" फ़ाइल खोलता हूं, बाकी सब कुछ बंद करता हूं और उससे छुटकारा पाता हूं। हम यहां एक नया config option add कर सकते हैं, जिसे मैंने बस training material website से copy किया है, और इसे outputDir कहा जाता है।

यह किसी scope के अंतर्गत नहीं है। यह params या कुछ भी के अंतर्गत नहीं है। यह top level है, और हम इसे एक string पर set कर सकते हैं। अब एक simple चीज बस इसे एक hard-coded string के रूप में results के अलावा कुछ और में बदलना है। लेकिन क्योंकि यह एक Nextflow config फ़ाइल में है, हम यहां थोड़ा clever हो सकते हैं और variables भी शामिल कर सकते हैं।

और तुम यहां देख सकते हो कि हमने एक params variable, params.batch शामिल किया है, जो इस string का हिस्सा है। इसका मतलब है कि हम variables को reuse कर सकते हैं जो अन्य स्थानों से आ रहे हैं। और इस मामले में, अगर हम --batch करते हैं, जब हम Nextflow Pipeline चलाते हैं, तो हमें अपने custom path में एक subdirectory मिलने वाली है जो batch name पर आधारित है।

ठीक है, तो चलो इसे try करते हैं और बस एक quick look डालते हैं कि results कैसे दिखते हैं। तो अगर मैं Nextflow run hello config और --batch my_run करता हूं। खुद को याद दिलाएं कि config कैसा दिखता था। तो यह custom-outdir-config है।

Tree custom-outdir-config। और तुम देख सकते हो कि batch को my_run कहा गया था। और फिर हमारे पास my_run नामक subdirectory है। तो वह dynamic file path काम कर गया।

और न केवल वह, यह अब एक डिफ़ॉल्ट results directory में नहीं गया, और मुझे base directory को बदलने के लिए command line पर कुछ भी specify करने की जरूरत नहीं थी। तो हमने डिफ़ॉल्ट outputDir के लिए डिफ़ॉल्ट value को सफलतापूर्वक reset कर दिया है।

## 2.2.2. batch और process names के साथ Subdirectories

ठीक है, चलो उसे थोड़ा आगे ले जाते हैं। यह config फ़ाइल के भीतर एक dynamic variable है। script के बारे में क्या? अब, अब तक हमारे पास यहां ये paths थे और ये भी dynamic हो सकते हैं। तो सिर्फ कुछ hard code करने के बजाय, हम कुछ squiggly brackets डाल सकते हैं और कुछ dynamic डाल सकते हैं।

तो उदाहरण के लिए, हमारे पास sayHello नामक हमारे processes हैं। हम sayHello.name कर सकते हैं, जो process का एक attribute है, जो थोड़ा boring है क्योंकि यह इस मामले में बस "sayHello" है। लेकिन यह variable है।

तो यह तुम्हें एक idea देता है। तो हम इसे यहां डाल सकते हैं और convertToUpper.name, collectGreetings.name, collectGreetings.name फिर से, और cowpy कह सकते हैं।

अब जब हम चलाएंगे, तो base directory अभी भी custom-outdir-config होने वाली है। और यह params.batch नामक एक subdirectory में होने वाली है, लेकिन उसके नीचे की sub directories को process name द्वारा organize किया जाना चाहिए।

चलो बस उसे try करते हैं और देखते हैं कि क्या यह काम करता है। तो मैं पिछली directory को हटाने वाला हूं ताकि हम confused न हों, और बिल्कुल same Nextflow Run command का उपयोग करें।

यह same तरीके से चलना चाहिए। मैं इन सभी पर dash resume का उपयोग कर सकता था ताकि इसे थोड़ा तेज बनाया जा सके और पहले calculated results का उपयोग किया जा सके। अब, अगर मैं tree custom-outdir-config करता हूं, तो तुम देख सकते हो कि यह results में नहीं है, यह हमारी base directory में है, batch name के साथ। और तुम देख सकते हो कि सभी results अब process के नाम पर named sub directories के भीतर organize किए गए हैं। तो हमारे पास दो अलग-अलग जगहें हैं जहां हम यहां dynamic output paths define कर रहे हैं।

ठीक है। अंतिम चीज, चलो उन intermediate folders को वापस जोड़ते हैं, जो हमारे पास पहले थे क्योंकि वे kind of अच्छे थे। Intermediates।

और हम इस params.batch के बारे में भी थोड़ा सोच सकते हैं, शायद एक pipeline developer के रूप में मुझे वास्तव में subdirectory में वह होना पसंद था, लेकिन अगर पाइपलाइन के end users CLI पर "-o" या -output-dir set कर रहे हैं, तो यह इस पूरे statement को पूरी तरह से overwrite कर रहा है, और हम उस subdirectory को खो देते हैं।

तो हम जो कर सकते हैं वह है कि हम उस dynamic path को outputDir config से बाहर ले सकते हैं, जिसे clobber किया जाएगा, और इसे output path में डाल सकते हैं, जो clobber नहीं किया गया है।

तो हम params.batch slash intermediates slash sayHello.name कर सकते हैं, और यह सब एक double quoted string में कर सकते हैं, ताकि यह Nextflow द्वारा interpolate किया जाए।

अब copy कर सकते हैं, whoops। इन्हें अन्य processes में copy करें। याद रखें कि इन सभी को quotes में डालें। और इन particular outputs से intermediates को हटाएं।

ठीक है? यह अब थोड़ा अधिक complex दिख रहा है, लेकिन तुम देख सकते हो कि हम वास्तव में अपने code में अच्छी organized output directory structure बनाना शुरू कर रहे हैं।

और जो वास्तव में अच्छा है वह यह है कि code में यह extra complexity CLI तक pass through नहीं होती है। तो हम अपनी command को -output-dir के साथ चला सकते हैं और जो भी batch variables हैं, बस पाइपलाइन को कैसे चलाना है इसके बारे में सोचते हुए और code में क्या है इसके बारे में वास्तव में बहुत अधिक नहीं सोचते हुए। और हमारी output फ़ाइलें एक बहुत अच्छी organized तरीके से बनाई जाने वाली हैं, जो मूल रूप से पाइपलाइन का उपयोग करने वाले लोगों के लिए अच्छा है।

बढ़िया। जैसे मैं यह लिख रहा हूं, मुझे एहसास होता है कि मैंने एक गलती की। देखें कि क्या किसी ने मुझे यहां पकड़ा। हमारे पास collectGreetings.name है, तो कुछ गलत हो गया है। और हां, निश्चित रूप से, मैं गलती से इन्हें squiggly brackets में डालना भूल गया।

तो याद रखें, जब तुम अपना code लिख रहे हो तो सावधान रहो और सुनिश्चित करो कि तुम Nextflow को बताओ कि क्या एक variable है और क्या बस एक string है। क्योंकि यह बिल्कुल वही करेगा जो तुम इसे करने के लिए कहते हो। और कुछ नहीं। सभी अच्छे computers की तरह। ठीक है, इससे यह ठीक हो जाना चाहिए।

## 2.3. workflow level पर publish mode सेट करें

इस script का एक bit है, जो मुझे अभी भी पसंद नहीं है, वह यह तथ्य है कि हम mode copy को बार-बार लिख रहे हैं, और अगर एक चीज है जो हमें पसंद नहीं है, तो वह है खुद को repeat करना।

तो हम इसे लेकर और config में move करके इसे थोड़ा साफ कर सकते हैं। और वास्तव में, हम इसे एक बार में पूरी पाइपलाइन के लिए सेट कर सकते हैं। तो हमें इसे कई बार कहने की जरूरत नहीं है।

हम अपनी config फ़ाइल पर जाते हैं और हमारे पास यहां workflow नामक एक नया scope है। और हम या तो squiggly brackets कर सकते हैं या हम dot notation कर सकते हैं। कोई फर्क नहीं पड़ता। Training material website dot notation का उपयोग करती है। मैं output कह सकता हूं और हम mix and match कर सकते हैं, तो mode equals copy। बढ़िया।

और अब हम यहां वापस जा सकते हैं और इन्हें delete कर सकते हैं। अब हम उन्हें place में छोड़ सकते थे। Config मूल रूप से यहां जो लिखा है उसे overwrite कर रहा है, लेकिन जैसा कि हमारे पास pipeline level config में है, और ये दोनों फ़ाइलें एक साथ ship होती हैं, वास्तव में इसे दो बार करने का कोई कारण नहीं है।

ठीक है। बस खुद को sanity check करें, क्योंकि जाहिर है हम गलतियां करते हैं। चलो, इसे फिर से चलाएं और बस check करें कि हम फ़ाइलों को publish करने के लिए copy mode का सही तरीके से उपयोग कर रहे हैं। तो हम script को फिर से चलाने वाले हैं और इस बार हमने results को config-output-mode नामक directory में डाल दिया है, देखें कि वहां फ़ाइलें कैसी दिखती हैं।

और फिर अगर मैं batch को देखने के लिए "ls -l" करता हूं, और हम उदाहरण के लिए cowpy को देख सकते हैं। और हमें देखना चाहिए, हां, कि यह यहां एक proper फ़ाइल है, जो एक soft link नहीं है, तो वह config attribute ठीक से applied किया गया है।

## 3. एक software packaging technology चुनें

ठीक है। अब तक हम inputs और outputs, उन फ़ाइलों पर focus कर रहे हैं जिनके साथ workflow चल रहा है। लेकिन infrastructure के बारे में क्या? मैंने शुरुआत में कहा था कि Nextflow तुम्हें विभिन्न computing setups पर same पाइपलाइन चलाने की अनुमति देता है। तो यह कैसा दिखता है?

यह दिखाने के लिए, हम cowpy चलाने के लिए Docker का उपयोग करने से switch करने जा रहे हैं, और इसके बजाय हम same काम करने के लिए Conda का उपयोग करेंगे।

मैं यह बहुत simply कर सकता हूं। अगर मैं code, "nextflow.config" पर जाता हूं। अगर तुम्हें याद है कि top पर, हमने पहले docker.enabled define किया था, और last chapter में ताकि हम cowpy के साथ container का उपयोग कर सकें।

मैं Nextflow को Docker का उपयोग न करने के लिए कहने वाला हूं। उसे false पर set करें। और मैं Conda enabled equals true कहने वाला हूं। तो Nextflow को बताएं, कृपया Conda का उपयोग करें।

अब सिर्फ Conda को enable करना अपने आप में पर्याप्त नहीं है। बस जैसा हमने Docker के साथ किया, हमें Nextflow को बताना होगा कि उसे जो सॉफ़्टवेयर चाहिए वह कहां से मिल सकता है।

तो अगर हम यहां modules में hop करें। और cowpy script खोलें। हम देख सकते हैं कि हमारे पास top पर एक container declaration है। और container को Docker द्वारा उपयोग किया जाता है, लेकिन Singularity, Apptainer, और कई अन्य software tools द्वारा भी।

लेकिन इसे Conda के लिए उपयोग नहीं किया जा सकता, इसलिए हमारे पास "conda" नामक एक अलग declaration है, और हम बस "cowpy" लिख सकते थे। और यह conda package resolution पर छोड़ देगा कि तुम्हारे local conda environment के अनुसार इसे solve करने का best तरीका figure out करे।

या यह अच्छा अभ्यास है कि training material website जो कहता है वह करें, जो एक specific conda channel को उसके double colon notation के साथ define करना है, और निश्चित रूप से software के एक specific version को define करें ताकि पाइपलाइन चलाने वाला हर व्यक्ति same version प्राप्त करे।

ध्यान दें कि containers इस संबंध में थोड़े superior हैं, क्योंकि जब तुम Conda के साथ कुछ install करते हो, तो यह अभी भी उस package के लिए सभी dependencies को figure out करने वाला है, और वे समय के साथ बदल सकते हैं। Dependency drift कहा जाता है।

तो containers, हालांकि, पूरी तरह से नीचे तक software dependencies के पूरे stack को lock करते हैं, इसलिए तुम थोड़ा अधिक confident हो सकते हो कि A, यह काम करने वाला है, और B, यह reproducible होगा।

तो अगर तुम Docker या Singularity या Apptainer का उपयोग करने में सक्षम हो, तो मैं निश्चित रूप से recommend करूंगा।

अब जो अच्छा है वह यह है कि module फ़ाइल, जो pipeline developer द्वारा लिखी गई है, अब दोनों Container और Conda है, और इसलिए हम इस पाइपलाइन को चलाने वाले व्यक्ति को बता रहे हैं, हमें कोई फर्क नहीं पड़ता कि तुम कौन सा software packaging solution उपयोग करते हो। यह Docker और Conda दोनों के साथ काम करेगा, और दोनों मामलों में software कहां से प्राप्त करना है यह है।

हम terminal को pull up कर सकते हैं और चलो इसे try करते हैं। तो Nextflow run hello config --batch conda। और पहली बार जब यह conda के साथ runs करता है, तो यह थोड़ा slow होने वाला है जब यह उस particular process पर पहुंचता है, क्योंकि इसे "conda install" चलाना होता है।

और यह सिर्फ इस एक process के लिए एक special conda environment बना रहा है। तो यह मेरे global conda environment का उपयोग नहीं कर रहा है, जो मेरे पास मेरे terminal पर है। यह सिर्फ उस एक process के लिए एक बना रहा है। यह अच्छा है क्योंकि यह तुम्हारे workflow में विभिन्न processes के बीच dependency clashes जैसी चीजों से बचता है। अगर तुम्हारे processes में tools हैं जिन्हें Python या ऐसी चीजों के विभिन्न versions की जरूरत है, तो यह ठीक है क्योंकि वे विभिन्न conda environments का उपयोग कर रहे हैं।

Nextflow इन conda environments को locally cache करता है, तुम देख सकते हो कि यह तुम्हें बताता है कि वह path कहां है, यह यहां work directory में है। और इसलिए अगली बार जब मैं इस script को Conda के साथ चलाता हूं, तो यह बहुत तेज होगा क्योंकि यह उस existing conda environment को find करेगा और बस इसे reuse करेगा। लेकिन पहली बार जब हम इसे करते हैं, तो इसे जाना होता है और इसे fetch करना होता है, resolve करना होता है, सभी dependencies को download करना होता है, और सब कुछ set up करना होता है।

ठीक है, बढ़िया, यह चला। हम बस खुद को याद दिला सकते हैं कि पाइपलाइन वर्तमान में क्या उपयोग करने के लिए configured है। अगर हम config फ़ाइल में देखें, तो यह अभी मेरे लिए "custom-outdir-config" था। देखें कि अगर मैं उस base directory तक जाता हूं। और मैंने --batch conda किया। वहां हमारी conda subdirectory है। तो यह काम कर गया और वहां हमारा cowpy output है।

तो इसने cowpy fetch किया, conda का उपयोग करके मेरे local system पर install किया, और process चलाया। और जो बढ़िया है वह यह है कि, उस end user के रूप में, मुझे वहां किसी भी software management के बारे में बिल्कुल भी सोचना नहीं पड़ा। Nextflow ने बस मेरे लिए सॉर्ट किया। मैंने कहा, मुझे इस system पर conda का उपयोग करना होगा। Pipeline developer ने कहा कि मुझे कौन से packages चाहिए। और Nextflow ने बाकी काम किया। बहुत powerful।

ध्यान दें कि तुम वास्तव में विभिन्न technologies के मिश्रण का उपयोग कर सकते हो। तो मैं specific processes के लिए Docker enable कर सकता हूं, और अन्य processes के लिए conda, या कह सकता हूं कि कुछ processes को बस जो भी local software मेरे पास installed था उसका उपयोग करना चाहिए। यह बहुत असामान्य है, लेकिन यह संभव है, और कुछ मामलों में, उदाहरण के लिए, अगर तुम कुछ ऐसे software का उपयोग कर रहे हो जिसे Docker में package करना मुश्किल हो सकता है, तो तुम्हारे पास एक escape है

## 4. एक execution platform चुनें

तो यह software packaging है। अन्य systems में portability का दूसरा भाग यह है कि jobs वास्तव में कहां चलती हैं। इस समय, मैं मूल रूप से अपने_laptop पर या इस Codespaces में चल रहा हूं, जो एक single computer है। कुछ भी fancy नहीं है। Nextflow jobs को जितना संभव हो सके parallelizing के बारे में थोड़ा clever हो रहा है, लेकिन यह सब एक system पर है।

अब, अगर तुम एक HPC पर चल रहे हो, तो तुम्हारे पास शायद किसी प्रकार का job scheduler है जैसे SLURM या PBS या कुछ, और तुम उस scheduler को jobs submit करोगे और यह सभी jobs को विभिन्न compute nodes में farm out करेगा।

चलाने का एक और तरीका cloud पर है। तो शायद तुम AWS Batch, या Azure Cloud, या Google का उपयोग कर रहे हो। और ये सभी एक समान system में काम करते हैं जहां तुम्हारे पास एक scheduler होता है और तुम jobs submit करते हो और उन्हें compute किए जाने के लिए विभिन्न स्थानों पर submit किया जाता है।

अब बहुत दूर के अतीत में जब मैंने bioinformatics शुरू किया, तो सभी के analysis चलाने के लिए software उनके computational infrastructure से बहुत tied था, जिसने इसे replicate करना लगभग असंभव बना दिया।

लेकिन Nextflow में इस config separation के साथ, और बहुत विभिन्न compute infrastructure backends के साथ interact करने की Nextflow की क्षमता के साथ, हमारी पाइपलाइन को बिल्कुल भी pipeline code को modify किए बिना लेना और बस उसे switch out करना बहुत simple है।

## 4.1. एक अलग backend को targeting करना

तो अगर हम अपनी "nextflow.config" फ़ाइल में जाते हैं, और हम अब कुछ process level config डाल सकते हैं। तो अगर मैं top पर process scope डालता हूं और मैं executor सेट कर सकता हूं, और यहां यह local पर set है, जो डिफ़ॉल्ट है।

ध्यान दें क्योंकि यह process level है, हम चीजों को विभिन्न processes पर target कर सकते हैं। और इसलिए तुम वास्तव में executors को process specific होने के लिए set up कर सकते हो और एक hybrid execution रख सकते हो, जहां कुछ jobs locally चल सकती हैं, जहां भी Nextflow job execute किया जा रहा है। कुछ को विभिन्न HPC में submit किया जाता है और कुछ को cloud में submit किया जा सकता है। तुम जितना चाहो उतना clever हो सकते हो।

अब, इस तरह के training environment में इसे demo करना बहुत मुश्किल है क्योंकि मेरे पास submit करने के लिए कोई HPC नहीं है। लेकिन मैं क्या कर सकता हूं अगर मैं slurm टाइप करता हूं, तो हम थोड़ा cheat कर सकते हैं और तुम इसका feel पा सकते हो।

और यह वास्तव में उन लोगों के लिए सबसे interesting है जो SLURM पर चलाने के आदी हैं और जानते हैं कि SLURM headers कैसे दिखते हैं। लेकिन अगर मैं Nextflow run, hello config करता हूं। यह fail होने वाला है क्योंकि यह एक cluster में jobs submit करने की कोशिश करने वाला है जो exist नहीं करता है। तो हमें sbatch उपलब्ध नहीं होने के बारे में कुछ प्रकार की error मिलेगी।

हां, written। वह tool है। वह CLI tool है जिसका उपयोग तुम slurm cluster में jobs submit करने के लिए करते हो। लेकिन हम जो कर सकते हैं वह यह है कि हम जा सकते हैं और command click करके अपनी work directory यहां देख सकते हैं, उस directory को खोल सकते हैं और .command.run को देख सकते हैं। और तुम .command.run फ़ाइल के top पर देख सकते हो, हमारे पास हमारे sbatch headers हैं, एक theoretical SLURM cluster को बताते हुए कि इस job submission को कैसे handle करना है।

और इसलिए तुम देख सकते हो कि Nextflow clever हो रहा है, यह सभी सही चीजें कर रहा है। बस हमारे पास submit करने के लिए कोई cluster नहीं था।

## 5. Compute resource allocations को control करें

विभिन्न computing infrastructures के बीच और क्या अलग है? एक और चीज यह है कि तुम्हारे पास कितने available resources हैं, और वास्तव में, कई compute environments में, यह एक requirement है कि तुम्हें specify करना होगा कि एक job को कितने CPUs और कितनी memory की जरूरत है।

फिर से, Nextflow हमारे लिए इसे abstract करता है, ताकि यह अब एक single compute environment type के लिए specific न हो, और हम यहां process level scope में टाइप कर सकें। CPUs equals one, memory equals two gigabytes। हमारी पाइपलाइन बहुत demanding नहीं है, तो वह ठीक होना चाहिए।

अब, मैंने यहां बस इन numbers का अनुमान लगाया है, लेकिन तुम कैसे जानते हो कि उपयोग करने के लिए resources की sensible मात्रा क्या है? एक बड़ी पाइपलाइन के इन सभी विभिन्न processes को कई samples के जाकर और यह समझना कि resource utilization क्या था, यह एक काफी मुश्किल काम है।

तो इसके लिए एक अच्छा approach यह है कि इन values को शुरू करने के लिए high numbers पर set करें, बस ताकि तुम्हारी पाइपलाइन बिना किसी errors के चले, और फिर Nextflow से तुम्हारे लिए एक usage report generate करने के लिए कहें।

यह करना बेहद आसान है, तो मैं एक terminal पर वापस जाने वाला हूं। ओह, मुझे याद रखना होगा कि उसे local पर वापस set करूं ताकि मेरी पाइपलाइन वास्तव में चले। और मैं Nextflow run कहने वाला हूं, और मैं एक command line flag -with-report का उपयोग करने वाला हूं।

और मैं उसे blank छोड़ सकता हूं और यह एक डिफ़ॉल्ट फ़ाइल name देगा, लेकिन मैं इसे एक specific फ़ाइल name देने वाला हूं ताकि, वह एक specific जगह पर save हो।

Enter दबाएं, और पाइपलाइन बिल्कुल सामान्य रूप से चलती है, लेकिन जब यह finish होती है, तो यह मेरे लिए एक अच्छी HTML report generate करने वाली है।

तो यहां sidebar में, मुझे यह HTML फ़ाइल मिली है। अगर मैं इसे locally चला रहा था, तो मैं इसे बस खोल देता। मैं हूं, क्योंकि मैं Codespaces में हूं, मैं उस पर right click करने वाला हूं और download पर क्लिक करने वाला हूं, जो इसे मेरे local computer पर download करने वाला है। और मैं इसे आसानी से web browser में खोल सकता हूं।

Nextflow किसी भी पाइपलाइन के लिए इस तरह की एक report generate कर सकता है और इसमें कुछ वास्तव में अच्छी जानकारी है। तो इन चीजों को हमेशा save करना अच्छा अभ्यास है। यह हमें बताता है कि हमने कब चलाया, हमने कहां चलाया, क्या यह successful था या नहीं, कौन से parameters उपयोग किए गए, CLI command क्या था, इस तरह की चीजें।

और resource usage के बारे में ये plots भी हैं। तो यह हमें बताता है कि प्रत्येक process के लिए CPU calls का कितना प्रतिशत उपयोग किया गया था एक box plot के रूप में यहां, क्योंकि प्रत्येक process के लिए कई tasks हैं, इसलिए हम distribution देख सकते हैं।

तुम यहां हमारे processes देख सकते हो, cowpy और collectGreetings में केवल एक single task था, इसलिए यह बस एक single line है। और हमारे पास CPU और memory और job duration दोनों हैं, और वे बहुत तेज थे।

अगर तुम Seqera Platform का उपयोग कर रहे हो, वैसे, तुम्हें बिना कुछ किए Platform interface में built-in same plots मिलते हैं। तो तुम्हें हमेशा यह जानकारी तुम्हारी fingertips पर मिलती है।

ठीक है, तो हम इस report का उपयोग कर सकते हैं और एक real run पर, और इसका feel पा सकते हैं कि हमारी पाइपलाइन द्वारा कितने CPUs और कितनी memory का उपयोग किया जा रहा है और वापस आएं और, और उन values को अपनी config फ़ाइल में वापस डालें, ताकि अगली बार शायद हम इतना अधिक request न करें। और हम थोड़ा अधिक lean हो सकते हैं।

अब तुम pipeline config फ़ाइलों को configure करने के बारे में वास्तव में clever हो सकते हो। और फिर से, अगर तुम Seqera Platform का उपयोग कर रहे हो, तो एक छोटे button की तलाश करो जो एक light bulb की तरह दिखता है। क्योंकि अगर तुम उस पर क्लिक करते हो, तो यह एक highly optimized config फ़ाइल generate करेगा, जो विशेष रूप से तुम्हारे data, तुम्हारे run और तुम्हारी पाइपलाइन के लिए tailored है। इसे सबसे कुशल तरीके से संभव चलाने के लिए।

लेकिन अभी के लिए, मैं कहने वाला हूं कि वास्तव में CPUs की डिफ़ॉल्ट संख्या जो Nextflow दे रहा था वह ठीक था और लेकिन केवल memory की एक gigabyte चाहिए।

## 5.3. एक specific process के लिए resource allocations सेट करें

अब, real life में, यह काफी असामान्य है कि तुम्हारी पाइपलाइन में सभी processes को same requirements की आवश्यकता होगी। तुम्हारे पास कुछ ऐसा हो सकता है जैसे MultiQC एक reporting tool के रूप में, जिसे resources के मामले में बहुत कम चाहिए और काफी जल्दी चलता है।

और फिर शायद तुम्हारे पास कुछ ऐसा है जो एक reference genome को indexing कर रहा है या कुछ alignment कर रहा है या कुछ अन्य job कर रहा है। कोई फर्क नहीं पड़ता कि यह क्या है, जो बहुत सारे resources लेता है। और इसलिए एक scheduler के लिए इन विभिन्न job submissions के लिए, तुम विभिन्न मात्रा में resources देना चाहते हो।

इस process scope के तहत, हम एक config define कर सकते हैं, जो विभिन्न तरीकों से specific processes को target करता है।

यहां हम withName का उपयोग कर रहे हैं, हम labels का भी उपयोग कर सकते हैं, और ये एक या कई processes को target करने के लिए एक pattern का उपयोग कर सकते हैं। यहां हम बस कह रहे हैं कि किसी भी processes जिनका नाम cowpy है उन्हें दो gigabytes memory और दो CPUs पर set करें, और क्योंकि यह top level process वाले की तुलना में अधिक specific selector है, यह इन मामलों में overwritten है, तो तुम यहां एक अच्छी config फ़ाइल बना सकते हो, जो वास्तव में उन्हें वास्तव में कुशल बनाने के लिए तुम्हारी पाइपलाइन में तुम्हारे सभी विभिन्न processes को tailor करती है।

## 5.5. Resource limits जोड़ें

अब एक pipeline developer के रूप में, मैं शायद tools को काफी अच्छी तरह से जानता हूं, और मैं चाहता हूं कि सब कुछ जितना संभव हो उतना तेज और अच्छी तरह से चले। तो यह हो सकता है कि मैं इनमें से कुछ के लिए काफी high numbers डालूं क्योंकि मुझे पता है कि अगर मैं cowpy को 20 CPUs देता हूं तो यह बहुत तेजी से चलेगा दो के बजाय।

यह तब तक ठीक है जब तक तुम अपने laptop पर या GitHub Actions Continuous Integration test पर, या किसी अन्य system पर चलाने नहीं जाते, जिसमें शायद 20 CPUs उपलब्ध नहीं हैं।

अब जब तुम पाइपलाइन चलाने की कोशिश करते हो, तो यह crash हो जाएगा क्योंकि Nextflow कहेगा, मैं इस job को कहीं भी submit नहीं कर सकता। मेरे पास available resources नहीं हैं।

अब उस hard crash से बचने के लिए, हम थोड़ा और config add कर सकते हैं, जो अब हमारे system के लिए specific है, जिसे resource limits कहा जाता है। और वह इस तरह दिखता है। यह फिर से process scope के अंतर्गत है।

और resource limits, तुम मूल रूप से ceiling specify कर सकते हो जो तुम्हारे पास उपलब्ध है। यहां एक map है, और तुम, इस map के भीतर, तुम memory, CPUs, और time सेट कर सकते हो।

अब क्या होता है जब Nextflow एक process से एक task submit करता है, तो यह देखता है कि क्या request किया गया है और यह मूल रूप से बस उसके और उसके बीच एक minimum करता है। तो अगर हमने 20 CPUs request किए, लेकिन केवल चार उपलब्ध हैं, तो यह चार request करेगा। पाइपलाइन crash नहीं होती है और यह pipeline developer द्वारा design किए गए के जितना करीब संभव हो उसका उपयोग करती है।

## 6. Preset configurations के बीच स्विच करने के लिए profiles का उपयोग करें

ठीक है। मैंने कहा कि यहां resource limits system specific हो सकते हैं, और शायद मेरी पाइपलाइन में मेरे पास एक Nextflow config फ़ाइल है, और मुझे पता है कि लोग इसका उपयोग विभिन्न स्थानों पर करने जा रहे हैं। अब, हर एक बार सभी को अपनी खुद की Nextflow config फ़ाइल बनाने के लिए मजबूर करने के बजाय, मैं जो कर सकता हूं वह यह है कि मैं configuration के विभिन्न presets को एक साथ config profiles में group कर सकता हूं।

मैं यहां थोड़ा नीचे scroll करने वाला हूं और फिर बस params के बाद, क्योंकि यहां config फ़ाइल का क्रम महत्वपूर्ण है, config फ़ाइल sequentially लोड की जाती है, इसलिए मैं इन profiles को बाकी सब के बाद डालने वाला हूं ताकि यह पहले defined params को override करे। और मैं training material से इन profiles को paste करने वाला हूं।

तो profiles नामक एक नया top, top level scope है। हमारे यहां arbitrary names हो सकते हैं। तो हमारे पास my_laptop और univ_hpc है। और यहां हम देख सकते हैं कि हम अन्य same config parameters सेट कर रहे हैं जो हम पहले थे। अब सिर्फ एक profile के भीतर। तो हमारे पास my_laptop पर चलाने के लिए एक local executor है और मैं HPC पर एक SLURM cluster में submit कर रहा हूं।

मैं locally Docker का उपयोग कर रहा हूं, HPC पर conda, और HPC system में बहुत अधिक resource limits हैं।

अब मैं -profile CLI option के साथ पाइपलाइन चला सकता हूं, कहें कि मैं कौन सा profile उपयोग करना चाहता हूं। तो मैं my_laptop का उपयोग करने वाला हूं, और Nextflow उस profile scope के भीतर सभी config को apply करेगा। तो मैं अब इसे try कर सकता हूं। यह पहले जैसी same command है। Nextflow run hello config, और मैं dash profile करता हूं, single dash क्योंकि यह core Nextflow option है, dash profile my_laptop।

यह अब उस सभी config option को batch apply करने वाला है। ओह, और तुम देख सकते हो, मैंने पहले कहा था कि यह हो सकता है कि process requirement, इसने चार CPUs मांगे और मेरे पास इस Codespaces instance पर केवल दो हैं।

तो यह process resource limits को try out करने का एक अच्छा अवसर है, और कहें कि मेरे my_laptop पर, या इस Codespaces में केवल दो CPUs हैं। अब अगर हम इसे फिर से चलाते हैं, तो इसे उस requirement को दो पर cap करना चाहिए और उम्मीद है कि पाइपलाइन चलेगी। बढ़िया।

## 6.2. test parameters का एक profile बनाएं

ध्यान दें कि इन profiles में केवल उनके infrastructure के बारे में configuration होना जरूरी नहीं है। तुम यहां किसी भी config की groupings रख सकते हो, parameters सहित।

तो एक और चीज जो तुम अक्सर लोगों की pipelines में देखोगे वह है एक test profile, जिसमें parameters शामिल होते हैं, जिन्हें तुम आमतौर पर per user basis पर submit करोगे। लेकिन यहां हमारे पास, मूल रूप से अलग-अलग sensible defaults हैं जब मैं test cases चलाना चाहता हूं।

और यह बढ़िया है क्योंकि मुझे इन सभी चीजों को जरूरी specify करने की जरूरत नहीं है, जो required parameters हो सकती हैं। अन्यथा मैं बस dash profile test कह सकता हूं और यह बस out of the box चलेगा।

अब ध्यान देने वाली बात यह है कि profiles को एक से अधिक भी combined किया जा सकता है। तो मैं यहां profile my_laptop कर सकता हूं, और फिर test भी add कर सकता हूं। मैं profile दो बार नहीं करता। मैं बस यहां बिना spaces के comma separated list करता हूं। और यह इन profiles को क्रम में apply करने वाला है। तो यह my_laptop profile से config लेगा, और फिर यह top पर test config apply करेगा।

वास्तव में convenient है और तुम देख सकते हो कि तुम यहां बहुत सारे sensible default groups कैसे set up कर सकते हो ताकि तुम्हारी पाइपलाइन को चलाना आसान हो।

## 6.3. resolved configuration देखने के लिए nextflow config का उपयोग करें

उम्मीद है, मैंने तुम्हें convince कर लिया है कि Nextflow config resolution powerful है, लेकिन मैं तुम्हें दोष नहीं दूंगा अगर तुम इस बिंदु पर थोड़ा cross-eyed हो रहे हो जब मैंने लगभग 20 विभिन्न तरीकों के बारे में कहा है config प्रदान करने और इन सभी विभिन्न layers को देने के बारे में जैसे एक onion skin।

तो अगर कभी तुम अनिश्चित महसूस कर रहे हो कि Nextflow के लिए final resolved config क्या है, तो जान लो कि "nextflow config" नामक एक command है, और हम इसे चला सकते हैं और यह हमें बताएगा कि हमारे current location पर resolved configuration क्या है।

तो जब मैं इसे यहां चलाता हूं, तो यह current working directory में "nextflow.config" फ़ाइल find करता है, और यह सभी विभिन्न config को process करता है, और यह मुझे resolved output देता है।

ध्यान दें कि Nextflow config फ़ाइल profile CLI option भी ले सकती है। तो अगर मैं इसे my_laptop और test profiles में resolve करने के लिए कहता हूं, और तुम देख सकते हो कि इसने my_laptop config option से resource limits भी apply की है और params भी set किए हैं जो test में थे।

तो यह सिर्फ explore करने का एक अच्छा तरीका है कि config resolution कैसे काम कर रहा है, अगर तुम बिल्कुल भी अनिश्चित हो।

## Wrap up

ठीक है, बस इतना ही। यह एक nutshell में Nextflow config है। तुम config के साथ बहुत कुछ कर सकते हो। यह वास्तव में powerful है। लेकिन ये अधिकांश सामान्य use cases हैं जिन्हें तुम खुद को करते हुए पाओगे, और ये concepts सभी विभिन्न options पर apply होते हैं।

अपनी पीठ थपथपाओ क्योंकि यह Hello Nextflow training course का अंत है। तुम अब उम्मीद से scratch से अपनी खुद की Nextflow पाइपलाइन लिखने, इसे configure करने और चलाने में confident हो, और तुम सभी ins और outs और ध्यान रखने वाली चीजें जानते हो।

एक और quiz है जिसे तुम config training page पर try कर सकते हो। तो नीचे जाओ और उसे try करो और सुनिश्चित करो कि तुमने config के बारे में इन सभी भागों को समझ लिया है।

और, अगले steps के बारे में एक quick wrap up के लिए अंतिम video में हमारे साथ शामिल हो जो इस training course के बाद करना अच्छा हो सकता है।

हमारे साथ बने रहने के लिए धन्यवाद। बहुत बढ़िया और मैं तुम्हें अगले video में मिलूंगा।
