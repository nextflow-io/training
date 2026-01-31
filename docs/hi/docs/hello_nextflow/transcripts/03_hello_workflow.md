# भाग 3: Hello Workflow - ट्रांसक्रिप्ट

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "महत्वपूर्ण नोट"

    यह पेज केवल ट्रांसक्रिप्ट दिखाता है। पूर्ण चरण-दर-चरण निर्देशों के लिए, [पाठ्यक्रम सामग्री](../03_hello_workflow.md) पर वापस जाएं।

    ट्रांसक्रिप्ट में दिखाए गए अनुभाग नंबर केवल संकेत के लिए प्रदान किए गए हैं और सामग्री में सभी अनुभाग नंबर शामिल नहीं हो सकते हैं।

## स्वागत

नमस्ते, "Hello Nextflow" प्रशिक्षण पाठ्यक्रम के भाग तीन में आपका स्वागत है।

इस अध्याय को "Hello Workflow" कहा जाता है।

अध्याय दो में, हमने एक process की एक सरल workflow बनाई थी, लेकिन वास्तविकता में, pipelines उपयोगी होती हैं क्योंकि वे विश्लेषण के कई चरणों को एक साथ जोड़ सकती हैं।

इस अध्याय में, हम उस प्रारंभिक उदाहरण को लेने जा रहे हैं और इसे थोड़ा और यथार्थवादी बनाने के लिए विस्तारित करेंगे।

हम कुछ अतिरिक्त चरण जोड़ने जा रहे हैं और हम देखेंगे कि उन चरणों को जोड़ने के लिए हम channels का उपयोग कैसे करते हैं।

हम कई tasks को देखने जा रहे हैं, जो एक single process में समा सकते हैं और हम processes को देखने जा रहे हैं जिनमें कई inputs और कई outputs हो सकते हैं।

ठीक है, चलिए शुरू करते हैं।

तो चलिए शुरू करते हैं। पहले की तरह ही। चलिए training.nextflow.io पर जाते हैं। Hello Nextflow, अध्याय तीन। Hello Workflow। और चलिए अपना workspace खोलें। मैंने अपने पिछले अध्यायों से अपनी सभी work files को साफ कर दिया है और मैं Hello Workflow को खोलने जा रहा हूं।

अब यह वही फ़ाइल है जिस पर हम अब तक काम कर रहे थे इसलिए यह परिचित लगनी चाहिए। हमारे पास हमारी say hello process है। हमारे पास हमारी params.greeting अपनी greetings CSV फ़ाइल के साथ है, और हमारे पास नीचे हमारी workflow है, जो उस CSV फ़ाइल को लोड करती है, channel बनाती है और इसे हमारी process में पास करती है।

## 0. वार्म-अप: hello-workflow.nf चलाएं

यदि आप चाहें, तो हम इसे आज़मा सकते हैं और दोबारा जांच सकते हैं कि यह हमारी उम्मीद के अनुसार काम कर रहा है। nextflow run hello workflow nf के लिए एक terminal लोड करें और enter पर क्लिक करें।

ठीक है, बढ़िया। हमारी तीन processes रन हुईं। हमारे पास अपनी तीन outputs के साथ हमारी results डायरेक्टरी है। Bonjour. Hello. Holà. तो चलिए उन files को बंद करें, terminal को बंद करें, script पर वापस जाएं।

## 1. workflow में दूसरा चरण जोड़ें

ठीक है। हमारे उदाहरण के लिए, हम बुनियादी बने रहेंगे और हम domain agnostic रहने की कोशिश कर रहे हैं। तो हमारी दूसरी process बस इन strings, इन शब्दों को एक सरल तरीके से manipulate करने जा रही है। हम इन files को लेने और उन्हें सभी uppercase बनाने के लिए translate Unix कमांड का उपयोग करने जा रहे हैं। हम "tr" कमांड के साथ ऐसा करते हैं।

## 1.1. uppercasing कमांड परिभाषित करें और इसे terminal में टेस्ट करें

हम इसे बस bash terminal में आज़मा सकते हैं, और देख सकते हैं कि क्या यह काम करता है। तो आप echo करें, Hello World, और फिर इसे pipe character के साथ tr पर पास करें, और हम इसे एक recognition pattern देते हैं, a से z और इसे क्या translate करना चाहिए। A से Z uppercase में।

यह बहुत सरल है क्योंकि यह शाब्दिक रूप से A से Z characters कर रहा है। तो यह किसी भी चीज़ पर काम नहीं करेगा जो accented है या ऐसा कुछ भी है। लेकिन उदाहरण के उद्देश्यों के लिए, आपको तस्वीर मिल जानी चाहिए।

enter दबाने जा रहा हूं और यह terminal पर प्रिंट करता है, HELLO WORLD capitals में। और बस पहले की तरह, हम चाहें तो इसे एक फ़ाइल में redirect कर सकते हैं। Outfile.

ठीक है। चलिए इसे साफ करें।

## 1.1. uppercasing step को Nextflow process के रूप में लिखें

चलिए अपनी script पर वापस जाएं और इस bash कमांड को handle करने के लिए एक नई process लिखें। मैं पिछली process को copy करने जा रहा हूं, इसे नीचे paste करूंगा, और इसे convert to upper कहूंगा। Uppercase के लिए। मैं वही publishDir results का उपयोग करने जा रहा हूं, लेकिन मैं यहां कुछ बदलाव करने जा रहा हूं। val लेने के बजाय, मैं एक path input फ़ाइल लेने जा रहा हूं, और मेरे पास यहां एक prefix upper होगा, ताकि हमारी output files आउटपुट को clobber न करें। और मैं input से variable name का उपयोग करने जा रहा हूं। और फिर मैं यहां नीचे एक script बदलने जा रहा हूं, और इसके बजाय मैं input फ़ाइल पर cat का उपयोग करने जा रहा हूं और जैसे हमने Bash TR में किया था, a-z, upper input file .txt। ठीक है, चलिए save पर क्लिक करें।

## 1.2. workflow block में नई process को कॉल जोड़ें

अब यदि मैं नीचे scroll करूं, तो हमें वास्तव में इस process को call करने की आवश्यकता है। बस script में process जोड़ना पर्याप्त नहीं है। हमें Nextflow को बताना होगा कि हमें इस process को चलाने की आवश्यकता है और इसे कहां करना है।

तो मैं यहीं जा रहा हूं, convert to upper और

ठीक है, हमें यहां एक error मिल रही है जो कहती है कि यह एक argument की उम्मीद करती है। निश्चित रूप से, हमें इस process को कुछ पास करने की आवश्यकता है ताकि इसके पास वास्तव में कुछ करने के लिए हो।

## 1.3. पहली process के output को दूसरी process में पास करें

हम जो करने जा रहे हैं वह है हम इस process से output लेने जा रहे हैं। तो मैं name लेता हूं, say hello, और जब मैं dot out करता हूं।

इस तरह के एक सरल उदाहरण के लिए, जहां हमारे पास एक process है जिसमें बस एक output है और हम इसे एक नई process में पास कर रहे हैं, तो इसमें एक input है जो हमें चाहिए होना चाहिए। तो मैं save पर क्लिक करने जा रहा हूं, terminal लाऊंगा, और चलिए इसे फिर से चलाने की कोशिश करते हैं।

## 1.4. workflow को फिर से चलाएं

अब, मैंने पिछली बार जब मैंने यह workflow चलाई थी तब से अपनी work डायरेक्टरी को साफ नहीं किया है। मैं इसे फिर से चलाने जा रहा हूं और मैं इसे एक अवसर के रूप में उपयोग करने जा रहा हूं यह दिखाने के लिए कि partial caching कैसे काम करती है। तो यदि मैं single dash resume करता हूं। उम्मीद है कि यह उस पहली process से outputs का पुन: उपयोग करना चाहिए, जो पिछली बार मैंने चलाए के बिल्कुल समान थे। लेकिन अब हमारे पास यहां एक नई process है जो पहले नहीं चली है, जो scratch से चलती है। और निश्चित रूप से, आप देख सकते हैं कि पहली process ने cache outputs का उपयोग किया, और दूसरे output ने तीन में से तीन चलाए। आप यह भी देख सकते हैं कि हमारे पास अब अपनी दोनों processes हैं, हमारी पहली process, say hello, तीन बार चली, और हमारी दूसरी process convert to upper तीन बार चली।

यदि मैं इसे फिर से चलाता हूं, एक reminder के रूप में, -ansi-log false के साथ, हमें देखना चाहिए कि छह अलग-अलग process tasks चलीं उनमें से प्रत्येक के लिए तीन। तो यह बिल्कुल वही कर रहा है जो हमने उम्मीद की थी। पहली process तीन बार चल रही है, उन outputs को दूसरी process पर पास कर रही है, जो फिर तीन बार चल रही है।

तो चलिए work डायरेक्टरी के अंदर देखते हैं और देखते हैं कि Nextflow इन file inputs को कैसे handle कर रहा है। यदि मैं दूसरी process से इस hash डायरेक्टरी को यहां लेता हूं तो हम इन files को देखने के लिए फिर से -a के साथ tree कमांड का उपयोग कर सकते हैं। आप यहां देख सकते हैं कि हमारे पास हमारी input फ़ाइल है, जो Bonjour-output.txt फ़ाइल है, और वह वास्तव में एक symlink है। यही वह है जो यह arrow हमें दिखा रहा है, और यह पिछली work डायरेक्टरी में फ़ाइल की ओर इशारा कर रहा है।

यह समझ में आता है। Nextflow अपनी स्वयं की encapsulated डायरेक्टरी में प्रत्येक कार्य के execution को handle करता है, इसलिए यह पूरी तरह से self enclosed है। हालांकि, इसे input के रूप में पिछले चरणों से files प्रदान करने की आवश्यकता है। उन files को प्राप्त करने के लिए work डायरेक्टरी के बाहर पहुंचने के बजाय, Nextflow उन्हें work डायरेक्टरी में stages करता है।

यदि हमारे पास यहां जैसी shared file system है, तो यह एक symlink का उपयोग करके ऐसा करता है ताकि यह किसी अतिरिक्त file space का उपयोग न करे। यदि हम विभिन्न स्थानों में buckets के साथ cloud storage का उपयोग करते हैं, तो यह उन files को fetch करेगा और वास्तव में उन्हें work डायरेक्टरी में copy करेगा।

चलिए command sh फ़ाइल पर एक नज़र डालते हैं। यदि मैं code work, command sh करता हूं, तो आप देख सकते हैं, निश्चित रूप से, यह local डायरेक्टरी से उस फ़ाइल तक पहुंच रहा है। तो सब कुछ बहुत self-contained और clean है।

हम results डायरेक्टरी की भी जांच कर सकते हैं और सुनिश्चित कर सकते हैं कि ये files ठीक से output हुईं। और निश्चित रूप से, results में, हम पहली process से सभी output files और दूसरी से सभी output files देख सकते हैं। और वे सभी uppercase में हैं जैसा हमने उम्मीद की थी।

यहीं पर Nextflow की शक्ति चमकने लगती है। कुछ बहुत न्यूनतम code के साथ और Nextflow ने इन tasks के parallel में execution को clean encapsulation के साथ अलग work directories के भीतर और input और output files की staging और file publishing सभी स्वचालित रूप से हमारे लिए बस out of the box में handle किया। तो आप देख सकते हैं कि जैसे-जैसे हम अपनी analysis workflows की जटिलता को scale करते हैं, यह functionality वास्तव में, वास्तव में मूल्यवान है।

## 2. सभी greetings को collect करने के लिए तीसरा चरण जोड़ें

ठीक है। ये steps one-to-one थे। हमारे पास पहली process से एक output था जो दूसरी process के लिए एक input में जा रहा था। इसके बाद, हम बात करने जा रहे हैं कि इन विभिन्न outputs को एक single process कार्य में कैसे collect किया जाए, जो फिर से, करने के लिए एक बहुत ही सामान्य बात है। तो चलिए जल्दी से terminal लाते हैं और इसका एक dry run करते हैं।

## 2.1. collection कमांड परिभाषित करें और इसे terminal में टेस्ट करें

मैं cheat करने जा रहा हूं और प्रशिक्षण सामग्री से example bash code को copy करूंगा और बस enter दबाऊंगा।

हम यहां देख सकते हैं कि हमने इस echo कमांड को तीन अलग-अलग output files के लिए तीन बार चलाया, जिसे मैं यहां देख सकता हूं। और फिर इन तीन अलग-अलग files में से प्रत्येक के output को print करने के लिए cat कमांड का उपयोग किया, और उसे एक single collected फ़ाइल में redirect किया।

और यदि मैं "cat COLLECTED-output" करता हूं, तो आप देख सकते हैं कि इसमें उन तीन अलग-अलग files की contents हैं, अब एक single फ़ाइल में।

## 2.2. collection step करने के लिए एक नई process बनाएं

तो चलिए देखते हैं कि क्या हम अपनी Nextflow pipeline के भीतर वही चीज़ replicate कर सकते हैं।

चलिए ऊपर scroll करें और एक तीसरी process बनाएं। मैं इस पिछले वाले को copy करने जा रहा हूं, और इस बार मैं इसे Collect Greetings कहने जा रहा हूं।

bash terminal में, हमने इसे collected output txt कहा था। तो मैं यहां same path output कहने जा रहा हूं। और मैं यहां redirection करने जा रहा हूं, इसलिए यह उसी तरह से save हो जाती है।

ठीक है। हमें उस कमांड की शुरुआत में जो होता है उसे बदलने की आवश्यकता है, और हमें सोचना होगा कि यहां input फ़ाइल क्या है। वास्तव में, यह process कई input files लेने जा रही है। मैं path रखने जा रहा हूं और मैं इसे एक नए variable input files में बदलने जा रहा हूं, plural।

फिर मैं फिर से, उन्हें cat करूंगा जैसे हमने अपनी bash script में किया था। और मैं यहां variable का उपयोग करने जा रहा हूं।

अब, आप सोच सकते हैं कि यह काम नहीं करेगा। हमने पहले failures देखी हैं जहां strings की एक array या paths की एक array को एक process में पास किया गया था और इससे error हुई थी। लेकिन वास्तव में, यहां Nextflow इसे स्वचालित रूप से सही तरीके से हमारे लिए handle करने जा रहा है। यह कई अलग-अलग input files लेने जा रहा है, और यह बस यहां विभिन्न file paths को print करने जा रहा है।

बेशक यह मदद करता है कि cat कमांड इस तरह file names की एक series ले सकती है। यदि मैं एक अलग कमांड का उपयोग कर रहा था जिसे प्रत्येक file path से पहले एक argument की आवश्यकता थी या कुछ और, तो हमें यहां थोड़ा और अधिक code और logic रखना होगा ताकि इन file paths की iteration को handle करने में सक्षम हो सकें। लेकिन इस मामले में, यह बस काम करना चाहिए।

## 2.3. workflow में collection step जोड़ें

ठीक है, चलिए workflow पर नीचे जाएं और अपनी नई process जोड़ें। Collect greetings। और फिर से, चलिए convert to upper out से output लेते हैं। चलिए इसे save करें।

इसे एक try दें। nextflow run hello workflow.

ठीक है, workflow चली, लेकिन कुछ अजीब है यहां। हमें पहले चरण के तीन executions मिले हैं, जिसकी हम उम्मीद करते हैं। दूसरे के लिए तीन tasks, लेकिन हमारे पास अंत में तीन tasks भी हैं जब हम यहां केवल एक single कार्य की उम्मीद करते थे जो सभी outputs को merge कर रहा था।

यदि हम अपनी results डायरेक्टरी में जाते हैं। हम यह भी देखते हैं कि collected output में तीनों के बजाय केवल एक single value है। ऐसा इसलिए है क्योंकि वह output फ़ाइल तीन अलग-अलग values के साथ तीन बार overwrite की गई थी।

यह समझ में आता है क्योंकि हमने यहां एक output को एक input में उसी तरह पास किया जैसे हमने पिछले चरण में किया था।

## 2.4. greetings को एक single input में collect करने के लिए एक operator का उपयोग करें

तो हमें यहां एक operator की आवश्यकता है जो इस channel को तीन elements के साथ लेगा और उन्हें एक single element में collapse करेगा, ताकि वह final process केवल एक बार चले।

ऐसा करने के लिए, हम collect operator का उपयोग करने जा रहे हैं। मैं workflow के भीतर सीधे ऐसा कर सकता हूं। मैं .out कर सकता हूं और यहां अंत में एक operator पर chain कर सकता हूं .collect.

save दबाएं। और फिर इस प्रशिक्षण के उद्देश्यों के लिए, मैं कुछ view operators भी करने जा रहा हूं जैसे हमने पहले किया था, ताकि हम collect operator का उपयोग करने से पहले और बाद में इस channel पर एक नज़र डाल सकें, ताकि हम समझ सकें कि क्या हो रहा है।

मैं इस channel को लेने जा रहा हूं, collect से छुटकारा पाऊंगा और dot view greetings करूंगा, और फिर मैं इस line को duplicate करूंगा, collect operator जोड़ूंगा। और उसे after में बदलूंगा।

यह अलग है जहां हम इसे call कर रहे हैं, लेकिन यह ठीक है क्योंकि हम same output channel पर same operator calls का उपयोग कर रहे हैं।

ठीक है, चलिए save दबाएं और चलिए इसे terminal में आज़माएं। lgoing nextflow run करने जा रहा हूं। Hello, workflow। हमारी script फिर से चलाएं।

ठीक है। यह बेहतर लग रहा है। पहले की तरह हम देख सकते हैं कि पहली दो processes तीन बार चलीं और अब हमारी final process केवल एक बार चली।

यदि हम देखें कि view operator द्वारा क्या print किया गया था, यहां नीचे, हमने collect से पहले कहा, जो यहां यह output है, और वह तीन बार print हुआ है। और आप देख सकते हैं कि उनमें से प्रत्येक के लिए एक single path है। और फिर collect के बाद, आप देख सकते हैं कि हमारे पास तीन paths की यह array है। तो यह जैसी हमने उम्मीद की थी वैसी है।

ठीक है, चलिए results फ़ाइल की जांच करें और देखें कि क्या यह इस बार हमारी उम्मीद के अनुसार है। निश्चित रूप से, अब फ़ाइल में तीन lines हैं - वह सफलतापूर्वक इन तीन outputs को एक single output फ़ाइल में concatenate कर दी। शानदार।

ठीक है, मैं clean up करने जा रहा हूं और चलिए अगले चरण पर जाते हैं। और मैं इन view statements को हटाने जा रहा हूं बस चीजों को clean रखने के लिए।

## 3. final output फ़ाइल को uniquely नाम देने के लिए एक process को एक से अधिक input पास करें

ठीक है। अब तक, हमारी सभी processes ने केवल एक single input लिया है। अब हम एक अभ्यास करने जा रहे हैं जहां हम एक process में एक से अधिक input जोड़ते हैं यह देखने के लिए कि यह कैसे काम करता है। ऐसा करने के लिए, हम इस collect greetings example का उपयोग करने जा रहे हैं।

हर बार जब मैंने workflow चलाई, तो यह results डायरेक्टरी में उस फ़ाइल को overwrite कर दी, जो शायद हम नहीं चाहते हैं।

## 3.1. output फ़ाइल के लिए user-defined name स्वीकार करने के लिए collector process को संशोधित करें

तो इस उदाहरण के लिए, हम एक अतिरिक्त पैरामीटर पास करने जा रहे हैं ताकि हम output फ़ाइल के name को customize कर सकें।

एक process में दूसरा input जोड़ना बहुत सरल है। मैं बस input block में दूसरी line जोड़ता हूं। इस बार यह एक value होने जा रहा है, path के बजाय, क्योंकि हम एक string पास करना चाहते हैं और मैं इसे batch underscore name कहने जा रहा हूं।

अब मैं इस variable का उपयोग script block में कर सकता हूं, और मैं collected dash dollar batch name कहने जा रहा हूं।

मैं यहां variable name के आसपास squiggly brackets का उपयोग कर रहा हूं। यह बस इसे string के बाकी हिस्सों से अलग रखने के लिए है, और शायद इस मामले में इसकी आवश्यकता नहीं है, लेकिन मुझे लगता है कि यह इसे पढ़ना आसान बनाता है।

ठीक है। अंत में, output path को update करना याद रखें क्योंकि अब फ़ाइल का name बदल गया है, तो मैं वही चीज़ करने जा रहा हूं और batch name को output of path में उम्मीद के अनुसार रखूंगा।

## 3.2. batch command-line पैरामीटर जोड़ें

अब हमें कहीं से एक batch name पास करने की आवश्यकता है, और मैं ऐसा करने के लिए एक दूसरा पैरामीटर बनाने जा रहा हूं ताकि हम इसे command line पर कर सकें जब हम workflow चलाएं।

तो मैं params batch name करने जा रहा हूं, और by default, चलिए इसे test batch कहते हैं। अब मैं इस special parameters variable को नीचे उपयोग कर सकता हूं, जहां हम process को call करते हैं।

और निश्चित रूप से VS Code हमें बता रहा है कि अब इस process के लिए पर्याप्त arguments नहीं हैं, और यह दूसरे input की उम्मीद करती है।

बस comma करें और हमारे नए variable को पास करें और error चली जाती है।

ध्यान दें कि यहां inputs का क्रम वास्तव में महत्वपूर्ण है। पहला process input path था, और दूसरा input name है। यदि मैं यहां क्रम बदलता हूं, तो मुझे भी क्रम बदलना होगा जब मैं process को call करूं। अन्यथा। इसके बाद, हम गलत input के लिए गलत channel पास करेंगे।

## 3.3. workflow चलाएं

ठीक है, चलिए इसे आज़माते हैं और देखते हैं कि क्या यह काम करता है। चलिए "nextflow run hello- workflow करते हैं। ठीक है, यह पहले की तरह चली। चलिए results डायरेक्टरी में देखते हैं।

निश्चित रूप से, हमारी फ़ाइल का name अब "collected test batch output txt" है। शानदार।

और अब देखते हैं कि क्या हम फिर से चलाकर उसे overwrite कर सकते हैं। इस बार मैं --batch_name करने जा रहा हूं ताकि यहां उस special parameter variable name से match हो। और मैं इसे demo output कहने जा रहा हूं।

workflow को फिर से चलाएं और हम देखेंगे कि क्या कुछ होता है।

ठीक है, अब हमारे पास एक collected demo output .txt है। और क्योंकि यह फ़ाइल का name उस एक से अलग है, इसने इसे overwrite नहीं किया है। दोनों अब results डायरेक्टरी में मौजूद हैं।

## 4. collector step में एक output जोड़ें

ठीक है, तो वहां हमने एक process को कई inputs देना दिखाया, लेकिन कई outputs के बारे में क्या? इस उदाहरण के लिए, हम greetings की संख्या की गणना करने जा रहे हैं जो processed हैं और इसे इस collect greeting step के लिए एक secondary output के रूप में output करें।

## 4.1. greetings की संख्या को count और output करने के लिए process को संशोधित करें

हम यहां थोड़ी trick करने जा रहे हैं। Nextflow processes में यह script block होता है एक multi-line string के साथ, और वह bash output के रूप में dot कमांड dot sh में पास किया जाता है। लेकिन हम वास्तव में उस से ऊपर कोई भी custom code लिख सकते हैं, और वह एक कार्य के हिस्से के रूप में execute होगा लेकिन bash script के भीतर शामिल नहीं होगा।

Nextflow syntax में built-in functions में से एक को size कहा जाता है। तो मैं path input लेने जा रहा हूं, और मैं count underscore greetings कहने जा रहा हूं, बस एक variable name define करने के लिए। मैं input files लेने जा रहा हूं और मैं इस पर "size" call करने जा रहा हूं।

यह function इस input channel के size को count करेगा और इसे एक variable को assign करेगा।

अब हम उस variable को output block के हिस्से के रूप में return कर सकते हैं। तो हम कहते हैं, val, क्योंकि यह value है, फ़ाइल नहीं। और count greetings.

अब यह अपने आप में पर्याप्त है, और हम अब इस process से इन विभिन्न outputs तक पहुंच सकते हैं। हालांकि, हमें उन्हें positional तरीके से access करना होगा। तो zero और one जैसी index key का उपयोग करके।

outputs तक पहुंचना थोड़ा आसान बनाने के लिए, हम उन्हें नाम दे सकते हैं और हम एक emit statement का उपयोग करके ऐसा करते हैं।

तो हम comma emit out file या जो कुछ भी मैं इसे call करना चाहता हूं करते हैं। और मैं यहां emit count करता हूं। यह मूल रूप से बस एक decorator है, जो हमें थोड़ा cleaner code लिखने में मदद करता है ताकि हम बाद में workflow block में specific outputs को आसानी से reference कर सकें।

## 4.2. workflow के अंत में output की रिपोर्ट करें

ठीक है। यदि मैं workflow block पर नीचे scroll करूं, तो अब मैं collect greetings के outputs ले सकता हूं, collect greetings, dot out कर सकता हूं, और हम देख सकते हैं कि हमारे दो named outputs यहां VS Code extension द्वारा सुझाए गए हैं। बहुत उपयोगी।

तो मैं dot count करने जा रहा हूं जो हमने अभी बनाया count value प्राप्त करने के लिए, और मैं view करने जा रहा हूं, ताकि यह command line में print हो। ताकि हम इसे देख सकें जब हम workflow चलाएं।

चलिए closure में यहां कुछ लिखते हैं बस इसे थोड़ा अच्छा बनाने के लिए। num greetings, there were greetings greetings.

और हम वास्तव में दूसरे output के बारे में परवाह नहीं करते क्योंकि हम इसे किसी अन्य processes के लिए input के रूप में उपयोग नहीं कर रहे हैं। लेकिन आप देख सकते हैं कि यदि हम चाहें तो इसे एक और process के लिए input के रूप में आसानी से कैसे पास कर सकते हैं, downstream।

## 4.3. workflow चलाएं

हम save पर क्लिक करने जा रहे हैं। चलिए terminal पर एक नज़र डालते हैं और इसे आज़माते हैं।

ठीक है, शानदार। यहां हम जाते हैं। तीन greetings हैं। यह बिल्कुल सही है।

ठीक है, बढ़िया सामान। यह इस अध्याय का अंत है। अब तक बनाने के लिए हम सब पूरा कर चुके हैं। अब आप काफी realistic workflow बनाना शुरू कर रहे हैं, जहां हम अपनी workflow के भीतर inputs और outputs और logic को handle करने में सक्षम हैं।

जैसे-जैसे ये workflow files लंबी होती जाती हैं, वे थोड़ी unwieldy होने लगती हैं। तो अगले अध्याय में, हम देखेंगे कि कैसे हम Nextflow code को अलग files में modularize कर सकते हैं ताकि workflow के भीतर code को ढूंढना और बनाए रखना आसान हो।

अध्याय चार के लिए अगले वीडियो में हमारे साथ शामिल हों। Hello Modules.

[अगला वीडियो ट्रांसक्रिप्ट :octicons-arrow-right-24:](04_hello_modules.md)
