# भाग 1: Hello World - वीडियो ट्रांसक्रिप्ट

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "महत्वपूर्ण नोट"

    यह पेज केवल ट्रांसक्रिप्ट दिखाता है। पूर्ण चरण-दर-चरण निर्देशों के लिए, [कोर्स सामग्री](../01_hello_world.md) पर वापस जाएं।

    ट्रांसक्रिप्ट में दिखाए गए सेक्शन नंबर केवल संकेत के उद्देश्य से दिए गए हैं और हो सकता है कि सामग्री में सभी सेक्शन नंबर शामिल न हों।

## स्वागत

नमस्ते, और वापस स्वागत है।

तुम अब "Hello Nextflow" कोर्स के भाग एक में हो जिसे "Hello World" कहा जाता है। इस अध्याय में, हम Nextflow की बुनियादी बातों की समझ बनाना शुरू करेंगे।

तो उम्मीद है कि तुम अब Codespaces या किसी समकक्ष जगह पर VS Code के साथ सेटअप हो, और तुम्हारे पास Explorer में workspace में Hello Nextflow फ़ोल्डर है जिसमें ये सभी अलग-अलग फ़ाइलें हैं।

हम terminal में Bash का उपयोग करके कुछ बहुत ही बुनियादी चीजें करके शुरू करेंगे, और फिर हम देखेंगे कि क्या हम Nextflow के भीतर वही चीजें कर सकते हैं ताकि तुम्हें syntax का अनुभव मिल सके।

## 0. वार्मअप

तो चलो वास्तव में सरल शुरू करते हैं। चलो बस "echo" से शुरू करते हैं, terminal पर कुछ प्रिंट करने के लिए। "Hello World"। मैं enter दबाता हूं और वह terminal पर जाता है। Hello World। उम्मीद है कि यह इस कोर्स को देख रहे किसी के लिए आश्चर्य की बात नहीं है।

ठीक है, चलो इसके साथ कुछ करते हैं। इसे सिर्फ terminal पर प्रिंट करने के बजाय, चलो इसे एक फ़ाइल में लिखते हैं। मैं अपने keyboard के up cursor को दबाने जा रहा हूं, जो Bash history के माध्यम से cycle करता है, तो यह मुझे मेरा आखिरी command देता है, और मैं इसके अंत में वहां जोड़ने जा रहा हूं, छोटा greater than symbol, जो इस command से output को एक फ़ाइल में redirect करता है, और मैं इसे output.txt कहने जा रहा हूं।

फिर से enter, उस command को चलाने के लिए, इस बार terminal में कुछ नहीं, लेकिन हम बाईं ओर देख सकते हैं, नई फ़ाइल यहां दिखाई दी है, जिसे output.txt कहा जाता है।

हम इसे terminal में cat जैसी किसी चीज़ से देख सकते हैं। तो cat output.txt और निश्चित रूप से यह "Hello World" कहता है। हम इसे double click भी कर सकते हैं और यह VS Code में code editor में खुल जाता है।

## 1.1. कोड की जांच करें

ठीक है। मैंने तुमसे कहा था कि यह सरल था। आगे क्या है? चलो इस प्रोसेस को लेने की कोशिश करते हैं और इसे फिर से करते हैं, लेकिन इस बार, चलो इसे Nextflow के अंदर करते हैं।

जैसा कि मैंने कहा, इस कोर्स के सभी अलग-अलग अध्याय एक script के साथ शुरू होते हैं और इसे Hello World कहा जाता है। तो मैं Hello World ढूंढने जा रहा हूं। यह इसे preview करता है अगर मैं इसे single click करता हूं, मैं इसे यहां editor में खोलने के लिए double click करने जा रहा हूं। और मैं बस जल्दी से terminal से छुटकारा पाने जा रहा हूं।

अब यह एक बहुत ही सरल script है, तो जितना सरल हो सकता है। यह केवल 22 लाइनों लंबा है, और यह मूल रूप से वही काम करता है। वास्तव में। इसमें से कुछ परिचित लगना चाहिए। यह वही है जो हमने अभी टाइप किया था। हम अपना bash command देख सकते हैं जो वहां एक फ़ाइल में redirect कर रहा है।

ठीक है। और क्या? इसके अलावा, इस फ़ाइल में, हम Nextflow की कुछ मुख्य अवधारणाओं को देखना शुरू कर सकते हैं। हमारे पास यहां लाल रंग में एक process है और एक workflow। ये Nextflow में विशेष keywords और विशेष terminology हैं।

## 1.1.1. प्रोसेस परिभाषा

एक workflow के भीतर अलग-अलग processes तुम्हारे workflow की अलग-अलग logical units को wrap करती हैं। प्रत्येक process एक काम करती है।

जब हम इसे चलाते हैं, तो यह एक task या कई tasks उत्पन्न करती है, जो pipeline के वास्तविक कार्य चरण हैं। सभी processes को फिर एक workflow block के भीतर orchestrate किया जाता है, जिसे हम नीचे देखते हैं, और इस मामले में बस उस एक process को चलाता है।

process का नाम यहां इस keyword के बाद आता है, और यह मूल रूप से कुछ भी हो सकता है। और फिर process की सामग्री इन curly brackets के भीतर होती है।

process के लिए वास्तव में केवल एक आवश्यकता है, जो यह है कि इसमें किसी प्रकार का script या exec block शामिल हो। यह यहां triple quotes में है, और यह bash script है जो working directory में लिखी जाती है जब हम pipeline चलाते हैं और यह वह चीज़ है जो वास्तव में तुम्हारे computer या server पर चलती है।

यह आमतौर पर bash होता है, लेकिन तुम यहां शीर्ष पर एक अलग hash bang भी डाल सकते हो, और यह Python script या R script हो सकती है। कोई फर्क नहीं पड़ता। इस script में जो कुछ भी है वह execute होगा।

एक और चीज़ है जो हमने इस process में जोड़ी है, जो output declaration है। यह Nextflow को बताता है कि यह process output.txt नामक एक output फ़ाइल की उम्मीद कर रही है। यह कहता है कि यह एक path है, तो इसे एक फ़ाइल की तरह handle किया जाना चाहिए, न कि कहें, अगर यह val था, तो यह कहेगा कि यह एक variable या value की तरह है।

ध्यान दें कि यह इस फ़ाइल को नहीं बना रहा है। यह वास्तव में इसे generate नहीं कर रहा है। वह यहां नीचे script द्वारा किया जाता है। यह सिर्फ Nextflow को इस filename के साथ एक output फ़ाइल की उम्मीद करने के लिए कह रहा है।

## 1.1.2. वर्कफ़्लो परिभाषा

ठीक है। और फिर नीचे हमारे पास यहां एक workflow है, और फिर से, हमारे पास एक declaration है। इसे Main कहा जाता है। यह workflow का script block के बराबर है, अगर तुम चाहो। यह workflow का वह हिस्सा है जो कुछ करता है। और इस मामले में, हम कह रहे हैं, sayHello नामक process को call करो।

सामान्य रूप से, निश्चित रूप से, तुम्हारी pipeline इससे बहुत अधिक जटिल दिखेगी। तुम्हारे पास शायद एक से अधिक process होंगी, और तुम उनके बीच data flow को orchestrate करने के लिए channels का उपयोग करोगे। हम इस कोर्स के अगले भागों में उस पर आने वाले हैं, लेकिन अभी के लिए, यह काफी है। यह एक valid pipeline है, जो काम करनी चाहिए।

मैं VS Code में यहां preview DAG भी click कर सकता हूं। DAG या DAG pipeline में data flow structure का एक representation है, और हम इसे side पर mermaid diagram के रूप में rendered देख सकते हैं। इस मामले में यह बहुत सरल है। एक box है, जो workflow है और एक process है, जिसे sayHello कहा जाता है, लेकिन जैसे-जैसे हम आगे बढ़ते हैं यह अधिक दिलचस्प लग सकता है।

## 1.2. वर्कफ़्लो चलाएं

ठीक है, चलो इस workflow को चलाने की कोशिश करते हैं और देखते हैं कि क्या होता है।

मैं नीचे फिर से terminal लाने जा रहा हूं, output को clear करूंगा, और मैं Nextflow Run टाइप करने जा रहा हूं। और फिर मैं बस script का नाम टाइप करने जा रहा हूं, जो hello-world.nf है। और मैं enter दबाने जा रहा हूं।

ठीक है, शीर्ष पर कुछ standard चीजें मिली हैं, जो हमें बताती हैं कि Nextflow चला और कौन सा version चल रहा था और script का नाम और सब कुछ क्या था।

और वास्तव में महत्वपूर्ण चीज़ जिसे हम यहां ढूंढ रहे हैं वह _यहां_ है, जो अलग-अलग tasks का सारांश है जो execute किए गए थे।

अगर तुम्हारा इस तरह दिखता है एक छोटे हरे tick के साथ, तो बहुत बढ़िया। तुमने अभी-अभी अपनी पहली pipeline चलाई है। शानदार।

यह हमें यहां process का नाम बताता है, जो चला, जिसे Say Hello कहा जाता था, और इसने हमें बताया कि यह एक बार चला और यह सफल रहा। यह जैसे-जैसे तुम आगे बढ़ते हो update होता है, तो जब तुम एक बड़ी pipeline चला रहे होते हो, तो तुम यहां progress को represented देखोगे। लेकिन क्योंकि यह इतना छोटा है, यह मूल रूप से तुरंत चलता है।

## 1.2.2. work directory में output और logs खोजें

अब जब तुम एक Nextflow pipeline चलाते हो, तो उन processes में से प्रत्येक को एक साथ जोड़ा जाता है, और प्रत्येक process, जैसा कि मैंने पहले कहा, tasks उत्पन्न कर सकती है एक या कई। तो इस मामले में, हमारे पास इस process से एक single task था। यह बस एक बार चला और वह इस task _hash_ के तहत किया गया था।

Nextflow सीधे तुम्हारी working directory में फ़ाइलों से नहीं निपटता है, यह work नामक एक विशेष फ़ोल्डर बनाता है। और अगर मैं "ls" करता हूं, तो हम देखेंगे कि यह यहां दिखाई दिया है: _work_, और इसके भीतर प्रत्येक single task के लिए sub directories हैं जो चलती हैं। और वह इस hash से मेल खाता है। तो तुम देख सकते हो अगर मैं "ls work/c4" पर जाता हूं, और फिर यह truncated है, लेकिन यह 203 से शुरू होता है, और वह working directory है, जो इस process द्वारा बनाई गई थी जब हमने pipeline चलाई थी। और तुम इसे side पर भी देख सकते हो।

जब मैं उन फ़ाइलों को list करता हूं, तो तुम देख सकते हो कि output.txt फ़ाइल generate हुई थी। तुम इसे यहां भी देख सकते हो। और कुछ hidden फ़ाइलें हैं, जो मेरे regular "ls" के साथ नहीं दिख रही हैं।

अगर मैं output.txt पर click करता हूं, तो निश्चित रूप से, हमारे पास हमारा output है। शानदार। तो pipeline काम कर गई।

यह एक line bash script चलाने के लिए काफी boilerplate लग सकता है, लेकिन जैसे-जैसे हमारी processes अधिक जटिल होती जाएंगी, यह अधिक समझ में आएगा। और Nextflow के साथ यह work directory और ये फ़ाइलें, जो बनाई जाती हैं, वास्तव में Nextflow को इतना शक्तिशाली बनाने की रीढ़ हैं।

प्रत्येक task, pipeline का प्रत्येक element हर दूसरे task से isolated है। यह reproducible है। वे एक दूसरे के साथ conflict नहीं करते हैं, और सब कुछ parallel में चल सकता है। यह वास्तव में एक अच्छा तरीका है जब तुम इसके आदी हो जाते हो क्योंकि इस isolation के कारण तुम अंदर जा सकते हो और देख सकते हो कि एक single task के लिए वास्तव में क्या हुआ और debug कर सकते हो।

चलो work directory में इन अन्य फ़ाइलों पर एक त्वरित नज़र डालते हैं। ऊपर से नीचे तक, हमारे पास _.command.begin_ नामक एक फ़ाइल है। यह खाली है। यह बस एक sentinel फ़ाइल है, जिसे Nextflow द्वारा बनाया गया है जो कहती है, ठीक है, मैं task शुरू कर रहा हूं। वहां कुछ भी दिलचस्प नहीं है।

फिर _.command.error_, _.command.log_ और _.command.out_ हैं। ये सभी bash command या इस script से outputs हैं जो चली। यह standard error है। यह standard out है, और यह दोनों का combined है जैसे वे बाहर आए। तो तुम्हें logical order मिलता है।

ठीक है, वे सभी इसके लिए भी खाली थे, तो बहुत दिलचस्प नहीं, लेकिन चीजें अधिक दिलचस्प हो जाती हैं जब तुम _.command.run_ पर पहुंचते हो।

यह आमतौर पर एक बहुत लंबी script होती है। और यह वह है जो Nextflow वास्तव में execute करता है। अगर तुम यहां अंदर जाते हो, तो तुम Nextflow की सभी inner logic को देखना शुरू करोगे और देखोगे कि यह क्या कर रहा है और यह तुम्हारी process को कैसे execute कर रहा है। यह इस बात पर निर्भर करेगा कि तुम कहां चला रहे हो, चाहे हम locally चला रहे हों या इसे SLURM को एक job के रूप में submit कर रहे हों, जिस स्थिति में हमारे पास शीर्ष पर SLURM headers होंगे। ये सभी अलग-अलग setups।

आम तौर पर, तुम्हें वास्तव में कभी भी इस फ़ाइल को देखने की आवश्यकता नहीं है। यह Nextflow द्वारा autogenerated है और इसमें तुम्हारी pipeline के लिए विशेष रूप से कुछ भी unique नहीं है। लेकिन वह वास्तव में वह core है जो चल रहा है।

अगला बहुत अधिक दिलचस्प है। _.command.sh_ generated script है, जो तुम्हारी process से आई, और यहां तुम देख सकते हो कि Nextflow ने Bash header जोड़ा, और फिर इसने हमारा command execute किया, जो हमारे script block में था।

और बस इतना ही _.command.run_ फ़ाइल करती है कि यह बस इस _.command.sh_ फ़ाइल को चलाती है।

यह वास्तव में एक उपयोगी है, जो आमतौर पर सबसे अधिक देखा जाता है जब तुम किसी चीज़ को debug करने की कोशिश कर रहे होते हो और यह check कर रहे होते हो कि तुम्हारी Nextflow pipeline की logic वह कर रही है जो तुम उम्मीद करते हो।

अंत में, हमारे पास _.exitcode_ नामक एक फ़ाइल है, और यह बस एक task से exit code को capture करती है, जो इस मामले में सफल था। तो exit code zero था।

अगर कुछ गलत हो जाता है, तुम memory से बाहर हो जाते हो या कुछ और और यह fail हो जाता है, तो यह समझने के लिए बहुत उपयोगी है कि क्या गलत हुआ।

## 1.3. वर्कफ़्लो को फिर से चलाएं

work directories के बारे में समझने के लिए एक और चीज़ यह है कि अगर मैं इस pipeline को बार-बार चलाता रहता हूं, तो अगर मैं _"nextflow run hello-world.nf"_ करता हूं, तो यह बिल्कुल वही काम करने जा रहा है, लेकिन इस बार इसका एक नया task id होगा। तुम देख सकते हो कि यहां यह hash अलग है, और अब अगर मैं work में देखता हूं, तो दो hash directories हैं। और ये फिर से एक दूसरे से अलग हैं।

तो हर बार जब तुम एक Nextflow workflow चलाते हो, जब तक कि तुम resume का उपयोग नहीं करते जो cache का उपयोग करता है, हम बाद में touch करेंगे, यह उन processes को नई work directories में फिर से चलाने जा रहा है, जो एक दूसरे से अलग हैं। तुम्हें कोई file name collisions नहीं मिलेंगे, तुम्हें उस तरह की कोई समस्या नहीं होगी। सब कुछ isolated और clean है।

और अगर हम इस directory में जाते हैं, तो तुम सभी समान फ़ाइलें और समान _output.txt_ देख सकते हो, जो scratch से फिर से बनाई गई है।

## 2. outputs प्रकाशित करें

ठीक है, यह Nextflow के लिए अपने आप में बहुत अच्छा है, जबकि यह तुम्हारी pipeline चला रहा है ताकि सभी चीजें एक दूसरे से अलग हों और clean हों और manage की जा सकें।

लेकिन यह बहुत उपयोगी नहीं है अगर तुम एक व्यक्ति हो जो अपने results को explore करने की कोशिश कर रहे हो। तुम वास्तव में हजारों और हजारों अलग-अलग work directories के माध्यम से अपनी result फ़ाइलों को खोजने की कोशिश नहीं करना चाहते हो। और तुम वास्तव में ऐसा करने के लिए नहीं हो। work directories का मतलब यह नहीं है कि वे तुम्हारी फ़ाइलों के बनने की अंतिम स्थिति हों।

हम अपनी फ़ाइलों को publish करके ऐसा करते हैं।

## 2.1.1. sayHello process के output को declare करें

तो अगर मैं अपनी script पर वापस जाता हूं, तो हम अपने workflow block में काम करने जा रहे हैं। हम इसे बताने जा रहे हैं कि किन फ़ाइलों की उम्मीद करनी है, कौन सी फ़ाइलें हमें परवाह करती हैं, और फिर हम नीचे output block नामक एक नया block बनाने जा रहे हैं।

यह नया syntax है, जो Nextflow के version 26.04 में syntax parser के साथ आया और default है। तो अगर तुमने पहले थोड़ा Nextflow का उपयोग किया है, तो यह उन चीजों में से एक है जो नई है।

तो हमें main block मिला है, और अगला मैं publish कहने जा रहा हूं और मैं Nextflow को बताने जा रहा हूं कि publishing से क्या उम्मीद करनी है। हम इसे _first_output_ कहने जा रहे हैं, और हम इसे _sayHello.out_ कहने जा रहे हैं।

मैंने गलती से वहां एक typo बनाया, लेकिन यह Nextflow VS Code extension की कुछ features को भी इंगित करने का एक अच्छा अवसर है। तुम देख सकते हो कि तुरंत इसने मुझे इसके नीचे एक छोटी wiggly red line दी जो कहती है कि कुछ गलत है। और अगर मैं इस पर hover करता हूं, तो यह मुझे बताने जा रहा है कि यह variable defined नहीं है। मुझे नहीं पता कि यह क्या है।

इस मामले में यह बहुत स्पष्ट है, मैंने एक typo बनाया। मेरा मतलब sayHello टाइप करना था, और फिर wiggly line चली जाती है।

अब यह purple है। Nextflow syntax parser जानता है कि यह एक process है और जब मैं इस पर hover करता हूं, तो यह मुझे एक reduced representation देता है कि यह process कैसी दिखती है। तो मैं बहुत जल्दी एक नज़र में देख सकता हूं कि यह कोई inputs नहीं लेती है और यह हमें यह output देती है। तो VS Code में इस extension के साथ काम करना तुम्हें code लिखते समय बहुत सारी contextual information देता है।

ध्यान दें कि हम इस process से output को _.out_ syntax के साथ refer कर सकते हैं। और इस समय हम इसे जो चाहें कह सकते हैं, यह सिर्फ एक arbitrary variable name है।

## 2.1.2. script में एक output: block जोड़ें

जहां यह महत्वपूर्ण हो जाता है जब हम यहां अपना नया block करते हैं, और यह अब workflow block के नीचे है, हम अब workflow के अंदर नहीं हैं। फिर से squiggly brackets। और यह वह जगह है जहां हम बस Nextflow को बताते हैं कि सभी फ़ाइलों को कहां रखना है, जो workflow द्वारा बनाई गई हैं।

अब मैं इस variable name को लेने जा रहा हूं, जो मैंने यहां बनाया था, और मैं इसे वहां रखने जा रहा हूं और इसके लिए कुछ squiggly brackets डालूंगा। और मैं Nextflow को एक path का उपयोग करने के लिए कहने जा रहा हूं। उफ़। Path, quote marks में। और मैं dot का उपयोग करने जा रहा हूं। यह बस Nextflow को results directory के root में फ़ाइल डालने के लिए कहता है। तो कोई sub directories या कुछ भी नहीं।

चलो अपनी workflow को फिर से चलाने की कोशिश करते हैं। अगर मैं _"nextflow run hello-world.nf"_ करता हूं, तो उम्मीद है कि यह मूल रूप से बिल्कुल वैसा ही दिखना चाहिए। यहां Nextflow के साथ वास्तव में कुछ भी नहीं बदला है। यह वही चीजें चला रहा है। यह बस उन्हें फिर से work directories में कर रहा है।

लेकिन अब अगर मैं _"ls results/"_ करता हूं, तो तुम देखोगे कि यहां results नामक एक नई directory बनाई गई है जो workflow publishing के लिए default base directory है। और उसमें _output.txt_ नामक एक फ़ाइल है।

अगर मैं _"ls -l results"_ करता हूं, तो तुम देखोगे कि यह वास्तव में work directory से soft linked है। तो यह एक real फ़ाइल नहीं है, यह work directory से linked है और इसने हमारे लिए वहां सभी फ़ाइलों को collect किया है।

## 2.2. एक custom location सेट करें

"Results" इस path के लिए default name है। अगर मैं workflow को फिर से चलाता हूं, और इस बार मैं _dash_ single hyphen करता हूं, यह है, क्योंकि यह एक core Nextflow option है। _"-Output-dir **my** results"_। short के लिए _"-o"_ भी कर सकता था। फिर यह एक अलग base directory सेट करने जा रहा है जहां फ़ाइलें stored हैं और एक बार फिर, यहां _myresults/_ में, अब हमारे पास एक _output.txt_ है।

यह बहुत अच्छा है, लेकिन हम शायद नहीं चाहते कि सभी फ़ाइलें बस root में हों। हम कुछ organization चाहते हैं, तो हम यहां एक subdirectory भी बना सकते हैं जिसे हम जो चाहें कह सकते हैं। चलो कहते हैं _"path 'hello_world'"_, और मैं बस इसे फिर से चलाता हूं। _"nextflow run hello-world.nf"_। यह results directory में एक subdirectory में जाना चाहिए और निश्चित रूप से, अब शीर्ष पर results के तहत हमारे पास _hello_world/_ है और हमारे पास _output.txt_ है।

ध्यान देने वाली महत्वपूर्ण बात, पुरानी _output.txt_ फ़ाइल अभी भी वहां है। results directory को wipe नहीं किया जाता है जब तुम ऐसा करते हो। बस नई फ़ाइलें वहां copy की जाती हैं। वे उन फ़ाइलों को overwrite कर देंगी जो पहले से वहां हैं अगर उनका same file name है, लेकिन वे पुरानी को clear नहीं करेंगी। तो तुम्हें थोड़ा सावधान रहने की आवश्यकता है जब तुम pipelines को फिर से चलाते हो। अगर तुम नहीं चाहते कि वे उन फ़ाइलों के ऊपर हों जो पहले से वहां हैं। सुनिश्चित करो कि तुम एक blank empty directory का उपयोग करते हो।

## 2.3. publish mode को copy पर सेट करें

ठीक है, मैंने उल्लेख किया कि ये फ़ाइलें soft links हैं, तो अगर मैं _"ls -l results/hello_world/"_ करता हूं, तो तुम देख सकते हो कि यह work directory से soft linking कर रहा है। यह आम तौर पर एक अच्छी बात है अगर तुम HPC जैसी किसी चीज़ पर काम कर रहे हो, और ये वास्तव में बहुत बड़ी फ़ाइलें हैं और तुम उन्हें duplicate नहीं करना चाहते, क्योंकि इसका मतलब है कि फ़ाइलें file system पर केवल एक बार stored हैं।

हालांकि, इसका मतलब यह है कि अगर तुम work directory को delete करते हो: अगर मैं _"rm -r work"_ करता हूं और उन सभी intermediate फ़ाइलों को clear करता हूं जो बनाई गई थीं। अब, अगर मैं इस फ़ाइल _"results/hello_world/"_ को पढ़ने की कोशिश करता हूं। यह एक soft link के रूप में एक फ़ाइल की ओर इशारा कर रहा होगा जो अब मौजूद नहीं है और data हमेशा के लिए चला गया है और irretrievable है, जो शायद बहुत अच्छा नहीं है।

तो आम तौर पर हम, मैं कहता हूं कि यह अच्छा practice है कि फ़ाइलों को soft linking के बजाय copy करें अगर तुम कर सकते हो, क्योंकि यह safer है। बस aware रहो कि यह दोगुना disk space का उपयोग करेगा जब तक कि तुम उन work directories को delete नहीं करते।

output block के साथ ऐसा करने के लिए, मैं यहां first output पर जाने जा रहा हूं। मैंने पहले path सेट किया और अब मैं mode सेट करने जा रहा हूं और तुम देख सकते हो जैसे मैं टाइप करता हूं, VS code extension, चीजें suggest कर रहा है यह जानता है कि यह यहां एक output directive है। और मैं copy कहने जा रहा हूं। मैं save hit करता हूं।

चलो workflow को फिर से चलाते हैं। यह फ़ाइलों को फिर से बनाने जा रहा है, नई work directory।

अब, अगर मैं _"ls -l results/hello_world/"_ पर जाता हूं तो तुम देख सकते हो कि यह एक real फ़ाइल है और यह अब soft link नहीं है, और Nextflow ने उसे copy किया। जानना अच्छा है। तो path और mode ऐसी चीजें हैं जो तुम खुद को काफी बार लिखते हुए पाओगे।

अब, निश्चित रूप से, यह बहुत सरल है। हम इसे अधिक complex और powerful बनाएंगे जैसे-जैसे हम आगे बढ़ते हैं, और तुम देखोगे कि इन चीजों को dynamic कैसे बनाया जाए और बहुत verbose नहीं।

## 2.4. process-level publishDir directives पर नोट

अब, मैंने कहा जैसे हमने इस पर शुरू किया, कि यह syntax का एक काफी नया रूप है। यह केवल Nextflow के नवीनतम versions में उपलब्ध है जैसा कि मैं इसे record करता हूं, और इसे Workflow Outputs कहा जाता है।

अगर तुम इसका उपयोग करते हो, तो यह बहुत अच्छा है। यह Nextflow के भीतर बहुत सारी अन्य cool features को unlock करता है, जैसे, Nextflow Lineage जो इन फ़ाइलों की heritage को track करने में मदद करता है जैसे वे बनाई जाती हैं, और जल्द ही 26.04 में default होगा। और भविष्य में बाद की तारीख में, यह तुम्हारे workflows लिखने का एकमात्र तरीका होगा।

हालांकि, जैसा कि हम अभी इस transition phase में हैं, तुम अच्छी तरह से wild में pipelines देख सकते हो, जो तुम publishDir नामक कुछ का उपयोग करते हो, जो इसे करने का पुराना तरीका है, और यह workflow और output level पर नहीं, बल्कि process level पर defined है।

और यह declaration मूल रूप से वही बात कहती है। यह कहती है, results नामक एक directory में results फ़ाइलों को publish करो, और copy mode का उपयोग करो। तो तुम देख सकते हो कि syntax बहुत समान है। लेकिन जब तुम अब नई pipelines लिख रहे हो, तो इस publishDir directive का उपयोग न करने की कोशिश करो, भले ही तुम इसे देखते हो, AI results में या documentation में या अन्य pipelines में, क्योंकि यह इसे करने का पुराना तरीका है।

2026 में हम सभी को workflow outputs का उपयोग करना चाहिए।

यह सब documented है, अगर तुम यह कर रहे हो और तुमने पहले Nextflow का उपयोग किया है, तो तुम यहां Nextflow docs पर जा सकते हो, nextflow.io/docs/। और अगर मैं tutorials तक scroll करता हूं, तो _Migrating to Workflow Outputs_ नामक एक tutorial है।

यह वास्तव में अच्छा है। यह सभी syntax के माध्यम से जाता है, यह पुराने syntax के बराबर कैसे है, हमने इसे क्यों बदला, और, एक timeline और सब कुछ है। और यह loads और loads of examples के साथ सभी अलग-अलग scenarios के माध्यम से जाता है। तो तुम आसानी से मौजूदा Nextflow code को नए syntax में convert कर सकते हो।

## 3.1. sayHello process को एक variable input की उम्मीद करने के लिए बदलें

ठीक है, तो हमारी simple script मिल गई है, जो एक process चला रही है, एक फ़ाइल बना रही है, Nextflow को बता रही है कि यह एक output है, और फिर हम Nextflow को बता रहे हैं कि उस फ़ाइल को कहां save करना है। यह एक अच्छी शुरुआत है।

लेकिन यह अधिक दिलचस्प होगा अगर यह सब hardcoded नहीं होता। तो अगला, चलो सोचते हैं कि Nextflow को कैसे बताया जाए कि यह process एक variable input ले सकती है, जो कुछ ऐसा है जिसे हम runtime पर control कर सकते हैं जब हम एक workflow launch करते हैं।

हमें ऐसा करने के लिए कुछ अलग-अलग चीजें करने की आवश्यकता है।

सबसे पहले, हमें इस process को बताने की आवश्यकता है कि यह एक input variable स्वीकार कर सकती है और हम यहां एक नए declaration block के रूप में _input_ टाइप करते हैं। और हम इसे _"val greeting"_ कहने जा रहे हैं।

val bit यहां नीचे path के बराबर है। यह Nextflow को बताता है कि यह एक variable है, इस मामले में एक string की तरह। और अगर तुम इस पर फिर से hover करते हो, तो यह तुम्हें extension से बताता है कि इसका क्या मतलब है।

अगला हम Nextflow को बताने जा रहे हैं कि इसके साथ क्या करना है। यह सिर्फ यह कहना काफी नहीं है कि एक variable है। तुम्हें script में कहना होगा कि उस variable का उपयोग कैसे करना है। और इसलिए मैं यहां इस hardcoded string से छुटकारा पाने जा रहा हूं, और मैं एक variable डालने जा रहा हूं।

मैं जल्दी से इसे squiggly brackets के बिना करने जा रहा हूं बस तुम्हें दिखाने के लिए कि यह, allowed है, और यह इसे करने का पुराना, style तरीका है। लेकिन अब नए syntax के साथ, हम वास्तव में इसे इस तरह squiggly brackets के अंदर डालने की सिफारिश करते हैं, और यह बहुत स्पष्ट करता है कि यह यहां Nextflow द्वारा interpolated किया जा रहा है।

बढ़िया। तो _"input greeting"_ _$\{greeting\}_ में जाता है। आखिरी चीज़ यह है कि हमें workflow level पर Nextflow को बताने की आवश्यकता है कि यह process अब एक input लेती है। और ऐसा करने के लिए, हम मूल रूप से इसे एक variable देने जा रहे हैं।

## 3.2. user input को capture करने के लिए एक command-line parameter सेट करें

हम इसे फिर से hard code कर सकते थे, जैसे Hello World, और वह ठीक काम करेगा, लेकिन स्पष्ट रूप से यह वास्तव में हमें कोई लाभ नहीं देता है। हम इसे run time पर configure करने में सक्षम होना चाहते थे, तो हम इसे CLI पर करने में सक्षम होना चाहते हैं, जब तुम Nextflow launch करते हो।

और जिस तरह से हम ऐसा करते हैं वह _params_ नामक एक विशेष Nextflow concept है। हम इसे _params.input_ कहने जा रहे हैं।

यह क्या करता है यह इस input variable को CLI पर expose करता है और वह वह जगह है जहां हम double dash का उपयोग करते हैं जब हम Nextflow launch करते हैं।

मैं इसे जो चाहूं कह सकता हूं, मैं इसे _hello, greeting_ कह सकता हूं। कोई फर्क नहीं पड़ता। मैं वहां जो भी करूंगा वह एक CLI option के रूप में exposed होगा जब हम एक pipeline launch करते हैं। और यह Nextflow द्वारा एक real magic trick है क्योंकि इसका मतलब है कि तुम इन parameters के साथ अपनी workflow script बहुत जल्दी बना सकते हो, और तुम अनिवार्य रूप से अपनी pipeline के लिए एक custom CLI बना रहे हो, जब तुम launch करते हो तो विभिन्न options को on the fly customize करना वास्तव में आसान बनाते हो।

तो। चलो इसे आज़माते हैं। अपने terminal पर वापस जाओ। हमारे पास यहां हमारा _"nextflow run"_ command है। और अब मैं _"--input"_ करने जा रहा हूं, जो _"params.input"_ से मेल खाता है जो हमने पहले देखा था। मुझे लगता है कि docs में यह French में है। Geraldine को French बोलना पसंद है। मैं इसे Swedish में करने जा रहा हूं क्योंकि मैं Sweden में रहता हूं। तो मैं कहने जा रहा हूं, "_Hej Världen_" और enter hit करूंगा।

Single quotes या double quotes का उपयोग कर सकते हो, यह सिर्फ प्रभावित करता है कि Bash इसे कैसे interpret करता है।

यह Nextflow pipeline को बिल्कुल उसी तरह चलाता है। तुम working directory और सब कुछ देख सकते हो वही है। लेकिन अब अगर मैं _"results/hello_world/output"_ पर जाता हूं। हम यहां अपनी अच्छी Swedish देख सकते हैं।

तो हमने dynamically रूप से एक CLI से एक parameter को एक input pass किया है। हमने उसे process को एक input के रूप में pass किया है और process ने उसे interpret किया और इसे एक script block में डाला, जिसने फिर dynamically रूप से उस script result के output को बदल दिया। बहुत cool।

यहां बहुत कम syntax के साथ काफी complex logic। और तुम उम्मीद कर सकते हो कि यह अब कैसे scale करना शुरू होता है। और यह वह तरीका है जिससे हम वास्तव में अपनी pipelines की logic और customizability को Nextflow script में बनाते हैं।

## 3.4. command line parameters के लिए default values का उपयोग करें

ठीक है, यह बहुत अच्छा है। हालांकि अब समस्या यह है, हर बार जब मैं इस pipeline को चलाता हूं, तो मुझे इसे चलाने के लिए dash, input करने की आवश्यकता होती है।

अगर मैं इस parameter के बिना चलाने की कोशिश करता हूं, तो अब Nextflow एक error throw करने जा रहा है जो कहता है कि इसे इस parameter की आवश्यकता थी और यह set नहीं था। और इसलिए यह नहीं जानता था कि क्या करना है।

वैसे, यह एक cool नई चीज़ है। अतीत में, Nextflow बस एक empty string के साथ चला होता, और तुम्हें सभी प्रकार की अजीब errors मिली होतीं, जिन्हें समझना मुश्किल होता। लेकिन नए Nextflow syntax parser में, यह थोड़ा अधिक सावधान है और यह तुम्हें तुरंत बताता है।

तो हम हमेशा हर single option को specify नहीं करना चाहते। sensible defaults को specify करना अच्छा practice है। तो हम अपनी script में ऐसा कैसे करते हैं?

तुम देखोगे कि जब हमने यह लिखा, तो हमने बस _params.input_ को सीधे वहां डाल दिया जहां हम इसका उपयोग कर रहे हैं। तो स्पष्ट solution यह है कि हम एक default define करते हैं, और हम यहां workflow में script के शीर्ष पर एक विशेष params block में ऐसा करते हैं। यह यहां workflow script में है।

फिर से, यहां कुछ नया syntax है, तो ध्यान दो। यह वास्तव में cool stuff है। हमारे पास parameter का नाम है, जिसकी यहां उम्मीद की जाएगी।

और फिर इस colon character के बाद, हम variable के type को define कर रहे हैं। तुम्हें ऐसा करने की आवश्यकता नहीं है, तुम इसे बस blank छोड़ सकते हो, लेकिन यह वास्तव में अच्छा है। यह Nextflow को बताता है कि हम एक string की उम्मीद कर रहे हैं और इसे वैसे ही treat करें।

अगर हम इसके बजाय एक number चाहते हैं, उदाहरण के लिए, हम float लिख सकते थे, और वह कहेगा कि हम एक floating point number चाहते हैं। और अगर हम उसके साथ चलाने की कोशिश करते हैं, तो यह एक error throw करेगा। अगर हम इसे एक string देते हैं, जो float नहीं है। और यह इसे वैसे ही pass भी करेगा। जैसे अगर हम string करते हैं, तो यह जानता है कि यह एक string है। और भले ही इसमें leading zeros हों और यह सभी numeric हो, यह अभी भी इसे एक actual string के रूप में pass करेगा।

तो वह type safety Nextflow की एक बहुत नई feature है, लेकिन तुम्हारे code को लिखने और चलाने के लिए safer बनाने के लिए वास्तव में शक्तिशाली है।

फिर उसके बाद हमारे पास एक equal symbol है और फिर यहां default value है। Nextflow मूल रूप से Barcelona में लिखा गया था, तो यह उपयुक्त लगता है कि हमारे पास यहां कुछ, Spanish है, _"Holà mundo!"_ एक default के रूप में।

ठीक है मैं उस script को save करने जा रहा हूं, वापस जाओ, script को फिर से _--input_ के बिना चलाओ। और इस बार यह चलना चाहिए और यह _results_ में हमारी नई फ़ाइल बनाएगा। और इस फ़ाइल में अब यह _"Holà mundo!"_ कहता है।

हालांकि यह सिर्फ एक default है, तो इसका मतलब यह नहीं है कि हम अभी भी पहले जैसा ही नहीं कर सकते। अगर मैं वापस जाता हूं और यहां अपनी पुरानी script ढूंढता हूं, _"Hej Världen"_, क्योंकि मैं command line पर _--input_ करता हूं, वह उस default को overwrite करेगा और output.txt फ़ाइल में फिर से उसका उपयोग करेगा।

तो script में यह केवल default value है जो मैं set कर रहा हूं।

जैसे-जैसे हम अपनी workflow को अधिक complex बनाते हैं और अधिक parameters शामिल करते हैं, script के शीर्ष पर यह params block उन सभी को एक जगह collect करना शुरू कर देगा।

और तुम अपनी script में इस काफी अच्छी symmetry के साथ समाप्त होते हो, जहां तुम्हारे पास effectively यहां तुम्हारे सभी workflow inputs हैं और नीचे तुम्हारे workflow outputs हैं। और यह बहुत स्पष्ट है कि बाहरी दुनिया के लिए तुम्हारी workflow का interface क्या है। तो तुम नए syntax के साथ बहुत जल्दी एक नई pipeline उठा सकते हो और समझ सकते हो कि इसका उपयोग कैसे करना है।

एक आखिरी cool चीज़। हमें इसके साथ एक default value set करने की आवश्यकता नहीं है। अगर हम params input करते हैं लेकिन default value set नहीं करते हैं, तो यह Nextflow को बताता है कि यह parameter required है, और फिर से, pipeline इसके बिना चलने में fail हो जाएगी, लेकिन यह तुम्हें इसके null होने के बारे में कुछ के बजाय एक अधिक उपयोगी error message देगा।

तो यह कहता है कि हम उम्मीद कर रहे हैं कि इसका input required है, लेकिन यह command line पर specified नहीं था। बहुत अच्छा।

ठीक है, तो उम्मीद है कि अब यह स्पष्ट है कि variable inputs और parameters के साथ अपनी Nextflow pipeline को कैसे set up करना है, default कैसे set करना है, types कैसे set करना है, यह एक Boolean true false flag या integer या यहां अलग-अलग types हो सकता है। उन्हें अपनी workflow में कैसे pass करना है, यह कहां से गुजरता है, और फिर अपनी process में interpolates करता है। और फिर तुम यह भी जानते हो कि जब तुम Nextflow launch करते हो तो command line पर उन्हें कैसे customize करना है। यह हमारे simple bash command से अधिक दिलचस्प दिखना शुरू हो रहा है।

## 4. workflow executions को manage करें

ठीक है। आगे क्या है? इस अध्याय के अंतिम भाग के लिए, हम सभी अलग-अलग workflow executions को manage करने के बारे में थोड़ी बात करने जा रहे हैं। अगर तुम मेरे sidebar में और work के नीचे Explorer में देखते हो, तो तुम देखोगे कि मैंने कुछ अलग-अलग pipelines चलाई हैं और ये work directories काफी लंबी हो रही हैं, इनमें से बहुत सारी हैं।

और दूसरी बात यह है, जैसा कि मैंने पहले कहा, हर बार जब मैं इस pipeline को फिर से चलाता हूं, तो यह work directories का एक नया set बना रहा है, और यह सभी processes को scratch से फिर से चला रहा है, जो एक अच्छी बात है। यह intended behavior है। यह reproducible है और यह सब कुछ fresh regenerate कर रहा है। लेकिन यह स्पष्ट रूप से, अगर तुम बहुत लंबे समय तक चलने वाली processes चला रहे हो, तो हमेशा अपनी pipeline को शुरुआत से शुरू करना annoying है अगर यह आधे रास्ते में crash हो गई, या अगर तुम pipeline के अंत में कुछ बदलते हो।

## 4.1. -resume के साथ एक workflow को फिर से launch करें

सौभाग्य से, Nextflow वास्तव में, यह जानने में अच्छा है कि पहले क्या चलाया गया है और क्या उपलब्ध है, और उन पुराने results को reuse करना बहुत, सरल है। हम बस command के अंत में एक, नया flag जोड़ते हैं _"-resume"_।

अब, ध्यान दें कि input पर दो hyphens हैं क्योंकि वह parameter है। resume पर केवल एक hyphen है क्योंकि वह एक core Nextflow option है।

यह लोगों को हर समय trip करता है, भले ही तुमने लंबे समय से Nextflow का उपयोग किया हो। तो हमेशा याद रखो एक या दो hyphens। निर्भर करता है कि यह एक core Nextflow option है या नहीं।

ठीक है, तो अब मैं _-resume_ करता हूं और मैं बिल्कुल वही workflow फिर से चलाता हूं। और इस बार यह एक key difference के साथ बहुत ज्यादा बिल्कुल वैसा ही दिखना चाहिए।

यहां output में, तुम देख सकते हो कि results cached थे। और वास्तव में, यहां यह task hash पिछले run के बिल्कुल समान है, और इसने बस उस work directory को इसकी संपूर्णता में reuse किया। inputs और outputs और script सभी unmodified थे। और इसलिए यह बस उससे उस फ़ाइल को लेता है और अगर pipeline में downstream steps हैं, तो यह उन्हें pipeline में अगले step पर pass करेगा।

तो यह अभी भी शुरू से अंत तक पूरी pipeline चला रहा है, लेकिन यह उन tasks में से प्रत्येक के लिए cached results का उपयोग कर रहा है, जहां यह कर सकता है।

अब, जब तुम _-resume_ करते हो, तो यह बस तुम्हारी working directory में आखिरी pipeline run को resume करता है, जो भी वह था। लेकिन तुम वास्तव में किसी भी पिछले run से resume कर सकते हो जो तुमने वहां किया है। और हमने अब काफी बहुत कुछ किया है।

## 4.2. past executions के log का निरीक्षण करें

उन सभी को देखने के लिए, हम _"nextflow run"_ के बजाय _"nextflow log"_ कर सकते हैं, और वह हमें इन सभी अलग-अलग का एक अच्छा output देगा.. मुझे अपनी screen को थोड़ा छोटा करने की आवश्यकता है ताकि हम इसे देख सकें, ये सभी अलग-अलग runs जब हमने उन्हें किया, session id, command और सब कुछ।

और हम यहां देख सकते हैं और हम इनमें से किसी का भी run name ले सकते हैं और फिर उन specific में से एक को resume कर सकते हैं। तो मैं वापस जा सकता हूं और मैं उस _hungry_ekeblad_ नामक को resume कर सकता हूं। और मैं बस उसे _resume_ के बाद डालता हूं।

अगर तुम curious हो, वैसे, ये सभी adjectives और scientists names Nextflow source code में हैं। यह Nextflow को अपना पहला ever pull request प्राप्त करने का एक वास्तव में अच्छा तरीका है इसे जाकर और ढूंढकर और अपने पसंदीदा scientists को जोड़कर।

और वैसे भी, तो मैंने वह किया और यह वापस गया और इसने इस workflow run से cached results को देखा, realize किया कि यह अभी भी उनका reuse कर सकता है, और इसने किया। तो मुझे फिर से cached results मिले।

## 4.3. पुरानी work directories को delete करें

यह बहुत अच्छा है। क्या होगा अगर मैं इन work directories को clean up करना चाहता हूं? यहां इनमें से loads हैं। loads of फ़ाइलें हैं। शायद मुझे fact के लिए पता है कि मैं आखिरी कुछ pipeline runs से resume करना चाहता हूं, लेकिन मुझे उससे पहले की सभी की परवाह नहीं है।

फिर मैं यहां एक pick कर सकता हूं और मैं एक और Nextflow command का उपयोग कर सकता हूं, जो _"nextflow clean"_ है, और मैं _"nextflow clean"_ कर सकता हूं, मैं _"-before"_ करने जा रहा हूं, और particular run name, जो इस मामले में _reverent_pike_ था और मैं _"-n"_ करने जा रहा हूं, जो Nextflow को बस एक dry run करने के लिए कहता है। तो यह बस मुझे बताता है कि यह क्या delete करेगा। वास्तव में कुछ भी किए बिना, तो यह इन work directories को remove करेगा।

वह sensible लगता है। तो मैं वही command फिर से करने जा रहा हूं, लेकिन _"-n"_ के बजाय मैं वास्तव में cleanup करने के लिए _"-f"_ करूंगा। और इस बार इसने वास्तव में इन सभी directories को removed किया है। और अगर मैं अंदर जाता हूं और work directories को देखता हूं, तो यह अब बहुत lighter दिख रहा है। शानदार।

तो यह है कि अपनी सभी local work directories को एक pretty safe तरीके से कैसे clean up करना है बिना cache को पूरी तरह से destroy किए। तो तुम अभी भी resume कर सकते हो अगर तुम चाहो।

अगर कभी तुम भूल जाते हो कि हर Nextflow command के लिए ये flags क्या हैं तो तुम _"nextflow help"_ कर सकते हो, और फिर command का नाम। तो अगर मैं _"nextflow help clean"_ करता हूं, तो तुम सभी अलग-अलग options देख सकते हो: _-after, -before, -but_, इस cleanup behavior को configure करने के सभी अलग-अलग तरीके। बहुत cool।

## सारांश

ठीक है, यह Hello Nextflow के भाग एक का अंत है। यह कोर्स की काफी intense शुरुआत है, लेकिन उम्मीद है कि अब तुम्हें एक Nextflow script कैसी दिखती है इसकी काफी अच्छी समझ है; अलग-अलग key parts के साथ, processes, workflows, outputs, और parameters। तुम जानते हो कि उन्हें command line से basic overrides के साथ कैसे configure करना है, dynamic input block को dynamic script के साथ कैसे बनाना है और तुम जानते हो कि अपने सभी workload executions को कैसे manage करना है: यह देखना कि तुमने पहले क्या चलाया है, resuming, cleaning up। बहुत सारी चीजें हैं। तुम एक लंबा रास्ता तय कर चुके हो। तो अगर तुम break लेना चाहते हो और एक quick walk around और एक cup of tea लेना चाहते हो, तो अब शायद एक अच्छा समय है। तुमने इसे earn किया है।

यहां से, हम मूल रूप से इस foundation पर building कर रहे हैं। हम इसे अधिक complex, अधिक powerful कैसे बना सकते हैं? हम इसे अधिक flexible कैसे बना सकते हैं? वे चीजें करें जो हम scale पर अपने analysis करना चाहते हैं।

## क्विज़

अब अगर तुम webpage पर भाग एक, hello world तक scroll करते हो तो तुम एक छोटा quiz देखोगे और यह कुछ नया है जो हमने Nextflow training के इस version के लिए किया है। और तुम जा सकते हो और खुद को quiz कर सकते हो यह check करने के लिए कि तुमने इस अध्याय में हमने जो सभी material किया है उसे समझ लिया है।

यह हमें भेजा नहीं जाता है या कुछ भी, यह बस तुम्हारे browser में stored है। तो हम नहीं जानते कि तुम्हारे answers क्या हैं, लेकिन यह बस एक छोटा self check है यह सुनिश्चित करने के लिए कि तुमने कुछ भी miss नहीं किया है या कुछ भी misunderstood नहीं किया है। और तुम इसे जितनी बार चाहो उतनी बार try कर सकते हो।

अगर तुम मेरे जैसे हो, तो शायद तुम अपने VS Code instance में terminal में रहना चाहते हो, जिस स्थिति में तुम _quiz_ command टाइप कर सकते हो और फिर बस इसे बता सकते हो कि तुम किस, chapter पर हो। तो हम _"Hello World"_ करते हैं, और फिर तुम बिल्कुल वही, quiz questions कर सकते हो, जो web browser में हैं, लेकिन बस अपने terminal में।

Cool। ठीक है। उम्मीद है कि तुम्हें वह पसंद आया। थोड़ा मज़ा करो और, हम तुम्हें अगले अध्याय में बस एक मिनट में देखेंगे Nextflow channels के बारे में बात करने के लिए।
