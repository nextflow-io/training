# भाग 6: Hello Config - ट्रांसक्रिप्ट

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "महत्वपूर्ण नोट्स"

    यह पेज केवल ट्रांसक्रिप्ट दिखाता है। पूर्ण चरण-दर-चरण निर्देशों के लिए, [कोर्स सामग्री](../06_hello_config.md) पर वापस जाएं।

    ट्रांसक्रिप्ट में दिखाए गए सेक्शन नंबर केवल संकेत के उद्देश्य से प्रदान किए गए हैं और सामग्री में सभी सेक्शन नंबर शामिल नहीं हो सकते हैं।

## स्वागत

नमस्ते, Hello Nextflow प्रशिक्षण कोर्स के भाग छह में आपका स्वागत है।

इस अध्याय को Hello Config कहा जाता है, और यह हमारे प्रशिक्षण कोर्स का अंतिम भाग है।

इस अध्याय में, हम Nextflow configuration के बारे में बात करने जा रहे हैं। Nextflow configuration वास्तव में शक्तिशाली है। यह हमें एक ही pipeline को कई विभिन्न compute infrastructures पर विभिन्न software provisioning के साथ चलाने की अनुमति देता है। और pipeline में विभिन्न options के साथ।

इसका मतलब है कि आप अन्य लोगों द्वारा बनाए गए Nextflow pipelines ले सकते हैं और उन्हें अपने system पर चला सकते हैं, भले ही वे पूरी तरह से अलग infrastructure के लिए बनाए गए हों। Nextflow को configure करने की यह क्षमता workflows को वास्तव में portable और shareable बनाती है।

इस अध्याय में, हम पिछले भागों में बनाए गए workflow का उपयोग करेंगे, लेकिन हम workflow code को बिल्कुल भी edit नहीं करेंगे। हम सिर्फ अपनी Nextflow config फ़ाइल को देखेंगे और देखेंगे कि config बदलने से Nextflow के चलने का तरीका कैसे बदलता है।

ठीक है, चलिए शुरू करते हैं।

पहले की तरह, चलिए training.nextflow.io पर जाकर शुरू करते हैं। बाईं ओर Hello Nextflow और अध्याय छह पर जाएं। Hello config। अब मैं अपने GitHub code Spaces environment में जा रहा हूं और उस script को दोबारा जांचूंगा जिसका हम उपयोग करेंगे।

## 0. Warmup: जांचें कि Docker सक्षम है और Hello Config workflow चलाएं

इसे Hello Config कहा जाता है, और यह वहीं से शुरू हो रहा है जहां हम पहले थे। तो हमारे तीन parameters के साथ बिल्कुल वैसा ही दिख रहा है। CSV फ़ाइल के लिए Greetings, आउटपुट batch name के लिए batch और cowpy name के लिए character। हमारे पास विभिन्न processes के चार imports हैं, और फिर हमारे पास एक workflow है जहां हम उन्हें एक साथ chain करते हैं।

मैं वास्तव में अब इस फ़ाइल को बंद करने जा रहा हूं क्योंकि हम इस अध्याय में Nextflow फ़ाइल को बिल्कुल भी touch नहीं करेंगे। हम पूरी तरह से configuration फ़ाइल के भीतर काम करेंगे। अगर मैं Nextflow dot config फ़ाइल में देखूं जिसे हमने पिछले अध्याय पांच में संक्षेप में देखा था, तो हम देख सकते हैं कि यहां हमारे पास एक single statement है: Docker enabled equals true, जो Nextflow को बता रहा है कि इस workflow को execute करते समय Docker का उपयोग करें।

मैं यहां pipeline root में Nextflow dot config का उपयोग कर रहा हूं, जो मेरे Nextflow चलाने पर automatically load हो जाता है। लेकिन याद रखें, Nextflow कई जगहों से config फ़ाइलें load कर सकता है।

अगर मैं Nextflow docs के साथ जांच करूं और Configuration पर जाऊं, तो आप इन जगहों की एक list और एक priority देख सकते हैं जिसमें वे load होते हैं।

ठीक है। चलिए जांचते हैं कि हमारा workflow हमारी अपेक्षा के अनुसार execute हो रहा है। एक terminal खोलें। Nextflow करें। Run। Hello, config। और enter दबाएं। हमारे पास वे चार processes चल रही होनी चाहिए, एक cowpy कमांड के साथ समाप्त होनी चाहिए। निश्चित रूप से, यह ठीक से काम किया। Docker enabled था, Docker को pull किया और मेरे लिए cowpy चलाया, बिल्कुल वैसे ही जैसे अध्याय पांच के अंत में किया था।

## 1. निर्धारित करें कि किस software packaging technology का उपयोग करना है

ठीक है। मान लीजिए मैं एक HPC पर चल रहा हूं और मेरे पास Docker installed नहीं है। इस परिदृश्य में करने के लिए सबसे अच्छी बात Singularity या Apptainer का उपयोग करना होगी। अगर मैं ऐसा करने जा रहा था, तो मैं module cowpy में जाऊंगा और इस container को singularity image का उपयोग करने के लिए बदलूंगा जैसा कि मैंने पिछले अध्याय में दिखाया था, एक oras:// के साथ, जो आप Seqera Containers से भी प्राप्त कर सकते हैं।

फिर मैं Nextflow dot config में जाकर Docker enabled को false पर सेट करूंगा और singularity enabled equals true करूंगा। या, अगर Apptainer का उपयोग कर रहे हैं, तो apptainer enabled equals true और यह काम करेगा।

Nextflow, यह containers के अलावा अन्य technologies को भी support करता है, जिससे आप परिचित हो सकते हैं वह है conda। यहां हम conda enabled equals true कर सकते हैं और Docker को false पर सेट कर सकते हैं। conda वही container directive का उपयोग नहीं करता है। इसके बजाय, हम यहां एक नया जोड़ सकते हैं जिसे conda कहा जाता है। फिर हम उस conda package को specify करते हैं जिसका हम उपयोग करना चाहते हैं। pipeline को यथासंभव reproducible बनाने की कोशिश करने के लिए जितना संभव हो उतना specific होना अच्छा अभ्यास है। तो मैं conda channel, conda-forge, और फिर cowpy, और exact version को specify करने जा रहा हूं, जो 1.1.5 था।

मैं सिर्फ cowpy भी लिख सकता था अगर मैं चाहता, लेकिन यह pipeline के विभिन्न executions पर cowpy के एक अलग version में resolve हो सकता है।

इसके बारे में अच्छी बात यह है कि मैंने docker directive को बिल्कुल भी touch नहीं किया है। यह Docker image अभी भी वहां है। मैं अब सिर्फ दो alternatives प्रदान कर रहा हूं, और इन्हें केवल एक config फ़ाइल का उपयोग करके on या off किया जा सकता है।

## 1.3. यह verify करने के लिए workflow चलाएं कि यह Conda का उपयोग कर सकता है

Conda अब enabled है, तो चलिए इसे आज़माते हैं।

बढ़िया। यह चल रहा है और आप देख सकते हैं कि यहां Nextflow से एक message है जो कह रहा है कि Nextflow मेरे लिए एक conda environment बना रहा है, और यह इस cache location का उपयोग कर रहा है।

Background में, Nextflow मेरे लिए "conda create" कमांड चला रहा है ताकि सिर्फ उन packages के साथ एक नया isolated conda environment बना सके जो मुझे चाहिए, और फिर उन conda packages को install और fetch कर रहा है ताकि यह process चला सके।

आप देख सकते हैं कि वहां थोड़ा समय लगा क्योंकि यह पहली बार environment बना रहा था और software install कर रहा था। हालांकि, इसने इस environment को cache किया है, इसलिए अगर मैं फिर से वही Nextflow कमांड चलाता हूं, तो यह बहुत जल्दी होना चाहिए क्योंकि यह उसी conda environment का पुन: उपयोग करेगा।

इसके बारे में अच्छी बातों में से एक यह है कि इन directives को process level पर specify किया जा सकता है, न कि सिर्फ पूरे workflow के लिए। इसलिए यदि आप चाहें, तो आप विभिन्न processes के लिए किस technology का उपयोग किया जाता है, इसे mix और match कर सकते हैं।

## 2. Process directives के साथ compute resources आवंटित करें

Nextflow configuration फ़ाइल सिर्फ software packaging से कहीं अधिक कर सकती है। हम Nextflow को यह भी बता सकते हैं कि वास्तव में pipeline में steps को कैसे चलाना है। एक उदाहरण host system को बताना है कि प्रत्येक executing task के लिए कौन से resources उपलब्ध कराए जाने चाहिए।

By default, Nextflow बहुत ज्यादा नहीं देता है। यह प्रत्येक process को एक single CPU और केवल दो gigabytes की memory देता है।

यह शायद कुछ ऐसा है जिसे हम बदलना चाहेंगे, ताकि जो processes चलने में लंबा समय लेती हैं वे अधिक resources प्राप्त कर सकें और अधिक तेज़ी से चल सकें, लेकिन यह जानना मुश्किल हो सकता है कि एक process को क्या allocate करना है। Nextflow के पास इसमें आपकी मदद करने के लिए कुछ अच्छी tricks हैं।

## 2.1. Resource utilization report generate करने के लिए workflow चलाएं

चलिए फिर से workflow चलाते हैं। इस बार, मैं एक अतिरिक्त argument जोड़ने जा रहा हूं, जो है dash with reports। यह एक core Nextflow option है, इसलिए यह एक single hyphen है। और फिर जो भी फ़ाइल name मुझे पसंद है। इस मामले में, मैं इसे report config one html कहने जा रहा हूं।

मैं फिर से workflow चलाने जा रहा हूं। यह बिल्कुल पहले की तरह चलने वाला है, लेकिन यह मुझे एक अतिरिक्त helper report देगा, जो आप देख सकते हैं कि अब sidebar में pop up हो गया है।

मैं इस फ़ाइल पर right click करने जा रहा हूं, download पर click करूंगा, जो इसे GitHub code Spaces से मेरे local system में download करता है, ताकि मैं आसानी से इसे यहां web browser में देख सकूं।

यह report किसी भी Nextflow run के लिए generate किया जा सकता है, और इसमें बहुत सारी information होती है। यह top पर कुछ metadata के साथ शुरू होता है कि कौन सी कमांड उपयोग की गई थी, workflow कब चला, इसमें कितना समय लगा, लेकिन जैसे-जैसे आप नीचे scroll करते हैं, हमें resources के बारे में अधिक detailed information मिलती है, जो pipeline में प्रत्येक step द्वारा उपयोग किए गए थे।

क्योंकि प्रत्येक process विभिन्न tasks के लिए कई बार चलता है। हमारे पास एक box plot है जो प्रत्येक process के लिए उपयोग किए गए resources की variation दिखा रहा है।

अगर मैं थोड़ा और नीचे scroll करूं, तो मुझे memory used और job duration के बारे में समान information दिखाई देती है। साथ ही disk read write।

आप कल्पना कर सकते हैं कि long running tasks के साथ एक बड़े pipeline के लिए, यह resources के configuration को fine tune करने के बारे में बहुत informative हो सकता है जो आप request कर रहे हैं ताकि आप over request न करें, बल्कि इतना प्रदान कर सकें कि यह जल्दी चले।

अगर मैं report को नीचे scroll करता रहूं, तो हम एक task table भी देखते हैं, जो हमें workflow में चलाए गए हर एक task के बारे में detailed information दिखाता है। इसमें resolved script जैसी information शामिल है, जो चलाई गई थी।

ठीक है, चलिए अपनी config फ़ाइल पर वापस जाते हैं। हमने देखा कि हमें वास्तव में अपने workflow के लिए बहुत ज्यादा की आवश्यकता नहीं है, तो चलिए Nextflow को बताते हैं कि हमें workflow में प्रत्येक process के लिए केवल एक gigabyte की memory की आवश्यकता है।

अब जब हम इसे इस तरह process level पर define करते हैं, तो यह pipeline में प्रत्येक single process पर apply होता है।

## 2.3. एक individual process के लिए resource allocations सेट करें

तर्क के लिए, चलिए मान लेते हैं कि cowpy वास्तव में बहुत heavy lifting कर रहा है और इसे अन्य tasks की तुलना में अधिक resources की आवश्यकता है। हम यहां config का एक extra block define कर सकते हैं, जो सिर्फ उस process पर apply होता है, with name cowpy का उपयोग करके।

इसे config selector कहा जाता है, और हम यहां विभिन्न processes को match करने के लिए अलग-अलग patterns define कर सकते हैं। उदाहरण के लिए, मैं cow star कर सकता हूं। फिर मैं उसके बाद कुछ curly brackets डालता हूं और चलिए इसे एक के बजाय दो gigabytes की memory दें और चलिए कहें दो CPUs।

अब Nextflow workflow में प्रत्येक process को एक gigabyte दे रहा होगा इस request के अलावा, जो अधिक specific है। तो यह इसे override कर देता है। और सिर्फ किसी भी processes के लिए जिन्हें cowpy कहा जाता है, दो gigs की memory और दो CPUs मिलेंगे।

ध्यान दें कि Nextflow resource utilization के बारे में clever है। इसलिए यदि आप इन संख्याओं को उच्च values पर रखना शुरू करते हैं, तो आप देखेंगे कि Nextflow job submissions को एक के बाद एक queue करना शुरू कर देता है, बजाय उन सभी को parallel में चलाने के, ताकि यह उपलब्ध resources को over request न करे।

## 2.4. Modified configuration के साथ workflow चलाएं

चलिए फिर से workflow चलाने की कोशिश करते हैं और चलिए इस बार एक नई report save करें।

ठीक है, हम इस फ़ाइल को download कर सकते हैं और एक नज़र डाल सकते हैं।

हां, आश्चर्यजनक रूप से नहीं, यह मूल रूप से बिल्कुल वैसा ही दिखता है क्योंकि यह एक dummy workflow है, जो कुछ भी real नहीं कर रहा है। लेकिन आप कल्पना कर सकते हैं कि इस तरह के reporting के साथ limits को define करने और real life workflows करने का यह iterative approach आपको appropriate configuration सेट करने के लिए एक evidence-based approach करने की अनुमति कैसे देता है और वास्तव में आपके लिए उपलब्ध computational resources का अधिकतम लाभ उठाता है।

आप इसके बारे में वास्तव में clever होना शुरू कर सकते हैं। Nextflow में failures को retry करने की built-in क्षमता है, और आप अपनी config फ़ाइल में इस तरह एक closure का उपयोग करके और dynamically उपलब्ध resources को सेट करके इसका लाभ उठा सकते हैं। तो यहां मैंने Nextflow को उस दो gigabyte को retry attempt से multiply करने के लिए कहा है। तो दूसरी retry को चार gigs मिलेंगे, तीसरी retry को छह gigs मिलेंगे और इसी तरह। यह इस प्रशिक्षण कोर्स के दायरे से थोड़ा बाहर है, लेकिन अगर आप रुचि रखते हैं, तो Nextflow docs को check करें, जिसमें dynamic retry logic के बारे में एक अच्छा section है।

## 2.5. Resource limits जोड़ें

अब, एक चीज जो आपको इसके बारे में पता चल सकती है वह यह है कि इस तरह की चीज़ से गलती से आपके system पर उपलब्ध resources से परे जाना काफी आसान हो सकता है। यदि आप उपलब्ध resources से अधिक resources की request करते हैं तो Nextflow आपके configuration के बारे में एक error throw करेगा और run को halt कर देगा। इससे बचने के लिए, आप resource limits नामक किसी चीज़ का उपयोग कर सकते हैं।

हमारे workflow में process scope के तहत, हम इस तरह resource limits define कर सकते हैं, जो एक array लेता है, और हम इस system पर उपलब्ध maximum memory CPUs और time specify कर सकते हैं।

यहां high values सेट करना उन resources की मात्रा को नहीं बढ़ाता है जिनकी request की जाती है। हम अभी भी अपनी requests में एक gigabyte का उपयोग करने जा रहे हैं, लेकिन इसका मतलब है कि यदि इनमें से कोई भी request 750 तक पहुंच जाती है, तो वे उस ceiling से टकराएंगी और उससे अधिक कुछ भी request नहीं किया जाएगा, जिसका अर्थ है कि Nextflow चलना जारी रखेगा और unavailable resources के कारण crash नहीं करेगा।

तो यह उपयोग करने के लिए एक अच्छा safeguard है, खासकर यदि आप अपने resource allocation के साथ dynamic logic का उपयोग कर रहे हैं।

दूसरी स्थिति जहां यह वास्तव में उपयोगी है वह यह है कि यदि आप public pipelines का उपयोग कर रहे हैं जो आपके द्वारा नियंत्रित नहीं हैं। वे configuration defaults के साथ आ सकते हैं, और Nextflow automatically आपके system पर चलाने के लिए किसी भी resource requests को thresholding करने का सही approach लेगा।

ठीक है, बढ़िया। हमने software के बारे में बात की है। हमने resource allocation के बारे में बात की है, और हमने सभी processes और specific processes दोनों के लिए config के विभिन्न scopes का वर्णन किया है।

## 3. Workflow parameters store करने के लिए एक parameter फ़ाइल का उपयोग करें

ठीक है, अब हम अपना ध्यान parameters की ओर मोड़ने जा रहे हैं। हम config फ़ाइल में parameters को define कर सकते हैं जैसा कि हमने पहले Nextflow script में किया था। तो params dot greeting equals hello या params scope का उपयोग करें और foo equals bar सेट करें।

और यह आपके workflow के लिए defaults सेट करने के लिए बहुत अच्छा है। हालांकि, जब आप pipelines चला रहे हैं, तो JSON या YAML फ़ाइल में parameters specify करना अच्छा हो सकता है।

इस तरह की फ़ाइल का उपयोग करना dash dash के साथ command line options specify करने से बेहतर है। जब आप एक workflow चलाते हैं, तो आपको कई parameters specify करने पड़ सकते हैं और उन सभी को एक single CLI पर लिखना tedious और error prone हो सकता है। साथ ही, यह संभावना नहीं है कि आप उन सभी parameters को याद रखेंगे जिनका आपने उपयोग किया था, इसलिए यदि आप उसे एक फ़ाइल में code करते हैं, तो भविष्य में समान parameters का उपयोग करके workflow को फिर से launch करना आसान है।

हमारे पास यहां test params नामक एक example फ़ाइल है और आप देख सकते हैं कि यह हमारे workflow में तीन अलग-अलग values के साथ तीन parameters specify करती है। व्यक्तिगत रूप से, मुझे YAML को JSON की तुलना में लिखना आसान लगता है। तो सिर्फ यह demonstrate करने के लिए कि यह काम करता है, मैं Test yaml नामक एक नई फ़ाइल बनाने जा रहा हूं और इन्हें copy करूंगा, quotes से छुटकारा पाऊंगा। और save दबाऊंगा।

ये JSON और YAML फ़ाइलें लिखने में आसान हो सकती हैं क्योंकि ये अधिक परिचित syntax हैं। लेकिन ध्यान दें कि ये केवल parameters के लिए हैं और वे केवल इस तरह key value syntax लेते हैं।

## 3.1. Parameter फ़ाइल का उपयोग करके workflow चलाएं

चलिए इसे आज़माते हैं। पहले जैसी ही कमांड करें। report से छुटकारा पाएं और मैं dash params file test params yaml करने जा रहा हूं।

नहीं, यह एक core Nextflow option है, इसलिए यह एक single hyphen है।

ठीक है। इसने workflow चलाया और इसने उस YAML फ़ाइल में parameters का उपयोग किया बजाय मैं उन सभी को command line पर specify करूं। सिर्फ इस simple example के लिए overkill लग सकता है, लेकिन आप कल्पना कर सकते हैं कि यदि आपके पास 10 या 20 विभिन्न parameters हैं, तो manually type करना एक pain हो सकता है, और यह code editor में edit करना और reproducibility के लिए पकड़े रखना बहुत आसान है।

## 3. निर्धारित करें कि काम करने के लिए किस executor(s) का उपयोग किया जाना चाहिए

ठीक है। हमने Docker और conda के साथ software packaging के बारे में बात की है। हमने CPUs और memory के साथ process resource requirements के बारे में बात की है। और हमने workflows चलाते समय parameters को specify करने के तरीके के बारे में थोड़ी बात की।

Configuration के अंतिम भाग वास्तव में execution हैं, underlying compute infrastructure ही, और यह Nextflow का real jewel in the crown है: कि हम इन same workflow को कई विभिन्न compute infrastructures में चला सकते हैं।

मैं वास्तव में एक second के लिए written training material पर switch करने जा रहा हूं। प्रशिक्षण के इस भाग के तहत, हम कुछ अलग-अलग examples देख सकते हैं कि विभिन्न executors, इस मामले में, HPC schedulers, एक job submit करने के लिए आवश्यक resource requirements को कैसे define करते हैं।

तो Slurm के लिए, आपके पास ये SBATCH headers हैं, जो dash dash mem और CPU number को define करते हैं। यदि आप PBS का उपयोग कर रहे हैं, तो आपके पास अलग headers हैं, और यदि आप Grid Engine का उपयोग करते हैं, तो आपके पास फिर से अलग headers हैं।

आप कल्पना कर सकते हैं कि यदि आप cloud पर चलाना चाहते हैं, तो यह और भी अलग है, चाहे वह AWS batch हो, Google Cloud, Azure, या अधिक।

इन underlying compute infrastructures में से प्रत्येक को executor कहा जाता है और Nextflow जानता है कि correct syntax के साथ jobs submit करने के लिए इन सभी विभिन्न executors से कैसे बात करनी है।

अच्छी खबर यह है कि आपको इसके बारे में जानने की ज़रूरत नहीं है। आपको बस Nextflow को बताना है कि किस executor का उपयोग करना है।

## 3.1. एक अलग backend को लक्षित करना

हम अपनी config फ़ाइल में वापस जाते हैं और process हम executor करते हैं, और मैं local type करने जा रहा हूं।

Local वास्तव में default है, यदि आप कोई अन्य executor specify नहीं करते हैं, तो local का उपयोग किया जाएगा, और इसका मतलब सिर्फ आपका host system है, जहां भी आपने Nextflow launch किया,

मैं इसके बजाय, Slurm specify कर सकता था। और वह Slurm jobs submit करेगा, या मैं AWS batch कह सकता था, और वह AWS batch को jobs submit करेगा।

कुछ मामलों में आपको कुछ अतिरिक्त configuration की आवश्यकता होती है, उदाहरण के लिए, cloud पर चलाने के लिए कुछ credentials की आवश्यकता होगी, लेकिन वास्तव में यह इसका core है, और पूरी तरह से अलग compute environment में अपना workflow चलाने के लिए यह config की एक या दो lines जितना सरल हो सकता है।

भले ही हम code spaces के भीतर एक simple system पर चल रहे हैं, मैं अभी भी इसके साथ थोड़ा खेल सकता हूं और दिखावा कर सकता हूं कि हम Slurm पर चल रहे हैं। अगर मैं फिर workflow launch करता हूं, Nextflow run, hello config। यह fail हो जाएगा क्योंकि यह Slurm को jobs submit करने में सक्षम नहीं होगा। लेकिन हम अभी भी work directories में जा सकते हैं और देख सकते हैं कि Nextflow ने क्या किया। तो अगर हम इस work directory में जाएं और Command Run को देखें। आप फ़ाइल के top पर देख सकते हैं, हमारे पास अब ये sbatch header lines हैं, जिन्होंने Slurm job के लिए आवश्यक resources को specify करने की कोशिश की।

## 4. Preset configurations select करने के लिए profiles का उपयोग करें

ठीक है, हम लगभग वहां हैं। इस अध्याय का अंतिम भाग configuration profiles के बारे में बात कर रहा है। यदि आप कई विभिन्न systems पर अपना pipeline चला रहे हैं, तो ये सभी विभिन्न Nextflow config फ़ाइलें होना annoying हो सकता है, जिन्हें आपको हर बार specify करने की आवश्यकता होती है।

इसके बजाय, आप अपनी Nextflow config फ़ाइल के भीतर configuration की groupings को encode कर सकते हैं, और उन groups को एक profile flag का उपयोग करके on और off switch कर सकते हैं। चलिए देखते हैं कि वह कैसा दिखता है।

## 4.1. Local development और HPC पर execution के बीच switching के लिए profiles बनाएं

हम यहां अपने example में दो profiles बनाने जा रहे हैं, एक मेरे laptop के लिए और एक भारी HPC system के लिए। मैं थोड़ा cheat करने जा रहा हूं और बस training material से code को copy करूंगा और इसे यहां रखूंगा।

हमारे पास profiles नामक एक नया scope है, और फिर हमारे पास प्रत्येक profile के लिए एक name है, जो कुछ भी हो सकता है। और उसके भीतर हमारे पास configuration है, जो बिल्कुल top level config की तरह दिखता है जो हमने पहले ही लिखा था। तो फिर से, हमारे पास process scope है। Docker scope।

My laptop नामक profile पर। मैं कह रहा हूं कि local executor का उपयोग करके चलाएं, तो मेरे host system पर और Docker का उपयोग करें।

यहां university HPC profile पर मैं कह रहा हूं कि jobs submit करने के लिए Slurm का उपयोग करें, Docker के बजाय conda का उपयोग करें, और मैं विभिन्न resource limits specify कर रहा हूं, जो मेरे द्वारा उपयोग किए जा रहे HPC पर nodes के system size से match हो सकते हैं।

By default, जब मैं Nextflow चलाता हूं तो इस configuration में से कोई भी उपयोग नहीं किया जाएगा, मुझे specify करना होगा कि मैं इनमें से किसी एक profile का उपयोग करना चाहता हूं।

## 4.2. एक profile के साथ workflow चलाएं

चलिए nextflow run hello config करते हैं। और मैं dash profile करने जा रहा हूं, single hyphen क्योंकि यह एक core Nextflow option है। और फिर जो name मैंने दिया, जो है my laptop। Nextflow को अब उस configuration profile के भीतर specify किए गए config के block का उपयोग करना चाहिए, और इसे apply करना चाहिए जब यह Nextflow चलाता है। यदि मैं दूसरे config block का उपयोग करना चाहता था, तो मुझे बस उस profile name को switch करना होगा। याद रखना बहुत आसान। उपयोग करना बहुत आसान।

## 4.3. एक test profile बनाएं

ध्यान दें, profiles में किसी भी तरह का configuration हो सकता है, इसलिए इसे आपके execution environment से संबंधित होना जरूरी नहीं है। उदाहरण के लिए, चलिए यहां एक नया profile बनाते हैं, जिसमें parameters का एक set है। हम इसे tux में बदल सकते हैं और my profile में बदल सकते हैं, और अब जब हम profile test करते हैं, तो यह इन parameters को specify करने जा रहा है, जो workflow के top level पर specify किए गए parameters को overwrite करेगा।

जब आप Nextflow चलाते हैं, तो आप कई profiles को chain कर सकते हैं और वे क्रम में apply होंगे।

## 4.4. Test profile के साथ locally workflow चलाएं

तो मैं पिछली कमांड ले सकता हूं और comma test कर सकता हूं। वह पहले my laptop config apply करेगा, और फिर यह test config apply करेगा। यदि कोई overlap है, तो right पर profile किसी भी पिछले profiles में configuration को overwrite कर देगा। अगर मैं enter दबाता हूं, तो देखते हैं क्या होता है।

ठीक है, हमें यहां एक नई results फ़ाइल मिली है। आप My Profile देख सकते हैं, जिसे मैंने options में से एक के रूप में specify किया था। और हम cowpy, my profile भी देख सकते हैं, और निश्चित रूप से, वहां tux है। तो यह काम किया।

## समापन

ठीक है! अद्भुत। बस इतना ही। आपने कोर्स के अंत तक पहुंच गए हैं। आपको थोड़ी celebration confetti मिलती है। इस अध्याय को समाप्त करने के लिए बधाई।

[अगली वीडियो ट्रांसक्रिप्ट :octicons-arrow-right-24:](07_next_steps.md)
