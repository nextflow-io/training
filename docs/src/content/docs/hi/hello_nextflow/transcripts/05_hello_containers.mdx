# भाग 5: Hello Containers - वीडियो ट्रांसक्रिप्ट

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता से अनुवादित - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "महत्वपूर्ण नोट"

    यह पेज केवल ट्रांसक्रिप्ट दिखाता है। पूर्ण चरण-दर-चरण निर्देशों के लिए, [कोर्स सामग्री](../05_hello_containers.md) पर वापस जाओ।

    ट्रांसक्रिप्ट में दिखाए गए सेक्शन नंबर केवल सांकेतिक उद्देश्यों के लिए हैं और सामग्री में सभी सेक्शन नंबर शामिल नहीं हो सकते।

## स्वागत और पृष्ठभूमि

नमस्ते, और Hello Nextflow में वापस स्वागत है। यह पार्ट 5 है जिसे Hello Containers कहा जाता है। और कोर्स के इस भाग में, हम इस बारे में बात करने जा रहे हैं कि पाइपलाइन के लिए सॉफ़्टवेयर आवश्यकताओं को कैसे encapsulate किया जाए ताकि पाइपलाइन चलाने वाले लोगों को सॉफ़्टवेयर इंस्टॉल करने के बारे में सोचना न पड़े।

अगर तुम मेरे जितने समय से बायोइन्फोर्मैटिक्स में काम कर रहे हो, तो तुम्हें शायद वो दिन याद होंगे जिन्हें मैं अक्सर बुरे पुराने दिन कहता हूँ, जहाँ जब तुम किसी और की पाइपलाइन चलाना चाहते थे या उनके काम को replicate करना चाहते थे, तो तुम घंटों या दिनों तक उन सभी विभिन्न सॉफ़्टवेयर टूल्स को इंस्टॉल करने की कोशिश में लगे रहते थे जो उन्होंने उपयोग किए थे, उसी versions पर, उन्हें अपनी मशीन पर compile करने की कोशिश करते थे, और यह एक दुःस्वप्न था। यह वाकई मुश्किल था।

अगर तुम HPC पर चला रहे हो, तो तुमने environment modules का उपयोग किया होगा जहाँ sysadmins तुम्हारे लिए सॉफ़्टवेयर इंस्टॉल करने की कोशिश करते थे, जो ठीक था, लेकिन फिर भी अपूर्ण।

लेकिन अब हमारे पास इसे करने के बेहतर तरीके हैं। Nextflow में विभिन्न सॉफ़्टवेयर container तकनीकों के लिए built-in सपोर्ट है। Docker सबसे आम है। वही हम आज उपयोग करने जा रहे हैं। यह Codespaces में अच्छी तरह से काम करता है। यह तुम्हारे local computer पर अच्छी तरह से काम करता है और यह cloud में भी अच्छी तरह से काम करता है।

लेकिन Singularity या Apptainer भी, जो HPC systems पर बहुत आम हैं और प्रभावी रूप से बिल्कुल उसी तरह काम करते हैं। या Podman, Shifter, और भी बहुत सारे हैं जो सभी बहुत समान हैं।

एक अतिरिक्त जो कुछ हद तक समान है लेकिन बिल्कुल नहीं, जिसे Nextflow सपोर्ट करता है वह है Conda। और Nextflow तुम्हारे लिए per process आधार पर Conda environments को manage कर सकता है, जो तुम्हारे अपने Conda environments करने से बहुत बेहतर है। और फिर से, पाइपलाइन के साथ ship किया जा सकता है।

हम इस चैप्टर की शुरुआत container technologies और Docker और वे कैसे काम करते हैं, इसके बारे में थोड़ी बात करके करने जा रहे हैं। और हम पहला आधा हिस्सा Docker में मैन्युअली करने जा रहे हैं ताकि तुम समझ सको कि hood के नीचे क्या हो रहा है और यह कैसे काम करता है। क्योंकि यह समझना वाकई महत्वपूर्ण है कि Nextflow क्या कर रहा है और यह समझना कि तुम्हारा workflow क्या कर रहा है जब इसे execute किया जा रहा है।

तो। चलो हमारे Codespaces में jump करते हैं। अब मैंने सब कुछ फिर से साफ़ कर दिया है, लेकिन अगर हम Hello Containers में जाते हैं, तो तुम्हें दिखना चाहिए कि हमारी सभी scripts और सब कुछ modules चैप्टर के अंत जैसा ही है। तो हमारे यहाँ modules directory में अलग-अलग modules हैं।

वे अभी भी वहाँ हैं। उन्हें वहाँ होना चाहिए ताकि यह चल सके। और workflow और output सभी समान हैं सिवाय इसके कि हमने output publishing path को Hello Containers में बदल दिया है, ताकि तुम्हारी files उस directory में समाप्त हों।

हम इसे अभी चला सकते हैं यह देखने के लिए कि यह काम करता है या नहीं, या हम terminal के साथ आगे बढ़ सकते हैं।

## 1. Container का 'मैन्युअली' उपयोग करो

हम अपने containers को manage करने के लिए Docker का उपयोग करने जा रहे हैं, और मैं "docker -v" करके चेक कर सकता हूँ कि यह मेरे Codespaces पर इंस्टॉल है, जो मुझे version दिखाता है जो इंस्टॉल है और सब कुछ, और यह ठीक से काम कर रहा है।

अब containers और Docker में दो concepts हैं जो वाकई महत्वपूर्ण हैं। एक को image कहा जाता है, और एक को container कहा जाता है। Image पूरे file system का snapshot है, अगर तुम चाहो तो, जिसे तुम उपयोग करोगे, और container चलता हुआ environment है। तो तुम एक image का उपयोग करके एक container बनाते हो।

एक बार जब तुम उस container में होते हो, तो यह आमतौर पर पूरे operating system की तरह काम करता है। यह बाहरी दुनिया से कट जाता है। यह बाकी सब चीज़ों से अलग हो जाता है, और यह एक अच्छी बात है। यही तरीका है जिससे हम Nextflow के साथ इतनी अच्छी reproducibility प्राप्त करते हैं।

क्योंकि container के अंदर चलने वाले tasks के लिए, वे तुम्हारे local system पर किसी भी config files से प्रभावित नहीं होते हैं। किसी भी अन्य external influences से, वे अपने छोटे sandbox में चलते हैं। फिर files बहुत ही reproducible तरीके से produce की जाती हैं क्योंकि तुम same underlying libraries, सभी same dependencies, हर व्यक्ति के लिए हर अलग computing environment पर चलने के लिए बिल्कुल same software का उपयोग कर रहे हो। जो ईमानदारी से मुझे लगता है कि शानदार और अद्भुत है कि यह काम करता है। और आज भी, मुझे अभी भी आश्चर्य होता है कि यह संभव है।

## 1.1. Container image को pull करो

तो हम कुछ Docker images और Docker का उपयोग करने की कोशिश करने जा रहे हैं, जब तुम इसे अपने system पर चलाते हो, तो तुम्हारे computer पर एक docker registry होती है, या इस मामले में, code space में, जो सभी विभिन्न images का track रखती है जो पहले download और उपयोग की गई हैं, और विभिन्न layers जिन पर वे बनाए गए हैं।

हम Docker के साथ locally हमारे पास कौन सी images हैं यह "docker image ls" करके देख सकते हैं। और इस मामले में तुम देख सकते हो कि यहाँ बहुत सारी Docker images हैं, जो सभी इस Codespaces को set up करने से संबंधित हैं। सभी dev containers और चीजों से संबंधित हैं। तो तुम्हें उनके बारे में ज्यादा चिंता करने की जरूरत नहीं है, लेकिन जैसे-जैसे हम और images add करते हैं और उन्हें download करते हैं, जैसे-जैसे यह कोर्स आगे बढ़ता है, तुम उस list को चेक कर सकते हो और तुम देखोगे कि local registry इन सभी चीजों का track रख रही है जो हमने pull की हैं।

लेकिन हम "docker pull" करके एक नई को grab करने जा रहे हैं। और यह Docker को web से एक नई image fetch करने के लिए कहता है।

फिर हम उस container के लिए URI डालते हैं। अब यह एक docker image हो सकती है जिसे तुमने locally build किया है और फिर internet पर push किया है। यह एक image हो सकती है जो किसी और ने बनाई है। Docker images बनाने के कई, कई, कई अलग-अलग तरीके हैं, लेकिन यकीनन सबसे सरल तरीकों में से एक है इसे outsource करना, और किसी और को तुम्हारे लिए यह करने देना।

और हम इस tutorial में जो उपयोग करने जा रहे हैं वह Seqera की एक service है जिसे Seqera Containers कहा जाता है।

अब, Seqera Containers पूरी तरह से free है, और यह एक open source software का उपयोग करती है जिसे हम Wave कहते हैं, जो Nextflow के complementary तरीके से containers को manage करने के लिए बनाया गया था। और यह कई common use cases को handle करता है जो हमें Nextflow के साथ deal करते हुए मिलते हैं।

यह बहुत आम है कि हमें जिस software की जरूरत है वह Conda में, Bioconda में, या conda-forge channels में या अन्य अधिक domain specific channels में packaged है। और Wave और Seqera Containers इससे images build करने में वाकई अच्छा है।

तो मैं इस web UI पर जा सकता हूँ और हम "cowpy" नामक package के साथ mess around करने जा रहे हैं। तो मैं उस package का नाम टाइप करता हूँ जो मुझे चाहिए। यह searches करता है, इसने इसे Python package index पर पाया है, तो मैं इसका उपयोग कर सकता हूँ। या अगर मैं थोड़ा और इंतजार करता हूँ, तो यह bioconda और conda-forge को search कर रहा है। और तुम देख सकते हो, मैं यहाँ कोई भी conda channel specify कर सकता हूँ। तो अगर तुम Nvidia channel या कुछ और ढूंढना चाहते हो, तो वह भी काम करना चाहिए।

और फिर मैं specify कर सकता हूँ कि मैं चाहता हूँ कि यह मेरे लिए एक docker image build करे या एक singularity image और यह भी कि मैं कौन सी CPU architecture चाहता हूँ। तो amd64 या arm64।

और एक बार bioconda results list हो जाने पर, मैं अब सभी अलग-अलग versions भी देख सकता हूँ जो available हैं। मैं इसे डालने जा रहा हूँ। और अब मैं searching जारी रख सकता हूँ और Conda से और packages प्राप्त कर सकता हूँ अगर मैं चाहूँ और इस container को जैसे मैं चाहता हूँ वैसे compose कर सकता हूँ, लेकिन मैं बस वही चाहता हूँ। तो मैं Get Container पर click करने जा रहा हूँ।

अब, किसी और ने पहले same container के लिए पहले ही request किया है और यह एक registry से return किया गया है, तो हमें यह तुरंत मिल जाता है। लेकिन अगर किसी और ने कभी भी इस software package या software packages के इस combination के लिए नहीं पूछा होता, तो Wave और Seqera Containers इसे हमारे लिए on the fly build करते।

हम इस URL को copy कर सकते हैं और हम view build details भी देख सकते हैं। और यह हमें दिखाता है कि service ने backend पर क्या किया। इसने एक conda environment file बनाई। एक docker file, और फिर यह है, docker build process चला रहा है। इसने एक scan भी चलाया, एक security scan, तो तुम कोई भी CVEs देख सकते हो। और यह तुम्हें बताता है कि यह कब बनाया गया था।

Wave और Seqera Containers इससे कहीं ज्यादा कर सकते हैं, लेकिन यह एक सरल use case है, जो सबसे आम है। और मुझे कहना चाहिए कि ये images कम से कम पाँच साल के लिए host की जाती हैं। तो तुम इन URLs को अपनी pipelines में build कर सकते हो और जान सकते हो कि वे जल्द ही कहीं नहीं जाने वाली हैं।

तो मेरे पास cowpy के लिए मेरी docker image के लिए मेरा URL है।

मैं अब "docker pull" उस URL को कर सकता हूँ, और यह सभी विभिन्न layers को fetch करेगा और इस image को download करेगा ताकि यह locally मेरे लिए available हो।

## 1.2. Container का उपयोग cowpy को one-off command के रूप में चलाने के लिए करो

ठीक है, अब चलो इसे वास्तव में उपयोग करने की कोशिश करते हैं। तो अब मैं अब "docker pull" के बजाय "docker run" command का उपयोग करने जा रहा हूँ, और मैं "--rm" flag का उपयोग करने जा रहा हूँ, जो बस Docker को बताता है कि एक बार जब यह मैंने जो पूछा है वह समाप्त हो जाए तो इस container को shut down कर दे। और फिर मैं container के लिए identifier डालता हूँ, जो बस एक URI है।

और फिर अंत में, मैं वह command specify करता हूँ जो मैं चाहता हूँ कि Docker container के अंदर चलाए जो इस image से generated है। मैं बस cowpy कहने जा रहा हूँ, जो tool का नाम है जो Conda Forge से इंस्टॉल किया गया है, जो image के अंदर available है।

मैं enter दबाने जा रहा हूँ और वहाँ तुम जाओ। हमने एक system पर cowpy चलाया है। हमारे पास एक छोटी गाय है जो हमें कुछ जानकारी दे रही है।

अब ध्यान दो कि cowpy मेरे local system पर इंस्टॉल नहीं है। तो अगर मैं इसे बस सभी Docker चीजों के बिना चलाता हूँ, तो यह कहता है, command नहीं मिला। तो इसने एक image pull की है। इसने Docker का उपयोग करके एक container बनाया है, और फिर यह उस container में गया है और हमारे लिए यह command चलाई है और हमें output हमारे terminal पर वापस दिया है। बहुत, बहुत बढ़िया।

## 1.3. Container का उपयोग cowpy को interactively चलाने के लिए करो

ठीक है, हम अब एक कदम और आगे जाने जा रहे हैं और इस container को interactively चलाने जा रहे हैं और थोड़ा poke around करने जा रहे हैं, ताकि हम देख सकें कि container के अंदर क्या हो रहा है।

तो अगर मैं वापस जाता हूँ और मैं अपनी run command लेता हूँ और मैं अंत में cowpy से छुटकारा पाने जा रहा हूँ, क्योंकि मैं वास्तव में cowpy नहीं चलाना चाहता। मैं एक Bash terminal चलाना चाहता हूँ।

और फिर मैं यहाँ वापस जाने जा रहा हूँ और मैं "-it" करने जा रहा हूँ, जो Interactive और Terminal या TTY के लिए है, और मैं enter दबाने जा रहा हूँ।

और अब तुम देख सकते हो कि prompt, जो मैं टाइप करने से पहले का हिस्सा है, बदल गया है। यह Codespaces prompt था जहाँ इसने directory कहा था, और अब यह base और root और tmp कहता है। तो मैं अब container के अंदर हूँ, और अगर मैं "ls" करता हूँ, तो तुम देखोगे कि इस directory में जो files मैं देखता हूँ वे मेरे workspace में मेरी files से अलग हैं।

और वास्तव में, मैं अपने local codespaces workspace या container के अंदर अपनी local drive से किसी भी files को नहीं देख सकता। Docker container runtime, पूरी तरह से isolated है और यह बाहर host file system से कोई भी files write या read नहीं कर सकता।

हालाँकि, मैं software को देख सकता हूँ जो container के अंदर इंस्टॉल है और इसे चला सकता हूँ। तो मैं cowpy चला सकता हूँ और हम cowpy का उपयोग कैसे करें इसके बारे में थोड़ा और देख सकते हैं। यहाँ मैं "cowpy 'Hello World'" कर सकता हूँ और यह, इसे बताता है कि वास्तव में मेरे quote को एक छोटे speech bubble के अंदर रखे। और तुम अलग-अलग types की गायों को भी चला सकते हो, तो यह एक गाय होना जरूरी नहीं है। तुम "-c" कर सकते हो। और मैं Sweden में हूँ, तो मैं एक moose चुनने जा रहा हूँ। बहुत बढ़िया। उसे कुछ antlers दिए।

और बहुत सारे अलग-अलग हैं जिनके साथ तुम खेल सकते हो, जो तुम training docs में described देख सकते हो।

## 1.3.4. Container में data mount करो

ठीक है। यह अच्छा होगा अगर हम अपने file system में files पर cowpy चला सकें।

बेशक, यह बहुत उपयोगी नहीं है कि बस container हो और कुछ भी access न हो। यह safe और reproducible हो सकता है, लेकिन यह बहुत उपयोगी नहीं है।

तो हम यह कैसे करते हैं? मैं exit टाइप करके इस Docker container से बाहर निकलने जा रहा हूँ, और तुम देख सकते हो कि prompt हमें बताता है कि हम अब फिर से अपने regular Codespaces में वापस हैं।

और मैं फिर से same command चलाने जा रहा हूँ। लेकिन इस बार मैं यहाँ कुछ additional flags add करने जा रहा हूँ। और महत्वपूर्ण एक है "-v", जो mounting a volume के लिए है, जो basically एक disk space का हिस्सा है।

"-v" दो भाग लेता है: एक string है और फिर एक colon और एक string। और पहला भाग local file system है, जिसे container में mount किया जाना चाहिए। और फिर दूसरा भाग यह है कि container के अंदर वह कहाँ समाप्त होना चाहिए।

अब मैं यहाँ बस अपना पूरा local file system load करना चाहता हूँ। तो "." current working directory है। तो मैं बस "." करने जा रहा हूँ और फिर ":", और फिर हम इसे container के अंदर एक नई directory में डालने जा रहे हैं जिसे "my_project" कहा जाता है। इसे वास्तव में कुछ भी कहा जा सकता है।

और फिर मैं फिर से चलाने जा रहा हूँ।

Working directory में जहाँ मैं dump किया गया हूँ, जो /tmp है, files वहाँ नहीं हैं। लेकिन अगर मैं "ls my_project" करता हूँ, तो वहाँ हमारे पास है: वही सभी files जो हमारे पास locally Codespaces पर थीं अब container के अंदर उस path पर available हैं।

यह read और write access है तो मैं इस directory में नई files बना सकता हूँ और वे मेरे host file system पर दिखाई देंगी। तो यह particular directory, फिर बिल्कुल वैसे ही behave करती है जैसे मैं container के बाहर था तो मैं अब read और write कर सकता हूँ और चीजें कर सकता हूँ।

## 1.3.5. Mounted data का उपयोग करो

ठीक है, चलो बस साबित करते हैं कि हम यह कर सकते हैं। मैं "cat /my_project/data/greetings.csv" करता हूँ। अगर तुम्हें याद है कि इस file की contents इस तरह दिखती हैं। मैं अब इसे cowpy में pipe कर सकता हूँ और गाय उस file के अलग-अलग outputs को अपने छोटे speech bubble में print करेगी, जो कि काफी fun है।

तो तुम देख सकते हो, हम अब अपने host system पर files के साथ interact करने के लिए container में software का उपयोग कर सकते हैं।

ठीक है, चलो वापस बाहर drop करते हैं और हम training material के बाकी हिस्से के साथ आगे बढ़ेंगे।

## 2. Nextflow में containers का उपयोग करो

तो यह containers का उपयोग करना वाकई बढ़िया है। उम्मीद है कि यह समझ में आता है। और तुम इन containers के value और यह analysis software चलाने के लिए क्यों उपयोगी है यह देख सकते हो।

लेकिन हम Nextflow के अंदर यह पूरी same process कैसे करते हैं? हम खुद बहुत सारे Docker commands नहीं चलाना चाहते। हम बस Nextflow को हमारे लिए यह सब handle करने देना चाहते हैं।

तो चलो इसके माध्यम से काम करते हैं। हम cowpy चलाने के लिए अपनी pipeline में एक नई process add करने जा रहे हैं। ठीक है, तो चलो अपनी नई process के लिए एक नया module बनाते हैं। तो modules में जाओ, चलो इसे cowpy.nf कहते हैं, और फिर मैं यहाँ training material से code copy करने जा रहा हूँ।

लेकिन तुम देख सकते हो कि process बहुत सरल है। यह अब तक हमने जो किए हैं उनकी तरह दिखता है, हमारे पास एक input block है जिसमें एक path है, जो हमारी input file है, और एक value भी है ताकि यह एक character हो, तो हम फिर से एक moose उपयोग कर सकते हैं अगर हम चाहें।

और फिर एक output, जो यहाँ एक single file है, एक path और फिर एक script। और हम वही कर रहे हैं जो हमने container के अंदर interactively किया था: हम input file को read करने के लिए "cat" कर रहे हैं। हम उस contents को cowpy पर pipe कर रहे हैं। हम उस input के आधार पर एक specific character चुन रहे हैं, हम cowpy input file नामक एक output file में लिख रहे हैं, जिसे फिर output पर echoed किया जाता है।

बढ़िया। चलो इसे include करते हैं। तो include \{ COWPY \} from "./modules/cowpy.nf", क्या मैंने इसे cowpy कहा था? हाँ।

और फिर चलो workflow के main block में नीचे अपनी नई process को call करते हैं। तो चलो cowpy चलाते हैं। और हम अपनी नई COWPY process लेंगे और हम collectGreetings.out कहने जा रहे हैं।

और फिर अगर तुम्हें याद है, इस module के लिए दो outputs थे। एक जिसे outfile कहा जाता है और एक जिसे report कहा जाता है। VS Code extension हमारे लिए इन्हें auto-suggest कर रहा है और हम .outfile चाहते हैं।

तुम हमेशा इस process में hop कर सकते हो। या तो इस पर hover करो और यह तुम्हें जल्दी से दिखाना चाहिए कि outputs क्या थे। और हम इसमें command click भी कर सकते हैं और यह module file खोलेगा अगर तुम अधिक detail में देखना चाहते हो।

तो यहाँ हम जाते हैं। वह वहाँ outfile है, और वह path है। तो अब यह हमारी cowpy process के लिए input file होगी। शानदार।

अब अगर तुम्हें याद है, एक cowpy process में दो inputs हैं। हमारे पास character के लिए value channel भी था। तो हम यहाँ "params.character" add कर सकते हैं। मैं इसे hard code कर सकता था अगर मैं चाहता, लेकिन चलो इसे एक CLI option बनाते हैं ताकि हम --character कर सकें।

सही। अब मुझे input parameter को define करने की जरूरत है जिसे हमने अभी call किया है और इसे एक default दें। तो character, String। और मुझे moose पसंद है, तो मैं इसे default रूप से moose पर set करने जा रहा हूँ।

सही, चलो इसे चलाने की कोशिश करते हैं। तो अगर मैं Nextflow run hello containers करता हूँ, तो हम देखेंगे कि क्या होता है।

मैं -resume का उपयोग कर सकता था अगर मेरे पास पुरानी work directories kicking around होतीं। और फिर से, ये पहली processes cached हो गई होतीं और यह थोड़ा तेज़ होता, लेकिन यह basically same होना चाहिए।

अब हम तुरंत देख सकते हैं कि जब यह हमारी नई process पर पहुँचा तो इसने एक error throw की है, यह हमें यहाँ बता रहा है कि cowpy process को execute करने में एक error था और यह exit status 127 के साथ exit हो गया। यह वह command है जिसे इसने चलाने की कोशिश की। यह, यह सही दिखता है, यह वैसा दिखता है जैसा हमने उम्मीद की थी। यह उस output file name ले रहा है, जो सही लगता है, यह इसे एक moose character के साथ चला रहा है और इसे save करने की कोशिश कर रहा है।

लेकिन तुम यहाँ command error देख सकते हो जो कह रहा है कि cowpy command नहीं मिला। और यह समझ में आता है क्योंकि हमने वास्तव में अभी तक Nextflow को एक container उपयोग करने के लिए नहीं कहा है। हमने बस इसे cowpy command दी है। और जैसा मैंने पहले कहा था, cowpy हमारे local system पर इंस्टॉल नहीं है। तो जब इसने इसे चलाने की कोशिश की, तो यह fail हो गया।

## 2.3.1. cowpy के लिए एक container specify करो

हमें Nextflow को बताने की जरूरत है कि एक container available है और यह इसका उपयोग कर सकता है। तो हम यह कैसे करते हैं?

अगर हम अपने module में pop करते हैं, तो हम top पर "container" नामक एक नया declaration add करने जा रहे हैं। और फिर हम इसे एक string पर set करने जा रहे हैं।

अब, अगर तुम्हें याद है, Seqera Containers में, मैं उस URL को copy कर सकता हूँ और मैं बस इसे यहाँ quotes में drop करता हूँ।

अब वापस जाओ और इसे फिर से चलाने की कोशिश करो।

देखते हैं कि यह इस बार काम करता है या नहीं।

दुर्भाग्य से, यह बिल्कुल उसी तरह fail हो जाता है, भले ही अब हमने process चलाने के लिए एक container define किया है। तो अपनी docker image का उपयोग करने के लिए, हमें Nextflow को workflow चलाते समय Docker usage enable करने के लिए कहने की जरूरत है।

और हम एक नई config file बनाकर ऐसा करने जा रहे हैं। तो मैं touch nextflow.config कहने जा रहा हूँ।

यह एक special file name है जहाँ अगर यह working directory में है जब मैं pipeline launch करता हूँ, तो यह automatically load हो जाएगी। तो अगर मैं इस nextflow.config file में जाता हूँ, तो तुम देख सकते हो कि यह वास्तव में पहले से मौजूद है, जिसे मैं भूल गया था। और हमारे पास यहाँ पहले से docker.enabled है, लेकिन यह false पर set है, जो default है।

तो अगर मैं इसे बदलकर equals True कर दूं, docker.enabled। और इन सभी config scopes के लिए Nextflow docs में reference docs हैं। और तुम यह भी देख सकते हो कि जब मैं VS Code extension के साथ hover करता हूँ, तो यह इसके specific docs में pull करता है और मुझे बताता है कि इसका क्या मतलब है और इसे कैसे set करना है।

तो अब हमने इसे true पर set कर दिया है, और अगर मैं Nextflow को फिर से चलाता हूँ, तो Nextflow अब जान जाएगा कि हमारे लिए उस docker image को fetch करना है अगर हमारे पास यह locally पहले से नहीं है, और फिर उस container environment के साथ उस process को execute करना है।

और इसलिए हम देख सकते हैं कि यह successfully चला है और हमारे पास cowpy के बगल में एक छोटा tick है। शानदार। अगर मैं ऊपर जाता हूँ और results directory में देखता हूँ, तो file अभी तक वहाँ नहीं है। और ऐसा इसलिए है क्योंकि हमें अभी भी इस output file को publish करने की जरूरत है जैसे कि, बाकी सभी की तरह।

तो हम workflow के भीतर publish block में जाते हैं, mycowpy equals cowpy.out कहते हैं।

और फिर यहाँ output block में नीचे, mycowpy, squiggly brackets path। ओह। Hello containers। Mode, copy।

अगर मैं अब फिर से चलाता हूँ, तो यह बिल्कुल उसी तरह चलना चाहिए। मैं -resume का उपयोग कर सकता था और मैं हर बार भूल जाता हूँ। और फिर मैं ऊपर जाता हूँ और अब हमारे पास cowpy-COLLECTED नामक एक नई file बनाई गई है, और वहाँ मेरा moose BONJOUR, HELLO, HOLà कह रहा है। शानदार।

अब बेशक मैं अब "--character" भी pass कर सकता हूँ। अलग-अलग options क्या हैं? मुझे लगता है कि एक Turkey है? तो मैं character Turkey का उपयोग कर सकता हूँ। यह बिल्कुल उसी तरह चलने वाला है। मैंने -resume का उपयोग करने का एक और अवसर चूक गया, और अब अगर हम अपनी file load करते हैं और अब हमारे पास एक Turkey है। शानदार।

## 2.3.4. Inspect करो कि Nextflow ने containerized task को कैसे launch किया

ठीक है। अंतिम छोटी बात। चलो बस जल्दी से यह command फिर से चलाते हैं, इस बार resume करते हैं, और work directory में एक quick look लेते हैं यह देखने के लिए कि Nextflow हमारे लिए यह सब काम करने के लिए hood के नीचे क्या कर रहा है।

इस बार यह बहुत तेज़ है, चलो इस work directory में जाते हैं, cd work/. अब अगर तुम्हें याद है कि हमारे पास यहाँ बहुत सारी dot files हैं और जिसमें हम इस मामले में interested हैं वह है जिसे मैंने कहा था कि हमें लगभग कभी देखने की जरूरत नहीं है, जिसे .command.run कहा जाता है।

अगर मैं code .command.run करता हूँ, तो यह इसे editor में खोलने वाला है। और मैं इस file में search कर सकता हूँ और अगर मैं scroll down करता हूँ तो मुझे docker run दिखना चाहिए। और तुम देख सकते हो कि Nextflow हमारे लिए docker run command कर रहा है, जब Docker config में enabled है। इसमें बहुत सारे अलग-अलग, flags और चीजें यहाँ हैं, लेकिन तुम "-v" flag देख सकते हो जिसका हमने खुद उपयोग किया था जब हम चला रहे थे। और तुम देख सकते हो कि यह local, workspace directory को container में mount कर रहा है, ताकि container हमारी input files को access कर सके और outputs को save कर सके। और फिर अंत में, यह .command.sh भी चला रहा है, जो generated script है, जिसमें cowpy command है।

और इसलिए तुम देख सकते हो कि Nextflow workflow logic ले रहा है, जो वह चीज़ है जिसकी हमें वास्तव में परवाह है, जो हमारे analysis के लिए specific है, और यह सभी clever behind the scenes चीज़ें कर रहा है ताकि Docker हमारे system पर काम करे।

और यह इसे वास्तव में portable तरीके से कर रहा है ताकि pipeline का एक end user technology को switch out कर सके जो वे उपयोग कर रहे हैं: Docker, Singularity, Apptainer, Conda। यह वास्तव में pipeline logic के लिए मायने नहीं रखता, लेकिन Nextflow सभी underlying infrastructure needs को handle करेगा, ताकि यह कहीं भी चले।

और यह वास्तव में Nextflow की superpower है। Reproducibility और portability। और Nextflow के साथ तुम वास्तव में अपना workflow share कर सकते हो और अन्य लोग इसे अपने systems पर चला सकते हैं और यह बस काम करेगा।

यह करना वास्तव में, वास्तव में मुश्किल चीज़ है, और अब तुम जानते हो कि अपने workflows के साथ इसे कैसे करना है।

ठीक है, इस chapter के लिए बस इतना ही। अगर तुम course के अंत तक जाते हो, तो तुम्हें, containers के बारे में फिर से एक quiz मिलेगी। उम्मीद है कि यह सब समझ में आया। यह analysis के साथ काम करने का वाकई बढ़िया तरीका है। और अगर तुम containers के लिए नए हो, तो मुझे उम्मीद है कि मैंने तुम्हें convince कर दिया है कि यह जाने का तरीका है, और तुम कभी पीछे मुड़कर नहीं देखोगे।

लेकिन उसके साथ, शायद थोड़ा break लो, और तुम कुछ मिनटों में मेरे साथ Hello Nextflow के अंतिम part 6 से गुजरने के लिए शामिल हो, जो configuration के बारे में है।

बहुत-बहुत धन्यवाद।
