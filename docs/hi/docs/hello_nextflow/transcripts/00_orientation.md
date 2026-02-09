# ओरिएंटेशन - वीडियो ट्रांसक्रिप्ट

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "महत्वपूर्ण नोट"

    यह पेज केवल ट्रांसक्रिप्ट दिखाता है। पूर्ण चरण-दर-चरण निर्देशों के लिए, [कोर्स सामग्री](../00_orientation.md) पर वापस जाओ।

## स्वागत

नमस्ते, और Hello Nextflow में स्वागत है। मेरा नाम Phil Ewels है। मैं Seqera में Open Source Software के लिए Product Manager हूँ, जो Nextflow के पीछे की कंपनी है।

यह कोर्स Nextflow के साथ workflows बनाने का एक hands-on परिचय है। यह उन लोगों के लिए डिज़ाइन किया गया है जो Nextflow में बिल्कुल नए हैं और अपनी खुद की pipelines विकसित करना चाहते हैं।

सभी उदाहरण सरल text processing के हैं, ताकि तुम domain expertise की आवश्यकता के बिना Nextflow concepts पर ध्यान केंद्रित कर सको, बस कुछ command line की जानकारी चाहिए।

हम Nextflow की मूल बातों से गुजरने वाले हैं: processes लिखना, उन्हें multi-step workflows में जोड़ना, containers के साथ software dependencies को manage करना, और विभिन्न computing environments के लिए pipelines को configure करना। अंत तक, तुमने शुरू से एक working pipeline बना ली होगी।

यह कोर्स pipelines _विकसित करने_ पर केंद्रित है। यदि तुम केवल मौजूदा pipelines को _चलाना_ चाहते हो बिना code में बहुत अधिक गोता लगाए, तो हमारे पास एक छोटा "Nextflow Run" कोर्स है जो तुम्हारे लिए बेहतर हो सकता है।

एक बार जब तुम यहाँ मूल बातें सीख लो, तो हमारे पास follow-on courses भी हैं जो इन concepts को वास्तविक scientific analysis में लागू करते हैं। हम तुम्हें सिखाएंगे कि nf-core community की pipelines और best practices का उपयोग कैसे करें।

यदि तुम फंस जाते हो तो community.seqera.io पर जाओ। वहाँ एक सक्रिय community forum है जिसमें training questions के लिए समर्पित एक section है। तुम इसे किसी भी समय उपयोग कर सकते हो, हालाँकि, हम quarterly training weeks भी चलाते हैं जिसमें विशेष रूप से मदद के लिए लोग उपलब्ध होते हैं। इसलिए यदि तुम उनमें से किसी एक के दौरान training कर रहे हो, तो निश्चित रूप से शर्माओ मत और मदद माँगो।

तुम Seqera AI से भी मदद माँग सकते हो। यह Nextflow code को समझाने और debugging में मदद करने में बहुत अच्छा है।

जब तुम Nextflow को scale पर चलाने के लिए तैयार हो, तो Seqera Platform ऐसा करने के लिए सबसे अच्छी जगह है। यह बिना किसी vendor lock-in के तुम्हारे infrastructure पर चलता है, pipeline launching से लेकर real-time monitoring तक, interactive analysis environments तक सब कुछ के साथ। लेकिन अभी के लिए, चलो बस fundamentals पर ध्यान केंद्रित करते हैं।

ठीक है, चलो शुरू करते हैं।

## training.nextflow.io

ठीक है। पहली बात यह ध्यान देने योग्य है कि training.nextflow.io पर सभी training courses बहुत interactive हैं। विचार यह है कि तुम training material और मेरे निर्देशों का पालन करो, और हम एक साथ training material से गुजरें। तो तुम्हें दो चीजों की आवश्यकता होगी: तुम्हें अपने laptop की आवश्यकता होगी और तुम्हें यह website खुली रखनी होगी। और बस इतना ही।

तो यह homepage है जैसा कि आज दिखता है जब मैं इसे record कर रहा हूँ। तुम देख सकते हो कि एक overview है, विभिन्न चीजों का, background का, और विभिन्न courses का जो हमारे पास हैं, जिसकी list हर समय बढ़ती रहती है।

Nextflow for newcomers वह जगह है जहाँ हम हैं। इसमें दो courses हैं, Nextflow Run, जो एक अलग course है, और, Hello Nextflow, जिसकी हमें परवाह है।

और तुम sidebar पर सभी विभिन्न courses भी देख सकते हो। मैं Hello Nextflow पर jump कर सकता हूँ, और हम सभी विभिन्न chapters देख सकते हैं जिन पर हम एक साथ काम करने वाले हैं।

यहाँ कुछ अन्य महत्वपूर्ण चीजें ध्यान देने योग्य हैं। सबसे पहले, training material versioned है, तो तुम यहाँ ऊपर देख सकते हो। यह कहता है 3.0 latest, जो record करते समय latest stable version है। यह समय के साथ बदलेगा। हम नए courses push करते हैं और समय के साथ material को update करते हैं। इसलिए यदि यह 3.1 या 3.2 है, तो बहुत चिंता मत करो। यदि यह 4.0 है, तो शायद एक नया video है, और तुम्हें शायद जाकर वह ढूंढना चाहिए क्योंकि शायद महत्वपूर्ण updates होंगे।

शीर्ष पर एक और dropdown यह language वाला है। अब यह version 3.0 के लिए बिल्कुल नया है। हमने पहले से translated material लिया है, जो Humans द्वारा, हाथ से किया गया था, और हमने उसे एक LLM में pass किया है और LLM translation का उपयोग करके training material के विभिन्न translations को maintain करने के लिए यह पूरा नया infrastructure set up किया है।

तो अब हमारे पास यहाँ ये सभी शानदार translations हैं। इसलिए यदि तुम Korean में सुनना चाहते हो, तो तुम पूरी website को Korean में load कर सकते हो। और, वहाँ follow कर सकते हो। इन सभी अन्य languages के लिए भी, Hindi और German और इसी तरह। मैं English में follow करने वाला हूँ। यह वह primary language है जिसमें हम material लिखते हैं।

कुछ अन्य buttons यदि तुम्हें light mode पसंद है। dark mode के बजाय, तुम यहाँ शीर्ष पर light mode में website को follow कर सकते हो।

और फिर हम जो कुछ भी देखते हैं वह एक single GitHub repository में है, जो open source है, जिसे nextflow-io/training कहा जाता है। और यदि तुम किसी भी समय इस button पर click करते हो, तो यह GitHub repository पर जाएगा। हम एक मिनट में उस पर वापस आएंगे।

## GitHub Codespaces सेट अप करना

ठीक है, तो अब तुमने इसे browser tab में खोल लिया है। चलो Hello Nextflow पर जाते हैं और click करते हैं। तुम intro page पर देख सकते हो, यह हमें कुछ requirements, overview, और lesson plan बताता है कि हम लगभग क्या cover करने वाले हैं, और फिर हम getting started में dive करने वाले हैं।

इस interactive tutorial को करने के विभिन्न तरीके हैं। यदि तुम खुश हो, तो तुम अपने खुद के computer पर अपनी खुद की Nextflow installation के साथ locally यह कर सकते हो। यदि हम Environment Options पर click करते हैं, तो तुम देख सकते हो कि local Devcontainers का उपयोग करके या तुम manual installation के साथ locally सभी software को install करके भी यह कैसे करना है, इसके बारे में अधिक विवरण हैं।

हम इसे Seqera Studios के साथ अच्छी तरह से काम करने पर काम कर रहे हैं, तो यह एक और option है। लेकिन अभी सबसे आम एक GitHub Codespaces का उपयोग करना है।

Codespaces, GitHub द्वारा चलाए जाने वाले remote server पर एक sandbox environment set up करता है। और यह एक निश्चित मात्रा के उपयोग के लिए free है, जो आमतौर पर training के लिए ठीक है। और यह तुम्हें एक VS Code instance, एक IDE के साथ set up करेगा जहाँ तुम repository से सभी files को access कर सकते हो, Nextflow चला सकते हो और सब कुछ। और हमने तुम्हारे लिए Codespaces को pre-configure किया है। तो इसमें वह सब कुछ है जो तुम्हें चाहिए।

इसकी खूबसूरती यह है कि Codespace set up करने के लिए बस एक click है। यह सभी के लिए समान है, और हम जानते हैं कि तुम्हारे पास सभी prerequisites पहले से installed हैं, इसलिए यह अच्छा और तेज़ है।

तो पहली चीज़ "Getting Started" पर जाना है। इस button को देखो, जो कहता है, _Open in Codespaces_। मैं इसे एक नए tab में खोलने के लिए command + click करने वाला हूँ, और यह हमें GitHub पर ले जाता है।

यह ऐसा दिखता है। हम देख सकते हैं, हमने यहाँ तुम्हारे लिए सभी options set किए हैं। यदि तुम चाहो, तो तुम change options पर click कर सकते हो। यहाँ कुछ चीजें तुम कर सकते हो। तुम एक बड़ा instance machine दे सकते हो, उदाहरण के लिए, यदि तुम पाते हो कि यह crash हो जाता है क्योंकि यह memory से बाहर हो जाता है या कुछ भी ऐसा। या training material के specific versions set कर सकते हो। लेकिन आमतौर पर तुम बस उसके साथ जा सकते हो जो हमने यहाँ set up किया है और तुम इसे देख सकते हो। इस मामले में यह, 3.0 release का उपयोग कर रहा है।

तो मैं create new Codespace पर click करने वाला हूँ। और वह मुझे अंदर ले जाता है।

ध्यान दें, यह भी कहता है no Codespace to resume वहाँ। यदि मैंने पहले एक Codespace बनाया है, तो training material पर उस button पर फिर से click करने से मुझे उसी page पर ले जाएगा और यह सभी Codespaces को list करेगा जो मेरे पास पहले से चल रहे हैं। फिर तुम सीधे उनमें वापस jump कर सकते हो और जहाँ तुमने छोड़ा था वहाँ से जारी रख सकते हो। तो कोई बात नहीं यदि तुमने अपना laptop बंद कर दिया।

वे कुछ मिनटों की inactivity के बाद automatically खुद को shut down कर लेते हैं, लेकिन कोई समस्या नहीं। तुम बस उन्हें restart कर सकते हो।

एक बार जब तुम एक नया Codespace start करते हो, तो यह इस page पर इस तरह बैठने वाला है और यह काफी देर तक load होने वाला है। तो अब एक quick break लेने का अच्छा समय है। शायद तुम toilet जाना भूल गए या तुम शुरू करने से पहले एक cup tea चाहते हो? अभी जाओ जब तुम इसके लिए wait कर रहे हो, क्योंकि यह थोड़ी देर के लिए वहाँ spin करने वाला है।

बस जल्दी से जब हम इसके load होने का wait कर रहे हैं, मैं github.com/codespaces पर भी जाने वाला हूँ और बस दिखाता हूँ कि यह overview page है जहाँ तुम सभी विभिन्न Codespaces देख सकते हो जो तुम्हारे पास वर्तमान में चल रहे हैं।

तुम देख सकते हो मेरे पास यहाँ nextflow-io/training के लिए एक है। कोई changes नहीं, क्योंकि मैंने अभी तक इसमें कुछ नहीं किया है। यह कितने resources का उपयोग कर रहा है, और तुम देख सकते हो कि इस समय यह setting up है। मैं यहाँ जा सकता हूँ, इस छोटे dropdown पर click कर सकता हूँ और delete पर click कर सकता हूँ। इसलिए यदि तुमने गलती से multiple Codespaces set up कर दिए हैं और तुम कुछ का उपयोग नहीं कर रहे हो, तो तुम पुराने को delete कर सकते हो और clean up कर सकते हो।

अंत में, इसमें जाने का एक और तरीका। यदि हम GitHub repository पर जाते हैं। और यह किसी भी GitHub repository के लिए काम करता है। code पर click करो। तुम्हारे पास locally repository को clone करने के लिए commands हैं। और Codespaces नाम का एक tab है। और फिर से, तुम एक नया बना सकते हो, और तुम उन्हें देख सकते हो जो पहले से चल रहे हैं।

तो फिर से, यदि तुम भूल जाते हो कि तुमने अपना Codespace कैसे बनाया, तो तुम हमेशा इस तरह से इसमें वापस आ सकते हो।

## VS Code interface

ठीक है, builders ने finish कर दिया है और अब यह GitHub Codespaces को load करना शुरू कर रहा है। यह हमेशा इतना लंबा नहीं लगता, तो चिंता मत करो। यह बस पहली बार है जब तुम Codespace बनाते हो। यदि तुम एक में वापस jump करते हो जो पहले से मौजूद है, तो यह बहुत तेज़ है।

बहुत अधिक अधीर मत बनो यदि यह पहली बार है, यह अभी तक finish नहीं हुआ है, भले ही यह हमें एक interface देना शुरू कर रहा है।

लेकिन जब हम final चीजों के set up होने का wait कर रहे हैं, मैं बस तुम्हें interface के माध्यम से ले जाऊंगा यदि तुम VS Code से थोड़ा unfamiliar हो।

सबसे पहले, AI stuff के लिए chat sidebar है, जिसकी हमें आवश्यकता नहीं है। तो मैं उसे close करने वाला हूँ, उससे छुटकारा पाने वाला हूँ और कुछ space free करने वाला हूँ।

बाईं ओर, हमारे पास file explorer है जो हमें Git repository में सभी files दिखाता है, जो workspace है जिसे हमने बनाया है। ध्यान दें, ये local files नहीं हैं। यह सब remote server पर है जहाँ हम काम कर रहे हैं। तुम local files को drag और drop कर सकते हो और चीजें, लेकिन अधिकांश भाग के लिए, हम आज उसके बारे में नहीं सोचने वाले हैं। हम बस purely remotely काम कर रहे हैं।

इस sidebar में अन्य tools हैं, उदाहरण के लिए, search। तो तुम एक बार में repository में सभी files को search कर सकते हो। और यदि हम training repo पर development work कर रहे होते, तो हम Git के साथ source control के साथ integration कर सकते थे और debugging और अन्य चीजें।

अन्य चीजें हैं, यहाँ ऊपर एक main kind of code editing window है, जिसने अभी readme का एक preview load किया है, जो training material के लिए है। तो इस मामले में यह markdown देख रहा है, लेकिन आमतौर पर यह एक code editor होगा।

और फिर उसके नीचे हमारे पास terminal है, जहाँ हम अपने सभी commands चलाने वाले हैं और सीधे Nextflow के साथ interact करने वाले हैं।

Codespace में सब कुछ pre-installed है, तो Nextflow command पहले से वहाँ है और इसी तरह।

ठीक है। जब तुम इतनी दूर पहुँच जाते हो, तो यह लगभग done होना चाहिए। तुम अब देख सकते हो कि इसने Nextflow language server को download कर लिया है और इसने हमारे लिए VS code में कुछ extensions set up किए हैं, जिसमें Nextflow extension शामिल है, जो उपयोगी होने वाला है। तो मैं उसे close कर सकता हूँ और मैं README.md को close कर सकता हूँ।

और अब तुम देख सकते हो कि मुझे बाईं ओर कुछ और मिला है। मैं यहाँ थोड़ा zoomed in हूँ, लेकिन यदि मैं zoom out करता हूँ तो तुम देख सकते हो कि एक button कहता है Nextflow with the Nextflow icon। और उसमें project को explore करने के लिए कुछ अच्छी चीजें हैं और चीजें, जिन पर हम बाद में course में वापस आएंगे।

ठीक है। यदि तुम कभी इनमें से कोई भी panels खो देते हो, तो ऊपर दाईं ओर ये buttons वास्तव में उपयोगी हैं और ये बस चीजों को show और hide करते हैं। तो वह Explorer को shows और hides करता है, नीचे terminal को shows और hides करता है। और इसी तरह।

मैं इनका बहुत उपयोग करने वाला हूँ क्योंकि मैं बहुत zoomed in हूँ, तो तुम्हें मेरी screen पर सभी text देखने में मदद करने की कोशिश करता हूँ, और इसलिए terminal के साथ full screen जाने में सक्षम होना उपयोगी है और फिर जब हम code देख रहे हों तो इसे hide करना। लेकिन अधिकांश समय तुम बस इस सभी चीजों को एक ही समय में खुली रख सकते हो।

ठीक है, और क्या देखना है? बहुत अधिक नहीं। ध्यान दें कि Nextflow, जैसा कि मैं कहता हूँ, installed है। तो मैं "nextflow -version" type कर सकता हूँ और यह कहते हुए आना चाहिए कि हमारे पास कौन सा version installed है।

यहाँ कुछ अन्य चीजें भी installed हैं। हर chapter के अंत में, हमारे पास quiz questions का एक set है, उदाहरण के लिए, website पर। और तुम उन्हें terminal में भी कर सकते हो यदि तुम चाहो तो quiz type करके।

कुछ अन्य keyboard shortcuts हैं जिनका मैं उपयोग करने वाला हूँ, बस यदि तुम curious हो। उदाहरण के लिए, अभी मैंने अपने Mac पर cmd+K दबाया और उसने terminal को clear कर दिया, सभी पिछले output से छुटकारा पाने के लिए। तो यह चीजों को clean रखने के लिए अच्छा है। यदि तुम मुझे ऐसा करते हुए देखते हो तो मैं इसे इस तरह कर रहा हूँ।

और यदि तुम terminal में नए हो, तो याद रखो कि तुम auto complete के लिए tab का उपयोग कर सकते हो, जो मैं paths को auto complete करने के लिए बहुत कर रहा हूँगा।

तो मैं यहाँ बाईं ओर देख सकता हूँ कि Hello Nextflow नाम का एक folder है, जिस पर हम काम करने वाले हैं। यदि मैं files को list करने के लिए "ls" करता हूँ, तो मैं "hel" कर सकता हूँ, tab hit कर सकता हूँ, auto completes। और तो यह paths को complete करने का एक बहुत तेज़ तरीका है।

## केवल Hello Nextflow folder खोलना

ठीक है। यह बहुत अच्छा है। हालाँकि इस repository में बहुत सारी चीजें हैं।

website generate करने के लिए सभी files हैं, और यहाँ multiple विभिन्न courses हैं, और तुम इस route से कर सकते हो और बस "Hello Nextflow" folder में click कर सकते हो। लेकिन वास्तव में बस इस पर purely focus करना अच्छा है।

तुम इसे अपने workspace के रूप में यहाँ चारों ओर clicking के एक bunch के साथ और एक project directory set करने और stuff के साथ set कर सकते हो। लेकिन सबसे आसान तरीका code type करना है, जो VS Code को launch करने के लिए CLI command है, और फिर "hello-nextflow"।

वह एक नया browser tab खोलेगा और तुम पुराने को close कर सकते हो। और यह बिल्कुल वैसा ही दिखता है। लेकिन अब तुम देख सकते हो कि हम इस subdirectory में हैं और सभी अन्य files invisible हैं, और हमारे पास एक cleaner kind of setup है।

तुम यहाँ देख सकते हो कि current working directory भी अब Hello Nextflow folder के भीतर है। तो nice और clean। हमें गलत जगह पर होने के बारे में चिंता करने की आवश्यकता नहीं है। ठीक है।

## 2026 के लिए नया Nextflow Syntax

इस बिंदु पर मुझे एक विशेष चीज़ का उल्लेख करना है। अभी, 2026 की शुरुआत में, हम Nextflow में विभिन्न features लाना शुरू कर रहे हैं, और बड़े नए में से एक Nextflow के अंदर एक नया language syntax parser है।

मूल रूप से वह engine जो तुम्हारी Nextflow files को पढ़ता है और उसे समझता है, runtime के लिए। syntax में कुछ changes हैं, और यह वास्तव में महत्वपूर्ण है कि तुम सही syntax parser enabled के साथ Nextflow का उपयोग करो।

इसके लिए तुम्हें दो चीजों की आवश्यकता है। तुम्हें Nextflow का एक up to date version चाहिए और तुम्हें यह सुनिश्चित करना होगा कि यह enabled है।

यदि मैं फिर से "nextflow -version" करता हूँ, तो तुम देखोगे कि Codespaces 25.10.2 के साथ चल रहा है और 25.10 इस stuff का उपयोग करने में सक्षम होने के लिए minimum version है।

यदि तुम 26.04 का उपयोग कर रहे हो, जो मेरे लिए अभी तक बाहर नहीं आया है, लेकिन जल्द ही आएगा। तो यह default रूप से नए syntax parser को चला रहा होगा, और तुम्हें कुछ और करने की आवश्यकता नहीं है।

लेकिन यदि तुम 25.10 चला रहे हो, तो तुम्हें strict syntax parser को enable करना होगा, जैसा कि इसे कहा जाता है, या v2 syntax parser।

यह एक environment variable के साथ किया जाता है। यह पहले से Codespaces में set है, तो तुम्हें कुछ भी करने की आवश्यकता नहीं है। लेकिन यदि तुम locally चला रहे हो, तो तुम्हें इसे set करना होगा, और मैं "echo $NXF_SYNTAX_PARSER" करके इसे verify कर सकता हूँ, और यह v2 पर set होना चाहिए।

तो यदि तुम locally चला रहे हो, तो बस "export NXF_SYNTAX_PARSER=v2" करो। बस इतना ही सरल। लेकिन याद रखो कि ऐसा करो। क्योंकि अन्यथा तुम कुछ अजीब discrepancies और, errors देखने वाले हो जैसे हम आगे बढ़ते हैं।

यदि तुम Nextflow version और syntax parser के आसपास इस सभी चीजों के बारे में बिल्कुल भी unsure हो, तो सबसे पहले, याद रखो, तुम्हें चिंता करने की आवश्यकता नहीं है यदि तुम Codespaces में हो। सब कुछ ठीक से set up होना चाहिए। लेकिन दूसरा, यदि तुम Nextflow training material पर जाते हो, यदि तुम नीचे जाते हो, version requirements के बारे में बात करते हो, तो यहाँ एक link है जो तुम्हें explore versions के आसपास help page पर ले जाता है, और यह इस सब के बारे में विस्तार से बताता है।

यदि तुम्हारे पास एक moment है तो इसे पढ़ना worth है। क्योंकि यह कुछ अलग-अलग terms को clarify करने में मदद करता है, जो तुम Nextflow का उपयोग करना शुरू करते समय सुन सकते हो। DSL1, DSL2, syntax parser one, syntax parser two, और इसी तरह की चीजें। तो बस उस पर एक check करना worth है और वह कुछ दोहराता है जो मैंने अभी कहा है।

यह वास्तव में उपयोगी भी है यदि तुमने पहले Nextflow code लिखा है और तुम एक refresher के लिए वापस आ रहे हो। यह तुम्हें कुछ चीजें बताता है जो changes करती हैं और तुम्हें Nextflow documentation के parts से link करती हैं, जो तुम्हें बताती हैं कि अपने Nextflow code को कैसे update करें।

## Course files

ठीक है। खुद को परिचित करने के लिए आखिरी चीज़ बस उन files को देखना है, जो इस directory में हैं। तुम या तो sidebar में देख सकते हो या अक्सर training material में, हम tree command का उपयोग करते हैं, -L, जो levels की संख्या है जिसमें देखना है। हम दो कहेंगे, और यदि मैं इसे full screen बनाता हूँ, तो तुम देखोगे कि यह मूल रूप से बिल्कुल वही mirrors करता है जो हम sidebar पर देखते हैं, लेकिन यह hidden files को exclude करता है, जो एक dot से शुरू होती हैं।

तो \*.nf files, यह Nextflow के लिए खड़ा है। तो ये Nextflow script files हैं, और यहाँ training material के विभिन्न chapters में से प्रत्येक के लिए एक starter file है, जिसे हम खोलेंगे और explore करेंगे और फिर edit करेंगे।

हम इन files को बदलेंगे जैसे हम आगे बढ़ते हैं, और इसलिए हर chapter के अंत तक, files लगभग वैसी ही दिखनी चाहिए जैसी अगले के लिए chapter की शुरुआत में थीं। लेकिन हम तुम्हें ये विभिन्न files देते हैं ताकि तुम हमेशा kind of fresh start कर सको और syntax को mess up करने के बारे में बहुत अधिक चिंता न करो।

यदि तुम्हें किसी ऐसी चीज़ से compare करने की आवश्यकता है जो निश्चित रूप से काम करनी चाहिए। तुम solutions folder में check कर सकते हो, और यह हर एक chapters के लिए एक final state की तरह है, तो तुम जो तुमने लिखा है उसकी तुलना वहाँ जो है उससे कर सकते हो।

एक data directory है। इसमें बस एक greetings.csv file है, जिसे हम course के कुछ में example, input data के रूप में उपयोग करेंगे, और एक config file और कुछ parameters जैसी चीजें, जिन्हें हम बाद में course में describe करेंगे।

## Wrap up

ठीक है, तो अब उम्मीद है कि सब कुछ चल रहा है। तुम्हारी screen मेरी जैसी दिखती है और तुम समझते हो कि सब कुछ कैसे प्राप्त करना है और सभी विभिन्न files क्या हैं।

यदि तुम getting started पर page के नीचे scroll करते हो, तो छोटा checkbox तुम्हें कहना चाहिए कि मैं समझता हूँ कि मैं क्या कर रहा हूँ। मेरा environment up और running है और तुम set हो, तुम्हारी working directory ठीक से "Hello Nextflow" folder पर है।

यदि तुमने उन सभी को tick किया है और वे green दिखते हैं। हम अगले video और अगले chapter पर carry on कर सकते हैं, जो part one है। Hello World। एक moment में मिलते हैं।
