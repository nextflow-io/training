# अगले कदम - वीडियो ट्रांसक्रिप्ट

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=Lp8hS8RdaMwbp5j5&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "महत्वपूर्ण नोट्स"

    यह पेज केवल ट्रांसक्रिप्ट दिखाता है। पूर्ण चरण-दर-चरण निर्देशों के लिए, [कोर्स मटेरियल](../next_steps.md) पर वापस जाओ।

## स्वागत

​

बधाई हो, तुमने यह कर दिखाया!

तुम अंत तक पहुँच गए हो और Hello Nextflow ट्रेनिंग कोर्स पूरा कर लिया है। वाकई उम्मीद है कि तुम्हें यह पसंद आया होगा। अंत तक हमारे साथ बने रहने के लिए धन्यवाद, और हम तुम्हारे समय और मेहनत की सराहना करते हैं जो तुमने Nextflow सीखने में लगाई है। हमें वाकई उम्मीद है कि यह तुम्हारे काम के लिए उपयोगी साबित होगा।

## training.nextflow.io पर अन्य कोर्स

training.nextflow.io पर वापस आते रहना मत भूलो। हम हर समय नए छोटे कोर्स जोड़ रहे हैं और बहुत सारे मटेरियल को भी रिफ्रेश करते रहते हैं जो पहले से यहाँ है। तो यह Hello Nextflow ट्रेनिंग कोर्स समय के साथ अपडेट होता रहेगा।

यह विशेष रूप से महत्वपूर्ण है क्योंकि हम Nextflow में syntax अपडेट कर रहे हैं, और 2026 में काफी सारे नए फीचर्स आएंगे, तो यह कोर्स अगली बार 2027 में थोड़ा अलग दिखेगा और महसूस होगा।

विशेष रूप से, मैं "Nextflow for Science" पेज के लिए एक call out देना चाहता हूँ। ये छोटे कोर्स हैं और इन्हें इस Hello Nextflow कोर्स के बाद फॉलो करने के लिए डिज़ाइन किया गया है। और ये दिखाते हैं कि Nextflow को विशिष्ट अलग-अलग use cases के साथ कैसे उपयोग करें, चाहे वह genomics हो या RNAseq, या कई अलग-अलग चीज़ें। हम हर समय और वैज्ञानिक use cases जोड़ने की कोशिश कर रहे हैं।

Side Quests भी हैं। जब हम Hello Nextflow जैसा कोर्स डेवलप करते हैं, तो बहुत कुछ होता है जिसे हम कवर कर सकते हैं, और सब कुछ scope में रखना मुश्किल होता है। तो अगर कोई विशेष टॉपिक है जो हमें लगता है कि लोगों के लिए दिलचस्प है, जिसे ज़्यादा गहराई से कवर करने लायक है, तो हम उसे Side Quest में डाल देते हैं।

जाकर देखो और अगर अलग-अलग चीज़ें हैं जो तुम्हारे काम के लिए relevant हो सकती हैं, जैसे nf-test या metadata के साथ अलग-अलग चीज़ें करना, और common scripting patterns, Side Quests चेक करो और देखो कि क्या वह और ज़्यादा सीखने के लिए उपयोगी हो सकता है।

nf-core पर भी एक कोर्स है। उम्मीद है कि तुम इस समय तक प्रोजेक्ट से परिचित हो गए होगे, लेकिन अगर नहीं हो, तो जाकर इसे चेक करो। वहाँ अलग-अलग प्रकार के analysis और अलग-अलग प्रकार के data के लिए लगभग 150 अलग-अलग pipelines हैं, तो पूरी तरह से संभव है कि जिस तरह के data analysis की तुम्हें ज़रूरत है उसके लिए एक pipeline out of the box तैयार होगा।

महत्वपूर्ण रूप से, nf-core में components भी हैं, लगभग 1700 अलग-अलग modules, अलग-अलग processes और tools के लिए wrappers। और nf-core के साथ आने वाली tooling से, तुम उन्हें mix और match कर सकते हो और अपनी खुद की pipeline बना सकते हो जैसे Lego bricks। बहुत तेज़ और ज़्यादा reproducible।

## Seqera Platform

जैसे-जैसे तुम Nextflow के साथ अपने उपयोग को बढ़ाते हो, Seqera Platform चेक करो, यह Nextflow चलाने का सबसे अच्छा तरीका है। तुम अपने खुद के infrastructure पर चला सकते हो, तो HPC या AWS, Azure, Google Cloud, Oracle और भी बहुत कुछ। तुम हमारे अपने Seqera Compute का भी उपयोग कर सकते हो अगर तुम बिल्कुल भी computing infrastructure मैनेज नहीं करना चाहते।

Seqera Platform वाकई इन जटिल cloud infrastructures के setup को सरल बनाता है features जैसे Batch Forge के साथ, जो तुम्हारे लिए environment बनाता है। और यह observability और audit logging और compliance के साथ भी वाकई मदद करता है।

यह objectively रूप से pipelines को सस्ता और तेज़ चलाता है technologies जैसे Fusion के साथ, जो disk access और data transfers को optimize करते हैं। और pipeline optimization भी है यह सुनिश्चित करने के लिए कि तुम्हारी pipelines का configuration जितना संभव हो सके tightly tuned है।

pipelines चलाने के अलावा भी पूरी तरह से अलग features हैं। हमारे पास Studios हैं जहाँ तुम interactive analyses चला सकते हो और किसी भी custom docker image से environments बना सकते हो जो तुम बनाते हो। और Data Explorer, जो तुम्हें अपने अलग-अलग file systems explore करने में मदद करता है चाहे वे कहीं भी हों।

Seqera Platform के लिए एक free tier है, तो तुम इन सभी features का उपयोग अभी मुफ्त में कर सकते हो। और हम तुम्हें Seqera Compute के साथ सौ डॉलर का मुफ्त compute credit भी देंगे अगर तुम अपने organizational email address से साइन अप करते हो। अंत में, एक academic program है, तो अगर तुम किसी university में काम कर रहे हो, pricing page चेक करो, वहाँ form ढूंढो और हमें बताओ, और हम तुम्हें Cloud Pro में मुफ्त में upgrade कर देंगे।

## Community सहायता और इवेंट्स

ठीक है। आगे बढ़ते हुए। अगर तुम्हें कभी Nextflow के साथ किसी सहायता की ज़रूरत है, community.seqera.io चेक करो। यह वाकई active है और हम उम्मीद करते हैं कि तुम्हें वहाँ देखेंगे और तुम्हारी अलग-अलग समस्याओं और use cases पर चर्चा करेंगे, और शायद अब तुम कुछ और लोगों की भी मदद कर सकते हो।

हमारे पास बहुत सारे events भी चल रहे हैं। हमारे पास nf-core और Nextflow से आने वाले community events हैं। हमारे पास मार्च में एक online और distributed nf-core hackathon है, पिछले साल हमारे पास एक हज़ार से ज़्यादा लोग शामिल हुए थे दुनिया भर की sites के साथ। तो अगर तुम कर सकते हो तो कृपया हमारे साथ शामिल हो।

और हमारे पास Nextflow Summit events भी हैं, एक Boston में, और फिर हमारे पास Barcelona में और online एक event है। शानदार talks जहाँ तुम लोगों को Nextflow का वाकई बड़े और wild और रोमांचक अलग-अलग तरीकों से उपयोग करते हुए सुन सकते हो। और उनसे जुड़े hackathons और in-person training भी हैं।

## Nextflow Podcast और blog

अगर तुम Nextflow ecosystem में चल रही चीज़ों के साथ updated रहना चाहते हो, तो seqera.io/blog ज़रूर चेक करो।

वहाँ Nextflow के लिए एक section है जहाँ तुम community में काम कर रहे लोगों के community blog posts सुन सकते हो, और Seqera से Nextflow और अन्य tools के updates के बारे में blog posts भी जो हम generate करते हैं।

मैं अपने pet project के लिए भी एक plug देना चाहूँगा, जो कि Nextflow Podcast है। इसे Spotify, या Apple Music, या YouTube पर चेक करो। हम periodically नए episodes डालते हैं जहाँ मैं दूसरे लोगों से बात करता हूँ, या तो Nextflow के साथ काम करने वालों या associated technologies, या community में लोगों से। और हम वाकई technical, deep dives करते हैं कि चीज़ें कैसे काम करती हैं और लोग क्या कर रहे हैं। तो अगर तुम interested हो, उन्हें चेक करो। वे वाकई मज़ेदार हैं।

## धन्यवाद

ठीक है, मैं धन्यवाद का एक set देना चाहूँगा। Seqera की training team इस मटेरियल के लिए ज़िम्मेदार है। मैं camera के सामने बैठा हूँ, लेकिन वाकई में सारी मेहनत उन दूसरे लोगों ने की है। विशेष shout out Geraldine के लिए, जिन्होंने लिखा है और इस training material course को Hello Nextflow और अन्य के लिए refresh किया है। और Jon के लिए भी, जिन्होंने वाकई मदद की है, खासकर नए Nextflow syntax के लिए syntax update करने में और कई कोर्स खुद लिखने में भी। scientific development team में अन्य जैसे Rike, Rob, Florian, और कई अन्य लोगों का भी हम जिस मटेरियल के साथ काम कर रहे हैं उसमें बहुत बड़ा input रहा है।

मैं community में लोगों को भी धन्यवाद देना चाहूँगा। नए translations, उदाहरण के लिए, जो बहुत नए हैं, ambassador program और अन्य जगहों के लोगों द्वारा बहुत प्रभावित किए गए हैं। और वाकई, training material की open source nature का मतलब है कि हमारे पास pull requests और issues काफी बार आते रहते हैं, जो वाकई हमारी मदद करते हैं।

## सर्वे

अब जब तुमने finish कर लिया है, अगर तुमने पहले से नहीं किया है, तो कृपया जल्दी से feedback survey करो। यह training.nextflow.io website पर Hello Nextflow section के ठीक नीचे है।

यह केवल पाँच questions हैं। यह वाकई, वाकई तेज़ है, लेकिन यह हमें roughly track करने की अनुमति देता है कि कितने लोग training कर रहे हैं और तुम हमें यह भी बता सकते हो कि training material को कैसे improve किया जाए। हम वाकई सभी responses चेक करते हैं, तो हम वहाँ तुम्हारे input को वाकई value देते हैं।

## अलविदा

एक बार फिर, इस course में, और इस यात्रा में हमारे साथ शामिल होने के लिए बहुत-बहुत धन्यवाद। अगर तुमने training material में कुछ भी देखा है जिसे तुम्हें लगता है कि improve किया जा सकता है तो एक GitHub issue या Pull Request drop करो। और मैं वाकई उम्मीद करता हूँ कि तुम्हें किसी और Nextflow training course में देखूँगा, या शायद किसी hackathon या event में। फिर से धन्यवाद।​
