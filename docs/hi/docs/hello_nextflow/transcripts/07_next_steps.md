# अगले कदम - वीडियो ट्रांसक्रिप्ट

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=Lp8hS8RdaMwbp5j5&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "महत्वपूर्ण नोट"

    यह पेज केवल ट्रांसक्रिप्ट दिखाता है। पूर्ण चरण-दर-चरण निर्देशों के लिए, [कोर्स सामग्री](../next_steps.md) पर वापस जाओ।

## स्वागत

​

बधाई हो, तुमने कर दिखाया!

तुम अंत तक पहुँच गए हो और Hello Nextflow प्रशिक्षण कोर्स पूरा कर लिया है। सच में उम्मीद है कि तुम्हें यह पसंद आया होगा। पूरे समय हमारे साथ बने रहने के लिए धन्यवाद, और हम वाकई में उस समय और मेहनत की सराहना करते हैं जो तुमने Nextflow सीखने में लगाई है। हमें सच में उम्मीद है कि यह तुम्हारे काम के लिए उपयोगी होगा।

## training.nextflow.io पर अन्य कोर्स

training.nextflow.io पर वापस आते रहना मत भूलना। हम हर समय नए छोटे कोर्स जोड़ रहे हैं और हम यहाँ पहले से मौजूद बहुत सारी सामग्री को भी रिफ्रेश करते हैं। तो यह Hello Nextflow प्रशिक्षण कोर्स समय के साथ अपडेट होता रहेगा।

यह विशेष रूप से महत्वपूर्ण है क्योंकि हम Nextflow में सिंटैक्स अपडेट कर रहे हैं, और 2026 में काफी सारे नए फीचर्स आएंगे, तो यह कोर्स अगली बार 2027 में थोड़ा अलग दिखेगा और महसूस होगा।

विशेष रूप से, मैं "Nextflow for Science" पेज के लिए एक कॉल आउट देना चाहता हूँ। ये छोटे कोर्स हैं और इन्हें इस Hello Nextflow कोर्स के बाद फॉलो करने के लिए डिज़ाइन किया गया है। और ये दिखाते हैं कि विशिष्ट अलग-अलग use cases के साथ Nextflow का उपयोग कैसे करें, चाहे वह genomics हो या RNAseq, या सभी प्रकार की अलग-अलग चीज़ें। हम हर समय अधिक वैज्ञानिक use cases जोड़ने की कोशिश कर रहे हैं।

Side Quests भी हैं। जब हम Hello Nextflow जैसा कोर्स विकसित करते हैं, तो बहुत कुछ होता है जिसे हम कवर कर सकते हैं, और सब कुछ scope में रखना मुश्किल होता है। तो अगर कोई विशेष विषय है जो हमें लगता है कि लोगों के लिए दिलचस्प है, जो अधिक गहराई से कवर करने योग्य है, तो हम उसे Side Quest में डालते हैं।

जाओ और देखो और अगर अलग-अलग चीज़ें हैं जो तुम्हारे काम के लिए प्रासंगिक हो सकती हैं, जैसे nf-test या metadata के साथ अलग-अलग चीज़ें करना, और सामान्य scripting patterns, Side Quests देखो और देखो कि क्या वह अधिक सीखने के लिए उपयोगी हो सकता है।

nf-core पर भी एक कोर्स है। उम्मीद है कि तुम इस समय तक प्रोजेक्ट से परिचित हो, लेकिन अगर नहीं हो, तो जाओ और इसे देखो। विभिन्न प्रकार के विश्लेषण और विभिन्न प्रकार के डेटा के लिए लगभग 150 अलग-अलग पाइपलाइन हैं, तो यह पूरी तरह से संभव है कि तुम्हें जिस प्रकार के डेटा विश्लेषण की आवश्यकता है, उसके लिए एक पाइपलाइन तैयार हो।

महत्वपूर्ण रूप से, nf-core में components भी हैं, लगभग 1700 अलग-अलग modules, अलग-अलग processes और tools के लिए wrappers। और nf-core के साथ आने वाली tooling के साथ, तुम उन्हें मिक्स और मैच कर सकते हो और Lego bricks की तरह अपनी खुद की पाइपलाइन बना सकते हो। बहुत तेज़ और अधिक reproducible।

## Seqera Platform

जैसे-जैसे तुम Nextflow के साथ अपने उपयोग को बढ़ाते हो, Seqera Platform देखो, यह Nextflow चलाने का सबसे अच्छा तरीका है। तुम अपने खुद के infrastructure पर चला सकते हो, तो HPC या AWS, Azure, Google Cloud, Oracle और भी बहुत कुछ। तुम हमारे अपने Seqera Compute का भी उपयोग कर सकते हो अगर तुम बिल्कुल भी computing infrastructure को manage नहीं करना चाहते।

Seqera Platform वास्तव में इन जटिल cloud infrastructures के सेटअप को सरल बनाता है Batch Forge जैसी features के साथ, जो तुम्हारे लिए environment बनाता है। और यह observability और audit logging और compliance में भी वास्तव में मदद करता है।

यह objectively रूप से पाइपलाइनों को सस्ता और तेज़ चलाता है Fusion जैसी technologies के साथ, जो disk access और data transfers को optimize करती हैं। और pipeline optimization भी है यह सुनिश्चित करने के लिए कि तुम्हारी पाइपलाइनों का configuration जितना संभव हो उतना tightly tuned है।

पाइपलाइन चलाने के अलावा भी पूरी तरह से अलग features हैं। हमारे पास Studios हैं जहाँ तुम interactive analyses चला सकते हो और किसी भी custom docker image से environments बना सकते हो जो तुम बनाते हो। और Data Explorer, जो तुम्हें अपने अलग-अलग file systems को explore करने में मदद करता है चाहे वे कहीं भी हों।

Seqera Platform के लिए एक free tier है, तो तुम अभी इन सभी features का उपयोग मुफ्त में कर सकते हो। और हम तुम्हें Seqera Compute के साथ सौ डॉलर का मुफ्त compute credit भी देंगे अगर तुम अपने organizational email address के साथ sign up करते हो। अंत में, एक academic program है, तो अगर तुम किसी university में काम कर रहे हो, pricing page देखो, वहाँ form ढूंढो और हमें बताओ, और हम तुम्हें मुफ्त में Cloud Pro में upgrade कर देंगे।

## Community help और events

ठीक है। आगे बढ़ते हुए। अगर तुम्हें कभी Nextflow के साथ किसी support की आवश्यकता है, तो community.seqera.io देखो। यह वास्तव में active है और हम उम्मीद करते हैं कि तुम्हें वहाँ देखेंगे और तुम्हारी अलग-अलग problems और use cases पर चर्चा करेंगे, और शायद अब तुम कुछ अन्य लोगों की भी मदद कर सकते हो।

हमारे पास बहुत सारे events भी चल रहे हैं। हमारे पास nf-core और Nextflow से आने वाले community events हैं। हमारे पास मार्च में एक online और distributed nf-core hackathon है, पिछले साल हमारे पास दुनिया भर की sites के साथ एक हज़ार से अधिक लोग शामिल हुए थे। तो कृपया अगर तुम कर सकते हो तो हमारे साथ शामिल हो।

और हमारे पास Nextflow Summit events भी हैं, एक Boston में, और फिर हमारे पास Barcelona में और online एक event है। शानदार talks जहाँ तुम लोगों के बारे में सुन सकते हो जो Nextflow का उपयोग वास्तव में बड़े और wild और रोमांचक अलग-अलग तरीकों से कर रहे हैं। और उनके साथ जुड़े hackathons और in-person training भी हैं।

## Nextflow Podcast और blog

अगर तुम Nextflow ecosystem में चल रही चीज़ों के साथ अपडेट रहना चाहते हो, तो seqera.io/blog ज़रूर देखो।

वहाँ Nextflow के लिए एक section है जहाँ तुम community में काम करने वाले लोगों से community blog posts सुन सकते हो, और Seqera से Nextflow और अन्य tools के updates के बारे में blog posts भी जो हम generate करते हैं।

मैं अपने pet project के लिए भी एक plug देना चाहूँगा, जो Nextflow Podcast है। इसे Spotify, या Apple Music, या YouTube पर देखो। हम समय-समय पर नए episodes डालते हैं जहाँ मैं अन्य लोगों से बात करता हूँ, या तो Nextflow या associated technologies के साथ काम करने वाले, या community में लोग। और हम वास्तविक technical, deep dives करते हैं कि चीज़ें कैसे काम करती हैं और लोग क्या कर रहे हैं। तो अगर तुम रुचि रखते हो, तो उन्हें देखो। वे वास्तव में मज़ेदार हैं।

## धन्यवाद

ठीक है, मैं धन्यवाद का एक सेट करना चाहूँगा। Seqera की training team इस सामग्री के लिए जिम्मेदार है। मैं camera के सामने बैठा हूँ, लेकिन वास्तव में सारी मेहनत उन अन्य लोगों ने की है। Geraldine के लिए विशेष shout out, जिन्होंने Hello Nextflow और अन्य के लिए इस प्रशिक्षण सामग्री कोर्स को लिखा है और refresh किया है। और Jon भी, जिन्होंने वास्तव में मदद की है, विशेष रूप से नए Nextflow syntax के लिए syntax अपडेट करने में और कई कोर्स खुद लिखने में भी। Scientific development team में अन्य जैसे Rike, Rob, Florian, और कई अन्य लोगों का उस सामग्री में बहुत बड़ा योगदान रहा है जिसके साथ हम काम कर रहे हैं।

मैं community में लोगों को भी धन्यवाद देना चाहूँगा। नए translations, उदाहरण के लिए, जो बहुत नए हैं, ambassador program और अन्य जगहों के लोगों से बहुत प्रभावित हुए हैं। और वास्तव में, training material की open source प्रकृति का मतलब है कि हमारे पास pull requests और issues काफी बार आते हैं, जो वास्तव में हमारी मदद करते हैं।

## सर्वेक्षण

अब जब तुमने समाप्त कर लिया है, अगर तुमने पहले से नहीं किया है, तो कृपया जल्दी से feedback survey करो। यह training.nextflow.io website पर Hello Nextflow section के ठीक नीचे है।

यह केवल पाँच प्रश्न हैं। यह वास्तव में, वास्तव में तेज़ है, लेकिन यह हमें मोटे तौर पर track करने की अनुमति देता है कि कितने लोग प्रशिक्षण कर रहे हैं और तुम हमें यह भी बता सकते हो कि प्रशिक्षण सामग्री को कैसे सुधारा जाए। हम वास्तव में सभी responses की जाँच करते हैं, तो हम वास्तव में वहाँ तुम्हारे input को महत्व देते हैं।

## अलविदा

एक बार फिर, इस कोर्स में, और इस यात्रा में हमारे साथ शामिल होने के लिए बहुत-बहुत धन्यवाद। एक GitHub issue या Pull Request drop करो अगर तुमने प्रशिक्षण सामग्री में कुछ भी देखा, जिसे तुम्हें लगता है कि सुधारा जा सकता है। और मैं वास्तव में उम्मीद करता हूँ कि तुम्हें किसी अन्य Nextflow प्रशिक्षण कोर्स में, या शायद किसी hackathon या event में देखूँगा। फिर से धन्यवाद।​
