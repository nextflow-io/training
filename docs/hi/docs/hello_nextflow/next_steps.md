# कोर्स सारांश

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello Nextflow training course complete करने पर बधाई! 🎉

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/xHOcx_4Ancg?si=Lp8hS8RdaMwbp5j5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [Nextflow YouTube channel पर पूरी playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) देखें।

:green_book: तुम video के साथ [video transcript](./transcripts/07_next_steps.md) पढ़ सकते हो।
///
-->

## तुम्हारी यात्रा

तुमने एक बहुत ही basic workflow से शुरू किया जो hardcoded command run करता था।
छह parts के दौरान, तुमने उस basic workflow को एक modular multi-step pipeline में transform किया जो Nextflow की key features exercise करता है जिसमें channels, operators, containers के लिए built-in support, और configuration options शामिल हैं।

### तुमने क्या बनाया

- Hello workflow का final form input के रूप में text greetings वाली CSV file लेता है।
- चार steps Nextflow processes (`sayHello`, `convertToUpper`, `collectGreetings`, और `cowpy`) के रूप में implement किए गए हैं जो separate module files में stored हैं।
- Results `results/` नामक directory में publish होते हैं।
- Pipeline का final output एक plain text file है जिसमें uppercased greetings बोलने वाले character की ASCII art है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** प्रत्येक greeting को उसकी अपनी output file में लिखता है (_जैसे_ "Hello-output.txt")
2. **`convertToUpper`:** प्रत्येक greeting को uppercase में convert करता है (_जैसे_ "HELLO")
3. **`collectGreetings`:** सभी uppercase greetings को एक single batch file में collect करता है
4. **`cowpy`:** `cowpy` tool का उपयोग करके ASCII art generate करता है

Workflow configuration flexible, reproducible तरीके से inputs और parameters provide करने को support करता है।

### प्राप्त कौशल

इस hands-on course के माध्यम से, तुमने सीखा कि कैसे:

- एक simple multi-step workflow बनाने के लिए पर्याप्त core Nextflow components describe और utilize करें
- Next-step concepts जैसे operators और channel factories describe करें
- Locally Nextflow workflow launch करें
- Nextflow द्वारा generate किए गए outputs (results) और log files खोजें और interpret करें
- Basic issues troubleshoot करें

तुम अब Nextflow में अपने खुद के pipelines develop करना शुरू करने के foundational knowledge से equipped हो।

## अपने skills बनाने के लिए next steps

यहाँ आगे क्या करना है इसके लिए हमारे top 3 suggestions हैं:

- [Nextflow for Science](../nf4_science/index.md) के साथ scientific analysis use case पर Nextflow apply करें
- [Hello nf-core](../hello_nf-core/index.md) के साथ nf-core शुरू करें
- [Side Quests](../side_quests/index.md) के साथ more advanced Nextflow features explore करें

Finally, हम recommend करते हैं कि तुम [**Seqera Platform**](https://seqera.io/) पर नज़र डालो, Nextflow के creators द्वारा develop किया गया एक cloud-based platform जो तुम्हारे workflows launch और manage करना, साथ ही तुम्हारा data manage करना और किसी भी environment में interactively analyses run करना और भी आसान बनाता है।

## Feedback survey

Move on करने से पहले, कृपया course survey complete करने के लिए एक minute लो! तुम्हारी feedback हमें सभी के लिए हमारी training materials improve करने में मदद करती है।

[Survey लें :material-arrow-right:](survey.md){ .md-button .md-button--primary }
