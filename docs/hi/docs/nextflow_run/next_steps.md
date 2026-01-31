# कोर्स सारांश

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow Run प्रशिक्षण कोर्स पूरा करने पर बधाई हो! 🎉

<!-- placeholder for video -->

## तुम्हारी यात्रा

तुमने एक बहुत ही basic workflow से शुरू किया, और इसे run करना, outputs खोजना, और इसके execution को manage करना सीखा।
फिर, तुमने उस workflow के increasingly अधिक complex versions के माध्यम से काम किया और essential concepts और mechanisms को पहचानना सीखा जो Nextflow pipelines को power करते हैं, जिसमें channels और operators, code modularization, और containers शामिल हैं।
अंत में, तुमने सीखा कि अपनी preferences और अपने computational infrastructure में fit करने के लिए pipeline की configuration को कैसे customize करें।

### तुमने क्या सीखा

तुम अब Hello pipeline के execution को manage करने, describe करने कि यह कैसे structured है, और involved code के main pieces identify करने में सक्षम हो।

- Hello workflow का final form input के रूप में text greetings containing एक CSV फ़ाइल लेता है।
- चार steps Nextflow processes (`sayHello`, `convertToUpper`, `collectGreetings`, और `cowpy`) के रूप में implemented हैं जो अलग module files में stored हैं।
- Results `results/` नाम की डायरेक्टरी में publish होते हैं।
- Pipeline का final output एक plain text फ़ाइल है जिसमें uppercased greetings बोलते हुए character की ASCII art है।

<figure class="excalidraw">
--8<-- "docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** प्रत्येक greeting को अपनी own output फ़ाइल में लिखता है (_e.g._ "Hello-output.txt")
2. **`convertToUpper`:** प्रत्येक greeting को uppercase में convert करता है (_e.g._ "HELLO")
3. **`collectGreetings`:** सभी uppercase greetings को एक single batch फ़ाइल में collect करता है
4. **`cowpy`:** `cowpy` tool का उपयोग करके ASCII art generate करता है

Workflow configuration flexible, reproducible तरीके से inputs और parameters provide करने का support करती है।

### अर्जित कौशल

इस hands-on course के माध्यम से, तुमने सीखा कि कैसे:

- Locally Nextflow workflow launch करें
- Nextflow द्वारा generated outputs (results) और log files खोजें और interpret करें
- एक simple multi-step workflow बनाने वाले core Nextflow components recognize करें
- Operators और channel factories जैसे next-step concepts describe करें
- Different computing environments के लिए pipelines configure करें

तुम अब अपने own work में existing Nextflow pipelines integrate करना शुरू करने के लिए foundational knowledge से equipped हो।

## अपने skills build करने के लिए next steps

यहां आगे क्या करना है इसके लिए हमारे top suggestions हैं:

- बस Nextflow run मत करो, इसे लिखो! [Hello Nextflow](../hello_nextflow/index.md) के साथ Nextflow developer बनो
- [Nextflow for Science](../nf4_science/index.md) के साथ scientific analysis use case पर Nextflow apply करो
- [Hello nf-core](../hello_nf-core/index.md) के साथ nf-core के साथ शुरू करो
- [Debugging Side Quest](../side_quests/debugging.md) के साथ troubleshooting techniques सीखो

अंत में, हम recommend करते हैं कि तुम [**Seqera Platform**](https://seqera.io/) पर एक नज़र डालो, Nextflow के creators द्वारा developed एक cloud-based platform जो तुम्हारी workflows launch और manage करना और भी आसान बनाता है, साथ ही तुम्हारे data manage करना और किसी भी environment में interactively analyses run करना।

## Help प्राप्त करना

Help resources और community support के लिए, [Help page](../help.md) देखो।

## Feedback survey

आगे बढ़ने से पहले, कृपया course survey complete करने के लिए एक मिनट लो! तुम्हारी feedback हमें सभी के लिए हमारी training materials improve करने में help करती है।

[Survey लें :material-arrow-right:](survey.md){ .md-button .md-button--primary }
