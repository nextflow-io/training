# कोर्स सारांश

Nextflow Run प्रशिक्षण कोर्स पूरा करने पर बधाई! 🎉

<!-- placeholder for video -->

## तुम्हारी यात्रा

तुमने एक बहुत ही बुनियादी workflow से शुरुआत की, और इसे चलाना, outputs खोजना, और इसके execution को manage करना सीखा।
फिर, तुमने उस workflow के क्रमशः अधिक जटिल versions पर काम किया और उन आवश्यक concepts और mechanisms को पहचानना सीखा जो Nextflow pipelines को शक्ति देते हैं, जिनमें channels और operators, code modularization, और containers शामिल हैं।
अंत में, तुमने सीखा कि किसी pipeline के configuration को अपनी preferences और अपने computational infrastructure के अनुसार कैसे customize करें।

### तुमने क्या सीखा

अब तुम Hello pipeline के execution को manage करने, यह वर्णन करने में सक्षम हो कि यह कैसे structured है, और इसमें शामिल मुख्य code pieces को identify करने में सक्षम हो।

- Hello workflow का अंतिम रूप input के रूप में एक CSV फ़ाइल लेता है जिसमें text greetings होती हैं।
- चार steps को Nextflow processes (`sayHello`, `convertToUpper`, `collectGreetings`, और `cowpy`) के रूप में implement किया गया है जो अलग-अलग module फ़ाइलों में stored हैं।
- परिणाम `results/` नामक directory में publish किए जाते हैं।
- Pipeline का अंतिम output एक plain text फ़ाइल है जिसमें एक character की ASCII art होती है जो uppercased greetings बोल रहा है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** प्रत्येक greeting को अपनी output फ़ाइल में लिखता है (_जैसे_ "Hello-output.txt")
2. **`convertToUpper`:** प्रत्येक greeting को uppercase में convert करता है (_जैसे_ "HELLO")
3. **`collectGreetings`:** सभी uppercase greetings को एक single batch फ़ाइल में collect करता है
4. **`cowpy`:** `cowpy` tool का उपयोग करके ASCII art generate करता है

Workflow configuration inputs और parameters को flexible, reproducible तरीके से provide करने का समर्थन करता है।

### प्राप्त कौशल

इस hands-on कोर्स के माध्यम से, तुमने सीखा है कि कैसे:

- Nextflow workflow को locally launch करें
- Nextflow द्वारा generate की गई outputs (results) और log फ़ाइलों को खोजें और interpret करें
- एक simple multi-step workflow बनाने वाले core Nextflow components को पहचानें
- Operators और channel factories जैसे next-step concepts का वर्णन करें
- विभिन्न computing environments के लिए pipelines को configure करें

अब तुम अपने काम में मौजूदा Nextflow pipelines को integrate करना शुरू करने के लिए foundational knowledge से लैस हो।

## अपने कौशल को बढ़ाने के लिए अगले कदम

यहाँ हमारे top सुझाव हैं कि आगे क्या करना है:

- सिर्फ़ Nextflow चलाओ नहीं, इसे लिखो! [Hello Nextflow](../hello_nextflow/index.md) के साथ Nextflow developer बनो
- [Nextflow for Science](../nf4_science/index.md) के साथ scientific analysis use case में Nextflow को apply करो
- [Hello nf-core](../hello_nf-core/index.md) के साथ nf-core से शुरुआत करो
- [Debugging Side Quest](../side_quests/debugging.md) के साथ troubleshooting techniques सीखो

अंत में, हम recommend करते हैं कि तुम [**Seqera Platform**](https://seqera.io/) पर एक नज़र डालो, जो Nextflow के creators द्वारा विकसित एक cloud-based platform है जो तुम्हारे workflows को launch और manage करना और भी आसान बनाता है, साथ ही तुम्हारे data को manage करने और किसी भी environment में interactively analyses चलाने में मदद करता है।

## मदद प्राप्त करना

मदद के संसाधनों और community support के लिए, [Help page](../help.md) देखो।

## Feedback survey

आगे बढ़ने से पहले, कृपया कोर्स survey पूरा करने के लिए एक मिनट निकालो! तुम्हारा feedback हमें सभी के लिए अपनी प्रशिक्षण सामग्री को बेहतर बनाने में मदद करता है।

[Survey लो :material-arrow-right:](survey.md){ .md-button .md-button--primary }
