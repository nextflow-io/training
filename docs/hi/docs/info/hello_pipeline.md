---
title: Hello pipeline
description: Hello pipeline क्या करती है और यह कैसे संरचित है, इसका पुनरावलोकन।
hide:
  - toc
  - footer
---

# Hello pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

हमारे अधिकांश प्रशिक्षण पाठ्यक्रम Nextflow अवधारणाओं और तंत्रों को प्रदर्शित करने के लिए एक सरल डोमेन-अज्ञेयवादी pipeline का उपयोग करते हैं।
Hello Nextflow पाठ्यक्रम दिखाता है कि इस pipeline को चरण-दर-चरण तरीके से कैसे विकसित किया जाए जो प्रत्येक डिज़ाइन और कार्यान्वयन निर्णय को समझाता है।
अन्य प्रशिक्षण इस pipeline, या इसके भागों को, प्रारंभिक बिंदु के रूप में उपयोग करते हैं।

यह पृष्ठ Hello Nextflow पाठ्यक्रम के पूरा होने पर pipeline की स्थिति का सारांश देता है।

### सारांश विवरण

Hello workflow एक CSV फ़ाइल लेता है जिसमें अभिवादन होते हैं, उन्हें अलग-अलग फ़ाइलों में लिखता है, प्रत्येक को uppercase में परिवर्तित करता है, उन्हें वापस एकत्र करता है और अभिवादन कहते हुए एक मज़ेदार चरित्र की ASCII तस्वीर वाली एक टेक्स्ट फ़ाइल आउटपुट करता है।

### Workflow चरण (processes)

चार चरणों को Nextflow processes (`sayHello`, `convertToUpper`, `collectGreetings`, और `cowpy`) के रूप में लागू किया गया है जो अलग मॉड्यूल फ़ाइलों में संग्रहीत हैं।

1. **`sayHello`:** प्रत्येक अभिवादन को उसकी अपनी आउटपुट फ़ाइल में लिखता है (जैसे, "Hello-output.txt")
2. **`convertToUpper`:** प्रत्येक अभिवादन को uppercase में परिवर्तित करता है (जैसे, "HELLO")
3. **`collectGreetings`:** सभी uppercase अभिवादनों को एक एकल batch फ़ाइल में एकत्र करता है
4. **`cowpy`:** `cowpy` टूल का उपयोग करके ASCII art उत्पन्न करता है

### आरेख

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### परिणाम

परिणाम `results/` नामक डायरेक्टरी में प्रकाशित किए जाते हैं, और pipeline का अंतिम आउटपुट (जब डिफ़ॉल्ट पैरामीटर के साथ चलाया जाता है) एक plain text फ़ाइल है जिसमें uppercase अभिवादन कहते हुए एक turkey की ASCII art है।

```txt title="results/cowpy-COLLECTED-test-batch-output.txt"
  _________
/ BONJOUR \
| HELLO   |
\ HOLà    /
---------
  \                                  ,+*^^*+___+++_
  \                           ,*^^^^              )
    \                       _+*                     ^**+_
    \                    +^       _ _++*+_+++_,         )
              _+^^*+_    (     ,+*^ ^          \+_        )
            {       )  (    ,(    ,_+--+--,      ^)      ^\
            { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
          {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
          ( /  (    (        ,___    ^*+_+* )   <    <      \
          U _/     )    *--<  ) ^\-----++__)   )    )       )
            (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
          (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
        (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
          *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
          \             \_)^)_)) ))^^^^^^^^^^))^^^^)
          (_             ^\__^^^^^^^^^^^^))^^^^^^^)
            ^\___            ^\__^^^^^^))^^^^^^^^)\\
                  ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                    ___) >____) >___   ^\_\_\_\_\_\_\)
                    ^^^//\\_^^//\\_^       ^(\_\_\_\)
                      ^^^ ^^ ^^^ ^
```

पाठ्यक्रम के आधार पर जिसमें pipeline प्रदर्शित है, तुम्हें विशिष्टताओं में कुछ भिन्नताएँ मिल सकती हैं।

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
