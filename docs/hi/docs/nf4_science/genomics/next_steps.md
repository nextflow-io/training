# कोर्स सारांश

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow for Genomics प्रशिक्षण कोर्स पूरा करने पर बधाई! 🎉

## तुम्हारी यात्रा

तुमने methodology को समझने के लिए टर्मिनल में variant calling टूल्स को manually चलाकर शुरुआत की।
फिर तुमने प्रोसेस को automate करने के लिए एक single-sample Nextflow पाइपलाइन बनाई, इसे parallel में कई samples को handle करने के लिए scale किया, और channel operators का उपयोग करके multi-sample joint genotyping जोड़ा।

### तुमने क्या बनाया

- एक variant calling पाइपलाइन जो BAM फ़ाइलों को इनपुट के रूप में लेती है और joint-called VCFs को आउटपुट के रूप में produce करती है।
- तीन processes (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER`, और `GATK_JOINTGENOTYPING`) जो अलग-अलग module फ़ाइलों में stored हैं।
- पाइपलाइन Nextflow के dataflow paradigm का उपयोग करके इनपुट samples की processing को automatically parallelize करती है।
- परिणाम `results/` नामक डायरेक्टरी में publish किए जाते हैं।

### प्राप्त कौशल

इस hands-on कोर्स के माध्यम से, तुमने सीखा है कि कैसे:

- एक single sample पर variant calling लागू करने के लिए linear workflow लिखें
- Index फ़ाइलों और reference genome resources जैसी accessory फ़ाइलों को उचित रूप से handle करें
- Per-sample variant calling को parallelize करने के लिए Nextflow के dataflow paradigm का लाभ उठाएं
- Relevant channel operators का उपयोग करके multi-sample joint calling को implement करें

अब तुम अपने काम में genomics analysis workflows पर Nextflow लागू करना शुरू करने के लिए तैयार हो।

## अपने कौशल को बढ़ाने के लिए अगले कदम

यहाँ हमारे शीर्ष सुझाव हैं कि आगे क्या करें:

- अन्य scientific analysis use cases पर Nextflow लागू करें [Nextflow for Science](../index.md) के साथ
- nf-core के साथ शुरुआत करें [Hello nf-core](../../hello_nf-core/index.md) के साथ
- अधिक advanced Nextflow features का पता लगाएं [Side Quests](../../side_quests/index.md) के साथ

अंत में, हम recommend करते हैं कि तुम [**Seqera Platform**](https://seqera.io/) पर एक नज़र डालो, जो Nextflow के creators द्वारा विकसित एक cloud-based platform है जो तुम्हारे workflows को launch और manage करना, साथ ही तुम्हारे data को manage करना और किसी भी environment में interactively analyses चलाना और भी आसान बनाता है।

## मदद प्राप्त करना

मदद के संसाधनों और community support के लिए, [Help page](../../help.md) देखें।

## फीडबैक सर्वे

आगे बढ़ने से पहले, कृपया कोर्स सर्वे को पूरा करने में एक मिनट का समय दें! तुम्हारा फीडबैक हमें सभी के लिए अपनी प्रशिक्षण सामग्री को बेहतर बनाने में मदद करता है।

[सर्वे लें :material-arrow-right:](survey.md){ .md-button .md-button--primary }
