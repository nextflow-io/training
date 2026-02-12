# कोर्स सारांश

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow For RNAseq प्रशिक्षण कोर्स पूरा करने पर बधाई!

## तुम्हारी यात्रा

तुमने RNAseq प्रोसेसिंग टूल्स को टर्मिनल में मैन्युअली चलाकर शुरुआत की ताकि methodology को समझ सको।
फिर तुमने प्रोसेस को automate करने के लिए एक single-sample Nextflow pipeline बनाई, इसे parallel में कई samples को हैंडल करने के लिए स्केल किया, और इसे paired-end डेटा को हैंडल करने और samples में QC reports को aggregate करने के लिए विस्तारित किया।

### तुमने क्या बनाया

- एक RNAseq प्रोसेसिंग pipeline जो FASTQ फ़ाइलों को input के रूप में लेती है और trimmed reads, alignments और aggregated QC reports को output के रूप में उत्पन्न करती है।
- Trimming (Trim Galore), alignment (HISAT2), quality control (FastQC) और report aggregation (MultiQC) के लिए processes जो अलग module फ़ाइलों में संग्रहीत हैं।
- Pipeline Nextflow के dataflow paradigm का उपयोग करके input samples की प्रोसेसिंग को स्वचालित रूप से parallelize करती है।
- अंतिम pipeline paired-end sequencing डेटा को हैंडल करती है।

### प्राप्त कौशल

इस hands-on कोर्स के माध्यम से, तुमने सीखा है कि कैसे:

- बुनियादी RNAseq प्रोसेसिंग और QC methods को लागू करने के लिए एक linear workflow लिखें
- FASTQ और reference genome resources जैसी domain-specific फ़ाइलों को उचित रूप से हैंडल करें
- Single-end और paired-end sequencing डेटा को हैंडल करें
- Per-sample RNAseq प्रोसेसिंग को parallelize करने के लिए Nextflow के dataflow paradigm का लाभ उठाएं
- संबंधित channel operators का उपयोग करके कई steps और samples में QC reports को aggregate करें

अब तुम अपने खुद के काम में RNAseq विश्लेषण workflows में Nextflow को लागू करना शुरू करने के लिए तैयार हो।

## अपने कौशल को बढ़ाने के लिए अगले कदम

यहां हमारी शीर्ष सिफारिशें हैं कि आगे क्या करना है:

- [Nextflow for Science](../index.md) के साथ अन्य वैज्ञानिक विश्लेषण उपयोग मामलों में Nextflow लागू करें
- [Hello nf-core](../../hello_nf-core/index.md) के साथ nf-core के साथ शुरुआत करें
- [Side Quests](../../side_quests/index.md) के साथ अधिक उन्नत Nextflow सुविधाओं का अन्वेषण करें

अंत में, हम अनुशंसा करते हैं कि तुम [**Seqera Platform**](https://seqera.io/) पर एक नज़र डालो, जो Nextflow के निर्माताओं द्वारा विकसित एक क्लाउड-आधारित प्लेटफ़ॉर्म है जो तुम्हारे workflows को लॉन्च और प्रबंधित करना, साथ ही तुम्हारे डेटा को प्रबंधित करना और किसी भी वातावरण में इंटरैक्टिव रूप से विश्लेषण चलाना और भी आसान बनाता है।

## सहायता प्राप्त करना

सहायता संसाधनों और community support के लिए, [Help पेज](../../help.md) देखें।

## फीडबैक सर्वेक्षण

आगे बढ़ने से पहले, कृपया कोर्स सर्वेक्षण को पूरा करने के लिए एक मिनट का समय निकालें! तुम्हारा फीडबैक हमें सभी के लिए अपनी प्रशिक्षण सामग्री को बेहतर बनाने में मदद करता है।

[सर्वेक्षण लें :material-arrow-right:](survey.md){ .md-button .md-button--primary }
