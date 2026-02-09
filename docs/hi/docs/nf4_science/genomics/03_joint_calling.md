# भाग 3: कोहोर्ट पर ज्वाइंट कॉलिंग

भाग 2 में, तुमने एक प्रति-नमूना वेरिएंट कॉलिंग पाइपलाइन बनाई जो प्रत्येक नमूने के डेटा को स्वतंत्र रूप से प्रोसेस करती थी।
अब हम इसे ज्वाइंट वेरिएंट कॉलिंग को लागू करने के लिए विस्तारित करने जा रहे हैं, जैसा कि [भाग 1](01_method.md) में बताया गया है।

## असाइनमेंट

इस कोर्स के इस भाग में, हम वर्कफ़्लो को निम्नलिखित करने के लिए विस्तारित करने जा रहे हैं:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Samtools का उपयोग करके प्रत्येक BAM इनपुट फ़ाइल के लिए एक इंडेक्स फ़ाइल जेनरेट करें
2. प्रत्येक BAM इनपुट फ़ाइल पर GATK HaplotypeCaller चलाएं ताकि प्रति-नमूना जीनोमिक वेरिएंट कॉल्स का एक GVCF जेनरेट हो सके
3. सभी GVCFs को इकट्ठा करें और उन्हें एक GenomicsDB डेटा स्टोर में संयोजित करें
4. कोहोर्ट-स्तर VCF उत्पन्न करने के लिए संयोजित GVCF डेटा स्टोर पर ज्वाइंट जीनोटाइपिंग चलाएं

यह भाग सीधे भाग 2 द्वारा उत्पादित वर्कफ़्लो पर आधारित है।

??? info "इस सेक्शन से कैसे शुरू करें"

    कोर्स का यह सेक्शन मानता है कि तुमने [भाग 2: प्रति-नमूना वेरिएंट कॉलिंग](./02_per_sample_variant_calling.md) पूरा कर लिया है और तुम्हारे पास एक कार्यशील `genomics.nf` पाइपलाइन है।

    यदि तुमने भाग 2 पूरा नहीं किया है या इस भाग के लिए नए सिरे से शुरू करना चाहते हो, तो तुम भाग 2 के समाधान को अपने शुरुआती बिंदु के रूप में उपयोग कर सकते हो।
    `nf4-science/genomics/` डायरेक्टरी के अंदर से ये कमांड चलाओ:

    ```bash
    cp solutions/part2/genomics-2.nf genomics.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    यह तुम्हें एक पूर्ण प्रति-नमूना वेरिएंट कॉलिंग वर्कफ़्लो देता है।
    तुम निम्नलिखित कमांड चलाकर परीक्षण कर सकते हो कि यह सफलतापूर्वक चलता है:

    ```bash
    nextflow run genomics.nf -profile test
    ```

## पाठ योजना

हमने इसे दो चरणों में विभाजित किया है:

1. **प्रति-नमूना वेरिएंट कॉलिंग चरण को GVCF उत्पन्न करने के लिए संशोधित करें।**
   यह प्रोसेस कमांड और आउटपुट को अपडेट करना कवर करता है।
2. **एक ज्वाइंट जीनोटाइपिंग चरण जोड़ें जो प्रति-नमूना GVCFs को संयोजित और जीनोटाइप करता है।**
   यह `collect()` ऑपरेटर, कमांड-लाइन निर्माण के लिए Groovy closures, और मल्टी-कमांड प्रोसेस का परिचय देता है।

!!! note

     सुनिश्चित करो कि तुम सही वर्किंग डायरेक्टरी में हो:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. प्रति-नमूना वेरिएंट कॉलिंग चरण को GVCF उत्पन्न करने के लिए संशोधित करें

भाग 2 की पाइपलाइन VCF फ़ाइलें उत्पन्न करती है, लेकिन ज्वाइंट कॉलिंग के लिए GVCF फ़ाइलों की आवश्यकता होती है।
हमें GVCF वेरिएंट कॉलिंग मोड को चालू करना होगा और आउटपुट फ़ाइल एक्सटेंशन को अपडेट करना होगा।

[भाग 1](01_method.md) से GVCF वेरिएंट कॉलिंग कमांड को याद करो:

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

भाग 2 में हमने जो बेस HaplotypeCaller कमांड रैप की थी, उसकी तुलना में अंतर `-ERC GVCF` पैरामीटर और `.g.vcf` आउटपुट एक्सटेंशन हैं।

### 1.1. HaplotypeCaller को GVCF emit करने के लिए कहें और आउटपुट एक्सटेंशन अपडेट करें

`modules/gatk_haplotypecaller.nf` मॉड्यूल फ़ाइल खोलो और दो बदलाव करो:

- GATK HaplotypeCaller कमांड में `-ERC GVCF` पैरामीटर जोड़ो;
- GATK कन्वेंशन के अनुसार, आउटपुट फ़ाइल पथ को संबंधित `.g.vcf` एक्सटेंशन का उपयोग करने के लिए अपडेट करो।

सुनिश्चित करो कि जब तुम `-ERC GVCF` जोड़ो तो पिछली लाइन के अंत में बैकस्लैश (`\`) जोड़ो।

=== "बाद में"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5 7"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "पहले"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

हमें नए फ़ाइल एक्सटेंशन से मेल खाने के लिए आउटपुट ब्लॉक को भी अपडेट करना होगा।
चूंकि हमने कमांड आउटपुट को `.vcf` से `.g.vcf` में बदल दिया है, प्रोसेस `output:` ब्लॉक को भी वही बदलाव दर्शाना चाहिए।

### 1.2. प्रोसेस आउटपुट ब्लॉक में आउटपुट फ़ाइल एक्सटेंशन अपडेट करें

=== "बाद में"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "पहले"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

हमें नए GVCF आउटपुट को दर्शाने के लिए वर्कफ़्लो के publish और output कॉन्फ़िगरेशन को भी अपडेट करना होगा।

### 1.3. नए GVCF आउटपुट के लिए publish टारगेट अपडेट करें

चूंकि अब हम VCFs के बजाय GVCFs उत्पन्न कर रहे हैं, हमें वर्कफ़्लो के `publish:` सेक्शन को अधिक वर्णनात्मक नामों का उपयोग करने के लिए अपडेट करना चाहिए।
हम स्पष्टता के लिए GVCF फ़ाइलों को उनकी अपनी सबडायरेक्टरी में भी व्यवस्थित करेंगे।

=== "बाद में"

    ```groovy title="genomics.nf" linenums="47" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="47"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

अब output ब्लॉक को मेल खाने के लिए अपडेट करो।

### 1.4. नई डायरेक्टरी संरचना के लिए output ब्लॉक अपडेट करें

हमें GVCF फ़ाइलों को `gvcf` सबडायरेक्टरी में रखने के लिए `output` ब्लॉक को भी अपडेट करना होगा।

=== "बाद में"

    ```groovy title="genomics.nf" linenums="53" hl_lines="3 5 6 8 9"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="53"
    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

मॉड्यूल, publish टारगेट, और output ब्लॉक सभी अपडेट होने के साथ, हम बदलावों का परीक्षण कर सकते हैं।

### 1.5. पाइपलाइन चलाओ

बदलावों को सत्यापित करने के लिए वर्कफ़्लो चलाओ।

```bash
nextflow run genomics.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (6)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Nextflow आउटपुट पहले जैसा ही दिखता है, लेकिन `.g.vcf` फ़ाइलें और उनकी इंडेक्स फ़ाइलें अब सबडायरेक्टरी में व्यवस्थित हैं।

??? abstract "डायरेक्टरी सामग्री (symlinks छोटे किए गए)"

    ```console
    results/
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

यदि तुम GVCF फ़ाइलों में से एक को खोलते हो और इसे स्क्रॉल करते हो, तो तुम सत्यापित कर सकते हो कि GATK HaplotypeCaller ने अनुरोध के अनुसार GVCF फ़ाइलें उत्पन्न कीं।

### सारांश

जब तुम किसी टूल कमांड की आउटपुट फ़ाइलनेम बदलते हो, तो प्रोसेस `output:` ब्लॉक और publish/output कॉन्फ़िगरेशन को मेल खाने के लिए अपडेट किया जाना चाहिए।

### आगे क्या है?

चैनल की सामग्री को इकट्ठा करना और उन्हें एकल इनपुट के रूप में अगली प्रोसेस में पास करना सीखो।

---

## 2. ज्वाइंट जीनोटाइपिंग चरण जोड़ें

अब हमें प्रति-नमूना GVCFs को इकट्ठा करना, उन्हें GenomicsDB डेटा स्टोर में संयोजित करना, और कोहोर्ट-स्तर VCF उत्पन्न करने के लिए ज्वाइंट जीनोटाइपिंग चलाना होगा।
जैसा कि [भाग 1](01_method.md) में बताया गया है, यह दो-टूल ऑपरेशन है: GenomicsDBImport GVCFs को संयोजित करता है, फिर GenotypeGVCFs अंतिम वेरिएंट कॉल्स उत्पन्न करता है।
हम दोनों टूल्स को `GATK_JOINTGENOTYPING` नामक एकल प्रोसेस में रैप करेंगे।

[भाग 1](01_method.md) से दोनों कमांड को याद करो:

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

पहली कमांड प्रति-नमूना GVCFs और एक intervals फ़ाइल लेती है, और एक GenomicsDB डेटा स्टोर उत्पन्न करती है।
दूसरी उस डेटा स्टोर, एक रेफरेंस जीनोम लेती है, और अंतिम कोहोर्ट-स्तर VCF उत्पन्न करती है।
कंटेनर URI HaplotypeCaller के समान है: `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`।

### 2.1. इनपुट सेट अप करें

ज्वाइंट जीनोटाइपिंग प्रोसेस को दो प्रकार के इनपुट की आवश्यकता है जो हमारे पास अभी तक नहीं हैं: एक मनमाना कोहोर्ट नाम, और सभी नमूनों से एकत्रित GVCF आउटपुट एक साथ बंडल किए गए।

#### 2.1.1. एक `cohort_name` पैरामीटर जोड़ें

हमें कोहोर्ट के लिए एक मनमाना नाम प्रदान करना होगा।
ट्रेनिंग सीरीज़ में बाद में तुम इस तरह की चीज़ों के लिए नमूना मेटाडेटा का उपयोग करना सीखोगे, लेकिन अभी के लिए हम सिर्फ `params` का उपयोग करके एक CLI पैरामीटर घोषित करते हैं और सुविधा के लिए इसे एक डिफ़ॉल्ट मान देते हैं।

=== "बाद में"

    ```groovy title="genomics.nf" linenums="14" hl_lines="3-4"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"

        // Base name for final output file
        cohort_name: String = "family_trio"
    }
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="14"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"
    }
    ```

#### 2.1.2. नमूनों में HaplotypeCaller आउटपुट इकट्ठा करें

यदि हम `GATK_HAPLOTYPECALLER` से आउटपुट चैनल को सीधे नई प्रोसेस में प्लग करते, तो Nextflow प्रत्येक नमूना GVCF पर अलग से प्रोसेस को कॉल करता।
हम सभी तीन GVCFs (और उनकी इंडेक्स फ़ाइलों) को बंडल करना चाहते हैं ताकि Nextflow उन सभी को एक साथ एकल प्रोसेस कॉल में सौंप दे।

हम `collect()` चैनल ऑपरेटर का उपयोग करके ऐसा कर सकते हैं।
GATK_HAPLOTYPECALLER की कॉल के ठीक बाद, `workflow` बॉडी में निम्नलिखित लाइनें जोड़ो:

=== "बाद में"

    ```groovy title="genomics.nf" hl_lines="4-6"
            intervals_file
        )

        // Collect variant calling outputs across samples
        all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

=== "पहले"

    ```groovy title="genomics.nf"
            intervals_file
        )
    ```

इसे तोड़ते हुए:

1. हम `.out` प्रॉपर्टी का उपयोग करके `GATK_HAPLOTYPECALLER` से आउटपुट चैनल लेते हैं।
2. क्योंकि हमने सेक्शन 1 में `emit:` का उपयोग करके आउटपुट को नाम दिया, हम `.vcf` के साथ GVCFs और `.idx` के साथ इंडेक्स फ़ाइलों का चयन कर सकते हैं। नामित आउटपुट के बिना, हमें `.out[0]` और `.out[1]` का उपयोग करना होगा।
3. `collect()` ऑपरेटर सभी फ़ाइलों को एकल एलिमेंट में बंडल करता है, इसलिए `all_gvcfs_ch` में सभी तीन GVCFs एक साथ हैं, और `all_idxs_ch` में सभी तीन इंडेक्स फ़ाइलें एक साथ हैं।

हम GVCFs और उनकी इंडेक्स फ़ाइलों को अलग से इकट्ठा कर सकते हैं (उन्हें tuples में एक साथ रखने के विपरीत) क्योंकि Nextflow निष्पादन के लिए सभी इनपुट फ़ाइलों को एक साथ stage करेगा, इसलिए इंडेक्स फ़ाइलें GVCFs के साथ मौजूद होंगी।

!!! tip

    तुम चैनल ऑपरेटर लागू करने से पहले और बाद में चैनलों की सामग्री का निरीक्षण करने के लिए `view()` ऑपरेटर का उपयोग कर सकते हो।

### 2.2. ज्वाइंट जीनोटाइपिंग प्रोसेस लिखें और इसे वर्कफ़्लो में कॉल करें

भाग 2 में हमने जो पैटर्न उपयोग किया था उसी का पालन करते हुए, हम मॉड्यूल फ़ाइल में प्रोसेस डेफिनिशन लिखेंगे, इसे वर्कफ़्लो में import करेंगे, और इसे हमने अभी तैयार किए गए इनपुट पर कॉल करेंगे।

#### 2.2.1. प्रत्येक GVCF को `-V` आर्गुमेंट देने के लिए एक स्ट्रिंग बनाएं

प्रोसेस डेफिनिशन भरना शुरू करने से पहले, एक चीज़ का पता लगाना है।
GenomicsDBImport कमांड प्रत्येक GVCF फ़ाइल के लिए एक अलग `-V` आर्गुमेंट की अपेक्षा करती है, इस तरह:

```bash
gatk GenomicsDBImport \
    -V reads_mother.bam.g.vcf \
    -V reads_father.bam.g.vcf \
    -V reads_son.bam.g.vcf \
    ...
```

यदि हम `-V ${all_gvcfs_ch}` लिखते, तो Nextflow बस फ़ाइलनेम को concatenate करता और कमांड का वह भाग इस तरह दिखता:

```groovy
-V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf
```

लेकिन हमें स्ट्रिंग को इस तरह दिखना चाहिए:

```groovy
-V reads_mother.bam.g.vcf -V reads_father.bam.g.vcf -V reads_son.bam.g.vcf
```

महत्वपूर्ण रूप से, हमें इस स्ट्रिंग को एकत्रित चैनल में जो भी फ़ाइलें हैं उनसे गतिशील रूप से बनाना होगा।
Nextflow (Groovy के माध्यम से) ऐसा करने का एक संक्षिप्त तरीका प्रदान करता है:

```groovy
def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
```

इसे तोड़ते हुए:

1. `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }` प्रत्येक फ़ाइल पथ पर iterate करता है और इसमें `-V ` prepend करता है, `["-V A.g.vcf", "-V B.g.vcf", "-V C.g.vcf"]` उत्पन्न करता है।
2. `.join(' ')` उन्हें स्पेस के साथ concatenate करता है: `"-V A.g.vcf -V B.g.vcf -V C.g.vcf"`।
3. परिणाम एक लोकल वेरिएबल `gvcfs_line` (जो `def` के साथ परिभाषित है) को असाइन किया जाता है, जिसे हम कमांड टेम्पलेट में interpolate कर सकते हैं।

यह लाइन प्रोसेस के `script:` ब्लॉक के अंदर, कमांड टेम्पलेट से पहले जाती है।
तुम `script:` और कमांड टेम्पलेट के शुरुआती `"""` के बीच मनमाना Groovy कोड रख सकते हो।

फिर तुम प्रोसेस के `script:` ब्लॉक में उस पूरी स्ट्रिंग को `gvcfs_line` के रूप में संदर्भित कर सकोगे।

#### 2.2.2. ज्वाइंट जीनोटाइपिंग प्रोसेस के लिए मॉड्यूल भरें

अब हम पूर्ण प्रोसेस लिखने का काम कर सकते हैं।

`modules/gatk_jointgenotyping.nf` खोलो और प्रोसेस डेफिनिशन की रूपरेखा की जांच करो।

आगे बढ़ो और ऊपर दी गई जानकारी का उपयोग करके प्रोसेस डेफिनिशन भरो, फिर नीचे "बाद में" टैब में समाधान के विरुद्ध अपने काम की जांच करो।

=== "पहले"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * GVCFs को GenomicsDB datastore में संयोजित करें और कोहोर्ट-स्तर कॉल्स उत्पन्न करने के लिए ज्वाइंट जीनोटाइपिंग चलाएं
     */
    process GATK_JOINTGENOTYPING {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "बाद में"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1" hl_lines="8 11-17 20-21 24-25 29-33"
    #!/usr/bin/env nextflow

    /*
     * GVCFs को GenomicsDB datastore में संयोजित करें और कोहोर्ट-स्तर कॉल्स उत्पन्न करने के लिए ज्वाइंट जीनोटाइपिंग चलाएं
     */
    process GATK_JOINTGENOTYPING {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx

        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb

        gatk GenotypeGVCFs \
            -R ${ref_fasta} \
            -V gendb://${cohort_name}_gdb \
            -L ${interval_list} \
            -O ${cohort_name}.joint.vcf
        """
    }
    ```

यहां कई चीजें ध्यान देने योग्य हैं।

पहले की तरह, कई इनपुट सूचीबद्ध हैं भले ही कमांड उन्हें सीधे संदर्भित नहीं करती: `all_idxs`, `ref_index`, और `ref_dict`।
उन्हें सूचीबद्ध करना सुनिश्चित करता है कि Nextflow इन फ़ाइलों को वर्किंग डायरेक्टरी में उन फ़ाइलों के साथ stage करता है जो कमांड में दिखाई देती हैं, जिन्हें GATK नामकरण कन्वेंशन के आधार पर खोजने की अपेक्षा करता है।

`gvcfs_line` वेरिएबल GenomicsDBImport के लिए `-V` आर्गुमेंट बनाने के लिए ऊपर वर्णित Groovy closure का उपयोग करता है।

यह प्रोसेस दो कमांड को सीरियल में चलाती है, जैसे तुम टर्मिनल में करोगे।
GenomicsDBImport प्रति-नमूना GVCFs को एक डेटा स्टोर में संयोजित करता है, फिर GenotypeGVCFs उस डेटा स्टोर को पढ़ता है और अंतिम कोहोर्ट-स्तर VCF उत्पन्न करता है।
GenomicsDB डेटा स्टोर (`${cohort_name}_gdb`) एक मध्यवर्ती artifact है जिसका उपयोग केवल प्रोसेस के भीतर किया जाता है; यह output ब्लॉक में दिखाई नहीं देता।

एक बार जब तुम यह पूरा कर लो, तो प्रोसेस उपयोग के लिए तैयार है।
इसे वर्कफ़्लो में उपयोग करने के लिए, तुम्हें मॉड्यूल को import करना होगा और एक प्रोसेस कॉल जोड़नी होगी।

#### 2.2.3. मॉड्यूल import करें

`genomics.nf` में import स्टेटमेंट जोड़ो, मौजूदा import स्टेटमेंट के नीचे:

=== "बाद में"

    ```groovy title="genomics.nf" linenums="21" hl_lines="3"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    include { GATK_JOINTGENOTYPING } from './modules/gatk_jointgenotyping.nf'
    ```

=== "पहले"

    ```groovy title="genomics.nf" linenums="21"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

प्रोसेस अब वर्कफ़्लो स्कोप में उपलब्ध है।

#### 2.2.4. प्रोसेस कॉल जोड़ें

`collect()` लाइनों के बाद, वर्कफ़्लो बॉडी में `GATK_JOINTGENOTYPING` की कॉल जोड़ो:

=== "बाद में"

    ```groovy title="genomics.nf" hl_lines="3-12"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

        // Combine GVCFs into a GenomicsDB data store and apply joint genotyping
        GATK_JOINTGENOTYPING(
            all_gvcfs_ch,
            all_idxs_ch,
            intervals_file,
            params.cohort_name,
            ref_file,
            ref_index_file,
            ref_dict_file
        )
    ```

=== "पहले"

    ```groovy title="genomics.nf"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

प्रोसेस अब पूरी तरह से wired up है।
अगला, हम कॉन्फ़िगर करते हैं कि आउटपुट कैसे publish किए जाते हैं।

### 2.3. आउटपुट हैंडलिंग कॉन्फ़िगर करें

हमें ज्वाइंट VCF आउटपुट publish करने की आवश्यकता है।
ज्वाइंट जीनोटाइपिंग परिणामों के लिए publish टारगेट और output ब्लॉक एंट्री जोड़ें।

#### 2.3.1. ज्वाइंट VCF के लिए publish टारगेट जोड़ें

वर्कफ़्लो के `publish:` सेक्शन में ज्वाइंट VCF और इसका इंडेक्स जोड़ो:

=== "बाद में"

    ```groovy title="genomics.nf" hl_lines="5-6"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
        joint_vcf = GATK_JOINTGENOTYPING.out.vcf
        joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
    ```

=== "पहले"

    ```groovy title="genomics.nf"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

अब output ब्लॉक को मेल खाने के लिए अपडेट करो।

#### 2.3.2. ज्वाइंट VCF के लिए output ब्लॉक एंट्री जोड़ें

ज्वाइंट VCF फ़ाइलों के लिए एंट्री जोड़ो।
हम उन्हें results डायरेक्टरी के रूट में रखेंगे क्योंकि यह अंतिम आउटपुट है।

=== "बाद में"

    ```groovy title="genomics.nf" hl_lines="11-16"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
        joint_vcf {
            path '.'
        }
        joint_vcf_idx {
            path '.'
        }
    }
    ```

=== "पहले"

    ```groovy title="genomics.nf"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

प्रोसेस, publish टारगेट, और output ब्लॉक सभी जगह पर होने के साथ, हम पूर्ण वर्कफ़्लो का परीक्षण कर सकते हैं।

### 2.4. वर्कफ़्लो चलाओ

सब कुछ काम करता है यह सत्यापित करने के लिए वर्कफ़्लो चलाओ।

```bash
nextflow run genomics.nf -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

पहले दो चरण पिछले रन से cached हैं, और नया `GATK_JOINTGENOTYPING` चरण सभी तीन नमूनों से एकत्रित इनपुट पर एक बार चलता है।
अंतिम आउटपुट फ़ाइल, `family_trio.joint.vcf` (और इसका इंडेक्स), results डायरेक्टरी में है।

??? abstract "डायरेक्टरी सामग्री (symlinks छोटे किए गए)"

    ```console
    results/
    ├── family_trio.joint.vcf -> */a6/7cc8ed*/family_trio.joint.vcf
    ├── family_trio.joint.vcf.idx -> */a6/7cc8ed*/family_trio.joint.vcf.idx
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

यदि तुम ज्वाइंट VCF फ़ाइल खोलते हो, तो तुम सत्यापित कर सकते हो कि वर्कफ़्लो ने अपेक्षित वेरिएंट कॉल्स उत्पन्न कीं।

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

अब तुम्हारे पास एक स्वचालित, पूरी तरह से reproducible ज्वाइंट वेरिएंट कॉलिंग वर्कफ़्लो है!

!!! note

    ध्यान रखो कि हमने तुम्हें जो डेटा फ़ाइलें दीं वे क्रोमोसोम 20 के केवल एक छोटे से हिस्से को कवर करती हैं।
    वेरिएंट callset का वास्तविक आकार लाखों वेरिएंट में गिना जाएगा।
    इसीलिए हम प्रशिक्षण उद्देश्यों के लिए केवल डेटा के छोटे सबसेट का उपयोग करते हैं!

### सारांश

तुम जानते हो कि चैनल से आउटपुट कैसे इकट्ठा करें और उन्हें दूसरी प्रोसेस में एकल इनपुट के रूप में कैसे बंडल करें।
तुम यह भी जानते हो कि Groovy closures का उपयोग करके कमांड लाइन कैसे बनाएं, और एकल प्रोसेस में कई कमांड कैसे चलाएं।

### आगे क्या है?

अपनी सफलता का जश्न मनाओ और एक अच्छा ब्रेक लो।

इस कोर्स के अगले भाग में, तुम nf-core से एक production-ready वेरिएंट कॉलिंग पाइपलाइन चलाना सीखोगे और इसकी तुलना तुमने मैन्युअली बनाई गई पाइपलाइन से करोगे।
