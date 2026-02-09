# भाग 3: कोहोर्ट पर जॉइंट कॉलिंग

भाग 2 में, तुमने एक per-sample वेरिएंट कॉलिंग पाइपलाइन बनाई थी जो हर सैंपल का डेटा स्वतंत्र रूप से प्रोसेस करती है।
अब हम इसे एक्सटेंड करके जॉइंट वेरिएंट कॉलिंग को इम्प्लीमेंट करने जा रहे हैं, जैसा कि [भाग 1](01_method.md) में बताया गया था।

## असाइनमेंट

कोर्स के इस भाग में, हम वर्कफ़्लो को निम्नलिखित कार्यों के लिए एक्सटेंड करने जा रहे हैं:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. हर BAM इनपुट फ़ाइल के लिए Samtools का उपयोग करके एक इंडेक्स फ़ाइल जेनरेट करो
2. हर BAM इनपुट फ़ाइल पर GATK HaplotypeCaller को रन करो ताकि per-sample जीनोमिक वेरिएंट कॉल्स का एक GVCF जेनरेट हो
3. सभी GVCFs को इकट्ठा करो और उन्हें एक GenomicsDB डेटा स्टोर में कंबाइन करो
4. कंबाइन किए गए GVCF डेटा स्टोर पर जॉइंट जेनोटाइपिंग रन करके कोहोर्ट-लेवल VCF बनाओ

यह भाग सीधे भाग 2 द्वारा बनाए गए वर्कफ़्लो पर आधारित है।

??? info "इस सेक्शन से कैसे शुरू करें"

    कोर्स का यह सेक्शन मानता है कि तुमने [भाग 2: Per-sample वेरिएंट कॉलिंग](./02_per_sample_variant_calling.md) पूरा कर लिया है और तुम्हारे पास एक काम करती हुई `genomics.nf` पाइपलाइन है।

    अगर तुमने भाग 2 पूरा नहीं किया है या इस भाग के लिए नए सिरे से शुरू करना चाहते हो, तो तुम भाग 2 के सॉल्यूशन को अपने शुरुआती बिंदु के रूप में उपयोग कर सकते हो।
    `nf4-science/genomics/` डायरेक्टरी के अंदर से ये कमांड्स रन करो:

    ```bash
    cp solutions/part2/genomics-2.nf genomics.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    इससे तुम्हें एक पूर्ण per-sample वेरिएंट कॉलिंग वर्कफ़्लो मिल जाएगा।
    तुम यह टेस्ट कर सकते हो कि यह सफलतापूर्वक रन हो रहा है, निम्नलिखित कमांड रन करके:

    ```bash
    nextflow run genomics.nf -profile test
    ```

## पाठ योजना

हमने इसे दो चरणों में विभाजित किया है:

1. **Per-sample वेरिएंट कॉलिंग स्टेप को संशोधित करके एक GVCF बनाओ।**
   इसमें प्रोसेस कमांड्स और आउटपुट्स को अपडेट करना शामिल है।
2. **एक जॉइंट जेनोटाइपिंग स्टेप जोड़ो जो per-sample GVCFs को कंबाइन और जेनोटाइप करे।**
   इसमें `collect()` ऑपरेटर, कमांड-लाइन कंस्ट्रक्शन के लिए Groovy closures, और मल्टी-कमांड प्रोसेसेस का परिचय होगा।

!!! note

     सुनिश्चित करो कि तुम सही वर्किंग डायरेक्टरी में हो:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. Per-sample वेरिएंट कॉलिंग स्टेप को GVCF बनाने के लिए संशोधित करो

भाग 2 की पाइपलाइन VCF फ़ाइलें बनाती है, लेकिन जॉइंट कॉलिंग के लिए GVCF फ़ाइलें चाहिए।
हमें GVCF वेरिएंट कॉलिंग मोड को स्विच ऑन करना होगा और आउटपुट फ़ाइल एक्सटेंशन को अपडेट करना होगा।

[भाग 1](01_method.md) से GVCF वेरिएंट कॉलिंग कमांड याद करो:

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

भाग 2 में हमने जिस बेस HaplotypeCaller कमांड को रैप किया था, उसकी तुलना में अंतर `-ERC GVCF` पैरामीटर और `.g.vcf` आउटपुट एक्सटेंशन हैं।

### 1.1. HaplotypeCaller को GVCF emit करने के लिए बताओ और आउटपुट एक्सटेंशन को अपडेट करो

`modules/gatk_haplotypecaller.nf` मॉड्यूल फ़ाइल खोलो और दो बदलाव करो:

- GATK HaplotypeCaller कमांड में `-ERC GVCF` पैरामीटर जोड़ो;
- आउटपुट फ़ाइल पाथ को अपडेट करके संबंधित `.g.vcf` एक्सटेंशन का उपयोग करो, GATK कन्वेंशन के अनुसार।

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

हमें आउटपुट ब्लॉक को भी नए फ़ाइल एक्सटेंशन से मैच करने के लिए अपडेट करना होगा।
चूंकि हमने कमांड आउटपुट को `.vcf` से `.g.vcf` में बदल दिया है, प्रोसेस `output:` ब्लॉक को भी वही बदलाव दर्शाना चाहिए।

### 1.2. प्रोसेस आउटपुट्स ब्लॉक में आउटपुट फ़ाइल एक्सटेंशन को अपडेट करो

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

हमें नए GVCF आउटपुट्स को दर्शाने के लिए वर्कफ़्लो के publish और output कॉन्फ़िगरेशन को भी अपडेट करना होगा।

### 1.3. नए GVCF आउटपुट्स के लिए publish टार्गेट्स को अपडेट करो

चूंकि अब हम VCFs के बजाय GVCFs बना रहे हैं, हमें वर्कफ़्लो के `publish:` सेक्शन को अधिक वर्णनात्मक नामों का उपयोग करने के लिए अपडेट करना चाहिए।
हम स्पष्टता के लिए GVCF फ़ाइलों को अपनी सबडायरेक्टरी में व्यवस्थित भी करेंगे।

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

अब आउटपुट ब्लॉक को मैच करने के लिए अपडेट करो।

### 1.4. नई डायरेक्टरी स्ट्रक्चर के लिए आउटपुट ब्लॉक को अपडेट करो

हमें GVCF फ़ाइलों को एक `gvcf` सबडायरेक्टरी में रखने के लिए `output` ब्लॉक को भी अपडेट करना होगा।

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

मॉड्यूल, publish टार्गेट्स, और आउटपुट ब्लॉक सभी अपडेट होने के साथ, हम बदलावों को टेस्ट कर सकते हैं।

### 1.5. पाइपलाइन को रन करो

बदलावों को वेरिफाई करने के लिए वर्कफ़्लो को रन करो।

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

Nextflow आउटपुट पहले जैसा दिखता है, लेकिन `.g.vcf` फ़ाइलें और उनकी इंडेक्स फ़ाइलें अब सबडायरेक्टरीज़ में व्यवस्थित हैं।

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

अगर तुम किसी एक GVCF फ़ाइल को खोलो और उसे स्क्रॉल करो, तो तुम वेरिफाई कर सकते हो कि GATK HaplotypeCaller ने अनुरोध के अनुसार GVCF फ़ाइलें बनाई हैं।

### सारांश

जब तुम किसी टूल कमांड की आउटपुट फ़ाइलनेम बदलते हो, तो प्रोसेस `output:` ब्लॉक और publish/output कॉन्फ़िगरेशन को मैच करने के लिए अपडेट करना होगा।

### आगे क्या है?

चैनल की सामग्री को इकट्ठा करना और उन्हें अगली प्रोसेस में एक सिंगल इनपुट के रूप में पास करना सीखो।

---

## 2. एक जॉइंट जेनोटाइपिंग स्टेप जोड़ो

अब हमें per-sample GVCFs को इकट्ठा करना होगा, उन्हें एक GenomicsDB डेटा स्टोर में कंबाइन करना होगा, और कोहोर्ट-लेवल VCF बनाने के लिए जॉइंट जेनोटाइपिंग रन करना होगा।
जैसा कि [भाग 1](01_method.md) में बताया गया है, यह दो-टूल ऑपरेशन है: GenomicsDBImport GVCFs को कंबाइन करता है, फिर GenotypeGVCFs अंतिम वेरिएंट कॉल्स बनाता है।
हम दोनों टूल्स को एक सिंगल प्रोसेस में रैप करेंगे जिसे `GATK_JOINTGENOTYPING` कहा जाएगा।

[भाग 1](01_method.md) से दोनों कमांड्स याद करो:

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

पहली कमांड per-sample GVCFs और एक intervals फ़ाइल लेती है, और एक GenomicsDB डेटा स्टोर बनाती है।
दूसरी उस डेटा स्टोर, एक रेफरेंस जीनोम लेती है, और अंतिम कोहोर्ट-लेवल VCF बनाती है।
कंटेनर URI HaplotypeCaller के समान है: `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`।

### 2.1. इनपुट्स सेट अप करो

जॉइंट जेनोटाइपिंग प्रोसेस को दो तरह के इनपुट्स चाहिए जो हमारे पास अभी नहीं हैं: एक arbitrary कोहोर्ट नाम, और सभी सैंपल्स से इकट्ठे किए गए GVCF आउटपुट्स एक साथ बंडल किए हुए।

#### 2.1.1. एक `cohort_name` पैरामीटर जोड़ो

हमें कोहोर्ट के लिए एक arbitrary नाम प्रदान करना होगा।
ट्रेनिंग सीरीज़ में बाद में तुम सीखोगे कि इस तरह की चीज़ों के लिए सैंपल मेटाडेटा का उपयोग कैसे करें, लेकिन अभी हम सिर्फ `params` का उपयोग करके एक CLI पैरामीटर डिक्लेयर करते हैं और सुविधा के लिए इसे एक डिफ़ॉल्ट वैल्यू देते हैं।

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

#### 2.1.2. सैंपल्स में HaplotypeCaller आउटपुट्स को इकट्ठा करो

अगर हम `GATK_HAPLOTYPECALLER` से आउटपुट चैनल को सीधे नई प्रोसेस में प्लग करें, तो Nextflow प्रोसेस को हर सैंपल GVCF पर अलग से कॉल करेगा।
हम सभी तीन GVCFs (और उनकी इंडेक्स फ़ाइलों) को बंडल करना चाहते हैं ताकि Nextflow उन सभी को एक साथ एक सिंगल प्रोसेस कॉल में दे।

हम `collect()` चैनल ऑपरेटर का उपयोग करके ऐसा कर सकते हैं।
GATK_HAPLOTYPECALLER की कॉल के ठीक बाद `workflow` बॉडी में निम्नलिखित लाइनें जोड़ो:

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
2. क्योंकि हमने सेक्शन 1 में `emit:` का उपयोग करके आउटपुट्स को नाम दिया था, हम `.vcf` से GVCFs और `.idx` से इंडेक्स फ़ाइलों को सिलेक्ट कर सकते हैं। नामित आउटपुट्स के बिना, हमें `.out[0]` और `.out[1]` का उपयोग करना पड़ता।
3. `collect()` ऑपरेटर सभी फ़ाइलों को एक सिंगल एलिमेंट में बंडल करता है, इसलिए `all_gvcfs_ch` में सभी तीन GVCFs एक साथ होते हैं, और `all_idxs_ch` में सभी तीन इंडेक्स फ़ाइलें एक साथ होती हैं।

हम GVCFs और उनकी इंडेक्स फ़ाइलों को अलग से इकट्ठा कर सकते हैं (tuples में उन्हें एक साथ रखने के विपरीत) क्योंकि Nextflow execution के लिए सभी इनपुट फ़ाइलों को एक साथ stage करेगा, इसलिए इंडेक्स फ़ाइलें GVCFs के साथ मौजूद होंगी।

!!! tip

    चैनल ऑपरेटर्स लागू करने से पहले और बाद में चैनलों की सामग्री को इंस्पेक्ट करने के लिए तुम `view()` ऑपरेटर का उपयोग कर सकते हो।

### 2.2. जॉइंट जेनोटाइपिंग प्रोसेस लिखो और इसे वर्कफ़्लो में कॉल करो

भाग 2 में उपयोग किए गए समान पैटर्न का पालन करते हुए, हम मॉड्यूल फ़ाइल में प्रोसेस डेफिनिशन लिखेंगे, इसे वर्कफ़्लो में इंपोर्ट करेंगे, और इसे जिन इनपुट्स को हमने अभी तैयार किया है, उन पर कॉल करेंगे।

#### 2.2.1. हर GVCF को `-V` आर्ग्युमेंट देने के लिए एक स्ट्रिंग कंस्ट्रक्ट करो

प्रोसेस डेफिनिशन को भरना शुरू करने से पहले, एक चीज़ काम करनी है।
GenomicsDBImport कमांड हर GVCF फ़ाइल के लिए एक अलग `-V` आर्ग्युमेंट की उम्मीद करता है, इस तरह:

```bash
gatk GenomicsDBImport \
    -V reads_mother.bam.g.vcf \
    -V reads_father.bam.g.vcf \
    -V reads_son.bam.g.vcf \
    ...
```

अगर हम `-V ${all_gvcfs_ch}` लिखें, तो Nextflow बस फ़ाइलनेम्स को concatenate कर देगा और कमांड का वह हिस्सा इस तरह दिखेगा:

```groovy
-V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf
```

लेकिन हमें स्ट्रिंग को इस तरह दिखना चाहिए:

```groovy
-V reads_mother.bam.g.vcf -V reads_father.bam.g.vcf -V reads_son.bam.g.vcf
```

महत्वपूर्ण रूप से, हमें इस स्ट्रिंग को collected चैनल में जो भी फ़ाइलें हैं, उनसे dynamically कंस्ट्रक्ट करनी होगी।
Nextflow (Groovy के माध्यम से) ऐसा करने का एक संक्षिप्त तरीका प्रदान करता है:

```groovy
def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
```

इसे तोड़ते हुए:

1. `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }` हर फ़ाइल पाथ पर iterate करता है और उसमें `-V ` prepend करता है, `["-V A.g.vcf", "-V B.g.vcf", "-V C.g.vcf"]` बनाता है।
2. `.join(' ')` उन्हें स्पेस के साथ concatenate करता है: `"-V A.g.vcf -V B.g.vcf -V C.g.vcf"`।
3. परिणाम एक लोकल वेरिएबल `gvcfs_line` को असाइन किया जाता है (`def` के साथ डिफाइन किया गया), जिसे हम कमांड टेम्पलेट में interpolate कर सकते हैं।

यह लाइन प्रोसेस के `script:` ब्लॉक के अंदर जाती है, कमांड टेम्पलेट से पहले।
तुम कमांड टेम्पलेट के opening `"""` से पहले `script:` के बीच arbitrary Groovy कोड रख सकते हो।

फिर तुम प्रोसेस के `script:` ब्लॉक में उस पूरी स्ट्रिंग को `gvcfs_line` के रूप में refer कर पाओगे।

#### 2.2.2. जॉइंट जेनोटाइपिंग प्रोसेस के लिए मॉड्यूल को भरो

अब हम पूर्ण प्रोसेस लिखने की ओर बढ़ सकते हैं।

`modules/gatk_jointgenotyping.nf` खोलो और प्रोसेस डेफिनिशन की outline को एग्ज़ामिन करो।

ऊपर दी गई जानकारी का उपयोग करके प्रोसेस डेफिनिशन को भरो, फिर नीचे "बाद में" टैब में सॉल्यूशन के खिलाफ अपना काम चेक करो।

=== "पहले"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * GVCFs को GenomicsDB datastore में कंबाइन करो और कोहोर्ट-लेवल कॉल्स बनाने के लिए जॉइंट जेनोटाइपिंग रन करो
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
     * GVCFs को GenomicsDB datastore में कंबाइन करो और कोहोर्ट-लेवल कॉल्स बनाने के लिए जॉइंट जेनोटाइपिंग रन करो
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

यहाँ कई चीज़ें उल्लेखनीय हैं।

पहले की तरह, कई इनपुट्स लिस्ट किए गए हैं भले ही कमांड्स उन्हें सीधे रेफरेंस नहीं करते: `all_idxs`, `ref_index`, और `ref_dict`।
उन्हें लिस्ट करना सुनिश्चित करता है कि Nextflow इन फ़ाइलों को वर्किंग डायरेक्टरी में उन फ़ाइलों के साथ stage करे जो कमांड्स में दिखाई देती हैं, जिन्हें GATK naming conventions के आधार पर खोजने की उम्मीद करता है।

`gvcfs_line` वेरिएबल GenomicsDBImport के लिए `-V` आर्ग्युमेंट्स कंस्ट्रक्ट करने के लिए ऊपर वर्णित Groovy closure का उपयोग करता है।

यह प्रोसेस serial में दो कमांड्स रन करता है, ठीक जैसे तुम टर्मिनल में करते। GenomicsDBImport per-sample GVCFs को एक डेटा स्टोर में कंबाइन करता है, फिर GenotypeGVCFs उस डेटा स्टोर को पढ़ता है और अंतिम कोहोर्ट-लेवल VCF बनाता है।
GenomicsDB डेटा स्टोर (`${cohort_name}_gdb`) एक intermediate artifact है जिसका उपयोग केवल प्रोसेस के अंदर किया जाता है; यह आउटपुट ब्लॉक में दिखाई नहीं देता।

एक बार जब तुम यह पूरा कर लो, तो प्रोसेस उपयोग के लिए तैयार है।
वर्कफ़्लो में इसका उपयोग करने के लिए, तुम्हें मॉड्यूल को इंपोर्ट करना होगा और एक प्रोसेस कॉल जोड़नी होगी।

#### 2.2.3. मॉड्यूल को इंपोर्ट करो

मौजूदा इंपोर्ट स्टेटमेंट्स के नीचे `genomics.nf` में इंपोर्ट स्टेटमेंट जोड़ो:

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

प्रोसेस अब वर्कफ़्लो scope में उपलब्ध है।

#### 2.2.4. प्रोसेस कॉल जोड़ो

`collect()` लाइनों के बाद वर्कफ़्लो बॉडी में `GATK_JOINTGENOTYPING` की कॉल जोड़ो:

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
अगला, हम configure करेंगे कि आउटपुट्स कैसे publish होते हैं।

### 2.3. आउटपुट हैंडलिंग को कॉन्फ़िगर करो

हमें जॉइंट VCF आउटपुट्स को publish करना होगा।
जॉइंट जेनोटाइपिंग परिणामों के लिए publish टार्गेट्स और आउटपुट ब्लॉक एंट्रीज़ जोड़ो।

#### 2.3.1. जॉइंट VCF के लिए publish टार्गेट्स जोड़ो

वर्कफ़्लो के `publish:` सेक्शन में जॉइंट VCF और उसका इंडेक्स जोड़ो:

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

अब मैच करने के लिए आउटपुट ब्लॉक को अपडेट करो।

#### 2.3.2. जॉइंट VCF के लिए आउटपुट ब्लॉक एंट्रीज़ जोड़ो

जॉइंट VCF फ़ाइलों के लिए एंट्रीज़ जोड़ो।
हम उन्हें results डायरेक्टरी की रूट पर रखेंगे क्योंकि यह अंतिम आउटपुट है।

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

प्रोसेस, publish टार्गेट्स, और आउटपुट ब्लॉक सभी के साथ, हम पूर्ण वर्कफ़्लो को टेस्ट कर सकते हैं।

### 2.4. वर्कफ़्लो को रन करो

सब कुछ काम करता है वेरिफाई करने के लिए वर्कफ़्लो को रन करो।

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

पहले दो स्टेप्स पिछले रन से cached हैं, और नया `GATK_JOINTGENOTYPING` स्टेप सभी तीन सैंपल्स से collected इनपुट्स पर एक बार रन होता है।
अंतिम आउटपुट फ़ाइल, `family_trio.joint.vcf` (और उसका इंडेक्स), results डायरेक्टरी में है।

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

अगर तुम जॉइंट VCF फ़ाइल खोलो, तो तुम वेरिफाई कर सकते हो कि वर्कफ़्लो ने अपेक्षित वेरिएंट कॉल्स बनाए।

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

अब तुम्हारे पास एक automated, पूर्ण रूप से reproducible जॉइंट वेरिएंट कॉलिंग वर्कफ़्लो है!

!!! note

    ध्यान रखो कि हमने तुम्हें जो डेटा फ़ाइलें दी हैं वे chromosome 20 के केवल एक छोटे से हिस्से को कवर करती हैं।
    वेरिएंट कॉलसेट का वास्तविक साइज़ लाखों वेरिएंट्स में गिना जाएगा।
    इसीलिए हम ट्रेनिंग उद्देश्यों के लिए केवल छोटे डेटा सबसेट्स का उपयोग करते हैं!

### सारांश

तुम जानते हो कि चैनल से आउटपुट्स कैसे इकट्ठा करें और उन्हें किसी अन्य प्रोसेस के लिए एक सिंगल इनपुट के रूप में बंडल करें।
तुम Groovy closures का उपयोग करके कमांड लाइन कैसे कंस्ट्रक्ट करें, और एक सिंगल प्रोसेस में मल्टिपल कमांड्स कैसे रन करें, यह भी जानते हो।

### आगे क्या है?

अपनी सफलता का जश्न मनाओ और एक अच्छा ब्रेक लो।

इस कोर्स के अगले भाग में, तुम सीखोगे कि nf-core से एक production-ready वेरिएंट कॉलिंग पाइपलाइन कैसे रन करें और इसकी तुलना उस पाइपलाइन से कैसे करें जो तुमने मैन्युअली बनाई है।
