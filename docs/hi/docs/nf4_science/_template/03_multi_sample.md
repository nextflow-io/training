# भाग 3: मल्टी-सैंपल एग्रीगेशन

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

भाग 2 में, तुमने एक per-sample प्रोसेसिंग पाइपलाइन बनाई जो प्रत्येक सैंपल को स्वतंत्र रूप से हैंडल करती थी।
अब हम इसे विस्तारित करके मल्टी-सैंपल {AGGREGATION_METHOD} को लागू करेंगे, जैसा कि [भाग 1](01_method.md) में बताया गया है।

## असाइनमेंट

इस कोर्स के इस भाग में, हम वर्कफ़्लो को निम्नलिखित कार्य करने के लिए विस्तारित करेंगे:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

1. {PER_SAMPLE_STEP_1}
2. {PER_SAMPLE_STEP_2}
3. {AGGREGATION_STEP_1}
4. {AGGREGATION_STEP_2}

यह भाग सीधे भाग 2 द्वारा तैयार की गई वर्कफ़्लो पर आधारित है।

??? info "इस खंड से कैसे शुरू करें"

    कोर्स का यह सेक्शन मानता है कि तुमने [भाग 2: सिंगल-सैंपल प्रोसेसिंग](./02_single_sample.md) पूरा कर लिया है और तुम्हारे पास एक कार्यशील `{DOMAIN_DIR}.nf` पाइपलाइन है।

    यदि तुमने भाग 2 पूरा नहीं किया है या इस भाग के लिए नए सिरे से शुरू करना चाहते हो, तो तुम भाग 2 के समाधान को अपने शुरुआती बिंदु के रूप में उपयोग कर सकते हो।
    `nf4-science/{DOMAIN_DIR}/` डायरेक्टरी के अंदर से ये कमांड चलाओ:

    ```bash
    cp solutions/part2/{DOMAIN_DIR}-2.nf {DOMAIN_DIR}.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    यह तुम्हें एक पूर्ण सिंगल-सैंपल प्रोसेसिंग वर्कफ़्लो देता है।
    तुम निम्नलिखित कमांड चलाकर परीक्षण कर सकते हो कि यह सफलतापूर्वक चलती है:

    ```bash
    nextflow run {DOMAIN_DIR}.nf -profile test
    ```

## पाठ योजना

हमने इसे दो चरणों में विभाजित किया है:

1. **{MODIFICATION_STEP_SUMMARY}.**
   यह प्रोसेस कमांड और आउटपुट को अपडेट करना कवर करता है।
2. **{AGGREGATION_STEP_SUMMARY}.**
   यह `collect()` ऑपरेटर {AND_OTHER_CONCEPTS} का परिचय देता है।

!!! note "नोट"

     सुनिश्चित करो कि तुम सही वर्किंग डायरेक्टरी में हो:
     `cd /workspaces/training/nf4-science/{DOMAIN_DIR}`

---

## 1. {MODIFICATION_STEP_TITLE}

{DESCRIPTION_OF_MODIFICATION_TO_EXISTING_PROCESSES}

[भाग 1](01_method.md) से संशोधित कमांड को याद करो:

```bash
{MODIFIED_COMMAND}
```

{EXPLAIN_DIFFERENCES_FROM_PART_2}

### 1.1. {MODIFICATION_SUBSTEP}

{INSTRUCTIONS_WITH_BEFORE_AFTER_TABS}

=== "बाद में"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {UPDATED_CODE}
    ```

=== "पहले"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {ORIGINAL_CODE}
    ```

### 1.2. संशोधन को सत्यापित करने के लिए वर्कफ़्लो चलाओ

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "कमांड आउटपुट"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### सारांश

तुम जानते हो कि वर्कफ़्लो के व्यवहार को अनुकूलित करने के लिए प्रोसेस कमांड और आउटपुट को कैसे संशोधित करें।

### आगे क्या है?

मल्टी-सैंपल एग्रीगेशन स्टेप जोड़ो।

---

## 2. {AGGREGATION_STEP_TITLE}

{DESCRIPTION_OF_AGGREGATION}

### 2.1. एग्रीगेशन मॉड्यूल लिखो

{MODULE_INSTRUCTIONS}

### 2.2. Per-sample आउटपुट को collect करो और उन्हें एग्रीगेशन प्रोसेस में फीड करो

{INSTRUCTIONS_USING_COLLECT_OPERATOR}

### 2.3. पूर्ण वर्कफ़्लो चलाओ

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "कमांड आउटपुट"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

---

### सारांश

तुम्हारे पास एक पूर्ण पाइपलाइन है जो सैंपल को व्यक्तिगत रूप से प्रोसेस करती है और सभी सैंपल के परिणामों को एग्रीगेट करती है।
तुम जानते हो कि मल्टी-सैंपल विश्लेषण के लिए per-sample आउटपुट को एग्रीगेट करने के लिए `collect()` जैसे चैनल ऑपरेटर का उपयोग कैसे करें।

### आगे क्या है?

इस कोर्स को पूरा करने पर बधाई! तुमने जो सीखा है उसकी समीक्षा करने और अगले चरणों का पता लगाने के लिए [कोर्स सारांश](next_steps.md) पर जाओ।
