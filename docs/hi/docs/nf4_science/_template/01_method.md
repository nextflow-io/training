# भाग 1: विधि का अवलोकन और मैनुअल परीक्षण

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

{BRIEF_METHOD_DESCRIPTION}

![पाइपलाइन अवलोकन](img/{PIPELINE_DIAGRAM}.png)

{OPTIONAL_EXTENDED_METHOD_DESCRIPTION}

### विधियाँ

{DESCRIBE_THE_TWO_APPROACHES: single-sample और multi-sample/aggregation}

इससे पहले कि हम कोई भी वर्कफ़्लो कोड लिखना शुरू करें, हम कुछ टेस्ट डेटा पर कमांड को मैनुअली आज़माने जा रहे हैं।

### डेटासेट

हम निम्नलिखित डेटा और संबंधित संसाधन प्रदान करते हैं:

- **{PRIMARY_INPUT_DESCRIPTION}**
- **{SAMPLE_DESCRIPTION}**
- **{ADDITIONAL_RESOURCES_DESCRIPTION}**

### सॉफ़्टवेयर

मुख्य टूल हैं [{TOOL_A}]({TOOL_A_URL}) और [{TOOL_B}]({TOOL_B_URL})।

ये टूल GitHub Codespaces वातावरण में इंस्टॉल नहीं हैं, इसलिए हम इन्हें कंटेनर के माध्यम से उपयोग करेंगे ([Hello Containers](../../hello_nextflow/05_hello_containers.md) देखें)।

!!! note "नोट"

    सुनिश्चित करो कि तुम `nf4-science/{DOMAIN_DIR}` डायरेक्टरी में हो ताकि जब तुम `pwd` टाइप करो तो पाथ का अंतिम भाग `{DOMAIN_DIR}` दिखे।

---

## 1. {SINGLE_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_SINGLE_SAMPLE_APPROACH}

इस सेक्शन में हम उन कमांड का परीक्षण करते हैं जो single-sample प्रोसेसिंग दृष्टिकोण को बनाते हैं।
ये वही कमांड हैं जिन्हें हम इस कोर्स के भाग 2 में Nextflow वर्कफ़्लो में wrap करेंगे।

1. {STEP_1_SUMMARY}
2. {STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

हम सिर्फ़ एक नमूने पर कमांड का परीक्षण करके शुरू करते हैं।

### 1.1. {FIRST_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.1.1. कंटेनर को Pull करो

कंटेनर इमेज डाउनलोड करने के लिए `docker pull` कमांड चलाओ:

```bash
docker pull {TOOL_A_CONTAINER_URI}
```

??? success "कमांड आउटपुट"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.1.2. कंटेनर को इंटरैक्टिवली स्पिन अप करो

कंटेनर को स्पिन अप करो और `data` डायरेक्टरी को माउंट करो ताकि टूल इनपुट फ़ाइलों तक पहुँच सकें:

```bash
docker run -it -v ./data:/data {TOOL_A_CONTAINER_URI}
```

तुम्हारा प्रॉम्प्ट बदल जाता है जो दर्शाता है कि तुम कंटेनर के अंदर हो।

#### 1.1.3. कमांड चलाओ

```bash
{TOOL_A_COMMAND}
```

??? success "कमांड आउटपुट"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.1.4. कंटेनर से बाहर निकलो

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करो।

```bash
exit
```

तुम्हारा प्रॉम्प्ट सामान्य हो जाना चाहिए।

### 1.2. {SECOND_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.2.1. कंटेनर को Pull करो

```bash
docker pull {TOOL_B_CONTAINER_URI}
```

??? success "कमांड आउटपुट"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.2.2. कंटेनर को इंटरैक्टिवली स्पिन अप करो

```bash
docker run -it -v ./data:/data {TOOL_B_CONTAINER_URI}
```

#### 1.2.3. कमांड चलाओ

```bash
{TOOL_B_COMMAND}
```

??? success "कमांड आउटपुट"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.2.4. कंटेनर से बाहर निकलो

कंटेनर से बाहर निकलने के लिए, `exit` टाइप करो।

```bash
exit
```

तुम्हारा प्रॉम्प्ट सामान्य हो जाना चाहिए।
यह single-sample प्रोसेसिंग टेस्ट को समाप्त करता है।

---

## 2. {MULTI_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_MULTI_SAMPLE_APPROACH}

{EXPLAIN_WHY_MULTI_SAMPLE_IS_NEEDED}

इस सेक्शन में हम multi-sample प्रोसेसिंग के लिए आवश्यक अतिरिक्त कमांड का परीक्षण करते हैं।
ये वही कमांड हैं जिन्हें हम इस कोर्स के भाग 3 में Nextflow वर्कफ़्लो में wrap करेंगे।

1. {AGGREGATION_STEP_1_SUMMARY}
2. {AGGREGATION_STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

### 2.1. {AGGREGATION_STEP_TITLE}

{STEP_INSTRUCTIONS_FOLLOWING_SAME_PATTERN_AS_ABOVE}

---

### सारांश

तुम जानते हो कि {TOOL_A} और {TOOL_B} कमांड को उनके संबंधित कंटेनर में कैसे टेस्ट करना है, जिसमें {MULTI_SAMPLE_SUMMARY} कैसे करना है यह भी शामिल है।

### आगे क्या है?

सीखो कि उन्हीं कमांड को वर्कफ़्लो में कैसे wrap करना है जो काम को execute करने के लिए कंटेनर का उपयोग करते हैं।
