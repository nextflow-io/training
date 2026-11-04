# भाग 2: कमांड लाइन से पाइपलाइन लॉन्च करना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[भाग 1](./01_run_with_seqera.md) में, तुमने Seqera के वेब इंटरफ़ेस से nf-core/rnaseq लॉन्च किया था।
अब हम यही काम `tw` CLI का उपयोग करके कमांड लाइन से करेंगे, और अपने workspace में एक नई पाइपलाइन जोड़ेंगे।

---

## 1. कमांड लाइन से पाइपलाइन लॉन्च करना

run view में, **Command line** टैब पर क्लिक करो।
तुम्हें वह सटीक `nextflow run` कमांड दिखेगी जो Platform ने तुम्हारी ओर से बनाई और सबमिट की — वही तरह की कमांड जो तुम Use nf-core कोर्स में मैन्युअली चला रहे थे।

Platform, Nextflow की जगह नहीं लेता; यह उसे ऑर्केस्ट्रेट करता है।
जो कुछ भी तुम वेब इंटरफ़ेस के ज़रिए कर सकते हो, वह सब `tw` CLI से टर्मिनल में भी किया जा सकता है — यह Platform API के साथ इंटरैक्ट करने का कमांड-लाइन टूल है।
यह स्क्रिप्ट या CI/CD पाइपलाइन से लॉन्च को ऑटोमेट करने के लिए उपयोगी है।

हम यह अभी उसी codespace से करेंगे जो तुमने पहले के कोर्स के लिए उपयोग किया था।

### 1.1. tw CLI इंस्टॉल करना

`tw` बाइनरी डाउनलोड और इंस्टॉल करने के लिए अपने Codespace टर्मिनल में निम्नलिखित कमांड चलाओ:

```bash
curl -fsSL https://github.com/seqeralabs/tower-cli/releases/latest/download/tw-linux-x86_64 -o tw
chmod +x tw
sudo mv tw /usr/local/bin/
```

इंस्टॉलेशन की जाँच करो:

```bash
tw --version
```

??? success "कमांड आउटपुट"

    ```console
    Tower CLI version 0.40.0 (build 26db579)
    ```

`tw` CLI इंस्टॉल हो गई है और कॉन्फ़िगर करने के लिए तैयार है।

### 1.2. एक्सेस टोकन प्राप्त करना

`tw` CLI, Seqera के साथ एक personal access token का उपयोग करके प्रमाणित होती है।

1. Seqera वेब इंटरफ़ेस में, ऊपर-दाईं ओर अपने अवतार पर क्लिक करो और **Your tokens** चुनो।
2. **Add token** पर क्लिक करो, उसे एक नाम दो (जैसे `training`), और **Add** पर क्लिक करो।
3. टोकन की वैल्यू कॉपी करो — यह केवल एक बार दिखाई जाएगी।
   अगर तुम इसे तुरंत कहीं सेव नहीं करते, तो तुम्हें एक नया टोकन जनरेट करना होगा।

### 1.3. CLI कॉन्फ़िगर करना

सुविधा के लिए, हम एक कॉन्फ़िगरेशन फ़ाइल सेट अप करेंगे जिसमें
तुम्हारा अभी जनरेट किया गया एक्सेस टोकन और workspace आइडेंटिफ़ायर होगा।

इस डायरेक्टरी में `.seqera_config` फ़ाइल को एडिटर में खोलो और दो वेरिएबल सेट करो:

- **`TOWER_ACCESS_TOKEN`**: वह टोकन जो तुमने सेक्शन 1.2 में जनरेट किया
- **`TOWER_WORKSPACE_ID`**: तुम्हारे workspace का न्यूमेरिक ID (`tw workspaces list` में `ID` कॉलम, जो तुम सेक्शन 1.4 में चलाओगे)

वैल्यू भरने के बाद, config लोड करो:

```bash
source .seqera_config
```

कनेक्शन की जाँच करो:

```bash
tw info
```

??? success "कमांड आउटपुट"

    ```console
        Details
    -------------------------+-----------------------------
     Tower API endpoint      | https://api.cloud.seqera.io
     Tower API version       | 1.150.0
     Tower version           | 26.1.0-cycle54
     CLI version             | 0.30.0 (fde9dec)
     CLI minimum API version | 1.148.0
     Authenticated user      | <your-name>

    System health status
    ---------------------------------------+----
     Remote API server connection check    | OK
     Tower API version check               | OK
     Authentication API credential's token | OK
    ```

`tw` CLI अब प्रमाणित है और तुम्हारे Seqera अकाउंट से कनेक्ट है।
config को रीलोड करने के लिए हर Codespace सेशन की शुरुआत में `source .seqera_config` चलाओ।

!!! tip "सुझाव"

    अगर तुम्हारे workspace में कोई primary compute environment सेट नहीं है, तो तुम अपनी config फ़ाइल में `export TOWER_COMPUTE_ENV=<compute-env-name>` जोड़ सकते हो।
    किसी भी config वैल्यू को कमांड लाइन पर फ्लैग पास करके ओवरराइड किया जा सकता है (जैसे `--compute-env other-env`)।
    सभी ऑप्शन और एनवायरनमेंट वेरिएबल की पूरी सूची के लिए [tw CLI reference](https://docs.seqera.io/platform/latest/cli/reference) देखो।

### 1.4. CLI से अपना workspace एक्सप्लोर करना

उन workspaces की सूची देखो जिन तक तुम्हारी पहुँच है:

```bash
tw workspaces list
```

??? success "कमांड आउटपुट"

    ```console
    Available workspaces:
    ID              | Name            | Full Name                | Visibility
    --------------- | --------------- | ------------------------ | ----------
    <workspace-id>  | my-workspace    | my-org/my-workspace      | PRIVATE
    ```

अपने workspace में runs देखो, जिसमें वह nf-core/rnaseq run भी शामिल है जो तुमने अभी लॉन्च की:

```bash
tw runs list
```

??? success "कमांड आउटपुट"

    ```console
    Pipeline runs for my-org/my-workspace workspace:
    ID        | Status   | Name          | Pipeline               | Run name
    --------- | -------- | ------------- | ---------------------- | --------
    <run-id>  | RUNNING  | ...           | nf-core/rnaseq         | happy_curie
    ```

वही run जो तुम वेब इंटरफ़ेस में मॉनिटर कर रहे हो, यहाँ भी दिखाई दे रही है।

!!! note "नोट"

    क्योंकि `.seqera_config` में `TOWER_WORKSPACE_ID` सेट है, तुम सभी `tw` कमांड से `--workspace` हटा सकते हो।
    config के बिना, तुम्हें इसे स्पष्ट रूप से पास करना होगा:

    ```bash
    tw runs list --workspace <org>/<workspace>
    ```

वेब इंटरफ़ेस में जो कुछ भी दिखता है, वह CLI से भी एक्सेस किया जा सकता है।

### 1.5. CLI से nf-core/rnaseq लॉन्च करना

[भाग 1](./01_run_with_seqera.md) में तुमने अपने workspace में जो पाइपलाइन जोड़ी थी, वह CLI में नाम से उपलब्ध है।
इसे `test` प्रोफ़ाइल के साथ लॉन्च करो:

```bash
tw launch nf-core-rnaseq --profile test
```

??? success "कमांड आउटपुट"

    ```console
    Launching pipeline nf-core-rnaseq
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

लिंक को अपने ब्राउज़र में खोलो और पुष्टि करो कि run **Runs** पैनल में दिखाई दे रही है।

एक बार जब तुम इसे चलते हुए देख लो, तो तुमने पुष्टि कर ली है कि CLI और वेब इंटरफ़ेस एक ही workspace के दो अलग-अलग दृश्य हैं।

!!! note "नोट"

    तुम पाइपलाइन को workspace में जोड़े बिना भी `tw launch` में सीधे GitHub URL पास कर सकते हो।
    हालाँकि, लॉन्च करने से पहले पाइपलाइन को स्पष्ट रूप से जोड़ना आमतौर पर बेहतर होता है: यह भविष्य की runs के लिए पाइपलाइन कॉन्फ़िगरेशन सेव करता है, इसे नाम से उपलब्ध कराता है, और Launchpad में सभी workspace सदस्यों को दिखाई देता है।

    `tw` का उपयोग करके कमांड लाइन से सीधे workspace में पाइपलाइन जोड़ना संभव है।
    अगला सेक्शन दिखाता है कि nf-core/demo पाइपलाइन के साथ यह कैसे करें।

### सारांश

तुम जानते हो कि `tw` CLI को कैसे प्रमाणित करें, अपना workspace कैसे देखें, और टर्मिनल से एक सेव की गई पाइपलाइन कैसे लॉन्च करें।

### आगे क्या है?

कमांड लाइन से अपने workspace में एक नई पाइपलाइन जोड़ो और उसे लॉन्च करो।

---

## 2. एक नई पाइपलाइन जोड़ना और चलाना

GitHub पर कोई भी Nextflow पाइपलाइन `tw pipelines add` से तुम्हारे workspace में जोड़ी जा सकती है, जब तक उसके root में `main.nf` एंट्री पॉइंट और `nextflow.config` हो।
nf-core/demo इसके साथ अभ्यास करने का एक अच्छा उदाहरण है: तुम इसे Use nf-core कोर्स में पहले ही चला चुके हो, इसलिए तुम जानते हो यह क्या करता है और क्या उम्मीद करनी है।

### 2.1. अपने workspace में nf-core/demo जोड़ना

workspace में पाइपलाइन रजिस्टर करने के लिए निम्नलिखित कमांड चलाओ:

```bash
tw pipelines add \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "कमांड आउटपुट"

    ```console
    New pipeline 'nf-core-demo' added at my-org/my-workspace workspace
    ```

पाइपलाइन अब रजिस्टर हो गई है और Launchpad में दिखाई देगी।

### 2.2. Launchpad में इसकी उपस्थिति की जाँच करना

यह पुष्टि करने के लिए कि पाइपलाइन जोड़ी गई, अपने workspace में पाइपलाइन की सूची देखो:

```bash
tw pipelines list
```

??? success "कमांड आउटपुट"

    ```console
    Pipelines at my-org/my-workspace:
    ID   | Name            | Repository
    ---- | --------------- | -----------------------------------------
    ...  | nf-core-demo    | https://github.com/nf-core/demo
    ...  | nf-core-rnaseq  | ...
    ```

अपने workspace को ब्राउज़र में खोलो और **Launchpad** पर क्लिक करो यह पुष्टि करने के लिए कि nf-core/demo अब nf-core/rnaseq के साथ दिखाई दे रही है।

!!! tip "सुझाव"

    तुम वेब इंटरफ़ेस के ज़रिए भी पाइपलाइन जोड़ सकते हो: बाईं साइडबार में **Launchpad** पर क्लिक करो, फिर **Add pipeline** पर क्लिक करो, और फ़ॉर्म भरो।

nf-core/demo एंट्री पर **Launch** बटन पर क्लिक करो उसका launch फ़ॉर्म खोलने के लिए।
तुम देखोगे कि `input` और `outdir` पैरामीटर लाल रंग में हाइलाइट हैं — ये आवश्यक फ़ील्ड हैं जिनकी कोई डिफ़ॉल्ट वैल्यू नहीं है, क्योंकि `tw pipelines add` केवल पाइपलाइन सोर्स रजिस्टर करता है बिना कोई पैरामीटर पूर्व-कॉन्फ़िगर किए।
अगले दो सेक्शन बताते हैं कि वे वैल्यू कैसे प्रदान करें: पहले वेब फ़ॉर्म के ज़रिए, फिर कमांड लाइन से।

### 2.3. वेब इंटरफ़ेस से nf-core/demo लॉन्च करना

launch फ़ॉर्म खुला होने पर, दो आवश्यक पैरामीटर भरो।

`input` के लिए, nf-core/demo test प्रोफ़ाइल से test samplesheet URL दर्ज करो।
तुम इसे पाइपलाइन रिपॉज़िटरी के अंदर `conf/test.config` में पा सकते हो, जिसे तुमने Use nf-core कोर्स में देखा था:

```
https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
```

`outdir` के लिए, एक cloud storage पाथ दर्ज करो जहाँ पाइपलाइन अपने परिणाम लिख सके।
अपने workspace के लिए कॉन्फ़िगर किए गए bucket का उपयोग करो, runs को व्यवस्थित रखने के लिए एक सबडायरेक्टरी के साथ:

```
s3://my-bucket/demo-results
```

दोनों फ़ील्ड भरने के बाद, नीले **Launch** बटन पर क्लिक करो।

run **Runs** पैनल में दिखाई देती है और test डेटासेट पर कुछ मिनटों में पूरी हो जानी चाहिए।
task टेबल और किसी भी execution रिपोर्ट को एक्सप्लोर करने के लिए run पर क्लिक करो।

### 2.4. CLI से nf-core/demo लॉन्च करना

`nextflow run` के विपरीत, `tw launch` कमांड `--input` या `--outdir` जैसे अलग-अलग पैरामीटर फ्लैग स्वीकार नहीं करती।
पैरामीटर YAML या JSON फ़ॉर्मेट की फ़ाइल के ज़रिए प्रदान किए जाने चाहिए, जिसे `--params-file` के साथ पास किया जाता है।
यह reproducibility को प्रोत्साहित करता है: एक सेव की गई पैरामीटर फ़ाइल दस्तावेज़ करती है कि किसी run के लिए कौन सी वैल्यू उपयोग की गईं, जिससे run कॉन्फ़िगरेशन को दोहराना या शेयर करना आसान हो जाता है।

अपनी working डायरेक्टरी में एक पैरामीटर फ़ाइल बनाओ:

```bash
touch params.yaml
```

इसे एडिटर में खोलो और आउटपुट पाथ जोड़ो:

```yaml title="params.yaml"
outdir: "s3://my-bucket/demo-results-cli"
```

अब तुम `test` प्रोफ़ाइल (जो `input` samplesheet प्रदान करती है) और params फ़ाइल (जो `outdir` प्रदान करती है) का उपयोग करके पाइपलाइन लॉन्च कर सकते हो:

```bash
tw launch nf-core-demo --profile test --params-file params.yaml
```

??? success "कमांड आउटपुट"

    ```console
    Launching pipeline nf-core-demo
    Run name: focused_feynman
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

यह पुष्टि करने के लिए लिंक खोलो कि run **Runs** पैनल में दिखाई दे रही है।

!!! tip "सुझाव"

    अगर तुम कुछ डिफ़ॉल्ट सेट करना चाहते हो, तो पैरामीटर फ़ाइल को initial setup स्टेप के दौरान शामिल कर सकते हो, साथ ही कुछ अतिरिक्त प्रॉपर्टी भी जो हमने पहले वेब फ़ॉर्म के ज़रिए की थीं:

    ```bash
    tw pipelines add \
      --name nf-core-demo-gg \
      --description "Demo pipeline with defaults for testing" \
      --params-file params.yaml \
      --profile test \
      https://github.com/nf-core/demo
    ```

### सारांश

तुम जानते हो कि GitHub पर होस्ट किसी भी Nextflow पाइपलाइन को अपने workspace में कैसे जोड़ें और लॉन्च करें — वेब इंटरफ़ेस से पैरामीटर मैन्युअली भरकर, और `tw` CLI से प्रोफ़ाइल और पैरामीटर फ़ाइल को मिलाकर।

---

## सारांश

इस भाग में तुमने सीखा:

- `tw` CLI को प्रमाणित करना और टर्मिनल से एक सेव की गई पाइपलाइन लॉन्च करना
- CLI का उपयोग करके GitHub से एक नई पाइपलाइन जोड़ना और यह पुष्टि करना कि वह Launchpad में दिखाई देती है
- Seqera वेब इंटरफ़ेस से आवश्यक पैरामीटर मैन्युअली भरकर पाइपलाइन लॉन्च करना
- Nextflow प्रोफ़ाइल और पैरामीटर फ़ाइल का उपयोग करके CLI से पाइपलाइन लॉन्च करना
