# भाग 3: इनपुट को व्यवस्थित करना

भाग 2 में, हमने कमांड लाइन पर कई पैरामीटर के साथ molkart चलाया।
अब हम इनपुट को मैनेज करने के दो बेहतर तरीके सीखेंगे: **पैरामीटर फ़ाइलें** और **सैंपलशीट**।

## 1. पैरामीटर फ़ाइलों का उपयोग करना

### 1.1. लंबी कमांड लाइन की समस्या

भाग 2 से हमारी कमांड याद करो:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results
```

यह काम करता है, लेकिन इसे दोबारा बनाना, शेयर करना या संशोधित करना मुश्किल है।
क्या होगा अगर तुम्हें अगले महीने फिर से यही विश्लेषण चलाना हो?
क्या होगा अगर कोई सहयोगी तुम्हारी सटीक सेटिंग्स का उपयोग करना चाहे?

### 1.2. समाधान: पैरामीटर फ़ाइल का उपयोग करो

`params.yaml` नाम की एक फ़ाइल बनाओ:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

अब तुम्हारी कमांड बन जाती है:

```bash
nextflow run ./molkart -params-file params.yaml -resume
```

बस इतना ही! पैरामीटर फ़ाइल तुम्हारी सटीक कॉन्फ़िगरेशन को डॉक्यूमेंट करती है और इसे दोबारा चलाना या शेयर करना आसान बनाती है।

### 1.3. पैरामीटर को ओवरराइड करना

तुम अभी भी कमांड लाइन से विशिष्ट पैरामीटर को ओवरराइड कर सकते हो:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

ऊपर की लाइन `segmentation_method` को `stardist` में बदल देती है और `--outdir` नाम को `params.yaml` फ़ाइल के पैरामीटर की जगह `stardist_results` कर देती है।
इसके अलावा, तुम देख सकते हो कि `-resume` फ्लैग ने हमें पिछले रन के प्री-प्रोसेसिंग परिणामों का पुन: उपयोग करने की अनुमति दी, जिससे समय की बचत हुई।
तुम इस पैटर्न का उपयोग पाइपलाइन के विभिन्न वेरिएशन को जल्दी से टेस्ट करने के लिए कर सकते हो।

### सारांश

पैरामीटर फ़ाइलें तुम्हारे विश्लेषण को पुनरुत्पादनीय और शेयर करने में आसान बनाती हैं।
किसी भी वास्तविक विश्लेषण कार्य के लिए इनका उपयोग करो।

### आगे क्या है?

सीखो कि सैंपलशीट कई नमूनों के बारे में जानकारी को कैसे व्यवस्थित करती हैं।

---

## 2. सैंपलशीट पैटर्न

### 2.1. सैंपलशीट क्या है?

सैंपलशीट एक CSV फ़ाइल है जो तुम्हारे इनपुट नमूनों का वर्णन करती है।
प्रत्येक पंक्ति एक नमूना है, और कॉलम उस नमूने के लिए फ़ाइलें और मेटाडेटा निर्दिष्ट करते हैं।

आओ उस सैंपलशीट को देखें जिसका हम उपयोग कर रहे हैं:

```bash
cat data/samplesheet.csv
```

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
```

कॉलम हैं:

- `sample`: यूनीक नमूना पहचानकर्ता
- `nuclear_image`: न्यूक्लियर स्टेनिंग इमेज (TIFF)
- `spot_table`: ट्रांसक्रिप्ट स्पॉट (TXT)
- `membrane_image`: मेम्ब्रेन स्टेनिंग इमेज (TIFF, वैकल्पिक)

### 2.2. फ़ाइल पाथ

सैंपलशीट कई पाथ टाइप स्वीकार करती हैं:

- **URL**: Nextflow स्वचालित रूप से डाउनलोड करता है (जैसा ऊपर दिखाया गया है)
- **लोकल पाथ**: `data/nuclear.tiff` या `/absolute/path/to/nuclear.tiff`
- **क्लाउड स्टोरेज**: `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

तुम एक ही सैंपलशीट में पाथ टाइप को मिक्स कर सकते हो।

### 2.3. अपनी खुद की सैंपलशीट बनाना

पहले, आओ टेस्ट डेटा फ़ाइलों को लोकली डाउनलोड करें:

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

अब आओ सैंपलशीट को इन लोकल फ़ाइलों का संदर्भ देने के लिए संशोधित करें:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! warning "चेतावनी"

    ध्यान दो कि सैंपलशीट में पाथ उस जगह के सापेक्ष हैं जहाँ तुम Nextflow **चलाते हो**, न कि जहाँ सैंपलशीट स्थित है।

अंत में, आओ लोकल फ़ाइल पाथ वाली सैंपलशीट के साथ nf-core/molkart को एक बार फिर से एक्जीक्यूट करें:

`nextflow run ./molkart -params-file params.yaml -resume`

जैसा कि तुम देख सकते हो, Nextflow इस रन को उसी तरह एक्जीक्यूट करता है जैसे जब फ़ाइलें Github से डाउनलोड की गई थीं। यह Nextflow की महान विशेषताओं में से एक है, यह तुम्हारे लिए डेटा को ठीक से स्टेज करता है, चाहे वह कहीं भी स्थित हो।

### सारांश

सैंपलशीट मल्टी-सैंपल डेटासेट को इस तरह व्यवस्थित करती हैं जो तुम्हें फ़ाइल पाथ के साथ अपने मेटाडेटा को स्पष्ट रूप से परिभाषित करने की अनुमति देती हैं।
अधिकांश nf-core पाइपलाइन इस पैटर्न का उपयोग करती हैं।

### आगे क्या है?

अब जब हमने इनपुट को कवर कर लिया है, आओ जानें कि विभिन्न कंप्यूटिंग वातावरण के लिए Nextflow पाइपलाइन को कैसे कॉन्फ़िगर करें।
