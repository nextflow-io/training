# भाग 3: इनपुट को व्यवस्थित करना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

भाग 2 में, हमने कमांड लाइन पर कई पैरामीटर के साथ molkart चलाया।
अब हम इनपुट को प्रबंधित करने के दो बेहतर तरीके सीखेंगे: **पैरामीटर फ़ाइलें** और **samplesheets**।

## 1. पैरामीटर फ़ाइलों का उपयोग

### 1.1. लंबी कमांड लाइनों की समस्या

भाग 2 से हमारी कमांड याद करें:

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

यह काम करता है, लेकिन इसे पुनः उत्पन्न करना, साझा करना या संशोधित करना कठिन है।
यदि आपको अगले महीने फिर से वही विश्लेषण चलाने की आवश्यकता हो तो क्या होगा?
यदि कोई सहयोगी आपकी सटीक सेटिंग्स का उपयोग करना चाहता है तो क्या होगा?

### 1.2. समाधान: पैरामीटर फ़ाइल का उपयोग करें

`params.yaml` नाम की एक फ़ाइल बनाएं:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

अब आपकी कमांड बन जाती है:

```bash
nextflow run ./molkart -params-file params.yaml -resume
```

बस इतना ही! पैरामीटर फ़ाइल आपके सटीक कॉन्फ़िगरेशन को दस्तावेज़ित करती है और इसे फिर से चलाना या साझा करना आसान बनाती है।

### 1.3. पैरामीटर को ओवरराइड करना

आप अभी भी कमांड लाइन से विशिष्ट पैरामीटर को ओवरराइड कर सकते हैं:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

उपरोक्त लाइन `segmentation_method` को `stardist` में बदल देती है और `--outdir` नाम को `params.yaml` फ़ाइल में दिए गए पैरामीटर के बजाय `stardist_results` में बदल देती है।
इसके अतिरिक्त, आप देख सकते हैं कि `-resume` फ़्लैग ने हमें पिछले रन से प्री-प्रोसेसिंग परिणामों का पुन: उपयोग करने की अनुमति दी, जिससे समय की बचत हुई।
आप इस पैटर्न का उपयोग pipeline के विभिन्न विविधताओं का शीघ्रता से परीक्षण करने के लिए कर सकते हैं।

### मुख्य बात

पैरामीटर फ़ाइलें आपके विश्लेषण को पुनरुत्पादनीय और साझा करने में आसान बनाती हैं।
किसी भी वास्तविक विश्लेषण कार्य के लिए इनका उपयोग करें।

### आगे क्या?

जानें कि samplesheets कैसे कई नमूनों के बारे में जानकारी को व्यवस्थित करती हैं।

---

## 2. Samplesheet पैटर्न

### 2.1. Samplesheet क्या है?

Samplesheet एक CSV फ़ाइल है जो आपके इनपुट नमूनों का वर्णन करती है।
प्रत्येक पंक्ति एक नमूना है, और कॉलम उस नमूने के लिए फ़ाइलों और मेटाडेटा को निर्दिष्ट करते हैं।

आइए देखें कि हम किस samplesheet का उपयोग कर रहे हैं:

```bash
cat data/samplesheet.csv
```

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
```

कॉलम हैं:

- `sample`: अद्वितीय नमूना पहचानकर्ता
- `nuclear_image`: परमाणु स्टेनिंग इमेज (TIFF)
- `spot_table`: ट्रांसक्रिप्ट स्पॉट्स (TXT)
- `membrane_image`: मेम्ब्रेन स्टेनिंग इमेज (TIFF, वैकल्पिक)

### 2.2. फ़ाइल पथ

Samplesheets कई पथ प्रकारों को स्वीकार करती हैं:

- **URLs**: Nextflow स्वचालित रूप से डाउनलोड करता है (जैसा कि ऊपर दिखाया गया है)
- **स्थानीय पथ**: `data/nuclear.tiff` या `/absolute/path/to/nuclear.tiff`
- **Cloud storage**: `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

आप एक ही samplesheet में पथ प्रकारों को मिक्स कर सकते हैं।

### 2.3. अपनी खुद की samplesheet बनाना

सबसे पहले, आइए परीक्षण डेटा फ़ाइलों को स्थानीय रूप से डाउनलोड करें:

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

अब आइए इन स्थानीय फ़ाइलों का संदर्भ देने के लिए samplesheet को संशोधित करें:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! warning "चेतावनी"

    ध्यान दें कि samplesheet में पथ वहाँ से सापेक्ष हैं जहाँ आप Nextflow **चलाते** हैं, न कि जहाँ samplesheet स्थित है।

अंत में, आइए स्थानीय फ़ाइल पथों के साथ samplesheet के साथ nf-core/molkart को एक बार फिर से निष्पादित करें:

`nextflow run ./molkart -params-file params.yaml -resume`

जैसा कि आप देख सकते हैं, Nextflow इस रन को उसी तरह निष्पादित करता है जैसे कि फ़ाइलें Github से डाउनलोड की गई थीं। यह Nextflow की महान विशेषताओं में से एक है, यह आपके लिए डेटा को ठीक से व्यवस्थित करता है, चाहे वह कहीं भी स्थित हो।

### मुख्य बात

Samplesheets बहु-नमूना डेटासेट को इस तरह व्यवस्थित करती हैं जो आपको फ़ाइल पथों के साथ अपने मेटाडेटा को स्पष्ट रूप से परिभाषित करने की अनुमति देती हैं।
अधिकांश nf-core pipelines इस पैटर्न का उपयोग करती हैं।

### आगे क्या?

अब जब हमने इनपुट को कवर कर लिया है, तो आइए जानें कि विभिन्न कंप्यूटिंग वातावरण के लिए Nextflow pipelines को कैसे कॉन्फ़िगर किया जाए।
