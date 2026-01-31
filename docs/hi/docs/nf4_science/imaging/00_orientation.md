# ओरिएंटेशन

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

यह ओरिएंटेशन मानता है कि आपने पहले से ही "Open in GitHub Codespaces" बटन पर क्लिक करके प्रशिक्षण वातावरण खोल लिया है।
यदि नहीं, तो कृपया अभी ऐसा करें, आदर्श रूप से दूसरे ब्राउज़र विंडो या टैब में ताकि आप इन निर्देशों को वापस देख सकें।

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "मशीन साइज़ की आवश्यकता"

    इस प्रशिक्षण कोर्स के लिए अपना Codespace बनाते समय **8-core मशीन** चुनना सुनिश्चित करें। bioimaging workflows के लिए अतिरिक्त कंप्यूट संसाधनों की आवश्यकता होती है।

## GitHub Codespaces

GitHub Codespaces वातावरण में इस प्रशिक्षण कोर्स के माध्यम से काम करने के लिए आवश्यक सभी सॉफ़्टवेयर, कोड और डेटा शामिल हैं, इसलिए आपको स्वयं कुछ भी इंस्टॉल करने की आवश्यकता नहीं है।
हालांकि, लॉग इन करने के लिए आपको एक (मुफ्त) GitHub अकाउंट की आवश्यकता है, और यदि आप इंटरफ़ेस से अपरिचित हैं तो आपको [GitHub Codespaces Orientation](../../envsetup/index.md) मिनी-कोर्स पूरा करके इसे जानने में कुछ मिनट देने चाहिए।

## Docker images को पहले से डाउनलोड करें

एक बार जब आप अपना Codespace खोल लें, तो आइए इस प्रशिक्षण कोर्स के लिए आवश्यक सभी Docker images को पहले से डाउनलोड करें।
यह बाद में समय बचाएगा और workflows के सुचारू निष्पादन को सुनिश्चित करेगा।

एक नया टर्मिनल टैब खोलें और निम्नलिखित कमांड चलाएं:

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

यह कमांड पृष्ठभूमि में सभी आवश्यक Docker images डाउनलोड करेगा।
जब यह चल रहा हो तब आप ओरिएंटेशन के बाकी हिस्से को जारी रख सकते हैं।

!!!tip "सुझाव"

    `-stub` फ्लैग pipeline को वास्तविक डेटा को प्रोसेस किए बिना तेज़ी से चलाने की अनुमति देता है, जो images डाउनलोड करने के लिए बिल्कुल सही है। आप टर्मिनल टैब में प्रगति की निगरानी कर सकते हैं।

## वर्किंग डायरेक्टरी

इस प्रशिक्षण कोर्स के दौरान, हम `nf4-science/imaging/` डायरेक्टरी में काम करेंगे।

टर्मिनल में यह कमांड चलाकर अभी डायरेक्टरी बदलें:

```bash
cd nf4-science/imaging/
```

!!!tip "सुझाव"

    यदि किसी कारणवश आप इस डायरेक्टरी से बाहर चले जाते हैं, तो आप हमेशा इसमें वापस आने के लिए पूरा path उपयोग कर सकते हैं, यह मानते हुए कि आप इसे GitHub Codespaces प्रशिक्षण वातावरण के भीतर चला रहे हैं:

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**अब, कोर्स शुरू करने के लिए, इस पेज के नीचे दाएं कोने में तीर पर क्लिक करें।**
