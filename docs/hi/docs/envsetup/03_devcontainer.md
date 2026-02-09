# लोकल Devcontainers

यदि तुम्हारे पास लोकल Docker इंस्टॉलेशन है या तुम एक इंस्टॉल करने के लिए तैयार हो, तो इन सामग्रियों के साथ लोकल रूप से काम करने का सबसे आसान तरीका Visual Studio Code की devcontainer सुविधा का उपयोग करना है। यह दृष्टिकोण मैनुअल इंस्टॉलेशन की आवश्यकता के बिना सभी आवश्यक टूल्स और डिपेंडेंसीज़ प्रदान करता है।

## आवश्यकताएं

लोकल devcontainer सेटअप का उपयोग करने के लिए, तुम्हें चाहिए:

- [Visual Studio Code](https://code.visualstudio.com/)
- एक लोकल Docker इंस्टॉलेशन, उदाहरण के लिए:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (Windows/macOS के लिए)
  - [Docker Engine](https://docs.docker.com/engine/install/) (Linux के लिए)
  - [Colima](https://github.com/abiosoft/colima) (macOS के लिए विकल्प)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (Docker Desktop में शामिल है, लेकिन अन्य Docker सेटअप के साथ अलग से इंस्टॉलेशन की आवश्यकता हो सकती है)
- VS Code के लिए [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

devcontainer खोलने का प्रयास करने से पहले तुम्हारा Docker इंस्टॉलेशन चल रहा होना चाहिए।

यह सत्यापित करने के लिए कि Docker buildx उपलब्ध है, चलाओ:

```bash
docker buildx version
```

यदि यह कमांड विफल होता है, तो आगे बढ़ने से पहले तुम्हें buildx एक्सटेंशन इंस्टॉल करना होगा।

## सेटअप निर्देश

VS Code devcontainers का उपयोग करके अपना लोकल वातावरण सेट अप करने के लिए इन चरणों का पालन करो:

### VS Code में "Dev Containers" एक्सटेंशन इंस्टॉल करो

- VS Code खोलो
- Extensions पर जाओ (Ctrl+Shift+X या macOS पर Cmd+Shift+X)
- "Dev Containers" खोजो
- "Install" पर क्लिक करो

![VS Code में Dev Containers एक्सटेंशन इंस्टॉल करना](img/install_extension.png)

### रिपॉज़िटरी क्लोन करो:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### रिपॉज़िटरी को VS Code में खोलो:

- VS Code लॉन्च करो
- मेनू से **File -> Open Folder** चुनो
- अभी क्लोन किए गए training रिपॉज़िटरी फ़ोल्डर को नेविगेट करके चुनो
- **Open** पर क्लिक करो

### Container में फिर से खोलो

यदि VS Code द्वारा "Reopen in Container" के लिए प्रॉम्प्ट किया जाता है, तो उस पर क्लिक करो। वैकल्पिक रूप से:

- F1 दबाओ (या Ctrl+Shift+P / macOS पर Cmd+Shift+P)
- "Dev Containers: Reopen in Container" टाइप करो
- **महत्वपूर्ण**: जब कॉन्फ़िगरेशन चुनने के लिए प्रॉम्प्ट किया जाए, तो **local-dev** devcontainer कॉन्फ़िगरेशन चुनो

![Reopen in Container प्रॉम्प्ट](img/reopen_prompt.png)

![लोकल कॉन्फ़िगरेशन चुनना](img/select_local_config.png)

कंटेनर के बिल्ड होने की प्रतीक्षा करो। पहली बार यह कुछ मिनट ले सकता है क्योंकि यह सभी आवश्यक कंपोनेंट्स को डाउनलोड और सेट अप करता है।

एक बार कंटेनर बिल्ड और चल रहा हो, तुम्हारे पास सभी आवश्यक टूल्स के साथ एक पूरी तरह से कॉन्फ़िगर किया गया वातावरण होगा, जिसमें शामिल हैं:

- Java
- Nextflow
- Docker
- Git
- और प्रशिक्षण के लिए आवश्यक सभी अन्य डिपेंडेंसीज़

![devcontainer चलाते हुए VS Code](img/running_container.png)

## Devcontainers का उपयोग करने के लाभ

devcontainer दृष्टिकोण का उपयोग करने से कई फायदे मिलते हैं:

- **स्थिरता**: विभिन्न मशीनों में एक सुसंगत विकास वातावरण सुनिश्चित करता है
- **सरलता**: सभी डिपेंडेंसीज़ पहले से इंस्टॉल और कॉन्फ़िगर की गई हैं
- **अलगाव**: विकास वातावरण तुम्हारे लोकल सिस्टम से अलग है
- **पुनरुत्पादनीयता**: devcontainer का उपयोग करने वाले सभी को समान सेटअप मिलता है
- **कोई मैनुअल इंस्टॉलेशन नहीं**: Java, Nextflow और अन्य टूल्स को मैनुअल रूप से इंस्टॉल करने की आवश्यकता नहीं

## अपने वातावरण की जांच करना

एक बार तुम्हारा devcontainer चल रहा हो, तुम यह सत्यापित कर सकते हो कि सब कुछ सही ढंग से सेट अप है:

```bash
nextflow info
```

यह Nextflow संस्करण और रनटाइम जानकारी प्रदर्शित करेगा, जो पुष्टि करता है कि तुम्हारा वातावरण ठीक से कॉन्फ़िगर किया गया है।

## समस्या निवारण

यदि तुम्हें devcontainer सेटअप के साथ समस्याओं का सामना करना पड़ता है:

1. सुनिश्चित करो कि devcontainer खोलने से पहले तुम्हारा Docker इंस्टॉलेशन (Docker Desktop, Colima, Docker Engine, आदि) चल रहा है
2. जांचो कि तुमने प्रॉम्प्ट किए जाने पर **local-dev** कॉन्फ़िगरेशन चुना है
3. सत्यापित करो कि Docker buildx इंस्टॉल है और `docker buildx version` चलाकर काम कर रहा है
4. यदि कंटेनर बिल्ड करने में विफल रहता है, तो "Dev Containers: Rebuild Container" कमांड चलाकर इसे फिर से बिल्ड करने का प्रयास करो
5. लगातार समस्याओं के लिए, [VS Code Dev Containers troubleshooting guide](https://code.visualstudio.com/docs/devcontainers/troubleshooting) देखो
