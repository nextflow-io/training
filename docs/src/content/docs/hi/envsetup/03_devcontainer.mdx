# Local Devcontainers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

अगर तुम्हारे पास local Docker installation है या install करने के लिए तैयार हो, तो इन materials के साथ स्थानीय रूप से काम करने का सबसे आसान तरीका Visual Studio Code की devcontainer feature का उपयोग करना है। यह approach मैन्युअल installation की आवश्यकता के बिना सभी आवश्यक tools और dependencies प्रदान करता है।

## आवश्यकताएँ

Local devcontainer setup का उपयोग करने के लिए, तुम्हें चाहिए:

- [Visual Studio Code](https://code.visualstudio.com/)
- एक local Docker installation, उदाहरण के लिए:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (Windows/macOS के लिए)
  - [Docker Engine](https://docs.docker.com/engine/install/) (Linux के लिए)
  - [Colima](https://github.com/abiosoft/colima) (macOS के लिए वैकल्पिक)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (Docker Desktop में शामिल, लेकिन अन्य Docker setups के साथ अलग installation की आवश्यकता हो सकती है)
- VS Code के लिए [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

Devcontainer खोलने का प्रयास करने से पहले तुम्हारी Docker installation चल रही होनी चाहिए।

यह verify करने के लिए कि Docker buildx उपलब्ध है, चलाओ:

```bash
docker buildx version
```

अगर यह command fail होता है, तो आगे बढ़ने से पहले तुम्हें buildx extension install करना होगा।

## Setup निर्देश

VS Code devcontainers का उपयोग करके अपना local वातावरण सेट करने के लिए इन steps का पालन करो:

### VS Code में "Dev Containers" extension install करो

- VS Code खोलो
- Extensions पर जाओ (Ctrl+Shift+X या macOS पर Cmd+Shift+X)
- "Dev Containers" search करो
- "Install" पर click करो

![Installing Dev Containers extension in VS Code](img/install_extension.png)

### Repository clone करो:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### VS Code में repository खोलो:

- VS Code launch करो
- Menu से **File -> Open Folder** चुनो
- तुमने अभी clone की training repository folder पर navigate करो और उसे चुनो
- **Open** पर click करो

### Container में फिर से खोलो

अगर VS Code "Reopen in Container" के लिए prompt करे, तो उस पर click करो। वैकल्पिक रूप से:

- F1 (या Ctrl+Shift+P / macOS पर Cmd+Shift+P) दबाओ
- "Dev Containers: Reopen in Container" type करो
- **महत्वपूर्ण**: जब configuration चुनने के लिए prompt हो, **local-dev** devcontainer configuration चुनो

![Reopen in Container prompt](img/reopen_prompt.png)

![Selecting local configuration](img/select_local_config.png)

Container build होने का इंतज़ार करो। पहली बार इसमें कुछ मिनट लग सकते हैं क्योंकि यह सभी आवश्यक components download और सेट करता है।

एक बार container build और run हो जाए, तुम्हारे पास सभी आवश्यक tools installed के साथ पूर्णतः configured वातावरण होगा, जिसमें शामिल हैं:

- Java
- Nextflow
- Docker
- Git
- और प्रशिक्षण के लिए आवश्यक अन्य सभी dependencies

![VS Code with devcontainer running](img/running_container.png)

## Devcontainers के उपयोग के लाभ

Devcontainer approach का उपयोग करने के कई फायदे हैं:

- **सुसंगतता**: विभिन्न machines पर एक सुसंगत development वातावरण सुनिश्चित करता है
- **सरलता**: सभी dependencies पूर्व-installed और configured हैं
- **अलगाव**: Development वातावरण तुम्हारे local system से अलग है
- **पुनरुत्पादनीयता**: Devcontainer का उपयोग करने वाले सभी को समान setup मिलता है
- **कोई मैन्युअल installation नहीं**: Java, Nextflow और अन्य tools मैन्युअल रूप से install करने की आवश्यकता नहीं

## अपने वातावरण की जाँच करना

एक बार तुम्हारा devcontainer चल रहा हो, तुम verify कर सकते हो कि सब कुछ सही ढंग से सेट है:

```bash
nextflow info
```

यह Nextflow version और runtime जानकारी दिखाएगा, जो confirm करेगा कि तुम्हारा वातावरण ठीक से configured है।

## समस्या निवारण

अगर तुम्हें devcontainer setup में समस्याएँ आती हैं:

1. Devcontainer खोलने से पहले सुनिश्चित करो कि तुम्हारी Docker installation (Docker Desktop, Colima, Docker Engine, आदि) चल रही है
2. जाँचो कि prompt होने पर तुमने **local-dev** configuration चुना है
3. `docker buildx version` चलाकर verify करो कि Docker buildx installed और काम कर रहा है
4. अगर container build fail हो जाए, तो "Dev Containers: Rebuild Container" command चलाकर इसे फिर से build करने का प्रयास करो
5. लगातार समस्याओं के लिए, [VS Code Dev Containers troubleshooting guide](https://code.visualstudio.com/docs/devcontainers/troubleshooting) देखो
