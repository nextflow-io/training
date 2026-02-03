# भाग 6: Hello Config

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube channel पर [पूरी playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) देखें।

:green_book: वीडियो transcript [यहाँ](./transcripts/06_hello_config.md) उपलब्ध है।
///
-->

यह section explore करेगा कि अपनी Nextflow pipeline का configuration कैसे set up और manage करें ताकि तुम इसके behavior को customize कर सको, इसे different environments में adapt कर सको, और resource usage optimize कर सको _बिना workflow code की single line को alter किए_।

ऐसा करने के multiple ways हैं, जो combination में use किए जा सकते हैं और [यहाँ](https://www.nextflow.io/docs/latest/config.html) described order of precedence के अनुसार interpret किए जाते हैं।

इस course के part में, हम तुम्हें सबसे simple और common configuration file mechanism, `nextflow.config` file दिखाएंगे, जो तुमने Part 5: Hello Containers में पहले ही encounter किया था।

हम Nextflow configuration के essential components जैसे process directives, executors, profiles, और parameter files को cover करेंगे।
इन configuration options को effectively utilize करना सीखकर, तुम अपनी pipelines की flexibility, scalability, और performance enhance कर सकते हो।

??? info "इस section से कैसे शुरू करें"

    Course का यह section मानता है कि तुमने [Hello Nextflow](./index.md) course के Parts 1-5 complete कर लिए हैं और एक complete working pipeline है।

    यदि तुम इस point से course शुरू कर रहे हो, तो तुम्हें solutions से `modules` directory और `nextflow.config` file copy करनी होगी:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    `nextflow.config` file में line `docker.enabled = true` है जो Docker containers का use enable करती है।

    यदि तुम Hello pipeline से familiar नहीं हो या तुम्हें reminder चाहिए, तो [यह info page](../info/hello_pipeline.md) देखो।

---

## 0. Warmup: `hello-config.nf` चलाएं

हम starting point के रूप में workflow script `hello-config.nf` use करेंगे।
यह इस training course के Part 5 में produce की गई script के equivalent है, except हमने output destinations change कर दिए हैं:

```groovy title="hello-config.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    uppercased {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    collected {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    batch_report {
        path 'hello_config'
        mode 'copy'
    }
    cowpy_art {
        path 'hello_config'
        mode 'copy'
    }
}
```

यह sure करने के लिए कि सब कुछ काम कर रहा है, कोई भी changes करने से पहले script को एक बार run करो:

```bash
nextflow run hello-config.nf
```

यदि यह तुम्हारे लिए काम किया, तो तुम अपनी pipelines configure करना सीखने के लिए ready हो।

---

## 1. Workflow input parameters manage करें

हम configuration के एक aspect से शुरू करेंगे जो हम अब तक काम कर रहे थे उसका simply extension है: input parameters का management।

Currently, हमारा workflow command-line के through कई parameter values accept करने के लिए set up है, जिनकी default values workflow script में `params` block में set हैं।
हालाँकि, तुम उन defaults को override करना चाहते हो बिना parameters को command line पर specify करने या original script file modify करने के।

ऐसा करने के multiple ways हैं; हम तुम्हें तीन basic ways दिखाएंगे जो बहुत commonly use होते हैं।

### 1.1. Default values को `nextflow.config` में move करें

यह सबसे simple approach है, हालाँकि यह possibly सबसे least flexible है क्योंकि main `nextflow.config` file कुछ ऐसा नहीं है जिसे तुम हर run के लिए edit करना चाहोगे।
लेकिन इसका advantage है कि यह parameters _declare_ करने (जो definitely workflow में belong करता है) versus _default values_ supply करने की concerns को separate करता है, जो configuration file में अधिक at home हैं।

#### 1.1.1. Configuration file में `params` block create करें

`nextflow.config` file में following code changes करो:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline पैरामीटर
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Note करो कि हमने simply workflow से configuration file में `params` block copy नहीं किया।
Syntax थोड़ा different है।
Workflow file में, वे typed declarations हैं।
Configuration में, वे value assignments हैं।

#### 1.1.2. Workflow file में `params` block से values remove करें

`hello-config.nf` workflow file में following code changes करो:

=== "After"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Pipeline पैरामीटर
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Before"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Pipeline पैरामीटर
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

अब workflow file itself इन parameters के लिए कोई default values set नहीं करती।

#### 1.1.3. Pipeline run करें

```bash
nextflow run hello-config.nf
```

यह अभी भी पहले जैसा same output produce करता है।

### 1.2. Run-specific configuration file use करें

यह great है, लेकिन sometimes तुम main configuration file के साथ mess किए बिना different default values के साथ कुछ temporary experiments run करना चाहते हो।
तुम ऐसा एक subdirectory में new `nextflow.config` file create करके कर सकते हो जिसे तुम अपने experiments के लिए working directory के रूप में use करोगे।

#### 1.2.1. Blank configuration के साथ working directory create करें

एक new directory create करके उसमें move करो:

```bash
mkdir -p tux-run
cd tux-run
```

फिर, उस directory में blank configuration file create करो:

```bash
touch nextflow.config
```

#### 1.2.2. Experimental configuration set up करें

अब new file open करो और जो parameters customize करना चाहते हो वो add करो:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Note करो कि input file का path directory structure reflect करना चाहिए।

#### 1.2.3. Pipeline run करें

```bash
nextflow run ../hello-config.nf
```

इस run में, Nextflow हमारी current directory में `nextflow.config` को pipeline की root directory में `nextflow.config` के साथ combine करता है, और thereby default character (turkey) को tux character के साथ override करता है।

!!! warning "चेतावनी"

    Next section पर move करने से पहले previous directory में वापस change करना sure करो!

    ```bash
    cd ..
    ```

### 1.3. Parameter file use करें

Subdirectory approach experimenting के लिए great काम करता है, लेकिन इसमें थोड़ी setup involve होती है और require होता है कि तुम paths accordingly adapt करो।
एक simpler approach है जब तुम अपनी pipeline को specific set of values के साथ run करना चाहते हो, या किसी और को minimal effort के साथ ऐसा करने enable करना चाहते हो।

Nextflow हमें YAML या JSON format में parameter file के through parameters specify करने allow करता है।

#### 1.3.1. Example parameter file examine करें

इसे demonstrate करने के लिए, हम current directory में एक example parameter file provide करते हैं, जिसका नाम `test-params.yaml` है:

```yaml title="test-params.yaml" linenums="1"
{
  input: "greetings.csv"
  batch: "yaml"
  character: "stegosaurus"
}
```

Note करो कि यदि तुम syntax को configuration file से compare करो तो equal signs (`=`) के बजाय colons (`:`) का use है।
Config file Groovy में लिखी है, जबकि parameter file YAML में लिखी है।

#### 1.3.2. Pipeline run करें

इस parameter file के साथ workflow run करने के लिए, simply base command में `-params-file <filename>` add करो।

```bash
nextflow run hello-config.nf -params-file test-params.yaml
```

Parameter file use करना overkill लग सकता है जब तुम्हारे पास specify करने के लिए only कुछ parameters हों, लेकिन कुछ pipelines दर्जनों parameters expect करती हैं।
उन cases में, parameter file use करना हमें massive command lines type किए बिना और workflow script modify किए बिना runtime पर parameter values provide करने allow करेगा।

### सीख

तुम जानते हो कि workflow inputs manage करने के लिए key configuration options का advantage कैसे लें।

### आगे क्या?

सीखो कि where और how तुम्हारे workflow outputs publish होते हैं यह कैसे manage करें।

---

## 2. Workflow outputs manage करें

अब तक हम workflow-level output declarations के लिए सभी paths hardcode कर रहे थे, और जैसा हमने note किया जब हमने multiple outputs add करना शुरू किया, इसमें थोड़ी repetition involve हो सकती है।

कुछ common ways देखते हैं जिनसे तुम इसे अधिक flexible बनाने के लिए configure कर सकते हो।

### 2.1. `outputDir` directory name customize करें

इस course के हर chapter के लिए, हम outputs को एक different hardcoded subdirectory में publish कर रहे थे।

इसे user-configurable parameter use करने के लिए change करते हैं।
हम इसके लिए एक whole new parameter create कर सकते हैं, लेकिन `batch` parameter use करते हैं क्योंकि वह right there है।

#### 2.1.1. Configuration file में `outputDir` के लिए value set करें

Nextflow जो path outputs publish करने के लिए use करता है वह `outputDir` option द्वारा controlled है।

`nextflow.config` file में following code add करो:

=== "After"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline पैरामीटर
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output सेटिंग्स
    */
    outputDir = "results/${params.batch}"
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline पैरामीटर
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

यह built-in default path, `results/`, को `results/` plus `batch` parameter की value as subdirectory के साथ replace करेगा।

#### 2.1.2. Hardcoded path का repeated part remove करें

हमारे पास अभी भी output options में hardcoded subdirectory है, तो इसे अब remove करते हैं।

Workflow file में following code changes करो:

=== "After"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_config'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_config'
            mode 'copy'
        }
    }
    ```

#### 2.1.3. Pipeline run करें

```bash
nextflow run hello-config.nf --batch outdir
```

यह अभी भी पहले जैसा same output produce करता है, except इस बार हम अपने outputs `results/outdir/` के under पाते हैं।

### 2.2. Process के अनुसार outputs organize करें

Outputs को further organize करने का एक popular way है इसे process के अनुसार करना, _i.e._ pipeline में run होने वाले प्रत्येक process के लिए subdirectories create करना।

#### 2.2.1. Output paths को process names के reference से replace करें

तुम्हें बस output path declaration में process का name `<task>.name` के रूप में reference करना है।

Workflow file में following changes करो:

=== "After"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

#### 2.2.2. Pipeline run करें

```bash
nextflow run hello-config.nf --batch pnames
```

इस बार हम अपने outputs `results/pnames/` के under पाते हैं, और वे process के अनुसार grouped हैं।

### 2.3. Workflow level पर publish mode set करें

Finally, repetitive code की amount reduce करने की spirit में, हम per-output `mode` declarations को configuration में single line से replace कर सकते हैं।

#### 2.3.1. Configuration file में `workflow.output.mode` add करें

`nextflow.config` file में following code add करो:

=== "After"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Output सेटिंग्स
    */
    outputDir = "results/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output सेटिंग्स
    */
    outputDir = "results/${params.batch}"
    ```

#### 2.3.2. Workflow file से output mode remove करें

Workflow file में following changes करो:

=== "After"

    ```groovy title="hello-config.nf" linenums="42"
    output {
        first_output {
            path { sayHello.process }
        }
        uppercased {
            path { convertToUpper.process }
        }
        collected {
            path { collectGreetings.process }
        }
        batch_report {
            path { collectGreetings.process }
        }
        cowpy_art {
            path { cowpy.process }
        }
    }
    ```

=== "Before"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.process }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.process }
            mode 'copy'
        }
        collected {
            path { collectGreetings.process }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.process }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.process }
            mode 'copy'
        }
    }
    ```

यह अधिक concise है।

#### 2.3.3. Pipeline run करें

```bash
nextflow run hello-config.nf --batch outmode
```

यह अभी भी पहले जैसा same output produce करता है, except इस बार हम अपने outputs `results/outmode/` के under पाते हैं।
वे सभी अभी भी proper copies हैं, symlinks नहीं।

### सीख

तुम जानते हो कि directories का naming और structure जहाँ तुम्हारे outputs publish होते हैं, साथ ही workflow output publishing mode कैसे control करें।

### आगे क्या?

सीखो कि अपने workflow configuration को अपने compute environment में कैसे adapt करें, software packaging technology से शुरू करके।

---

## 3. Software packaging technology select करें

अब तक हम configuration elements देख रहे थे जो control करते हैं कि inputs कैसे जाते हैं और outputs कहाँ से आते हैं।
अब specifically अपने workflow configuration को अपने compute environment में adapt करने पर focus करने का time है।

उस path पर पहला step है यह specify करना कि software packages जो प्रत्येक step में run होंगे वे कहाँ से आएंगे।
क्या वे पहले से local compute environment में installed हैं?
क्या हमें images retrieve करनी और उन्हें container system के through run करना है?
या हमें Conda packages retrieve करने और local Conda environment build करना है?

इस training course के बहुत पहले part में (Parts 1-4) हमने अपने workflow में बस locally installed software use किया।
फिर Part 5 में, हमने Docker containers और `nextflow.config` file introduce की, जिसे हमने Docker containers का use enable करने के लिए use किया।

अब देखते हैं कि हम `nextflow.config` file के through एक alternative software packaging option कैसे configure कर सकते हैं।

### 3.1. Config file में Docker disable और Conda enable करें

मान लो हम एक HPC cluster पर काम कर रहे हैं और admin security reasons के लिए Docker का use allow नहीं करता।
Fortunately हमारे लिए, Nextflow कई अन्य container technologies support करता है जिसमें Singularity (जो HPC पर अधिक widely use होती है) शामिल है, और software package managers जैसे Conda।

हम अपनी configuration file को Docker के बजाय Conda use करने के लिए change कर सकते हैं।

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

यह Nextflow को उन processes के लिए Conda environments create और utilize करने allow करेगा जिनके पास Conda packages specified हैं।

### 3.2. Process definition में Conda package specify करें

हम पहले से ही `cowpy` tool contain करने वाले Conda package के लिए URI retrieve कर चुके हैं: `conda-forge::cowpy==1.1.5`

अब हम `conda` directive का use करके URI को `cowpy` process definition में add करते हैं:

=== "After"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "Before"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

Clear करने के लिए, हम `docker` directive _replace_ नहीं कर रहे, हम एक alternative option _add_ कर रहे हैं।

### 3.3. Workflow run करें verify करने के लिए कि यह Conda use कर सकता है

```bash
nextflow run hello-config.nf --batch conda
```

यह बिना issue के काम करना चाहिए और पहले जैसे same outputs `results/conda` के under produce करना चाहिए।

Behind the scenes, Nextflow ने Conda packages retrieve किए और environment create किया।

??? info "Docker और Conda mix और match करना"

    चूंकि ये directives per process assign किए जाते हैं, 'mix और match' करना possible है, _i.e._ अपने workflow में कुछ processes को Docker के साथ और अन्य को Conda के साथ run करने के लिए configure करना, example के लिए, यदि तुम जो compute infrastructure use कर रहे हो वह दोनों support करता है।
    उस case में, तुम अपनी configuration file में Docker और Conda दोनों enable करोगे।
    यदि किसी given process के लिए दोनों available हैं, Nextflow containers को prioritize करेगा।

### सीख

तुम जानते हो कि प्रत्येक process को कौन सा software package use करना चाहिए यह कैसे configure करें, और technologies के बीच switch कैसे करें।

### आगे क्या?

सीखो कि Nextflow द्वारा actually work करने के लिए use किए जाने वाले execution platform को कैसे change करें।

---

## 4. Execution platform select करें

अब तक, हम अपनी pipeline को local executor के साथ run कर रहे थे।
यह प्रत्येक task को उस machine पर execute करता है जिस पर Nextflow run हो रहा है।

Local executor convenient और efficient है, लेकिन यह उस single machine तक limited है।
बहुत large workloads के लिए, तुम discover कर सकते हो कि तुम्हारी local machine bottleneck है।

Nextflow [कई different execution backends](https://www.nextflow.io/docs/latest/executor.html) support करता है, जिसमें HPC schedulers (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor और अन्य) साथ ही cloud execution backends जैसे (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes और अधिक) शामिल हैं।

### 4.1. Different backend target करना

Executor की choice एक process directive द्वारा set होती है जिसे `executor` कहते हैं।
By default यह `local` पर set है, तो following configuration implied है:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Different backend target करने के लिए executor set करने के लिए, तुम simply वह executor specify करोगे जो तुम चाहते हो।

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "चेतावनी"

    हम actually training environment में इसे test नहीं कर सकते क्योंकि यह HPC से connect करने के लिए set up नहीं है।

### 4.2. Execution parameters के लिए backend-specific syntax deal करना

Most high-performance computing platforms allow (और sometimes require) करते हैं कि तुम certain parameters specify करो जैसे resource allocation requests और limitations (e.g. number of CPUs और memory) और use करने के लिए job queue का name।

Unfortunately, इनमें से प्रत्येक system different technologies, syntaxes और configurations use करता है यह define करने के लिए कि job कैसे define और relevant scheduler को submit किया जाना चाहिए।

Fortunately, Nextflow यह सब simplify करता है।
यह एक standardized syntax provide करता है ताकि तुम relevant properties जैसे `cpus`, `memory` और `queue` सिर्फ एक बार specify कर सको।
फिर, runtime पर, Nextflow उन settings को executor setting के based पर appropriate backend-specific scripts generate करने के लिए use करेगा।

### सीख

अब तुम जानते हो कि different kinds of computing infrastructure use करने के लिए executor कैसे change करें।

### आगे क्या?

सीखो कि Nextflow में resource allocations और limitations कैसे evaluate और express करें।

---

## 5. Compute resource allocations control करें

Most high-performance computing platforms allow (और sometimes require) करते हैं कि तुम certain resource allocation parameters specify करो जैसे number of CPUs और memory।

By default, Nextflow प्रत्येक process के लिए single CPU और 2GB memory use करेगा।
Corresponding process directives को `cpus` और `memory` कहा जाता है, तो following configuration implied है:

```groovy title="Built-in configuration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

तुम इन values को modify कर सकते हो, either सभी processes के लिए या specific named processes के लिए, अपनी configuration file में additional process directives use करके।

लेकिन तुम कैसे जानते हो कि कौन सी values use करनी हैं?

### 5.1. Resource utilization report generate करने के लिए workflow run करें

यदि तुम up front नहीं जानते कि तुम्हारे processes को कितनी CPU और memory की likely need होगी, तुम कुछ resource profiling कर सकते हो।

Conveniently, Nextflow में इसके लिए built-in tools included हैं, और request पर तुम्हारे लिए report happily generate करेगा।

ऐसा करने के लिए, अपनी command line में `-with-report <filename>.html` add करो।

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

Report एक html file है, जिसे तुम download करके अपने browser में open कर सकते हो।

Report को देखने के लिए कुछ minutes लो और identify करो कि resources adjust करने के लिए कुछ opportunities हैं या नहीं।

### 5.2. सभी processes के लिए resource allocations set करें

Profiling show करती है कि हमारी training workflow में processes बहुत lightweight हैं, तो default memory allocation को 1GB per process तक reduce करते हैं।

अपनी `nextflow.config` file में following add करो:

```groovy title="nextflow.config" linenums="4"
/*
* Process सेटिंग्स
*/
process {
    memory = 1.GB
}
```

### 5.3. Specific process के लिए resource allocations set करें

साथ ही, हम pretend करेंगे कि `cowpy` process को दूसरों से अधिक resources require होती हैं, बस demonstrate करने के लिए कि individual process के लिए allocations कैसे adjust करें।

=== "After"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }
    ```

इस configuration के साथ, सभी processes 1GB memory और single CPU (implied default) request करेंगे, except `cowpy` process, जो 2GB और 2 CPUs request करेगा।

### 5.4. Updated configuration के साथ workflow run करें

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

तुम probably कोई real difference notice नहीं करोगे चूंकि यह इतना small workload है, लेकिन यह approach है जो तुम real-world workflow की performance और resource requirements analyze करने के लिए use करोगे।

### 5.5. Resource limits add करें

Depending on तुम कौन सा computing executor और compute infrastructure use कर रहे हो, कुछ constraints हो सकते हैं कि तुम क्या allocate कर सकते हो (या must)।

तुम `resourceLimits` directive use कर सकते हो relevant limitations set करने के लिए:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow इन values को उस executor के depending appropriate instructions में translate करेगा जो तुमने specify किया।

### सीख

तुम जानते हो कि resource utilization assess करने के लिए profiling report कैसे generate करें और सभी processes के लिए और/या individual processes के लिए resource allocations कैसे modify करें, साथ ही HPC पर running के लिए resource limitations set करें।

### आगे क्या?

सीखो कि preset configuration profiles कैसे set up करें और runtime पर उनके बीच switch करें।

---

## 6. Preset configurations के बीच switch करने के लिए profiles use करें

हमने तुम्हें कई ways दिखाए हैं जिनसे तुम अपनी pipeline configuration customize कर सकते हो depending on तुम किस project पर काम कर रहे हो या तुम कौन सा compute environment use कर रहे हो।

तुम alternative settings के बीच switch करना चाहते हो depending on तुम कौन सी computing infrastructure use कर रहे हो।
Example के लिए, तुम अपने laptop पर locally develop और small-scale tests run करना चाहते हो, फिर HPC या cloud पर full-scale workloads run करना चाहते हो।

Nextflow तुम्हें कितनी भी profiles set up करने देता है जो different configurations describe करती हैं, जिन्हें तुम फिर runtime पर command-line argument use करके select कर सकते हो, बजाय configuration file itself modify करने के।

### 6.1. Local development और HPC पर execution के बीच switch करने के लिए profiles create करें

दो alternative profiles set up करते हैं; एक regular computer पर small scale loads run करने के लिए, जहाँ हम Docker containers use करेंगे, और एक Slurm scheduler के साथ university HPC पर running के लिए, जहाँ हम Conda packages use करेंगे।

#### 6.1.1. Profiles set up करें

अपनी `nextflow.config` file में following add करो:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
}
```

तुम देखते हो कि university HPC के लिए, हम resource limitations भी specify कर रहे हैं।

#### 6.1.2. Profile के साथ workflow run करें

अपनी Nextflow command line में profile specify करने के लिए, हम `-profile` argument use करते हैं।

`my_laptop` configuration के साथ workflow run करने की try करते हैं।

```bash
nextflow run hello-config.nf -profile my_laptop
```

जैसा तुम देख सकते हो, यह हमें runtime पर configurations के बीच बहुत conveniently toggle करने allow करता है।

!!! warning "चेतावनी"

    `univ_hpc` profile training environment में properly run नहीं होगी चूंकि हमारे पास Slurm scheduler तक access नहीं है।

### 6.2. Test parameters की profile create करें

Profiles सिर्फ infrastructure configuration के लिए नहीं हैं।
हम उन्हें workflow parameters के लिए default values set करने के लिए भी use कर सकते हैं, ताकि दूसरों के लिए workflow को try out करना easier हो बिना appropriate input values खुद gather किए।

#### 6.2.1. Profile set up करें

यदि हम अपने workflow के लिए test profile add करें, तो `profiles` block बन जाता है:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
    test {
        params.greeting = 'greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

#### 6.2.2. Test profile के साथ workflow locally run करें

Conveniently, profiles mutually exclusive नहीं हैं, तो हम following syntax `-profile <profile1>,<profile2>` use करके अपनी command line में multiple profiles specify कर सकते हैं।

अपने previous command में test profile add करने की try करते हैं:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

यह Docker use करेगा जहाँ possible हो और `results/test` के under outputs produce करेगा, और इस बार character comedic duo `dragonandcow` है।

### 6.3. Resolved configuration देखने के लिए `nextflow config` use करें

जैसा ऊपर noted है, sometimes same parameter profiles में different values पर set हो सकता है जिन्हें तुम combine करना चाहते हो।
और more generally, कई places हैं जहाँ configuration के elements stored हो सकते हैं, और sometimes same properties different places में different values पर set हो सकती हैं।

Nextflow किसी भी conflicts को resolve करने के लिए set [order of precedence](https://www.nextflow.io/docs/latest/config.html) apply करता है।

Fortunately, Nextflow में एक convenient utility tool included है जिसे `config` कहते हैं जो तुम्हारे लिए वह whole process automate कर सकता है।

#### 6.3.1. Default configuration resolve करें

यह command run करो default द्वारा apply होने वाली configuration resolve करने के लिए।

```bash
nextflow config
```

यह तुम्हें base configuration show करता है जो तुम्हें मिलती है यदि तुम command line में कुछ extra specify नहीं करते।

#### 6.3.2. Specific settings activated के साथ configuration resolve करें

यदि तुम command-line parameters provide करते हो, e.g. एक या अधिक profiles enable करना या parameter file load करना, command additionally उन्हें account में लेगा।

```bash
nextflow config -profile my_laptop,test
```

यह complex projects के लिए especially useful हो जाता है जिनमें configuration की multiple layers involve होती हैं।

### सीख

तुम जानते हो कि minimal hassle के साथ runtime पर preset configuration select करने के लिए profiles कैसे use करें।
More generally, तुम जानते हो कि अपने workflow executions को different compute platforms suit करने के लिए कैसे configure करें और अपनी analyses की reproducibility enhance करें।

### आगे क्या?

Celebrate करो और खुद को एक big pat on the back दो! तुमने अपना बहुत पहला Nextflow developer course complete कर लिया है।

Final [course summary](./next_steps.md) पर जाओ review करने के लिए कि तुमने क्या सीखा और पता लगाओ कि आगे क्या आता है।

---

## Quiz

<quiz>
उस configuration file का क्या नाम है जो Nextflow automatically load करता है?
- [ ] `config.nf`
- [ ] `pipeline.config`
- [x] `nextflow.config`
- [ ] `workflow.config`
</quiz>

<quiz>
जब same parameter config file और command line दोनों में set हो तो क्या precedence लेता है?
- [ ] Config file value
- [x] Command line value
- [ ] पहली encountered value
- [ ] Neither; यह error cause करता है

और जानें: [1.1. Default values को `nextflow.config` में move करें](#11-default-values-को-nextflowconfig-में-move-करें)
</quiz>

<quiz>
क्या same configuration में Docker और Conda दोनों enabled हो सकते हैं?
- [x] हाँ, Nextflow process directives के depending दोनों use कर सकता है
- [ ] नहीं, एक time पर सिर्फ एक enabled हो सकता है
- [ ] हाँ, लेकिन सिर्फ profiles में
- [ ] नहीं, वे mutually exclusive हैं
</quiz>

<quiz>
यदि Docker और Conda दोनों enabled हैं और process के पास दोनों directives हैं, तो कौन prioritized होता है?
- [x] Docker (containers)
- [ ] Conda
- [ ] पहला defined
- [ ] यह error cause करता है

और जानें: [3. Software packaging technology select करें](#3-software-packaging-technology-select-करें)
</quiz>

<quiz>
Nextflow processes के लिए default memory allocation क्या है?
- [ ] 1 GB
- [x] 2 GB
- [ ] 4 GB
- [ ] No limit
</quiz>

<quiz>
Config file में specific process के लिए resource requirements कैसे set करते हो?
- [ ] `#!groovy processName.memory = '4 GB'`
- [ ] `#!groovy process.memory.processName = '4 GB'`
- [x] `#!groovy process { withName: 'processName' { memory = '4 GB' } }`
- [ ] `#!groovy resources.processName.memory = '4 GB'`

और जानें: [5.3. Specific process के लिए resource allocations set करें](#53-specific-process-के-लिए-resource-allocations-set-करें)
</quiz>

<quiz>
कौन सा command line option resource utilization report generate करता है?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

और जानें: [5.1. Resource utilization report generate करने के लिए workflow run करें](#51-resource-utilization-report-generate-करने-के-लिए-workflow-run-करें)
</quiz>

<quiz>
`resourceLimits` directive क्या करता है?
- [ ] Minimum resource requirements set करता है
- [ ] Processes को resources allocate करता है
- [x] Maximum resources को cap करता है जो request किए जा सकते हैं
- [ ] Resource usage monitor करता है

और जानें: [5.5. Resource limits add करें](#55-resource-limits-add-करें)
</quiz>

<quiz>
Nextflow में default executor क्या है?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

और जानें: [4. Execution platform select करें](#4-execution-platform-select-करें)
</quiz>

<quiz>
Nextflow run करते समय parameter file कैसे specify करते हो?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

और जानें: [1.3. Parameter file use करें](#13-parameter-file-use-करें)
</quiz>

<quiz>
Profiles किसके लिए use की जा सकती हैं? (सभी लागू select करें)
- [x] Infrastructure-specific settings define करने के लिए
- [x] Different environments के लिए resource limits set करने के लिए
- [x] Test parameters provide करने के लिए
- [ ] New processes define करने के लिए

और जानें: [6. Preset configurations के बीच switch करने के लिए profiles use करें](#6-preset-configurations-के-बीच-switch-करने-के-लिए-profiles-use-करें)
</quiz>

<quiz>
Single command में multiple profiles कैसे specify करते हो?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

और जानें: [6. Preset configurations के बीच switch करने के लिए profiles use करें](#6-preset-configurations-के-बीच-switch-करने-के-लिए-profiles-use-करें)
</quiz>
