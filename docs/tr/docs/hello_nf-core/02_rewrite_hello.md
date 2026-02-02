# Bölüm 2: Hello'yu nf-core için yeniden yazma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu Hello nf-core eğitim kursunun ikinci bölümünde, [Hello Nextflow](../hello_nextflow/index.md) başlangıç kursu tarafından üretilen pipeline'ın nf-core uyumlu bir versiyonunu nasıl oluşturacağınızı gösteriyoruz.

Eğitimin ilk bölümünde nf-core pipeline'larının çok sayıda yardımcı dosya ile oldukça detaylı bir yapıyı takip ettiğini fark etmişsinizdir.
Tüm bunları sıfırdan oluşturmak çok yorucu olurdu, bu yüzden nf-core topluluğu, süreci başlatmak için bunun yerine bir şablondan yapılmasını sağlayan araçlar geliştirmiştir.

Bu araçları kullanarak bir pipeline iskeleti oluşturmayı, ardından mevcut 'normal' pipeline kodunu nf-core iskeleti üzerine nasıl uyarlayacağınızı göstereceğiz.

Eğer Hello pipeline'ına aşina değilseniz veya hatırlatmaya ihtiyacınız varsa, [bu bilgi sayfasına](../info/hello_pipeline.md) bakın.

---

## 1. Yeni bir pipeline projesi oluşturma

İlk olarak, yeni pipeline için iskeleyi oluşturuyoruz.

!!! note

    Terminalinizde `hello-nf-core` dizininde olduğunuzdan emin olun.

### 1.1. Şablon tabanlı pipeline oluşturma aracını çalıştırma

Yeni bir pipeline oluşturmak için `nf-core pipelines create` komutu ile başlayalım.
Bu, nf-core temel şablonunu kullanarak, bir pipeline adı, açıklaması ve yazarı ile özelleştirilmiş yeni bir pipeline iskeleti oluşturacaktır.

```bash
nf-core pipelines create
```

Bu komutu çalıştırmak, pipeline oluşturma için bir Metin Kullanıcı Arayüzü (TUI) açacaktır:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/VwjXNXONHlY?si=d0HkFSISnKn76TeI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

Bu TUI, pipeline'ınız hakkında temel bilgiler sağlamanızı isteyecek ve pipeline iskelenize dahil edilecek veya hariç tutulacak özellikler seçeneği sunacaktır.

- Hoş geldiniz ekranında **Let's go!** butonuna tıklayın.
- `Choose pipeline type` ekranında **Custom**'a tıklayın.
- Pipeline bilgilerinizi aşağıdaki gibi girin (`< YOUR NAME >` kısmını kendi adınızla değiştirerek), ardından **Next**'e tıklayın.

```
[ ] GitHub organisation: core
[ ] Workflow name: hello
[ ] A short description of your pipeline: A basic nf-core style version of Hello Nextflow
[ ] Name of the main author(s): < YOUR NAME >
```

- Template features ekranında, `Toggle all features`'ı **kapalı** olarak ayarlayın, ardından aşağıdakileri seçerek **etkinleştirin**. Seçimlerinizi kontrol edin ve **Continue**'ye tıklayın.

```
[ ] Add testing profiles
[ ] Use nf-core components
[ ] Use nf-schema
[ ] Add configuration files
[ ] Add documentation
```

- `Final details` ekranında **Finish**'e tıklayın. Pipeline'ın oluşturulmasını bekleyin, ardından **Continue**'ye tıklayın.
- Create GitHub repository ekranında **Finish without creating a repo**'ya tıklayın. Bu, daha sonra bir GitHub repository'si oluşturmak için talimatları görüntüleyecektir. Bunları görmezden gelin ve **Close**'a tıklayın.

TUI kapandığında, aşağıdaki konsol çıktısını görmelisiniz.

??? success "Komut çıktısı"

    ```console
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    INFO     Launching interactive nf-core pipeline creation tool.
    ```

Pipeline oluşturmanın çalıştığına dair konsol çıktısında açık bir onay yoktur, ancak `core-hello` adında yeni bir dizin görmelisiniz.

Şablonu kullanarak kendinize ne kadar iş kazandırdığınızı görmek için yeni dizinin içeriğini görüntüleyin.

```bash
tree core-hello
```

??? abstract "Dizin içeriği"

    ```console
    core-hello/
    ├── assets
    │   ├── samplesheet.csv
    │   └── schema_input.json
    ├── conf
    │   ├── base.config
    │   ├── modules.config
    │   ├── test.config
    │   └── test_full.config
    ├── docs
    │   ├── output.md
    │   ├── README.md
    │   └── usage.md
    ├── main.nf
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── README.md
    ├── subworkflows
    │   ├── local
    │   │   └── utils_nfcore_hello_pipeline
    │   │       └── main.nf
    │   └── nf-core
    │       ├── utils_nextflow_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       └── nextflow.config
    │       ├── utils_nfcore_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       ├── main.workflow.nf.test.snap
    │       │       └── nextflow.config
    │       └── utils_nfschema_plugin
    │           ├── main.nf
    │           ├── meta.yml
    │           └── tests
    │               ├── main.nf.test
    │               ├── nextflow.config
    │               └── nextflow_schema.json
    └── workflows
        └── hello.nf

    14 directories, 34 files
    ```

Bu çok fazla dosya!

Umarız bunların çoğunu `nf-core/demo` pipeline yapısını incelediğimizde karşılaştığımız dosyalarla aynı olarak tanırsınız.
Ancak hala biraz kaybolmuş hissediyorsanız endişelenmeyin; bu eğitim boyunca önemli kısımları birlikte inceleyeceğiz.

!!! note

    Bu eğitimin ilk bölümünde incelediğimiz `nf-core/demo` pipeline'ına kıyasla önemli bir fark, `modules` dizininin olmamasıdır.
    Bunun nedeni, varsayılan nf-core modüllerinden herhangi birini dahil etmeyi seçmemiş olmamızdır.

### 1.2. İskeletin işlevsel olduğunu test etme

İnanın ya da inanmayın, gerçek iş yapmak için henüz herhangi bir modül eklememiş olsanız bile, pipeline iskeleti aslında test profilini kullanarak çalıştırılabilir, tıpkı `nf-core/demo` pipeline'ını çalıştırdığımız gibi.

```bash
nextflow run ./core-hello -profile docker,test --outdir core-hello-results
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `./core-hello/main.nf` [scruffy_marconi] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-47-18

    Core Nextflow options
      runName                   : scruffy_marconi
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    -[core/hello] Pipeline completed successfully-
    ```

Bu, tüm temel bağlantıların yerinde olduğunu gösterir.
Peki çıktılar nerede? Var mı?

Aslında, standart yürütme raporlarını içeren `core-hello-results` adında yeni bir sonuç dizini oluşturuldu:

```bash
tree core-hello-results
```

??? abstract "Dizin içeriği"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-18.json
        └── pipeline_dag_2025-11-21_04-47-18.html

    1 directory, 6 files
    ```

Ne çalıştırıldığını görmek için raporlara bir göz atabilirsiniz ve cevap şu: hiçbir şey!

![boş yürütme zaman çizelgesi raporu](./img/execution_timeline_empty.png)

Kodda gerçekte ne olduğuna bir göz atalım.

### 1.3. Yer tutucu workflow'u inceleme

`main.nf` dosyasının içine bakarsanız, `workflows/hello` içinden `HELLO` adlı bir workflow import ettiğini göreceksiniz.

Bu, Bölüm 1'de karşılaştığımız `workflows/demo.nf` workflow'una eşdeğerdir ve bazı nf-core işlevleri zaten yerinde olan ilgilendiğimiz workflow için yer tutucu olarak hizmet eder.

```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="15 17 19 35"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

[Hello Nextflow](../hello_nextflow/index.md)'da geliştirilen temel bir Nextflow workflow'una kıyasla, burada yeni olan birkaç şey fark edeceksiniz (yukarıdaki vurgulanan satırlar):

- Workflow bloğunun bir adı var
- Workflow girdileri `take:` anahtar kelimesi kullanılarak bildirilir ve kanal oluşturma üst workflow'a taşınır
- Workflow içeriği bir `main:` bloğunun içine yerleştirilir
- Çıktılar `emit:` anahtar kelimesi kullanılarak bildirilir

Bunlar, workflow'u **birleştirilebilir** yapan Nextflow'un isteğe bağlı özellikleridir, yani başka bir workflow içinden çağrılabilir.

!!! note "Derinlemesine birleştirilebilir workflow'lar"

    [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Side Quest, workflow kompozisyonunu çok daha derinlemesine inceliyor, birden fazla workflow'u birlikte nasıl oluşturacağınızı ve aralarındaki karmaşık veri akışlarını nasıl yöneteceğinizi içeriyor. Burada birleştirilebilirliği tanıtıyoruz çünkü bu, pipeline başlatma, ana analiz workflow'u ve tamamlama görevlerini ayrı, yeniden kullanılabilir bileşenlere organize etmek için iç içe workflow'lar kullanan nf-core şablon mimarisinin temel bir gereksinimidir.

İlgilendiğimiz workflow'dan ilgili mantığı bu yapıya eklemek zorunda kalacağız.
Bunun için ilk adım orijinal workflow'umuzu birleştirilebilir hale getirmektir.

### Çıkarım

Artık nf-core araçlarını kullanarak bir pipeline iskeleti oluşturmayı biliyorsunuz.

### Sırada ne var?

Basit bir workflow'u nf-core uyumlu hale getirmenin bir ön adımı olarak nasıl birleştirilebilir yapacağınızı öğrenin.

---

## 2. Orijinal Hello Nextflow workflow'unu birleştirilebilir hale getirme

Şimdi workflow'umuzu nf-core iskeletine entegre etmek için çalışma zamanı. Hatırlatma olarak, [Hello Nextflow](../hello_nextflow/index.md) eğitim kursumuzdaki workflow ile çalışıyoruz.

!!! tip

    O pipeline'a aşina değilseniz veya hatırlatmaya ihtiyacınız varsa, [The Hello pipeline](../info/hello_pipeline.md) sayfasına bakın.

Size, tamamlanmış Hello Nextflow workflow'unun temiz, tamamen işlevsel bir kopyasını `original-hello` dizininde modülleri ve girdi olarak kullanmasını beklediği varsayılan CSV dosyası ile birlikte sağlıyoruz.

```bash
tree original-hello/
```

??? abstract "Dizin içeriği"

    ```console
    original-hello/
    ├── hello.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config

    1 directory, 6 files
    ```

Çalıştığından emin olmak için çekinmeden çalıştırın:

```bash
nextflow run original-hello/hello.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/hello.nf` [goofy_babbage] DSL2 - revision: e9e72441e9

    executor >  local (8)
    [a4/081cec] sayHello (1)       | 3 of 3 ✔
    [e7/7e9058] convertToUpper (3) | 3 of 3 ✔
    [0c/17263b] collectGreetings   | 1 of 1 ✔
    [94/542280] cowpy              | 1 of 1 ✔
    ```

Kodu incelemek için `hello.nf` workflow dosyasını açalım, aşağıda tam olarak gösterilmiştir (modüller içinde olan process'ler sayılmaz):

```groovy title="original-hello/hello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
* Pipeline parametreleri
*/
params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Modülleri dahil et
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow {

  // bir CSV dosyasından girdiler için bir kanal oluştur
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // bir selamlama yayınla
  sayHello(greeting_ch)

  // selamlamayı büyük harfe dönüştür
  convertToUpper(sayHello.out)

  // tüm selamlamaları tek bir dosyada topla
  collectGreetings(convertToUpper.out.collect(), params.batch)

  // cowpy ile selamlamaların ASCII sanatını oluştur
  cowpy(collectGreetings.out.outfile, params.character)
}
```

Gördüğünüz gibi, bu workflow kendi başına çalıştırılabilen basit bir adsız workflow olarak yazılmıştır.
Bunu nf-core şablonunun gerektirdiği gibi bir üst workflow içinden çalıştırılabilir hale getirmek için **birleştirilebilir** yapmamız gerekiyor.

Gerekli değişiklikleri tek tek inceleyelim.

### 2.1. Workflow'a ad verme

İlk olarak, workflow'a bir üst workflow'dan ona başvurabilmek için bir ad verelim.

=== "Sonra"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow HELLO {
    ```

=== "Önce"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow {
    ```

Workflow adları için modül adlarında olduğu gibi aynı kurallar geçerlidir.

### 2.2. Kanal oluşturmayı `take` ile değiştirme

Şimdi, kanal oluşturmayı beklenen girdileri bildiren basit bir `take` ifadesi ile değiştirin.

=== "Sonra"

    ```groovy title="original-hello/hello.nf" linenums="18"
        take:
        // channel of greetings
        greeting_ch
    ```

=== "Önce"

    ```groovy title="original-hello/hello.nf" linenums="18"
        // bir CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.greeting)
                            .splitCsv()
                            .map { line -> line[0] }
    ```

Bu, girdilerin nasıl sağlandığının ayrıntılarını üst workflow'a bırakır.

Bu arada, `params.greeting = 'greetings.csv'` satırını da yorum satırı yapabiliriz

=== "Sonra"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parametreleri
        */
        //params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

=== "Önce"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parametreleri
        */
        params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

!!! note

    Nextflow dil sunucusu uzantısı yüklüyse, sözdizimi denetleyicisi kodunuzu kırmızı dalgalı çizgilerle işaretleyecektir.
    Bunun nedeni, bir `take:` ifadesi koyarsanız, aynı zamanda bir `main:` de olması gerektiğidir.

    Bunu bir sonraki adımda ekleyeceğiz.

### 2.3. Workflow işlemlerinin önüne `main` ifadesi ekleme

Ardından, workflow gövdesinde çağrılan işlemlerin geri kalanından önce bir `main` ifadesi ekleyin.

=== "Sonra"

    ```groovy title="original-hello/hello.nf" linenums="22" hl_lines="1"
        main:

        // bir selamlama yayınla
        sayHello(greeting_ch)

        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)

        // tüm selamlamaları tek bir dosyada topla
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // cowpy ile ASCII art oluştur
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Önce"

    ```groovy title="original-hello/hello.nf" linenums="21"
        // bir selamlama yayınla
        sayHello(greeting_ch)

        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)

        // tüm selamlamaları tek bir dosyada topla
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // cowpy ile ASCII art oluştur
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Bu temel olarak 'bu workflow'un _yaptığı_ budur' der.

### 2.4. `emit` ifadesi ekleme

Son olarak, workflow'un son çıktılarının ne olduğunu bildiren bir `emit` ifadesi ekleyin.

```groovy title="original-hello/hello.nf" linenums="35"
    emit:
    cowpy_hellos = cowpy.out
```

Bu, orijinal workflow'a kıyasla koda tamamen yeni bir eklentidir.

### 2.5. Tamamlanan değişikliklerin özeti

Tüm değişiklikleri açıklandığı gibi yaptıysanız, workflow'unuz şimdi şöyle görünmelidir:

```groovy title="original-hello/hello.nf" linenums="1" hl_lines="16 18-20 22 36-37"
#!/usr/bin/env nextflow

/*
* Pipeline parametreleri
*/
// params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Modülleri dahil et
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow HELLO {

    take:
    // channel of greetings
    greeting_ch

    main:

    // bir selamlama yayınla
    sayHello(greeting_ch)

    // selamlamayı büyük harfe dönüştür
    convertToUpper(sayHello.out)

    // tüm selamlamaları tek bir dosyada topla
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // cowpy ile selamlamaların ASCII sanatını oluştur
    cowpy(collectGreetings.out.outfile, params.character)

    emit:
    cowpy_hellos = cowpy.out
}
```

Bu, Nextflow'un ihtiyaç duyduğu her şeyi açıklar, girdi kanalına ne besleyeceğimiz HARİÇ.
Bu, **giriş noktası** workflow'u olarak da adlandırılan üst workflow'da tanımlanacaktır.

### 2.6. Kukla bir giriş noktası workflow'u yapma

Birleştirilebilir workflow'umuzu karmaşık nf-core iskeletine entegre etmeden önce, doğru çalıştığını doğrulayalım.
Birleştirilebilir workflow'u izole bir şekilde test etmek için basit bir kukla giriş noktası workflow'u yapabiliriz.

Aynı `original-hello` dizininde `main.nf` adında boş bir dosya oluşturun.

```bash
touch original-hello/main.nf
```

Aşağıdaki kodu `main.nf` dosyasına kopyalayın.

```groovy title="original-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow

// import the workflow code from the hello.nf file
include { HELLO } from './hello.nf'

// declare input parameter
params.greeting = 'greetings.csv'

workflow {
  // bir CSV dosyasından girdiler için bir kanal oluştur
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // call the imported workflow on the channel of greetings
  HELLO(greeting_ch)

  // view the outputs emitted by the workflow
  HELLO.out.view { output -> "Output: $output" }
}
```

Burada yapılacak iki önemli gözlem var:

- İçe aktarılan workflow'u çağırma sözdizimi, modülleri çağırma sözdizimi ile esasen aynıdır.
- Girdileri workflow'a çekmeyle ilgili her şey (girdi parametresi ve kanal oluşturma) artık bu üst workflow'da bildirilir.

!!! note

    Giriş noktası workflow dosyasını `main.nf` olarak adlandırmak bir kural, bir gereklilik değildir.

    Bu kurala uyarsanız, `nextflow run` komutunuzda workflow dosya adını belirtmeyi atlayabilirsiniz.
    Nextflow, yürütme dizininde `main.nf` adlı bir dosyayı otomatik olarak arayacaktır.

    Ancak, isterseniz giriş noktası workflow dosyasını başka bir şey olarak adlandırabilirsiniz.
    Bu durumda, `nextflow run` komutunuzda workflow dosya adını belirttiğinizden emin olun.

### 2.7. Workflow'un çalıştığını test etme

Sonunda birleştirilebilir workflow'un çalıştığını doğrulamak için ihtiyacımız olan tüm parçalara sahibiz.
Hadi çalıştıralım!

```bash
nextflow run ./original-hello
```

Burada `main.nf` adlandırma kuralını kullanmanın avantajını görüyorsunuz.
Giriş noktası workflow'unu `something_else.nf` olarak adlandırsaydık, `nextflow run original-hello/something_else.nf` yapmak zorunda kalacaktık.

Tüm değişiklikleri doğru yaptıysanız, bu tamamlanana kadar çalışmalıdır.

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/main.nf` [friendly_wright] DSL2 - revision: 1ecd2d9c0a

    executor >  local (8)
    [24/c6c0d8] HELLO:sayHello (3)       | 3 of 3 ✔
    [dc/721042] HELLO:convertToUpper (3) | 3 of 3 ✔
    [48/5ab2df] HELLO:collectGreetings   | 1 of 1 ✔
    [e3/693b7e] HELLO:cowpy              | 1 of 1 ✔
    Output: /workspaces/training/hello-nf-core/work/e3/693b7e48dc119d0c54543e0634c2e7/cowpy-COLLECTED-test-batch-output.txt
    ```

Bu, HELLO workflow'umuzu başarıyla birleştirilebilir hale getirdiğimiz anlamına gelir.

### Çıkarım

Bir workflow'u adlandırarak ve `take`, `main` ve `emit` ifadeleri ekleyerek nasıl birleştirilebilir yapacağınızı ve bir giriş noktası workflow'undan nasıl çağıracağınızı biliyorsunuz.

### Sırada ne var?

Temel bir birleştirilebilir workflow'u nf-core iskeletine nasıl aşılayacağınızı öğrenin.

---

## 3. Güncellenmiş workflow mantığını yer tutucu workflow'a uydurma

Artık birleştirilebilir workflow'umuzun doğru çalıştığını doğruladığımıza göre, bölüm 1'de oluşturduğumuz nf-core pipeline iskeletine geri dönelim.
Az önce geliştirdiğimiz birleştirilebilir workflow'u nf-core şablon yapısına entegre etmek istiyoruz, böylece sonuç şuna benzer bir şey olmalıdır.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/core-hello.svg"
</figure>

Peki bunu nasıl gerçekleştiririz? `core-hello/workflows/hello.nf` içindeki (nf-core iskeleti) `HELLO` workflow'unun mevcut içeriğine bir göz atalım.

```groovy title="core-hello/workflows/hello.nf" linenums="1"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Genel olarak bu kod, workflow'da çalıştırılan herhangi bir yazılım aracının sürümünü yakalamakla ilgisi olan bazı temizlik dışında çok az şey yapar.

Bölüm 2'de geliştirdiğimiz orijinal workflow'un birleştirilebilir versiyonundan ilgili kodu eklememiz gerekiyor.

Bunu şu aşamalarda ele alacağız:

1. Modülleri kopyalama ve modül import'larını kurma
2. `take` bildirimini olduğu gibi bırakma
3. Workflow mantığını `main` bloğuna ekleme
4. `emit` bloğunu güncelleme

!!! note

    Bu ilk geçiş için versiyon yakalamayı görmezden geliyoruz ve bunun nasıl bağlanacağını bu eğitimin daha sonraki bir bölümünde inceleyeceğiz.

### 3.1. Modülleri kopyalama ve modül import'larını kurma

Hello Nextflow workflow'umuzdaki dört process, `original-hello/modules/` içinde modül olarak saklanır.
Bu modülleri nf-core proje yapısına (`core-hello/modules/local/` altına) kopyalamamız ve nf-core workflow dosyasına import ifadeleri eklememiz gerekiyor.

İlk olarak modül dosyalarını `original-hello/`'dan `core-hello/`'ya kopyalayalım:

```bash
mkdir -p core-hello/modules/local/
cp original-hello/modules/* core-hello/modules/local/.
```

Şimdi modül dizininin `core-hello/` altında listelendiğini görmelisiniz.

```bash
tree core-hello/modules
```

??? abstract "Dizin içeriği"

    ```console
    core-hello/modules
    └── local
        ├── collectGreetings.nf
        ├── convertToUpper.nf
        ├── cowpy.nf
        └── sayHello.nf

    1 directory, 4 files
    ```

Şimdi modül import ifadelerini kuralım.

Bunlar `original-hello/hello.nf` workflow'undaki import ifadeleriydi:

```groovy title="original-hello/hello.nf" linenums="9"
// Modülleri dahil et
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'
```

`core-hello/workflows/hello.nf` dosyasını açın ve bu import ifadelerini aşağıda gösterildiği gibi ona aktarın.

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="8-11"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    ```

Burada iki ilginç gözlem daha:

- Import ifadelerinin biçimlendirmesini nf-core stil kuralına uyacak şekilde uyarladık.
- Modüllerin göreceli yollarını, artık farklı bir iç içe yerleşme seviyesinde saklandıklarını yansıtacak şekilde güncelledik.

### 3.2. `take` bildirimini olduğu gibi bırakma

nf-core projesinin, tipik olarak sütunsal veri içeren bir CSV dosyası olan samplesheet kavramı etrafında çok sayıda önceden oluşturulmuş işlevselliği vardır.
`greetings.csv` dosyamız esasen bu olduğundan, mevcut `take` bildirimini olduğu gibi tutacağız ve bir sonraki adımda sadece girdi kanalının adını güncelleyeceğiz.

```groovy title="core-hello/workflows/hello.nf" linenums="21"
    take:
    ch_samplesheet // channel: samplesheet read in from --input
```

Girdi işleme bu workflow'un yukarısında yapılacaktır (bu kod dosyasında değil).

### 3.3. Workflow mantığını `main` bloğuna ekleme

Şimdi modüllerimiz workflow'a açık olduğuna göre, workflow mantığını `main` bloğuna ekleyebiliriz.

Hatırlatma olarak, orijinal workflow'daki ilgili kod şudur, birleştirilebilir yaptığımızda pek değişmedi (sadece `main:` satırını ekledik):

```groovy title="original-hello/hello.nf" linenums="22"
    main:

    // bir selamlama yayınla
    sayHello(greeting_ch)

    // selamlamayı büyük harfe dönüştür
    convertToUpper(sayHello.out)

    // tüm selamlamaları tek bir dosyada topla
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // cowpy ile selamlamaların ASCII sanatını oluştur
    cowpy(collectGreetings.out.outfile, params.character)
```

`main:`'den sonra gelen kodu workflow'un yeni versiyonuna kopyalamamız gerekiyor.

Orada zaten workflow tarafından çalıştırılan araçların sürümlerini yakalamakla ilgili bazı kodlar var. Şimdilik bunu olduğu gibi bırakacağız (araç versiyonlarıyla daha sonra ilgileneceğiz).
`ch_versions = channel.empty()` başlatmasını en üstte tutacağız, ardından workflow mantığımızı ekleyeceğiz, versiyon harmanlama kodunu sonda tutacağız.
Bu sıralama mantıklıdır çünkü gerçek bir pipeline'da, process'ler workflow çalışırken `ch_versions` kanalına eklenecek sürüm bilgisi yayar.

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19" hl_lines="10-20"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input

        main:

        ch_versions = Channel.empty()

        // bir selamlama yayınla
        sayHello(greeting_ch)

        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)

        // tüm selamlamaları tek bir dosyada topla
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // cowpy ile ASCII art oluştur
        cowpy(collectGreetings.out.outfile, params.character)

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input
        main:

        ch_versions = Channel.empty()

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

Ayrıca kodu daha okunabilir yapmak için `main:`'den önce boş bir satır eklediğimizi fark edeceksiniz.

Bu harika görünüyor, ancak `sayHello()` process'ine ilettiğimiz kanalın adını aşağıda gösterildiği gibi `greeting_ch`'den `ch_samplesheet`'e, `take:` anahtar kelimesi altında yazılanla eşleşecek şekilde güncellememiz gerekiyor.

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // bir selamlama yayınla (updated to use the nf-core convention for samplesheets)
        sayHello(ch_samplesheet)
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // bir selamlama yayınla
        sayHello(greeting_ch)
    ```

Şimdi workflow mantığı doğru şekilde bağlandı.

### 3.4. `emit` bloğunu güncelleme

Son olarak, workflow'un son çıktılarının bildirimini içerecek şekilde `emit` bloğunu güncellememiz gerekiyor.

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55"
        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

Bu, HELLO workflow'unun kendisinde yapmamız gereken değişiklikleri tamamlar.
Bu noktada, uygulamaya koymayı amaçladığımız genel kod yapısına ulaştık.

### Çıkarım

Birleştirilebilir bir workflow'un temel parçalarını bir nf-core yer tutucu workflow'una nasıl yerleştireceğinizi biliyorsunuz.

### Sırada ne var?

nf-core pipeline iskeletinde girdilerin nasıl işlendiğini nasıl uyarlayacağınızı öğrenin.

---

## 4. Girdi işlemeyi uyarlama

Artık workflow mantığımızı nf-core iskeletine başarıyla entegre ettiğimize göre, bir kritik parçayı daha ele almamız gerekiyor: girdi verilerimizin doğru şekilde işlendiğinden emin olmak.
nf-core şablonu, karmaşık genomik veri kümeleri için tasarlanmış sofistike girdi işleme ile birlikte gelir, bu yüzden onu daha basit `greetings.csv` dosyamızla çalışacak şekilde uyarlamamız gerekiyor.

### 4.1. Girdilerin nerede işlendiğini belirleme

İlk adım, girdi işlemenin nerede yapıldığını bulmaktır.

Hello Nextflow workflow'unu birleştirilebilir olacak şekilde yeniden yazdığımızda, girdi parametresi bildirimini bir seviye yukarı, `main.nf` giriş noktası workflow'una taşıdığımızı hatırlayabilirsiniz.
O halde, pipeline iskeleti kapsamında oluşturulan üst düzey `main.nf` giriş noktası workflow'una bir göz atalım:

```groovy title="core-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core/hello
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/core/hello
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HELLO  } from './workflows/hello'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hello_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hello_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CORE_HELLO {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    HELLO (
        samplesheet
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    CORE_HELLO (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

nf-core projesi iç içe subworkflow'ları yoğun bir şekilde kullanır, bu yüzden bu kısım ilk yaklaşımda biraz kafa karıştırıcı olabilir.

Burada önemli olan iki workflow tanımlanmıştır:

- `CORE_HELLO`, `core-hello/workflows/hello.nf` içinde az önce uyarlamayı bitirdiğimiz HELLO workflow'unu çalıştırmak için ince bir sarmalayıcıdır.
- `CORE_HELLO`'yu ve iki diğer subworkflow'u, `PIPELINE_INITIALISATION` ve `PIPELINE_COMPLETION`'ı çağıran adsız bir workflow.

İşte birbirleriyle nasıl ilişkili olduklarının bir diyagramı:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/hello-nested-workflows.svg"
</figure>

Önemlisi, bu seviyede bir girdi kanalı oluşturan herhangi bir kod bulamıyoruz, yalnızca `--input` parametresi aracılığıyla sağlanan bir samplesheet'e referanslar.

Biraz araştırma, girdi işlemenin uygun bir şekilde `core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf`'den import edilen `PIPELINE_INITIALISATION` subworkflow'u tarafından yapıldığını ortaya çıkarır.

Bu dosyayı açıp aşağı kaydırırsak, bu kod bloğuna geliriz:

```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76"
    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
```

Bu, samplesheet'i ayrıştıran ve HELLO workflow'u tarafından tüketilmeye hazır bir biçimde ileten kanal fabrikasıdır.

!!! note

    Yukarıdaki sözdizimi daha önce kullandığımızdan biraz farklı, ancak temelde bu:

    ```groovy
    channel.<...>.set { ch_samplesheet }
    ```

    şuna eşdeğerdir:

    ```groovy
    ch_samplesheet = channel.<...>
    ```

Bu kod, yazım zamanında çok alana özgü olan ve basit pipeline projemiz için uygun olmayan nf-core pipeline şablonuna dahil edilen örnek samplesheet'e oldukça özgü bazı ayrıştırma ve doğrulama adımları içerir.

### 4.2. Şablonlu girdi kanalı kodunu değiştirme

İyi haber, pipeline'ımızın ihtiyaçlarının çok daha basit olması, bu yüzden tüm bunları orijinal Hello Nextflow workflow'unda geliştirdiğimiz kanal oluşturma koduyla değiştirebiliriz.

Hatırlatma olarak, kanal oluşturma şöyle görünüyordu (çözümler dizininde görüldüğü gibi):

```groovy title="solutions/composable-hello/main.nf" linenums="10" hl_lines="4"
    // bir CSV dosyasından girdiler için bir kanal oluştur
    greeting_ch = channel.fromPath(params.greeting)
        .splitCsv()
        .map { line -> line[0] }
```

Yani sadece bunu küçük değişikliklerle başlatma workflow'una eklememiz gerekiyor: kanal adını `greeting_ch`'den `ch_samplesheet`'e ve parametre adını `params.greeting`'den `params.input`'a güncelliyoruz (vurgulanan satıra bakın).

=== "Sonra"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-7"
        //
        // Create channel from input file provided through params.input
        //

        ch_samplesheet = channel.fromPath(params.input)
            .splitCsv()
            .map { line -> line[0] }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

=== "Önce"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-23"
        //
        // Create channel from input file provided through params.input
        //

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

Bu, girdi işlemenin çalışmasını sağlamak için yapmamız gereken değişiklikleri tamamlar.

Mevcut haliyle, bu bizi şema doğrulaması için nf-core'un yerleşik yeteneklerinden yararlanmamıza izin vermez, ancak bunu daha sonra ekleyebiliriz.
Şimdilik, test verileri üzerinde başarıyla çalıştırabileceğimiz bir şeye ulaşmak için mümkün olduğunca basit tutmaya odaklanıyoruz.

### 4.3. Test profilini güncelleme

Test verileri ve parametrelerden bahsetmişken, bu pipeline'ın test profilini şablonda sağlanan örnek samplesheet yerine `greetings.csv` mini-samplesheet'i kullanacak şekilde güncelleyelim.

`core-hello/conf` altında, küçük bir veri örneğini ve tam boyutlu birini test etmek için tasarlanan iki şablonlu test profili buluyoruz: `test.config` ve `test_full.config`.
Pipeline'ımızın amacı göz önüne alındığında, tam boyutlu bir test profili kurmanın gerçekten bir anlamı yok, bu yüzden `test_full.config`'i görmezden gelmekten veya silmekten çekinmeyin.
Birkaç varsayılan parametreyle `greetings.csv` dosyamızda çalışacak şekilde `test.config`'i kurmaya odaklanacağız.

#### 4.3.1. `greetings.csv` dosyasını kopyalama

Öncelikle `greetings.csv` dosyasını pipeline projemizde uygun bir yere kopyalamamız gerekiyor.
Tipik olarak küçük test dosyaları `assets` dizininde saklanır, o halde dosyayı çalışma dizinimizden kopyalayalım.

```bash
cp greetings.csv core-hello/assets/.
```

Şimdi `greetings.csv` dosyası test girdisi olarak kullanılmaya hazır.

#### 4.3.2. `test.config` dosyasını güncelleme

Şimdi `test.config` dosyasını aşağıdaki gibi güncelleyebiliriz:

=== "Sonra"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-10"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Girdi verileri
        input  = "${projectDir}/assets/greetings.csv"

        // Diğer parametreler
        batch     = 'test'
        character = 'tux'
    }
    ```

=== "Önce"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-8"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Girdi verileri
        // TODO nf-core: Test verilerinizin nf-core/test-datasets üzerindeki yollarını belirtin
        // TODO nf-core: Test için gerekli parametreleri verin, böylece komut satırı bayrakları gerekmez
        input  = params.pipelines_testdata_base_path + 'viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
    }
    ```

Önemli noktalar:

- **`${projectDir}` kullanımı**: Bu, ana workflow betiğinin bulunduğu dizine (pipeline kökü) işaret eden bir Nextflow örtük değişkenidir. Bunu kullanmak, yolun pipeline'ın nereden çalıştırıldığından bağımsız olarak çalışmasını sağlar.
- **Mutlak yollar**: `${projectDir}` kullanarak, küçük test dosyaları için önemli olan mutlak bir yol oluşturuyoruz.
- **Test verisi konumu**: nf-core pipeline'ları tipik olarak küçük test dosyaları için pipeline repository'si içindeki `assets/` dizininde test verilerini saklar veya daha büyük dosyalar için harici test veri kümelerine referans verir.

Ve bu arada, bunun çok basit makinelerde (Github Codespaces'teki minimal VM'ler gibi) çalışacağından emin olmak için varsayılan kaynak limitlerini sıkılaştıralım:

=== "Sonra"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 2,
            memory: '4.GB',
            time: '1.h'
        ]
    }
    ```

=== "Önce"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 4,
            memory: '15.GB',
            time: '1.h'
        ]
    }
    ```

Bu, yapmamız gereken kod değişikliklerini tamamlar.

### 4.4. Pipeline'ı test profiliyle çalıştırma

Bu çok şeydi, ancak sonunda pipeline'ı çalıştırmayı deneyebiliriz!
Henüz doğrulamayı kurmadığımız için komut satırına `--validate_params false` eklememiz gerektiğini unutmayın (bu daha sonra gelecek).

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

Tüm değişiklikleri doğru yaptıysanız, tamamlanana kadar çalışmalıdır.

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `core-hello/main.nf` [condescending_allen] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-11-21_07-29-37

    Core Nextflow options
      runName                   : condescending_allen
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (1)
    [ed/727b7e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [45/bb6096] CORE_
