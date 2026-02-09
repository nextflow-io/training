# Bölüm 2: Hello'yu nf-core için Yeniden Yazın

Bu Hello nf-core eğitim kursunun ikinci bölümünde, [Hello Nextflow](../hello_nextflow/index.md) başlangıç kursunda üretilen pipeline'ın nf-core uyumlu bir versiyonunu nasıl oluşturacağınızı gösteriyoruz.

Eğitimin ilk bölümünde nf-core pipeline'larının oldukça detaylı bir yapıya ve birçok yardımcı dosyaya sahip olduğunu fark etmişsinizdir.
Tüm bunları sıfırdan oluşturmak çok yorucu olacağından, nf-core topluluğu süreci başlatmak için bunun yerine bir şablondan yapılmasını sağlayan araçlar geliştirmiştir.

Size bu araçları kullanarak bir pipeline iskeleti oluşturmayı, ardından mevcut 'normal' pipeline kodunu nf-core iskeleti üzerine nasıl uyarlayacağınızı göstereceğiz.

Hello pipeline'ına aşina değilseniz veya hatırlatmaya ihtiyacınız varsa, [bu bilgi sayfasına](../info/hello_pipeline.md) bakın.

---

## 1. Yeni bir pipeline projesi oluşturun

İlk olarak, yeni pipeline için iskelet oluşturuyoruz.

!!! note "Not"

    Terminalinizde `hello-nf-core` dizininde olduğunuzdan emin olun.

### 1.1. Şablon tabanlı pipeline oluşturma aracını çalıştırın

`nf-core pipelines create` komutuyla yeni bir pipeline oluşturarak başlayalım.
Bu, nf-core temel şablonunu kullanarak, pipeline adı, açıklaması ve yazarıyla özelleştirilmiş yeni bir pipeline iskeleti oluşturacaktır.

```bash
nf-core pipelines create
```

Bu komutu çalıştırmak, pipeline oluşturma için bir Metin Kullanıcı Arayüzü (TUI) açacaktır:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/VwjXNXONHlY?si=d0HkFSISnKn76TeI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

Bu TUI, pipeline'ınız hakkında temel bilgiler sağlamanızı isteyecek ve pipeline iskeleti içine dahil edilecek veya hariç tutulacak özellikleri seçme imkanı sunacaktır.

- Hoş geldiniz ekranında **Let's go!** düğmesine tıklayın.
- `Choose pipeline type` ekranında **Custom** seçeneğine tıklayın.
- Pipeline bilgilerinizi aşağıdaki gibi girin (`< YOUR NAME >` yerine kendi adınızı yazın), ardından **Next** düğmesine tıklayın.

```
[ ] GitHub organisation: core
[ ] Workflow name: hello
[ ] A short description of your pipeline: A basic nf-core style version of Hello Nextflow
[ ] Name of the main author(s): < YOUR NAME >
```

- Template features ekranında, `Toggle all features` seçeneğini **kapalı** konuma getirin, ardından aşağıdakileri seçerek **etkinleştirin**. Seçimlerinizi kontrol edin ve **Continue** düğmesine tıklayın.

```
[ ] Add testing profiles
[ ] Use nf-core components
[ ] Use nf-schema
[ ] Add configuration files
[ ] Add documentation
```

- `Final details` ekranında **Finish** düğmesine tıklayın. Pipeline'ın oluşturulmasını bekleyin, ardından **Continue** düğmesine tıklayın.
- Create GitHub repository ekranında **Finish without creating a repo** seçeneğine tıklayın. Bu, daha sonra bir GitHub deposu oluşturmak için talimatları gösterecektir. Bunları görmezden gelin ve **Close** düğmesine tıklayın.

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

Konsol çıktısında pipeline oluşturmanın başarılı olduğuna dair açık bir onay yoktur, ancak `core-hello` adında yeni bir dizin görmelisiniz.

Şablonu kullanarak kendinize ne kadar iş tasarrufu sağladığınızı görmek için yeni dizinin içeriğini görüntüleyin.

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

Umarım bunların çoğunu `nf-core/demo` pipeline yapısını incelediğimizde karşılaştığımız dosyalarla aynı olarak tanırsınız.
Ancak hala biraz kaybolmuş hissediyorsanız endişelenmeyin; bu eğitim boyunca önemli kısımları birlikte inceleyeceğiz.

!!! note "Not"

    Bu eğitimin ilk bölümünde incelediğimiz `nf-core/demo` pipeline'ına kıyasla önemli bir fark, `modules` dizininin olmamasıdır.
    Bunun nedeni, varsayılan nf-core modüllerinden herhangi birini dahil etmeyi seçmemiş olmamızdır.

### 1.2. İskeletin işlevsel olduğunu test edin

İster inanın ister inanmayın, gerçek iş yapması için henüz herhangi bir modül eklememiş olsanız bile, pipeline iskeleti test profili kullanılarak çalıştırılabilir, tıpkı `nf-core/demo` pipeline'ını çalıştırdığımız gibi.

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

Neyin çalıştırıldığını görmek için raporlara göz atabilirsiniz ve cevap: hiçbir şey!

![boş yürütme zaman çizelgesi raporu](./img/execution_timeline_empty.png)

Kodda gerçekte ne olduğuna bir göz atalım.

### 1.3. Yer tutucu iş akışını inceleyin

`main.nf` dosyasının içine bakarsanız, `workflows/hello` dosyasından `HELLO` adlı bir iş akışı içe aktardığını göreceksiniz.

Bu, Bölüm 1'de karşılaştığımız `workflows/demo.nf` iş akışına eşdeğerdir ve bazı nf-core işlevleri zaten yerinde olan ilgilendiğimiz iş akışı için bir yer tutucu görevi görür.

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

[Hello Nextflow](../hello_nextflow/index.md)'da geliştirilen gibi temel bir Nextflow iş akışıyla karşılaştırıldığında, burada yeni olan birkaç şey fark edeceksiniz (yukarıda vurgulanan satırlar):

- İş akışı bloğunun bir adı var
- İş akışı girdileri `take:` anahtar kelimesi kullanılarak bildirilir ve kanal oluşturma üst iş akışına taşınır
- İş akışı içeriği bir `main:` bloğunun içine yerleştirilir
- Çıktılar `emit:` anahtar kelimesi kullanılarak bildirilir

Bunlar, iş akışını **birleştirilebilir** yapan Nextflow'un isteğe bağlı özellikleridir, yani başka bir iş akışı içinden çağrılabilir.

!!! note "Birleştirilebilir iş akışları derinlemesine"

    [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Yan Görevi, iş akışı kompozisyonunu çok daha derinlemesine incelemektedir; birden fazla iş akışını birlikte nasıl oluşturacağınızı ve aralarında karmaşık veri akışlarını nasıl yöneteceğinizi içerir. Birleştirilebilirliği burada tanıtıyoruz çünkü bu, nf-core şablon mimarisinin temel bir gereksinimidir; bu mimari, pipeline başlatma, ana analiz iş akışı ve tamamlama görevlerini ayrı, yeniden kullanılabilir bileşenler halinde düzenlemek için iç içe iş akışları kullanır.

İlgili mantığı ilgilendiğimiz iş akışından bu yapıya takmamız gerekecek.
Bunun için ilk adım, orijinal iş akışımızı birleştirilebilir hale getirmektir.

### Özet

Artık nf-core araçlarını kullanarak bir pipeline iskeleti oluşturmayı biliyorsunuz.

### Sırada ne var?

Basit bir iş akışını nf-core uyumlu hale getirmenin bir ön hazırlığı olarak nasıl birleştirilebilir yapacağınızı öğrenin.

---

## 2. Orijinal Hello Nextflow iş akışını birleştirilebilir yapın

Şimdi iş akışımızı nf-core iskeletine entegre etmeye başlama zamanı.
Hatırlatma olarak, [Hello Nextflow](../hello_nextflow/index.md) eğitim kursumuzdaki iş akışıyla çalışıyoruz.

!!! tip "İpucu"

    Bu pipeline'a aşina değilseniz veya hatırlatmaya ihtiyacınız varsa, [The Hello pipeline](../info/hello_pipeline.md) sayfasına bakın.

Size, tamamlanmış Hello Nextflow iş akışının temiz, tamamen işlevsel bir kopyasını `original-hello` dizininde modülleri ve girdi olarak kullanmayı beklediği varsayılan CSV dosyasıyla birlikte sağlıyoruz.

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

Çalıştığından emin olmak için çalıştırmaktan çekinmeyin:

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

Kodu incelemek için `hello.nf` iş akışı dosyasını açalım, aşağıda tam olarak gösterilmiştir (modüllerde olan süreçler sayılmaz):

```groovy title="original-hello/hello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow {

  // create a channel for inputs from a CSV file
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // emit a greeting
  sayHello(greeting_ch)

  // convert the greeting to uppercase
  convertToUpper(sayHello.out)

  // collect all the greetings into one file
  collectGreetings(convertToUpper.out.collect(), params.batch)

  // generate ASCII art of the greetings with cowpy
  cowpy(collectGreetings.out.outfile, params.character)
}
```

Gördüğünüz gibi, bu iş akışı kendi başına çalıştırılabilen basit bir isimsiz iş akışı olarak yazılmıştır.
nf-core şablonunun gerektirdiği gibi bir üst iş akışı içinden çalıştırılabilir hale getirmek için **birleştirilebilir** yapmamız gerekiyor.

Gerekli değişiklikleri tek tek inceleyelim.

### 2.1. İş akışını adlandırın

İlk olarak, iş akışına bir üst iş akışından referans verebilmemiz için bir ad verelim.

=== "Sonra"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow HELLO {
    ```

=== "Önce"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow {
    ```

İş akışı adları için modül adlarıyla aynı kurallar geçerlidir.

### 2.2. Kanal oluşturmayı `take` ile değiştirin

Şimdi, kanal oluşturmayı beklenen girdileri bildiren basit bir `take` ifadesiyle değiştirin.

=== "Sonra"

    ```groovy title="original-hello/hello.nf" linenums="18"
        take:
        // channel of greetings
        greeting_ch
    ```

=== "Önce"

    ```groovy title="original-hello/hello.nf" linenums="18"
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.greeting)
                            .splitCsv()
                            .map { line -> line[0] }
    ```

Bu, girdilerin nasıl sağlandığının ayrıntılarını üst iş akışına bırakır.

Bu arada, `params.greeting = 'greetings.csv'` satırını da yorumlayabiliriz

=== "Sonra"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        //params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

=== "Önce"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

!!! note "Not"

    Nextflow dil sunucusu uzantısı yüklüyse, sözdizimi denetleyicisi kodunuzu kırmızı dalgalı çizgilerle işaretleyecektir.
    Bunun nedeni, bir `take:` ifadesi koyarsanız, aynı zamanda bir `main:` ifadesine de sahip olmanız gerektiğidir.

    Bunu bir sonraki adımda ekleyeceğiz.

### 2.3. İş akışı işlemlerinin önüne `main` ifadesi ekleyin

Ardından, iş akışının gövdesinde çağrılan işlemlerin geri kalanından önce bir `main` ifadesi ekleyin.

=== "Sonra"

    ```groovy title="original-hello/hello.nf" linenums="22" hl_lines="1"
        main:

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Önce"

    ```groovy title="original-hello/hello.nf" linenums="21"
        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Bu temelde 'bu iş akışının _yaptığı şey_ budur' der.

### 2.4. `emit` ifadesi ekleyin

Son olarak, iş akışının nihai çıktılarının ne olduğunu bildiren bir `emit` ifadesi ekleyin.

```groovy title="original-hello/hello.nf" linenums="35"
    emit:
    cowpy_hellos = cowpy.out
```

Bu, orijinal iş akışına kıyasla koda tamamen yeni bir eklentidir.

### 2.5. Tamamlanan değişikliklerin özeti

Tüm değişiklikleri açıklandığı gibi yaptıysanız, iş akışınız şimdi şöyle görünmelidir:

```groovy title="original-hello/hello.nf" linenums="1" hl_lines="16 18-20 22 36-37"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
// params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow HELLO {

    take:
    // channel of greetings
    greeting_ch

    main:

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // generate ASCII art of the greetings with cowpy
    cowpy(collectGreetings.out.outfile, params.character)

    emit:
    cowpy_hellos = cowpy.out
}
```

Bu, Nextflow'un ihtiyaç duyduğu her şeyi tanımlar, ANCAK girdi kanalına ne besleyeceğimizi tanımlamaz.
Bu, **giriş noktası** iş akışı olarak da adlandırılan üst iş akışında tanımlanacaktır.

### 2.6. Bir kukla giriş noktası iş akışı yapın

Birleştirilebilir iş akışımızı karmaşık nf-core iskeletine entegre etmeden önce, doğru çalıştığını doğrulayalım.
Birleştirilebilir iş akışını izole bir şekilde test etmek için basit bir kukla giriş noktası iş akışı yapabiliriz.

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
  // create a channel for inputs from a CSV file
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // call the imported workflow on the channel of greetings
  HELLO(greeting_ch)

  // view the outputs emitted by the workflow
  HELLO.out.view { output -> "Output: $output" }
}
```

Burada yapılması gereken iki önemli gözlem var:

- İçe aktarılan iş akışını çağırma sözdizimi, modülleri çağırma sözdizimi ile temelde aynıdır.
- Girdileri iş akışına çekmeyle ilgili her şey (girdi parametresi ve kanal oluşturma) artık bu üst iş akışında bildirilir.

!!! note "Not"

    Giriş noktası iş akışı dosyasını `main.nf` olarak adlandırmak bir kural, gereklilik değildir.

    Bu kuralı takip ederseniz, `nextflow run` komutunuzda iş akışı dosya adını belirtmeyi atlayabilirsiniz.
    Nextflow, yürütme dizininde otomatik olarak `main.nf` adlı bir dosya arayacaktır.

    Ancak, isterseniz giriş noktası iş akışı dosyasını başka bir şey olarak adlandırabilirsiniz.
    Bu durumda, `nextflow run` komutunuzda iş akışı dosya adını belirttiğinizden emin olun.

### 2.7. İş akışının çalıştığını test edin

Sonunda birleştirilebilir iş akışının çalıştığını doğrulamak için ihtiyacımız olan tüm parçalara sahibiz.
Hadi çalıştıralım!

```bash
nextflow run ./original-hello
```

Burada `main.nf` adlandırma kuralını kullanmanın avantajını görüyorsunuz.
Giriş noktası iş akışını `something_else.nf` olarak adlandırmış olsaydık, `nextflow run original-hello/something_else.nf` yapmamız gerekirdi.

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

Bu, HELLO iş akışımızı başarıyla birleştirilebilir hale getirdiğimiz anlamına gelir.

### Özet

Bir iş akışına ad vererek ve `take`, `main` ve `emit` ifadeleri ekleyerek nasıl birleştirilebilir yapacağınızı ve bir giriş noktası iş akışından nasıl çağıracağınızı biliyorsunuz.

### Sırada ne var?

Temel bir birleştirilebilir iş akışını nf-core iskeletine nasıl aşılayacağınızı öğrenin.

---

## 3. Güncellenmiş iş akışı mantığını yer tutucu iş akışına yerleştirin

Artık birleştirilebilir iş akışımızın doğru çalıştığını doğruladığımıza göre, bölüm 1'de oluşturduğumuz nf-core pipeline iskeletine geri dönelim.
Az önce geliştirdiğimiz birleştirilebilir iş akışını nf-core şablon yapısına entegre etmek istiyoruz, böylece sonuç şuna benzer görünmelidir.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/core-hello.svg"
</figure>

Peki bunu nasıl gerçekleştiririz? `core-hello/workflows/hello.nf` dosyasındaki (nf-core iskeleti) `HELLO` iş akışının mevcut içeriğine bir göz atalım.

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

Genel olarak bu kod, pipeline'da çalıştırılan herhangi bir yazılım aracının sürümünü yakalamakla ilgili bazı işlemler dışında çok az şey yapar.

Bölüm 2'de geliştirdiğimiz orijinal iş akışının birleştirilebilir versiyonundan ilgili kodu eklememiz gerekiyor.

Bunu aşağıdaki aşamalarda ele alacağız:

1. Modülleri kopyalayın ve modül içe aktarmalarını ayarlayın
2. `take` bildirimini olduğu gibi bırakın
3. İş akışı mantığını `main` bloğuna ekleyin
4. `emit` bloğunu güncelleyin

!!! note "Not"

    Bu ilk geçiş için sürüm yakalamayı görmezden geleceğiz ve bunu bu eğitimin sonraki bir bölümünde nasıl bağlayacağımıza bakacağız.

### 3.1. Modülleri kopyalayın ve modül içe aktarmalarını ayarlayın

Hello Nextflow iş akışımızdaki dört süreç, `original-hello/modules/` altında modül olarak saklanır.
Bu modülleri nf-core proje yapısına (`core-hello/modules/local/` altına) kopyalamamız ve nf-core iş akışı dosyasına içe aktarma ifadeleri eklememiz gerekiyor.

İlk olarak modül dosyalarını `original-hello/`'dan `core-hello/`'ya kopyalayalım:

```bash
mkdir -p core-hello/modules/local/
cp original-hello/modules/* core-hello/modules/local/.
```

Artık `core-hello/` altında listelenen modüller dizinini görmelisiniz.

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

Şimdi modül içe aktarma ifadelerini ayarlayalım.

Bunlar `original-hello/hello.nf` iş akışındaki içe aktarma ifadeleriydi:

```groovy title="original-hello/hello.nf" linenums="9"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'
```

`core-hello/workflows/hello.nf` dosyasını açın ve bu içe aktarma ifadelerini aşağıda gösterildiği gibi aktarın.

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

Burada iki ilginç gözlem daha var:

- İçe aktarma ifadelerinin biçimlendirmesini nf-core stil kuralını takip edecek şekilde uyarladık.
- Modüllere giden göreli yolları, artık farklı bir iç içe geçme seviyesinde saklandıklarını yansıtacak şekilde güncelledik.

### 3.2. `take` bildirimini olduğu gibi bırakın

nf-core projesinin, tipik olarak sütunlu veri içeren bir CSV dosyası olan samplesheet kavramı etrafında önceden oluşturulmuş birçok işlevselliği vardır.
`greetings.csv` dosyamız esasen bu olduğundan, mevcut `take` bildirimini olduğu gibi tutacağız ve bir sonraki adımda sadece girdi kanalının adını güncelleyeceğiz.

```groovy title="core-hello/workflows/hello.nf" linenums="21"
    take:
    ch_samplesheet // channel: samplesheet read in from --input
```

Girdi işleme bu iş akışının yukarısında yapılacaktır (bu kod dosyasında değil).

### 3.3. İş akışı mantığını `main` bloğuna ekleyin

Artık modüllerimiz iş akışı için kullanılabilir olduğuna göre, iş akışı mantığını `main` bloğuna takabiliriz.

Hatırlatma olarak, bu orijinal iş akışındaki ilgili koddur, birleştirilebilir yaptığımızda pek değişmedi (sadece `main:` satırını ekledik):

```groovy title="original-hello/hello.nf" linenums="22"
    main:

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // generate ASCII art of the greetings with cowpy
    cowpy(collectGreetings.out.outfile, params.character)
```

`main:` ifadesinden sonra gelen kodu iş akışının yeni versiyonuna kopyalamamız gerekiyor.

Orada zaten iş akışı tarafından çalıştırılan araçların sürümlerini yakalamakla ilgili bazı kodlar var. Şimdilik bunu olduğu gibi bırakacağız (araç sürümleriyle daha sonra ilgileneceğiz).
En üstte `ch_versions = channel.empty()` başlatmasını tutacağız, ardından iş akışı mantığımızı ekleyeceğiz, sürüm toplama kodunu sonda tutacağız.
Bu sıralama mantıklıdır çünkü gerçek bir pipeline'da, süreçler iş akışı çalışırken `ch_versions` kanalına eklenecek sürüm bilgilerini yayınlardı.

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19" hl_lines="10-20"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input

        main:

        ch_versions = Channel.empty()

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generate ASCII art of the greetings with cowpy
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

Ayrıca kodu daha okunabilir hale getirmek için `main:` ifadesinden önce boş bir satır eklediğimizi fark edeceksiniz.

Bu harika görünüyor, ancak `sayHello()` sürecine geçirdiğimiz kanalın adını `greeting_ch`'den `ch_samplesheet`'e aşağıda gösterildiği gibi güncellememiz gerekiyor, `take:` anahtar kelimesi altında yazılanla eşleşmesi için.

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emit a greeting (updated to use the nf-core convention for samplesheets)
        sayHello(ch_samplesheet)
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emit a greeting
        sayHello(greeting_ch)
    ```

Şimdi iş akışı mantığı doğru şekilde bağlanmıştır.

### 3.4. `emit` bloğunu güncelleyin

Son olarak, `emit` bloğunu iş akışının nihai çıktılarının bildirimini içerecek şekilde güncellememiz gerekiyor.

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

Bu, HELLO iş akışının kendisine yapmamız gereken değişiklikleri tamamlar.
Bu noktada, uygulamayı hedeflediğimiz genel kod yapısını elde ettik.

### Özet

Birleştirilebilir bir iş akışının temel parçalarını bir nf-core yer tutucu iş akışına nasıl yerleştireceğinizi biliyorsunuz.

### Sırada ne var?

nf-core pipeline iskeletinde girdilerin nasıl işlendiğini nasıl uyarlayacağınızı öğrenin.

---

## 4. Girdi işlemeyi uyarlayın

Artık iş akışı mantığımızı nf-core iskeletine başarıyla entegre ettiğimize göre, bir kritik parçayı daha ele almamız gerekiyor: girdi verilerimizin doğru şekilde işlendiğinden emin olmak.
nf-core şablonu, karmaşık genomik veri kümeleri için tasarlanmış sofistike girdi işlemeyle birlikte gelir, bu nedenle basit `greetings.csv` dosyamızla çalışacak şekilde uyarlamamız gerekiyor.

### 4.1. Girdilerin nerede işlendiğini belirleyin

İlk adım, girdi işlemenin nerede yapıldığını bulmaktır.

Hello Nextflow iş akışını birleştirilebilir olacak şekilde yeniden yazdığımızda, girdi parametresi bildirimini bir seviye yukarı, `main.nf` giriş noktası iş akışına taşıdığımızı hatırlayabilirsiniz.
O halde pipeline iskeleti kapsamında oluşturulan üst seviye `main.nf` giriş noktası iş akışına bir göz atalım:

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

nf-core projesi iç içe geçmiş alt iş akışlarını yoğun bir şekilde kullanır, bu nedenle bu kısım ilk yaklaşımda biraz kafa karıştırıcı olabilir.

Burada önemli olan iki iş akışının tanımlanmış olmasıdır:

- `CORE_HELLO`, `core-hello/workflows/hello.nf` dosyasında az önce uyarlamayı bitirdiğimiz HELLO iş akışını çalıştırmak için ince bir sarmalayıcıdır.
- `CORE_HELLO`'yu ve ayrıca `PIPELINE_INITIALISATION` ve `PIPELINE_COMPLETION` olmak üzere iki alt iş akışını çağıran isimsiz bir iş akışı.

İşte birbirleriyle nasıl ilişkili olduklarının bir diyagramı:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/hello-nested-workflows.svg"
</figure>

Önemlisi, bu seviyede bir girdi kanalı oluşturan herhangi bir kod bulamıyoruz, sadece `--input` parametresi aracılığıyla sağlanan bir samplesheet'e referanslar var.

Biraz araştırma, girdi işlemenin uygun bir şekilde `PIPELINE_INITIALISATION` alt iş akışı tarafından yapıldığını ortaya çıkarır; bu, `core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf` dosyasından içe aktarılır.

Bu dosyayı açıp aşağı kaydırırsak, şu kod parçasına geliyoruz:

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

Bu, samplesheet'i ayrıştıran ve HELLO iş akışı tarafından tüketilmeye hazır bir biçimde ileten kanal fabrikasıdır.

!!! note "Not"

    Yukarıdaki sözdizimi daha önce kullandığımızdan biraz farklı, ancak temelde bu:

    ```groovy
    channel.<...>.set { ch_samplesheet }
    ```

    buna eşdeğerdir:

    ```groovy
    ch_samplesheet = channel.<...>
    ```

Bu kod, yazım sırasında çok alana özgü olan ve basit pipeline projemiz için uygun olmayan nf-core pipeline şablonuna dahil edilen örnek samplesheet'e oldukça özgü bazı ayrıştırma ve doğrulama adımlarını içerir.

### 4.2. Şablonlu girdi kanalı kodunu değiştirin

İyi haber şu ki pipeline'ımızın ihtiyaçları çok daha basit, bu nedenle tüm bunları orijinal Hello Nextflow iş akışında geliştirdiğimiz kanal oluşturma koduyla değiştirebiliriz.

Hatırlatma olarak, kanal oluşturma şöyle görünüyordu (çözümler dizininde görüldüğü gibi):

```groovy title="solutions/composable-hello/main.nf" linenums="10" hl_lines="4"
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)
        .splitCsv()
        .map { line -> line[0] }
```

Bu yüzden bunu başlatma iş akışına küçük değişikliklerle takmamız yeterli: kanal adını `greeting_ch`'den `ch_samplesheet`'e ve parametre adını `params.greeting`'den `params.input`'a güncelliyoruz (vurgulanan satıra bakın).

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

Bu, girdi işlemenin çalışması için yapmamız gereken değişiklikleri tamamlar.

Mevcut haliyle, bu bize şema doğrulaması için nf-core'un yerleşik yeteneklerinden yararlanma imkanı vermeyecek, ancak bunu daha sonra ekleyebiliriz.
Şimdilik, test verileri üzerinde başarıyla çalıştırabileceğimiz bir şeye ulaşmak için mümkün olduğunca basit tutmaya odaklanıyoruz.

### 4.3. Test profilini güncelleyin

Test verileri ve parametrelerden bahsetmişken, bu pipeline için test profilini şablon tarafından sağlanan örnek samplesheet yerine `greetings.csv` mini-samplesheet'ini kullanacak şekilde güncelleyelim.

`core-hello/conf` altında, küçük bir veri örneğini ve tam boyutlu bir örneği test etmeyi amaçlayan iki şablonlu test profili buluyoruz: `test.config` ve `test_full.config`.
Pipeline'ımızın amacı göz önüne alındığında, tam boyutlu bir test profili kurmanın gerçekten bir anlamı yok, bu nedenle `test_full.config`'i görmezden gelmekten veya silmekten çekinmeyin.
`test.config`'i birkaç varsayılan parametreyle `greetings.csv` dosyamız üzerinde çalışacak şekilde ayarlamaya odaklanacağız.

#### 4.3.1. `greetings.csv` dosyasını kopyalayın

İlk olarak `greetings.csv` dosyasını pipeline projemizde uygun bir yere kopyalamamız gerekiyor.
Tipik olarak küçük test dosyaları `assets` dizininde saklanır, bu nedenle dosyayı çalışma dizinimizden kopyalayalım.

```bash
cp greetings.csv core-hello/assets/.
```

Artık `greetings.csv` dosyası test girdisi olarak kullanılmaya hazır.

#### 4.3.2. `test.config` dosyasını güncelleyin

Şimdi `test.config` dosyasını aşağıdaki gibi güncelleyebiliriz:

=== "Sonra"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-10"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Input data
        input  = "${projectDir}/assets/greetings.csv"

        // Other parameters
        batch     = 'test'
        character = 'tux'
    }
    ```

=== "Önce"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-8"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Input data
        // TODO nf-core: Specify the paths to your test data on nf-core/test-datasets
        // TODO nf-core: Give any required params for the test so that command line flags are not needed
        input  = params.pipelines_testdata_base_path + 'viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
    }
    ```

Önemli noktalar:

- **`${projectDir}` kullanımı**: Bu, ana iş akışı betiğinin bulunduğu dizine (pipeline kökü) işaret eden bir Nextflow örtük değişkenidir. Bunu kullanmak, yolun pipeline'ın nereden çalıştırıldığından bağımsız olarak çalışmasını sağlar.
- **Mutlak yollar**: `${projectDir}` kullanarak, pipeline ile birlikte gelen test verileri için önemli olan mutlak bir yol oluşturuyoruz.
- **Test verisi konumu**: nf-core pipeline'ları tipik olarak test verilerini küçük test dosyaları için pipeline deposu içindeki `assets/` dizininde saklar veya daha büyük dosyalar için harici test veri kümelerine referans verir.

Ve bu arada, bunun çok temel makinelerde (Github Codespaces'teki minimal VM'ler gibi) çalışacağından emin olmak için varsayılan kaynak limitlerini sıkılaştıralım:

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

### 4.4. Pipeline'ı test profiliyle çalıştırın

Bu çok şeydi, ama sonunda pipeline'ı çalıştırmayı deneyebiliriz!
Doğrulamayı henüz kurmadığımız için komut satırına `--validate_params false` eklememiz gerektiğini unutmayın (bu daha sonra gelecek).

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
    [45/bb6096] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [81/7e2e34] CORE_HELLO:HELLO:collectGreetings   [100%] 1 of 1 ✔
    [96/9442a1] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Gördüğünüz gibi, bu başlatma alt iş akışı sayesinde başlangıçta tipik nf-core özetini üretti ve her modül için satırlar artık tam PIPELINE:WORKFLOW:module adlarını gösteriyor.

### 4.5. Pipeline çıktılarını bulun

Şimdi soru şu: pipeline'ın çıktıları nerede?
Ve cevap oldukça ilginç: sonuçları aramak için şimdi iki farklı yer var.

Daha önce hatırlayabileceğiniz gibi, yeni oluşturulan iş akışının ilk çalıştırması çeşitli yürütme raporları ve meta verileri içeren `core-hello-results/` adlı bir dizin üretti.

```bash
tree core-hello-results
```

??? abstract "Dizin içeriği"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_report_2025-11-21_07-29-37.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_07-29-37.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── execution_trace_2025-11-21_07-29-37.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-13.json
        ├── params_2025-11-21_07-29-41.json
        └── pipeline_dag_2025-11-21_04-47-18.html
        └── pipeline_dag_2025-11-21_07-29-37.html

    1 directory, 12 files
    ```

İş akışı hala sadece bir yer tutucu olduğunda ilk çalıştırmadan aldığımız raporlara ek olarak başka bir yürütme raporu seti aldığımızı görüyorsunuz.
Bu sefer beklendiği gibi çalıştırılan tüm görevleri görüyorsunuz.

![Hello pipeline için yürütme zaman çizelgesi raporu](./img/execution_timeline_hello.png)

!!! note "Not"

    Bir kez daha görevler paralel olarak çalıştırılmadı çünkü Github Codespaces'te minimalist bir makine üzerinde çalışıyoruz.
    Bunların paralel çalıştığını görmek için, codespace'inizin CPU tahsisini ve test yapılandırmasındaki kaynak limitlerini artırmayı deneyin.

Bu harika, ama gerçek pipeline sonuçlarımız orada değil!

İşte olan şey: modüllerin kendisinde hiçbir şeyi değiştirmedik, bu nedenle modül seviyesindeki `publishDir` yönergeleri tarafından işlenen çıktılar hala orijinal pipeline'da belirtildiği gibi bir `results` dizinine gidiyor.

```bash
tree results
```

??? abstract "Dizin içeriği"

    ```console
    results
    ├── Bonjour-output.txt
    ├── COLLECTED-test-batch-output.txt
    ├── COLLECTED-test-output.txt
    ├── cowpy-COLLECTED-test-batch-output.txt
    ├── cowpy-COLLECTED-test-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt

    0 directories, 10 files
    ```

Ah, işte oradalar, orijinal Hello pipeline'ının önceki çalıştırmalarının çıktılarıyla karışık.

Demo pipeline'ının çıktıları gibi düzgün bir şekilde organize olmalarını istiyorsak, çıktıların yayınlanacak şekilde nasıl ayarlandığını değiştirmemiz gerekecek.
Bunu bu eğitim kursunda daha sonra nasıl yapacağınızı göstereceğiz.

<!-- TODO: Update this once we've updated Hello Nextflow to use workflow-level outputs -->

Ve işte bu! Orijinal pipeline ile aynı sonucu elde etmek için çok fazla iş gibi görünebilir, ancak otomatik olarak oluşturulan tüm bu güzel raporları alıyorsunuz ve artık girdi doğrulaması ve sonraki bir bölümde ele alacağımız bazı düzgün meta veri işleme yetenekleri dahil olmak üzere nf-core'un ek özelliklerinden yararlanmak için sağlam bir temele sahipsiniz.

---

### Özet

nf-core şablonunu kullanarak normal bir Nextflow pipeline'ını nf-core tarzı bir pipeline'a nasıl dönüştüreceğinizi biliyorsunuz.
Bunun bir parçası olarak, bir iş akışını nasıl birleştirilebilir yapacağınızı ve özel bir nf-core tarzı pipeline geliştirirken en yaygın olarak uyarlanması gereken nf-core şablonunun öğelerini nasıl tanımlayacağınızı öğrendiniz.

### Sırada ne var?

Bir mola verin, bu zor bir işti! Hazır olduğunuzda, nf-core/modules deposundan topluluk tarafından sürdürülen modüllerden nasıl yararlanacağınızı öğrenmek için [Bölüm 3: Bir nf-core modülü kullanın](./03_use_module.md) bölümüne geçin.
