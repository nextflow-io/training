# Bölüm 4: Bir nf-core modülü oluşturun

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } YZ destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core eğitim kursunun bu dördüncü bölümünde, modülleri taşınabilir ve sürdürülebilir kılan temel kuralları uygulayarak bir nf-core modülünün nasıl oluşturulacağını gösteriyoruz.

nf-core projesi, Bölüm 2'de iş akışı için kullandığımıza benzer şekilde, düzgün yapılandırılmış modül şablonlarını otomatik olarak oluşturan bir komut (`nf-core modules create`) sağlar.
Ancak, öğretim amaçları için, bunu manuel olarak yaparak başlayacağız: `core-hello` pipeline'ınızdaki yerel `cowpy` modülünü adım adım nf-core tarzı bir modüle dönüştüreceğiz.
Bundan sonra, gelecekte daha verimli çalışmak için şablon tabanlı modül oluşturmayı nasıl kullanacağınızı göstereceğiz.

??? info "Bu bölümden nasıl başlanır"

    Bu bölüm, [Bölüm 3: Bir nf-core modülü kullanın](./03_use_module.md) bölümünü tamamladığınızı ve `CAT_CAT` modülünü pipeline'ınıza entegre ettiğinizi varsayar.

    Bölüm 3'ü tamamlamadıysanız veya bu bölüm için yeni başlamak istiyorsanız, başlangıç noktanız olarak `core-hello-part3` çözümünü kullanabilirsiniz.
    Bu komutları `hello-nf-core/` dizininin içinden çalıştırın:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    Bu size `CAT_CAT` modülünün zaten entegre edildiği bir pipeline verir.
    Aşağıdaki komutu çalıştırarak başarıyla çalıştığını test edebilirsiniz:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. `cowpy`'yi bir nf-core modülüne dönüştürün

Bu bölümde, `core-hello` pipeline'ınızdaki yerel `cowpy` modülüne nf-core kurallarını uygulayacak ve onu nf-core topluluk standartlarını takip eden bir modüle dönüştüreceğiz.

`cowpy` süreç modülü için mevcut kod şudur:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
process cowpy {

    publishDir 'results', mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """
}
```

Aşağıdaki nf-core kurallarını aşamalı olarak uygulayacağız:

1. **Süreç adını `COWPY` olarak büyük harfe çevirin** - kurala uymak için.
2. **`COWPY`'yi metadata demetlerini kullanacak şekilde güncelleyin** - örnek metadata'sını iş akışı boyunca yaymak için.
3. **Araç argüman yapılandırmasını `ext.args` ile merkezileştirin** - arayüzü minimal tutarken modül çok yönlülüğünü artırmak için.
4. **Çıktı adlandırmasını `ext.prefix` ile standartlaştırın** - tutarlılığı teşvik etmek için.
5. **Yayınlama yapılandırmasını merkezileştirin** - tutarlılığı teşvik etmek için.

Her adımdan sonra, her şeyin beklendiği gibi çalıştığını test etmek için pipeline'ı çalıştıracağız.

!!! warning "Çalışma dizini"

    Bu bölümdeki tüm dosya düzenlemeleri ve komut yürütmeleri için `core-hello` dizininde (pipeline kök dizininizde) olduğunuzdan emin olun.

    ```bash
    cd core-hello
    ```

### 1.1. Süreç adını büyük harfe çevirin

Bu tamamen stilistik bir kuraldır (teknik bir gerekçe yoktur) ancak nf-core modülleri için norm olduğundan, uyalım.

Üç grup değişiklik yapmamız gerekiyor:

1. Modüldeki süreç adını güncelleyin
2. İş akışı başlığındaki modül import ifadesini güncelleyin
3. İş akışı gövdesindeki süreç çağrısını ve emit bildirimini güncelleyin

Hadi başlayalım!

#### 1.1.1. Modüldeki süreç adını güncelleyin

`cowpy.nf` modül dosyasını açın (`core-hello/modules/local/` altında) ve süreç adını büyük harfe değiştirin:

=== "Sonra"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "Önce"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

Bu durumda büyük harfe çevirme tamamen basittir.

Süreç adı birkaç kelimeden oluşuyorsa, örneğin başlangıçta camel case olarak MyCowpyTool adında bir sürecimiz olsaydı, nf-core kuralı bunları ayırmak için alt çizgi kullanmak olurdu, bu da MY_COWPY_TOOL'u verir.

#### 1.1.2. Modül import ifadesini güncelleyin

Süreç adları büyük/küçük harfe duyarlıdır, bu nedenle süreç adını değiştirdiğimize göre, `hello.nf` iş akışı başlığındaki modül import ifadesini buna göre güncellememiz gerekiyor:

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { cowpy                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Sürece yapılan çağrıları güncellemeyi önlemek için import ifadesinde bir takma ad kullanabiliriz, ancak bu büyük harfe çevirme kuralını benimsemenin amacını bir şekilde boşa çıkarır.

#### 1.1.3. Süreç çağrısını ve emit bildirimini güncelleyin

Şimdi `hello.nf` iş akışı bloğundaki sürece yapılan iki referansı güncelleyelim:

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // generate ASCII art of the greetings with cowpy
    COWPY(CAT_CAT.out.file_out)

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
    cowpy_hellos   = COWPY.out.cowpy_output
    versions       = ch_versions
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // generate ASCII art of the greetings with cowpy
    cowpy(CAT_CAT.out.file_out)

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
    cowpy_hellos   = cowpy.out.cowpy_output
    versions       = ch_versions
    ```

**Her iki** değişikliği de yaptığınızdan emin olun, aksi takdirde bunu çalıştırdığınızda bir hata alırsınız.

#### 1.1.4. Test etmek için pipeline'ı çalıştırın

Bu değişikliklerden sonra her şeyin doğru çalıştığını test etmek için iş akışını çalıştıralım.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [elegant_plateau] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2026-01-06_04-51-29

    Core Nextflow options
      runName                   : elegant_plateau
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [7b/66ceb5] CORE_HELLO:HELLO:sayHello (3)       | 3 of 3 ✔
    [8e/1bafb9] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [bb/203575] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [39/715489] CORE_HELLO:HELLO:COWPY              | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Pekala, bu çalışıyor! Şimdi daha önemli değişiklikler yapmaya geçelim.

### 1.2. `COWPY`'yi metadata demetlerini kullanacak şekilde güncelleyin

`core-hello` pipeline'ının mevcut sürümünde, aşağıdaki diyagramın üst yarısında gösterildiği gibi, `COWPY`'ye geçmek için `CAT_CAT`'in çıktı demetinden dosyayı çıkarıyoruz.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

`COWPY`'nin metadata demetlerini doğrudan kabul etmesi daha iyi olurdu, bu da metadata'nın diyagramın alt yarısında gösterildiği gibi iş akışı boyunca akmasına izin verir.

Bu amaçla, aşağıdaki değişiklikleri yapmamız gerekecek:

1. Girdi ve çıktı tanımlarını güncelleyin
2. İş akışındaki süreç çağrısını güncelleyin
3. İş akışındaki emit bloğunu güncelleyin

Tüm bunları yaptıktan sonra, her şeyin daha önce olduğu gibi çalıştığını test etmek için pipeline'ı çalıştıracağız.

#### 1.2.1. Girdi ve çıktı tanımlarını güncelleyin

`cowpy.nf` modül dosyasına dönün ve aşağıda gösterildiği gibi metadata demetlerini kabul edecek şekilde değiştirin.

=== "Sonra"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output
    ```

=== "Önce"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"
    ```

Gördüğünüz gibi, hem **ana girdiyi** hem de **çıktıyı** Bölüm 3'te tanıtılan `tuple val(meta), path(input_file)` desenini takip eden bir demete değiştirdik.
Çıktı için, çıktı kanalına açıklayıcı bir ad vermek için `emit: cowpy_output` ekledik.

Artık sürecin ne beklediğini değiştirdiğimize göre, süreç çağrısında ona ne sağladığımızı güncellememiz gerekiyor.

#### 1.2.2. İş akışındaki süreç çağrısını güncelleyin

İyi haber şu ki bu değişiklik süreç çağrısını basitleştirecek.
Artık `CAT_CAT`'in çıktısı ve `COWPY`'nin girdisi aynı 'şekilde', yani her ikisi de `tuple val(meta), path(input_file)` yapısından oluştuğundan, `CAT_CAT` sürecinin çıktısından dosyayı açıkça çıkarmak zorunda kalmak yerine onları doğrudan bağlayabiliriz.

`hello.nf` iş akışı dosyasını açın (`core-hello/workflows/` altında) ve `COWPY` çağrısını aşağıda gösterildiği gibi güncelleyin.

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generate ASCII art of the greetings with cowpy
        COWPY(ch_for_cowpy, params.character)
    ```

Artık `COWPY`'yi doğrudan `CAT_CAT.out.file_out` üzerinde çağırıyoruz.

Sonuç olarak, artık `ch_for_cowpy` kanalını oluşturmamıza gerek yok, bu nedenle bu satır (ve yorum satırı) tamamen silinebilir.

#### 1.2.3. İş akışındaki emit bloğunu güncelleyin

`COWPY` artık adlandırılmış bir çıktı olan `cowpy_output`'u yayınladığından, `hello.nf` iş akışının `emit:` bloğunu bunu kullanacak şekilde güncelleyebiliriz.

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out.cowpy_output
        versions       = ch_versions
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out
        versions       = ch_versions
    ```

Bu teknik olarak gerekli değildir, ancak mümkün olduğunda adlandırılmış çıktılara başvurmak iyi bir uygulamadır.

#### 1.2.4. Test etmek için pipeline'ı çalıştırın

Bu değişikliklerden sonra her şeyin doğru çalıştığını test etmek için iş akışını çalıştıralım.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [modest_saha] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-16-55

    Core Nextflow options
      runName                   : modest_saha
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [a8/447993] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [00/1fc59c] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [57/ac800d] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [b7/092f2b] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Pipeline başarıyla çalışmalı, metadata artık `CAT_CAT`'ten `COWPY`'ye akıyor.

Bu, `COWPY`'nin metadata demetlerini işlemesi için yapmamız gerekenleri tamamlıyor.
Şimdi, nf-core modül desenlerinden yararlanmak için başka neler yapabileceğimize bakalım.

### 1.3. Araç argüman yapılandırmasını `ext.args` ile merkezileştirin

Mevcut durumunda, `COWPY` süreci `character` parametresi için bir değer almayı bekler.
Sonuç olarak, araç tarafından belirlenen varsayılanlardan memnun olsak bile, süreci her çağırdığımızda bir değer sağlamamız gerekir.
`COWPY` için bu kabul edilebilir bir sorun değildir, ancak birçok isteğe bağlı parametreye sahip araçlar için oldukça zahmetli olabilir.

nf-core projesi, araç argümanlarını yapılandırma dosyaları aracılığıyla daha rahat yönetmek için [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) adlı bir Nextflow özelliğini kullanmayı önerir.

Her araç seçeneği için süreç girdileri bildirmek yerine, modülü komut satırının oluşturulmasında `ext.args`'a başvuracak şekilde yazarsınız.
Ardından, tüm modüller için yapılandırma ayrıntılarını birleştiren `modules.config` dosyasında kullanmak istediğiniz argümanları ve değerleri tutacak şekilde `ext.args` değişkenini ayarlamak yeterlidir.
Nextflow, bu argümanları değerleriyle birlikte çalışma zamanında araç komut satırına ekleyecektir.

Bu yaklaşımı `COWPY` modülüne uygulayalım.
Aşağıdaki değişiklikleri yapmamız gerekecek:

1. `COWPY` modülünü güncelleyin
2. `modules.config` dosyasında `ext.args`'ı yapılandırın
3. `hello.nf` iş akışını güncelleyin

Tüm bunları yaptıktan sonra, her şeyin daha önce olduğu gibi çalıştığını test etmek için pipeline'ı çalıştıracağız.

#### 1.3.1. `COWPY` modülünü güncelleyin

Hadi yapalım.
`cowpy.nf` modül dosyasını açın (`core-hello/modules/local/` altında) ve aşağıda gösterildiği gibi `ext.args`'a başvuracak şekilde değiştirin.

=== "Sonra"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

=== "Önce"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

Üç değişiklik yaptığımızı görebilirsiniz.

1. **`input:` bloğunda, `val character` girdisini kaldırdık.**
   İleriye dönük olarak, bu argümanı aşağıda açıklandığı gibi `ext.args` yapılandırması aracılığıyla sağlayacağız.

2. **`script:` bloğunda, `def args = task.ext.args ?: ''` satırını ekledik.**
   Bu satır, `args` değişkeninin değerini belirlemek için `?:` operatörünü kullanır: boş değilse `task.ext.args`'ın içeriği veya boşsa boş bir dize.
   Genellikle `ext.args`'a atıfta bulunsak da, bu kod modül düzeyindeki `ext.args` yapılandırmasını çıkarmak için `task.ext.args`'a başvurmalıdır.

3. **Komut satırında, `-c "$character"`'ı `$args` ile değiştirdik.**
   Nextflow'un `modules.config` dosyasındaki `ext.args`'ta ayarlanan araç argümanlarını enjekte edeceği yer burasıdır.

Sonuç olarak, modül arayüzü artık daha basittir: yalnızca temel metadata ve dosya girdilerini bekler.

!!! note

    `?:` operatörüne genellikle 'Elvis operatörü' denir çünkü `?` karakteri saçındaki dalgayı simgeleyen, yana dönük bir Elvis Presley yüzüne benzer.

#### 1.3.2. `modules.config` dosyasında `ext.args`'ı yapılandırın

Artık `character` bildirimini modülden çıkardığımıza göre, onu `modules.config` yapılandırma dosyasındaki `ext.args`'a eklememiz gerekiyor.

Özellikle, `process {}` bloğuna şu küçük kod parçasını ekleyeceğiz:

```groovy title="Eklenecek kod"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

`withName:` sözdizimi bu yapılandırmayı yalnızca `COWPY` sürecine atar ve `ext.args = { "-c ${params.character}" }` basitçe `character` parametresinin değerini içerecek bir dize oluşturur.
Nextflow'a parametrenin değerini çalışma zamanında değerlendirmesini söyleyen süslü parantezlerin kullanımına dikkat edin.

Mantıklı mı? Hadi ekleyelim.

`conf/modules.config` dosyasını açın ve yapılandırma kodunu aşağıda gösterildiği gibi `process {}` bloğunun içine ekleyin.

=== "Sonra"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]

        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "Önce"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    ```

Umarım bir pipeline'daki tüm modüllerin `ext.args`'larının bu dosyada belirtildiğini hayal edebilirsiniz, aşağıdaki faydalarla:

- **Modül arayüzü basit kalır** - Yalnızca temel metadata ve dosya girdilerini kabul eder
- **Pipeline hala `params.character`'ı açığa çıkarır** - Son kullanıcılar bunu daha önce olduğu gibi yapılandırabilir
- **Modül artık taşınabilir** - Belirli bir parametre adı beklemeden diğer pipeline'larda yeniden kullanılabilir
- **Yapılandırma merkezileştirilmiştir** - `modules.config` dosyasında, iş akışı mantığını temiz tutar

`modules.config` dosyasını tüm pipeline'ların modül başına yapılandırmayı merkezileştirdiği yer olarak kullanarak, modüllerimizi farklı pipeline'lar arasında daha yeniden kullanılabilir hale getiririz.

#### 1.3.3. `hello.nf` iş akışını güncelleyin

`COWPY` modülü artık `character` parametresini girdi olarak gerektirmediğinden, iş akışı çağrısını buna göre güncellememiz gerekiyor.

`hello.nf` iş akışı dosyasını açın (`core-hello/workflows/` altında) ve `COWPY` çağrısını aşağıda gösterildiği gibi güncelleyin.

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

İş akışı kodu artık daha temiz: `params.character`'ı doğrudan sürece geçirmemize gerek yok.
Modül arayüzü minimal tutulur, bu da onu daha taşınabilir hale getirir, pipeline ise yapılandırma yoluyla açık seçeneği sağlamaya devam eder.

#### 1.3.4. Test etmek için pipeline'ı çalıştırın

İş akışının hala beklendiği gibi çalıştığını test edelim, `ext.args` yapılandırmasının çalıştığını doğrulamak için farklı bir karakter belirterek.

Daha... esrarengiz seçeneklerden biri olan `kosh`'u kullanarak bu komutu çalıştırın:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [exotic_planck] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-23-13

    Core Nextflow options
      runName                   : exotic_planck
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [13/9e3c0e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [e2/5b0ee5] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b6/4fb569] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [38/eb29ea] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Bu daha önce olduğu gibi başarıyla çalışmalıdır.

`ext.args` yapılandırmasının çalıştığını doğrulamak için çıktıyı kontrol edelim.
Dosya tarayıcısında çıktıyı bulun veya çıktı dosyasına bakmak için görev hash'ini (yukarıdaki örnekte `38/eb29ea` kısmı) kullanın:

```bash
cat work/38/eb29ea*/cowpy-test.txt
```

??? success "Komut çıktısı"

    ```console
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
        \
        \
          \
      ___       _____     ___
    /   \     /    /|   /   \
    |     |   /    / |  |     |
    |     |  /____/  |  |     |
    |     |  |    |  |  |     |
    |     |  | {} | /   |     |
    |     |  |____|/    |     |
    |     |    |==|     |     |
    |      \___________/      |
    |                         |
    |                         |
    ```

ASCII sanatının `kosh` karakteriyle görüntülendiğini görmelisiniz, bu da `ext.args` yapılandırmasının çalıştığını doğrular!

??? info "(İsteğe bağlı) Komut dosyasını inceleyin"

    Yapılandırmanın tam olarak nasıl uygulandığını görmek istiyorsanız, `.command.sh` dosyasını inceleyebilirsiniz:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    `cowpy` komutunu `-c kosh` argümanıyla göreceksiniz:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    Bu, `.command.sh` dosyasının `ext.args` yapılandırmasına göre doğru şekilde oluşturulduğunu gösterir.

Burada neyi başardığımızı düşünmek için bir dakikanızı ayırın.
Bu yaklaşım, modül arayüzünü temel veriler (dosyalar, metadata ve zorunlu örnek başına parametreler) üzerinde odaklanmış tutarken, aracın davranışını kontrol eden seçenekler yapılandırma yoluyla ayrı olarak ele alınır.

Bu, `cowpy` gibi basit bir araç için gereksiz görünebilir, ancak çok sayıda isteğe bağlı argümana sahip veri analizi araçları için büyük bir fark yaratabilir.

Bu yaklaşımın faydalarını özetlemek gerekirse:

- **Temiz arayüz**: Modül temel veri girdilerine (metadata ve dosyalar) odaklanır
- **Esneklik**: Kullanıcılar yapılandırma yoluyla araç argümanlarını belirtebilir, örneğe özgü değerler dahil
- **Tutarlılık**: Tüm nf-core modülleri bu deseni takip eder
- **Taşınabilirlik**: Modüller sabit kodlanmış araç seçenekleri olmadan yeniden kullanılabilir
- **İş akışı değişikliği yok**: Araç seçeneklerini eklemek veya değiştirmek iş akışı kodunu güncellemeyi gerektirmez

!!! note

    `ext.args` sistemi, burada ele alınmayan, metadata'ya göre argüman değerlerini dinamik olarak değiştirme dahil olmak üzere güçlü ek yeteneklere sahiptir. Daha fazla ayrıntı için [nf-core modül spesifikasyonlarına](https://nf-co.re/docs/guidelines/components/modules) bakın.

### 1.4. Çıktı adlandırmasını `ext.prefix` ile standartlaştırın

Artık `COWPY` sürecine metamap'e erişim verdiğimize göre, başka bir yararlı nf-core deseninden yararlanmaya başlayabiliriz: metadata'ya dayalı çıktı dosyalarını adlandırma.

Burada, `meta.id` (metamap'e dahil edilen tanımlayıcı) kullanarak modüller arasında çıktı dosyası adlandırmasını standartlaştırmamıza izin verecek `ext.prefix` adlı bir Nextflow özelliğini kullanacağız, aynı zamanda istenirse modülleri ayrı ayrı yapılandırabilme yeteneğini koruyacağız.

Bu, `ext.args` ile yaptığımıza benzer olacak, ilerledikçe detaylandıracağımız birkaç farkla.

Bu yaklaşımı `COWPY` modülüne uygulayalım.
Aşağıdaki değişiklikleri yapmamız gerekecek:

1. `COWPY` modülünü güncelleyin
2. `modules.config` dosyasında `ext.prefix`'i yapılandırın

(İş akışında değişiklik gerekmez.)

Bunu yaptıktan sonra, her şeyin daha önce olduğu gibi çalıştığını test etmek için pipeline'ı çalıştıracağız.

#### 1.4.1. `COWPY` modülünü güncelleyin

`cowpy.nf` modül dosyasını açın (`core-hello/modules/local/` altında) ve aşağıda gösterildiği gibi `ext.prefix`'e başvuracak şekilde değiştirin.

=== "Sonra"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 6 8"
        output:
            tuple val(meta), path("${prefix}.txt"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat $input_file | cowpy $args > ${prefix}.txt
        """
    }
    ```

=== "Önce"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 7"
        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

Üç değişiklik yaptığımızı görebilirsiniz.

1. **`script:` bloğunda, `prefix = task.ext.prefix ?: "${meta.id}"` satırını ekledik.**
   Bu satır, `prefix` değişkeninin değerini belirlemek için `?:` operatörünü kullanır: boş değilse `task.ext.prefix`'in içeriği veya boşsa metamap'ten tanımlayıcı (`meta.id`).
   Genellikle `ext.prefix`'e atıfta bulunsak da, bu kod modül düzeyindeki `ext.prefix` yapılandırmasını çıkarmak için `task.ext.prefix`'e başvurmalıdır.

2. **Komut satırında, `cowpy-${input_file}`'ı `${prefix}.txt` ile değiştirdik.**
   Nextflow'un yukarıdaki satır tarafından belirlenen `prefix` değerini enjekte edeceği yer burasıdır.

3. **`output:` bloğunda, `path("cowpy-${input_file}")`'ı `path("${prefix}.txt")` ile değiştirdik.**
   Bu, komut satırında yazılana göre dosya yolunun ne olacağını yeniden belirtir.

Sonuç olarak, çıktı dosyası adı artık mantıklı bir varsayılan (metamap'ten tanımlayıcı) ve uygun dosya formatı uzantısı kullanılarak oluşturulur.

#### 1.4.2. `modules.config` dosyasında `ext.prefix`'i yapılandırın

Bu durumda mantıklı varsayılan bizim zevkimiz için yeterince açıklayıcı değil; daha önce olduğu gibi araç adını içeren özel bir adlandırma deseni kullanmak istiyoruz, `cowpy-<id>.txt`.

Bunu, `ext.args` ile `character` parametresi için yaptığımız gibi `modules.config` dosyasında `ext.prefix`'i yapılandırarak yapacağız, bu sefer `withName: 'COWPY' {}` bloğu zaten mevcut ve sadece aşağıdaki satırı eklememiz gerekiyor:

```groovy title="Eklenecek kod"
ext.prefix = { "cowpy-${meta.id}" }
```

Bu istediğimiz dizeyi oluşturacaktır.
Bir kez daha, bu sefer Nextflow'a `meta.id` değerini çalışma zamanında değerlendirmesini söylemek için süslü parantezler kullandığımıza dikkat edin.

Hadi ekleyelim.

`conf/modules.config` dosyasını açın ve yapılandırma kodunu aşağıda gösterildiği gibi `process {}` bloğunun içine ekleyin.

=== "Sonra"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Önce"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    ```

Merak ediyorsanız, `ext.prefix` closure'ı doğru metadata parçasına erişebilir çünkü yapılandırma, metadata'nın mevcut olduğu süreç yürütme bağlamında değerlendirilir.

#### 1.4.3. Test etmek için pipeline'ı çalıştırın

İş akışının hala beklendiği gibi çalıştığını test edelim.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [admiring_turing] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-29-02

    Core Nextflow options
      runName                   : admiring_turing
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [b2/e08524] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [13/88939f] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [23/4554e1] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [a3/c6cbe9] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Sonuçlar dizinindeki çıktıya bakın.
Cowpy çıktı dosyasını daha önceki gibi aynı adlandırmayla görmelisiniz: varsayılan batch adına dayalı `cowpy-test.txt`.

??? abstract "Dizin içeriği"

    ```console hl_lines="3"
    results
    ├── Bonjour-output.txt
    ├── cowpy-test.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Modül veya iş akışı kodunda herhangi bir değişiklik yapmak zorunda kalmadan adlandırma desenini değiştirebileceğinizden emin olmak için `conf/modules.config` dosyasındaki `ext.prefix` yapılandırmasını değiştirmekten çekinmeyin.

Alternatif olarak, bu kısmın hala anında özelleştirilebilir olduğundan emin olmak için komut satırında farklı bir `--batch` parametresi belirterek bunu tekrar çalıştırmayı da deneyebilirsiniz.

Bu, `ext.prefix`'in modül arayüzünü esnek tutarken tercih ettiğiniz adlandırma kuralını korumanıza nasıl izin verdiğini gösterir.

Bu yaklaşımın faydalarını özetlemek gerekirse:

- **Standartlaştırılmış adlandırma**: Çıktı dosyaları genellikle metadata'dan örnek kimlikleri kullanılarak adlandırılır
- **Yapılandırılabilir**: Kullanıcılar gerekirse varsayılan adlandırmayı geçersiz kılabilir
- **Tutarlı**: Tüm nf-core modülleri bu deseni takip eder
- **Öngörülebilir**: Çıktı dosyalarının ne olarak adlandırılacağını bilmek kolaydır

Oldukça iyi, değil mi?
Peki, modülümüzü nf-core yönergelerine uyacak şekilde geliştirmek için yapmamız gereken bir önemli değişiklik daha var.

### 1.5. Yayınlama yapılandırmasını merkezileştirin

Çıktıları iki farklı dizine yayınladığımızı fark etmiş olabilirsiniz:

- **`results`** — Yerel modüllerimiz için başlangıçtan beri kullandığımız, modül başına `publishDir` yönergeleri kullanılarak ayrı ayrı ayarlanan orijinal çıktı dizini;
- **`core-hello-results`** — Komut satırında `--outdir` ile ayarlanan, nf-core günlüklerini ve `CAT_CAT` tarafından yayınlanan sonuçları alan çıktı dizini.

Bu dağınık ve optimal değil; her şey için tek bir konuma sahip olmak daha iyi olurdu.
Elbette, yerel modüllerimizin her birine gidip `publishDir` yönergesini `core-hello-results` dizinini kullanacak şekilde manuel olarak güncelleyebiliriz, ancak bir dahaki sefere çıktı dizinini değiştirmeye karar verirsek ne olacak?

Bireysel modüllerin yayınlama kararları alması açıkça gidilecek yol değil, özellikle aynı modülün birçok farklı pipeline'da, farklı ihtiyaçları veya tercihleri olan insanlar tarafından kullanılabileceği bir dünyada.
Çıktıların nerede yayınlanacağını iş akışı yapılandırması düzeyinde kontrol edebilmek istiyoruz.

"Hey," diyebilirsiniz, "`CAT_CAT` çıktılarını `--outdir`'e gönderiyor. Belki `publishDir` yönergesini kopyalamalıyız?"

Evet, bu harika bir fikir.

Ancak bir `publishDir` yönergesi yok. (Devam edin, modül koduna bakın.)

Bunun nedeni, nf-core pipeline'larının bireysel modüllerde değil `conf/modules.config` dosyasında `publishDir`'i yapılandırarak iş akışı düzeyinde kontrolü merkezileştirmesidir.
Özellikle, nf-core şablonu, geçersiz kılan bir yönerge sağlanmadıkça tüm modüllere uygulanan varsayılan bir `publishDir` yönergesi (önceden tanımlanmış bir dizin yapısıyla) bildirir.

Bu harika görünmüyor mu? Bu varsayılan yönergeden yararlanmak için tek yapmamız gereken yerel modüllerimizden mevcut `publishDir` yönergesini kaldırmak olabilir mi?

Bunu `COWPY` üzerinde deneyelim ve ne olduğunu görelim, ardından nasıl çalıştığını anlamak için varsayılan yapılandırmanın koduna bakacağız.

Son olarak, istenirse varsayılan davranışı nasıl geçersiz kılacağımızı göstereceğiz.

#### 1.5.1. `COWPY`'den `publishDir` yönergesini kaldırın

Hadi bunu yapalım.
`cowpy.nf` modül dosyasını açın (`core-hello/modules/local/` altında) ve aşağıda gösterildiği gibi `publishDir` yönergesini kaldırın.

=== "Sonra"

    ```groovy title="core-hello/modules/local/cowpy.nf (alıntı)" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "Önce"

    ```groovy title="core-hello/modules/local/cowpy.nf (alıntı)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

Bu kadar!

#### 1.5.2. Test etmek için pipeline'ı çalıştırın

Şimdi pipeline'ı çalıştırırsak ne olduğuna bakalım.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [silly_caravaggio] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-35-56

    Core Nextflow options
      runName                   : silly_caravaggio
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [db/39978e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [b5/bf6a8d] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b7/c61842] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [46/5839d6] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Mevcut çalışma dizininize bakın.
Şimdi `core-hello-results` ayrıca `COWPY` modülünün çıktılarını içeriyor.

??? abstract "Dizin içeriği"

    ```console hl_lines="4-5"
    core-hello-results/
    ├── cat
    │   └── test.txt
    ├── cowpy
    │   └── cowpy-test.txt
    └── pipeline_info
        ├── execution_report_2025-12-27_06-16-55.html
        ├── execution_report_2025-12-27_06-23-13.html
        ├── execution_report_2025-12-27_06-29-02.html
        ├── execution_report_2025-12-27_06-35-56.html
        ├── execution_timeline_2025-12-27_06-16-55.html
        ├── execution_timeline_2025-12-27_06-23-13.html
        ├── execution_timeline_2025-12-27_06-29-02.html
        ├── execution_timeline_2025-12-27_06-35-56.html
        ├── execution_trace_2025-12-27_06-16-55.txt
        ├── execution_trace_2025-12-27_06-23-13.txt
        ├── execution_trace_2025-12-27_06-29-02.txt
        ├── execution_trace_2025-12-27_06-35-56.txt
        ├── hello_software_versions.yml
        ├── params_2025-12-27_06-17-00.json
        ├── params_2025-12-27_06-23-17.json
        ├── params_2025-12-27_06-29-07.json
        ├── params_2025-12-27_06-36-01.json
        ├── pipeline_dag_2025-12-27_06-16-55.html
        ├── pipeline_dag_2025-12-27_06-23-13.html
        ├── pipeline_dag_2025-12-27_06-29-02.html
        └── pipeline_dag_2025-12-27_06-35-56.html
    ```

Nextflow'un iş akışının ve modülün adlarına göre bu dizin hiyerarşisini oluşturduğunu görebilirsiniz.

Sorumlu kod `conf/modules.config` dosyasında yaşar.
Bu, nf-core şablonunun bir parçası olan ve tüm süreçlere uygulanan varsayılan `publishDir` yapılandırmasıdır:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

Bu karmaşık görünebilir, bu yüzden üç bileşenin her birine bakalım:

- **`path:`** Süreç adına göre çıktı dizinini belirler.
  `task.process` içinde bulunan bir sürecin tam adı, iş akışı ve modül içe aktarmalarının hiyerarşisini içerir (`CORE_HELLO:HELLO:CAT_CAT` gibi).
  `tokenize` işlemleri bu hiyerarşiyi kaldırarak sadece süreç adını alır, ardından herhangi bir alt çizgiden önceki ilk kısmı alır (varsa) ve küçük harfe dönüştürür.
  Bu, `CAT_CAT`'in sonuçlarının `${params.outdir}/cat/` dizinine yayınlanmasını belirleyen şeydir.
- **`mode:`** Dosyaların nasıl yayınlanacağını kontrol eder (kopyalama, sembolik bağlantı vb.).
  Bu, `params.publish_dir_mode` parametresi aracılığıyla yapılandırılabilir.
- **`saveAs:`** Hangi dosyaların yayınlanacağını filtreler.
  Bu örnek, `versions.yml` dosyaları için `null` döndürerek bunların yayınlanmasını önler.

Bu, çıktıları düzenlemek için tutarlı bir mantık sağlar.

Bir pipeline'daki tüm modüller bu kuralı benimsediğinde çıktı daha da iyi görünür, bu nedenle pipeline'ınızdaki diğer modüllerden `publishDir` yönergelerini silmekten çekinmeyin.
Bu varsayılan, nf-core yönergelerini takip etmek için açıkça değiştirmediğimiz modüllere bile uygulanacaktır.

Bununla birlikte, girdilerinizi farklı şekilde düzenlemeye karar verebilirsiniz ve iyi haber şu ki bunu yapmak kolaydır.

#### 1.5.3. Varsayılanı geçersiz kılın

Varsayılan `publishDir` yönergesini geçersiz kılmak için, `conf/modules.config` dosyasına kendi yönergelerinizi ekleyebilirsiniz.

Örneğin, 'COWPY' süreci için özel bir `publishDir` yönergesi eklediğimiz bu örnekte olduğu gibi, `withName:` seçicisini kullanarak tek bir süreç için varsayılanı geçersiz kılabilirsiniz.

```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    withName: 'COWPY' {
        ext.args = { "-c ${params.character}" }
        publishDir = [
            path: 'my_custom_results'
        ]
    }
}
```

Aslında bu değişikliği yapmayacağız, ancak bununla oynamaktan ve hangi mantığı uygulayabileceğinizi görmekten çekinmeyin.

Mesele şu ki, bu sistem size her iki dünyanın da en iyisini verir: varsayılan olarak tutarlılık ve talep üzerine yapılandırmayı özelleştirme esnekliği.

Özetlemek gerekirse, şunları elde edersiniz:

- **Tek doğruluk kaynağı**: Tüm yayınlama yapılandırması `modules.config` dosyasında yaşar
- **Yararlı varsayılan**: Süreçler modül başına yapılandırma olmadan kutudan çıkar çıkmaz çalışır
- **Kolay özelleştirme**: Yayınlama davranışını modül kodunda değil yapılandırmada geçersiz kılın
- **Taşınabilir modüller**: Modüller çıktı konumlarını sabit kodlamaz

Bu, kesinlikle kullanmayı öğrenmeniz gereken nf-core modül özelliklerini tamamlar, ancak [nf-core modül spesifikasyonlarında](https://nf-co.re/docs/guidelines/components/modules) okuyabileceğiniz başkaları da vardır.

### Özet

Artık yerel modülleri nf-core kurallarını takip edecek şekilde nasıl uyarlayacağınızı biliyorsunuz:

- Modüllerinizi metadata demetlerini kabul edecek ve yayacak şekilde tasarlayın;
- Modül arayüzlerini minimal ve taşınabilir tutmak için `ext.args` kullanın;
- Yapılandırılabilir, standartlaştırılmış çıktı dosyası adlandırması için `ext.prefix` kullanın;
- Tutarlı bir sonuçlar dizin yapısı için varsayılan merkezileştirilmiş `publishDir` yönergesini benimseyin.

### Sırada ne var?

Modülleri kolay yoldan oluşturmak için nf-core'un yerleşik şablon tabanlı araçlarını nasıl kullanacağınızı öğrenin.

---

## 2. nf-core araçlarıyla bir modül oluşturun

Artık nf-core modül desenlerini manuel olarak uygulayarak öğrendiğinize göre, pratikte modülleri nasıl oluşturacağınıza bakalım.

### 2.1. Bir şablondan modül iskeleti oluşturun

Pipeline'lar oluşturmak için var olana benzer şekilde, nf-core projesi, tüm bu desenler başlangıçtan itibaren yerleşik olarak, bir şablona dayalı düzgün yapılandırılmış modüller oluşturmak için araçlar sağlar.

#### 2.1.1. Modül oluşturma komutunu çalıştırın

`nf-core modules create` komutu, öğrendiğiniz tüm kuralları zaten takip eden bir modül şablonu oluşturur.

Bu komutu çalıştırarak minimal bir şablonla `COWPY` modülünün yeni bir sürümünü oluşturalım:

```bash
nf-core modules create --empty-template COWPY
```

`--empty-template` bayrağı, temel yapıyı görmeyi kolaylaştıran, ekstra kod olmadan temiz bir başlangıç şablonu oluşturur.

Komut etkileşimli olarak çalışır ve kurulum boyunca size rehberlik eder.
Metadata'yı önceden doldurmak için Bioconda ve bio.tools gibi paket depolarından araç bilgilerini otomatik olarak arar.

Birkaç yapılandırma seçeneği için isteneceksiniz:

- **Yazar bilgileri**: Atıf için GitHub kullanıcı adınız
- **Kaynak etiketi**: Önceden tanımlanmış bir hesaplama gereksinimleri seti.
  nf-core projesi, hafif araçlar için `process_single` ve zorlu araçlar için `process_high` gibi standart etiketler sağlar.
  Bu etiketler, farklı yürütme ortamlarında kaynak tahsisini yönetmeye yardımcı olur.
- **Metadata gereksinimi**: Modülün bir `meta` haritası aracılığıyla örneğe özgü bilgilere ihtiyacı olup olmadığı (veri işleme modülleri için genellikle evet).

Araç, paket bilgilerini bulma ve yapıyı kurma karmaşıklığını ele alır, aracın belirli mantığını uygulamaya odaklanmanıza olanak tanır.

#### 2.1.2. Modül iskeletini inceleyin

Araç, `modules/local/` dizininde (veya nf-core/modules deposundaysanız `modules/nf-core/` dizininde) tam bir modül yapısı oluşturur:

??? abstract "Dizin içeriği"

    ```console
    modules/local/cowpy
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
    ```

Her dosya belirli bir amaca hizmet eder:

- **`main.nf`**: Tüm nf-core desenleri yerleşik olarak süreç tanımı
- **`meta.yml`**: Girdileri, çıktıları ve aracı açıklayan modül belgeleri
- **`environment.yml`**: Bağımlılıklar için Conda ortam spesifikasyonu
- **`tests/main.nf.test`**: Modülün çalıştığını doğrulamak için nf-test test durumları

!!! tip "Test hakkında daha fazla bilgi edinin"

    Oluşturulan test dosyası, Nextflow pipeline'ları ve modülleri için bir test çerçevesi olan nf-test kullanır. Bu testleri nasıl yazacağınızı ve çalıştıracağınızı öğrenmek için [nf-test yan görevine](../side_quests/nf-test.md) bakın.

Oluşturulan `main.nf`, az önce öğrendiğiniz tüm desenleri ve bazı ek özellikleri içerir:

```groovy title="modules/local/cowpy/main.nf" hl_lines="11 21 22"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(input)        // Desen 1: Metadata demetleri ✓

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                // Desen 2: ext.args ✓
    def prefix = task.ext.prefix ?: "${meta.id}"  // Desen 3: ext.prefix ✓

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
}
```

Yukarıda manuel olarak uyguladığınız tüm desenlerin zaten orada olduğuna dikkat edin!

Şablon ayrıca birkaç ek nf-core kuralı içerir.
Bunlardan bazıları kutudan çıkar çıkmaz çalışırken, diğerleri aşağıda açıklandığı gibi doldurmamız gereken yer tutuculardır.

**Olduğu gibi çalışan özellikler:**

- **`tag "$meta.id"`**: Daha kolay izleme için günlüklerdeki süreç adlarına örnek kimliği ekler
- **`label 'process_single'`**: CPU/bellek gereksinimlerini yapılandırmak için kaynak etiketi
- **`when:` bloğu**: `task.ext.when` yapılandırması aracılığıyla koşullu yürütmeye izin verir

Bu özellikler zaten işlevseldir ve modülleri daha sürdürülebilir hale getirir.

**Aşağıda özelleştireceğimiz yer tutucular:**

- **`input:` ve `output:` blokları**: Aracımıza uyacak şekilde güncelleyeceğimiz genel bildirimler
- **`script:` bloğu**: `cowpy` komutunu ekleyeceğimiz bir yorum içerir
- **`stub:` bloğu**: Doğru çıktıları üretecek şekilde güncelleyeceğimiz şablon
- **Container ve ortam**: Paket bilgileriyle dolduracağımız yer tutucular

Sonraki bölümler bu özelleştirmeleri tamamlamayı anlatır.

### 2.2. Container ve conda ortamını ayarlayın

nf-core yönergeleri, modülün bir parçası olarak hem bir container hem de bir Conda ortamı belirtmemizi gerektirir.

#### 2.2.1. Container

Container için, conda-forge paketleri dahil olmak üzere herhangi bir Conda paketinden otomatik olarak bir container oluşturmak için [Seqera Containers](https://seqera.io/containers/) kullanabilirsiniz.
Bu durumda daha önce olduğu gibi aynı önceden oluşturulmuş container'ı kullanıyoruz.

Varsayılan kod Docker ve Singularity arasında geçiş yapmayı teklif eder, ancak bu satırı basitleştireceğiz ve yukarıda Seqera Containers'dan aldığımız Docker container'ını belirteceğiz.

=== "Sonra"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273"
```

=== "Önce"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"
```

#### 2.2.2. Conda ortamı

Conda ortamı için, modül kodu `conda "${moduleDir}/environment.yml"` belirtir, bu da `environment.yml` dosyasında yapılandırılması gerektiği anlamına gelir.

Modül oluşturma aracı bizi `cowpy` paketini Bioconda'da (biyoinformatik araçları için birincil kanal) bulamadığı konusunda uyardı.
Ancak, `cowpy` conda-forge'da mevcuttur, bu nedenle `environment.yml` dosyasını şu şekilde tamamlayabilirsiniz:

=== "Sonra"

    ```yaml title="modules/local/cowpy/environment.yml"  linenums="1" hl_lines="1 3 5"
    name: COWPY
    channels:
      - conda-forge
    dependencies:
      - cowpy=1.1.5
    ```

=== "Önce"

    ```yaml title="modules/local/cowpy/environment.yml" linenums="1"
    ---
    # yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
    channels:
      - conda-forge
      - bioconda
    dependencies:
      # TODO nf-core: List required Conda package(s).
      #               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
      #               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
      - "YOUR-TOOL-HERE"
    ```

nf-core'a göndermek için varsayılanları daha yakından takip etmemiz gerekir, ancak kendi kullanımımız için kodu bu şekilde basitleştirebiliriz.

!!! tip "Bioconda vs conda-forge paketleri"

    - **Bioconda paketleri**: Otomatik olarak BioContainers oluşturulur, kullanıma hazır container'lar sağlar
    - **conda-forge paketleri**: Conda tarifinden talep üzerine container'lar oluşturmak için Seqera Containers kullanabilir

    Çoğu biyoinformatik araç Bioconda'dadır, ancak conda-forge araçları için Seqera Containers, konteynerleştirme için kolay bir çözüm sağlar.

### 2.3. `COWPY` mantığını ekleyin

Şimdi `COWPY` sürecinin ne yaptığına özgü kod öğelerini güncelleyelim: girdiler ve çıktılar ve script bloğu.

#### 2.3.1. Girdiler ve çıktılar

Oluşturulan şablon, belirli aracınız için özelleştirmeniz gereken genel girdi ve çıktı bildirimleri içerir.
Bölüm 1'deki manuel `COWPY` modülümüze bakarak, bunu bir kılavuz olarak kullanabiliriz.

Girdi ve çıktı bloklarını güncelleyin:

=== "Sonra"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: cowpy_output
    path "versions.yml"           , emit: versions
    ```

=== "Önce"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions
    ```

Bu şunları belirtir:

- Girdi dosyası parametre adı (genel `input` yerine `input_file`)
- Yapılandırılabilir önek desenini kullanan çıktı dosya adı (joker karakter `*` yerine `${prefix}.txt`)
- Açıklayıcı bir emit adı (genel `output` yerine `cowpy_output`)

Sözdizimini doğrulamak için Nextflow dil sunucusunu kullanıyorsanız, `${prefix}` kısmı bu aşamada bir hata olarak işaretlenecektir çünkü henüz script bloğuna eklemedik.
Şimdi buna geçelim.

#### 2.3.2. Script bloğu

Şablon, gerçek araç komutunu eklemeniz gereken script bloğunda bir yorum yer tutucusu sağlar.

Daha önce manuel olarak yazdığımız modüle dayanarak, aşağıdaki düzenlemeleri yapmalıyız:

=== "Sonra"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="3 6"
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat $input_file | cowpy $args > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Önce"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="6"
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Temel değişiklikler:

- `def prefix`'i sadece `prefix` olarak değiştirin (`def` olmadan) çıktı bloğunda erişilebilir hale getirmek için
- Yorumu hem `$args` hem de `${prefix}.txt` kullanan gerçek `cowpy` komutuyla değiştirin

`ext.args` ve `ext.prefix` yapılandırmasını `COWPY` süreci için `modules.config` dosyasına ekleme işini zaten yapmamış olsaydık, şimdi yapmamız gerekirdi.

#### 2.3.3. Stub bloğunu uygulama

Nextflow bağlamında, bir [stub](https://www.nextflow.io/docs/latest/process.html#stub) bloğu, gerçek komutu yürütmeden bir pipeline'ın mantığının hızlı prototiplenmesi ve test edilmesi için kullanılan hafif, sahte bir script tanımlamanıza olanak tanır.

<!-- TODO (gelecek) Bu çok yüzeysel ama gerçekten açıklanmalı veya en azından stub'lar hakkında bir açıklamaya bağlantı verilmeli (referans belgesi de pek yardımcı değil). Şu anda bu, stub'lar hakkında zaten bilmeyenler için büyük olasılıkla anlamsız olacaktır. -->

Bu gizemli görünüyorsa çok endişelenmeyin; bunu tamlık için dahil ediyoruz ancak tamamen isteğe bağlı olduğu için uğraşmak istemiyorsanız stub bölümünü de silebilirsiniz.

=== "Sonra"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Önce"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Temel değişiklikler:

- Script bloğuyla eşleşmesi için `def prefix`'i sadece `prefix` olarak değiştirin
- `echo $args` satırını kaldırın (bu sadece şablon yer tutucu koduydu)
- Stub, script bloğunun ürettiğiyle eşleşen boş bir `${prefix}.txt` dosyası oluşturur

Bu, gerçek aracın çalışmasını beklemeden iş akışı mantığını ve dosya işlemeyi test etmenizi sağlar.

Ortam kurulumunu (bölüm 2.2), girdileri/çıktıları (bölüm 2.3.1), script bloğunu (bölüm 2.3.2) ve stub bloğunu (bölüm 2.3.3) tamamladıktan sonra, modül test edilmeye hazırdır!

### 2.4. Yeni `COWPY` modülünü değiştirin ve pipeline'ı çalıştırın

Bu yeni `COWPY` modülü sürümünü denemek için tek yapmamız gereken, `hello.nf` iş akışı dosyasındaki import ifadesini yeni dosyaya işaret edecek şekilde değiştirmektir.

=== "Sonra"

    ```groovy title="workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Önce"

    ```groovy title="modules/local/cowpy/main.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Pipeline'ı test etmek için çalıştıralım.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Komut çıktısı"

    ```console hl_lines="33"
      N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [prickly_neumann] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_08-23-51

    Core Nextflow options
      runName                   : prickly_neumann
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [e9/008ede] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [f0/d70cfe] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [be/0ecc58] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [11/8e082f] CORE_HELLO:HELLO:COWPY (test)       [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Bu, daha önce olduğu gibi aynı sonuçları üretir.

### Özet

Artık her şeyi sıfırdan yazmak yerine şablonları kullanarak modülleri verimli bir şekilde oluşturmak için yerleşik nf-core araçlarını nasıl kullanacağınızı biliyorsunuz.

### Sırada ne var?

Modülleri nf-core'a katkıda bulunmanın faydalarının neler olduğunu ve ilgili ana adımların ve gereksinimlerin neler olduğunu öğrenin.

---

## 3. Modülleri nf-core'a geri katkıda bulunma

[nf-core/modules](https://github.com/nf-core/modules) deposu, iyi test edilmiş, standartlaştırılmış modüllerin katkılarını memnuniyetle karşılar.

### 3.1. Neden katkıda bulunmalı?

Modüllerinizi nf-core'a katkıda bulunmak:

- Araçlarınızı [nf-co.re/modules](https://nf-co.re/modules) adresindeki modüller kataloğu aracılığıyla tüm nf-core topluluğuna kullanılabilir hale getirir
- Devam eden topluluk bakımını ve iyileştirmelerini sağlar
- Kod incelemesi ve otomatik test yoluyla kalite güvencesi sağlar
- Çalışmanıza görünürlük ve tanınma kazandırır

### 3.2. Katkıda bulunanın kontrol listesi

nf-core'a bir modül katkıda bulunmak için aşağıdaki adımlardan geçmeniz gerekecektir:

1. [nf-co.re/modules](https://nf-co.re/modules) adresinde zaten var olup olmadığını kontrol edin
2. [nf-core/modules](https://github.com/nf-core/modules) deposunu fork edin
3. Şablonu oluşturmak için `nf-core modules create` kullanın
4. Modül mantığını ve testleri doldurun
5. `nf-core modules test tool/subtool` ile test edin
6. `nf-core modules lint tool/subtool` ile lint yapın
7. Bir pull request gönderin

Ayrıntılı talimatlar için [nf-core bileşenleri eğitimine](https://nf-co.re/docs/tutorials/nf-core_components/components) bakın.

### 3.3. Kaynaklar

- **Bileşenler eğitimi**: [Modül oluşturma ve katkıda bulunma için tam kılavuz](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Modül spesifikasyonları**: [Teknik gereksinimler ve yönergeler](https://nf-co.re/docs/guidelines/components/modules)
- **Topluluk desteği**: [nf-core Slack](https://nf-co.re/join) - `#modules` kanalına katılın

### Özet

Artık nf-core modüllerini nasıl oluşturacağınızı biliyorsunuz! Modülleri taşınabilir ve sürdürülebilir kılan dört temel deseni öğrendiniz:

- **Metadata demetleri** metadata'yı iş akışı boyunca yayar
- **`ext.args`** isteğe bağlı argümanları yapılandırma yoluyla ele alarak modül arayüzlerini basitleştirir
- **`ext.prefix`** çıktı dosyası adlandırmasını standartlaştırır
- **Merkezileştirilmiş yayınlama** modüllerde sabit kodlanmış yerine `modules.config` dosyasında yapılandırılan `publishDir` aracılığıyla

`COWPY`'yi adım adım dönüştürerek, bu desenlerin derin bir anlayışını geliştirdiniz, bu da sizi nf-core modülleriyle çalışmak, hata ayıklamak ve oluşturmak için donanımlı hale getirdi.
Pratikte, bu desenler başlangıçtan itibaren yerleşik olarak düzgün yapılandırılmış modüller oluşturmak için `nf-core modules create` kullanacaksınız.

Son olarak, modülleri nf-core topluluğuna nasıl katkıda bulunacağınızı öğrendiniz, araçları dünya çapındaki araştırmacılara kullanılabilir hale getirirken devam eden topluluk bakımından yararlanıyorsunuz.

### Sırada ne var?

Hazır olduğunuzda, pipeline'ınıza şema tabanlı girdi doğrulamasını nasıl ekleyeceğinizi öğrenmek için [Bölüm 5: Girdi doğrulaması](./05_input_validation.md) bölümüne devam edin.
