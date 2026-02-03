# Bölüm 4: Bir nf-core modülü oluşturma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core eğitim kursunun bu dördüncü bölümünde, modülleri taşınabilir ve sürdürülebilir kılan temel kuralları uygulayarak bir nf-core modülünün nasıl oluşturulacağını göstereceğiz.

nf-core projesi, Bölüm 2'de iş akışı için kullandığımıza benzer şekilde, düzgün yapılandırılmış modül şablonlarını otomatik olarak oluşturan bir komut (`nf-core modules create`) sağlar.
Ancak öğretim amaçları için, manuel olarak başlayacağız: `core-hello` iş hattınızdaki yerel `cowpy` modülünü adım adım nf-core tarzı bir modüle dönüştüreceğiz.
Bundan sonra, gelecekte daha verimli çalışmak için şablon tabanlı modül oluşturma yöntemini size göstereceğiz.

??? info "Bu bölüme nasıl başlanır"

    Bu bölüm, [Bölüm 3: Bir nf-core modülü kullanma](./03_use_module.md) kısmını tamamladığınızı ve `CAT_CAT` modülünü iş hattınıza entegre ettiğinizi varsayar.

    Bölüm 3'ü tamamlamadıysanız veya bu bölüm için yeni başlamak istiyorsanız, başlangıç noktanız olarak `core-hello-part3` çözümünü kullanabilirsiniz.
    Bu komutları `hello-nf-core/` dizini içinden çalıştırın:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    Bu size `CAT_CAT` modülünün zaten entegre edildiği bir iş hattı verir.
    Başarıyla çalıştığını test etmek için aşağıdaki komutu çalıştırabilirsiniz:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. `cowpy`'yi bir nf-core modülüne dönüştürme

Bu bölümde, `core-hello` iş hattınızdaki yerel `cowpy` modülüne nf-core kurallarını uygulayarak onu nf-core topluluk standartlarını izleyen bir modüle dönüştüreceğiz.

Mevcut `cowpy` process modülünün kodu şu şekildedir:

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

Aşağıdaki nf-core kurallarını adım adım uygulayacağız:

1. **Process adını `COWPY` olarak büyük harfe çevir** - Kurala uymak için.
2. **`COWPY`'yi metadata demetlerini kullanacak şekilde güncelle** - Örnek metadatasını iş akışı boyunca yaymak için.
3. **Araç argüman yapılandırmasını `ext.args` ile merkezileştir** - Arayüzü minimal tutarken modül çok yönlülüğünü artırmak için.
4. **Çıktı adlandırmasını `ext.prefix` ile standartlaştır** - Tutarlılığı teşvik etmek için.
5. **Yayınlama yapılandırmasını merkezileştir** - Tutarlılığı teşvik etmek için.

Her adımdan sonra, her şeyin beklendiği gibi çalıştığını test etmek için iş hattını çalıştıracağız.

!!! warning "Çalışma dizini"

    Bu bölümdeki tüm dosya düzenlemeleri ve komut çalıştırmaları için `core-hello` dizininde (iş hattı kök dizininizde) olduğunuzdan emin olun.

    ```bash
    cd core-hello
    ```

### 1.1. Process adını büyük harfe çevirme

Bu tamamen stilistik bir kuraldır (teknik bir gerekçesi yoktur) ancak nf-core modülleri için norm olduğu için uyalım.

Üç dizi değişiklik yapmamız gerekiyor:

1. Modüldeki process adını güncelle
2. İş akışı başlığındaki modül import ifadesini güncelle
3. İş akışı gövdesindeki process çağrısını ve emit bildirimini güncelle

Başlayalım!

#### 1.1.1. Modüldeki process adını güncelleme

`cowpy.nf` modül dosyasını açın (`core-hello/modules/local/` altında) ve process adını büyük harfe çevirin:

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

Eğer process adı birkaç kelimeden oluşsaydı, örneğin orijinalde camel case ile yazılmış MyCowpyTool gibi bir process'imiz olsaydı, nf-core kuralı onları ayırmak için alt çizgi kullanmak olurdu, bu da MY_COWPY_TOOL sonucunu verir.

#### 1.1.2. Modül import ifadesini güncelleme

Process adları büyük/küçük harfe duyarlıdır, bu nedenle process adını değiştirdiğimize göre, modül import ifadesini `hello.nf` iş akışı başlığında buna göre güncellememiz gerekiyor:

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

Process çağrılarını güncellemek zorunda kalmamak için import ifadesinde bir alias kullanabilirdik, ancak bu büyük harfe çevirme kuralını benimsemenin amacını biraz bozardı.

#### 1.1.3. Process çağrısını ve emit bildirimini güncelleme

Şimdi `hello.nf` workflow bloğundaki process'e yapılan iki referansı güncelleyelim:

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // cowpy ile selamlamaların ASCII sanatını oluştur
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
    // cowpy ile selamlamaların ASCII sanatını oluştur
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

**Her iki** değişikliği de yaptığınızdan emin olun, aksi takdirde bunu çalıştırdığınızda hata alırsınız.

#### 1.1.4. Test etmek için iş hattını çalıştırma

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

### 1.2. `COWPY`'yi metadata demetlerini kullanacak şekilde güncelleme

`core-hello` iş hattının mevcut versiyonunda, aşağıdaki diyagramın üst yarısında gösterildiği gibi, `CAT_CAT`'in çıktı demetinden dosyayı çıkarıp `COWPY`'ye aktarıyoruz.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

`COWPY`'nin metadata demetlerini doğrudan kabul etmesi daha iyi olurdu, bu da metadatanın diyagramın alt yarısında gösterildiği gibi iş akışı boyunca akmasına izin verir.

Bu amaçla, aşağıdaki değişiklikleri yapmamız gerekecek:

1. Girdi ve çıktı tanımlarını güncelle
2. İş akışındaki process çağrısını güncelle
3. İş akışındaki emit bloğunu güncelle

Tüm bunları yaptıktan sonra, her şeyin daha önce olduğu gibi çalıştığını test etmek için iş hattını çalıştıracağız.

#### 1.2.1. Girdi ve çıktı tanımlarını güncelleme

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

Gördüğünüz gibi, hem **ana girdide** hem de **çıktıda** Bölüm 3'te tanıtılan `tuple val(meta), path(input_file)` desenini izleyen bir demet kullanacak şekilde değiştirdik.
Çıktı için, bu fırsattan yararlanarak çıktı kanalına açıklayıcı bir ad vermek için `emit: cowpy_output` ekledik.

Artık process'in ne beklediğini değiştirdiğimize göre, process çağrısında ona sağladığımız şeyi güncellememiz gerekiyor.

#### 1.2.2. İş akışındaki process çağrısını güncelleme

İyi haber şu ki bu değişiklik process çağrısını basitleştirecek.
Artık `CAT_CAT` çıktısı ve `COWPY` girdisi aynı 'şekle' sahip, yani her ikisi de bir `tuple val(meta), path(input_file)` yapısından oluşuyor, bu nedenle `CAT_CAT` process'inin çıktısından dosyayı açıkça çıkarmak yerine onları doğrudan bağlayabiliriz.

`hello.nf` iş akışı dosyasını açın (`core-hello/workflows/` altında) ve `COWPY` çağrısını aşağıda gösterildiği gibi güncelleyin.

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // cowpy ile selamlamaların ASCII sanatını oluştur
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // cowpy ile selamlamaların ASCII sanatını oluştur
        COWPY(ch_for_cowpy, params.character)
    ```

Artık `COWPY`'yi doğrudan `CAT_CAT.out.file_out` üzerinde çağırıyoruz.

Sonuç olarak, artık `ch_for_cowpy` kanalını oluşturmamıza gerek yok, bu nedenle o satır (ve yorum satırı) tamamen silinebilir.

#### 1.2.3. İş akışındaki emit bloğunu güncelleme

`COWPY` artık adlandırılmış bir çıktı, `cowpy_output`, yayınladığından, `hello.nf` iş akışının `emit:` bloğunu bunu kullanacak şekilde güncelleyebiliriz.

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

Bu teknik olarak gerekli değil, ancak mümkün olduğunda adlandırılmış çıktılara başvurmak iyi bir uygulamadır.

#### 1.2.4. Test etmek için iş hattını çalıştırma

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

İş hattı başarıyla çalışmalı, metadata artık `CAT_CAT`'ten `COWPY`'ye akıyor.

`COWPY`'yi metadata demetlerini işleyecek şekilde yapmak için yapmamız gerekenleri tamamladık.
Şimdi, nf-core modül desenlerinden yararlanmak için başka neler yapabileceğimize bakalım.

### 1.3. Araç argüman yapılandırmasını `ext.args` ile merkezileştirme

Mevcut durumunda, `COWPY` process'i `character` parametresi için bir değer almayı bekler.
Sonuç olarak, process'i her çağırdığımızda, araç tarafından belirlenen varsayılanlarla mutlu olsak bile bir değer sağlamamız gerekir.
`COWPY` için bu kabul edilebilir bir sorun değil, ancak birçok opsiyonel parametresi olan araçlar için oldukça zahmetli olabilir.

nf-core projesi, araç argümanlarını yapılandırma dosyaları aracılığıyla daha rahat yönetmek için [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) adlı bir Nextflow özelliği kullanmayı önerir.

Her araç seçeneği için process girdileri bildirmek yerine, modülü komut satırının oluşturulmasında `ext.args`'a başvuracak şekilde yazarsınız.
Sonra sadece tüm modüller için yapılandırma detaylarını merkezileştiren `modules.config` dosyasında kullanmak istediğiniz argümanları ve değerleri tutacak şekilde `ext.args` değişkenini ayarlamanız yeterlidir.
Nextflow bu argümanları değerleriyle birlikte çalışma zamanında araç komut satırına ekleyecektir.

Bu yaklaşımı `COWPY` modülüne uygulayalım.
Aşağıdaki değişiklikleri yapmamız gerekecek:

1. `COWPY` modülünü güncelle
2. `modules.config` dosyasında `ext.args`'ı yapılandır
3. `hello.nf` iş akışını güncelle

Tüm bunları yaptıktan sonra, her şeyin daha önce olduğu gibi çalıştığını test etmek için iş hattını çalıştıracağız.

#### 1.3.1. `COWPY` modülünü güncelleme

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
   Bundan sonra, bu argümanı aşağıda açıklandığı gibi `ext.args` yapılandırması aracılığıyla sağlayacağız.

2. **`script:` bloğunda, `def args = task.ext.args ?: ''` satırını ekledik.**
   Bu satır, `args` değişkeninin değerini belirlemek için `?:` operatörünü kullanır: boş değilse `task.ext.args` içeriği, boşsa boş bir string.
   Genel olarak `ext.args`'a atıfta bulunsak da, bu kodun modül düzeyindeki `ext.args` yapılandırmasını çekmek için `task.ext.args`'a başvurması gerektiğini unutmayın.

3. **Komut satırında, `-c "$character"` yerine `$args` kullandık.**
   Nextflow'un `modules.config` dosyasındaki `ext.args` içinde ayarlanan tüm araç argümanlarını enjekte edeceği yer burasıdır.

Sonuç olarak, modül arayüzü artık daha basit: sadece temel metadata ve dosya girdilerini bekliyor.

!!! note

    `?:` operatörü genellikle 'Elvis operatörü' olarak adlandırılır çünkü yanal olarak Elvis Presley yüzüne benzer, `?` karakteri saçındaki dalgayı simgeler.

#### 1.3.2. `modules.config` dosyasında `ext.args`'ı yapılandırma

Artık `character` bildirimini modülden çıkardığımıza göre, onu `modules.config` yapılandırma dosyasındaki `ext.args`'a eklememiz gerekiyor.

Özellikle, `process {}` bloğuna şu kod parçasını ekleyeceğiz:

```groovy title="Eklenecek kod"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

`withName:` sözdizimi bu yapılandırmayı yalnızca `COWPY` process'ine atar ve `ext.args = { "-c ${params.character}" }` basitçe `character` parametresinin değerini içerecek bir string oluşturur.
Nextflow'a parametre değerini çalışma zamanında değerlendirmesini söyleyen küme parantezlerinin kullanımına dikkat edin.

Anlamlı mı? Hadi ekleyelim.

`conf/modules.config` dosyasını açın ve yapılandırma kodunu aşağıda gösterildiği gibi `process {}` bloğu içine ekleyin.

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

Umarım bir iş hattındaki tüm modüllerin `ext.args` değerlerinin bu dosyada belirtildiğini hayal edebilirsiniz, bu da aşağıdaki faydaları sağlar:

- **Modül arayüzü basit kalır** - Sadece temel metadata ve dosya girdilerini kabul eder
- **İş hattı hala `params.character` gösterir** - Son kullanıcılar daha önce olduğu gibi yapılandırabilir
- **Modül artık taşınabilir** - Belirli bir parametre adı beklemeden diğer iş hatlarında yeniden kullanılabilir
- **Yapılandırma merkezileştirilmiş** - `modules.config` içinde, iş akışı mantığını temiz tutar

`modules.config` dosyasını tüm iş hatlarının modül başına yapılandırmayı merkezileştirdiği yer olarak kullanarak, modüllerimizi farklı iş hatlarında daha yeniden kullanılabilir hale getiriyoruz.

#### 1.3.3. `hello.nf` iş akışını güncelleme

`COWPY` modülü artık girdi olarak `character` parametresini gerektirmediğinden, iş akışı çağrısını buna göre güncellememiz gerekiyor.

`hello.nf` iş akışı dosyasını açın (`core-hello/workflows/` altında) ve `COWPY` çağrısını aşağıda gösterildiği gibi güncelleyin.

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // cowpy ile selamlamaların ASCII sanatını oluştur
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // cowpy ile selamlamaların ASCII sanatını oluştur
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

İş akışı kodu artık daha temiz: `params.character`'ı doğrudan process'e geçirmemize gerek yok.
Modül arayüzü minimal tutularak daha taşınabilir hale getirilirken, iş hattı hala yapılandırma aracılığıyla açık seçeneği sağlar.

#### 1.3.4. Test etmek için iş hattını çalıştırma

İş akışının hala beklendiği gibi çalıştığını test edelim, `ext.args` yapılandırmasının çalıştığını doğrulamak için farklı bir karakter belirterek.

Daha... gizemli seçeneklerden biri olan `kosh` kullanarak bu komutu çalıştırın:

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
Çıktıyı dosya tarayıcıda bulun veya çıktı dosyasına bakmak için görev hash'ini (yukarıdaki örnekte `38/eb29ea` kısmı) kullanın:

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

`kosh` karakteriyle gösterilen ASCII sanatını görmelisiniz, `ext.args` yapılandırmasının çalıştığını doğrulayarak!

??? info "(Opsiyonel) Komut dosyasını inceleyin"

    Yapılandırmanın tam olarak nasıl uygulandığını görmek istiyorsanız, `.command.sh` dosyasını inceleyebilirsiniz:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    `-c kosh` argümanıyla `cowpy` komutunu göreceksiniz:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    Bu, `.command.sh` dosyasının `ext.args` yapılandırmasına göre doğru şekilde oluşturulduğunu gösterir.

Burada ne başardığımızı düşünmek için bir dakikanızı ayırın.
Bu yaklaşım, modül arayüzünü temel verilere (dosyalar, metadata ve zorunlu örnek başına parametreler) odaklı tutarken, aracın davranışını kontrol eden seçenekler ayrı olarak yapılandırma yoluyla ele alınır.

Bu, `cowpy` gibi basit bir araç için gereksiz görünebilir, ancak çok sayıda opsiyonel argümanı olan veri analiz araçları için büyük fark yaratabilir.

Bu yaklaşımın faydalarını özetlemek gerekirse:

- **Temiz arayüz**: Modül temel veri girdilerine (metadata ve dosyalar) odaklanır
- **Esneklik**: Kullanıcılar, örneğe özel değerler dahil olmak üzere, yapılandırma yoluyla araç argümanlarını belirtebilir
- **Tutarlılık**: Tüm nf-core modülleri bu deseni takip eder
- **Taşınabilirlik**: Modüller sabit kodlanmış araç seçenekleri olmadan yeniden kullanılabilir
- **İş akışı değişikliği yok**: Araç seçeneklerini eklemek veya değiştirmek iş akışı kodunu güncellemeyi gerektirmez

!!! note

    `ext.args` sistemi, burada ele alınmayan, metadata'ya göre argüman değerlerini dinamik olarak değiştirme dahil olmak üzere güçlü ek yeteneklere sahiptir. Daha fazla ayrıntı için [nf-core modül spesifikasyonlarına](https://nf-co.re/docs/guidelines/components/modules) bakın.

### 1.4. Çıktı adlandırmasını `ext.prefix` ile standartlaştırma

Artık `COWPY` process'ine metamap'e erişim verdiğimize göre, başka bir yararlı nf-core deseninden yararlanmaya başlayabiliriz: çıktı dosyalarını metadata'ya göre adlandırma.

Burada `ext.prefix` adlı bir Nextflow özelliğini kullanacağız, bu da `meta.id` kullanarak (metamap'e dahil edilen tanımlayıcı) modüller arasında çıktı dosya adlandırmasını standartlaştırmamıza ve yine de istenirse modülleri ayrı ayrı yapılandırabilmemize olanak tanıyacak.

Bu, `ext.args` ile yaptığımıza benzer olacak, ilerledikçe detaylandıracağımız birkaç farkla birlikte.

Bu yaklaşımı `COWPY` modülüne uygulayalım.
Aşağıdaki değişiklikleri yapmamız gerekecek:

1. `COWPY` modülünü güncelle
2. `modules.config` dosyasında `ext.prefix` yapılandır

(İş akışında değişiklik gerekmez.)

Bunları yaptıktan sonra, her şeyin daha önce olduğu gibi çalıştığını test etmek için iş hattını çalıştıracağız.

#### 1.4.1. `COWPY` modülünü güncelleme

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
   Bu satır, `prefix` değişkeninin değerini belirlemek için `?:` operatörünü kullanır: boş değilse `task.ext.prefix` içeriği, boşsa metamap'ten (`meta.id`) tanımlayıcı.
   Genel olarak `ext.prefix`'e atıfta bulunsak da, bu kodun modül düzeyindeki `ext.prefix` yapılandırmasını çekmek için `task.ext.prefix`'e başvurması gerektiğini unutmayın.

2. **Komut satırında, `cowpy-${input_file}` yerine `${prefix}.txt` kullandık.**
   Nextflow'un yukarıdaki satır tarafından belirlenen `prefix` değerini enjekte edeceği yer burasıdır.

3. **`output:` bloğunda, `path("cowpy-${input_file}")` yerine `path("${prefix}.txt")` kullandık.**
   Bu sadece komut satırında yazılana göre dosya yolunun ne olacağını yeniden ifade eder.

Sonuç olarak, çıktı dosya adı artık mantıklı bir varsayılan (metamap'ten tanımlayıcı) ve uygun dosya formatı uzantısı kullanılarak oluşturuluyor.

#### 1.4.2. `modules.config` dosyasında `ext.prefix` yapılandırma

Bu durumda mantıklı varsayılan bizim zevkimiz için yeterince ifade edici değil; daha önce olduğu gibi `cowpy-<id>.txt` araç adını içeren özel bir adlandırma deseni kullanmak istiyoruz.

Bunu, `ext.args` ile `character` parametresi için yaptığımız gibi `modules.config` dosyasında `ext.prefix` yapılandırarak yapacağız, bu sefer `withName: 'COWPY' {}` bloğu zaten mevcut ve sadece aşağıdaki satırı eklememiz gerekiyor:

```groovy title="Eklenecek kod"
ext.prefix = { "cowpy-${meta.id}" }
```

Bu istediğimiz string'i oluşturacak.
Bu sefer Nextflow'a `meta.id` değerini çalışma zamanında değerlendirmesini söylemek için bir kez daha küme parantezleri kullandığımıza dikkat edin.

Hadi ekleyelim.

`conf/modules.config` dosyasını açın ve yapılandırma kodunu aşağıda gösterildiği gibi `process {}` bloğu içine ekleyin.

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

Merak ediyorsanız, `ext.prefix` closure'ı doğru metadata parçasına erişime sahiptir çünkü yapılandırma, metadata'nın mevcut olduğu process yürütmesi bağlamında değerlendirilir.

#### 1.4.3. Test etmek için iş hattını çalıştırma

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

Sonuçlar dizinindeki çıktıya bir göz atın.
Varsayılan batch adına göre, cowpy çıktı dosyasını daha önce olduğu gibi aynı adlandırmayla görmelisiniz: `cowpy-test.txt`.

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

Modül veya iş akışı kodunda herhangi bir değişiklik yapmak zorunda kalmadan adlandırma desenini değiştirebileceğinizi kendinize kanıtlamak için `conf/modules.config` dosyasındaki `ext.prefix` yapılandırmasını değiştirmekten çekinmeyin.

Alternatif olarak, bu kısmın hala anında özelleştirilebilir olduğunu kendinize kanıtlamak için komut satırında farklı bir `--batch` parametresi belirterek bunu tekrar çalıştırmayı da deneyebilirsiniz.

Bu, `ext.prefix`'in modül arayüzünü esnek tutarken tercih ettiğiniz adlandırma kuralını korumanıza nasıl izin verdiğini gösterir.

Bu yaklaşımın faydalarını özetlemek gerekirse:

- **Standartlaştırılmış adlandırma**: Çıktı dosyaları genellikle metadata'dan örnek ID'leri kullanılarak adlandırılır
- **Yapılandırılabilir**: Kullanıcılar gerekirse varsayılan adlandırmayı geçersiz kılabilir
- **Tutarlı**: Tüm nf-core modülleri bu deseni takip eder
- **Öngörülebilir**: Çıktı dosyalarının ne olarak adlandırılacağını bilmek kolay

Oldukça iyi, değil mi?
Peki, modülümüzü nf-core yönergelerine uyacak şekilde iyileştirmek için yapmamız gereken bir önemli değişiklik daha var.

### 1.5. Yayınlama yapılandırmasını merkezileştirme

Çıktıları iki farklı dizine yayınladığımızı fark etmiş olabilirsiniz:

- **`results`** — Başlangıçtan beri yerel modüllerimiz için kullandığımız, modül başına `publishDir` yönergeleri kullanılarak ayrı ayrı ayarlanan orijinal çıktı dizini;
- **`core-hello-results`** — Komut satırında `--outdir` ile ayarlanan, nf-core loglarını ve `CAT_CAT` tarafından yayınlanan sonuçları alan çıktı dizini.

Bu dağınık ve optimal değil; her şey için tek bir konuma sahip olmak daha iyi olurdu.
Tabii ki, yerel modüllerimizin her birine gidip `publishDir` yönergesini `core-hello-results` dizinini kullanacak şekilde manuel olarak güncelleyebiliriz, ancak bir dahaki sefere çıktı dizinini değiştirmeye karar verirsek ne olur?

Bireysel modüllerin yayınlama kararları vermesi açıkça gidilecek yol değil, özellikle aynı modülün birçok farklı iş hattında, farklı ihtiyaçları veya tercihleri olan insanlar tarafından kullanılabileceği bir dünyada.
Çıktıların nereye yayınlanacağını iş akışı yapılandırma düzeyinde kontrol edebilmek istiyoruz.

"Hey," diyebilirsiniz, "`CAT_CAT` çıktılarını `--outdir`'e gönderiyor. Belki onun `publishDir` yönergesini kopyalamalıyız?"

Evet, bu harika bir fikir.

Ama onun bir `publishDir` yönergesi yok. (Devam edin, modül koduna bakın.)

Bunun nedeni, nf-core iş hatlarının, bireysel modüllerde değil `conf/modules.config` dosyasında `publishDir` yapılandırarak iş akışı düzeyinde kontrolü merkezileştirmesidir.
Özellikle, nf-core şablonu, geçersiz kılıcı bir yönerge sağlanmadıkça tüm modüllere uygulanan varsayılan bir `publishDir` yönergesi (önceden tanımlanmış bir dizin yapısıyla) bildirir.

Bu harika gelmiyor mu? Bu varsayılan yönergeden yararlanmak için yapmamız gereken tek şey yerel modüllerimizden mevcut `publishDir` yönergesini kaldırmak olabilir mi?

Ne olacağını görmek için bunu `COWPY` üzerinde deneyelim, sonra nasıl çalıştığını anlamak için varsayılan yapılandırma koduna bakacağız.

Son olarak, istenirse varsayılan davranışın nasıl geçersiz kılınacağını göstereceğiz.

#### 1.5.1. `COWPY`'den `publishDir` yönergesini kaldırma

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

#### 1.5.2. Test etmek için iş hattını çalıştırma

Şimdi iş hattını çalıştırırsak ne olacağına bakalım.

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

Mevcut çalışma dizininize bir göz atın.
Şimdi `core-hello-results` ayrıca `COWPY` modülünün çıktılarını da içeriyor.

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

Nextflow'un iş akışı ve modül adlarına göre bu dizin hiyerarşisini oluşturduğunu görebilirsiniz.

Sorumlu kod `conf/modules.config` dosyasında bulunur.
Bu, nf-core şablonunun bir parçası olan ve tüm process'lere uygulanan varsayılan `publishDir` yapılandırmasıdır:

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

- **`path:`** Process adına göre çıktı dizinini belirler.
  `task.process` içinde yer alan bir process'in tam adı, iş akışı ve modül import hiyerarşisini içerir (`CORE_HELLO:HELLO:CAT_CAT` gibi).
  `tokenize` işlemleri bu hiyerarşiyi soyarak sadece process adını alır, sonra herhangi bir alt çizgiden önceki ilk kısmı alır (varsa) ve küçük harfe dönüştürür.
  Bu, `CAT_CAT
