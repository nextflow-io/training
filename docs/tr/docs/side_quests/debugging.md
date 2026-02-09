# İş Akışlarında Hata Ayıklama

Hata ayıklama, saatlerce süren hayal kırıklığından sizi kurtarabilecek ve daha etkili bir Nextflow geliştiricisi olmanıza yardımcı olabilecek kritik bir beceridir. Kariyeriniz boyunca, özellikle yeni başladığınızda, iş akışlarınızı oluştururken ve sürdürürken hatalarla karşılaşacaksınız. Sistematik hata ayıklama yaklaşımlarını öğrenmek, sorunları hızlı bir şekilde belirlemenize ve çözmenize yardımcı olacaktır.

### Öğrenme hedefleri

Bu yan görevde, Nextflow iş akışları için **sistematik hata ayıklama tekniklerini** keşfedeceğiz:

- **Sözdizimi hatası ayıklama**: IDE özelliklerini ve Nextflow hata mesajlarını etkili bir şekilde kullanma
- **Kanal hata ayıklama**: Veri akışı sorunlarını ve kanal yapısı problemlerini teşhis etme
- **Süreç hata ayıklama**: Yürütme hatalarını ve kaynak sorunlarını araştırma
- **Yerleşik hata ayıklama araçları**: Nextflow'un önizleme modunu, stub çalıştırmasını ve çalışma dizinlerini kullanma
- **Sistematik yaklaşımlar**: Verimli hata ayıklama için dört aşamalı bir metodoloji

Sonunda, can sıkıcı hata mesajlarını çözümler için net yol haritalarına dönüştüren sağlam bir hata ayıklama metodolojisine sahip olacaksınız.

### Ön koşullar

Bu yan görevi üstlenmeden önce şunları yapmalısınız:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmak.
- Temel Nextflow kavramlarını ve mekanizmalarını (süreçler, kanallar, operatörler) rahatça kullanabiliyor olmak

**İsteğe bağlı:** Önce [Nextflow Geliştirme için IDE Özellikleri](./ide_features.md) yan görevini tamamlamanızı öneririz.
Bu, burada yoğun olarak kullanacağımız hata ayıklamayı destekleyen IDE özelliklerinin (sözdizimi vurgulama, hata algılama vb.) kapsamlı bir şekilde ele alınmasını içerir.

---

## 0. Başlayın

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine gidin

Bu eğitim için dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/debugging
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri gözden geçirin

Pratik için kullanacağımız çeşitli hata türlerine sahip örnek iş akışları bulacaksınız:

??? abstract "Dizin içeriği"

    ```console
    .
    ├── bad_bash_var.nf
    ├── bad_channel_shape.nf
    ├── bad_channel_shape_viewed_debug.nf
    ├── bad_channel_shape_viewed.nf
    ├── bad_number_inputs.nf
    ├── badpractice_syntax.nf
    ├── bad_resources.nf
    ├── bad_syntax.nf
    ├── buggy_workflow.nf
    ├── data
    │   ├── sample_001.fastq.gz
    │   ├── sample_002.fastq.gz
    │   ├── sample_003.fastq.gz
    │   ├── sample_004.fastq.gz
    │   ├── sample_005.fastq.gz
    │   └── sample_data.csv
    ├── exhausted.nf
    ├── invalid_process.nf
    ├── missing_output.nf
    ├── missing_software.nf
    ├── missing_software_with_stub.nf
    ├── nextflow.config
    └── no_such_var.nf
    ```

Bu dosyalar, gerçek dünya geliştirmede karşılaşacağınız yaygın hata ayıklama senaryolarını temsil eder.

#### Görevi gözden geçirin

Göreviniz her iş akışını çalıştırmak, hata(ları) belirlemek ve düzeltmektir.

Her hatalı iş akışı için:

1. **İş akışını çalıştırın** ve hatayı gözlemleyin
2. **Hata mesajını analiz edin**: Nextflow size ne söylüyor?
3. **Sorunu bulun** - sağlanan ipuçlarını kullanarak koddaki problemi tespit edin
4. **Hatayı düzeltin** ve çözümünüzün çalıştığını doğrulayın
5. Bir sonraki bölüme geçmeden önce **dosyayı sıfırlayın** (`git checkout <dosyaadı>` kullanın)

Alıştırmalar basit sözdizimi hatalarından daha ince çalışma zamanı sorunlarına doğru ilerler.
Çözümler satır içinde tartışılır, ancak ilerlemeden önce her birini kendiniz çözmeye çalışın.

#### Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Sözdizimi Hataları

Sözdizimi hataları, Nextflow kodu yazarken karşılaşacağınız en yaygın hata türüdür. Kod, Nextflow DSL'nin beklenen sözdizimi kurallarına uymadığında ortaya çıkarlar. Bu hatalar iş akışınızın hiç çalışmasını engellerler, bu nedenle bunları hızlı bir şekilde nasıl belirleyip düzelteceğinizi öğrenmek önemlidir.

### 1.1. Eksik parantezler

En yaygın sözdizimi hatalarından biri ve bazen hata ayıklaması en karmaşık olanlardan biri **eksik veya eşleşmeyen parantezlerdir**.

Pratik bir örnekle başlayalım.

#### Pipeline'ı çalıştırın

```bash
nextflow run bad_syntax.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [stupefied_bhabha] DSL2 - revision: ca6327fad2

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

**Sözdizimi hata mesajlarının temel öğeleri:**

- **Dosya ve konum**: Hatayı içeren dosyayı ve satır/sütunu gösterir (`bad_syntax.nf:24:1`)
- **Hata açıklaması**: Ayrıştırıcının beklemediği şeyi bulduğunu açıklar (`Unexpected input: '<EOF>'`)
- **EOF göstergesi**: `<EOF>` (Dosya Sonu) mesajı, ayrıştırıcının hala daha fazla içerik beklerken dosyanın sonuna ulaştığını gösterir - kapatılmamış parantezlerin klasik bir işareti

#### Kodu kontrol edin

Şimdi, hataya neyin sebep olduğunu anlamak için `bad_syntax.nf` dosyasını inceleyelim:

```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
// Missing closing brace for the process

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
```

Bu örneğin amacı doğrultusunda, hatanın nerede olduğunu göstermek için bir yorum bıraktık. Nextflow VSCode uzantısı da size neyin yanlış olabileceği konusunda bazı ipuçları vermeli, eşleşmeyen parantezi kırmızıyla işaretlemeli ve dosyanın erken bitişini vurgulamalıdır:

![Kötü sözdizimi](img/bad_syntax.png)

**Parantez hataları için hata ayıklama stratejisi:**

1. VS Code'un parantez eşleştirmesini kullanın (imleci bir parantezin yanına yerleştirin)
2. Parantezle ilgili mesajlar için Sorunlar panelini kontrol edin
3. Her açılış `{` parantezinin karşılık gelen bir kapanış `}` parantezine sahip olduğundan emin olun

#### Kodu düzeltin

Yorumu eksik kapanış parantezi ile değiştirin:

=== "Sonra"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }  // Add the missing closing brace

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

=== "Önce"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    // Missing closing brace for the process

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

#### Pipeline'ı çalıştırın

Şimdi çalıştığını doğrulamak için iş akışını tekrar çalıştırın:

```bash
nextflow run bad_syntax.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [insane_faggin] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [48/cd7f54] PROCESS_FILES (1) | 3 of 3 ✔
    ```

### 1.2. Yanlış süreç anahtar kelimeleri veya yönergeleri kullanma

Bir diğer yaygın sözdizimi hatası **geçersiz süreç tanımıdır**. Bu, gerekli blokları tanımlamayı unutursanız veya süreç tanımında yanlış yönergeler kullanırsanız olabilir.

#### Pipeline'ı çalıştırın

```bash
nextflow run invalid_process.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [nasty_jepsen] DSL2 - revision: da9758d614

    Error invalid_process.nf:3:1: Invalid process definition -- check for missing or out-of-order section labels
    │   3 | process PROCESS_FILES {
    │     | ^^^^^^^^^^^^^^^^^^^^^^^
    │   4 |     inputs:
    │   5 |     val sample_name
    │   6 |
    ╰   7 |     output:

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Kodu kontrol edin

Hata "Geçersiz süreç tanımı" belirtiyor ve sorunun etrafındaki bağlamı gösteriyor. 3-7. satırlara bakıldığında, 4. satırda `inputs:` görebiliriz, bu da sorundur. `invalid_process.nf` dosyasını inceleyelim:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // ERROR: Should be 'input' not 'inputs'
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
```

Hata bağlamındaki 4. satıra bakıldığında, sorunu fark edebiliriz: doğru `input` yönergesi yerine `inputs` kullanıyoruz. Nextflow VSCode uzantısı da bunu işaretleyecektir:

![Geçersiz süreç mesajı](img/invalid_process_message.png)

#### Kodu düzeltin

[Dokümantasyona](https://www.nextflow.io/docs/latest/process.html#) başvurarak yanlış anahtar kelimeyi doğrusuyla değiştirin:

=== "Sonra"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Fixed: Changed 'inputs' to 'input'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

=== "Önce"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // ERROR: Should be 'input' not 'inputs'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

#### Pipeline'ı çalıştırın

Şimdi çalıştığını doğrulamak için iş akışını tekrar çalıştırın:

```bash
nextflow run invalid_process.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [silly_fermi] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [b7/76cd9d] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.3. Hatalı değişken adları kullanma

Script bloklarınızda kullandığınız değişken adları geçerli olmalı, ya girdilerden ya da script'ten önce eklenen groovy kodundan türetilmelidir. Ancak pipeline geliştirmenin başında karmaşıklıkla uğraşırken, değişken adlandırmada hata yapmak kolaydır ve Nextflow bunu size hızlıca bildirecektir.

#### Pipeline'ı çalıştırın

```bash
nextflow run no_such_var.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [gloomy_meninsky] DSL2 - revision: 0c4d3bc28c

    Error no_such_var.nf:17:39: `undefined_var` is not defined
    │  17 |     echo "Using undefined variable: ${undefined_var}" >> ${output_pref
    ╰     |                                       ^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Hata derleme zamanında yakalanır ve 17. satırdaki tanımsız değişkeni doğrudan işaret eder, sorunun tam olarak nerede olduğunu gösteren bir işaretle.

#### Kodu kontrol edin

`no_such_var.nf` dosyasını inceleyelim:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

Hata mesajı, değişkenin script şablonunda tanınmadığını belirtiyor ve işte orada- script bloğunda kullanılan ancak başka bir yerde tanımlanmamış `${undefined_var}` değişkenini görmelisiniz.

#### Kodu düzeltin

'Böyle bir değişken yok' hatası alırsanız, değişkeni tanımlayarak (girdi değişken adlarını düzelterek veya script'ten önce groovy kodunu düzenleyerek) veya gerekli değilse script bloğundan kaldırarak düzeltebilirsiniz:

=== "Sonra"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Removed the line with undefined_var
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Önce"

    ```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Pipeline'ı çalıştırın

Şimdi çalıştığını doğrulamak için iş akışını tekrar çalıştırın:

```bash
nextflow run no_such_var.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [suspicious_venter] DSL2 - revision: 6ba490f7c5

    executor >  local (3)
    [21/237300] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.4. Bash değişkenlerinin hatalı kullanımı

Nextflow'da başlarken, Nextflow (Groovy) ve Bash değişkenleri arasındaki farkı anlamak zor olabilir. Bu, script bloğunun Bash içeriğinde değişkenleri kullanmaya çalışırken ortaya çıkan başka bir hatalı değişken hatası biçimi oluşturabilir.

#### Pipeline'ı çalıştırın

```bash
nextflow run bad_bash_var.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [infallible_mandelbrot] DSL2 - revision: 0853c11080

    Error bad_bash_var.nf:13:42: `prefix` is not defined
    │  13 |     echo "Processing ${sample_name}" > ${prefix}.txt
    ╰     |                                          ^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Kodu kontrol edin

Hata, `${prefix}` değişkeninin kullanıldığı 13. satırı işaret ediyor. Soruna neyin sebep olduğunu görmek için `bad_bash_var.nf` dosyasını inceleyelim:

```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    prefix="${sample_name}_output"
    echo "Processing ${sample_name}" > ${prefix}.txt  # HATA: ${prefix} Groovy sözdizimi, Bash değil
    """
}
```

Bu örnekte, `prefix` değişkenini Bash'te tanımlıyoruz, ancak bir Nextflow sürecinde ona başvurmak için kullandığımız `$` sözdizimi (`${prefix}`) Bash değil Groovy değişkeni olarak yorumlanır. Değişken Groovy bağlamında mevcut olmadığı için 'böyle bir değişken yok' hatası alırız.

#### Kodu düzeltin

Bir Bash değişkeni kullanmak istiyorsanız, dolar işaretini şu şekilde kaçırmalısınız:

=== "Sonra"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > \${prefix}.txt  # Düzeltildi: Dolar işareti kaçırıldı
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Önce"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > ${prefix}.txt  # HATA: ${prefix} Groovy sözdizimi, Bash değil
        """
    }
    ```

Bu, Nextflow'a bunu bir Bash değişkeni olarak yorumlamasını söyler.

#### Pipeline'ı çalıştırın

Şimdi çalıştığını doğrulamak için iş akışını tekrar çalıştırın:

```bash
nextflow run bad_bash_var.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [naughty_franklin] DSL2 - revision: 58c1c83709

    executor >  local (3)
    [4e/560285] PROCESS_FILES (2) | 3 of 3 ✔
    ```

!!! tip "Groovy vs Bash Değişkenleri"

    String birleştirme veya önek/sonek işlemleri gibi basit değişken manipülasyonları için, script bloğundaki Bash değişkenleri yerine script bölümündeki Groovy değişkenlerini kullanmak genellikle daha okunabilirdir:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    Bu yaklaşım, dolar işaretlerini kaçırma ihtiyacını ortadan kaldırır ve kodu okumayı ve sürdürmeyi kolaylaştırır.

### 1.5. Workflow Bloğu Dışındaki İfadeler

Nextflow VSCode uzantısı, hatalara neden olacak kod yapısıyla ilgili sorunları vurgular. Yaygın bir örnek, kanalları `workflow {}` bloğunun dışında tanımlamaktır - bu artık bir sözdizimi hatası olarak uygulanmaktadır.

#### Pipeline'ı çalıştırın

```bash
nextflow run badpractice_syntax.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [intergalactic_colden] DSL2 - revision: 5e4b291bde

    Error badpractice_syntax.nf:3:1: Statements cannot be mixed with script declarations -- move statements into a process or workflow
    │   3 | input_ch = channel.of('sample1', 'sample2', 'sample3')
    ╰     | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Hata mesajı sorunu açıkça belirtiyor: ifadeler (kanal tanımları gibi) bir workflow veya process bloğunun dışında script bildirimleriyle karıştırılamaz.

#### Kodu kontrol edin

Hataya neyin sebep olduğunu görmek için `badpractice_syntax.nf` dosyasını inceleyelim:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Channel defined outside workflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {
    PROCESS_FILES(input_ch)
}
```

VSCode uzantısı ayrıca `input_ch` değişkeninin workflow bloğunun dışında tanımlandığını vurgulayacaktır:

![Ölümcül olmayan sözdizimi hatası](img/nonlethal.png)

#### Kodu düzeltin

Kanal tanımını workflow bloğunun içine taşıyın:

=== "Sonra"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Moved inside workflow block
        PROCESS_FILES(input_ch)
    }
    ```

=== "Önce"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Channel defined outside workflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        PROCESS_FILES(input_ch)
    }
    ```

#### Pipeline'ı çalıştırın

Düzeltmenin çalıştığını doğrulamak için iş akışını tekrar çalıştırın:

```bash
nextflow run badpractice_syntax.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [naughty_ochoa] DSL2 - revision: 5e4b291bde

    executor >  local (3)
    [6a/84a608] PROCESS_FILES (2) | 3 of 3 ✔
    ```

Girdi kanallarınızı workflow bloğu içinde tanımlı tutun ve genel olarak uzantının yaptığı diğer önerileri takip edin.

### Özet

Nextflow hata mesajlarını ve IDE görsel göstergelerini kullanarak sözdizimi hatalarını sistematik olarak belirleyebilir ve düzeltebilirsiniz. Yaygın sözdizimi hataları arasında eksik parantezler, yanlış süreç anahtar kelimeleri, tanımsız değişkenler ve Bash ile Nextflow değişkenlerinin uygunsuz kullanımı yer alır. VSCode uzantısı, bunların çoğunu çalışma zamanından önce yakalamaya yardımcı olur. Araç setinizdeki bu sözdizimi hata ayıklama becerileriyle, en yaygın Nextflow sözdizimi hatalarını hızlı bir şekilde çözebilecek ve daha karmaşık çalışma zamanı sorunlarıyla uğraşmaya geçebileceksiniz.

### Sırada ne var?

Sözdizimi doğru olsa bile ortaya çıkan daha karmaşık kanal yapısı hatalarını ayıklamayı öğrenin.

---

## 2. Kanal Yapısı Hataları

Kanal yapısı hataları sözdizimi hatalarından daha incedir çünkü kod sözdizimsel olarak doğrudur, ancak veri şekilleri süreçlerin beklediğiyle eşleşmez. Nextflow pipeline'ı çalıştırmaya çalışacaktır, ancak girdi sayısının beklediğiyle eşleşmediğini bulabilir ve başarısız olabilir. Bu hatalar genellikle yalnızca çalışma zamanında görünür ve iş akışınızdan akan verileri anlamayı gerektirir.

!!! tip "`.view()` ile Kanalları Hata Ayıklama"

    Bu bölüm boyunca, iş akışınızın herhangi bir noktasında kanal içeriğini incelemek için `.view()` operatörünü kullanabileceğinizi unutmayın. Bu, kanal yapısı sorunlarını anlamak için en güçlü hata ayıklama araçlarından biridir. Bu tekniği bölüm 2.4'te ayrıntılı olarak inceleyeceğiz, ancak örnekler üzerinde çalışırken kullanmaktan çekinmeyin.

    ```groovy
    my_channel.view()  // Shows what's flowing through the channel
    ```

### 2.1. Yanlış Sayıda Girdi Kanalı

Bu hata, bir sürecin beklediğinden farklı sayıda kanal geçirdiğinizde oluşur.

#### Pipeline'ı çalıştırın

```bash
nextflow run bad_number_inputs.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [happy_swartz] DSL2 - revision: d83e58dcd3

    Error bad_number_inputs.nf:23:5: Incorrect number of call arguments, expected 1 but received 2
    │  23 |     PROCESS_FILES(samples_ch, files_ch)
    ╰     |     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Kodu kontrol edin

Hata mesajı, çağrının 1 argüman beklediğini ancak 2 aldığını açıkça belirtiyor ve 23. satırı işaret ediyor. `bad_number_inputs.nf` dosyasını inceleyelim:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Process expects only 1 input

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create two separate channels
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: Passing 2 channels but process expects only 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

Süreç yalnızca bir tane tanımlarken birden fazla girdi kanalı sağlayan eşleşmeyen `PROCESS_FILES` çağrısını görmelisiniz. VSCode uzantısı da süreç çağrısının altını kırmızıyla çizecek ve fare ile üzerine geldiğinizde bir tanı mesajı sağlayacaktır:

![Yanlış sayıda argüman mesajı](img/incorrect_num_args.png)

#### Kodu düzeltin

Bu özel örnek için, süreç tek bir kanal bekliyor ve ikinci kanala ihtiyaç duymuyor, bu nedenle yalnızca `samples_ch` kanalını geçirerek düzeltebiliriz:

=== "Sonra"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Fixed: Pass only the channel the process expects
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Önce"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERROR: Passing 2 channels but process expects only 1
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### Pipeline'ı çalıştırın

```bash
nextflow run bad_number_inputs.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [big_euler] DSL2 - revision: e302bd87be

    executor >  local (3)
    [48/497f7b] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Bu örnekten daha yaygın olarak, bir sürece ek girdiler ekleyebilir ve workflow çağrısını buna göre güncellemeyi unutabilirsiniz, bu da bu tür bir hataya yol açabilir. Neyse ki, hata mesajı uyumsuzluk hakkında oldukça net olduğu için bu, anlaşılması ve düzeltilmesi daha kolay hatalardan biridir.

### 2.2. Kanal Tükenmesi (Süreç Beklenenden Daha Az Çalışır)

Bazı kanal yapısı hataları çok daha incedir ve hiç hata üretmez. Muhtemelen bunların en yaygını, yeni Nextflow kullanıcılarının kuyruk kanallarının tükenebileceğini ve öğelerin biteceğini anlamadaki zorluğunu yansıtır, bu da iş akışının erken bitmesi anlamına gelir.

#### Pipeline'ı çalıştırın

```bash
nextflow run exhausted.nf
```

??? success "Komut çıktısı"

```console title="Exhausted channel output"
 N E X T F L O W   ~  version 25.10.2

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

Bu iş akışı hatasız tamamlanır, ancak yalnızca tek bir örneği işler!

#### Kodu kontrol edin

Bunun doğru olup olmadığını görmek için `exhausted.nf` dosyasını inceleyelim:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Define variables in Groovy code before the script
    output_prefix = "${reference}_${sample_name}"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {

    reference_ch = channel.of('baseline_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

Süreç üç kez yerine yalnızca bir kez çalışır çünkü `reference_ch` kanalı, ilk süreç yürütmesinden sonra tükenen bir kuyruk kanalıdır. Bir kanal tükendiğinde, diğer kanalların hala öğeleri olsa bile tüm süreç durur.

Bu, birden fazla örnek üzerinde yeniden kullanılması gereken tek bir referans dosyanız olduğu yaygın bir desendir. Çözüm, referans kanalını süresiz olarak yeniden kullanılabilecek bir değer kanalına dönüştürmektir.

#### Kodu düzeltin

Kaç dosyanın etkilendiğine bağlı olarak bunu ele almanın birkaç yolu vardır.

**Seçenek 1**: Çok fazla yeniden kullandığınız tek bir referans dosyanız var. Basitçe tekrar tekrar kullanılabilecek bir değer kanal türü oluşturabilirsiniz. Bunu yapmanın üç yolu vardır:

**1a** `channel.value()` kullanın:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Value channel can be reused
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** `first()` [operatörünü](https://www.nextflow.io/docs/latest/reference/operator.html#first) kullanın:

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** `collect()` [operatörünü](https://www.nextflow.io/docs/latest/reference/operator.html#collect) kullanın:

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Seçenek 2**: Daha karmaşık senaryolarda, belki de örnek kanalındaki tüm örnekler için birden fazla referans dosyanız olduğunda, iki kanalı demetler halinde birleştiren yeni bir kanal oluşturmak için `combine` operatörünü kullanabilirsiniz:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Creates cartesian product

    PROCESS_FILES(combined_ch)
}
```

`.combine()` operatörü iki kanalın kartezyen çarpımını oluşturur, böylece `reference_ch` içindeki her öğe `input_ch` içindeki her öğeyle eşleştirilir. Bu, sürecin referansı kullanırken her örnek için çalışmasına olanak tanır.

Bu, süreç girdisinin ayarlanmasını gerektirir. Örneğimizde, süreç tanımının başlangıcının aşağıdaki gibi ayarlanması gerekir:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

Bu yaklaşım tüm durumlarda uygun olmayabilir.

#### Pipeline'ı çalıştırın

Yukarıdaki düzeltmelerden birini deneyin ve iş akışını tekrar çalıştırın:

```bash
nextflow run exhausted.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `exhausted.nf` [maniac_leavitt] DSL2 - revision: f372a56a7d

    executor >  local (3)
    [80/0779e9] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Artık yalnızca bir yerine üç örneğin de işlendiğini görmelisiniz.

### 2.3. Yanlış Kanal İçerik Yapısı

İş akışları belirli bir karmaşıklık seviyesine ulaştığında, her kanalın iç yapılarını takip etmek biraz zor olabilir ve insanlar yaygın olarak sürecin beklediği ile kanalın gerçekte içerdiği arasında uyumsuzluklar oluştururlar. Bu, kanal sayısının yanlış olduğu daha önce tartıştığımız sorundan daha incedir. Bu durumda, doğru sayıda girdi kanalınız olabilir, ancak bu kanallardan bir veya daha fazlasının iç yapısı sürecin beklediğiyle eşleşmez.

#### Pipeline'ı çalıştırın

```bash
nextflow run bad_channel_shape.nf
```

??? failure "Komut çıktısı"

    ```console
    Launching `bad_channel_shape.nf` [hopeful_pare] DSL2 - revision: ffd66071a1

    executor >  local (3)
    executor >  local (3)
    [3f/c2dcb3] PROCESS_FILES (3) [  0%] 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `[sample1, file1.txt]_output.txt` expected by process `PROCESS_FILES (1)`


    Command executed:

      echo "Processing [sample1, file1.txt]" > [sample1, file1.txt]_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/d6/1fb69d1d93300bbc9d42f1875b981e

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Kodu kontrol edin

Hata mesajındaki köşeli parantezler burada ipucunu sağlar - süreç demeti tek bir değer olarak ele alıyor, bu istediğimiz şey değil. `bad_channel_shape.nf` dosyasını inceleyelim:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Expects single value, gets tuple

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

Demetlerden oluşan bir kanal oluşturduğumuzu görebilirsiniz: `['sample1', 'file1.txt']`, ancak süreç tek bir değer bekliyor, `val sample_name`. Yürütülen komut, sürecin `[sample3, file3.txt]_output.txt` adlı bir dosya oluşturmaya çalıştığını gösteriyor, bu amaçlanan çıktı değil.

#### Kodu düzeltin

Bunu düzeltmek için, süreç her iki girdiyi de gerektiriyorsa süreci bir demet kabul edecek şekilde ayarlayabiliriz:

=== "Seçenek 1: Süreçte demet kabul et"

    === "Sonra"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Fixed: Accept tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

    === "Önce"

        ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                val sample_name  // Expects single value, gets tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "Seçenek 2: İlk öğeyi çıkar"

    === "Sonra"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Fixed: Extract first element
        }
        ```

    === "Önce"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### Pipeline'ı çalıştırın

Çözümlerden birini seçin ve iş akışını yeniden çalıştırın:

```bash
nextflow run bad_channel_shape.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape.nf` [clever_thompson] DSL2 - revision: 8cbcae3746

    executor >  local (3)
    [bb/80a958] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 2.4. Kanal Hata Ayıklama Teknikleri

#### Kanal İncelemesi için `.view()` Kullanımı

Kanallar için en güçlü hata ayıklama aracı `.view()` operatörüdür. `.view()` ile, hata ayıklamaya yardımcı olmak için kanallarınızın şeklini tüm aşamalarda anlayabilirsiniz.

#### Pipeline'ı çalıştırın

Bunu çalışırken görmek için `bad_channel_shape_viewed.nf` dosyasını çalıştırın:

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [maniac_poisson] DSL2 - revision: b4f24dc9da

    executor >  local (3)
    [c0/db76b3] PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

#### Kodu kontrol edin

`.view()` öğesinin nasıl kullanıldığını görmek için `bad_channel_shape_viewed.nf` dosyasını inceleyelim:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debug: Show original channel content
    .map { tuple -> tuple[0] }        // Transform: Extract first element
    .view { "After mapping: $it" }    // Debug: Show transformed channel content

    PROCESS_FILES(input_ch)
}
```

#### Kodu düzeltin

Gelecekte kanal içeriğini anlamak için `.view()` işlemlerini aşırı kullanmaktan sizi kurtarmak için, yardımcı olacak bazı yorumlar eklemeniz önerilir:

```groovy title="bad_channel_shape_viewed.nf (with comments)" linenums="16" hl_lines="8 9"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

Bu, iş akışlarınız karmaşıklık açısından büyüdükçe ve kanal yapısı daha opak hale geldikçe daha önemli hale gelecektir.

#### Pipeline'ı çalıştırın

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [marvelous_koch] DSL2 - revision: 03e79cdbad

    executor >  local (3)
    [ff/d67cec] PROCESS_FILES (2) | 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

### Özet

Birçok kanal yapısı hatası geçerli Nextflow sözdizimi ile oluşturulabilir. Veri akışını anlayarak, inceleme için `.view()` operatörlerini kullanarak ve beklenmeyen demet yapılarını gösteren köşeli parantezler gibi hata mesajı desenlerini tanıyarak kanal yapısı hatalarını ayıklayabilirsiniz.

### Sırada ne var?

Süreç tanımlarının oluşturduğu hatalar hakkında bilgi edinin.

---

## 3. Süreç Yapısı Hataları

Süreçlerle ilgili karşılaşacağınız hataların çoğu, komutu oluştururken yaptığınız hatalarla veya temel yazılımla ilgili sorunlarla ilgili olacaktır. Bununla birlikte, yukarıdaki kanal sorunlarına benzer şekilde, sözdizimi hatası olarak nitelendirilmeyen ancak çalışma zamanında hatalara neden olacak süreç tanımında hatalar yapabilirsiniz.

### 3.1. Eksik Çıktı Dosyaları

Süreç yazarken yaygın bir hata, sürecin beklediği ile üretilenler arasında uyumsuzluk yaratan bir şey yapmaktır.

#### Pipeline'ı çalıştırın

```bash
nextflow run missing_output.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [zen_stone] DSL2 - revision: 37ff61f926

    executor >  local (3)
    executor >  local (3)
    [fd/2642e9] process > PROCESS_FILES (2) [ 66%] 2 of 3, failed: 2
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Missing output file(s) `sample3.txt` expected by process `PROCESS_FILES (3)`


    Command executed:

      echo "Processing sample3" > sample3_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/02/9604d49fb8200a74d737c72a6c98ed

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Kodu kontrol edin

Hata mesajı, sürecin `sample3.txt` adlı bir çıktı dosyası üretmesini beklediğini, ancak script'in aslında `sample3_output.txt` oluşturduğunu gösteriyor. `missing_output.nf` dosyasındaki süreç tanımını inceleyelim:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Expects: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
    """
}
```

`output:` bloğundaki çıktı dosya adı ile script'te kullanılan arasında bir uyumsuzluk olduğunu görmelisiniz. Bu uyumsuzluk sürecin başarısız olmasına neden olur. Bu tür bir hatayla karşılaşırsanız, geri dönün ve çıktıların süreç tanımınız ile çıktı bloğunuz arasında eşleştiğini kontrol edin.

Sorun hala net değilse, gerçek oluşturulan çıktı dosyalarını belirlemek için çalışma dizininin kendisini kontrol edin:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

Bu örnek için bu, `output:` tanımımızın aksine çıktı dosya adına bir `_output` sonekinin dahil edildiğini vurgulayacaktır.

#### Kodu düzeltin

Çıktı dosya adını tutarlı hale getirerek uyumsuzluğu düzeltin:

=== "Sonra"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Fixed: Match the script output

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }
    ```

=== "Önce"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}.txt"  // Expects: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
        """
    }
    ```

#### Pipeline'ı çalıştırın

```bash
nextflow run missing_output.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [elated_hamilton] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [16/1c437c] PROCESS_FILES (3) | 3 of 3 ✔
    ```

### 3.2. Eksik yazılım

Başka bir hata sınıfı, yazılım sağlamadaki hatalardan kaynaklanır. `missing_software.nf`, sözdizimsel olarak geçerli bir iş akışıdır, ancak kullandığı `cowpy` komutunu sağlamak için bazı harici yazılımlara bağlıdır.

#### Pipeline'ı çalıştırın

```bash
nextflow run missing_software.nf
```

??? failure "Komut çıktısı"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Süreç, belirttiğimiz komuta erişemiyor. Bazen bu, bir script'in iş akışı `bin` dizininde mevcut olması ancak çalıştırılabilir hale getirilmemiş olmasından kaynaklanır. Diğer zamanlarda, yazılımın iş akışının çalıştığı container veya ortamda yüklü olmamasından kaynaklanır.

#### Kodu kontrol edin

O `127` çıkış koduna dikkat edin - size tam olarak sorunu söyler. `missing_software.nf` dosyasını inceleyelim:

```groovy title="missing_software.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

#### Kodu düzeltin

Burada biraz samimiyetsiz davrandık ve aslında kodda yanlış bir şey yok. Sadece süreci söz konusu komuta erişimi olacak şekilde çalıştırmak için gerekli yapılandırmayı belirtmemiz gerekiyor. Bu durumda sürecin bir container tanımı var, bu yüzden tek yapmamız gereken iş akışını Docker etkinleştirilmiş olarak çalıştırmak.

#### Pipeline'ı çalıştırın

`nextflow.config` dosyasında sizin için bir Docker profili ayarladık, böylece iş akışını şununla çalıştırabilirsiniz:

```bash
nextflow run missing_software.nf -profile docker
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software.nf` [awesome_stonebraker] DSL2 - revision: 0296d12839

    executor >  local (3)
    [38/ab20d1] PROCESS_FILES (1) | 3 of 3 ✔
    ```

!!! note

    Nextflow'un container'ları nasıl kullandığı hakkında daha fazla bilgi edinmek için [Hello Nextflow](../hello_nextflow/05_hello_containers.md) bölümüne bakın

### 3.3. Hatalı kaynak yapılandırması

Üretim kullanımında, süreçlerinizde kaynakları yapılandırıyor olacaksınız. Örneğin `memory`, süreciniz için mevcut maksimum bellek miktarını tanımlar ve süreç bunu aşarsa, zamanlayıcınız genellikle süreci sonlandırır ve `137` çıkış kodu döndürür. Bunu burada gösteremeyiz çünkü `local` yürütücüsünü kullanıyoruz, ancak `time` ile benzer bir şey gösterebiliriz.

#### Pipeline'ı çalıştırın

`bad_resources.nf`, 1 milisaniyelik gerçekçi olmayan bir zaman sınırına sahip süreç yapılandırmasına sahiptir:

```bash
nextflow run bad_resources.nf -profile docker
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [disturbed_elion] DSL2 - revision: 27d2066e86

    executor >  local (3)
    [c0/ded8e1] PROCESS_FILES (3) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (2)'

    Caused by:
      Process exceeded running time limit (1ms)

    Command executed:

      cowpy sample2 > sample2_output.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/53/f0a4cc56d6b3dc2a6754ff326f1349

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Kodu kontrol edin

`bad_resources.nf` dosyasını inceleyelim:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // ERROR: Unrealistic time limit

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // Takes 1 second, but time limit is 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

Sürecin bir saniyeden uzun süreceğini biliyoruz (emin olmak için içine bir sleep ekledik), ancak süreç 1 milisaniye sonra zaman aşımına uğrayacak şekilde ayarlanmış. Birisi yapılandırmasıyla biraz gerçekçi olmamış!

#### Kodu düzeltin

Zaman sınırını gerçekçi bir değere yükseltin:

=== "Sonra"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // Fixed: Realistic time limit

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

=== "Önce"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '1 ms'  // ERROR: Unrealistic time limit

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // Takes 1 second, but time limit is 1ms
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### Pipeline'ı çalıştırın

```bash
nextflow run bad_resources.nf -profile docker
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [friendly_mcclintock] DSL2 - revision: 381567d2c1

    executor >  local (3)
    [c2/9b4c41] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Hata mesajlarınızı dikkatlice okuduğunuzdan emin olursanız, bunun gibi başarısızlıklar sizi çok uzun süre şaşırtmamalıdır. Ancak kaynak yönergelerinizi uygun şekilde yapılandırabilmeniz için çalıştırdığınız komutların kaynak gereksinimlerini anladığınızdan emin olun.

### 3.4. Süreç Hata Ayıklama Teknikleri

Süreçler başarısız olduğunda veya beklenmedik şekilde davrandığında, neyin yanlış gittiğini araştırmak için sistematik tekniklere ihtiyacınız vardır. Çalışma dizini, süreç yürütmesini hata ayıklamak için ihtiyacınız olan tüm bilgileri içerir.

#### Çalışma Dizini İncelemesini Kullanma

Süreçler için en güçlü hata ayıklama aracı, çalışma dizinini incelemektir. Bir süreç başarısız olduğunda, Nextflow o belirli süreç yürütmesi için ne olduğunu anlamak için gereken tüm dosyaları içeren bir çalışma dizini oluşturur.

#### Pipeline'ı çalıştırın

Çalışma dizini incelemesini göstermek için önceki `missing_output.nf` örneğini kullanalım (gerekirse bir çıktı adlandırma uyumsuzluğunu yeniden oluşturun):

```bash
nextflow run missing_output.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [irreverent_payne] DSL2 - revision: 3d5117f7e2

    executor >  local (3)
    [5d/d544a4] PROCESS_FILES (2) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `sample1.txt` expected by process `PROCESS_FILES (1)`

    Command executed:

      echo "Processing sample1" > sample1_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/1e/2011154d0b0f001cd383d7364b5244

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Çalışma dizinini kontrol edin

Bu hatayı aldığınızda, çalışma dizini tüm hata ayıklama bilgilerini içerir. Hata mesajından çalışma dizini yolunu bulun ve içeriğini inceleyin:

```bash
# Hata mesajından çalışma dizinini bulun
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

Ardından anahtar dosyaları inceleyebilirsiniz:

##### Komut Script'ini Kontrol Edin

`.command.sh` dosyası tam olarak hangi komutun yürütüldüğünü gösterir:

```bash
# Yürütülen komutu görüntüle
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

Bu şunları ortaya çıkarır:

- **Değişken ikamesi**: Nextflow değişkenlerinin düzgün şekilde genişletilip genişletilmediği
- **Dosya yolları**: Girdi dosyalarının doğru şekilde bulunup bulunmadığı
- **Komut yapısı**: Script sözdiziminin doğru olup olmadığı

Aranacak yaygın sorunlar:

- **Eksik tırnaklar**: Boşluk içeren değişkenlerin düzgün tırnak içine alınması gerekir
- **Yanlış dosya yolları**: Var olmayan veya yanlış konumlardaki girdi dosyaları
- **Yanlış değişken adları**: Değişken referanslarındaki yazım hataları
- **Eksik ortam kurulumu**: Belirli ortamlara bağlı komutlar

##### Hata Çıktısını Kontrol Edin

`.command.err` dosyası gerçek hata mesajlarını içerir:

```bash
# Hata çıktısını görüntüle
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

Bu dosya şunları gösterecektir:

- **Çıkış kodları**: 127 (komut bulunamadı), 137 (sonlandırıldı), vb.
- **İzin hataları**: Dosya erişim sorunları
- **Yazılım hataları**: Uygulamaya özgü hata mesajları
- **Kaynak hataları**: Bellek/zaman sınırı aşıldı

##### Standart Çıktıyı Kontrol Edin

`.command.out` dosyası komutunuzun ne ürettiğini gösterir:

```bash
# Standart çıktıyı görüntüle
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

Bu şunları doğrulamaya yardımcı olur:

- **Beklenen çıktı**: Komutun doğru sonuçları üretip üretmediği
- **Kısmi yürütme**: Komutun başlayıp yarıda başarısız olup olmadığı
- **Hata ayıklama bilgisi**: Script'inizden herhangi bir tanı çıktısı

##### Çıkış Kodunu Kontrol Edin

`.exitcode` dosyası süreç için çıkış kodunu içerir:

```bash
# Çıkış kodunu görüntüle
cat work/*/*/.exitcode
```

Yaygın çıkış kodları ve anlamları:

- **Çıkış kodu 127**: Komut bulunamadı - yazılım kurulumunu kontrol edin
- **Çıkış kodu 137**: Süreç sonlandırıldı - bellek/zaman sınırlarını kontrol edin

##### Dosya Varlığını Kontrol Edin

Süreçler eksik çıktı dosyaları nedeniyle başarısız olduğunda, gerçekte hangi dosyaların oluşturulduğunu kontrol edin:

```bash
# Çalışma dizinindeki tüm dosyaları listele
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

Bu şunları belirlemeye yardımcı olur:

- **Dosya adlandırma uyumsuzlukları**: Beklenenden farklı adlara sahip çıktı dosyaları
- **İzin sorunları**: Oluşturulamayan dosyalar
- **Yol sorunları**: Yanlış dizinlerde oluşturulan dosyalar

Önceki örneğimizde, bu beklenen `sample3.txt` dosyamızın mevcut olmadığını, ancak `sample3_output.txt` dosyasının olduğunu doğruladı:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### Özet

Süreç hata ayıklama, neyin yanlış gittiğini anlamak için çalışma dizinlerini incelemeyi gerektirir. Anahtar dosyalar arasında `.command.sh` (yürütülen script), `.command.err` (hata mesajları) ve `.command.out` (standart çıktı) bulunur. 127 (komut bulunamadı) ve 137 (süreç sonlandırıldı) gibi çıkış kodları, başarısızlık türü hakkında anında tanı ipuçları sağlar.

### Sırada ne var?

Nextflow'un yerleşik hata ayıklama araçları ve sorun giderme için sistematik yaklaşımlar hakkında bilgi edinin.

---

## 4. Yerleşik Hata Ayıklama Araçları ve İleri Teknikler

Nextflow, iş akışı yürütmesini hata ayıklamak ve analiz etmek için birkaç güçlü yerleşik araç sağlar. Bu araçlar, neyin yanlış gittiğini, nerede yanlış gittiğini ve bunu verimli bir şekilde nasıl düzelteceğinizi anlamanıza yardımcı olur.

### 4.1. Gerçek Zamanlı Süreç Çıktısı

Bazen çalışan süreçlerin içinde neler olduğunu görmeniz gerekir. Gerçek zamanlı süreç çıktısını etkinleştirebilirsiniz, bu size her görevin yürütülürken tam olarak ne yaptığını gösterir.

#### Pipeline'ı çalıştırın

Önceki örneklerimizden `bad_channel_shape_viewed.nf`, `.view()` kullanarak kanal içeriğini yazdırdı, ancak `bad_channel_shape_viewed_debug.nf` dosyasında gösterdiğimiz gibi, sürecin içinden değişkenleri yankılamak için `debug` yönergesini de kullanabiliriz. İş akışını çalıştırın:

```bash
nextflow run bad_channel_shape_viewed_debug.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed_debug.nf` [agitated_crick] DSL2 - revision: ea3676d9ec

    executor >  local (3)
    [c6/2dac51] process > PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    Sample name inside process is sample2

    Sample name inside process is sample1

    Sample name inside process is sample3
    ```

#### Kodu kontrol edin

`debug` yönergesinin nasıl çalıştığını görmek için `bad_channel_shape_viewed_debug.nf` dosyasını inceleyelim:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Enable real-time output

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Sample name inside process is ${sample_name}"
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}
```

`debug` yönergesi, bir sürecin ortamını anlamak için hızlı ve kullanışlı bir yol olabilir.

### 4.2. Önizleme Modu

Bazen herhangi bir süreç çalışmadan önce sorunları yakalamak istersiniz. Nextflow bu tür proaktif hata ayıklama için bir bayrak sağlar: `-preview`.

#### Pipeline'ı çalıştırın

Önizleme modu, komutları yürütmeden iş akışı mantığını test etmenizi sağlar. Bu, iş akışınızın yapısını hızlı bir şekilde kontrol etmek ve gerçek komutları çalıştırmadan süreçlerin doğru şekilde bağlandığından emin olmak için oldukça yararlı olabilir.

!!! note

    `bad_syntax.nf` dosyasını daha önce düzelttiyseniz, bu komutu çalıştırmadan önce script bloğundan sonra kapanış parantezini kaldırarak sözdizimi hatasını yeniden ekleyin.

Bu komutu çalıştırın:

```bash
nextflow run bad_syntax.nf -preview
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [magical_mercator] DSL2 - revision: 550b9a8873

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Önizleme modu, herhangi bir süreç çalıştırmadan sözdizimi hatalarını erken yakalamak için özellikle yararlıdır. Yürütmeden önce iş akışı yapısını ve süreç bağlantılarını doğrular.

### 4.3. Mantık Testi için Stub Çalıştırma

Bazen hatalar hata ayıklaması zordur çünkü komutlar çok uzun sürer, özel yazılım gerektirir veya karmaşık nedenlerle başarısız olur. Stub çalıştırma, gerçek komutları yürütmeden iş akışı mantığını test etmenizi sağlar.

#### Pipeline'ı çalıştırın

Bir Nextflow süreci geliştirirken, gerçek komutu çalıştırmadan doğru biçimde çıktılar üreten 'sahte' komutları tanımlamak için `stub` yönergesini kullanabilirsiniz. Bu yaklaşım, gerçek yazılımın karmaşıklıklarıyla uğraşmadan önce iş akışı mantığınızın doğru olduğunu doğrulamak istediğinizde özellikle değerlidir.

Örneğin, önceki `missing_software.nf` dosyamızı hatırlıyor musunuz? `-profile docker` ekleyene kadar iş akışının çalışmasını engelleyen eksik yazılımın olduğu? `missing_software_with_stub.nf` çok benzer bir iş akışıdır. Aynı şekilde çalıştırırsak, aynı hatayı üreteceğiz:

```bash
nextflow run missing_software_with_stub.nf
```

??? failure "Komut çıktısı"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Ancak, bu iş akışı `-stub-run` ile çalıştırırsak, `docker` profili olmadan bile hata üretmeyecektir:

```bash
nextflow run missing_software_with_stub.nf -stub-run
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software_with_stub.nf` [astonishing_shockley] DSL2 - revision: f1f4f05d7d

    executor >  local (3)
    [b5/2517a3] PROCESS_FILES (3) | 3 of 3 ✔
    ```

#### Kodu kontrol edin

`missing_software_with_stub.nf` dosyasını inceleyelim:

```groovy title="missing_software.nf (with stub)" hl_lines="16-19" linenums="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """

    stub:
    """
    touch ${sample_name}_output.txt
    """
}
```

`missing_software.nf` dosyasına göre, bu süreç, Nextflow stub modunda çalıştırıldığında `script:` içinde belirtilen yerine kullanılacak bir komut belirten bir `stub:` yönergesine sahiptir.

Burada kullandığımız `touch` komutu herhangi bir yazılıma veya uygun girdilere bağlı değildir ve tüm durumlarda çalışacaktır, bu da süreç iç yapıları hakkında endişelenmeden iş akışı mantığını hata ayıklamamıza olanak tanır.

**Stub çalıştırma şunları hata ayıklamaya yardımcı olur:**

- Kanal yapısı ve veri akışı
- Süreç bağlantıları ve bağımlılıkları
- Parametre yayılımı
- Yazılım bağımlılıkları olmadan iş akışı mantığı

### 4.4. Sistematik Hata Ayıklama Yaklaşımı

Artık bireysel hata ayıklama tekniklerini öğrendiğinize göre - trace dosyalarından ve çalışma dizinlerinden önizleme moduna, stub çalıştırmaya ve kaynak izlemeye kadar - bunları sistematik bir metodolojide bir araya getirelim. Yapılandırılmış bir yaklaşıma sahip olmak, karmaşık hatalar tarafından bunalmaktan kaçınmanızı ve önemli ipuçlarını kaçırmamanızı sağlar.

Bu metodoloji, ele aldığımız tüm araçları verimli bir iş akışında birleştirir:

**Dört Aşamalı Hata Ayıklama Yöntemi:**

**Aşama 1: Sözdizimi Hatası Çözümü (5 dakika)**

1. VSCode veya IDE'nizde kırmızı alt çizgileri kontrol edin
2. Sözdizimi sorunlarını belirlemek için `nextflow run workflow.nf -preview` komutunu çalıştırın
3. Tüm sözdizimi hatalarını düzeltin (eksik parantezler, sondaki virgüller, vb.)
4. Devam etmeden önce iş akışının başarıyla ayrıştırıldığından emin olun

**Aşama 2: Hızlı Değerlendirme (5 dakika)**

1. Çalışma zamanı hata mesajlarını dikkatlice okuyun
2. Bunun bir çalışma zamanı, mantık veya kaynak hatası olup olmadığını kontrol edin
3. Temel iş akışı mantığını test etmek için önizleme modunu kullanın

**Aşama 3: Ayrıntılı Araştırma (15-30 dakika)**

1. Başarısız görevin çalışma dizinini bulun
2. Log dosyalarını inceleyin
3. Kanalları incelemek için `.view()` operatörleri ekleyin
4. Yürütme olmadan iş akışı mantığını test etmek için `-stub-run` kullanın

**Aşama 4: Düzeltme ve Doğrulama (15 dakika)**

1. Minimal hedefli düzeltmeler yapın
2. Resume ile test edin: `nextflow run workflow.nf -resume`
3. Tam iş akışı yürütmesini doğrulayın

!!! tip "Verimli Hata Ayıklama için Resume Kullanımı"

    Bir sorunu belirlediğinizde, başarılı iş akışı bölümlerini yeniden çalıştırarak zaman kaybetmeden düzeltmelerinizi test etmek için verimli bir yola ihtiyacınız vardır. Nextflow'un `-resume` işlevselliği hata ayıklama için paha biçilmezdir.

    [Hello Nextflow](../hello_nextflow/) üzerinde çalıştıysanız `-resume` ile karşılaşmış olacaksınız ve sorun sürecinizden önceki süreçler çalışırken beklerken kendinizi kurtarmak için hata ayıklarken bundan iyi kullanım yapmanız önemlidir.

    **Resume hata ayıklama stratejisi:**

    1. Başarısızlığa kadar iş akışını çalıştırın
    2. Başarısız görev için çalışma dizinini inceleyin
    3. Belirli sorunu düzeltin
    4. Yalnızca düzeltmeyi test etmek için resume edin
    5. İş akışı tamamlanana kadar tekrarlayın

#### Hata Ayıklama Yapılandırma Profili

Bu sistematik yaklaşımı daha da verimli hale getirmek için, ihtiyacınız olan tüm araçları otomatik olarak etkinleştiren özel bir hata ayıklama yapılandırması oluşturabilirsiniz:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Conservative resources for debugging
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

Ardından pipeline'ı bu profil etkinleştirilmiş olarak çalıştırabilirsiniz:

```bash
nextflow run workflow.nf -profile debug
```

Bu profil gerçek zamanlı çıktıyı etkinleştirir, çalışma dizinlerini korur ve daha kolay hata ayıklama için paralelleştirmeyi sınırlar.

### 4.5. Pratik Hata Ayıklama Alıştırması

Şimdi sistematik hata ayıklama yaklaşımını pratiğe dökme zamanı. `buggy_workflow.nf` iş akışı, gerçek dünya geliştirmede karşılaşacağınız sorun türlerini temsil eden birkaç yaygın hata içerir.

!!! exercise

    `buggy_workflow.nf` dosyasındaki tüm hataları belirlemek ve düzeltmek için sistematik hata ayıklama yaklaşımını kullanın. Bu iş akışı bir CSV dosyasından örnek verileri işlemeye çalışır ancak yaygın hata ayıklama senaryolarını temsil eden birden fazla kasıtlı hata içerir.

    İlk hatayı görmek için iş akışını çalıştırarak başlayın:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "Komut çıktısı"

        ```console
        N E X T F L O W   ~  version 25.10.2

        Launching `buggy_workflow.nf` [wise_ramanujan] DSL2 - revision: d51a8e83fd

        ERROR ~ Range [11, 12) out of bounds for length 11

         -- Check '.nextflow.log' file for details
        ```

        Bu şifreli hata, `params{}` bloğunda 11-12. satırlar civarında bir ayrıştırma sorunu olduğunu gösterir. v2 ayrıştırıcısı yapısal sorunları erken yakalar.

    Öğrendiğiniz dört aşamalı hata ayıklama yöntemini uygulayın:

    **Aşama 1: Sözdizimi Hatası Çözümü**
    - VSCode veya IDE'nizde kırmızı alt çizgileri kontrol edin
    - Sözdizimi sorunlarını belirlemek için `nextflow run workflow.nf -preview` komutunu çalıştırın
    - Tüm sözdizimi hatalarını düzeltin (eksik parantezler, sondaki virgüller, vb.)
    - Devam etmeden önce iş akışının başarıyla ayrıştırıldığından emin olun

    **Aşama 2: Hızlı Değerlendirme**
    - Çalışma zamanı hata mesajlarını dikkatlice okuyun
    - Hataların çalışma zamanı, mantık veya kaynakla ilgili olup olmadığını belirleyin
    - Temel iş akışı mantığını test etmek için `-preview` modunu kullanın

    **Aşama 3: Ayrıntılı Araştırma**
    - Başarısız görevler için çalışma dizinlerini inceleyin
    - Kanalları incelemek için `.view()` operatörleri ekleyin
    - Çalışma dizinlerindeki log dosyalarını kontrol edin
    - Yürütme olmadan iş akışı mantığını test etmek için `-stub-run` kullanın

    **Aşama 4: Düzeltme ve Doğrulama**
    - Hedefli düzeltmeler yapın
    - Düzeltmeleri verimli bir şekilde test etmek için `-resume` kullanın
    - Tam iş akışı yürütmesini doğrulayın

    **Emrinizde Hata Ayıklama Araçları:**
    ```bash
    # Sözdizimi kontrolü için önizleme modu
    nextflow run buggy_workflow.nf -preview

    # Ayrıntılı çıktı için debug profili
    nextflow run buggy_workflow.nf -profile debug

    # Mantık testi için stub çalıştırma
    nextflow run buggy_workflow.nf -stub-run

    # Düzeltmelerden sonra resume
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution
        `buggy_workflow.nf`, tüm ana hata ayıklama kategorilerini kapsayan 9 veya 10 farklı hata içerir (nasıl saydığınıza bağlı olarak). İşte her hatanın sistematik bir dökümü ve nasıl düzeltileceği

        O sözdizimi hatalarıyla başlayalım:

        **Hata 1: Sözdizimi Hatası - Sondaki Virgül**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // ERROR: Trailing comma
        ```
        **Düzeltme:** Sondaki virgülü kaldırın
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **Hata 2: Sözdizimi Hatası - Eksik Kapanış Parantezi**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // ERROR: Missing closing brace for processFiles process
        ```
        **Düzeltme:** Eksik kapanış parantezini ekleyin
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Add missing closing brace
        ```

        **Hata 3: Değişken Adı Hatası**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // ERROR: should be sample_id
        cat ${input_file} > ${sample}_result.txt  // ERROR: should be sample_id
        ```
        **Düzeltme:** Doğru girdi değişken adını kullanın
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Hata 4: Tanımsız Değişken Hatası**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // ERROR: sample_ids undefined
        ```
        **Düzeltme:** Doğru kanalı kullanın ve örnek ID'lerini çıkarın
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        Bu noktada iş akışı çalışacaktır, ancak hala hatalar alacağız (örneğin `processFiles` içinde `Path value cannot be null`), bunlar kötü kanal yapısından kaynaklanır.

        **Hata 5: Kanal Yapısı Hatası - Yanlış Map Çıktısı**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // ERROR: processFiles expects tuple
        ```
        **Düzeltme:** processFiles'ın beklediği demet yapısını döndürün
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        Ancak bu, yukarıdaki `heavyProcess()` çalıştırma düzeltmemizi bozacaktır, bu nedenle yalnızca örnek ID'lerini o sürece geçirmek için bir map kullanmamız gerekecektir:

        **Hata 6: heavyProcess için kötü kanal yapısı**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // ERROR: input_ch now has 2 elements per emission- heavyProcess only needs 1 (the first)
        ```
        **Düzeltme:** Doğru kanalı kullanın ve örnek ID'lerini çıkarın
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Şimdi biraz daha ilerliyoruz ancak `No such variable: i` hakkında bir hata alıyoruz, çünkü bir Bash değişkenini kaçırmadık.

        **Hata 7: Bash Değişkeni Kaçırma Hatası**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // ERROR: $i not escaped
        ```
        **Düzeltme:** Bash değişkenini kaçırın
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Şimdi `Process exceeded running time limit (1ms)` alıyoruz, bu nedenle ilgili süreç için çalışma zamanı sınırını düzeltiriz:

        **Hata 8: Kaynak Yapılandırma Hatası**
        ```groovy linenums="36"
        time '1 ms'  // ERROR: Unrealistic time limit
        ```
        **Düzeltme:** Gerçekçi bir zaman sınırına yükseltin
        ```groovy linenums="36"
        time '100 s'
        ```

        Sonra çözülecek bir `Missing output file(s)` hatamız var:

        **Hata 9: Çıktı Dosya Adı Uyumsuzluğu**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // ERROR: Wrong filename, should match output declaration
        ```
        **Düzeltme:** Çıktı bildirimiyle eşleştirin
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        İlk iki süreç çalıştı, ancak üçüncüsü çalışmadı.

        **Hata 10: Çıktı Dosya Adı Uyumsuzluğu**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Error: attempting to take input from the pwd rather than a process
        handleFiles(file_ch)
        ```
        **Düzeltme:** Önceki süreçten çıktıyı alın
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        Bununla, tüm iş akışı çalışmalıdır.

        **Tam Düzeltilmiş İş Akışı:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Hata ayıklama alıştırmaları için hatalı iş akışı
        * Bu iş akışı öğrenme amaçlı birkaç kasıtlı hata içerir
        */

        params{
            // Parameters with missing validation
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Girdi/çıktı uyumsuzluğu olan süreç
        */
        process processFiles {
            publishDir "${params.output}/processed", mode: 'copy'

            input:
                tuple val(sample_id), path(input_file)

            output:
                path "${sample_id}_result.txt"

            script:
            """
            echo "Processing: ${sample_id}"
            cat ${input_file} > ${sample_id}_result.txt
            """
        }

        /*
        * Kaynak sorunları olan süreç
        */
        process heavyProcess {
            publishDir "${params.output}/heavy", mode: 'copy'

            time '100 s'

            input:
                val sample_id

            output:
                path "${sample_id}_heavy.txt"

            script:
            """
            # Ağır hesaplamayı simüle et
            for i in {1..1000000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * Dosya işleme sorunları olan süreç
        */
        process handleFiles {
            publishDir "${params.output}/files", mode: 'copy'

            input:
                path input_file

            output:
                path "processed_${input_file}"

            script:
            """
            if [ -f "${input_file}" ]; then
                cp ${input_file} processed_${input_file}
            fi
            """
        }

        /*
        * Kanal sorunları olan ana iş akışı
        */
        workflow {

            // Channel with incorrect usage
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map{it[0]})

            handleFiles(heavyProcess.out)
        }
        ```

**Kapsanan Hata Kategorileri:**

- **Sözdizimi hataları**: Eksik parantezler, sondaki virgüller, tanımsız değişkenler
- **Kanal yapısı hataları**: Yanlış veri şekilleri, tanımsız kanallar
- **Süreç hataları**: Çıktı dosya uyumsuzlukları, değişken kaçırma
- **Kaynak hataları**: Gerçekçi olmayan zaman sınırları

**Temel Hata Ayıklama Dersleri:**

1. **Hata mesajlarını dikkatlice okuyun** - genellikle doğrudan soruna işaret ederler
2. **Sistematik yaklaşımlar kullanın** - bir seferde bir hatayı düzeltin ve `-resume` ile test edin
3. **Veri akışını anlayın** - kanal yapısı hataları genellikle en incedir
4. **Çalışma dizinlerini kontrol edin** - süreçler başarısız olduğunda, loglar size tam olarak neyin yanlış gittiğini söyler

---

## Özet

Bu yan görevde, Nextflow iş akışlarını hata ayıklamak için bir dizi sistematik teknik öğrendiniz.
Bu teknikleri kendi çalışmanızda uygulamak, bilgisayarınızla savaşarak daha az zaman harcamanızı, sorunları daha hızlı çözmenizi ve gelecekteki sorunlardan kendinizi korumanızı sağlayacaktır.

### Temel desenler

**1. Sözdizimi hatalarını nasıl belirleyip düzeltirsiniz**:

- Nextflow hata mesajlarını yorumlama ve sorunları bulma
- Yaygın sözdizimi hataları: eksik parantezler, yanlış anahtar kelimeler, tanımsız değişkenler
- Nextflow (Groovy) ve Bash değişkenleri arasında ayrım yapma
- Erken hata algılama için VS Code uzantı özelliklerini kullanma

```groovy
// Missing brace - look for red underlines in IDE
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- missing!

// Wrong keyword
inputs:  // Should be 'input:'

// Undefined variable - escape with backslash for Bash variables
echo "${undefined_var}"      // Nextflow variable (error if not defined)
echo "\${bash_var}"          // Bash variable (escaped)
```

**2. Kanal yapısı sorunlarını nasıl hata ayıklarsınız**:

- Kanal kardinalitesi ve tükenme sorunlarını anlama
- Kanal içerik yapısı uyumsuzluklarını hata ayıklama
- Kanal incelemesi için `.view()` operatörlerini kullanma
- Çıktıda köşeli parantezler gibi hata mesajı desenlerini tanıma

```groovy
// Inspect channel content
my_channel.view { "Content: $it" }

// Convert queue to value channel (prevents exhaustion)
reference_ch = channel.value('ref.fa')
// or
reference_ch = channel.of('ref.fa').first()
```

**3. Süreç yürütme sorunlarını nasıl giderirsiniz**:

- Eksik çıktı dosya hatalarını teşhis etme
- Çıkış kodlarını anlama (eksik yazılım için 127, bellek sorunları için 137)
- Çalışma dizinlerini ve komut dosyalarını araştırma
- Kaynakları uygun şekilde yapılandırma

```bash
# Gerçekte ne yürütüldüğünü kontrol et
cat work/ab/cdef12/.command.sh

# Hata çıktısını kontrol et
cat work/ab/cdef12/.command.err

# Çıkış kodu 127 = komut bulunamadı
# Çıkış kodu 137 = sonlandırıldı (bellek/zaman sınırı)
```

**4. Nextflow'un yerleşik hata ayıklama araçlarını nasıl kullanırsınız**:

- Önizleme modundan ve gerçek zamanlı hata ayıklamadan yararlanma
- Mantık testi için stub çalıştırma uygulama
- Verimli hata ayıklama döngüleri için resume uygulama
- Dört aşamalı sistematik hata ayıklama metodolojisini takip etme

!!! tip "Hızlı Hata Ayıklama Referansı"

    **Sözdizimi hataları?** → VSCode uyarılarını kontrol edin, `nextflow run workflow.nf -preview` çalıştırın

    **Kanal sorunları?** → İçeriği incelemek için `.view()` kullanın: `my_channel.view()`

    **Süreç başarısızlıkları?** → Çalışma dizini dosyalarını kontrol edin:

    - `.command.sh` - yürütülen script
    - `.command.err` - hata mesajları
    - `.exitcode` - çıkış durumu (127 = komut bulunamadı, 137 = sonlandırıldı)

    **Gizemli davranış?** → İş akışı mantığını test etmek için `-stub-run` ile çalıştırın

    **Düzeltmeler yaptınız mı?** → Test ederken zaman kazanmak için `-resume` kullanın: `nextflow run workflow.nf -resume`

---

### Ek kaynaklar

- [Nextflow sorun giderme kılavuzu](https://www.nextflow.io/docs/latest/troubleshooting.html): Resmi sorun giderme dokümantasyonu
- [Nextflow kanallarını anlama](https://www.nextflow.io/docs/latest/channel.html): Kanal türleri ve davranışına derinlemesine bakış
- [Süreç yönergeleri referansı](https://www.nextflow.io/docs/latest/process.html#directives): Mevcut tüm süreç yapılandırma seçenekleri
- [nf-test](https://www.nf-test.com/): Nextflow pipeline'ları için test çerçevesi
- [Nextflow Slack topluluğu](https://www.nextflow.io/slack-invite.html): Topluluktan yardım alın

Üretim iş akışları için şunları düşünün:

- Ölçekte izleme ve hata ayıklama için [Seqera Platform](https://seqera.io/platform/) kurulumu
- Tekrarlanabilir yazılım ortamları için [Wave container'ları](https://seqera.io/wave/) kullanımı

**Unutmayın:** Etkili hata ayıklama, pratikle gelişen bir beceridir. Burada edindiğiniz sistematik metodoloji ve kapsamlı araç seti, Nextflow geliştirme yolculuğunuz boyunca size iyi hizmet edecektir.

---

## Sırada ne var?

[Yan Görevler menüsüne](./index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
