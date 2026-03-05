# İş Akışlarında Hata Ayıklama

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hata ayıklama, saatlerce süren hayal kırıklığından sizi kurtarabilecek ve sizi daha etkili bir Nextflow geliştiricisi haline getirebilecek kritik bir beceridir. Kariyeriniz boyunca, özellikle yeni başlarken, iş akışlarınızı oluştururken ve sürdürürken hatalarla karşılaşacaksınız. Sistematik hata ayıklama yaklaşımlarını öğrenmek, sorunları hızla belirlemenize ve çözmenize yardımcı olacaktır.

### Öğrenme hedefleri

Bu yan görevde, Nextflow iş akışları için **sistematik hata ayıklama tekniklerini** keşfedeceğiz:

- **Sözdizimi hatası ayıklama**: IDE özelliklerini ve Nextflow hata mesajlarını etkili şekilde kullanma
- **Kanal hata ayıklama**: Veri akışı sorunlarını ve kanal yapısı problemlerini teşhis etme
- **Process hata ayıklama**: Çalıştırma hatalarını ve kaynak sorunlarını araştırma
- **Yerleşik hata ayıklama araçları**: Nextflow'un önizleme modu, stub çalıştırma ve çalışma dizinlerinden yararlanma
- **Sistematik yaklaşımlar**: Verimli hata ayıklama için dört aşamalı bir metodoloji

Sonunda, sinir bozucu hata mesajlarını çözümler için net yol haritalarına dönüştüren sağlam bir hata ayıklama metodolojisine sahip olacaksınız.

### Ön koşullar

Bu yan görevi üstlenmeden önce şunları yapmalısınız:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olun.
- Temel Nextflow kavramlarını ve mekanizmalarını (process'ler, kanal'lar, operatör'ler) rahatça kullanabilin

**İsteğe bağlı:** [Nextflow Geliştirme için IDE Özellikleri](./ide_features.md) yan görevini önce tamamlamanızı öneririz.
Bu, burada yoğun olarak kullanacağımız hata ayıklamayı destekleyen IDE özelliklerinin (sözdizimi vurgulama, hata algılama vb.) kapsamlı bir şekilde ele alınmasını sağlar.

---

## 0. Başlangıç

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine gidin

Bu eğitim için dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/debugging
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri inceleyin

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

#### Ödevi inceleyin

Göreviniz her bir iş akışını çalıştırmak, hataları belirlemek ve düzeltmektir.

Her hatalı iş akışı için:

1. **İş akışını çalıştırın** ve hatayı gözlemleyin
2. **Hata mesajını analiz edin**: Nextflow size ne söylüyor?
3. **Sorunu belirleyin** sağlanan ipuçlarını kullanarak kodda
4. **Hatayı düzeltin** ve çözümünüzün çalıştığını doğrulayın
5. Bir sonraki bölüme geçmeden önce **dosyayı sıfırlayın** (`git checkout <dosyaadı>` kullanın)

Alıştırmalar basit sözdizimi hatalarından daha ince çalışma zamanı sorunlarına doğru ilerler.
Çözümler satır içinde tartışılır, ancak ilerlemeden önce her birini kendiniz çözmeye çalışın.

#### Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım
- [ ] Ödevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Sözdizimi Hataları

Sözdizimi hataları, Nextflow kodu yazarken karşılaşacağınız en yaygın hata türüdür. Kod, Nextflow DSL'nin beklenen sözdizimi kurallarına uymadığında meydana gelirler. Bu hatalar iş akışınızın hiç çalışmasını engellerler, bu nedenle bunları nasıl hızlı bir şekilde belirleyip düzelteceğinizi öğrenmek önemlidir.

### 1.1. Eksik parantezler

En yaygın sözdizimi hatalarından biri ve bazen hata ayıklaması en karmaşık olanlardan biri **eksik veya uyumsuz parantezlerdir**.

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

**Sözdizimi hata mesajlarının ana öğeleri:**

- **Dosya ve konum**: Hatayı içeren dosya ve satır/sütunu gösterir (`bad_syntax.nf:24:1`)
- **Hata açıklaması**: Ayrıştırıcının beklemediği şeyi bulduğunu açıklar (`Unexpected input: '<EOF>'`)
- **EOF göstergesi**: `<EOF>` (End Of File) mesajı, ayrıştırıcının hala daha fazla içerik beklerken dosyanın sonuna ulaştığını gösterir - kapatılmamış parantezlerin klasik bir işareti

#### Kodu kontrol edin

Şimdi, hatanın nedenini anlamak için `bad_syntax.nf` dosyasını inceleyelim:

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

    // Girdi kanalı oluştur
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // İşlemi girdi kanalıyla çağır
    PROCESS_FILES(input_ch)
}
```

Bu örneğin amacı için, hatanın nerede olduğunu göstermek için size bir yorum bıraktık. Nextflow VSCode uzantısı da size neyin yanlış olabileceği konusunda bazı ipuçları vermelidir, uyumsuz parantezi kırmızıya boyar ve dosyanın erken bitişini vurgular:

![Kötü sözdizimi](img/bad_syntax.png)

**Parantez hataları için hata ayıklama stratejisi:**

1. VS Code'un parantez eşleştirmesini kullanın (imleci bir parantezin yanına yerleştirin)
2. Parantezle ilgili mesajlar için Sorunlar panelini kontrol edin
3. Her açılan `{` için karşılık gelen kapatan `}` olduğundan emin olun

#### Kodu düzeltin

Yorumu eksik kapatan parantez ile değiştirin:

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

        // Girdi kanalı oluştur
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // İşlemi girdi kanalıyla çağır
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

        // Girdi kanalı oluştur
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // İşlemi girdi kanalıyla çağır
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

### 1.2. Yanlış process anahtar kelimelerini veya yönergelerini kullanma

Bir diğer yaygın sözdizimi hatası **geçersiz bir process tanımıdır**. Bu, gerekli blokları tanımlamayı unutursanız veya process tanımında yanlış yönergeler kullanırsanız meydana gelebilir.

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

Hata "Geçersiz process tanımı" belirtiyor ve sorunun etrafındaki bağlamı gösteriyor. 3-7. satırlara bakıldığında, 4. satırda `inputs:` görebiliriz, bu sorun. `invalid_process.nf` dosyasını inceleyelim:

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

    // Girdi kanalı oluştur
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // İşlemi girdi kanalıyla çağır
    PROCESS_FILES(input_ch)
}
```

Hata bağlamındaki 4. satıra bakıldığında, sorunu tespit edebiliriz: doğru `input` yönergesi yerine `inputs` kullanıyoruz. Nextflow VSCode uzantısı bunu da işaretleyecektir:

![Geçersiz process mesajı](img/invalid_process_message.png)

#### Kodu düzeltin

[Belgelere](https://www.nextflow.io/docs/latest/process.html#) başvurarak yanlış anahtar kelimeyi doğru olanla değiştirin:

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

        // Girdi kanalı oluştur
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // İşlemi girdi kanalıyla çağır
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

        // Girdi kanalı oluştur
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // İşlemi girdi kanalıyla çağır
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

### 1.3. Kötü değişken isimleri kullanma

Script bloklarınızda kullandığınız değişken isimleri geçerli olmalıdır, girdilerden veya script'ten önce eklenen groovy kodundan türetilmelidir. Ancak pipeline geliştirmenin başında karmaşıklıkla boğuşurken değişken adlandırmada hata yapmak kolaydır ve Nextflow bunu size hızlıca bildirecektir.

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

Hata derleme zamanında yakalanır ve doğrudan 17. satırdaki tanımsız değişkeni gösterir, sorunun tam olarak nerede olduğunu gösteren bir şapka ile.

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

Hata mesajı değişkenin script şablonunda tanınmadığını belirtir ve işte orada- script bloğunda kullanılan ancak başka yerde tanımlanmamış `${undefined_var}` görmelisiniz.

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

### 1.4. Bash değişkenlerinin kötü kullanımı

Nextflow'a başlarken, Nextflow (Groovy) ve Bash değişkenleri arasındaki farkı anlamak zor olabilir. Bu, script bloğunun Bash içeriğinde değişkenleri kullanmaya çalışırken görünen başka bir kötü değişken hatası biçimi oluşturabilir.

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

Hata, `${prefix}` kullanılan 13. satırı gösterir. Neyin soruna neden olduğunu görmek için `bad_bash_var.nf` dosyasını inceleyelim:

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
    echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
    """
}
```

Bu örnekte, `prefix` değişkenini Bash'te tanımlıyoruz, ancak bir Nextflow process'inde ona atıfta bulunmak için kullandığımız `$` sözdizimi (`${prefix}`) Bash değil Groovy değişkeni olarak yorumlanır. Değişken Groovy bağlamında mevcut olmadığından, 'böyle bir değişken yok' hatası alırız.

#### Kodu düzeltin

Bash değişkeni kullanmak istiyorsanız, dolar işaretini şu şekilde kaçırmalısınız:

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
        echo "Processing ${sample_name}" > \${prefix}.txt  # Fixed: Escaped the dollar sign
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
        echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
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

    String birleştirme veya önek/sonek işlemleri gibi basit değişken manipülasyonları için, script bloğundaki Bash değişkenleri yerine script bölümünde Groovy değişkenleri kullanmak genellikle daha okunabilir:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    Bu yaklaşım dolar işaretlerini kaçırmaya gerek bırakmaz ve kodu okumayı ve sürdürmeyi kolaylaştırır.

### 1.5. Workflow Bloğu Dışında İfadeler

Nextflow VSCode uzantısı, hatalara neden olacak kod yapısıyla ilgili sorunları vurgular. Yaygın bir örnek, `workflow {}` bloğunun dışında kanal tanımlamaktır - bu artık bir sözdizimi hatası olarak zorlanmaktadır.

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

Hata mesajı sorunu açıkça belirtir: ifadeler (kanal tanımları gibi) bir workflow veya process bloğunun dışında script bildirimleriyle karıştırılamaz.

#### Kodu kontrol edin

Hataya neyin neden olduğunu görmek için `badpractice_syntax.nf` dosyasını inceleyelim:

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

VSCode uzantısı ayrıca `input_ch` değişkenini workflow bloğunun dışında tanımlandığı için vurgulayacaktır:

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

Girdi kanallarınızı workflow bloğu içinde tanımlanmış tutun ve genel olarak uzantının yaptığı diğer önerileri takip edin.

### Çıkarımlar

Nextflow hata mesajlarını ve IDE görsel göstergelerini kullanarak sözdizimi hatalarını sistematik olarak belirleyebilir ve düzeltebilirsiniz. Yaygın sözdizimi hataları arasında eksik parantezler, yanlış process anahtar kelimeleri, tanımsız değişkenler ve Bash ile Nextflow değişkenlerinin uygunsuz kullanımı bulunur. VSCode uzantısı bunların çoğunu çalışma zamanından önce yakalamanıza yardımcı olur. Bu sözdizimi hata ayıklama becerileri araç setinizde olduğunda, en yaygın Nextflow sözdizimi hatalarını hızlı bir şekilde çözebilecek ve daha karmaşık çalışma zamanı sorunlarıyla başa çıkmaya geçebileceksiniz.

### Sırada ne var?

Sözdizimi doğru olsa bile ortaya çıkan daha karmaşık kanal yapısı hatalarında hata ayıklamayı öğrenin.

---

## 2. Kanal Yapısı Hataları

Kanal yapısı hataları sözdizimi hatalarından daha inceliklidir çünkü kod sözdizimsel olarak doğrudur, ancak veri şekilleri process'lerin beklediği ile eşleşmez. Nextflow pipeline'ı çalıştırmaya çalışacaktır, ancak girdi sayısının beklediği ile eşleşmediğini bulabilir ve başarısız olabilir. Bu hatalar genellikle yalnızca çalışma zamanında görünür ve iş akışınız boyunca akan verilerin anlaşılmasını gerektirir.

!!! tip "`.view()` ile Kanallarda Hata Ayıklama"

    Bu bölüm boyunca, iş akışınızdaki herhangi bir noktada kanal içeriğini incelemek için `.view()` operatörünü kullanabileceğinizi unutmayın. Bu, kanal yapısı sorunlarını anlamak için en güçlü hata ayıklama araçlarından biridir. Bu tekniği 2.4. bölümünde ayrıntılı olarak keşfedeceğiz, ancak örnekler üzerinde çalışırken kullanmaktan çekinmeyin.

    ```groovy
    my_channel.view()  // Kanaldan ne aktığını gösterir
    ```

### 2.1. Yanlış Sayıda Girdi Kanalı

Bu hata, bir process'in beklediğinden farklı sayıda kanal geçirdiğinizde oluşur.

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

Hata mesajı, çağrının 1 argüman beklediğini ancak 2 aldığını açıkça belirtir ve 23. satırı gösterir. `bad_number_inputs.nf` dosyasını inceleyelim:

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

Process yalnızca bir tane tanımlarken birden fazla girdi kanalı sağlayan uyumsuz `PROCESS_FILES` çağrısını görmelisiniz. VSCode uzantısı ayrıca process çağrısının altını kırmızıya çizer ve fare ile üzerine geldiğinizde bir tanı mesajı sağlar:

![Yanlış sayıda argüman mesajı](img/incorrect_num_args.png)

#### Kodu düzeltin

Bu özel örnek için, process tek bir kanal bekler ve ikinci kanala ihtiyaç duymaz, bu nedenle yalnızca `samples_ch` kanalını geçirerek düzeltebiliriz:

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

Bu örnekten daha yaygın olarak, bir process'e ek girdiler ekleyebilir ve workflow çağrısını buna göre güncellemeyi unutabilirsiniz, bu da bu tür bir hataya yol açabilir. Neyse ki, hata mesajı uyumsuzluk hakkında oldukça açık olduğundan, bu anlaşılması ve düzeltilmesi daha kolay hatalardan biridir.

### 2.2. Kanal Tükenmesi (Process Beklenenden Daha Az Çalışır)

Bazı kanal yapısı hataları çok daha inceliklidir ve hiç hata üretmezler. Muhtemelen bunların en yaygını, yeni Nextflow kullanıcılarının queue kanallarının tükenebileceğini ve öğelerin bitebileceğini anlamalarındaki zorluğu yansıtır, bu da iş akışının erken bitmesi anlamına gelir.

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

Process üç kez yerine yalnızca bir kez çalışır çünkü `reference_ch` kanalı, ilk process yürütmesinden sonra tükenen bir queue kanalıdır. Bir kanal tükendiğinde, diğer kanalların hala öğeleri olsa bile tüm process durur.

Bu, birden çok örnek arasında yeniden kullanılması gereken tek bir referans dosyanız olduğu yaygın bir desendir. Çözüm, referans kanalını süresiz olarak yeniden kullanılabilecek bir value kanalına dönüştürmektir.

#### Kodu düzeltin

Kaç dosyanın etkilendiğine bağlı olarak bunu ele almanın birkaç yolu vardır.

**Seçenek 1**: Çok fazla yeniden kullandığınız tek bir referans dosyanız var. Defalarca kullanılabilecek bir value kanal türü oluşturabilirsiniz. Bunu yapmanın üç yolu vardır:

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

**Seçenek 2**: Daha karmaşık senaryolarda, belki de sample kanalındaki tüm örnekler için birden fazla referans dosyanız olduğunda, iki kanalı tuple'lara birleştiren yeni bir kanal oluşturmak için `combine` operatörünü kullanabilirsiniz:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Creates cartesian product

    PROCESS_FILES(combined_ch)
}
```

`.combine()` operatörü iki kanalın kartezyen çarpımını oluşturur, bu nedenle `reference_ch` içindeki her öğe `input_ch` içindeki her öğe ile eşleştirilecektir. Bu, process'in her örnek için çalışmasına izin verirken yine de referansı kullanmasını sağlar.

Bu, process girdisinin ayarlanmasını gerektirir. Örneğimizde, process tanımının başlangıcının aşağıdaki gibi ayarlanması gerekecektir:

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

Artık yalnızca bir tane yerine üç örneğin de işlendiğini görmelisiniz.

### 2.3. Yanlış Kanal İçerik Yapısı

İş akışları belirli bir karmaşıklık düzeyine ulaştığında, her kanalın iç yapılarını takip etmek biraz zor olabilir ve insanlar yaygın olarak process'in beklediği ile kanalın gerçekte içerdiği arasında uyumsuzluklar yaratırlar. Bu, daha önce tartıştığımız, kanal sayısının yanlış olduğu sorundan daha inceliklidir. Bu durumda, doğru sayıda girdi kanalına sahip olabilirsiniz, ancak bu kanallardan birinin veya daha fazlasının iç yapısı process'in beklediği ile eşleşmez.

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

Hata mesajındaki köşeli parantezler burada ipucunu sağlar - process tuple'ı tek bir değer olarak ele alır, bu bizim istediğimiz şey değil. `bad_channel_shape.nf` dosyasını inceleyelim:

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

Tuple'lardan oluşan bir kanal oluşturduğumuzu görebilirsiniz: `['sample1', 'file1.txt']`, ancak process tek bir değer bekliyor, `val sample_name`. Yürütülen komut, process'in `[sample3, file3.txt]_output.txt` adında bir dosya oluşturmaya çalıştığını gösterir, bu amaçlanan çıktı değildir.

#### Kodu düzeltin

Bunu düzeltmek için, process her iki girdiyi de gerektiriyorsa process'i bir tuple kabul edecek şekilde ayarlayabiliriz:

=== "Seçenek 1: Process'te tuple kabul et"

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

Kanallar için en güçlü hata ayıklama aracı `.view()` operatörüdür. `.view()` ile hata ayıklamaya yardımcı olmak için kanallarınızın şeklini tüm aşamalarda anlayabilirsiniz.

#### Pipeline'ı çalıştırın

Bunu eylemde görmek için `bad_channel_shape_viewed.nf` dosyasını çalıştırın:

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

Gelecekte kanal içeriğini anlamak için `.view()` işlemlerini aşırı kullanmaktan sizi kurtarmak için yardımcı olması için bazı yorumlar eklemeniz tavsiye edilir:

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

Bu, iş akışlarınız karmaşıklık olarak büyüdükçe ve kanal yapısı daha opak hale geldikçe daha önemli hale gelecektir.

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

### Çıkarımlar

Birçok kanal yapısı hatası geçerli Nextflow sözdizimi ile oluşturulabilir. Veri akışını anlayarak, inceleme için `.view()` operatörlerini kullanarak ve beklenmedik tuple yapılarını gösteren köşeli parantezler gibi hata mesajı desenlerini tanıyarak kanal yapısı hatalarında hata ayıklayabilirsiniz.

### Sırada ne var?

Process tanımlarından kaynaklanan hatalar hakkında bilgi edinin.

---

## 3. Process Yapısı Hataları

Process'lerle ilgili karşılaşacağınız hataların çoğu, komutu oluştururken yaptığınız hatalara veya temel yazılımla ilgili sorunlara bağlı olacaktır. Bununla birlikte, yukarıdaki kanal sorunlarına benzer şekilde, sözdizimi hatası olarak nitelendirilmeyen ancak çalışma zamanında hatalara neden olacak process tanımında hatalar yapabilirsiniz.

### 3.1. Eksik Çıktı Dosyaları

Process'leri yazarken yaygın bir hata, process'in beklediği ile oluşturulan arasında uyumsuzluk yaratan bir şey yapmaktır.

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

Hata mesajı, process'in `sample3.txt` adında bir çıktı dosyası üretmesini beklediğini, ancak script'in aslında `sample3_output.txt` oluşturduğunu gösterir. `missing_output.nf` içindeki process tanımını inceleyelim:

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

`output:` bloğundaki çıktı dosya adı ile script'te kullanılan dosya adı arasında bir uyumsuzluk olduğunu görmelisiniz. Bu uyumsuzluk process'in başarısız olmasına neden olur. Bu tür bir hatayla karşılaşırsanız, geri dönüp process tanımınız ile çıktı bloğunuz arasındaki çıktıların eşleşip eşleşmediğini kontrol edin.

Sorun hâlâ net değilse, oluşturulan gerçek çıktı dosyalarını belirlemek için çalışma dizininin kendisini kontrol edin:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

Bu örnekte bu, `output:` tanımımızın aksine çıktı dosya adına bir `_output` son ekinin dahil edildiğini bize gösterecektir.

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

### 3.2. Eksik Yazılım

Bir diğer hata sınıfı, yazılım sağlama hatalarından kaynaklanır. `missing_software.nf` sözdizimsel olarak geçerli bir iş akışıdır, ancak kullandığı `cowpy` komutunu sağlamak için bazı harici yazılımlara bağımlıdır.

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

Process, belirttiğimiz komuta erişemiyor. Bazen bunun nedeni, bir script'in iş akışının `bin` dizininde bulunması ancak çalıştırılabilir yapılmamış olmasıdır. Diğer zamanlarda ise yazılımın, iş akışının çalıştığı container veya ortamda yüklü olmamasıdır.

#### Kodu kontrol edin

O `127` çıkış koduna dikkat edin - size sorunu tam olarak söyler. `missing_software.nf` dosyasını inceleyelim:

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

Burada biraz sahtekârlık yaptık ve aslında kodda yanlış bir şey yok. Sadece process'i söz konusu komuta erişebilecek şekilde çalıştırmak için gerekli yapılandırmayı belirtmemiz gerekiyor. Bu durumda process'in bir container tanımı var, bu yüzden tek yapmamız gereken iş akışını Docker etkin şekilde çalıştırmak.

#### Pipeline'ı çalıştırın

`nextflow.config` içinde sizin için bir Docker profili oluşturduk, böylece iş akışını şu şekilde çalıştırabilirsiniz:

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

### 3.3. Kötü Kaynak Yapılandırması

Üretim kullanımında, process'leriniz üzerinde kaynakları yapılandırıyor olacaksınız. Örneğin `memory`, process'inize sunulan maksimum bellek miktarını tanımlar ve process bunu aşarsa, zamanlayıcınız genellikle process'i sonlandırır ve `137` çıkış kodu döndürür. `local` yürütücüyü kullandığımız için bunu burada gösteremeyiz, ancak `time` ile benzer bir şey gösterebiliriz.

#### Pipeline'ı çalıştırın

`bad_resources.nf`, 1 milisaniyelik gerçekçi olmayan bir zaman sınırına sahip process yapılandırmasına sahiptir:

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

    time '1 ms'  // HATA: Gerçekçi olmayan zaman sınırı

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // 1 saniye sürer, ancak zaman sınırı 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

Process'in bir saniyeden fazla süreceğini biliyoruz (emin olmak için bir sleep ekledik), ancak process 1 milisaniye sonra zaman aşımına uğrayacak şekilde ayarlanmış. Birisi yapılandırmasında biraz gerçekçi olmayan davranmış!

#### Kodu düzeltin

Zaman sınırını gerçekçi bir değere yükseltin:

=== "Sonra"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // Düzeltildi: Gerçekçi zaman sınırı

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

        time '1 ms'  // HATA: Gerçekçi olmayan zaman sınırı

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // 1 saniye sürer, ancak zaman sınırı 1ms
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

Hata mesajlarınızı okuduğunuzdan emin olursanız, bu tür başarısızlıklar sizi çok uzun süre şaşırtmamalıdır. Ancak kaynak direktiflerinizi uygun şekilde yapılandırabilmeniz için çalıştırdığınız komutların kaynak gereksinimlerini anladığınızdan emin olun.

### 3.4. Process Hata Ayıklama Teknikleri

Process'ler başarısız olduğunda veya beklenmedik şekilde davrandığında, neyin yanlış gittiğini araştırmak için sistematik tekniklere ihtiyacınız var. Çalışma dizini, process yürütmesini hata ayıklamak için gereken tüm bilgileri içerir.

#### Çalışma Dizini İncelemesini Kullanma

Process'ler için en güçlü hata ayıklama aracı, çalışma dizinini incelemektir. Bir process başarısız olduğunda, Nextflow o belirli process yürütmesi için ne olduğunu anlamak için gereken tüm dosyaları içeren bir çalışma dizini oluşturur.

#### Pipeline'ı çalıştırın

Çalışma dizini incelemesini göstermek için daha önceki `missing_output.nf` örneğini kullanalım (gerekirse bir çıktı adlandırma uyumsuzluğunu yeniden oluşturun):

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

`.command.sh` dosyası tam olarak hangi komutun çalıştırıldığını gösterir:

```bash
# Çalıştırılan komutu görüntüleyin
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

Bu şunları ortaya çıkarır:

- **Değişken ikamesi**: Nextflow değişkenlerinin düzgün şekilde genişletilip genişletilmediği
- **Dosya yolları**: Girdi dosyalarının doğru şekilde konumlandırılıp konumlandırılmadığı
- **Komut yapısı**: Script sözdiziminin doğru olup olmadığı

Dikkat edilmesi gereken yaygın sorunlar:

- **Eksik tırnak işaretleri**: Boşluk içeren değişkenler uygun tırnak işaretleri gerektirir
- **Yanlış dosya yolları**: Var olmayan veya yanlış konumlardaki girdi dosyaları
- **Hatalı değişken adları**: Değişken referanslarındaki yazım hataları
- **Eksik ortam kurulumu**: Belirli ortamlara bağımlı komutlar

##### Hata Çıktısını Kontrol Edin

`.command.err` dosyası gerçek hata mesajlarını içerir:

```bash
# Hata çıktısını görüntüleyin
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

Bu dosya şunları gösterecektir:

- **Çıkış kodları**: 127 (komut bulunamadı), 137 (sonlandırıldı), vb.
- **İzin hataları**: Dosya erişim sorunları
- **Yazılım hataları**: Uygulamaya özgü hata mesajları
- **Kaynak hataları**: Bellek/zaman sınırı aşımı

##### Standart Çıktıyı Kontrol Edin

`.command.out` dosyası komutunuzun ürettiğini gösterir:

```bash
# Standart çıktıyı görüntüleyin
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

Bu şunları doğrulamaya yardımcı olur:

- **Beklenen çıktı**: Komutun doğru sonuçları üretip üretmediği
- **Kısmi yürütme**: Komutun başlayıp yarıda başarısız olup olmadığı
- **Hata ayıklama bilgisi**: Script'inizden gelen tanılama çıktısı

##### Çıkış Kodunu Kontrol Edin

`.exitcode` dosyası process için çıkış kodunu içerir:

```bash
# Çıkış kodunu görüntüleyin
cat work/*/*/.exitcode
```

Yaygın çıkış kodları ve anlamları:

- **Çıkış kodu 127**: Komut bulunamadı - yazılım kurulumunu kontrol edin
- **Çıkış kodu 137**: Process sonlandırıldı - bellek/zaman sınırlarını kontrol edin

##### Dosya Varlığını Kontrol Edin

Process'ler eksik çıktı dosyaları nedeniyle başarısız olduğunda, gerçekte hangi dosyaların oluşturulduğunu kontrol edin:

```bash
# Çalışma dizinindeki tüm dosyaları listeleyin
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

Bu şunları belirlemeye yardımcı olur:

- **Dosya adı uyumsuzlukları**: Beklenenden farklı adlara sahip çıktı dosyaları
- **İzin sorunları**: Oluşturulamayan dosyalar
- **Yol sorunları**: Yanlış dizinlerde oluşturulan dosyalar

Daha önceki örneğimizde bu, beklenen `sample3.txt` dosyası mevcut olmasa da `sample3_output.txt` dosyasının mevcut olduğunu bize doğruladı:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### Çıkarımlar

Process hata ayıklaması, neyin yanlış gittiğini anlamak için çalışma dizinlerini incelemeyi gerektirir. Anahtar dosyalar `.command.sh` (çalıştırılan script), `.command.err` (hata mesajları) ve `.command.out` (standart çıktı) dosyalarını içerir. 127 (komut bulunamadı) ve 137 (process sonlandırıldı) gibi çıkış kodları, başarısızlık türü hakkında anında tanılama ipuçları sağlar.

### Sırada ne var?

Nextflow'un yerleşik hata ayıklama araçları ve sorun giderme için sistematik yaklaşımlar hakkında bilgi edinin.

---

## 4. Yerleşik Hata Ayıklama Araçları ve İleri Teknikler

Nextflow, iş akışı yürütmesini hata ayıklamak ve analiz etmek için çeşitli güçlü yerleşik araçlar sağlar. Bu araçlar neyin yanlış gittiğini, nerede yanlış gittiğini ve nasıl verimli bir şekilde düzeltileceğini anlamanıza yardımcı olur.

### 4.1. Gerçek Zamanlı Process Çıktısı

Bazen çalışan process'lerin içinde neler olduğunu görmeniz gerekir. Her görevin yürütülürken tam olarak ne yaptığını gösteren gerçek zamanlı process çıktısını etkinleştirebilirsiniz.

#### Pipeline'ı çalıştırın

Daha önceki örneklerimizden `bad_channel_shape_viewed.nf`, `.view()` kullanarak kanal içeriğini yazdırıyordu, ancak process'in kendisi içinden değişkenleri yansıtmak için `debug` direktifini de kullanabiliriz; bunu `bad_channel_shape_viewed_debug.nf` dosyasında gösteriyoruz. İş akışını çalıştırın:

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

`debug` direktifinin nasıl çalıştığını görmek için `bad_channel_shape_viewed_debug.nf` dosyasını inceleyelim:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Gerçek zamanlı çıktıyı etkinleştir

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

`debug` direktifi, bir process'in ortamını anlamak için hızlı ve pratik bir yol olabilir.

### 4.2. Önizleme Modu

Bazen herhangi bir process çalışmadan önce sorunları yakalamak istersiniz. Nextflow bu tür proaktif hata ayıklama için bir bayrak sağlar: `-preview`.

#### Pipeline'ı çalıştırın

Önizleme modu, komutları çalıştırmadan iş akışı mantığını test etmenizi sağlar. Bu, iş akışınızın yapısını hızlıca kontrol etmek ve gerçek komutları çalıştırmadan process'lerin doğru şekilde bağlandığından emin olmak için oldukça faydalı olabilir.

!!! note

    Daha önce `bad_syntax.nf` dosyasını düzelttiyseniz, bu komutu çalıştırmadan önce script bloğundan sonra kapanış süslü parantezini kaldırarak sözdizimi hatasını yeniden ekleyin.

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

Önizleme modu, herhangi bir process çalıştırmadan sözdizimi hatalarını erken yakalamak için özellikle faydalıdır. Yürütmeden önce iş akışı yapısını ve process bağlantılarını doğrular.

### 4.3. Mantık Testi için Stub Çalıştırma

Bazen komutlar çok uzun sürdüğü, özel yazılım gerektirdiği veya karmaşık nedenlerle başarısız olduğu için hataların ayıklanması zordur. Stub çalıştırma, gerçek komutları yürütmeden iş akışı mantığını test etmenizi sağlar.

#### Pipeline'ı çalıştırın

Bir Nextflow process'i geliştirirken, gerçek komutu çalıştırmadan doğru biçimde çıktılar üreten 'sahte' komutlar tanımlamak için `stub` direktifini kullanabilirsiniz. Bu yaklaşım, gerçek yazılımın karmaşıklıklarıyla uğraşmadan önce iş akışı mantığınızın doğru olduğunu doğrulamak istediğinizde özellikle değerlidir.

Örneğin, daha önceki `missing_software.nf` dosyamızı hatırlıyor musunuz? `-profile docker` ekleyene kadar iş akışının çalışmasını engelleyen eksik yazılıma sahip olan? `missing_software_with_stub.nf` çok benzer bir iş akışıdır. Aynı şekilde çalıştırırsak, aynı hatayı üreteceğiz:

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

Ancak, `docker` profili olmadan bile `-stub-run` ile çalıştırırsak bu iş akışı hata üretmeyecektir:

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

`missing_software.nf` ile karşılaştırıldığında, bu process'te Nextflow stub modunda çalıştırıldığında `script:` bloğunda belirtilen komut yerine kullanılacak bir komutu belirten bir `stub:` direktifi bulunur.

Burada kullandığımız `touch` komutu herhangi bir yazılıma veya uygun girdilere bağlı değildir ve tüm durumlarda çalışır; bu sayede process iç yapısı hakkında endişelenmeden iş akışı mantığını ayıklamamıza olanak tanır.

**Stub çalıştırma şunların hata ayıklamasına yardımcı olur:**

- Kanal yapısı ve veri akışı
- Process bağlantıları ve bağımlılıkları
- Parametre yayılımı
- Yazılım bağımlılıkları olmadan iş akışı mantığı

### 4.4. Sistematik Hata Ayıklama Yaklaşımı

Artık bireysel hata ayıklama tekniklerini öğrendiğinize göre - izleme dosyaları ve çalışma dizinlerinden önizleme moduna, stub çalıştırmaya ve kaynak izlemeye kadar - bunları sistematik bir metodoloji içinde bir araya getirelim. Yapılandırılmış bir yaklaşıma sahip olmak, karmaşık hatalar karşısında bunalmanızı önler ve önemli ipuçlarını kaçırmamanızı sağlar.

Bu metodoloji, ele aldığımız tüm araçları verimli bir iş akışında birleştirir:

**Dört Aşamalı Hata Ayıklama Yöntemi:**

**Aşama 1: Sözdizimi Hatası Çözümü (5 dakika)**

1. VSCode veya IDE'nizde kırmızı alt çizgileri kontrol edin
2. Sözdizimi sorunlarını belirlemek için `nextflow run workflow.nf -preview` komutunu çalıştırın
3. Tüm sözdizimi hatalarını düzeltin (eksik süslü parantezler, sondaki virgüller, vb.)
4. Devam etmeden önce iş akışının başarıyla ayrıştırıldığından emin olun

**Aşama 2: Hızlı Değerlendirme (5 dakika)**

1. Çalışma zamanı hata mesajlarını dikkatlice okuyun
2. Çalışma zamanı, mantık veya kaynak hatası olup olmadığını kontrol edin
3. Temel iş akışı mantığını test etmek için önizleme modunu kullanın

**Aşama 3: Ayrıntılı İnceleme (15-30 dakika)**

1. Başarısız görevin çalışma dizinini bulun
2. Günlük dosyalarını inceleyin
3. Kanalları incelemek için `.view()` operatörleri ekleyin
4. Yürütme olmadan iş akışı mantığını test etmek için `-stub-run` kullanın

**Aşama 4: Düzeltme ve Doğrulama (15 dakika)**

1. Minimal hedefli düzeltmeler yapın
2. Resume ile test edin: `nextflow run workflow.nf -resume`
3. Tam iş akışı yürütmesini doğrulayın

!!! tip "Verimli Hata Ayıklama için Resume Kullanımı"

    Bir sorun belirledikten sonra, iş akışınızın başarılı bölümlerini yeniden çalıştırmadan zaman kaybetmeden düzeltmelerinizi test etmenin verimli bir yoluna ihtiyacınız var. Nextflow'un `-resume` işlevselliği hata ayıklama için çok değerlidir.

    [Hello Nextflow](../hello_nextflow/) üzerinde çalıştıysanız `-resume` ile karşılaşmış olacaksınız ve sorun process'inizden önceki process'lerin çalışmasını beklerken kendinizi kurtarmak için hata ayıklama sırasında bundan iyi yararlanmanız önemlidir.

    **Resume hata ayıklama stratejisi:**

    1. İş akışını başarısız olana kadar çalıştırın
    2. Başarısız görev için çalışma dizinini inceleyin
    3. Belirli sorunu düzeltin
    4. Yalnızca düzeltmeyi test etmek için resume kullanın
    5. İş akışı tamamlanana kadar tekrarlayın

#### Hata Ayıklama Yapılandırma Profili

Bu sistematik yaklaşımı daha da verimli hale getirmek için, ihtiyacınız olan tüm araçları otomatik olarak etkinleştiren özel bir hata ayıklama yapılandırması oluşturabilirsiniz:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Hata ayıklama için tutucu kaynaklar
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

Ardından pipeline'ı bu profil etkinleştirilerek çalıştırabilirsiniz:

```bash
nextflow run workflow.nf -profile debug
```

Bu profil gerçek zamanlı çıktıyı etkinleştirir, çalışma dizinlerini korur ve daha kolay hata ayıklama için paralelleştirmeyi sınırlar.

### 4.5. Pratik Hata Ayıklama Alıştırması

Şimdi sistematik hata ayıklama yaklaşımını pratiğe koyma zamanı. `buggy_workflow.nf` iş akışı, gerçek dünya geliştirmede karşılaşacağınız hata türlerini temsil eden birkaç yaygın hata içerir.

!!! exercise

    `buggy_workflow.nf` dosyasındaki tüm hataları belirlemek ve düzeltmek için sistematik hata ayıklama yaklaşımını kullanın. Bu iş akışı, bir CSV dosyasından örnek verileri işlemeye çalışır ancak yaygın hata ayıklama senaryolarını temsil eden birden fazla kasıtlı hata içerir.

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

        Bu şifreli hata, `params{}` bloğundaki 11-12. satır civarında bir ayrıştırma sorununa işaret eder. v2 ayrıştırıcı yapısal sorunları erken yakalar.

    Öğrendiğiniz dört aşamalı hata ayıklama yöntemini uygulayın:

    **Aşama 1: Sözdizimi Hatası Çözümü**
    - VSCode veya IDE'nizde kırmızı alt çizgileri kontrol edin
    - Sözdizimi sorunlarını belirlemek için `nextflow run workflow.nf -preview` komutunu çalıştırın
    - Tüm sözdizimi hatalarını düzeltin (eksik süslü parantezler, sondaki virgüller, vb.)
    - Devam etmeden önce iş akışının başarıyla ayrıştırıldığından emin olun

    **Aşama 2: Hızlı Değerlendirme**
    - Çalışma zamanı hata mesajlarını dikkatlice okuyun
    - Hataların çalışma zamanı, mantık veya kaynak ile ilgili olup olmadığını belirleyin
    - Temel iş akışı mantığını test etmek için `-preview` modunu kullanın

    **Aşama 3: Ayrıntılı İnceleme**
    - Başarısız görevler için çalışma dizinlerini inceleyin
    - Kanalları incelemek için `.view()` operatörleri ekleyin
    - Çalışma dizinlerindeki günlük dosyalarını kontrol edin
    - Yürütme olmadan iş akışı mantığını test etmek için `-stub-run` kullanın

    **Aşama 4: Düzeltme ve Doğrulama**
    - Hedefli düzeltmeler yapın
    - Düzeltmeleri verimli şekilde test etmek için `-resume` kullanın
    - Tam iş akışı yürütmesini doğrulayın

    **Kullanabileceğiniz Hata Ayıklama Araçları:**
    ```bash
    # Sözdizimi kontrolü için önizleme modu
    nextflow run buggy_workflow.nf -preview

    # Ayrıntılı çıktı için hata ayıklama profili
    nextflow run buggy_workflow.nf -profile debug

    # Mantık testi için stub çalıştırma
    nextflow run buggy_workflow.nf -stub-run

    # Düzeltmelerden sonra resume
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution
        `buggy_workflow.nf` dosyası, tüm önemli hata ayıklama kategorilerini kapsayan 9 veya 10 farklı hata içerir (nasıl saydığınıza bağlı). Her hatanın ve nasıl düzeltileceğinin sistematik bir dökümü aşağıdadır

        Sözdizimi hatalarıyla başlayalım:

        **Hata 1: Sözdizimi Hatası - Sondaki Virgül**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // HATA: Sondaki virgül
        ```
        **Düzeltme:** Sondaki virgülü kaldırın
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **Hata 2: Sözdizimi Hatası - Eksik Kapanış Süslü Parantezi**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // HATA: processFiles process'i için eksik kapanış süslü parantezi
        ```
        **Düzeltme:** Eksik kapanış süslü parantezini ekleyin
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Eksik kapanış süslü parantezini ekleyin
        ```

        **Hata 3: Değişken Adı Hatası**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // HATA: sample_id olmalı
        cat ${input_file} > ${sample}_result.txt  // HATA: sample_id olmalı
        ```
        **Düzeltme:** Doğru girdi değişken adını kullanın
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Hata 4: Tanımsız Değişken Hatası**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // HATA: sample_ids tanımsız
        ```
        **Düzeltme:** Doğru kanalı kullanın ve örnek kimliklerini çıkarın
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        Bu noktada iş akışı çalışacak, ancak hâlâ hatalar alacağız (ör. `processFiles` içinde `Path value cannot be null`), kötü kanal yapısından kaynaklanan.

        **Hata 5: Kanal Yapısı Hatası - Yanlış Map Çıktısı**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // HATA: processFiles tuple bekliyor
        ```
        **Düzeltme:** processFiles'ın beklediği tuple yapısını döndürün
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        Ancak bu, yukarıda `heavyProcess()` çalıştırmamızı bozacak, bu yüzden o process'e sadece örnek kimliklerini geçirmek için bir map kullanmamız gerekecek:

        **Hata 6: heavyProcess için kötü kanal yapısı**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // HATA: input_ch artık emisyon başına 2 öğeye sahip - heavyProcess sadece 1'e ihtiyaç duyuyor (ilki)
        ```
        **Düzeltme:** Doğru kanalı kullanın ve örnek kimliklerini çıkarın
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Şimdi biraz daha ilerliyoruz ama `No such variable: i` hatası alıyoruz, çünkü bir Bash değişkenini kaçış karakteriyle yazmadık.

        **Hata 7: Bash Değişken Kaçış Hatası**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // HATA: $i kaçış karakterli değil
        ```
        **Düzeltme:** Bash değişkenini kaçış karakteriyle yazın
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Şimdi `Process exceeded running time limit (1ms)` alıyoruz, bu yüzden ilgili process için çalışma süresi sınırını düzeltiyoruz:

        **Hata 8: Kaynak Yapılandırma Hatası**
        ```groovy linenums="36"
        time '1 ms'  // HATA: Gerçekçi olmayan zaman sınırı
        ```
        **Düzeltme:** Gerçekçi bir zaman sınırına yükseltin
        ```groovy linenums="36"
        time '100 s'
        ```

        Ardından çözülecek bir `Missing output file(s)` hatası var:

        **Hata 9: Çıktı Dosya Adı Uyumsuzluğu**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // HATA: Yanlış dosya adı, çıktı tanımıyla eşleşmeli
        ```
        **Düzeltme:** Çıktı tanımıyla eşleştirin
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        İlk iki process çalıştı, ancak üçüncüsü çalışmadı.

        **Hata 10: Çıktı Dosya Adı Uyumsuzluğu**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Hata: process'ten değil çalışma dizininden girdi almaya çalışıyor
        handleFiles(file_ch)
        ```
        **Düzeltme:** Önceki process'in çıktısını alın
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        Bununla birlikte, tüm iş akışı çalışmalıdır.

        **Tam Düzeltilmiş İş Akışı:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Hata ayıklama alıştırmaları için hatalı iş akışı
        * Bu iş akışı öğrenme amacıyla birkaç kasıtlı hata içerir
        */

        params{
            // Eksik doğrulamalı parametreler
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Girdi/çıktı uyumsuzluğu olan process
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
        * Kaynak sorunları olan process
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
        * Dosya işleme sorunları olan process
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

            // Hatalı kullanımlı kanal
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

- **Sözdizimi hataları**: Eksik süslü parantezler, sondaki virgüller, tanımsız değişkenler
- **Kanal yapısı hataları**: Yanlış veri şekilleri, tanımsız kanallar
- **Process hataları**: Çıktı dosya uyumsuzlukları, değişken kaçış
- **Kaynak hataları**: Gerçekçi olmayan zaman sınırları

**Temel Hata Ayıklama Dersleri:**

1. **Hata mesajlarını dikkatlice okuyun** - genellikle doğrudan soruna işaret ederler
2. **Sistematik yaklaşımlar kullanın** - bir seferde bir hatayı düzeltin ve `-resume` ile test edin
3. **Veri akışını anlayın** - kanal yapısı hataları genellikle en inceliklileridir
4. **Çalışma dizinlerini kontrol edin** - process'ler başarısız olduğunda, günlükler tam olarak neyin yanlış gittiğini söyler

---

## Özet

Bu yan görevde, Nextflow iş akışlarında hata ayıklamak için bir dizi sistematik teknik öğrendiniz.
Bu teknikleri kendi çalışmalarınızda uygulamak, bilgisayarınızla mücadele ederek daha az zaman harcamanızı, sorunları daha hızlı çözmenizi ve gelecekteki sorunlardan kendinizi korumanızı sağlayacaktır.

### Temel kalıplar

**1. Sözdizimi hatalarını belirleme ve düzeltme**:

- Nextflow hata mesajlarını yorumlama ve sorunları bulma
- Yaygın sözdizimi hataları: eksik süslü parantezler, yanlış anahtar kelimeler, tanımsız değişkenler
- Nextflow (Groovy) ve Bash değişkenlerini ayırt etme
- Erken hata tespiti için VS Code uzantı özelliklerini kullanma

```groovy
// Eksik süslü parantez - IDE'de kırmızı alt çizgilere bakın
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- eksik!

// Yanlış anahtar kelime
inputs:  // 'input:' olmalı

// Tanımsız değişken - Bash değişkenleri için ters eğik çizgi ile kaçış yapın
echo "${undefined_var}"      // Nextflow değişkeni (tanımlı değilse hata)
echo "\${bash_var}"          // Bash değişkeni (kaçış karakterli)
```

**2. Kanal yapısı sorunlarını ayıklama**:

- Kanal kardinalitesi ve tükenme sorunlarını anlama
- Kanal içerik yapısı uyumsuzluklarını ayıklama
- Kanal incelemesi için `.view()` operatörlerini kullanma
- Çıktıda köşeli parantez gibi hata kalıplarını tanıma

```groovy
// Kanal içeriğini inceleyin
my_channel.view { "Content: $it" }

// Queue kanalını value kanalına dönüştürün (tükenmeyi önler)
reference_ch = channel.value('ref.fa')
// veya
reference_ch = channel.of('ref.fa').first()
```

**3. Process yürütme sorunlarını giderme**:

- Eksik çıktı dosyası hatalarını teşhis etme
- Çıkış kodlarını anlama (127 eksik yazılım için, 137 bellek sorunları için)
- Çalışma dizinlerini ve komut dosyalarını inceleme
- Kaynakları uygun şekilde yapılandırma

```bash
# Gerçekte neyin çalıştırıldığını kontrol edin
cat work/ab/cdef12/.command.sh

# Hata çıktısını kontrol edin
cat work/ab/cdef12/.command.err

# Çıkış kodu 127 = komut bulunamadı
# Çıkış kodu 137 = sonlandırıldı (bellek/zaman sınırı)
```

**4. Nextflow'un yerleşik hata ayıklama araçlarını kullanma**:

- Önizleme modu ve gerçek zamanlı hata ayıklamadan yararlanma
- Mantık testi için stub çalıştırma uygulama
- Verimli hata ayıklama döngüleri için resume uygulama
- Dört aşamalı sistematik hata ayıklama metodolojisini takip etme

!!! tip "Hızlı Hata Ayıklama Referansı"

    **Sözdizimi hataları mı?** → VSCode uyarılarını kontrol edin, `nextflow run workflow.nf -preview` çalıştırın

    **Kanal sorunları mı?** → İçeriği incelemek için `.view()` kullanın: `my_channel.view()`

    **Process başarısızlıkları mı?** → Çalışma dizini dosyalarını kontrol edin:

    - `.command.sh` - çalıştırılan script
    - `.command.err` - hata mesajları
    - `.exitcode` - çıkış durumu (127 = komut bulunamadı, 137 = sonlandırıldı)

    **Gizemli davranış mı?** → İş akışı mantığını test etmek için `-stub-run` ile çalıştırın

    **Düzeltmeler mi yaptınız?** → Test ederken zaman kazanmak için `-resume` kullanın: `nextflow run workflow.nf -resume`

---

### Ek kaynaklar

- [Nextflow sorun giderme kılavuzu](https://www.nextflow.io/docs/latest/troubleshooting.html): Resmi sorun giderme belgeleri
- [Nextflow kanallarını anlama](https://www.nextflow.io/docs/latest/channel.html): Kanal türleri ve davranışları hakkında derinlemesine bilgi
- [Process direktifleri referansı](https://www.nextflow.io/docs/latest/process.html#directives): Tüm mevcut process yapılandırma seçenekleri
- [nf-test](https://www.nf-test.com/): Nextflow pipeline'ları için test çerçevesi
- [Nextflow Slack topluluğu](https://www.nextflow.io/slack-invite.html): Topluluktan yardım alın

Üretim iş akışları için şunları düşünün:

- Ölçekte izleme ve hata ayıklama için [Seqera Platform](https://seqera.io/platform/) kurulumu
- Tekrarlanabilir yazılım ortamları için [Wave containers](https://seqera.io/wave/) kullanımı

**Unutmayın:** Etkili hata ayıklama, pratikle gelişen bir beceridir. Burada edindiğiniz sistematik metodoloji ve kapsamlı araç seti, Nextflow geliştirme yolculuğunuz boyunca size iyi hizmet edecektir.

---

## Sırada ne var?

[Yan Görevler menüsüne](./index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
