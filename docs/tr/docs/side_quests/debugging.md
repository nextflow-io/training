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

Hata mesajı, process'in `sample3.txt` adında bir çıktı dosyası üretmesini beklediğini, ancak script'in aslında `sample3_output.txt` oluşturduğunu gösterir. `missing_output.nf` içindeki process tanımını inc
