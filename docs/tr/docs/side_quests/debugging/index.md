# İş Akışlarında Hata Ayıklama

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hata ayıklama, saatlerce süren hayal kırıklığından sizi kurtarabilecek ve daha etkili bir Nextflow geliştiricisi olmanızı sağlayacak kritik bir beceridir. Kariyeriniz boyunca, özellikle başlangıç aşamasında, iş akışlarınızı oluştururken ve bakımını yaparken hatalarla karşılaşacaksınız. Sistematik hata ayıklama yaklaşımlarını öğrenmek, sorunları hızla tespit edip çözmenize yardımcı olacaktır.

### Öğrenme hedefleri

Bu yan görevde, Nextflow iş akışları için **sistematik hata ayıklama teknikleri** keşfedeceğiz:

- **Sözdizimi hatası ayıklama**: IDE özelliklerini ve Nextflow hata mesajlarını etkin biçimde kullanma
- **Kanal hata ayıklama**: Veri akışı sorunlarını ve kanal yapısı problemlerini teşhis etme
- **Süreç hata ayıklama**: Yürütme hatalarını ve kaynak sorunlarını araştırma
- **Yerleşik hata ayıklama araçları**: Nextflow'un önizleme modu, stub çalıştırma ve work dizinlerinden yararlanma
- **Sistematik yaklaşımlar**: Verimli hata ayıklama için dört aşamalı bir metodoloji

Sonunda, sinir bozucu hata mesajlarını çözüm için net yol haritalarına dönüştüren sağlam bir hata ayıklama metodolojisine sahip olacaksınız.

### Ön koşullar

Bu yan göreve başlamadan önce:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmalısınız.
- Temel Nextflow kavramları ve mekanizmalarını (süreçler, kanallar, operatörler) rahatça kullanabiliyor olmalısınız.

**İsteğe bağlı:** Önce [IDE Features for Nextflow Development](../dev_environment/) yan görevini tamamlamanızı öneririz.
Bu görev, burada yoğun biçimde kullanacağımız hata ayıklamayı destekleyen IDE özelliklerini (sözdizimi vurgulama, hata tespiti vb.) kapsamlı şekilde ele almaktadır.

---

## 0. Başlarken

#### Eğitim kod alanını açın

Henüz yapmadıysanız, eğitim ortamını [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı şekilde açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Bu eğitim için gerekli dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/debugging
```

VSCode'u bu dizine odaklamak için şu komutu çalıştırabilirsiniz:

```bash
code .
```

#### Materyalleri inceleyin

Alıştırma için kullanacağımız çeşitli hata türleri içeren örnek iş akışları bulacaksınız:

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

Bu dosyalar, gerçek dünya geliştirme sürecinde karşılaşacağınız yaygın hata ayıklama senaryolarını temsil etmektedir.

#### Görevi inceleyin

Göreviniz her iş akışını çalıştırmak, hataları tespit etmek ve düzeltmektir.

Her hatalı iş akışı için:

1. **İş akışını çalıştırın** ve hatayı gözlemleyin
2. **Hata mesajını analiz edin**: Nextflow size ne söylüyor?
3. **Sağlanan ipuçlarını kullanarak** koddaki sorunu bulun
4. **Hatayı düzeltin** ve çözümünüzün çalıştığını doğrulayın
5. **Bir sonraki bölüme geçmeden önce dosyayı sıfırlayın** (`git checkout <dosyaadı>` kullanın)

Alıştırmalar basit sözdizimi hatalarından daha ince çalışma zamanı sorunlarına doğru ilerlemektedir.
Çözümler satır içinde ele alınmaktadır; ancak ileriye okumadan önce her birini kendiniz çözmeye çalışın.

#### Hazırlık kontrol listesi

Başlamaya hazır mısınız?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Kod alanım çalışıyor
- [ ] Çalışma dizinimi uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebildiyseniz, başlayabilirsiniz.

---

## 1. Sözdizimi Hataları

Sözdizimi hataları, Nextflow kodu yazarken karşılaşacağınız en yaygın hata türüdür. Kod, Nextflow DSL'nin beklenen sözdizimi kurallarına uymadığında ortaya çıkarlar. Bu hatalar iş akışınızın hiç çalışmamasına neden olur; bu nedenle bunları hızla tespit edip düzeltmeyi öğrenmek önemlidir.

### 1.1. Eksik süslü parantezler

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

**Sözdizimi hata mesajlarının temel unsurları:**

- **Dosya ve konum**: Hatanın hangi dosya ve satır/sütunda bulunduğunu gösterir (`bad_syntax.nf:24:1`)
- **Hata açıklaması**: Ayrıştırıcının beklemediği şeyi açıklar (`Unexpected input: '<EOF>'`)
- **EOF göstergesi**: `<EOF>` (Dosya Sonu) mesajı, ayrıştırıcının daha fazla içerik beklerken dosyanın sonuna ulaştığını gösterir; bu, kapatılmamış süslü parantezlerin klasik bir işaretidir

#### Kodu inceleyin

Şimdi hataya neyin neden olduğunu anlamak için `bad_syntax.nf` dosyasını inceleyelim:

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
// Süreç için kapanış süslü parantezi eksik

workflow {

    // Girdi kanalı oluştur
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Süreci girdi kanalıyla çağır
    PROCESS_FILES(input_ch)
}
```

Bu örnek için hatanın nerede olduğunu göstermek amacıyla bir yorum bıraktık. Nextflow VSCode eklentisi de eşleşmeyen süslü parantezi kırmızıyla işaretleyerek ve dosyanın erken sonlandığını vurgulayarak size ipuçları vermelidir:

![Hatalı sözdizimi](img/bad_syntax.png)

**Parantez hataları için hata ayıklama stratejisi:**

1. VS Code'un parantez eşleştirme özelliğini kullanın (imleci bir parantezin yanına getirin)
2. Parantezle ilgili mesajlar için Sorunlar panelini kontrol edin
3. Her açılan `{` için karşılık gelen bir kapanış `}` olduğundan emin olun

#### Kodu düzeltin

Yorumu eksik kapanış süslü parantezi ile değiştirin:

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
    }  // Eksik kapanış süslü parantezini ekle

    workflow {

        // Girdi kanalı oluştur
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Süreci girdi kanalıyla çağır
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
    // Süreç için kapanış süslü parantezi eksik

    workflow {

        // Girdi kanalı oluştur
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Süreci girdi kanalıyla çağır
        PROCESS_FILES(input_ch)
    }
    ```

#### Pipeline'ı çalıştırın

Çalıştığını doğrulamak için iş akışını tekrar çalıştırın:

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

### 1.2. Yanlış süreç anahtar sözcükleri veya yönergeleri kullanma

Bir diğer yaygın sözdizimi hatası **geçersiz süreç tanımıdır**. Bu durum, gerekli blokları tanımlamayı unutursanız veya süreç tanımında yanlış yönergeler kullanırsanız ortaya çıkabilir.

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

#### Kodu inceleyin

Hata, "Geçersiz süreç tanımı" olduğunu belirtir ve sorunun çevresindeki bağlamı gösterir. 3-7. satırlara bakıldığında, 4. satırda `inputs:` ifadesinin sorun olduğu görülmektedir. `invalid_process.nf` dosyasını inceleyelim:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // HATA: 'inputs' değil 'input' olmalı
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

    // Süreci girdi kanalıyla çağır
    PROCESS_FILES(input_ch)
}
```

Hata bağlamındaki 4. satıra bakıldığında sorun fark edilebilir: doğru `input` yönergesi yerine `inputs` kullanılmış. Nextflow VSCode eklentisi de bunu işaretleyecektir:

![Geçersiz süreç mesajı](img/invalid_process_message.png)

#### Kodu düzeltin

[Belgelere](https://www.nextflow.io/docs/latest/process.html#) başvurarak yanlış anahtar sözcüğü doğrusuyla değiştirin:

=== "Sonra"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Düzeltildi: 'inputs' yerine 'input' kullanıldı
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

        // Süreci girdi kanalıyla çağır
        PROCESS_FILES(input_ch)
    }
    ```

=== "Önce"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // HATA: 'inputs' değil 'input' olmalı
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

        // Süreci girdi kanalıyla çağır
        PROCESS_FILES(input_ch)
    }
    ```

#### Pipeline'ı çalıştırın

Çalıştığını doğrulamak için iş akışını tekrar çalıştırın:

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

Betik bloklarında kullandığınız değişken adları geçerli olmalı; ya girdilerden türetilmeli ya da betikten önce eklenen Groovy kodundan gelmelidir. Ancak pipeline geliştirmenin başında karmaşıklıkla boğuşurken değişken adlandırmada hata yapmak kolaydır ve Nextflow bunu hemen size bildirir.

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

Hata derleme zamanında yakalanır ve doğrudan 17. satırdaki tanımsız değişkene işaret eder; bir şapka işareti sorunun tam olarak nerede olduğunu gösterir.

#### Kodu inceleyin

`no_such_var.nf` dosyasını inceleyelim:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Betikten önce Groovy kodunda değişkenleri tanımla
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // HATA: undefined_var tanımlı değil
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

Hata mesajı, değişkenin betik şablonunda tanınmadığını belirtir; betik bloğunda `${undefined_var}` kullanıldığını ancak başka bir yerde tanımlanmadığını görebilirsiniz.

#### Kodu düzeltin

'No such variable' (böyle bir değişken yok) hatası alırsanız, değişkeni tanımlayarak (girdi değişken adlarını düzelterek veya betikten önce Groovy kodunu düzenleyerek) ya da gerekli değilse betik bloğundan kaldırarak düzeltebilirsiniz:

=== "Sonra"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Betikten önce Groovy kodunda değişkenleri tanımla
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // undefined_var içeren satır kaldırıldı
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
        // Betikten önce Groovy kodunda değişkenleri tanımla
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // HATA: undefined_var tanımlı değil
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Pipeline'ı çalıştırın

Çalıştığını doğrulamak için iş akışını tekrar çalıştırın:

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

Nextflow'a yeni başlarken, Nextflow (Groovy) ve Bash değişkenleri arasındaki farkı anlamak güç olabilir. Bu durum, betik bloğunun Bash içeriğinde değişken kullanmaya çalışırken ortaya çıkan başka bir hatalı değişken hatasına yol açabilir.

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

#### Kodu inceleyin

Hata, `${prefix}` ifadesinin kullanıldığı 13. satıra işaret etmektedir. Soruna neyin neden olduğunu görmek için `bad_bash_var.nf` dosyasını inceleyelim:

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

Bu örnekte `prefix` değişkenini Bash'te tanımlıyoruz; ancak bir Nextflow sürecinde buna başvurmak için kullandığımız `$` sözdizimi (`${prefix}`), Groovy değişkeni olarak yorumlanır, Bash değişkeni olarak değil. Değişken Groovy bağlamında mevcut olmadığından 'no such variable' (böyle bir değişken yok) hatası alırız.

#### Kodu düzeltin

Bir Bash değişkeni kullanmak istiyorsanız, dolar işaretini şu şekilde kaçırmanız gerekir:

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

Çalıştığını doğrulamak için iş akışını tekrar çalıştırın:

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

!!! tip "Groovy ve Bash Değişkenleri"

    Dize birleştirme veya önek/sonek işlemleri gibi basit değişken manipülasyonları için, betik bloğundaki Bash değişkenleri yerine betik bölümünde Groovy değişkenleri kullanmak genellikle daha okunabilirdir:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    Bu yaklaşım, dolar işaretlerini kaçırma ihtiyacını ortadan kaldırır ve kodu okumayı ve bakımını yapmayı kolaylaştırır.

### 1.5. İş Akışı Bloğu Dışındaki İfadeler

Nextflow VSCode eklentisi, hatalara neden olacak kod yapısı sorunlarını vurgular. Yaygın bir örnek, `workflow {}` bloğu dışında kanal tanımlamaktır; bu artık sözdizimi hatası olarak uygulanmaktadır.

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

Hata mesajı sorunu açıkça belirtmektedir: ifadeler (kanal tanımları gibi) bir iş akışı veya süreç bloğu dışında betik bildirimleriyle karıştırılamaz.

#### Kodu inceleyin

Hataya neyin neden olduğunu görmek için `badpractice_syntax.nf` dosyasını inceleyelim:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // HATA: Kanal iş akışı dışında tanımlanmış

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Betikten önce Groovy kodunda değişkenleri tanımla
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

VSCode eklentisi de `input_ch` değişkeninin iş akışı bloğu dışında tanımlandığını vurgulayacaktır:

![Önemsiz sözdizimi hatası](img/nonlethal.png)

#### Kodu düzeltin

Kanal tanımını iş akışı bloğunun içine taşıyın:

=== "Sonra"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Betikten önce Groovy kodunda değişkenleri tanımla
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // İş akışı bloğunun içine taşındı
        PROCESS_FILES(input_ch)
    }
    ```

=== "Önce"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // HATA: Kanal iş akışı dışında tanımlanmış

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Betikten önce Groovy kodunda değişkenleri tanımla
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

Girdi kanallarınızı iş akışı bloğu içinde tanımlayın ve genel olarak eklentinin yaptığı diğer önerilere uyun.

### Özetle

Nextflow hata mesajlarını ve IDE görsel göstergelerini kullanarak sözdizimi hatalarını sistematik biçimde tespit edip düzeltebilirsiniz. Yaygın sözdizimi hataları arasında eksik süslü parantezler, yanlış süreç anahtar sözcükleri, tanımsız değişkenler ve Bash ile Nextflow değişkenlerinin hatalı kullanımı yer almaktadır. VSCode eklentisi, bunların çoğunu çalışma zamanından önce yakalamaya yardımcı olur. Bu sözdizimi hata ayıklama becerileriyle, en yaygın Nextflow sözdizimi hatalarını hızla çözebilir ve daha karmaşık çalışma zamanı sorunlarıyla başa çıkmaya geçebilirsiniz.

### Sırada ne var?

Sözdizimi doğru olsa bile ortaya çıkan daha karmaşık kanal yapısı hatalarını nasıl ayıklayacağınızı öğrenin.

---

## 2. Kanal Yapısı Hataları

Kanal yapısı hataları, sözdizimi hatalarından daha ince bir yapıya sahiptir; çünkü kod sözdizimsel olarak doğrudur, ancak veri şekilleri süreçlerin beklediğiyle eşleşmez. Nextflow pipeline'ı çalıştırmaya çalışır; ancak girdi sayısının beklediğiyle eşleşmediğini fark edip başarısız olabilir. Bu hatalar genellikle yalnızca çalışma zamanında ortaya çıkar ve iş akışınızda akan verilerin anlaşılmasını gerektirir.

!!! tip "`.view()` ile Kanal Hata Ayıklama"

    Bu bölüm boyunca, iş akışınızın herhangi bir noktasında kanal içeriğini incelemek için `.view()` operatörünü kullanabileceğinizi unutmayın. Bu, kanal yapısı sorunlarını anlamak için en güçlü hata ayıklama araçlarından biridir. Bu tekniği 2.4. bölümde ayrıntılı olarak ele alacağız; ancak örnekler üzerinde çalışırken kullanmaktan çekinmeyin.

    ```groovy
    my_channel.view()  // Kanaldan geçen verileri gösterir
    ```

### 2.1. Yanlış Sayıda Girdi Kanalı

Bu hata, bir sürecin beklediğinden farklı sayıda kanal ilettiğinizde ortaya çıkar.

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

#### Kodu inceleyin

Hata mesajı, çağrının 1 bağımsız değişken beklediğini ancak 2 aldığını açıkça belirtir ve 23. satıra işaret eder. `bad_number_inputs.nf` dosyasını inceleyelim:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Süreç yalnızca 1 girdi bekliyor

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // İki ayrı kanal oluştur
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // HATA: 2 kanal iletiliyor ancak süreç yalnızca 1 bekliyor
    PROCESS_FILES(samples_ch, files_ch)
}
```

Süreç yalnızca bir girdi kanalı tanımlarken birden fazla girdi kanalı sağlayan eşleşmeyen `PROCESS_FILES` çağrısını görebilirsiniz. VSCode eklentisi de süreç çağrısının altını kırmızıyla çizer ve üzerine geldiğinizde bir tanı mesajı gösterir:

![Yanlış bağımsız değişken sayısı mesajı](img/incorrect_num_args.png)

#### Kodu düzeltin

Bu özel örnek için süreç tek bir kanal beklemekte ve ikinci kanala ihtiyaç duymamaktadır; bu nedenle yalnızca `samples_ch` kanalını ileterek düzeltebiliriz:

=== "Sonra"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Süreç yalnızca 1 girdi bekliyor

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // İki ayrı kanal oluştur
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Düzeltildi: Yalnızca sürecin beklediği kanal iletildi
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Önce"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Süreç yalnızca 1 girdi bekliyor

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // İki ayrı kanal oluştur
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // HATA: 2 kanal iletiliyor ancak süreç yalnızca 1 bekliyor
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

Bu örnekten daha yaygın olarak, bir sürece ek girdiler ekleyip iş akışı çağrısını buna göre güncellemeyi unutabilirsiniz; bu da bu tür bir hataya yol açabilir. Neyse ki bu, anlaşılması ve düzeltilmesi en kolay hatalardan biridir; çünkü hata mesajı uyuşmazlık hakkında oldukça açıktır.

### 2.2. Kanal Tükenmesi (Süreç Beklenenden Az Çalışıyor)

Bazı kanal yapısı hataları çok daha ince bir yapıya sahiptir ve hiçbir hata üretmez. Bunların en yaygını, yeni Nextflow kullanıcılarının queue channel'larının tükenebileceğini ve öğelerinin biteceğini anlamakta yaşadığı güçlüğü yansıtır; bu da iş akışının erken bitmesine neden olur.

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

Bu iş akışı hatasız tamamlanır; ancak yalnızca tek bir örneği işler!

#### Kodu inceleyin

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
    // Betikten önce Groovy kodunda değişkenleri tanımla
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

`reference_ch` kanalı, ilk süreç yürütmesinden sonra tükenen bir queue channel olduğundan süreç üç kez yerine yalnızca bir kez çalışır. Bir kanal tükendiğinde, diğer kanallarda hâlâ öğeler olsa bile tüm süreç durur.

Bu, birden fazla örnek için yeniden kullanılması gereken tek bir referans dosyanızın olduğu yaygın bir kalıptır. Çözüm, referans kanalını süresiz olarak yeniden kullanılabilecek bir value channel'a dönüştürmektir.

#### Kodu düzeltin

Kaç dosyanın etkilendiğine bağlı olarak bunu ele almanın birkaç yolu vardır.

**Seçenek 1**: Çok sık yeniden kullandığınız tek bir referans dosyanız var. Defalarca kullanılabilen bir value channel türü oluşturabilirsiniz. Bunu yapmanın üç yolu vardır:

**1a** `channel.value()` kullanın:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Value channel yeniden kullanılabilir
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** `first()` [operatörünü](https://www.nextflow.io/docs/latest/reference/operator.html#first) kullanın:

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Value channel'a dönüştür
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** `collect()` [operatörünü](https://www.nextflow.io/docs/latest/reference/operator.html#collect) kullanın:

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Value channel'a dönüştür
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Seçenek 2**: Örnek kanalındaki tüm örnekler için birden fazla referans dosyanızın olduğu daha karmaşık senaryolarda, iki kanalı demetler halinde birleştiren yeni bir kanal oluşturmak için `combine` operatörünü kullanabilirsiniz:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Kartezyen çarpım oluşturur

    PROCESS_FILES(combined_ch)
}
```

`.combine()` operatörü iki kanalın kartezyen çarpımını üretir; böylece `reference_ch` içindeki her öğe `input_ch` içindeki her öğeyle eşleştirilir. Bu, sürecin referansı kullanmaya devam ederken her örnek için çalışmasına olanak tanır.

Bu, süreç girdisinin ayarlanmasını gerektirir. Örneğimizde, süreç tanımının başı şu şekilde ayarlanmalıdır:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

Bu yaklaşım her durumda uygun olmayabilir.

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

Artık yalnızca bir örnek yerine üç örneğin de işlendiğini görmelisiniz.

### 2.3. Yanlış Kanal İçeriği Yapısı

İş akışları belirli bir karmaşıklık düzeyine ulaştığında, her kanalın iç yapısını takip etmek biraz güçleşebilir ve insanlar genellikle sürecin beklediği ile kanalın gerçekte içerdiği arasında uyuşmazlıklar oluşturur. Bu, kanal sayısının yanlış olduğu daha önce ele aldığımız sorundan daha ince bir yapıya sahiptir. Bu durumda, doğru sayıda girdi kanalına sahip olabilirsiniz; ancak bir veya daha fazla kanalın iç yapısı sürecin beklediğiyle eşleşmez.

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

#### Kodu inceleyin

Hata mesajındaki köşeli parantezler burada ipucu vermektedir; süreç demeti tek bir değer olarak ele almaktadır ki bu istediğimiz şey değildir. `bad_channel_shape.nf` dosyasını inceleyelim:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Tek değer bekliyor, demet alıyor

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Kanal demet yayınlıyor, ancak süreç tek değer bekliyor
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

Demetlerden oluşan bir kanal ürettiğimizi görebilirsiniz: `['sample1', 'file1.txt']`; ancak süreç tek bir değer olan `val sample_name` beklemektedir. Yürütülen komut, sürecin `[sample3, file3.txt]_output.txt` adında bir dosya oluşturmaya çalıştığını göstermektedir; bu, istenen çıktı değildir.

#### Kodu düzeltin

Bunu düzeltmek için, süreç her iki girdiyi de gerektiriyorsa süreci bir demet kabul edecek şekilde ayarlayabiliriz:

=== "Seçenek 1: Süreçte demet kabul etme"

    === "Sonra"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Düzeltildi: Demet kabul edildi

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Kanal demet yayınlıyor, ancak süreç tek değer bekliyor
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
                val sample_name  // Tek değer bekliyor, demet alıyor

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Kanal demet yayınlıyor, ancak süreç tek değer bekliyor
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "Seçenek 2: İlk öğeyi çıkarma"

    === "Sonra"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Kanal demet yayınlıyor, ancak süreç tek değer bekliyor
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Düzeltildi: İlk öğe çıkarıldı
        }
        ```

    === "Önce"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Kanal demet yayınlıyor, ancak süreç tek değer bekliyor
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

Kanallar için en güçlü hata ayıklama aracı `.view()` operatörüdür. `.view()` ile hata ayıklamaya yardımcı olmak amacıyla kanallarınızın şeklini tüm aşamalarda anlayabilirsiniz.

#### Pipeline'ı çalıştırın

Bunu uygulamada görmek için `bad_channel_shape_viewed.nf` dosyasını çalıştırın:

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

#### Kodu inceleyin

`.view()` operatörünün nasıl kullanıldığını görmek için `bad_channel_shape_viewed.nf` dosyasını inceleyelim:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Kanal demet yayınlıyor, ancak süreç tek değer bekliyor
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Hata ayıklama: Orijinal kanal içeriğini göster
    .map { tuple -> tuple[0] }        // Dönüşüm: İlk öğeyi çıkar
    .view { "After mapping: $it" }    // Hata ayıklama: Dönüştürülmüş kanal içeriğini göster

    PROCESS_FILES(input_ch)
}
```

#### Kodu düzeltin

Kanal içeriğini anlamak için gelecekte aşırı `.view()` işlemleri kullanmaktan kaçınmak adına yorum eklemek faydalıdır:

```groovy title="bad_channel_shape_viewed.nf (with comments)" linenums="16" hl_lines="8 9"
workflow {

    // Kanal demet yayınlıyor, ancak süreç tek değer bekliyor
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

İş akışlarınız karmaşıklık kazandıkça ve kanal yapısı daha az şeffaf hale geldikçe bu daha da önemli hale gelecektir.

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

### Özetle

Geçerli Nextflow sözdizimi ile pek çok kanal yapısı hatası oluşturulabilir. Veri akışını anlayarak, inceleme için `.view()` operatörlerini kullanarak ve beklenmedik demet yapılarını gösteren köşeli parantezler gibi hata mesajı kalıplarını tanıyarak kanal yapısı hatalarını ayıklayabilirsiniz.

### Sırada ne var?

Süreç tanımlarından kaynaklanan hatalar hakkında bilgi edinin.

---

## 3. Süreç Yapısı Hataları

Süreçlerle ilgili karşılaşacağınız hataların çoğu, komutu oluştururken yaptığınız hatalarla veya temel yazılımla ilgili sorunlarla bağlantılı olacaktır. Bununla birlikte, yukarıdaki kanal sorunlarına benzer şekilde, sözdizimi hatası olarak nitelendirilmeyen ancak çalışma zamanında hatalara neden olan süreç tanımı hataları da yapabilirsiniz.

### 3.1. Eksik Çıktı Dosyaları

Süreç yazarken yaygın bir hata, sürecin beklediği ile üretilen arasında uyuşmazlığa yol açan bir şey yapmaktır.

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

#### Kodu inceleyin

Hata mesajı, sürecin `sample3.txt` adında bir çıktı dosyası üretmesini beklediğini; ancak betiğin aslında `sample3_output.txt` oluşturduğunu belirtmektedir. `missing_output.nf` dosyasındaki süreç tanımını inceleyelim:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Beklenen: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Oluşturulan: sample3_output.txt
    """
}
```

`output:` bloğundaki çıktı dosyası adı ile betikteki ad arasında uyuşmazlık olduğunu görmelisiniz. Bu uyuşmazlık sürecin başarısız olmasına neden olur. Bu tür bir hatayla karşılaşırsanız, süreç tanımınız ile çıktı bloğunuz arasındaki çıktıların eşleşip eşleşmediğini kontrol edin.

Sorun hâlâ net değilse, gerçekte oluşturulan çıktı dosyalarını belirlemek için work dizinini kontrol edin:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

Bu örnek için bu, `output:` tanımımızın aksine çıktı dosyası adına `_output` sonekinin eklendiğini bize gösterecektir.

#### Kodu düzeltin

Çıktı dosyası adını tutarlı hale getirerek uyuşmazlığı düzeltin:

=== "Sonra"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Düzeltildi: Betik çıktısıyla eşleştirildi

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
        path "${sample_name}.txt"  // Beklenen: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Oluşturulan: sample3_output.txt
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

Bir diğer hata sınıfı, yazılım sağlamadaki hatalardan kaynaklanır. `missing_software.nf`, sözdizimsel olarak geçerli bir iş akışıdır; ancak kullandığı `cowpy` komutunu sağlamak için bazı harici yazılımlara bağımlıdır.

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

Süreç, belirttiğimiz komuta erişemiyor. Bazen bu, iş akışının `bin` dizininde bir betiğin mevcut olmasına rağmen çalıştırılabilir yapılmamış olmasından kaynaklanır. Diğer durumlarda ise yazılımın iş akışının çalıştığı konteyner veya ortamda yüklü olmamasından kaynaklanır.

#### Kodu inceleyin

O `127` çıkış koduna dikkat edin; tam olarak sorunu size söyler. `missing_software.nf` dosyasını inceleyelim:

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

Burada biraz yanıltıcı davrandık; aslında kodda yanlış bir şey yok. Sadece süreci, söz konusu komuta erişimi olacak şekilde çalıştırmak için gerekli yapılandırmayı belirtmemiz gerekiyor. Bu durumda sürecin bir konteyner tanımı var; bu nedenle tek yapmamız gereken iş akışını Docker etkin olarak çalıştırmaktır.

#### Pipeline'ı çalıştırın

`nextflow.config` dosyasında sizin için bir Docker profili ayarladık; bu nedenle iş akışını şu şekilde çalıştırabilirsiniz:

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

!!! note "Not"

    Nextflow'un konteynerleri nasıl kullandığı hakkında daha fazla bilgi edinmek için [Hello Nextflow](../hello_nextflow/05_hello_containers.md) bölümüne bakın.

### 3.3. Hatalı kaynak yapılandırması

Üretim kullanımında, süreçlerinizde kaynakları yapılandıracaksınız. Örneğin `memory`, sürecinizin kullanabileceği maksimum bellek miktarını tanımlar; süreç bunu aşarsa zamanlayıcı genellikle süreci sonlandırır ve `137` çıkış kodu döndürür. Bunu burada gösteremiyoruz çünkü `local` yürütücüyü kullanıyoruz; ancak `time` ile benzer bir şey gösterebiliriz.

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

#### Kodu inceleyin

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

Sürecin bir saniyeden uzun süreceğini biliyoruz (bunu garantilemek için bir sleep ekledik); ancak süreç 1 milisaniye sonra zaman aşımına uğrayacak şekilde ayarlanmış. Birisi yapılandırmasında biraz gerçekçi olmayan bir değer kullanmış!

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

Hata mesajlarınızı dikkatle okuduğunuzda bu tür hatalar sizi uzun süre şaşırtmamalıdır. Ancak kaynak yönergelerinizi uygun şekilde yapılandırabilmek için çalıştırdığınız komutların kaynak gereksinimlerini anladığınızdan emin olun.

### 3.4. Süreç Hata Ayıklama Teknikleri

Süreçler başarısız olduğunda veya beklenmedik şekilde davrandığında, neyin yanlış gittiğini araştırmak için sistematik tekniklere ihtiyaç duyarsınız. Work dizini, süreç yürütmesini hata ayıklamak için ihtiyacınız olan tüm bilgileri içerir.

#### Work Dizini İncelemesini Kullanma

Süreçler için en güçlü hata ayıklama aracı, work dizinini incelemektir. Bir süreç başarısız olduğunda, Nextflow o belirli süreç yürütmesi için ne olduğunu anlamak amacıyla gereken tüm dosyaları içeren bir work dizini oluşturur.

#### Pipeline'ı çalıştırın

Work dizini incelemesini göstermek için daha önceki `missing_output.nf` örneğini kullanalım (gerekirse bir çıktı adlandırma uyuşmazlığını yeniden oluşturun):

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

#### Work dizinini kontrol edin

Bu hatayı aldığınızda, work dizini tüm hata ayıklama bilgilerini içerir. Hata mesajından work dizini yolunu bulun ve içeriğini inceleyin:

```bash
# Hata mesajından work dizinini bulun
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

Ardından temel dosyaları inceleyebilirsiniz:

##### Komut Betiğini Kontrol Edin

`.command.sh` dosyası tam olarak hangi komutun yürütüldüğünü gösterir:

```bash
# Yürütülen komutu görüntüle
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

Bu şunları ortaya koyar:

- **Değişken genişletme**: Nextflow değişkenlerinin doğru şekilde genişletilip genişletilmediği
- **Dosya yolları**: Girdi dosyalarının doğru konumda olup olmadığı
- **Komut yapısı**: Betik sözdiziminin doğru olup olmadığı

Dikkat edilmesi gereken yaygın sorunlar:

- **Eksik tırnak işaretleri**: Boşluk içeren değişkenlerin uygun şekilde tırnak içine alınması gerekir
- **Yanlış dosya yolları**: Mevcut olmayan veya yanlış konumdaki girdi dosyaları
- **Yanlış değişken adları**: Değişken referanslarındaki yazım hataları
- **Eksik ortam kurulumu**: Belirli ortamlara bağımlı komutlar

##### Hata Çıktısını Kontrol Edin

`.command.err` dosyası gerçek hata mesajlarını içerir:

```bash
# Hata çıktısını görüntüle
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

Bu dosya şunları gösterir:

- **Çıkış kodları**: 127 (komut bulunamadı), 137 (sonlandırıldı) vb.
- **İzin hataları**: Dosya erişim sorunları
- **Yazılım hataları**: Uygulamaya özgü hata mesajları
- **Kaynak hataları**: Bellek/zaman sınırı aşıldı

##### Standart Çıktıyı Kontrol Edin

`.command.out` dosyası komutunuzun ürettiklerini gösterir:

```bash
# Standart çıktıyı görüntüle
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

Bu şunları doğrulamaya yardımcı olur:

- **Beklenen çıktı**: Komutun doğru sonuçları üretip üretmediği
- **Kısmi yürütme**: Komutun başlayıp başlamadığı ancak yarıda başarısız olup olmadığı
- **Hata ayıklama bilgisi**: Betiğinizden gelen tanı çıktısı

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
# Work dizinindeki tüm dosyaları listele
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

Bu şunları belirlemeye yardımcı olur:

- **Dosya adlandırma uyuşmazlıkları**: Beklenenden farklı adlara sahip çıktı dosyaları
- **İzin sorunları**: Oluşturulamayan dosyalar
- **Yol sorunları**: Yanlış dizinlerde oluşturulan dosyalar

Daha önceki örneğimizde bu, beklenen `sample3.txt` dosyasının mevcut olmadığını; ancak `sample3_output.txt` dosyasının var olduğunu bize doğruladı:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### Özetle

Süreç hata ayıklaması, neyin yanlış gittiğini anlamak için work dizinlerini incelemeyi gerektirir. Temel dosyalar arasında `.command.sh` (yürütülen betik), `.command.err` (hata mesajları) ve `.command.out` (standart çıktı) yer almaktadır. 127 (komut bulunamadı) ve 137 (süreç sonlandırıldı) gibi çıkış kodları, hata türü hakkında anında tanı ipuçları sağlar.

### Sırada ne var?

Nextflow'un yerleşik hata ayıklama araçları ve sorun gidermeye yönelik sistematik yaklaşımlar hakkında bilgi edinin.

---

## 4. Yerleşik Hata Ayıklama Araçları ve Gelişmiş Teknikler

Nextflow, iş akışı yürütmesini hata ayıklamak ve analiz etmek için çeşitli güçlü yerleşik araçlar sunar. Bu araçlar, neyin yanlış gittiğini, nerede yanlış gittiğini ve bunu nasıl verimli şekilde düzelteceğinizi anlamanıza yardımcı olur.

### 4.1. Gerçek Zamanlı Süreç Çıktısı

Bazen çalışan süreçlerin içinde neler olduğunu görmeniz gerekir. Her görevin yürütülürken tam olarak ne yaptığını gösteren gerçek zamanlı süreç çıktısını etkinleştirebilirsiniz.

#### Pipeline'ı çalıştırın

Daha önceki örneklerimizden `bad_channel_shape_viewed.nf`, `.view()` kullanarak kanal içeriğini yazdırdı; ancak `bad_channel_shape_viewed_debug.nf` dosyasında gösterdiğimiz gibi, sürecin kendisinden değişkenleri yankılamak için `debug` yönergesini de kullanabiliriz. İş akışını çalıştırın:

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

#### Kodu inceleyin

`debug` yönergesinin nasıl çalıştığını görmek için `bad_channel_shape_viewed_debug.nf` dosyasını inceleyelim:

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

`debug` yönergesi, bir sürecin ortamını anlamak için hızlı ve kullanışlı bir yol olabilir.

### 4.2. Önizleme Modu

Bazen herhangi bir süreç çalışmadan önce sorunları yakalamak istersiniz. Nextflow bu tür proaktif hata ayıklama için bir bayrak sağlar: `-preview`.

#### Pipeline'ı çalıştırın

Önizleme modu, komutları yürütmeden iş akışı mantığını test etmenizi sağlar. Bu, herhangi bir gerçek komut çalıştırmadan iş akışınızın yapısını hızla kontrol etmek ve süreçlerin doğru şekilde bağlandığından emin olmak için oldukça kullanışlı olabilir.

!!! note "Not"

    Daha önce `bad_syntax.nf` dosyasını düzelttiyseniz, bu komutu çalıştırmadan önce betik bloğundan sonra kapanış süslü parantezini kaldırarak sözdizimi hatasını yeniden ekleyin.

Şu komutu çalıştırın:

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

Önizleme modu, herhangi bir süreç çalıştırmadan sözdizimi hatalarını erken yakalamak için özellikle kullanışlıdır. Yürütmeden önce iş akışı yapısını ve süreç bağlantılarını doğrular.

### 4.3. Mantık Testi için Stub Çalıştırma

Bazen hatalar, komutların çok uzun sürmesi, özel yazılım gerektirmesi veya karmaşık nedenlerle başarısız olması nedeniyle hata ayıklaması güç olabilir. Stub çalıştırma, gerçek komutları yürütmeden iş akışı mantığını test etmenizi sağlar.

#### Pipeline'ı çalıştırın

Bir Nextflow süreci geliştirirken, gerçek komutu çalıştırmadan doğru biçimde çıktılar üreten 'sahte' komutlar tanımlamak için `stub` yönergesini kullanabilirsiniz. Bu yaklaşım, gerçek yazılımın karmaşıklıklarıyla uğraşmadan önce iş akışı mantığınızın doğru olduğunu doğrulamak istediğinizde özellikle değerlidir.

Örneğin, daha önceki `missing_software.nf` dosyamızı hatırlıyor musunuz? `-profile docker` ekleyene kadar iş akışının çalışmasını engelleyen eksik yazılımın olduğu dosyayı? `missing_software_with_stub.nf` çok benzer bir iş akışıdır. Aynı şekilde çalıştırırsak aynı hatayı alırız:

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

Ancak bu iş akışı, `docker` profili olmadan bile `-stub-run` ile çalıştırıldığında hata üretmez:

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

#### Kodu inceleyin

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

`missing_software.nf` dosyasına kıyasla, bu süreçte Nextflow'un stub modunda çalıştırıldığında `script:` bölümünde belirtilen komut yerine kullanılacak bir komutu belirten `stub:` yönergesi bulunmaktadır.

Burada kullandığımız `touch` komutu herhangi bir yazılıma veya uygun girdilere bağımlı değildir ve tüm durumlarda çalışır; bu da süreç iç yapısıyla uğraşmadan iş akışı mantığını hata ayıklamamıza olanak tanır.

**Stub çalıştırma şunların hata ayıklamasına yardımcı olur:**

- Kanal yapısı ve veri akışı
- Süreç bağlantıları ve bağımlılıkları
- Parametre yayılımı
- Yazılım bağımlılıkları olmadan iş akışı mantığı

### 4.4. Sistematik Hata Ayıklama Yaklaşımı

İzleme dosyaları ve work dizinlerinden önizleme moduna, stub çalıştırmaya ve kaynak izlemeye kadar bireysel hata ayıklama tekniklerini öğrendiğinize göre, bunları sistematik bir metodoloji halinde bir araya getirelim. Yapılandırılmış bir yaklaşıma sahip olmak, karmaşık hatalar karşısında bunalmaktan sizi korur ve önemli ipuçlarını kaçırmamanızı sağlar.

Bu metodoloji, ele aldığımız tüm araçları verimli bir iş akışında birleştirir:

**Dört Aşamalı Hata Ayıklama Yöntemi:**

**Aşama 1: Sözdizimi Hatası Çözümü (5 dakika)**

1. VSCode veya IDE'nizde kırmızı alt çizgileri kontrol edin
2. Sözdizimi sorunlarını belirlemek için `nextflow run workflow.nf -preview` komutunu çalıştırın
3. Tüm sözdizimi hatalarını düzeltin (eksik süslü parantezler, sondaki virgüller vb.)
4. Devam etmeden önce iş akışının başarıyla ayrıştırıldığından emin olun

**Aşama 2: Hızlı Değerlendirme (5 dakika)**

1. Çalışma zamanı hata mesajlarını dikkatlice okuyun
2. Çalışma zamanı, mantık veya kaynak hatası olup olmadığını kontrol edin
3. Temel iş akışı mantığını test etmek için önizleme modunu kullanın

**Aşama 3: Ayrıntılı Araştırma (15-30 dakika)**

1. Başarısız görevin work dizinini bulun
2. Log dosyalarını inceleyin
3. Kanalları incelemek için `.view()` operatörleri ekleyin
4. Yürütme olmadan iş akışı mantığını test etmek için `-stub-run` kullanın

**Aşama 4: Düzeltme ve Doğrulama (15 dakika)**

1. Minimal ve hedefli düzeltmeler yapın
2. Resume ile test edin: `nextflow run workflow.nf -resume`
3. Tam iş akışı yürütmesini doğrulayın

!!! tip "Verimli Hata Ayıklama için Resume Kullanımı"

    Bir sorunu tespit ettikten sonra, iş akışınızın başarılı bölümlerini yeniden çalıştırarak zaman kaybetmeden düzeltmelerinizi test etmek için verimli bir yola ihtiyacınız vardır. Nextflow'un `-resume` işlevi, hata ayıklama için son derece değerlidir.

    [Hello Nextflow](../hello_nextflow/) üzerinde çalıştıysanız `-resume` ile karşılaşmış olmalısınız; sorunlu sürecinizden önceki süreçlerin çalışmasını beklerken kendinizi bekletmemek için hata ayıklama sırasında bundan iyi yararlanmanız önemlidir.

    **Resume hata ayıklama stratejisi:**

    1. İş akışını başarısız olana kadar çalıştırın
    2. Başarısız görev için work dizinini inceleyin
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

            // Hata ayıklama için muhafazakâr kaynaklar
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

Ardından pipeline'ı bu profil etkin olarak çalıştırabilirsiniz:

```bash
nextflow run workflow.nf -profile debug
```

Bu profil gerçek zamanlı çıktıyı etkinleştirir, work dizinlerini korur ve daha kolay hata ayıklama için paralelleştirmeyi sınırlar.

### 4.5. Pratik Hata Ayıklama Alıştırması

Artık sistematik hata ayıklama yaklaşımını pratiğe dökme zamanı. `buggy_workflow.nf` iş akışı, gerçek dünya geliştirme sürecinde karşılaşacağınız sorun türlerini temsil eden çeşitli yaygın hatalar içermektedir.

!!! exercise "Alıştırma"

    `buggy_workflow.nf` dosyasındaki tüm hataları tespit edip düzeltmek için sistematik hata ayıklama yaklaşımını kullanın. Bu iş akışı, bir CSV dosyasından örnek verileri işlemeye çalışır; ancak yaygın hata ayıklama senaryolarını temsil eden birden fazla kasıtlı hata içerir.

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

        Bu şifreli hata, `params{}` bloğundaki 11-12. satırlar civarında bir ayrıştırma sorununa işaret etmektedir. v2 ayrıştırıcısı yapısal sorunları erken yakalar.

    Öğrendiğiniz dört aşamalı hata ayıklama yöntemini uygulayın:

    **Aşama 1: Sözdizimi Hatası Çözümü**
    - VSCode veya IDE'nizde kırmızı alt çizgileri kontrol edin
    - Sözdizimi sorunlarını belirlemek için `nextflow run workflow.nf -preview` komutunu çalıştırın
    - Tüm sözdizimi hatalarını düzeltin (eksik süslü parantezler, sondaki virgüller vb.)
    - Devam etmeden önce iş akışının başarıyla ayrıştırıldığından emin olun

    **Aşama 2: Hızlı Değerlendirme**
    - Çalışma zamanı hata mesajlarını dikkatlice okuyun
    - Hataların çalışma zamanı, mantık veya kaynakla ilgili olup olmadığını belirleyin
    - Temel iş akışı mantığını test etmek için `-preview` modunu kullanın

    **Aşama 3: Ayrıntılı Araştırma**
    - Başarısız görevler için work dizinlerini inceleyin
    - Kanalları incelemek için `.view()` operatörleri ekleyin
    - Work dizinlerindeki log dosyalarını kontrol edin
    - Yürütme olmadan iş akışı mantığını test etmek için `-stub-run` kullanın

    **Aşama 4: Düzeltme ve Doğrulama**
    - Hedefli düzeltmeler yapın
    - Düzeltmeleri verimli şekilde test etmek için `-resume` kullanın
    - Tam iş akışı yürütmesini doğrulayın

    **Kullanabileceğiniz Hata Ayıklama Araçları:**
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

    ??? solution "Çözüm"
        `buggy_workflow.nf` dosyası, tüm önemli hata ayıklama kategorilerini kapsayan 9 veya 10 farklı hata içermektedir (nasıl saydığınıza bağlı olarak). Her hatanın ve nasıl düzeltileceğinin sistematik bir dökümü aşağıda verilmiştir.

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
        // HATA: processFiles süreci için kapanış süslü parantezi eksik
        ```
        **Düzeltme:** Eksik kapanış süslü parantezini ekleyin
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Eksik kapanış süslü parantezini ekle
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

        Bu noktada iş akışı çalışacak; ancak hâlâ hatalar alacağız (örneğin `processFiles` içinde `Path value cannot be null`), bunun nedeni hatalı kanal yapısıdır.

        **Hata 5: Kanal Yapısı Hatası - Yanlış Map Çıktısı**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // HATA: processFiles demet bekliyor
        ```
        **Düzeltme:** processFiles'ın beklediği demet yapısını döndürün
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        Ancak bu, yukarıdaki `heavyProcess()` çalıştırma düzeltmemizi bozacak; bu nedenle o sürece yalnızca örnek kimliklerini iletmek için bir map kullanmamız gerekecek:

        **Hata 6: heavyProcess için Hatalı Kanal Yapısı**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // HATA: input_ch artık emisyon başına 2 öğeye sahip - heavyProcess yalnızca 1'e (ilkine) ihtiyaç duyuyor
        ```
        **Düzeltme:** Doğru kanalı kullanın ve örnek kimliklerini çıkarın
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Artık biraz daha ilerliyoruz; ancak bir Bash değişkenini kaçırmadığımız için `No such variable: i` hatası alıyoruz.

        **Hata 7: Bash Değişkeni Kaçırma Hatası**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // HATA: $i kaçırılmamış
        ```
        **Düzeltme:** Bash değişkenini kaçırın
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Artık `Process exceeded running time limit (1ms)` hatası alıyoruz; bu nedenle ilgili süreç için çalışma süresi sınırını düzeltiyoruz:

        **Hata 8: Kaynak Yapılandırma Hatası**
        ```groovy linenums="36"
        time '1 ms'  // HATA: Gerçekçi olmayan zaman sınırı
        ```
        **Düzeltme:** Gerçekçi bir zaman sınırına yükseltin
        ```groovy linenums="36"
        time '100 s'
        ```

        Ardından çözülmesi gereken bir `Missing output file(s)` hatası var:

        **Hata 9: Çıktı Dosyası Adı Uyuşmazlığı**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // HATA: Yanlış dosya adı, çıktı bildirimiyle eşleşmeli
        ```
        **Düzeltme:** Çıktı bildirimiyle eşleştirin
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        İlk iki süreç çalıştı; ancak üçüncüsü çalışmadı.

        **Hata 10: Çıktı Dosyası Adı Uyuşmazlığı**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Hata: Bir süreçten değil, mevcut çalışma dizininden girdi almaya çalışıyor
        handleFiles(file_ch)
        ```
        **Düzeltme:** Önceki sürecin çıktısını alın
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        Bununla birlikte tüm iş akışı çalışmalıdır.

        **Tam Düzeltilmiş İş Akışı:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Hata ayıklama alıştırmaları için hatalı iş akışı
        * Bu iş akışı öğrenme amacıyla birkaç kasıtlı hata içermektedir
        */

        params{
            // Doğrulama eksik parametreler
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Girdi/çıktı uyuşmazlığı olan süreç
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

            // Yanlış kullanımlı kanal
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
- **Süreç hataları**: Çıktı dosyası uyuşmazlıkları, değişken kaçırma
- **Kaynak hataları**: Gerçekçi olmayan zaman sınırları

**Temel Hata Ayıklama Dersleri:**

1. **Hata mesajlarını dikkatlice okuyun** - genellikle doğrudan soruna işaret ederler
2. **Sistematik yaklaşımlar kullanın** - bir seferde bir hatayı düzeltin ve `-resume` ile test edin
3. **Veri akışını anlayın** - kanal yapısı hataları genellikle en ince olanlardır
4. **Work dizinlerini kontrol edin** - süreçler başarısız olduğunda, loglar tam olarak neyin yanlış gittiğini söyler

---

## Özet

Bu yan görevde, Nextflow iş akışlarını hata ayıklamak için bir dizi sistematik teknik öğrendiniz.
Bu teknikleri kendi çalışmalarınızda uygulamak, bilgisayarınızla daha az zaman harcamanızı, sorunları daha hızlı çözmenizi ve kendinizi gelecekteki sorunlardan korumanızı sağlayacaktır.

### Temel kalıplar

**1. Sözdizimi hatalarını nasıl tespit edip düzeltirsiniz**:

- Nextflow hata mesajlarını yorumlama ve sorunları bulma
- Yaygın sözdizimi hataları: eksik süslü parantezler, yanlış anahtar sözcükler, tanımsız değişkenler
- Nextflow (Groovy) ve Bash değişkenleri arasındaki farkı ayırt etme
- Erken hata tespiti için VS Code eklenti özelliklerini kullanma

```groovy
// Eksik süslü parantez - IDE'de kırmızı alt çizgileri arayın
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- eksik!

// Yanlış anahtar sözcük
inputs:  // 'input:' olmalı

// Tanımsız değişken - Bash değişkenleri için ters eğik çizgiyle kaçırın
echo "${undefined_var}"      // Nextflow değişkeni (tanımlı değilse hata)
echo "\${bash_var}"          // Bash değişkeni (kaçırılmış)
```

**2. Kanal yapısı sorunlarını nasıl hata ayıklarsınız**:

- Kanal kardinalitesini ve tükenme sorunlarını anlama
- Kanal içeriği yapısı uyuşmazlıklarını hata ayıklama
- Kanal incelemesi için `.view()` operatörlerini kullanma
- Çıktıdaki köşeli parantezler gibi hata kalıplarını tanıma

```groovy
// Kanal içeriğini incele
my_channel.view { "Content: $it" }

// Queue channel'ı value channel'a dönüştür (tükenmeyi önler)
reference_ch = channel.value('ref.fa')
// veya
reference_ch = channel.of('ref.fa').first()
```

**3. Süreç yürütme sorunlarını nasıl giderirsiniz**:

- Eksik çıktı dosyası hatalarını teşhis etme
- Çıkış kodlarını anlama (eksik yazılım için 127, bellek sorunları için 137)
- Work dizinlerini ve komut dosyalarını araştırma
- Kaynakları uygun şekilde yapılandırma

```bash
# Gerçekte neyin yürütüldüğünü kontrol edin
cat work/ab/cdef12/.command.sh

# Hata çıktısını kontrol edin
cat work/ab/cdef12/.command.err

# Çıkış kodu 127 = komut bulunamadı
# Çıkış kodu 137 = sonlandırıldı (bellek/zaman sınırı)
```

**4. Nextflow'un yerleşik hata ayıklama araçlarını nasıl kullanırsınız**:

- Önizleme modundan ve gerçek zamanlı hata ayıklamadan yararlanma
- Mantık testi için stub çalıştırmayı uygulama
- Verimli hata ayıklama döngüleri için resume kullanma
- Dört aşamalı sistematik hata ayıklama metodolojisini takip etme

!!! tip "Hızlı Hata Ayıklama Başvurusu"

    **Sözdizimi hataları mı?** → VSCode uyarılarını kontrol edin, `nextflow run workflow.nf -preview` komutunu çalıştırın

    **Kanal sorunları mı?** → İçeriği incelemek için `.view()` kullanın: `my_channel.view()`

    **Süreç hataları mı?** → Work dizini dosyalarını kontrol edin:

    - `.command.sh` - yürütülen betik
    - `.command.err` - hata mesajları
    - `.exitcode` - çıkış durumu (127 = komut bulunamadı, 137 = sonlandırıldı)

    **Gizemli davranış mı?** → İş akışı mantığını test etmek için `-stub-run` ile çalıştırın

    **Düzeltmeler yaptınız mı?** → Test ederken zaman kazanmak için `-resume` kullanın: `nextflow run workflow.nf -resume`

---

### Ek kaynaklar

- [Nextflow sorun giderme kılavuzu](https://www.nextflow.io/docs/latest/troubleshooting.html): Resmi sorun giderme belgeleri
- [Nextflow kanallarını anlama](https://www.nextflow.io/docs/latest/channel.html): Kanal türleri ve davranışlarına derinlemesine bakış
- [Süreç yönergeleri başvurusu](https://www.nextflow.io/docs/latest/process.html#directives): Mevcut tüm süreç yapılandırma seçenekleri
- [nf-test](https://www.nf-test.com/): Nextflow pipeline'ları için test çerçevesi
- [Nextflow Slack topluluğu](https://www.nextflow.io/slack-invite.html): Topluluktan yardım alın

Üretim iş akışları için şunları göz önünde bulundurun:

- Ölçekte izleme ve hata ayıklama için [Seqera Platform](https://seqera.io/platform/) kurulumu
- Tekrarlanabilir yazılım ortamları için [Wave containers](https://seqera.io/wave/) kullanımı

**Unutmayın:** Etkili hata ayıklama, pratikle gelişen bir beceridir. Burada edindiğiniz sistematik metodoloji ve kapsamlı araç seti, Nextflow geliştirme yolculuğunuz boyunca size iyi hizmet edecektir.

---

## Sırada ne var?

[Yan Görevler menüsüne](../) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
