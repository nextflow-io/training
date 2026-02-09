# Bölüm 2: Örnek bazında varyant çağırma

Bölüm 1'de Samtools ve GATK komutlarını ilgili konteynırlarında manuel olarak test ettiniz.
Şimdi aynı komutları bir Nextflow iş akışına dönüştüreceğiz.

## Görev

Bu bölümde, aşağıdaki işlemleri yapan bir iş akışı geliştireceğiz:

1. [Samtools](https://www.htslib.org/) kullanarak her BAM girdi dosyası için bir dizin dosyası oluşturma
2. Her BAM girdi dosyası üzerinde GATK HaplotypeCaller'ı çalıştırarak örnek bazında varyant çağrılarını VCF (Variant Call Format) formatında üretme

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Bu, Bölüm 1'deki adımları tekrar üretir; burada bu komutları konteynırlarında manuel olarak çalıştırmıştınız.

Başlangıç noktası olarak, iş akışının ana bölümlerini özetleyen bir `genomics.nf` iş akışı dosyası ve modüllerin yapısını özetleyen iki modül dosyası (samtools_index.nf ve gatk_haplotypecaller.nf) sunuyoruz.
Bu dosyalar işlevsel değildir; amaçları sadece kodun ilginç kısımlarını doldurmanız için bir iskelet görevi görmektir.

## Ders planı

Geliştirme sürecini daha eğitici hale getirmek için bunu dört adıma ayırdık:

1. **Bir BAM dosyası üzerinde Samtools index çalıştıran tek aşamalı bir iş akışı yazma.**
   Bu, bir modül oluşturma, içe aktarma ve bir iş akışında çağırma konularını kapsar.
2. **Dizinlenmiş BAM dosyası üzerinde GATK HaplotypeCaller'ı çalıştıran ikinci bir süreç ekleme.**
   Bu, süreç çıktılarını girdilere zincirleme ve yardımcı dosyaları işleme konularını tanıtır.
3. **İş akışını bir örnek grubu üzerinde çalışacak şekilde uyarlama.**
   Bu, paralel yürütmeyi kapsar ve ilişkili dosyaları bir arada tutmak için demetleri tanıtır.
4. **İş akışının toplu olarak girdi dosyaları içeren bir metin dosyasını kabul etmesini sağlama.**
   Bu, toplu olarak girdi sağlamak için yaygın bir deseni gösterir.

Her adım, iş akışı geliştirmenin belirli bir yönüne odaklanır.

---

## 1. Bir BAM dosyası üzerinde Samtools index çalıştıran tek aşamalı bir iş akışı yazma

Bu ilk adım temellere odaklanır: bir BAM dosyası yükleme ve bunun için bir dizin oluşturma.

[Bölüm 1](01_method.md)'deki `samtools index` komutunu hatırlayın:

```bash
samtools index '<input_bam>'
```

Komut girdi olarak bir BAM dosyası alır ve yanında bir `.bai` dizin dosyası üretir.
Konteyner URI'si `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464` idi.

Bu bilgiyi alıp Nextflow'da üç aşamada sarmalayacağız:

1. Girdiyi ayarlama
2. Dizinleme sürecini yazma ve iş akışında çağırma
3. Çıktı yönetimini yapılandırma

### 1.1. Girdiyi ayarlama

Bir girdi parametresi bildirmemiz, uygun bir varsayılan değer sağlamak için bir test profili oluşturmamız ve bir girdi kanalı oluşturmamız gerekiyor.

#### 1.1.1. Bir girdi parametresi bildirimi ekleme

Ana iş akışı dosyası `genomics.nf`'de, `Pipeline parameters` bölümü altında, `reads_bam` adında bir CLI parametresi bildirin.

=== "Sonra"

    ```groovy title="genomics.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        reads_bam: Path
    }
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Bu, CLI parametresini ayarlar, ancak geliştirme sırasında iş akışını her çalıştırdığımızda dosya yolunu yazmak istemiyoruz.
Varsayılan değer sağlamak için birden fazla seçenek vardır; burada bir test profili kullanıyoruz.

#### 1.1.2. `nextflow.config`'de varsayılan değerli bir test profili oluşturma

Bir test profili, komut satırında girdileri belirtmeden bir iş akışını denemek için uygun varsayılan değerler sağlar.
Bu, Nextflow ekosisteminde yaygın bir kuraldır (daha fazla ayrıntı için [Hello Config](../../hello_nextflow/06_hello_config.md)'e bakın).

`nextflow.config`'e `reads_bam` parametresini test BAM dosyalarından birine ayarlayan bir `test` profili içeren bir `profiles` bloğu ekleyin.

=== "Sonra"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        }
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Burada, iş akışı betiğinin bulunduğu dizine işaret eden yerleşik bir Nextflow değişkeni olan `${projectDir}` kullanıyoruz.
Bu, mutlak yolları sabit kodlamadan veri dosyalarına ve diğer kaynaklara referans vermeyi kolaylaştırır.

#### 1.1.3. Girdi kanalını ayarlama

İş akışı bloğunda, parametre değerinden `.fromPath` kanal fabrikasını kullanarak bir girdi kanalı oluşturun ([Hello Channels](../../hello_nextflow/02_hello_channels.md)'da kullanıldığı gibi).

=== "Sonra"

    ```groovy title="genomics.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="13"
    workflow {

        main:
        // Create input channel
    ```

Şimdi bu girdi üzerinde dizinleme işlemini çalıştıracak süreci oluşturmamız gerekiyor.

### 1.2. Dizinleme sürecini yazma ve iş akışında çağırma

Modül dosyasında süreç tanımını yazmamız, bir include ifadesi kullanarak iş akışına içe aktarmamız ve girdi üzerinde çağırmamız gerekiyor.

#### 1.2.1. Dizinleme süreci için modülü doldurma

`modules/samtools_index.nf` dosyasını açın ve süreç tanımının ana hatlarını inceleyin.
Ana yapısal öğeleri tanımanız gerekir; aksi takdirde, bir hatırlatma için [Hello Nextflow](../../hello_nextflow/01_hello_world.md)'u okumayı düşünün.

Devam edin ve yukarıda sağlanan bilgileri kullanarak süreç tanımını kendiniz doldurun, ardından çalışmanızı aşağıdaki "Sonra" sekmesindeki çözümle kontrol edin.

=== "Önce"

    ```groovy title="modules/samtools_index.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * BAM dizin dosyası oluştur
     */
    process SAMTOOLS_INDEX {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Sonra"

    ```groovy title="modules/samtools_index.nf" linenums="1" hl_lines="8 11 14 18"
    #!/usr/bin/env nextflow

    /*
     * BAM dizin dosyası oluştur
     */
    process SAMTOOLS_INDEX {

        container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

        input:
        path input_bam

        output:
        path "${input_bam}.bai"

        script:
        """
        samtools index '$input_bam'
        """
    }
    ```

Bunu tamamladıktan sonra süreç tamamlanmış olur.
İş akışında kullanmak için modülü içe aktarmanız ve bir süreç çağrısı eklemeniz gerekir.

#### 1.2.2. Modülü dahil etme

`genomics.nf`'de, süreci iş akışına kullanılabilir hale getirmek için bir `include` ifadesi ekleyin:

=== "Sonra"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    ```

Süreç artık iş akışı kapsamında kullanılabilir.

#### 1.2.3. Girdi üzerinde dizinleme sürecini çağırma

Şimdi, iş akışı bloğuna girdi kanalını argüman olarak geçirerek `SAMTOOLS_INDEX`'e bir çağrı ekleyelim.

=== "Sonra"

    ```groovy title="genomics.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="14"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Call processes
    ```

İş akışı artık girdiyi yükler ve üzerinde dizinleme sürecini çalıştırır.
Şimdi, çıktının nasıl yayınlanacağını yapılandırmamız gerekiyor.

### 1.3. Çıktı yönetimini yapılandırma

Hangi süreç çıktılarının yayınlanacağını bildirmemiz ve nereye gideceklerini belirtmemiz gerekiyor.

#### 1.3.1. `publish:` bölümünde bir çıktı bildirme

İş akışı bloğunun içindeki `publish:` bölümü hangi süreç çıktılarının yayınlanması gerektiğini bildirir.
`SAMTOOLS_INDEX` çıktısını `bam_index` adlı bir hedefe atayın.

=== "Sonra"

    ```groovy title="genomics.nf" linenums="22" hl_lines="2"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="22"
        publish:
        // Declare outputs to publish
    }
    ```

Şimdi Nextflow'a yayınlanan çıktıyı nereye koyacağını söylememiz gerekiyor.

#### 1.3.2. `output {}` bloğunda çıktı hedefini yapılandırma

`output {}` bloğu iş akışının dışında bulunur ve her adlandırılmış hedefin nereye yayınlandığını belirtir.
`bam_index` için bir `bam/` alt dizinine yayınlayan bir hedef ekleyelim.

=== "Sonra"

    ```groovy title="genomics.nf" linenums="26" hl_lines="2-4"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="26"
    output {
        // Configure publish targets
    }
    ```

!!! note

    Varsayılan olarak, Nextflow çıktı dosyalarını sembolik bağlantılar olarak yayınlar, bu da gereksiz kopyalamayı önler.
    Burada kullandığımız veri dosyaları çok küçük olsa da, genomik alanında çok büyük olabilirler.
    Sembolik bağlantılar `work` dizininizi temizlediğinizde bozulur, bu nedenle üretim iş akışları için varsayılan yayınlama modunu `'copy'` olarak geçersiz kılmak isteyebilirsiniz.

### 1.4. İş akışını çalıştırma

Bu noktada, tam olarak işlevsel olması gereken tek adımlı bir dizinleme iş akışımız var. Çalıştığını test edelim!

Test profilinde ayarlanan varsayılan değeri kullanmak ve komut satırına yolu yazmak zorunda kalmamak için `-profile test` ile çalıştırabiliriz.

```bash
nextflow run genomics.nf -profile test
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Dizin dosyasının doğru şekilde oluşturulduğunu çalışma dizinine veya sonuçlar dizinine bakarak kontrol edebilirsiniz.

??? abstract "Çalışma dizini içeriği"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Sonuçlar dizini içeriği"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

İşte orada!

### Özet

Bir süreç içeren bir modül oluşturmayı, bir iş akışına içe aktarmayı, bir girdi kanalıyla çağırmayı ve sonuçları yayınlamayı biliyorsunuz.

### Sırada ne var?

Dizinleme sürecinin çıktısını alan ve varyant çağırma işlemini çalıştırmak için kullanan ikinci bir adım ekleyin.

---

## 2. Dizinlenmiş BAM dosyası üzerinde GATK HaplotypeCaller'ı çalıştıran ikinci bir süreç ekleme

Artık girdi dosyamız için bir dizinimiz olduğuna göre, varyant çağırma adımını ayarlamaya geçebiliriz.

[Bölüm 1](01_method.md)'deki `gatk HaplotypeCaller` komutunu hatırlayın:

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

Komut bir BAM dosyası (`-I`), bir referans genom (`-R`) ve bir aralık dosyası (`-L`) alır ve diziniyle birlikte bir VCF dosyası (`-O`) üretir.
Araç ayrıca BAM dizininin, referans dizininin ve referans sözlüğünün ilgili dosyalarıyla birlikte bulunmasını bekler.
Konteyner URI'si `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867` idi.

Daha önce olduğu gibi aynı üç aşamayı izliyoruz:

1. Girdileri ayarlama
2. Varyant çağırma sürecini yazma ve iş akışında çağırma
3. Çıktı yönetimini yapılandırma

### 2.1. Girdileri ayarlama

Varyant çağırma adımı birkaç ek girdi dosyası gerektirir.
Bunlar için parametreler bildirmemiz, test profiline varsayılan değerler eklememiz ve bunları yüklemek için değişkenler oluşturmamız gerekiyor.

#### 2.1.1. Yardımcı girdiler için parametre bildirimleri ekleme

Yeni sürecimiz sağlanması gereken bir avuç ek dosya beklediğinden, `genomics.nf`'de `Pipeline parameters` bölümü altında bunlar için parametre bildirimleri ekleyin:

=== "Sonra"

    ```groovy title="genomics.nf" linenums="9" hl_lines="5-9"
    params {
        // Primary input
        reads_bam: Path

        // Accessory files
        reference: Path
        reference_index: Path
        reference_dict: Path
        intervals: Path
    }
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="9"
    params {
        // Primary input
        reads_bam: Path
    }
    ```

Daha önce olduğu gibi, varsayılan değerleri satır içi yerine test profili aracılığıyla sağlıyoruz.

#### 2.1.2. Test profiline yardımcı dosya varsayılanlarını ekleme

Bölüm 1.1.2'de `reads_bam` için yaptığımız gibi, `nextflow.config`'deki test profiline yardımcı dosyalar için varsayılan değerler ekleyin:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
    }
    ```

Şimdi bu dosya yollarını iş akışında kullanmak için yükleyen değişkenler oluşturmamız gerekiyor.

#### 2.1.3. Yardımcı dosyalar için değişkenler oluşturma

İş akışı bloğunun içine yardımcı dosya yolları için değişkenler ekleyin:

=== "Sonra"

    ```groovy title="genomics.nf" linenums="21" hl_lines="7-11"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Load the file paths for the accessory files (reference and intervals)
        ref_file        = file(params.reference)
        ref_index_file  = file(params.reference_index)
        ref_dict_file   = file(params.reference_dict)
        intervals_file  = file(params.intervals)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="21"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

`file()` sözdizimi, Nextflow'a bu girdileri açıkça dosya yolları olarak işlemesini söyler.
Bu konuda daha fazla bilgiyi [Working with files](../../side_quests/working_with_files.md) Yan Görevinde bulabilirsiniz.

### 2.2. Varyant çağırma sürecini yazma ve iş akışında çağırma

Modül dosyasında süreç tanımını yazmamız, bir include ifadesi kullanarak iş akışına içe aktarmamız ve girdi okumaları artı dizinleme adımının çıktısı ve yardımcı dosyalar üzerinde çağırmamız gerekiyor.

#### 2.2.1. Varyant çağırma süreci için modülü doldurma

`modules/gatk_haplotypecaller.nf` dosyasını açın ve süreç tanımının ana hatlarını inceleyin.

Yukarıda sağlanan bilgileri kullanarak süreç tanımını kendiniz doldurun, ardından çalışmanızı aşağıdaki "Sonra" sekmesindeki çözümle kontrol edin.

=== "Önce"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * GATK HaplotypeCaller ile varyant çağır
     */
    process GATK_HAPLOTYPECALLER {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Sonra"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1" hl_lines="8 11-16 19-20 24-28"
    #!/usr/bin/env nextflow

    /*
     * GATK HaplotypeCaller ile varyant çağır
     */
    process GATK_HAPLOTYPECALLER {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx

        script:
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    }
    ```

Bu sürecin GATK komutunun kendisinin gerektirdiğinden daha fazla girdiye sahip olduğunu fark edeceksiniz.
GATK, adlandırma kurallarına göre BAM dizin dosyasını ve referans genomun yardımcı dosyalarını aramayı bilir, ancak Nextflow alan-agnostiktir ve bu kuralları bilmez.
Nextflow'un bunları çalışma dizininde çalışma zamanında aşamalandırması için açıkça listelememiz gerekir; aksi takdirde GATK eksik dosyalar hakkında bir hata verecektir.

Benzer şekilde, Nextflow'un sonraki adımlar için takip etmesi amacıyla çıktı VCF'nin dizin dosyasını (`"${input_bam}.vcf.idx"`) açıkça listeliyoruz.
Her çıktı kanalına bir ad atamak için `emit:` sözdizimini kullanıyoruz; bu, çıktıları yayınlama bloğuna bağladığımızda faydalı olacaktır.

Bunu tamamladıktan sonra süreç tamamlanmış olur.
İş akışında kullanmak için modülü içe aktarmanız ve bir süreç çağrısı eklemeniz gerekir.

#### 2.2.2. Yeni modülü içe aktarma

`genomics.nf`'yi yeni modülü içe aktarmak için güncelleyin:

=== "Sonra"

    ```groovy title="genomics.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

Süreç artık iş akışı kapsamında kullanılabilir.

#### 2.2.3. Süreç çağrısını ekleme

İş akışı gövdesine, `main:` altına süreç çağrısını ekleyin:

=== "Sonra"

    ```groovy title="genomics.nf" linenums="33" hl_lines="4-12"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
            ref_file,
            ref_index_file,
            ref_dict_file,
            intervals_file
        )
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="33"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

Hello Nextflow eğitim serisinden `*.out` sözdizimini tanımanız gerekir; Nextflow'a `SAMTOOLS_INDEX` tarafından çıkarılan kanalı almasını ve bunu `GATK_HAPLOTYPECALLER` süreç çağrısına bağlamasını söylüyoruz.

!!! note

    Girdilerin, sürecin girdi bloğunda listelendikleri sırayla tamamen aynı sırada süreç çağrısında sağlandığına dikkat edin.
    Nextflow'da girdiler konumsaldır, yani aynı sırayı _takip etmelisiniz_; ve tabii ki aynı sayıda eleman olmalıdır.

### 2.3. Çıktı yönetimini yapılandırma

Yeni çıktıları yayınlama bildirimine eklememiz ve nereye gideceklerini yapılandırmamız gerekiyor.

#### 2.3.1. Varyant çağırma çıktıları için yayınlama hedefleri ekleme

`publish:` bölümüne VCF ve dizin çıktılarını ekleyin:

=== "Sonra"

    ```groovy title="genomics.nf" linenums="45" hl_lines="3-4"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="45"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

Şimdi Nextflow'a yeni çıktıları nereye koyacağını söylememiz gerekiyor.

#### 2.3.2. Yeni çıktı hedeflerini yapılandırma

`output {}` bloğuna `vcf` ve `vcf_idx` hedefleri için girdiler ekleyin, her ikisini de bir `vcf/` alt dizinine yayınlayın:

=== "Sonra"

    ```groovy title="genomics.nf" linenums="51" hl_lines="5-10"
    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="49"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

VCF ve dizini, her ikisi de `vcf/` alt dizinine giden ayrı hedefler olarak yayınlanır.

### 2.4. İş akışını çalıştırma

Genişletilmiş iş akışını çalıştırın, bu sefer dizinleme adımını tekrar çalıştırmak zorunda kalmamak için `-resume` ekleyin.

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Şimdi konsol çıktısına bakarsak, iki sürecin listelendiğini görüyoruz.

İlk süreç beklendiği gibi önbellekleme sayesinde atlandı, ikinci süreç ise yepyeni olduğu için çalıştırıldı.

Çıktı dosyalarını sonuçlar dizininde bulacaksınız (çalışma dizinine sembolik bağlantılar olarak).

??? abstract "Dizin içeriği"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

VCF dosyasını açarsanız, GATK komutunu doğrudan konteynırda çalıştırarak oluşturduğunuz dosyadaki içerikle aynı içeriği görmelisiniz.

??? abstract "Dosya içeriği"

    ```console title="reads_mother.bam.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Bu, çalışmamızdaki her örnek için üretmeyi önemsediğimiz çıktıdır.

### Özet

Gerçek analiz çalışması yapan ve yardımcı dosyalar gibi genomik dosya formatı özelliklerini ele alabilen iki adımlı modüler bir iş akışını nasıl yapacağınızı biliyorsunuz.

### Sırada ne var?

İş akışının birden fazla örneği toplu olarak işlemesini sağlayın.

---

## 3. İş akışını bir örnek grubu üzerinde çalışacak şekilde uyarlama

Tek bir örnek üzerinde işlemi otomatikleştirebilen bir iş akışına sahip olmak güzel, ama ya 1000 örneğiniz varsa?
Tüm örneklerinizi döngüye alan bir bash betiği yazmanız mı gerekiyor?

Hayır, şükürler olsun! Sadece kodda küçük bir değişiklik yapın ve Nextflow bunu da sizin için halledecek.

### 3.1. Girdiyi üç örneği listeleyecek şekilde güncelleme

Birden fazla örnek üzerinde çalıştırmak için test profilini tek bir dosya yolu yerine bir dosya yolları dizisi sağlayacak şekilde güncelleyin.
Bu, çok örnekli yürütmeyi test etmenin hızlı bir yoludur; bir sonraki adımda bir girdi dosyası kullanarak daha ölçeklenebilir bir yaklaşıma geçeceğiz.

İlk olarak, diziler tür açıklamaları kullanamayacağından, parametre bildirimindeki tür açıklamasını yorumlayın:

=== "Sonra"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (array of three samples)
        reads_bam //: Path
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input
        reads_bam: Path
    ```

Ardından test profilini üç örneğin hepsini listeleyecek şekilde güncelleyin:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2-6"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

İş akışı gövdesindeki kanal fabrikası (`.fromPath`) tek bir dosya yolu kadar birden fazla dosya yolunu da kabul eder, bu nedenle başka değişiklik gerekmez.

### 3.2. İş akışını çalıştırma

Artık tesisat üç test örneğinin tamamı üzerinde çalışacak şekilde ayarlandığına göre iş akışını çalıştırmayı deneyin.

```bash
nextflow run genomics.nf -profile test -resume
```

İlginç olan: bu _çalışabilir_ VEYA _başarısız olabilir_. Örneğin, işte başarılı olan bir çalıştırma:

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

İş akışı çalıştırmanız başarılı olduysa, şu gibi bir hata alana kadar tekrar çalıştırın:

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

    executor >  local (4)
    [01/eea165] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 1 ✔
    [a5/fa9fd0] GATK_HAPLOTYPECALLER (3) | 1 of 3, cached: 1
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Process `GATK_HAPLOTYPECALLER (2)` terminated with an error exit status (2)

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_father.bam         -O reads_father.bam.vcf         -L intervals.bed

    Command exit status:
      2

    Command error:
      ...
      A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
      ...
    ```

GATK komut hata çıktısına bakarsanız, şu gibi bir satır olacaktır:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Bu garip, çünkü iş akışının ilk adımında BAM dosyalarını açıkça dizinledik. Tesisatta bir sorun olabilir mi?

### 3.3. Sorunu giderme

Neyin yanlış gittiğini anlamak için çalışma dizinlerini inceleyeceğiz ve `view()` operatörünü kullanacağız.

#### 3.3.1. İlgili çağrılar için çalışma dizinlerini kontrol etme

Konsol çıktısında listelenen başarısız `GATK_HAPLOTYPECALLER` süreç çağrısının çalışma dizinine bakın.

??? abstract "Dizin içeriği"

    ```console
    work/a5/fa9fd0994b6beede5fb9ea073596c2
    ├── intervals.bed -> /workspaces/training/nf4-science/genomics/data/ref/intervals.bed
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/01/eea16597bd6e810fb4cf89e60f8c2d/reads_father.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    ├── reads_son.bam.vcf
    ├── reads_son.bam.vcf.idx
    ├── ref.dict -> /workspaces/training/nf4-science/genomics/data/ref/ref.dict
    ├── ref.fasta -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta
    └── ref.fasta.fai -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta.fai
    ```

Bu dizinde listelenen BAM dosyasının ve BAM dizininin adlarına özellikle dikkat edin: `reads_son.bam` ve `reads_father.bam.bai`.

Ne? Nextflow bu süreç çağrısının çalışma dizininde bir dizin dosyası aşamalandırmış, ancak yanlış olanı. Bu nasıl oldu?

#### 3.3.2. Kanal içeriğini incelemek için [view() operatörünü](https://www.nextflow.io/docs/latest/reference/operator.html#view) kullanma

Kanalın içeriğini görüntülemek için `GATK_HAPLOTYPECALLER` süreç çağrısından önce iş akışı gövdesine şu iki satırı ekleyin:

=== "Sonra"

    ```groovy title="genomics.nf" hl_lines="3-5"
        SAMTOOLS_INDEX(reads_ch)

        // temporary diagnostics
        reads_ch.view()
        SAMTOOLS_INDEX.out.view()

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

=== "Önce"

    ```groovy title="genomics.nf"
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

Ardından iş akışı komutunu tekrar çalıştırın.

```bash
nextflow run genomics.nf -profile test
```

Bir kez daha, bu başarılı olabilir veya başarısız olabilir. İşte başarısız bir çalıştırma için iki `.view()` çağrısının çıktısının neye benzediği:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

İlk üç satır girdi kanalına, ikincisi çıktı kanalına karşılık gelir.
Üç örnek için BAM dosyalarının ve dizin dosyalarının aynı sırada listelenmediğini görebilirsiniz!

!!! note

    Birden fazla öğe içeren bir kanal üzerinde bir Nextflow süreci çağırdığınızda, Nextflow yürütmeyi mümkün olduğunca paralelleştirmeye çalışacak ve çıktıları kullanılabilir oldukları sırayla toplayacaktır.
    Sonuç olarak, karşılık gelen çıktılar orijinal girdilerin verildiği farklı bir sırada toplanabilir.

Şu anda yazıldığı gibi, iş akışı betiğimiz dizin dosyalarının dizinleme adımından girdilerin verildiği anne/baba/oğul sırasıyla aynı sırada listelenerek çıkacağını varsayıyor.
Ancak bunun böyle olacağı garanti değildir, bu yüzden bazen (her zaman olmasa da) yanlış dosyalar ikinci adımda eşleştirilir.

Bunu düzeltmek için BAM dosyalarının ve dizin dosyalarının kanallar aracılığıyla birlikte seyahat etmesini sağlamamız gerekiyor.

!!! tip

    İş akışı kodundaki `view()` ifadeleri hiçbir şey yapmaz, bu nedenle bunları bırakmak sorun değildir.
    Ancak konsol çıktınızı karmaşıklaştıracaklar, bu nedenle sorunu gidermeyi bitirdiğinizde bunları kaldırmanızı öneririz.

### 3.4. İş akışını dizin dosyalarını doğru şekilde işleyecek şekilde güncelleme

Düzeltme, her BAM dosyasını diziniyle birlikte bir demete paketlemek, ardından akış aşağısındaki süreci ve iş akışı tesisatını buna uyacak şekilde güncellemektir.

#### 3.4.1. SAMTOOLS_INDEX modülünün çıktısını bir demete değiştirme

Bir BAM dosyasının ve dizininin yakından ilişkili kalmasını sağlamanın en basit yolu, bunları dizin görevinden çıkan bir demete birlikte paketlemektir.

!!! note

    **Demet**, bir fonksiyondan birden fazla değer döndürmek için yaygın olarak kullanılan sonlu, sıralı bir öğe listesidir. Demetler, birden fazla girdi veya çıktıyı ilişkilerini ve sıralarını koruyarak süreçler arasında geçirmek için özellikle yararlıdır.

BAM dosyasını içerecek şekilde `modules/samtools_index.nf`'deki çıktıyı güncelleyin:

=== "Sonra"

    ```groovy title="modules/samtools_index.nf" linenums="14" hl_lines="2"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Önce"

    ```groovy title="modules/samtools_index.nf" linenums="14"
        output:
        path "${input_bam}.bai"
    ```

Bu şekilde, her dizin dosyası orijinal BAM dosyasıyla sıkı bir şekilde eşleştirilecek ve dizinleme adımının genel çıktısı dosya çiftleri içeren tek bir kanal olacaktır.

#### 3.4.2. GATK_HAPLOTYPECALLER modülünün girdisini bir demet kabul edecek şekilde değiştirme

İlk sürecin çıktısının 'şeklini' değiştirdiğimiz için, ikinci sürecin girdi tanımını buna uyacak şekilde güncellememiz gerekiyor.

`modules/gatk_haplotypecaller.nf`'yi güncelleyin:

=== "Sonra"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11" hl_lines="2"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Önce"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11"
        input:
        path input_bam
        path input_bam_index
    ```

Şimdi süreç çağrısında ve yayınlama hedeflerinde yeni demet yapısını yansıtmak için iş akışını güncellememiz gerekiyor.

#### 3.4.3. İş akışında GATK_HAPLOTYPECALLER çağrısını güncelleme

BAM dosyası artık `SAMTOOLS_INDEX` tarafından çıktı kanalına paketlendiği için, artık `GATK_HAPLOTYPECALLER` sürecine orijinal `reads_ch`'yi sağlamamıza gerek yok.

`genomics.nf`'deki çağrıyı güncelleyin:

=== "Sonra"

    ```groovy title="genomics.nf" linenums="42" hl_lines="2"
        GATK_HAPLOTYPECALLER(
            SAMTOOLS_INDEX.out,
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="42"
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
    ```

Son olarak, yeni çıktı yapısını yansıtmak için yayınlama hedeflerini güncellememiz gerekiyor.

#### 3.4.4. Dizinlenmiş BAM çıktısı için yayınlama hedefini güncelleme

SAMTOOLS_INDEX çıktısı artık hem BAM dosyasını hem de dizinini içeren bir demet olduğundan, yayınlama hedefini `bam_index`'ten `indexed_bam`'e yeniden adlandırarak içeriğini daha iyi yansıtın:

=== "Sonra"

    ```groovy title="genomics.nf" linenums="46" hl_lines="2 8"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="46"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

Bu değişikliklerle, BAM ve dizini birlikte seyahat etmeleri garanti edilir, bu nedenle eşleştirme her zaman doğru olacaktır.

### 3.5. Düzeltilmiş iş akışını çalıştırma

Bunun ileriye dönük güvenilir şekilde çalışacağından emin olmak için iş akışını tekrar çalıştırın.

```bash
nextflow run genomics.nf -profile test
```

Bu sefer (ve her seferinde) her şey doğru şekilde çalışmalıdır:

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Sonuçlar dizini artık her örnek için (demetten) hem BAM hem de BAI dosyalarını, VCF çıktılarıyla birlikte içerir:

??? abstract "Sonuçlar dizini içeriği"

    ```console
    results/
    ├── bam/
    │   ├── reads_father.bam -> ...
    │   ├── reads_father.bam.bai -> ...
    │   ├── reads_mother.bam -> ...
    │   ├── reads_mother.bam.bai -> ...
    │   ├── reads_son.bam -> ...
    │   └── reads_son.bam.bai -> ...
    └── vcf/
        ├── reads_father.bam.vcf -> ...
        ├── reads_father.bam.vcf.idx -> ...
        ├── reads_mother.bam.vcf -> ...
        ├── reads_mother.bam.vcf.idx -> ...
        ├── reads_son.bam.vcf -> ...
        └── reads_son.bam.vcf.idx -> ...
    ```

İlişkili dosyaları demetlere paketleyerek, doğru dosyaların her zaman iş akışı boyunca birlikte seyahat etmesini sağladık.
İş akışı artık herhangi bir sayıda örneği güvenilir bir şekilde işliyor, ancak bunları yapılandırmada tek tek listelemek çok ölçeklenebilir değil.
Bir sonraki adımda, girdileri bir dosyadan okumaya geçeceğiz.

### Özet

İş akışınızın birden fazla örnek üzerinde (bağımsız olarak) çalışmasını nasıl sağlayacağınızı biliyorsunuz.

### Sırada ne var?

Örnekleri toplu olarak işlemeyi kolaylaştırın.

---

## 4. İş akışının toplu girdi dosyaları içeren bir metin dosyasını kabul etmesini sağlama

Bir iş akışına birden fazla veri girdi dosyası sağlamanın çok yaygın bir yolu, dosya yollarını içeren bir metin dosyasıyla yapmaktır.
Satır başına bir dosya yolu listeleyen basit bir metin dosyası kadar basit olabilir veya dosya ek meta veriler içerebilir, bu durumda genellikle örnek listesi olarak adlandırılır.

Burada size basit durumu nasıl yapacağınızı göstereceğiz.

### 4.1. Sağlanan girdi dosya yollarını listeleyen metin dosyasını inceleme

Zaten `data/` dizininde bulabileceğiniz `sample_bams.txt` adlı girdi dosya yollarını listeleyen bir metin dosyası yaptık.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Gördüğünüz gibi, satır başına bir dosya yolu listeledik ve bunlar mutlak yollar.

!!! note

    Burada kullandığımız dosyalar sadece GitHub Codespaces'inizin yerel dosya sistemindedir, ancak bulut depolamasındaki dosyalara da işaret edebiliriz.
    Sağlanan Codespaces ortamını kullanmıyorsanız, dosya yollarını yerel kurulumunuza uyacak şekilde uyarlamanız gerekebilir.

### 4.2. Parametreyi ve test profilini güncelleme

`reads_bam` parametresini tek tek örnekleri listelemek yerine `sample_bams.txt` dosyasına işaret edecek şekilde değiştirin.

Params bloğundaki tür açıklamasını geri yükleyin (çünkü tekrar tek bir yoldur):

=== "Sonra"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (file of input files, one per line)
        reads_bam: Path
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input (array of three samples)
        reads_bam
    ```

Ardından test profilini metin dosyasına işaret edecek şekilde güncelleyin:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2"
    test {
        params.reads_bam = "${projectDir}/data/sample_bams.txt"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

Dosya listesi artık kodda hiç yaşamıyor, bu doğru yönde büyük bir adım.

### 4.3. Kanal fabrikasını bir dosyadan satırları okuyacak şekilde güncelleme

Şu anda, girdi kanal fabrikamız kendisine verdiğimiz tüm dosyaları dizinleme sürecine beslemek istediğimiz veri girdileri olarak ele alıyor.
Ona artık dosya yollarını listeleyen bir dosya verdiğimize göre, davranışını dosyayı ayrıştıracak ve içerdiği dosya yollarını veri girdileri olarak ele alacak şekilde değiştirmemiz gerekiyor.

Bunu [Hello Nextflow'un Bölüm 2](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file)'sinde kullandığımız aynı deseni kullanarak yapabiliriz: dosyayı ayrıştırmak için [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) operatörünü uygulama, ardından her satırın ilk alanını seçmek için bir `map` işlemi.

=== "Sonra"

    ```groovy title="genomics.nf" linenums="24" hl_lines="1-4"
        // Create input channel from a CSV file listing input file paths
        reads_ch = Channel.fromPath(params.reads_bam)
                .splitCsv()
                .map { line -> file(line[0]) }
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="24"
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

Teknik olarak bunu [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) operatörünü kullanarak daha basit bir şekilde yapabilirdik, çünkü girdi dosyamız şu anda yalnızca dosya yolları içeriyor.
Ancak, daha çok yönlü `splitCsv` operatörünü (`map` ile desteklenen) kullanarak, dosya yollarını içeren dosyaya meta veri eklemeye karar verme durumunda iş akışımızı geleceğe hazır hale getirebiliriz.

!!! tip

    Operatörlerin burada ne yaptığını anladığınızdan emin değilseniz, bu, kanal içeriklerinin bunları uygulamadan önce ve sonra neye benzediğine bakmak için `.view()` operatörünü kullanmak için başka bir harika fırsattır.

### 4.4. İş akışını çalıştırma

İş akışını bir kez daha çalıştırın. Bu, daha önce olduğu gibi aynı sonucu üretmeli, değil mi?

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

Evet! Aslında, Nextflow süreç çağrılarının tamamen aynı olduğunu doğru bir şekilde algılıyor ve `-resume` ile çalıştırdığımız için her şeyi yeniden çalıştırmaya bile zahmet etmiyor.

Ve işte bu kadar! Basit varyant çağırma iş akışımız istediğimiz tüm temel özelliklere sahip.

### Özet

Bir BAM dosyasını dizinlemek ve GATK kullanarak örnek bazında varyant çağırma uygulamak için çok adımlı modüler bir iş akışının nasıl yapılacağını biliyorsunuz.

Daha genel olarak, gerçek iş yapan basit bir genomik boru hattı oluşturmak için temel Nextflow bileşenlerini ve mantığını nasıl kullanacağınızı öğrendiniz, genomik dosya formatlarının özelliklerini ve araç gereksinimlerini dikkate alarak.

### Sırada ne var?

Başarınızı kutlayın ve ekstra uzun bir mola verin!

Bu kursun bir sonraki bölümünde, bu basit örnek bazında varyant çağırma iş akışını verilere ortak varyant çağırma uygulamak için nasıl dönüştüreceğinizi öğreneceksiniz.
