# Bölüm 3: Çok-örnekli paired-end uygulaması

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Daha önce, her örneğin verilerini bağımsız olarak işleyen örnek başına bir varyant çağırma boru hattı oluşturdunuz.
Bu kursun bu bölümünde, basit iş akışımızı bir üst seviyeye taşıyarak keyfi sayıda örneği işleyebilen güçlü bir toplu otomasyon aracına dönüştüreceğiz.
Bunu yaparken, aynı zamanda yeni çalışmalarda daha yaygın olan paired-end verileri kabul edecek şekilde de güncelleyeceğiz.

??? info "Bu bölümden nasıl başlanır"

    Kursun bu bölümü, [Bölüm 1: Yönteme Genel Bakış](./01_method.md) ve [Bölüm 2: Tek-örnekli uygulama](./02_single-sample.md) bölümlerini tamamladığınızı ve doldurulmuş modül dosyalarına sahip çalışan bir `rnaseq.nf` boru hattınız olduğunu varsayar.

    Bölüm 2'yi tamamlamadıysanız veya bu bölüm için yeni başlamak istiyorsanız, Bölüm 2 çözümünü başlangıç noktanız olarak kullanabilirsiniz.
    Bu komutları `nf4-science/rnaseq/` dizininin içinden çalıştırın:

    ```bash
    cp solutions/part2/rnaseq-2.nf rnaseq.nf
    cp solutions/part2/modules/fastqc.nf modules/
    cp solutions/part2/modules/trim_galore.nf modules/
    cp solutions/part2/modules/hisat2_align.nf modules/
    cp solutions/part2/nextflow.config .
    ```

    Bu size eksiksiz bir tek-örnekli işleme iş akışı verir.
    Başarıyla çalıştığını test edebilirsiniz:

    ```bash
    nextflow run rnaseq.nf -profile test
    ```

## Görev

Kursun bu bölümünde, iş akışını aşağıdakileri yapacak şekilde genişleteceğiz:

1. Örnek bilgilerini bir CSV örnek tablosundan okuma
2. Tüm örnekler üzerinde örnek başına QC, kırpma ve hizalama işlemlerini paralel olarak çalıştırma
3. Tüm QC raporlarını kapsamlı bir MultiQC raporunda toplama

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-03.svg"
</figure>

Bu, [Bölüm 1: Yönteme genel bakış](./01_method.md#2-multi-sample-qc-aggregation) bölümünün ikinci kısmındaki adımları otomatikleştirir; burada bu komutları konteynırlarında manuel olarak çalıştırmıştınız.

## Ders planı

Bunu üç aşamaya ayırdık:

1. **İş akışının birden fazla girdi örneğini kabul etmesini sağlama.**
   Bu, tek bir dosya yolundan CSV örnek tablosuna geçişi, `splitCsv()` ile ayrıştırmayı ve mevcut tüm süreçleri birden fazla örnek üzerinde çalıştırmayı kapsar.
2. **Kapsamlı QC raporu oluşturma ekleme.**
   Bu, örnekler arasında çıktıları toplamak için `collect()` operatörünü tanıtır ve birleşik bir rapor üretmek için bir MultiQC süreci ekler.
3. **Paired-end RNAseq verilerine geçiş.**
   Bu, süreçleri paired-end girdiler için uyarlamayı (demetler kullanarak), paired-end modüller oluşturmayı ve ayrı bir test profili ayarlamayı kapsar.

Bu, [Bölüm 1: Yönteme Genel Bakış](./01_method.md)'ta açıklanan yöntemi uygular (çok-örnekli kullanım durumunu kapsayan ikinci bölüm) ve doğrudan Bölüm 2 tarafından üretilen iş akışı üzerine inşa edilir.

!!! tip "İpucu"

     Doğru çalışma dizininde olduğunuzdan emin olun:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. İş akışının birden fazla girdi örneğini kabul etmesini sağlama

Birden fazla örnek üzerinde çalıştırmak için, girdiyi nasıl yönettiğimizi değiştirmemiz gerekiyor: tek bir dosya yolu sağlamak yerine, örnek bilgilerini bir CSV dosyasından okuyacağız.

`data/` dizininde örnek ID'leri ve FASTQ dosya yollarını içeren bir CSV dosyası sağlıyoruz.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Bu CSV dosyası, sütunları adlandıran bir başlık satırı içerir.

Bunun hala single-end okuma verisi olduğuna dikkat edin.

!!! warning "Uyarı"

    CSV'deki dosya yolları, ortamınızla eşleşmesi gereken mutlak yollardır.
    Bunu sağladığımız eğitim ortamında çalıştırmıyorsanız, yolları sisteminize uyacak şekilde güncellemeniz gerekecektir.

### 1.1. Test profilinde birincil girdiyi dosya yollarının CSV'si olacak şekilde değiştirme

İlk olarak, tek FASTQ yolu yerine CSV dosya yolunu sağlamak için `nextflow.config` dosyasındaki test profilini güncellememiz gerekiyor.

=== "Sonra"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Ardından, bu CSV'den okumak için kanal oluşturmayı güncellememiz gerekecek.

### 1.2. CSV girdisini ayrıştırmak için kanal fabrikasını güncelleme

Dosyanın içeriğini sadece dosya yolu yerine kanala yüklememiz gerekiyor.

Bunu [Hello Nextflow'un Bölüm 2'sinde](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file) kullandığımız aynı kalıbı kullanarak yapabiliriz: dosyayı ayrıştırmak için [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) operatörünü uygulayarak, ardından her satırdan FASTQ dosya yolunu çıkarmak için bir `map` işlemi.

=== "Sonra"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="4-7"
    workflow {

        main:
        // Bir CSV dosyasının içeriğinden girdi kanalı oluşturma
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

=== "Önce"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Bir dosya yolundan girdi kanalı oluşturma
        read_ch = channel.fromPath(params.input)
    ```

Hello Nextflow kursunda karşılaştığınız şeye kıyasla yeni olan bir şey, bu CSV'nin bir başlık satırına sahip olmasıdır, bu yüzden `splitCsv()` çağrısına `#!groovy header: true` ekliyoruz.
Bu, `map` işleminde sütunlara ada göre başvurmamızı sağlar: `#!groovy row.fastq_path`, her satırın `fastq_path` sütunundan dosya yolunu çıkarır.

Girdi işleme güncellendi ve iş akışı test edilmeye hazır.

### 1.3. İş akışını çalıştırma

İş akışı artık örnek bilgilerini bir CSV dosyasından okur ve tüm örnekleri paralel olarak işler.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [golden_curry] DSL2 - revision: 2a5ba5be1e

    executor >  local (18)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6 ✔
    [cc/16859f] TRIM_GALORE (6)  [100%] 6 of 6 ✔
    [68/4c27b5] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
    ```

Bu sefer her adım 6 kez çalıştırılıyor, CSV dosyasındaki her örnek için bir kez.

İş akışını birden fazla dosya üzerinde çalıştırmak için gereken tek şey buydu.
Nextflow tüm paralelliği bizim için yönetiyor.

### Çıkarım

Tek dosyalı girdiden, Nextflow'un paralel olarak işlediği CSV tabanlı çok-örnekli girdiye nasıl geçeceğinizi biliyorsunuz.

### Sırada ne var?

Tüm örneklerden metrikleri birleştiren bir QC raporu toplama adımı ekleyin.

---

## 2. Ön işleme QC metriklerini tek bir MultiQC raporunda toplama

Tüm bunlar çok sayıda QC raporu üretiyor ve bireysel raporları incelemek zorunda kalmak istemiyoruz.
Bu, bir MultiQC raporu toplama adımı eklemek için mükemmel bir nokta.

[Bölüm 1'den](01_method.md) `multiqc` komutunu hatırlayın:

```bash
multiqc . -n <output_name>.html
```

Komut, tanınan QC çıktı dosyaları için geçerli dizini tarar ve bunları tek bir HTML raporunda toplar.
Konteyner URI'si `community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c` idi.

Ek bir parametre ayarlamamız, girdileri hazırlamamız, süreci yazmamız, bağlamamız ve çıktı işlemeyi güncellememiz gerekiyor.

### 2.1. Girdileri ayarlama

MultiQC süreci bir rapor adı parametresine ve önceki tüm adımlardan toplanan QC çıktılarına ihtiyaç duyar.

#### 2.1.1. Bir `report_id` parametresi ekleme

Çıktı raporunu adlandırmak için bir parametre ekleyin.

=== "Sonra"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="8-9"
    params {
        // Birincil girdi
        input: Path

        // Referans genom arşivi
        hisat2_index_zip: Path

        // Rapor ID'si
        report_id: String
    }
    ```

=== "Önce"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Birincil girdi
        input: Path

        // Referans genom arşivi
        hisat2_index_zip: Path
    }
    ```

Rapor ID'si varsayılanını test profiline ekleyin:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="1" hl_lines="7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Ardından, MultiQC süreci için girdileri hazırlamamız gerekecek.

#### 2.1.2. Önceki adımlardan QC çıktılarını toplama ve birleştirme

`MULTIQC` sürecine önceki adımlardan gelen tüm QC ile ilgili çıktıları bir arada vermemiz gerekiyor.

Bunun için, birden fazla kanalı tek bir kanalda toplayan `.mix()` operatörünü kullanıyoruz.
`channel.empty()`'den başlıyoruz ve birleştirmek istediğimiz tüm çıktı kanallarını karıştırıyoruz.
Bu, `.mix()`'i doğrudan çıktı kanallarından birine zincirleme yapmaktan daha temizdir, çünkü tüm girdileri simetrik olarak ele alır.

İş akışımızda, toplanması gereken QC ile ilgili çıktılar şunlardır:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Bunları tek bir kanalda karıştırıyoruz, ardından raporları tüm örnekler arasında tek bir listeye toplamak için `.collect()` kullanıyoruz.

Bu satırları `HISAT2_ALIGN` çağrısından sonra iş akışı gövdesine ekleyin:

=== "Sonra"

    ```groovy title="rnaseq.nf" linenums="38" hl_lines="4-13"
        // Bir referans genoma hizalama
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

        // Kapsamlı QC raporu oluşturma
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
        multiqc_files_list = multiqc_files_ch.collect()
    ```

=== "Önce"

    ```groovy title="rnaseq.nf" linenums="38"
        // Bir referans genoma hizalama
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

Ara değişkenler kullanmak her adımı netleştirir: `multiqc_files_ch` tüm bireysel QC dosyalarını tek bir kanalda karıştırılmış olarak içerir ve `multiqc_files_list`, MultiQC'ye geçmeye hazır toplanan pakettir.

### 2.2. QC toplama sürecini yazma ve iş akışında çağırma

Daha önce olduğu gibi, süreç tanımını doldurmamız, modülü içe aktarmamız ve süreç çağrısını eklememiz gerekiyor.

#### 2.2.1. QC toplama süreci için modülü doldurma

`modules/multiqc.nf` dosyasını açın ve süreç tanımının ana hatlarını inceleyin.

Devam edin ve yukarıda sağlanan bilgileri kullanarak süreç tanımını kendiniz doldurun, ardından çalışmanızı aşağıdaki "Sonra" sekmesindeki çözümle karşılaştırın.

=== "Önce"

    ```groovy title="modules/multiqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * MultiQC ile QC raporlarını toplama
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Sonra"

    ```groovy title="modules/multiqc.nf" linenums="1" hl_lines="8 11 12 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * MultiQC ile QC raporlarını toplama
     */
    process MULTIQC {

        container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"

        input:
        path '*'
        val output_name

        output:
        path "${output_name}.html", emit: report
        path "${output_name}_data", emit: data

        script:
        """
        multiqc . -n ${output_name}.html
        """
    }
    ```

Bu süreç, QC dosyaları için girdi niteleyicisi olarak `#!groovy path '*'` kullanır.
`'*'` joker karakteri, Nextflow'a toplanan tüm dosyaları belirli adlar gerektirmeden çalışma dizinine yerleştirmesini söyler.
`val output_name` girdisi, rapor dosya adını kontrol eden bir dizedir.

`multiqc .` komutu geçerli dizini (tüm yerleştirilmiş QC dosyalarının bulunduğu yer) tarar ve raporu oluşturur.

Bunu tamamladığınızda, süreç kullanıma hazırdır.

#### 2.2.2. Modülü dahil etme

İçe aktarma ifadesini `rnaseq.nf` dosyasına ekleyin:

=== "Sonra"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="5"
    // Modül INCLUDE ifadeleri
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    include { MULTIQC } from './modules/multiqc.nf'
    ```

=== "Önce"

    ```groovy title="rnaseq.nf" linenums="3"
    // Modül INCLUDE ifadeleri
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Şimdi süreç çağrısını iş akışına ekleyin.

#### 2.2.3. Süreç çağrısını ekleme

Toplanan QC dosyalarını ve rapor ID'sini `MULTIQC` sürecine geçirin:

=== "Sonra"

    ```groovy title="rnaseq.nf" linenums="49" hl_lines="2"
        multiqc_files_list = multiqc_files_ch.collect()
        MULTIQC(multiqc_files_list, params.report_id)
    ```

=== "Önce"

    ```groovy title="rnaseq.nf" linenums="49"
        multiqc_files_list = multiqc_files_ch.collect()
    ```

MultiQC süreci artık iş akışına bağlandı.

### 2.3. Çıktı işlemeyi güncelleme

MultiQC çıktılarını yayınlama bildirimine eklememiz ve nereye gideceklerini yapılandırmamız gerekiyor.

#### 2.3.1. MultiQC çıktıları için yayınlama hedefleri ekleme

MultiQC çıktılarını `publish:` bölümüne ekleyin:

=== "Sonra"

    ```groovy title="rnaseq.nf" linenums="52" hl_lines="9-10"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Önce"

    ```groovy title="rnaseq.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

Ardından, Nextflow'a bu çıktıları nereye koyacağını söylememiz gerekecek.

#### 2.3.2. Yeni çıktı hedeflerini yapılandırma

MultiQC hedefleri için `output {}` bloğuna girdiler ekleyin, bunları bir `multiqc/` alt dizinine yayınlayın:

=== "Sonra"

    ```groovy title="rnaseq.nf" linenums="82" hl_lines="4-9"
        align_log {
            path 'align'
        }
        multiqc_report {
            path 'multiqc'
        }
        multiqc_data {
            path 'multiqc'
        }
    }
    ```

=== "Önce"

    ```groovy title="rnaseq.nf" linenums="82"
        align_log {
            path 'align'
        }
    }
    ```

Çıktı yapılandırması tamamlandı.

### 2.4. İş akışını çalıştırma

Önceki işleme adımlarının önbelleğe alınması ve yalnızca yeni MultiQC adımının çalışması için `-resume` kullanıyoruz.

```bash
nextflow run rnaseq.nf -profile test -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [modest_pare] DSL2 - revision: fc724d3b49

    executor >  local (1)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6, cached: 6 ✔
    [2c/8d8e1e] TRIM_GALORE (5)  [100%] 6 of 6, cached: 6 ✔
    [a4/7f9c44] HISAT2_ALIGN (6) [100%] 6 of 6, cached: 6 ✔
    [56/e1f102] MULTIQC          [100%] 1 of 1 ✔
    ```

Önbelleğe alınmış süreç çağrılarından sonra MULTIQC'ye tek bir çağrı eklenmiştir.

MultiQC çıktılarını sonuçlar dizininde bulabilirsiniz.

```bash
tree -L 2 results/multiqc
```

```console title="Çıktı"
results/multiqc
├── all_single-end_data
│   ├── cutadapt_filtered_reads_plot.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Counts.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt
│   ├── fastqc_adapter_content_plot.txt
│   ├── fastqc_overrepresented_sequences_plot.txt
│   ├── fastqc_per_base_n_content_plot.txt
│   ├── fastqc_per_base_sequence_quality_plot.txt
│   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   ├── fastqc_per_sequence_quality_scores_plot.txt
│   ├── fastqc_sequence_counts_plot.txt
│   ├── fastqc_sequence_duplication_levels_plot.txt
│   ├── fastqc_sequence_length_distribution_plot.txt
│   ├── fastqc-status-check-heatmap.txt
│   ├── fastqc_top_overrepresented_sequences_table.txt
│   ├── hisat2_se_plot.txt
│   ├── multiqc_citations.txt
│   ├── multiqc_cutadapt.txt
│   ├── multiqc_data.json
│   ├── multiqc_fastqc.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc_hisat2.txt
│   ├── multiqc.log
│   ├── multiqc_software_versions.txt
│   └── multiqc_sources.txt
└── all_single-end.html
```

Son `all_single-end.html` dosyası, göz atmayı kolaylaştıran tek bir HTML dosyasında pratik bir şekilde paketlenmiş tam toplu rapordur.

### Çıkarım

Birden fazla kanaldan çıktıları nasıl toplayacağınızı, bunları `.mix()` ve `.collect()` ile nasıl paketleyeceğinizi ve bir toplama sürecine nasıl geçireceğinizi biliyorsunuz.

### Sırada ne var?

İş akışını paired-end RNAseq verilerini işleyecek şekilde uyarlayın.

---

## 3. Paired-end RNAseq verilerinin işlenmesini etkinleştirme

Şu anda iş akışımız yalnızca single-end RNAseq verilerini işleyebiliyor.
Paired-end RNAseq verilerini görmek giderek daha yaygın hale geliyor, bu nedenle bunu da işleyebilmek istiyoruz.

İş akışını veri türünden tamamen bağımsız hale getirmek biraz daha gelişmiş Nextflow dil özelliklerini kullanmayı gerektirecektir, bu yüzden bunu burada yapmayacağız, ancak neyin uyarlanması gerektiğini göstermek için paired-end işleme versiyonu yapabiliriz.

### 3.1. İş akışını kopyalama ve girdileri güncelleme

Single-end iş akışı dosyasını kopyalayarak başlıyoruz ve bunu paired-end veriler için güncelliyoruz.

#### 3.1.1. İş akışı dosyasını kopyalama

Paired-end versiyonu için başlangıç noktası olarak kullanmak üzere iş akışı dosyasının bir kopyasını oluşturun.

```bash
cp rnaseq.nf rnaseq_pe.nf
```

Şimdi yeni dosyadaki parametreleri ve girdi işlemeyi güncelleyin.

#### 3.1.2. Paired-end test profili ekleme

`data/` dizininde örnek ID'leri ve paired FASTQ dosya yollarını içeren ikinci bir CSV dosyası sağlıyoruz.

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Bu dosyaya işaret eden ve paired-end rapor ID'si kullanan `nextflow.config` dosyasına bir `test_pe` profili ekleyin.

=== "Sonra"

    ```groovy title="nextflow.config" linenums="1" hl_lines="9-13"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
        test_pe {
            params.input = "${projectDir}/data/paired-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_paired-end"
        }
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

Paired-end veriler için test profili hazır.

#### 3.1.3. Kanal fabrikasını güncelleme

`.map()` operatörünün her iki FASTQ dosya yolunu da alması ve bunları bir liste olarak döndürmesi gerekiyor.

=== "Sonra"

    ```groovy title="rnaseq_pe.nf" linenums="25" hl_lines="4"
        // Bir CSV dosyasının içeriğinden girdi kanalı oluşturma
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
    ```

=== "Önce"

    ```groovy title="rnaseq_pe.nf" linenums="25"
        // Bir CSV dosyasının içeriğinden girdi kanalı oluşturma
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

Girdi işleme paired-end veriler için yapılandırıldı.

### 3.2. FASTQC modülünü paired-end veriler için uyarlama

Paired-end versiyonu oluşturmak için modülü kopyalayın:

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

FASTQC süreci girdisinin değişmesi gerekmiyor — Nextflow iki dosyadan oluşan bir liste aldığında, her ikisini de yerleştirir ve `reads` her iki dosya adına da genişler.
Gereken tek değişiklik çıktı bloğundadır: artık örnek başına iki FastQC raporu aldığımız için, `simpleName` tabanlı kalıplardan joker karakterlere geçiyoruz.

=== "Sonra"

    ```groovy title="modules/fastqc_pe.nf" linenums="10" hl_lines="2 3"
        output:
        path "*_fastqc.zip", emit: zip
        path "*_fastqc.html", emit: html
    ```

=== "Önce"

    ```groovy title="modules/fastqc_pe.nf" linenums="10"
        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html
    ```

Bu, süreci single-end veya paired-end verilerden birini işleyebilecek şekilde genelleştirir.

Paired-end versiyonunu kullanmak için `rnaseq_pe.nf` dosyasındaki içe aktarmayı güncelleyin:

=== "Sonra"

    ```groovy title="rnaseq_pe.nf" linenums="4" hl_lines="1"
    include { FASTQC } from './modules/fastqc_pe.nf'
    ```

=== "Önce"

    ```groovy title="rnaseq_pe.nf" linenums="4"
    include { FASTQC } from './modules/fastqc.nf'
    ```

FASTQC modülü ve içe aktarması paired-end veriler için güncellendi.

### 3.3. TRIM_GALORE modülünü paired-end veriler için uyarlama

Paired-end versiyonu oluşturmak için modülü kopyalayın:

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Bu modül daha önemli değişiklikler gerektirir:

- Girdi tek bir yoldan iki yoldan oluşan bir demete değişir
- Komut `--paired` bayrağını ekler ve her iki okuma dosyasını da alır
- Çıktı, Trim Galore'un paired-end adlandırma kurallarını yansıtacak şekilde değişir ve her okuma dosyası için ayrı FastQC raporları üretir

=== "Sonra"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8" hl_lines="2 5 7 8 12"
        input:
        tuple path(read1), path(read2)

        output:
        tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
        path "*_trimming_report.txt", emit: trimming_reports
        path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
        path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

        script:
        """
        trim_galore --fastqc --paired ${read1} ${read2}
        """
    ```

=== "Önce"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8"
        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    ```

`rnaseq_pe.nf` dosyasındaki içe aktarmayı güncelleyin:

=== "Sonra"

    ```groovy title="rnaseq_pe.nf" linenums="5" hl_lines="1"
    include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
    ```

=== "Önce"

    ```groovy title="rnaseq_pe.nf" linenums="5"
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

TRIM_GALORE modülü ve içe aktarması paired-end veriler için güncellendi.

### 3.4. HISAT2_ALIGN modülünü paired-end veriler için uyarlama

Paired-end versiyonu oluşturmak için modülü kopyalayın:

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Bu modül benzer değişiklikler gerektirir:

- Girdi tek bir yoldan iki yoldan oluşan bir demete değişir
- HISAT2 komutu `-U` (eşleşmemiş) yerine `-1` ve `-2` (eşleşmiş) okuma argümanlarına değişir
- Komuttaki ve çıktı bildirimlerindeki `reads.simpleName`'in tüm kullanımları, artık çiftin belirli bir üyesine başvurduğumuz için `read1.simpleName` ile değişir

=== "Sonra"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8" hl_lines="2 6 7 12 13 14"
        input:
        tuple path(read1), path(read2)
        path index_zip

        output:
        path "${read1.simpleName}.bam", emit: bam
        path "${read1.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
            --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
            samtools view -bS -o ${read1.simpleName}.bam
        """
    ```

=== "Önce"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8"
        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    ```

`rnaseq_pe.nf` dosyasındaki içe aktarmayı güncelleyin:

=== "Sonra"

    ```groovy title="rnaseq_pe.nf" linenums="6" hl_lines="1"
    include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
    ```

=== "Önce"

    ```groovy title="rnaseq_pe.nf" linenums="6"
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

HISAT2_ALIGN modülü ve içe aktarması paired-end veriler için güncellendi.

### 3.5. MultiQC toplamasını paired-end çıktıları için güncelleme

Paired-end `TRIM_GALORE` süreci artık bir yerine iki ayrı FastQC rapor kanalı (`fastqc_reports_1` ve `fastqc_reports_2`) üretiyor.
Her ikisini de dahil etmek için `rnaseq_pe.nf` dosyasındaki `.mix()` bloğunu güncelleyin:

=== "Sonra"

    ```groovy title="rnaseq_pe.nf" linenums="40" hl_lines="5 6"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            HISAT2_ALIGN.out.log,
        )
    ```

=== "Önce"

    ```groovy title="rnaseq_pe.nf" linenums="40"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
    ```

MultiQC toplaması artık her iki paired-end FastQC rapor setini de içeriyor.

### 3.6. Çıktı işlemeyi paired-end çıktıları için güncelleme

`publish:` bölümü ve `output {}` bloğu da paired-end `TRIM_GALORE` sürecinden gelen iki ayrı FastQC rapor kanalını yansıtmalıdır.

`rnaseq_pe.nf` dosyasındaki `publish:` bölümünü güncelleyin:

=== "Sonra"

    ```groovy title="rnaseq_pe.nf" linenums="52" hl_lines="6-7"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc_1 = TRIM_GALORE.out.fastqc_reports_1
        trimming_fastqc_2 = TRIM_GALORE.out.fastqc_reports_2
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Önce"

    ```groovy title="rnaseq_pe.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

`output {}` bloğundaki ilgili girdileri güncelleyin:

=== "Sonra"

    ```groovy title="rnaseq_pe.nf" linenums="77" hl_lines="4-9"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc_1 {
            path 'trimming'
        }
        trimming_fastqc_2 {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

=== "Önce"

    ```groovy title="rnaseq_pe.nf" linenums="77"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

Paired-end iş akışı artık tamamen güncellendi ve çalıştırılmaya hazır.

### 3.7. İş akışını çalıştırma

Bu önbelleğe almayacağı ve işlenecek verilerin öncekinden iki kat daha fazla olduğu için `-resume` kullanmıyoruz, ancak yine de bir dakikadan kısa sürede tamamlanması gerekiyor.

```bash
nextflow run rnaseq_pe.nf -profile test_pe
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq_pe.nf` [reverent_kare] DSL2 - revision: 9c376cc219

    executor >  local (19)
    [c5/cbde15] FASTQC (5)       [100%] 6 of 6 ✔
    [e4/fa2784] TRIM_GALORE (5)  [100%] 6 of 6 ✔
    [3a/e23049] HISAT2_ALIGN (5) [100%] 6 of 6 ✔
    [e6/a3ccd9] MULTIQC          [100%] 1 of 1 ✔
    ```

Şimdi iş akışımızın biri single-end okuma verileri için, diğeri paired-end veriler için olmak üzere biraz farklılaşan iki versiyonuna sahibiz.
Bir sonraki mantıklı adım, iş akışının her iki veri türünü de anında kabul etmesini sağlamak olacaktır; bu, bu kursun kapsamı dışındadır, ancak bunu bir takip kursunda ele alabiliriz.

---

### Çıkarım

Tek örnekli bir iş akışını birden fazla örneğin işlenmesini paralelleştirmek, kapsamlı bir QC raporu oluşturmak ve iş akışını paired-end okuma verilerini kullanacak şekilde uyarlamak için nasıl uyarlayacağınızı biliyorsunuz.

### Sırada ne var?

Kendinize büyük bir alkış verin! Nextflow for RNAseq kursunu tamamladınız.

Öğrendiklerinizi gözden geçirmek ve sırada ne olduğunu öğrenmek için son [kurs özetine](./next_steps.md) gidin.
