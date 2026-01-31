# Bölüm 3: Çok-örnekli paired-end uygulaması

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu kursun son bölümünde, basit iş akışımızı bir üst seviyeye taşıyarak keyfi sayıda örneği işleyebilen güçlü bir toplu otomasyon aracına dönüştüreceğiz.
Bunu yaparken, aynı zamanda yeni çalışmalarda daha yaygın olan paired-end verileri kabul edecek şekilde de değiştireceğiz.

Bunu üç aşamada yapacağız:

1. İş akışının birden fazla girdi örneğini kabul etmesini ve yürütmeyi paralelleştirmesini sağlama
2. Kapsamlı QC raporu oluşturma ekleme
3. Paired-end RNAseq verilerine geçiş

---

## 1. İş akışının birden fazla girdi örneğini kabul etmesini ve yürütmeyi paralelleştirmesini sağlama

Girdileri nasıl yönettiğimizi değiştirmemiz gerekecek.

### 1.1. Birincil girdiyi tek bir dosya yerine dosya yollarının CSV'si olacak şekilde değiştirme

`data/` dizininde örnek ID'leri ve FASTQ dosya yollarını içeren bir CSV dosyası sağlıyoruz.
Bu CSV dosyası bir başlık satırı içerir.
FASTQ dosya yollarının mutlak yollar olduğuna dikkat edin.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Birincil girdi parametresini `input_csv` olarak yeniden adlandıralım ve varsayılanı `single-end.csv` dosyasının yolu olacak şekilde değiştirelim.

```groovy title="rnaseq.nf" linenums="13"
params {
    // Birincil girdi
    input_csv: Path = "data/single-end.csv"

    // Referans genom arşivi
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 1.2. Girdi olarak bir CSV'yi işlemek için girdi kanalı fabrikasını güncelleme

Dosyanın içeriğini sadece dosya yolu yerine kanala yüklemek isteyeceğiz, bu nedenle CSV formatını ayrıştırmak için `.splitCsv()` operatörünü, ardından istediğimiz belirli bilgi parçasını (FASTQ dosya yolunu) almak için `.map()` operatörünü kullanıyoruz.

```groovy title="rnaseq.nf" linenums="16"
    // Bir CSV dosyasının içeriğinden girdi kanalı oluşturma
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }
```

### 1.3. İş akışının çalıştığını test etmek için çalıştırma

```bash
nextflow run rnaseq.nf
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

Bu sefer her adımın sağladığımız 6 veri dosyasının her birinde 6 kez çalıştırıldığını görüyoruz.

İş akışını birden fazla dosya üzerinde çalıştırmak için gereken tek şey buydu!
Nextflow tüm paralelliği bizim için yönetiyor.

---

## 2. Ön işleme QC metriklerini tek bir MultiQC raporunda toplama

Tüm bunlar çok sayıda QC raporu üretiyor ve bireysel raporları incelemek zorunda kalmak istemiyoruz.
Bu, bir MultiQC raporu toplama adımı eklemek için mükemmel bir nokta!

### 2.1. QC toplama işlemi için bir modül oluşturma

`MULTIQC` işlemini barındırmak için `modules/multiqc.nf` adlı bir modül dosyası oluşturalım:

```bash
touch modules/multiqc.nf
```

Dosyayı kod düzenleyicide açın ve aşağıdaki kodu içine kopyalayın:

```groovy title="modules/multiqc.nf" linenums="1"
#!/usr/bin/env nextflow

process MULTIQC {

    container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"
    publishDir "results/multiqc", mode: 'symlink'

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

### 2.2. Modülü iş akışı dosyasına aktarma

`rnaseq.nf` dosyasına `include { MULTIQC } from './modules/multiqc.nf'` ifadesini ekleyin:

```groovy title="rnaseq.nf" linenums="3"
// Modül INCLUDE ifadeleri
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
```

### 2.3. Bir `report_id` parametresi ekleme ve mantıklı bir varsayılan değer verme

```groovy title="rnaseq.nf" linenums="9"
params {
    // Birincil girdi
    input_csv: Path = "data/single-end.csv"

    // Referans genom arşivi
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Rapor ID'si
    report_id: String = "all_single-end"
}
```

### 2.4. İşlemi önceki adımların çıktıları üzerinde çağırma

`MULTIQC` işlemine önceki adımlardan gelen tüm QC ile ilgili çıktıları vermemiz gerekiyor.

Bunun için, birden fazla kanalı tek bir kanalda toplayan `.mix()` operatörünü kullanacağız.

A, B, C ve D adlı dört işlemimiz olsaydı ve her birinin basit bir `.out` kanalı olsaydı, sözdizimi şöyle görünürdü: `A.out.mix( B.out, C.out, D.out )`. Gördüğünüz gibi, onu birleştirmek istediğiniz kanalların ilkine (hangisi olduğu önemli değil) uyguluyorsunuz ve virgülle ayrılmış diğerlerini takip eden parantez içine ekliyorsunuz.

İş akışımızda, toplanması gereken şu çıktılar var:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Yani sözdizimi örneği şu hale geliyor:

```groovy title="MULTIQC çağrısında .mix() uygulaması"
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        )
```

Bu, örnek başına QC raporlarını toplayacak.
Ancak bunları tüm örnekler arasında toplamak istediğimiz için, tüm örneklere ait raporları tek bir `MULTIQC` çağrısında toplamak amacıyla `collect()` operatörünü eklememiz gerekiyor.
Ayrıca ona `report_id` parametresini de vermemiz gerekiyor.

Bu bize şunu verir:

```groovy title="Tamamlanmış MULTIQC çağrısı" linenums="33"
    // Kapsamlı QC raporu oluşturma
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

Tam iş akışı bloğu bağlamında, sonuç şöyle görünüyor:

```groovy title="rnaseq.nf" linenums="18"
workflow {
    // Bir CSV dosyasının içeriğinden girdi kanalı oluşturma
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }

    /// İlk kalite kontrolü
    FASTQC(read_ch)

    // Adaptör kırpma ve kırpma sonrası QC
    TRIM_GALORE(read_ch)

    // Bir referans genoma hizalama
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))

    // Kapsamlı QC raporu oluşturma
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
}
```

### 2.5. İş akışının çalıştığını test etmek için çalıştırma

```bash
nextflow run rnaseq.nf -resume
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

Bu sefer önbelleğe alınmış işlem çağrılarından sonra eklenen tek bir MULTIQC çağrısı görüyoruz:

Çıktıları `TRIM_GALORE` işleminde `publishDir` yönergesi tarafından belirtildiği gibi `results/trimming` altında bulabilirsiniz.

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

---

## 3. Paired-end RNAseq verilerinin işlenmesini etkinleştirme

Şu anda iş akışımız yalnızca single-end RNAseq verilerini işleyebiliyor.
Paired-end RNAseq verilerini görmek giderek daha yaygın hale geliyor, bu nedenle bunu da işleyebilmek istiyoruz.

İş akışını veri türünden tamamen bağımsız hale getirmek biraz daha gelişmiş Nextflow dil özelliklerini kullanmayı gerektirecektir, bu yüzden bunu burada yapmayacağız, ancak neyin uyarlanması gerektiğini göstermek için paired-end işleme versiyonu yapabiliriz.

### 3.1. `rnaseq_pe.nf` adlı bir iş akışı kopyası oluşturma

```bash
cp rnaseq.nf rnaseq_pe.nf
```

### 3.2. Varsayılan `input_csv`'yi paired-end verilere işaret edecek şekilde değiştirme

`data/` dizininde örnek ID'leri ve paired FASTQ dosya yollarını içeren ikinci bir CSV dosyası sağlıyoruz

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

`input_csv` varsayılanını `paired-end.csv` dosyasının yolu olacak şekilde değiştirelim.

```groovy title="rnaseq_pe.nf" linenums="15"
params {
    // Birincil girdi
    input_csv: Path = "data/paired-end.csv"

    // Referans genom arşivi
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Rapor ID'si
    report_id: String = "all_single-end"
}
```

### 3.3. Kanal fabrikasını güncelleme

`.map()` operatörüne artık her iki FASTQ dosya yolunu da almasını söylememiz gerekiyor.

Yani `row -> file(row.fastq_path)`, `row -> [file(row.fastq_1), file(row.fastq_2)]` olur

```groovy title="rnaseq_pe.nf" linenums="19"
    // Bir CSV dosyasının içeriğinden girdi kanalı oluşturma
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
```

### 3.4. FASTQC işleminin paired-end versiyonunu oluşturma

Her iki versiyonu da elimizde bulundurmak için modülün bir kopyasını oluşturalım.

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

Yeni `fastqc_pe.nf` modül dosyasını kod düzenleyicide açın ve aşağıdaki kod değişikliklerini yapın:

- `script` bloğunda (satır 17) `fastqc $reads`'i `fastqc ${reads}` olarak değiştirin, böylece `reads` girdisi artık tek bir yol yerine iki yoldan oluşan bir demet olduğu için açılacaktır.
- Çıktı dosyalarını ayrı ayrı işlemek zorunda kalmamak için `${reads.simpleName}`'i bir joker karakterle (`*`) değiştirin.

```groovy title="modules/fastqc_pe.nf" linenums="8"
    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads}
    """
```

Teknik olarak bu, `FASTQC` işlemini single-end veya paired-end RNAseq verilerinden birini işleyebilecek şekilde genelleştirir.

Son olarak, modülün paired-end versiyonunu kullanmak için modül içe aktarma ifadesini güncelleyin.

```groovy title="rnaseq_pe.nf" linenums="4"
include { FASTQC } from './modules/fastqc_pe.nf'
```

### 3.5. TRIM_GALORE işleminin paired-end versiyonunu oluşturma

Her iki versiyonu da elimizde bulundurmak için modülün bir kopyasını oluşturun.

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Yeni `trim_galore_pe.nf` modül dosyasını kod düzenleyicide açın ve aşağıdaki kod değişikliklerini yapın:

- Girdi bildirimini `path reads`'ten `tuple path(read1), path(read2)` olarak değiştirin
- `script` bloğundaki komutu güncelleyin, `$reads`'i `--paired ${read1} ${read2}` ile değiştirin
- Eklenen dosyaları ve farklı adlandırma kurallarını yansıtacak şekilde çıktı bildirimlerini güncelleyin, her şeyi listelemek zorunda kalmamak için joker karakterler kullanın.

```groovy title="modules/trim_galore_pe.nf" linenums="8"
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

Son olarak, modülün paired-end versiyonunu kullanmak için modül içe aktarma ifadesini güncelleyin.

```groovy title="rnaseq_pe.nf" linenums="5"
include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
```

### 3.6. MULTIQC işlemine yapılan çağrıyı TRIM_GALORE'dan iki rapor bekleyecek şekilde güncelleme

`TRIM_GALORE` işlemi artık ek bir çıktı kanalı üretiyor, bu yüzden bunu MultiQC'ye beslememiz gerekiyor.

`TRIM_GALORE.out.fastqc_reports,`'i `TRIM_GALORE.out.fastqc_reports_1,` artı `TRIM_GALORE.out.fastqc_reports_2,` ile değiştirin:

```groovy title="rnaseq_pe.nf" linenums="33"
    // Kapsamlı QC raporu oluşturma
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports_1,
        TRIM_GALORE.out.fastqc_reports_2,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

MultiQC üzerindeyken, `report_id` parametresinin varsayılanını da `"all_single-end"`'den `"all_paired-end"`'e güncelleyelim.

```groovy title="rnaseq_pe.nf" linenums="9"
params {
    // Birincil girdi
    input_csv: Path = "data/paired-end.csv"

    // Referans genom arşivi
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Rapor ID'si
    report_id: String = "all_paired-end"
}
```

### 3.7. HISAT2_ALIGN işleminin paired-end versiyonunu oluşturma

Her iki versiyonu da elimizde bulundurmak için modülün bir kopyasını oluşturun.

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Yeni `hisat2_align_pe.nf` modül dosyasını kod düzenleyicide açın ve aşağıdaki kod değişikliklerini yapın:

- Girdi bildirimini `path reads`'ten `tuple path(read1), path(read2)` olarak değiştirin
- `script` bloğundaki komutu güncelleyin, `-U $reads`'i `-1 ${read1} -2 ${read2}` ile değiştirin
- `${reads.simpleName}`'in tüm örneklerini `script` bloğundaki komutta ve çıktı bildirimlerinde `${read1.simpleName}` ile değiştirin.

```groovy title="modules/hisat2_align_pe.nf" linenums="8"
    input:
    tuple path(read1), path(read2)
    path index_zip

    output:
    path "${read1.simpleName}.bam", emit: bam
    path "${read1.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
        --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
        samtools view -bS -o ${read1.simpleName}.bam
    """
```

Son olarak, modülün paired-end versiyonunu kullanmak için modül içe aktarma ifadesini güncelleyin.

```groovy title="rnaseq_pe.nf" linenums="5"
include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
```

### 3.8. İş akışının çalıştığını test etmek için çalıştırma

Bu önbelleğe almayacağı ve işlenecek verilerin öncekinden iki kat daha fazla olduğu için `-resume` kullanmıyoruz, ancak yine de bir dakikadan kısa sürede tamamlanması gerekiyor.

```bash
nextflow run rnaseq_pe.nf
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

İşte bu kadar! Şimdi iş akışımızın biri single-end okuma verileri için, diğeri paired-end veriler için olmak üzere biraz farklılaşan iki versiyonuna sahibiz.
Bir sonraki mantıklı adım, iş akışının her iki veri türünü de anında kabul etmesini sağlamak olacaktır; bu, bu kursun kapsamı dışındadır, ancak bunu bir takip kursunda ele alabiliriz.

---

### Çıkarımlar

Tek örnekli bir iş akışını birden fazla örneğin işlenmesini paralelleştirmek, kapsamlı bir QC raporu oluşturmak ve gerekirse iş akışını paired-end okuma verilerini kullanacak şekilde uyarlamak için nasıl uyarlayacağınızı biliyorsunuz.

### Sırada ne var?

Tebrikler, Nextflow For RNAseq mini-kursunu tamamladınız! Başarınızı kutlayın ve hak ettiğiniz bir mola verin!

Ardından, bu eğitim kursu hakkındaki deneyiminizle ilgili çok kısa bir anketi tamamlamanızı rica ediyoruz, sonra sizi daha fazla eğitim kaynağı ve yararlı bağlantılara sahip bir sayfaya götüreceğiz.
