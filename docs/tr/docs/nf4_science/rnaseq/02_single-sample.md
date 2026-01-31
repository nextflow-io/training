# Bölüm 2: Tek örnek implementasyonu

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu kurumun bu bölümünde, Bölüm 1'de çalıştırdığımız tüm komutları otomatikleştirmek için mümkün olan en basit iş akışını yazacağız ve her seferinde yalnızca bir örnek işlemeyi hedefleyeceğiz.

Bunu üç aşamada yapacağız:

1. İlk kalite kontrol adımını çalıştıran tek aşamalı bir iş akışı yazın
2. Adaptör kırpma ve kırpma sonrası kalite kontrolünü ekleyin
3. Referans genoma hizalamayı ekleyin

!!! warning "Ön Koşul"

    Bu derse başlamadan önce kursun Bölüm 1'ini tamamlamanız gerekmektedir.
    Özellikle, 2.1-3 bölümlerini tamamlamak bu dersteki hizalama adımı için gereken genom dizin dosyasını (`data/genome_index.tar.gz`) oluşturur.

---

## 1. İlk kalite kontrolünü çalıştıran tek aşamalı bir iş akışı yazın

Tek uçlu RNAseq okumalarını içeren bir FASTQ dosyası üzerinde FastQC aracını çalıştıran basit bir iş akışı yazarak başlayalım.

Size iş akışının ana bölümlerini özetleyen `rnaseq.nf` adında bir iş akışı dosyası sağlıyoruz.

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Modül INCLUDE ifadeleri

/*
 * Pipeline parametreleri
 */

// Birincil girdi

workflow {

    // Girdi kanalı oluştur

    // Process'leri çağır

}
```

Bu iş akışı kodunun doğru ancak işlevsel olmadığını unutmayın; amacı yalnızca gerçek iş akışını yazmak için kullanacağınız bir iskelet görevi görmektir.

### 1.1. Modülleri depolamak için bir dizin oluşturun

Her bir process için bağımsız modüller oluşturacağız, böylece bunları yönetmeyi ve yeniden kullanmayı kolaylaştıracağız. Hadi bunları depolamak için bir dizin oluşturalım.

```bash
mkdir modules
```

### 1.2. Kalite kontrol metrikleri toplama process'i için bir modül oluşturun

`FASTQC` process'ini barındırmak için `modules/fastqc.nf` adında bir modül dosyası oluşturalım:

```bash
touch modules/fastqc.nf
```

Dosyayı kod düzenleyicide açın ve içine aşağıdaki kodu kopyalayın:

```groovy title="modules/fastqc.nf" linenums="1"
#!/usr/bin/env nextflow

process FASTQC {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/fastqc", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_fastqc.zip", emit: zip
    path "${reads.simpleName}_fastqc.html", emit: html

    script:
    """
    fastqc $reads
    """
}
```

Bu eğitim serisinin Bölüm 1 ve Bölüm 2'de öğrendiklerinizden tüm parçaları tanıyor olmalısınız; dikkate değer tek değişiklik bu sefer `publishDir` yönergesi için `mode: symlink` kullanıyor olmamız ve `publishDir` tanımlamak için bir parametre kullanıyor olmamızdır.

!!! note

    Burada kullandığımız veri dosyaları çok küçük olsa da, genomikte çok büyük olabilirler. Öğretim ortamında gösterim amacıyla, gereksiz dosya kopyalarından kaçınmak için 'symlink' yayınlama modunu kullanıyoruz. Bunu son iş akışlarınızda yapmamalısınız, çünkü `work` dizininizi temizlediğinizde sonuçları kaybedersiniz.

### 1.3. Modülü iş akışı dosyasına import edin

`rnaseq.nf` dosyasına `include { FASTQC } from './modules/fastqc.nf'` ifadesini ekleyin:

```groovy title="rnaseq.nf" linenums="3"
// Modül INCLUDE ifadeleri
include { FASTQC } from './modules/fastqc.nf'
```

### 1.4. Bir girdi bildirimi ekleyin

Varsayılan değere sahip bir girdi parametresi bildirin:

```groovy title="rnaseq.nf" linenums="10"
params {
    // Birincil girdi
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
}
```

### 1.5. İş akışı bloğunda bir girdi kanalı oluşturun

Girdi kanalını oluşturmak için temel bir `.fromPath()` kanal fabrikası kullanın:

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Girdi kanalı oluştur from a file path
    read_ch = channel.fromPath(params.reads)

    // Process'leri çağır

}
```

### 1.6. Girdi kanalı üzerinde `FASTQC` process'ini çağırın

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Girdi kanalı oluştur from a file path
    read_ch = channel.fromPath(params.reads)

    // İlk kalite kontrolü
    FASTQC(read_ch)

}
```

### 1.7. Çalıştığını test etmek için iş akışını çalıştırın

Komut satırından bir girdi belirtmek için `--reads` parametresini kullanabiliriz, ancak geliştirme sırasında tembel olabilir ve kurduğumuz test varsayılanını kullanabiliriz.

```bash
nextflow run rnaseq.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    ```

Bölüm 1'i tamamladıysanız ve konteyneri zaten çektiyseniz bu çok hızlı çalışmalıdır.
Eğer atladıysanız, Nextflow konteyneri sizin için çekecektir; bunun gerçekleşmesi için bir şey yapmanıza gerek yok, ancak bir dakika kadar beklemeniz gerekebilir.

Çıktıları `FASTQC` process'inde `publishDir` yönergesi tarafından belirtildiği şekilde `results/fastqc` altında bulabilirsiniz.

```bash
ls results/fastqc
```

```console title="Çıktı"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

## 2. Adaptör kırpma ve kırpma sonrası kalite kontrolünü ekleyin

Kırpma işlemi için Cutadapt'i ve kırpma sonrası kalite kontrolü için FastQC'yi bir araya getiren Trim_Galore sarmalayıcısını kullanacağız.

### 2.1. Kırpma ve kalite kontrol process'i için bir modül oluşturun

`TRIM_GALORE` process'ini barındırmak için `modules/trim_galore.nf` adında bir modül dosyası oluşturalım:

```bash
touch modules/trim_galore.nf
```

Dosyayı kod düzenleyicide açın ve içine aşağıdaki kodu kopyalayın:

```groovy title="modules/trim_galore.nf" linenums="1"
#!/usr/bin/env nextflow

process TRIM_GALORE {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/trimming", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
    path "${reads}_trimming_report.txt", emit: trimming_reports
    path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

    script:
    """
    trim_galore --fastqc $reads
    """
}
```

### 2.2. Modülü iş akışı dosyasına import edin

`rnaseq.nf` dosyasına `include { TRIM_GALORE } from './modules/trim_galore.nf'` ifadesini ekleyin:

```groovy title="rnaseq.nf" linenums="3"
// Modül INCLUDE ifadeleri
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. Girdi kanalı üzerinde process'i çağırın

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // Girdi kanalı oluştur from a file path
    read_ch = channel.fromPath(params.reads)

    // İlk kalite kontrolü
    FASTQC(read_ch)

    // Adaptör kırpma ve kırpma sonrası QC
    TRIM_GALORE(read_ch)
}
```

### 2.4. Çalıştığını test etmek için iş akışını çalıştırın

```bash
nextflow run rnaseq.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    [c2/e4a9bb] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    ```

Bu da çok hızlı çalışmalıdır, çünkü çok küçük bir girdi dosyası üzerinde çalışıyoruz.

Çıktıları `TRIM_GALORE` process'inde `publishDir` yönergesi tarafından belirtildiği şekilde `results/trimming` altında bulabilirsiniz.

```bash
ls results/trimming
```

```console title="Çıktı"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. Okumaları referans genoma hizalayın

Son olarak, FastQC tarzı kalite kontrol metriklerini de yayınlayacak olan Hisat2 kullanarak genom hizalama adımını çalıştırabiliriz.

### 3.1. HiSat2 process'i için bir modül oluşturun

`HISAT2_ALIGN` process'ini barındırmak için `modules/hisat2_align.nf` adında bir modül dosyası oluşturalım:

```bash
touch modules/hisat2_align.nf
```

Dosyayı kod düzenleyicide açın ve içine aşağıdaki kodu kopyalayın:

```groovy title="modules/hisat2_align.nf" linenums="1"
#!/usr/bin/env nextflow

process HISAT2_ALIGN {

    container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"
    publishDir "results/align", mode: 'symlink'

    input:
    path reads
    path index_zip

    output:
    path "${reads.simpleName}.bam", emit: bam
    path "${reads.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -U $reads \
        --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
        samtools view -bS -o ${reads.simpleName}.bam
    """
}
```

### 3.2. Modülü iş akışı dosyasına import edin

`rnaseq.nf` dosyasına `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` ifadesini ekleyin:

```groovy title="rnaseq.nf" linenums="3"
// Modül INCLUDE ifadeleri
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. Genom dizinini sağlamak için bir parametre bildirimi ekleyin

Varsayılan değere sahip bir girdi parametresi bildirin:

```groovy title="rnaseq.nf" linenums="8"
params {
    // Birincil girdi
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

    // Referans genom arşivi
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 3.4. `TRIM_GALORE` tarafından çıktılanan kırpılmış okumalar üzerinde `HISAT2_ALIGN` process'ini çağırın

Kırpılmış okumalar, önceki adım tarafından çıktılanan `TRIM_GALORE.out.trimmed_reads` kanalındadır.

Ek olarak, Hisat2 aracına sıkıştırılmış genom dizini tarball'ını sağlamak için `file (params.hisat2_index_zip)` kullanıyoruz.

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // Girdi kanalı oluştur from a file path
    read_ch = channel.fromPath(params.reads)

    // İlk kalite kontrolü
    FASTQC(read_ch)

    // Adaptör kırpma ve kırpma sonrası QC
    TRIM_GALORE(read_ch)

    // Referans genoma hizalama
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.5. Çalıştığını test etmek için iş akışını çalıştırın

```bash
nextflow run rnaseq.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [extravagant_khorana] DSL2 - revision: 701b41bd16

    executor >  local (3)
    [e4/d15ad4] FASTQC (1)       [100%] 1 of 1 ✔
    [c6/12b2be] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    [c6/7a9f13] HISAT2_ALIGN (1) [100%] 1 of 1 ✔
    ```

Çıktıları `HISAT2_ALIGN` process'inde `publishDir` yönergesi tarafından belirtildiği şekilde `results/align` altında bulabilirsiniz.

```bash
ls results/align
```

```console title="Çıktı"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Bu, her bir örneğe uygulamamız gereken temel işlemi tamamlar.

_MultiQC rapor toplama işlemini Bölüm 2'de, iş akışını aynı anda birden fazla örneği kabul edecek şekilde değiştirdikten sonra ekleyeceğiz._

---

### Çıkarımlar

Tek uçlu RNAseq örneklerini tek tek işlemek için tüm temel adımları nasıl sarmalayacağınızı biliyorsunuz.

### Sırada ne var?

İş akışını birden fazla örneği paralel olarak işleyecek şekilde nasıl değiştireceğinizi, tüm örnekler için tüm adımlarda kalite kontrol raporlarını nasıl toplayacağınızı ve iş akışını çift uçlu RNAseq verileri üzerinde çalışacak şekilde nasıl etkinleştireceğinizi öğrenin.
