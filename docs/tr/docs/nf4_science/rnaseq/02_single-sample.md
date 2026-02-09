# Bölüm 2: Tek örnekli uygulama

Bu kursun bu bölümünde, Bölüm 1'de çalıştırdığımız tüm komutları otomatikleştirmek için mümkün olan en basit iş akışını yazacağız ve her seferinde yalnızca bir örneği işlemeyi hedefleyeceğiz.

Bunu üç aşamada yapacağız:

1. İlk QC adımını çalıştıran tek aşamalı bir iş akışı yazma
2. Adaptör kırpma ve kırpma sonrası QC ekleme
3. Referans genoma hizalama ekleme

!!! warning "Ön Koşul"

    Bu derse başlamadan önce kursun Bölüm 1'ini tamamlamış olmanız gerekir.
    Özellikle, 2.1-3 bölümlerinde çalışmak, bu dersteki hizalama adımı için gerekli olan genom dizin dosyasını (`data/genome_index.tar.gz`) oluşturur.

---

## 1. İlk QC'yi çalıştıran tek aşamalı bir iş akışı yazma

Tek uçlu RNAseq okumalarını içeren bir FASTQ dosyası üzerinde FastQC aracını çalıştıran basit bir iş akışı yazarak başlayalım.

Size iş akışının ana bölümlerini özetleyen bir iş akışı dosyası olan `rnaseq.nf` sağlıyoruz.

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Module INCLUDE statements

/*
 * Pipeline parametreleri
 */

// Primary input

workflow {

    // Create input channel

    // Call processes

}
```

Bu iş akışı kodunun doğru olduğunu ancak işlevsel olmadığını unutmayın; amacı sadece gerçek iş akışını yazmak için kullanacağınız bir iskelet görevi görmektir.

### 1.1. Modülleri depolamak için bir dizin oluşturma

Her süreç için bağımsız modüller oluşturacağız, böylece bunları yönetmek ve yeniden kullanmak daha kolay olacak, bu yüzden bunları depolamak için bir dizin oluşturalım.

```bash
mkdir modules
```

### 1.2. QC metrikleri toplama süreci için bir modül oluşturma

`FASTQC` sürecini barındırmak için `modules/fastqc.nf` adında bir modül dosyası oluşturalım:

```bash
touch modules/fastqc.nf
```

Dosyayı kod düzenleyicide açın ve aşağıdaki kodu içine kopyalayın:

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

Bu eğitim serisinin Bölüm 1 ve Bölüm 2'sinde öğrendiklerinizden tüm parçaları tanıyor olmalısınız; dikkate değer tek değişiklik, bu sefer `publishDir` yönergesi için `mode: symlink` kullanıyor olmamız ve `publishDir`'i tanımlamak için bir parametre kullanıyor olmamızdır.

!!! note

    Burada kullandığımız veri dosyaları çok küçük olsa da, genomik alanında çok büyük olabilirler. Öğretim ortamında gösterim amacıyla, gereksiz dosya kopyalarından kaçınmak için 'symlink' yayınlama modunu kullanıyoruz. Bunu nihai iş akışlarınızda yapmamalısınız, çünkü `work` dizininizi temizlediğinizde sonuçları kaybedersiniz.

### 1.3. Modülü iş akışı dosyasına aktarma

`rnaseq.nf` dosyasına `include { FASTQC } from './modules/fastqc.nf'` ifadesini ekleyin:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
```

### 1.4. Bir girdi bildirimi ekleme

Varsayılan değere sahip bir girdi parametresi bildirin:

```groovy title="rnaseq.nf" linenums="10"
params {
    // Primary input
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
}
```

### 1.5. İş akışı bloğunda bir girdi kanalı oluşturma

Girdi kanalını oluşturmak için temel bir `.fromPath()` kanal fabrikası kullanın:

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Call processes

}
```

### 1.6. Girdi kanalı üzerinde `FASTQC` sürecini çağırma

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

}
```

### 1.7. Çalıştığını test etmek için iş akışını çalıştırma

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

Bölüm 1'de çalıştıysanız ve konteyneri zaten çektiyseniz bu çok hızlı çalışmalıdır.
Eğer atladıysanız, Nextflow konteyneri sizin için çekecektir; bunun gerçekleşmesi için bir şey yapmanıza gerek yok, ancak bir dakikaya kadar beklemeniz gerekebilir.

Çıktıları `FASTQC` sürecinde `publishDir` yönergesi tarafından belirtildiği gibi `results/fastqc` altında bulabilirsiniz.

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

## 2. Adaptör kırpma ve kırpma sonrası kalite kontrolü ekleme

Kırpmanın kendisi için Cutadapt'ı ve kırpma sonrası kalite kontrolü için FastQC'yi bir araya getiren Trim_Galore sarmalayıcısını kullanacağız.

### 2.1. Kırpma ve QC süreci için bir modül oluşturma

`TRIM_GALORE` sürecini barındırmak için `modules/trim_galore.nf` adında bir modül dosyası oluşturalım:

```bash
touch modules/trim_galore.nf
```

Dosyayı kod düzenleyicide açın ve aşağıdaki kodu içine kopyalayın:

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

### 2.2. Modülü iş akışı dosyasına aktarma

`rnaseq.nf` dosyasına `include { TRIM_GALORE } from './modules/trim_galore.nf'` ifadesini ekleyin:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. Girdi kanalı üzerinde süreci çağırma

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)
}
```

### 2.4. Çalıştığını test etmek için iş akışını çalıştırma

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

Çıktıları `TRIM_GALORE` sürecinde `publishDir` yönergesi tarafından belirtildiği gibi `results/trimming` altında bulabilirsiniz.

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. Okumaları referans genoma hizalama

Son olarak, FastQC tarzı kalite kontrol metriklerini de yayınlayacak olan Hisat2 kullanarak genom hizalama adımını çalıştırabiliriz.

### 3.1. HiSat2 süreci için bir modül oluşturma

`HISAT2_ALIGN` sürecini barındırmak için `modules/hisat2_align.nf` adında bir modül dosyası oluşturalım:

```bash
touch modules/hisat2_align.nf
```

Dosyayı kod düzenleyicide açın ve aşağıdaki kodu içine kopyalayın:

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

### 3.2. Modülü iş akışı dosyasına aktarma

`rnaseq.nf` dosyasına `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` ifadesini ekleyin:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. Genom dizinini sağlamak için bir parametre bildirimi ekleme

Varsayılan değere sahip bir girdi parametresi bildirin:

```groovy title="rnaseq.nf" linenums="8"
params {
    // Primary input
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 3.4. `TRIM_GALORE` tarafından çıktılanan kırpılmış okumalar üzerinde `HISAT2_ALIGN` sürecini çağırma

Kırpılmış okumalar, önceki adım tarafından çıktılanan `TRIM_GALORE.out.trimmed_reads` kanalındadır.

Ek olarak, Hisat2 aracına sıkıştırılmış genom dizin tarball'ını sağlamak için `file (params.hisat2_index_zip)` kullanıyoruz.

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.5. Çalıştığını test etmek için iş akışını çalıştırma

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

Çıktıları `HISAT2_ALIGN` sürecinde `publishDir` yönergesi tarafından belirtildiği gibi `results/align` altında bulabilirsiniz.

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Bu, her örneğe uygulamamız gereken temel işlemeyi tamamlar.

_MultiQC rapor toplama işlemini Bölüm 2'de, iş akışını aynı anda birden fazla örneği kabul edecek şekilde değiştirdikten sonra ekleyeceğiz._

---

### Özet

Tek uçlu RNAseq örneklerini ayrı ayrı işlemek için tüm temel adımları nasıl saracağınızı biliyorsunuz.

### Sırada ne var?

İş akışını birden fazla örneği paralel olarak işleyecek şekilde nasıl değiştireceğinizi, tüm örnekler için tüm adımlarda QC raporlarını nasıl toplayacağınızı ve iş akışının çift uçlu RNAseq verileri üzerinde çalışmasını nasıl etkinleştireceğinizi öğrenin.
