# Bölüm 1: Yönteme genel bakış ve manuel test

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Toplu RNAseq verilerini işlemek ve analiz etmek için birden fazla geçerli yöntem bulunmaktadır.
Bu kurs için, [Babraham Enstitüsü](https://www.babraham.ac.uk/)'nden Dr. Simon Andrews ve Laura Biggins tarafından [burada](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) açıklanan yöntemi takip ediyoruz.

Amacımız, aşağıdaki işleme adımlarını uygulayan bir iş akışı geliştirmektir: toplu RNAseq örneğindeki okumalar üzerinde ilk kalite kontrolünü çalıştırmak, okumalardan adaptör dizilerini kırpmak, okumaları referans genoma hizalamak ve kapsamlı bir kalite kontrol (QC) raporu oluşturmak.

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/preprocess.svg"
</figure>

- **FASTQC:** Kırpma öncesi okuma verileri üzerinde FastQC kullanarak QC gerçekleştirin
- **TRIM_GALORE:** Trim Galore kullanarak adaptör dizilerini kırpın ve kırpma sonrası QC gerçekleştirin (Cutadapt ve FastQC'yi bir araya getirir)
- **HISAT2_ALIGN:** Hisat2 kullanarak okumaları referans genoma hizalayın
- **MULTIQC:** MultiQC kullanarak kapsamlı bir QC raporu oluşturun

Ancak, herhangi bir iş akışı kodu yazmaya başlamadan önce, komutları bazı test verileri üzerinde manuel olarak deneyeceğiz.
İhtiyacımız olan araçlar GitHub Codespaces ortamında yüklü olmadığından, bunları konteynerler aracılığıyla kullanacağız (bkz. [Merhaba Konteynerler](../../hello_nextflow/05_hello_containers.md)).

!!! note "Not"

     `nf4-science/rnaseq` dizininde olduğunuzdan emin olun. `pwd` yazdığınızda gösterilen yolun son kısmı `rnaseq` olmalıdır.

---

## 1. İlk QC ve adaptör kırpma

Hem `fastqc` hem de `trim_galore` yüklü bir konteyner imajını çekeceğiz, etkileşimli olarak başlatacağız ve örnek veri dosyalarından biri üzerinde kırpma ve QC komutlarını çalıştıracağız.

### 1.1. Konteyneri çekin

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Bu, sistem imajı indirirken aşağıdaki konsol çıktısını verir:

??? success "Komut çıktısı"

    ```console
    0.6.10--1bf8ca4e1967cd18: Pulling from library/trim-galore
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    32ec762be2d0: Pull complete
    d2cb90387285: Pull complete
    Digest: sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907
    Status: Downloaded newer image for community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    ```

### 1.2. Konteyneri etkileşimli olarak başlatın

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

İsteminiz `(base) root@b645838b3314:/tmp#` gibi bir şeye dönüşecek, bu da artık konteynerin içinde olduğunuzu gösterir.

Komutun `-v ./data:/data` kısmı, `data/` dizininin içeriğine konteyner içinden erişmemizi sağlayacaktır.

```bash
ls /data/reads
```

??? success "Komut çıktısı"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gzO
    ```

### 1.3. İlk `fastqc` komutunu çalıştırın

Okuma verileri üzerinde kalite kontrol metriklerini toplamak için `fastqc`'yi çalıştıralım.

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Komut çıktısı"

    ```console
    application/gzip
    Started analysis of ENCSR000COQ1_1.fastq.gz
    Approx 5% complete for ENCSR000COQ1_1.fastq.gz
    Approx 10% complete for ENCSR000COQ1_1.fastq.gz
    Approx 15% complete for ENCSR000COQ1_1.fastq.gz
    Approx 20% complete for ENCSR000COQ1_1.fastq.gz
    Approx 25% complete for ENCSR000COQ1_1.fastq.gz
    Approx 30% complete for ENCSR000COQ1_1.fastq.gz
    Approx 35% complete for ENCSR000COQ1_1.fastq.gz
    Approx 40% complete for ENCSR000COQ1_1.fastq.gz
    Approx 45% complete for ENCSR000COQ1_1.fastq.gz
    Approx 50% complete for ENCSR000COQ1_1.fastq.gz
    Approx 55% complete for ENCSR000COQ1_1.fastq.gz
    Approx 60% complete for ENCSR000COQ1_1.fastq.gz
    Approx 65% complete for ENCSR000COQ1_1.fastq.gz
    Approx 70% complete for ENCSR000COQ1_1.fastq.gz
    Approx 75% complete for ENCSR000COQ1_1.fastq.gz
    Approx 80% complete for ENCSR000COQ1_1.fastq.gz
    Approx 85% complete for ENCSR000COQ1_1.fastq.gz
    Approx 90% complete for ENCSR000COQ1_1.fastq.gz
    Approx 95% complete for ENCSR000COQ1_1.fastq.gz
    Analysis complete for ENCSR000COQ1_1.fastq.gz
    ```

Bu çok hızlı çalışmalıdır.
Çıktı dosyalarını orijinal verilerle aynı dizinde bulabilirsiniz:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

<!-- switch to tree -->

```console title="Çıktı"
/data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
```

### 1.4. `trim_galore` ile adaptör dizilerini kırpın

Şimdi adaptör dizilerini kırpmak ve kırpma sonrası QC metriklerini toplamak için Cutadapt ve FastQC'yi bir araya getiren `trim_galore`'ı çalıştıralım.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

`--fastqc` bayrağı, kırpma tamamlandıktan sonra komutun otomatik olarak bir QC toplama adımı çalıştırmasını sağlar.

_Çıktı çok ayrıntılı olduğu için aşağıdaki kısaltılmıştır._

??? success "Komut çıktısı"

    ```console
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    <...>

    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

Çıktı dosyalarını çalışma dizininde bulabilirsiniz:

```bash
ls ENCSR000COQ1_1*
```

```console title="Çıktı"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
```

### 1.5. Çıktı dosyalarını konteyner dışındaki dosya sistemine taşıyın

Konteyner içinde kalan her şey gelecekteki çalışmalara erişilemez olacağından bunları yeni bir dizine taşıyalım.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

### 1.6. Konteynerden çıkın

```bash
exit
```

---

## 2. Okumaları referans genoma hizalayın

`hisat2`'nin yüklü olduğu bir konteyner imajını çekeceğiz, etkileşimli olarak başlatacağız ve RNAseq verilerini referans genoma hizalamak için hizalama komutunu çalıştıracağız.

### 2.1. `hisat2` konteynerini çekin

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "Komut çıktısı"

    ```console
    Unable to find image 'community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e' locally
    5e49f68a37dc010e: Pulling from library/hisat2_samtools
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    e74ed5dd390b: Pull complete
    abfcf0185e51: Pull complete
    Digest: sha256:29d8e1a3172a2bdde7be813f7ebec22d331388194a7c0de872b4ccca4bed8f45
    Status: Downloaded newer image for community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
    ```

### 2.2. `hisat2` konteynerini etkileşimli olarak başlatın

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Komut daha öncekiyle aynıdır, ilgili konteyner URI'si değiştirilmiştir.

### 2.3. Hisat2 genom dizin dosyalarını oluşturun

Hisat2, genom referansının çok özel bir biçimde sağlanmasını gerektirir ve sağladığımız `genome.fa` FASTA dosyasını doğrudan kullanamaz, bu nedenle bu fırsatı ilgili kaynakları oluşturmak için kullanacağız.

```bash
hisat2-build /data/genome.fa genome_index
```

Çıktı çok ayrıntılı olduğundan aşağıdaki kısaltılmıştır:

<!-- TODO: switch to full output -->

??? success "Komut çıktısı"

    ```console
    Settings:
      Output files: "genome_index.*.ht2"
    <...>
    Total time for call to driver() for forward index: 00:00:16
    ```

Bu, çalışma dizininde bulabileceğiniz birden fazla genom dizin dosyası oluşturur.

```bash
ls genome_index.*
```

```console title="Çıktı"
genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
```

Bunları birazdan kullanacağız, ancak önce bu genom dizin dosyalarıyla sıkıştırılmış bir tarball oluşturalım; bunlara daha sonra ihtiyacımız olacak ve bunları oluşturmak tipik olarak bir iş akışının parçası olarak yapmak istediğimiz bir şey değildir.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

Bu, genom dizin dosyalarını içeren bir `genome_index.tar.gz` tarball'ını dosya sistemimizde `data/` dizininde saklar, bu da bu kursun 2. Bölümünde işimize yarayacaktır.

### 2.4. `hisat2` komutunu çalıştırın

Şimdi hizalama komutunu çalıştırabiliriz, bu komut `hisat2` ile hizalama adımını gerçekleştirir ve ardından çıktıyı bir BAM dosyası olarak yazmak için `samtools`'a aktarır.

Okuma verisi girdisi, önceki adımda `trim_galore` ile oluşturduğumuz `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` dosyasıdır.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Komut çıktısı"

    ```console
    HISAT2 summary stats:
            Total reads: 27816
                    Aligned 0 time: 1550 (5.57%)
                    Aligned 1 time: 25410 (91.35%)
                    Aligned >1 times: 856 (3.08%)
            Overall alignment rate: 94.43%
    ```

Bu, çok küçük bir test dosyası olduğu için neredeyse anında çalışır.
Gerçek ölçekte bu çok daha uzun sürebilir.

Bir kez daha çıktı dosyalarını çalışma dizininde bulabilirsiniz:

```bash
ls ENCSR000COQ1_1*
```

```console title="Çıktı"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

### 2.5. Çıktı dosyalarını konteyner dışındaki dosya sistemine taşıyın

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

### 2.6. Konteynerden çıkın

```bash
exit
```

---

## 3. Kapsamlı bir QC raporu oluşturun

`multiqc`'nin yüklü olduğu bir konteyner imajını çekeceğiz, etkileşimli olarak başlatacağız ve öncesi/sonrası FastQC rapor dosyaları üzerinde bir rapor oluşturma komutu çalıştıracağız.

### 3.1. `multiqc` konteynerini çekin

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Komut çıktısı"

    ```console
    ad8f247edb55897c: Pulling from library/pip_multiqc
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    3f229294c69a: Pull complete
    5a5ad47fd84c: Pull complete
    Digest: sha256:0ebb1d9605395a7df49ad0eb366b21f46afd96a5090376b0d8941cf5294a895a
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

### 3.2. `multiqc` konteynerini etkileşimli olarak başlatın

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

### 3.3. `multiqc` komutunu çalıştırın

```bash
multiqc /data/reads /data/trimmed /data/aligned -n ENCSR000COQ1_1_QC
```

??? success "Komut çıktısı"

    ```console

    /// MultiQC 🔍 v1.27.1

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 20/20
                hisat2 | Found 1 reports
              cutadapt | Found 1 reports
                fastqc | Found 1 reports
        write_results | Data        : ENCSR000COQ1_1_QC_data
        write_results | Report      : ENCSR000COQ1_1_QC.html
              multiqc | MultiQC complete
    ```

MultiQC, uyumlu QC raporları için dizinlerde arama yapabilir ve bulduğu her şeyi birleştirir.

Burada aracın oluşturduğumuz üç QC raporunu da bulduğunu görüyoruz: `fastqc` ile yaptığımız ilk QC, `cutadapt`'ten (via `trim_galore`) kırpma sonrası rapor ve `hisat2` tarafından üretilen hizalama sonrası QC.

Çıktı dosyaları bir kez daha çalışma dizinindedir:

```bash
ls ENCSR000COQ1_1_QC*
```

```console title="Çıktı"
ENCSR000COQ1_1_QC.html

ENCSR000COQ1_1_QC_data:
cutadapt_filtered_reads_plot.txt                     fastqc_top_overrepresented_sequences_table.txt
cutadapt_trimmed_sequences_plot_3_Counts.txt         hisat2_se_plot.txt
cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc.log
fastqc-status-check-heatmap.txt                      multiqc_citations.txt
fastqc_adapter_content_plot.txt                      multiqc_cutadapt.txt
fastqc_per_base_n_content_plot.txt                   multiqc_data.json
fastqc_per_base_sequence_quality_plot.txt            multiqc_fastqc.txt
fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_general_stats.txt
fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_hisat2.txt
fastqc_per_sequence_quality_scores_plot.txt          multiqc_software_versions.txt
fastqc_sequence_counts_plot.txt                      multiqc_sources.txt
fastqc_sequence_duplication_levels_plot.txt
```

### 3.4. Çıktı dosyalarını konteyner dışındaki dosya sistemine taşıyın

```bash
mkdir /data/final_qc
mv ENCSR000COQ1_1_QC** /data/final_qc
```

### 3.5. Konteynerden çıkın

```bash
exit
```

---

### Çıkarımlar

Tüm bireysel komutları ilgili konteynerlerde etkileşimli olarak test ettiniz.

### Sırada ne var?

Aynı komutları, çalışmayı yürütmek için konteynerler kullanan çok adımlı bir iş akışına nasıl sarmalayacağınızı öğrenin.
