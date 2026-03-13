# Bölüm 1: Yönteme genel bakış

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Toplu RNAseq verilerini işlemek ve analiz etmek için birden fazla geçerli yöntem bulunmaktadır.
Bu kurs için, [Babraham Enstitüsü](https://www.babraham.ac.uk/)'nden Dr. Simon Andrews ve Laura Biggins tarafından [burada](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) açıklanan yöntemi takip ediyoruz.

Amacımız, aşağıdaki işleme adımlarını uygulayan bir iş akışı geliştirmektir: toplu RNAseq örneğindeki okumalar üzerinde ilk kalite kontrolünü çalıştırmak, okumalardan adaptör dizilerini kırpmak, okumaları referans genoma hizalamak ve kapsamlı bir kalite kontrol (QC) raporu oluşturmak.

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

- **FASTQC:** Kırpma öncesi okuma verileri üzerinde FastQC kullanarak QC gerçekleştirin
- **TRIM_GALORE:** Trim Galore kullanarak adaptör dizilerini kırpın ve kırpma sonrası QC gerçekleştirin (Cutadapt ve FastQC'yi bir araya getirir)
- **HISAT2_ALIGN:** Hisat2 kullanarak okumaları referans genoma hizalayın
- **MULTIQC:** MultiQC kullanarak kapsamlı bir QC raporu oluşturun

### Yöntemler

Bu işleme adımlarını iki aşamada nasıl uygulayacağınızı göstereceğiz.
İlk olarak, QC, kırpma ve hizalama araçlarını bir örnek üzerinde çalıştıran **tek örnekli işleme** ile başlayacağız.
Ardından, aynı araçları birden fazla örnek üzerinde çalıştıran ve toplu bir kalite kontrol raporu oluşturan **çok örnekli işleme**'ye geçeceğiz.

Her iki yaklaşım için de herhangi bir iş akışı kodu yazmaya başlamadan önce, komutları bazı test verileri üzerinde manuel olarak deneyeceğiz.

### Veri Seti

Aşağıdaki verileri ve ilgili kaynakları sağlıyoruz:

- **RNAseq verileri** (`reads/`): Altı örnekten FASTQ dosyaları, dosya boyutlarını küçük tutmak için küçük bir bölgeye indirilmiştir. Her örneğin eşleştirilmiş uç okumaları vardır (örnek başına iki dosya), ancak önce yalnızca tek uçlu okumalarla çalışmaya başlıyoruz.
- **Bir referans genom** (`genome.fa`): insan kromozom 20'sinin küçük bir bölgesi (hg19/b37'den).
- **CSV örnek listeleri** (`single-end.csv` ve `paired-end.csv`): örnek veri dosyalarının kimliklerini ve yollarını listeleyen dosyalar.

### Yazılım

İlgili dört ana araç, kalite kontrol metriklerini toplamak için [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), adaptör kırpma için [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) (Cutadapt ve FastQC'yi kırpma sonrası QC için bir araya getirir), referans genoma eklenmiş hizalama için [HISAT2](http://daehwankimlab.github.io/hisat2/) ve toplu QC raporu oluşturma için [MultiQC](https://multiqc.info/)'dir.

Bu araçlar GitHub Codespaces ortamında yüklü olmadığından, bunları Seqera Containers hizmeti aracılığıyla alınan konteynerler üzerinden kullanacağız (bkz. [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! tip "İpucu"

     `nf4-science/rnaseq` dizininde olduğunuzdan emin olun. `pwd` yazdığınızda gösterilen yolun son kısmı `rnaseq` olmalıdır.

---

## 1. Tek örnekli işleme

Bu bölümde, tek bir RNAseq örneğini işleyen komutları test ediyoruz: kalite kontrol, adaptör kırpma ve referans genoma hizalama.
Bunlar, bu kursun 2. Bölümünde bir Nextflow iş akışına saracağımız komutlardır.

1. FastQC kullanarak bir FASTQ dosyası üzerinde ilk QC'yi çalıştırın
2. Trim Galore kullanarak adaptör dizilerini kırpın ve kırpma sonrası QC'yi çalıştırın
3. Kırpılmış okumaları HISAT2 kullanarak referans genoma hizalayın

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-02.svg"
</figure>

Bu komutları yalnızca bir örnek üzerinde test ederek başlıyoruz.

### 1.1. QC ve adaptör kırpma

İlk olarak, örnek veri dosyalarından biri üzerinde QC ve kırpma komutlarını çalıştırmak istiyoruz.

#### 1.1.1. Konteyneri çekin

Hem `fastqc` hem de `trim_galore` yüklü bir konteyner imajını çekelim:

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

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

Bu imajı daha önce indirmediyseniz, tamamlanması bir dakika sürebilir.
Tamamlandığında, konteyner imajının yerel bir kopyasına sahip olursunuz.

#### 1.1.2. Konteyneri etkileşimli olarak başlatın

Konteyneri etkileşimli olarak çalıştırmak için `-it` bayraklarıyla `docker run` kullanın.
`-v ./data:/data` seçeneği, yerel `data/` dizinimizi bağlayarak girdi dosyalarına konteyner içinden erişmemizi sağlar.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Komut çıktısı"

    ```console
    (base) root@b645838b3314:/tmp#
    ```

İsteminiz `(base) root@b645838b3314:/tmp#` gibi bir şeye dönüşecek; bu, artık konteynerin içinde olduğunuzu gösterir.

Dizi veri dosyalarını `/data/reads` altında görebildiğinizi doğrulayın:

```bash
ls /data/reads
```

??? abstract "Dizin içeriği"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gz
    ```

Bununla, ilk komutunuzu denemeye hazırsınız.

#### 1.1.3. FastQC komutunu çalıştırın

Yukarıda referans verilen yöntem bize tek bir dosya üzerinde QC çalıştırmak için komut satırını verir.
Yalnızca girdi dosyasını sağlamamız gerekir; araç otomatik olarak orijinal verilerle aynı dizinde çıktı dosyaları oluşturacaktır.

Bir veri dosyası üzerinde `fastqc` komutunu çalıştırın:

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

??? abstract "Dizin içeriği"

    ```console
    /data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
    ```

Bir HTML raporu ve QC metriklerini içeren bir ZIP arşivi görmelisiniz.
Bu, ilk adımın testini tamamlar.

#### 1.1.4. Trim Galore ile adaptör dizilerini kırpın

Şimdi adaptör dizilerini kırpmak ve kırpma sonrası QC metriklerini toplamak için Cutadapt ve FastQC'yi bir araya getiren `trim_galore`'ı çalıştıralım.
Yukarıda belirtildiği gibi, yazılım aynı konteynere dahildir; bu nedenle orada herhangi bir değişiklik gerekmez.

Komut basittir; kırpma tamamlandıktan sonra otomatik olarak bir QC toplama adımı çalıştırmak için `--fastqc` bayrağını eklememiz yeterlidir.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Komut çıktısı"

    ```console hl_lines="54 55 56 58 59 60"
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)



    AUTO-DETECTING ADAPTER TYPE
    ===========================
    Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> /data/reads/ENCSR000COQ1_1.fastq.gz <<)

    Found perfect matches for the following adapter sequences:
    Adapter type	Count	Sequence	Sequences analysed	Percentage
    Illumina	9	AGATCGGAAGAGC	27816	0.03
    smallRNA	0	TGGAATTCTCGG	27816	0.00
    Nextera	0	CTGTCTCTTATA	27816	0.00
    Using Illumina adapter for trimming (count: 9). Second best hit was smallRNA (count: 0)

    Writing report to 'ENCSR000COQ1_1.fastq.gz_trimming_report.txt'

    SUMMARISING RUN PARAMETERS
    ==========================
    Input filename: /data/reads/ENCSR000COQ1_1.fastq.gz
    Trimming mode: single-end
    Trim Galore version: 0.6.10
    Cutadapt version: 4.9
    Number of cores used for trimming: 1
    Quality Phred score cutoff: 20
    Quality encoding type selected: ASCII+33
    Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    Maximum trimming error rate: 0.1 (default)
    Minimum required adapter overlap (stringency): 1 bp
    Minimum required sequence length before a sequence gets removed: 20 bp
    Running FastQC on the data once trimming has completed
    Output file(s) will be GZIP compressed

    Cutadapt seems to be fairly up-to-date (version 4.9). Setting -j 1
    Writing final adapter and quality trimmed output to ENCSR000COQ1_1_trimmed.fq.gz


      >>> Now performing quality (cutoff '-q 20') and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file /data/reads/ENCSR000COQ1_1.fastq.gz <<<
    This is cutadapt 4.9 with Python 3.12.7
    Command line parameters: -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC /data/reads/ENCSR000COQ1_1.fastq.gz
    Processing single-end reads on 1 core ...
    Finished in 0.373 s (13.399 µs/read; 4.48 M reads/minute).

    === Summary ===

    Total reads processed:                  27,816
    Reads with adapters:                     9,173 (33.0%)
    Reads written (passing filters):        27,816 (100.0%)

    Total basepairs processed:     2,114,016 bp
    Quality-trimmed:                       0 bp (0.0%)
    Total written (filtered):      2,100,697 bp (99.4%)

    === Adapter 1 ===

    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 9173 times

    Minimum overlap: 1
    No. of allowed errors:
    1-9 bp: 0; 10-13 bp: 1

    Bases preceding removed adapters:
      A: 27.4%
      C: 37.4%
      G: 20.9%
      T: 14.3%
      none/other: 0.0%

    Overview of removed sequences
    length	count	expect	max.err	error counts
    1	6229	6954.0	0	6229
    2	2221	1738.5	0	2221
    3	581	434.6	0	581
    4	88	108.7	0	88
    5	33	27.2	0	33
    6	2	6.8	0	2
    7	1	1.7	0	1
    9	1	0.1	0	1
    10	2	0.0	1	2
    12	1	0.0	1	0 1
    14	4	0.0	1	3 1
    16	1	0.0	1	1
    19	1	0.0	1	1
    22	1	0.0	1	1
    29	4	0.0	1	0 4
    33	3	0.0	1	3

    RUN STATISTICS FOR INPUT FILE: /data/reads/ENCSR000COQ1_1.fastq.gz
    =============================================
    27816 sequences processed in total
    Sequences removed because they became shorter than the length cutoff of 20 bp:	0 (0.0%)


      >>> Now running FastQC on the data <<<

    application/gzip
    Started analysis of ENCSR000COQ1_1_trimmed.fq.gz
    Approx 5% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 10% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 15% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 20% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 25% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 30% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 35% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 40% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 45% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 50% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 55% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 60% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 65% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 70% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 75% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 80% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 85% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 90% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 95% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

Çıktı çok ayrıntılı olduğundan, yukarıdaki örnekte en ilgili satırları vurguladık.
Çıktı dosyalarını çalışma dizininde bulabilirsiniz:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Dizin içeriği"

    ```console
    ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
    ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Bu, kırpılmış okumaları, kırpma raporunu ve kırpma sonrası QC dosyalarını içerir.

#### 1.1.5. Çıktı dosyalarını taşıyın

Konteyner içinde kalan her şey gelecekteki çalışmalara erişilemez olacağından, bu dosyaları bağlı dosya sistemindeki bir dizine taşımamız gerekir.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

??? abstract "Dizin içeriği"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    └── ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Dosyalar artık normal dosya sisteminizde erişilebilir.

#### 1.1.6. Konteynerden çıkın

Konteynerden çıkmak için `exit` yazın.

```bash
exit
```

İsteminiz normale dönmelidir; bu, ilk iki adımın testini tamamlar.

### 1.2. Okumaları referans genoma hizalayın

Ardından, kırpılmış RNAseq okumalarını bir referans genoma hizalamak için hizalama komutunu çalıştırmak istiyoruz.

#### 1.2.1. Konteyneri çekin

`hisat2` ve `samtools` yüklü bir konteyner imajını çekelim:

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

Bazı katmanların `Already exists` gösterdiğini fark edeceksiniz; bunun nedeni, bu katmanların daha önce çektiğimiz Trim Galore konteyner imajıyla paylaşılmasıdır.
Sonuç olarak, bu çekme işlemi ilkinden daha hızlı olmalıdır.

#### 1.2.2. Konteyneri etkileşimli olarak başlatın

Konteyneri etkileşimli olarak başlatın; ilgili konteyner URI'si değiştirilerek daha öncekiyle aynı yaklaşımı kullanın.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

İsteminiz, konteynerin içinde olduğunuzu belirtmek için tekrar değişecektir.

#### 1.2.3. Genom dizin dosyalarını oluşturun

HISAT2, genom referansının çok özel bir biçimde sağlanmasını gerektirir ve sağladığımız `genome.fa` FASTA dosyasını doğrudan kullanamaz; bu nedenle bu fırsatı ilgili kaynakları oluşturmak için kullanacağız.

```bash
hisat2-build /data/genome.fa genome_index
```

??? success "Komut çıktısı"

    ```console hl_lines="1 2 218"
    Settings:
      Output files: "genome_index.*.ht2"
      Line rate: 6 (line is 64 bytes)
      Lines per side: 1 (side is 64 bytes)
      Offset rate: 4 (one in 16)
      FTable chars: 10
      Strings: unpacked
      Local offset rate: 3 (one in 8)
      Local fTable chars: 6
      Local sequence length: 57344
      Local sequence overlap between two consecutive indexes: 1024
      Endianness: little
      Actual local endianness: little
      Sanity checking: disabled
      Assertions: disabled
      Random seed: 0
      Sizeofs: void*:8, int:4, long:8, size_t:8
    Input files DNA, FASTA:
      /data/genome.fa
    Reading reference sizes
      Time reading reference sizes: 00:00:00
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:00
      Time to read SNPs and splice sites: 00:00:00
    Using parameters --bmax 6542727 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 6542727 --dcv 1024
    Constructing suffix-array element generator
    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:00:01
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:00
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:00:00
      Sanity-checking and returning
    Building samples
    Reserving space for 12 sample suffixes
    Generating random suffixes
    QSorting 12 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 12 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 1, merged 7; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 4.98493e+06 (target: 6542726)
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering GFM loop
    Getting block 1 of 7
      Reserving size (6542727) for bucket 1
      Calculating Z arrays for bucket 1
      Entering block accumulator loop for bucket 1:
      bucket 1: 10%
      bucket 1: 20%
      bucket 1: 30%
      bucket 1: 40%
      bucket 1: 50%
      bucket 1: 60%
      bucket 1: 70%
      bucket 1: 80%
      bucket 1: 90%
      bucket 1: 100%
      Sorting block of length 3540952 for bucket 1
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 3540953 for bucket 1
    Getting block 2 of 7
      Reserving size (6542727) for bucket 2
      Calculating Z arrays for bucket 2
      Entering block accumulator loop for bucket 2:
      bucket 2: 10%
      bucket 2: 20%
      bucket 2: 30%
      bucket 2: 40%
      bucket 2: 50%
      bucket 2: 60%
      bucket 2: 70%
      bucket 2: 80%
      bucket 2: 90%
      bucket 2: 100%
      Sorting block of length 6195795 for bucket 2
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6195796 for bucket 2
    Getting block 3 of 7
      Reserving size (6542727) for bucket 3
      Calculating Z arrays for bucket 3
      Entering block accumulator loop for bucket 3:
      bucket 3: 10%
      bucket 3: 20%
      bucket 3: 30%
      bucket 3: 40%
      bucket 3: 50%
      bucket 3: 60%
      bucket 3: 70%
      bucket 3: 80%
      bucket 3: 90%
      bucket 3: 100%
      Sorting block of length 6199288 for bucket 3
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6199289 for bucket 3
    Getting block 4 of 7
      Reserving size (6542727) for bucket 4
      Calculating Z arrays for bucket 4
      Entering block accumulator loop for bucket 4:
      bucket 4: 10%
      bucket 4: 20%
      bucket 4: 30%
      bucket 4: 40%
      bucket 4: 50%
      bucket 4: 60%
      bucket 4: 70%
      bucket 4: 80%
      bucket 4: 90%
      bucket 4: 100%
      Sorting block of length 6454986 for bucket 4
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 6454987 for bucket 4
    Getting block 5 of 7
      Reserving size (6542727) for bucket 5
      Calculating Z arrays for bucket 5
      Entering block accumulator loop for bucket 5:
      bucket 5: 10%
      bucket 5: 20%
      bucket 5: 30%
      bucket 5: 40%
      bucket 5: 50%
      bucket 5: 60%
      bucket 5: 70%
      bucket 5: 80%
      bucket 5: 90%
      bucket 5: 100%
      Sorting block of length 3493181 for bucket 5
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3493182 for bucket 5
    Getting block 6 of 7
      Reserving size (6542727) for bucket 6
      Calculating Z arrays for bucket 6
      Entering block accumulator loop for bucket 6:
      bucket 6: 10%
      bucket 6: 20%
      bucket 6: 30%
      bucket 6: 40%
      bucket 6: 50%
      bucket 6: 60%
      bucket 6: 70%
      bucket 6: 80%
      bucket 6: 90%
      bucket 6: 100%
      Sorting block of length 5875908 for bucket 6
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 5875909 for bucket 6
    Getting block 7 of 7
      Reserving size (6542727) for bucket 7
      Calculating Z arrays for bucket 7
      Entering block accumulator loop for bucket 7:
      bucket 7: 10%
      bucket 7: 20%
      bucket 7: 30%
      bucket 7: 40%
      bucket 7: 50%
      bucket 7: 60%
      bucket 7: 70%
      bucket 7: 80%
      bucket 7: 90%
      bucket 7: 100%
      Sorting block of length 3134429 for bucket 7
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3134430 for bucket 7
    Exited GFM loop
    fchr[A]: 0
    fchr[C]: 9094775
    fchr[G]: 17470759
    fchr[T]: 25839994
    fchr[$]: 34894545
    Exiting GFM::buildToDisk()
    Returning from initFromVector
    Wrote 15826295 bytes to primary GFM file: genome_index.1.ht2
    Wrote 8723644 bytes to secondary GFM file: genome_index.2.ht2
    Re-opening _in1 and _in2 as input streams
    Returning from GFM constructor
    Returning from initFromVector
    Wrote 15353415 bytes to primary GFM file: genome_index.5.ht2
    Wrote 8883598 bytes to secondary GFM file: genome_index.6.ht2
    Re-opening _in5 and _in5 as input streams
    Returning from HGFM constructor
    Headers:
        len: 34894545
        gbwtLen: 34894546
        nodes: 34894546
        sz: 8723637
        gbwtSz: 8723637
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 0
        eftabSz: 0
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 2180910
        offsSz: 8723640
        lineSz: 64
        sideSz: 64
        sideGbwtSz: 48
        sideGbwtLen: 192
        numSides: 181743
        numLines: 181743
        gbwtTotLen: 11631552
        gbwtTotSz: 11631552
        reverse: 0
        linearFM: Yes
    Total time for call to driver() for forward index: 00:00:12
    ```

Çıktı çok ayrıntılı olduğundan, yukarıdaki örnekte bazı ilgili satırları vurguladık.

Bu, çalışma dizininde bulabileceğiniz birden fazla genom dizin dosyası oluşturur.

```bash
ls genome_index.*
```

??? abstract "Dizin içeriği"

    ```console
    genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
    genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
    ```

Bu dosyalara daha sonra ihtiyacımız olacak; bunları oluşturmak tipik olarak bir iş akışının parçası olarak yapmak istediğimiz bir şey değildir. Bu nedenle, gerektiğinde kolayca aktarabileceğimiz genom dizin dosyalarını içeren sıkıştırılmış bir tarball oluşturacağız.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

??? success "Komut çıktısı"

    ```console
    genome_index.1.ht2
    genome_index.2.ht2
    genome_index.3.ht2
    genome_index.4.ht2
    genome_index.5.ht2
    genome_index.6.ht2
    genome_index.7.ht2
    genome_index.8.ht2
    ```

Genom dizin dosyalarını içeren ortaya çıkan `genome_index.tar.gz` tarball'ını birkaç dakika içinde dosya sistemimizde `data/` dizinine taşıyacağız.
Bu, bu kursun 2. Bölümünde işimize yarayacaktır.

#### 1.2.4. Hizalama komutunu çalıştırın

Şimdi hizalama komutunu çalıştırabiliriz; bu komut `hisat2` ile hizalama adımını gerçekleştirir ve ardından çıktıyı bir BAM dosyası olarak yazmak için `samtools`'a aktarır.

Okuma verisi girdisi, önceki adımda `trim_galore` ile oluşturduğumuz `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` dosyasıdır.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Komut çıktısı"

    ```console hl_lines="6"
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

??? abstract "Dizin içeriği"

    ```console title="Çıktı"
    ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
    ```

Hizalama bir BAM dosyası ve hizalama istatistiklerini içeren bir log dosyası üretti.

#### 1.2.5. Çıktı dosyalarını taşıyın

Daha önce olduğu gibi, çıktı dosyalarını bağlı dosya sistemindeki bir dizine taşıyın; böylece konteynerden çıktıktan sonra erişilebilir kalırlar.

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

Bu yapıldıktan sonra, ihtiyacımız olan her şeye sahibiz.

#### 1.2.6. Konteynerden çıkın

Konteynerden çıkmak için `exit` yazın.

```bash
exit
```

İsteminiz normale dönmelidir.
Bu, tek örnekli işleme test çalışmasını sonlandırır.

!!! example "Bunu bir iş akışı olarak yazın!"

    Bu analizi bir Nextflow iş akışı olarak uygulamaya başlamak isterseniz hemen [Bölüm 2](./02_single-sample.md)'ye geçebilirsiniz.
    Sadece Bölüm 3'e geçmeden önce ikinci test turunu tamamlamak için geri dönmeniz gerekecek.

---

## 2. Çok örnekli QC toplama

Az önce test ettiğimiz komutlar bir seferde bir örneği işler.
Pratikte, tipik olarak birçok örneği işlememiz ve ardından genel veri setinin kalitesini değerlendirmek için tüm bunlar arasında QC sonuçlarını toplamamız gerekir.

[MultiQC](https://multiqc.info/), birçok yaygın biyoinformatik aracından QC raporları için dizinlerde arama yapan ve bunları tek bir kapsamlı HTML raporunda toplayan bir araçtır.
FastQC, Cutadapt (Trim Galore aracılığıyla) ve HISAT2'den çıktıları diğerleri arasında tanıyabilir.

Burada aynı örnek başına araçlarla iki ek örneği işliyoruz, ardından üç örnek arasında QC raporlarını toplamak için MultiQC kullanıyoruz.
Bunlar, bu kursun 3. Bölümünde bir Nextflow iş akışına saracağımız komutlardır.

1. Trim Galore kullanarak ek örnekler üzerinde QC ve kırpma çalıştırın
2. HISAT2 kullanarak ek örnekler üzerinde hizalama çalıştırın
3. MultiQC kullanarak tüm QC raporlarını kapsamlı bir raporda toplayın

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

### 2.1. Ek örnekleri QC ve kırpma

Örnek başına QC ve kırpma komutları, bölüm 1.1'de çalıştırdığımızla aynıdır.
Konteyner imajını zaten çektik; bu nedenle doğrudan başlatabiliriz.

#### 2.1.1. Konteyneri başlatın

Bu konteyner imajını bölüm 1.1'de zaten çektik; bu nedenle doğrudan başlatabiliriz:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

İsteminiz, konteynerin içinde olduğunuzu belirtmek için değişir.

#### 2.1.2. Ek örnekler üzerinde QC ve kırpma çalıştırın

İki örnek daha üzerinde FastQC ve Trim Galore'ı birbiri ardına çalıştırın.

```bash
fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
trim_galore --fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

Bu tamamlandığında, çalışma dizininde her iki örnek için de Trim Galore çıktı dosyalarına sahip olmalısınız.

#### 2.1.3. Çıktı dosyalarını taşıyın

Trim Galore çıktı dosyalarını bölüm 1'de kullandığımız aynı dizine taşıyın.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/trimmed
```

??? abstract "Dizin içeriği"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    ├── ENCSR000COQ1_1_trimmed_fastqc.zip
    ├── ENCSR000COQ2_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ2_1_trimmed.fq.gz
    ├── ENCSR000COQ2_1_trimmed_fastqc.html
    ├── ENCSR000COQ2_1_trimmed_fastqc.zip
    ├── ENCSR000COR1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COR1_1_trimmed.fq.gz
    ├── ENCSR000COR1_1_trimmed_fastqc.html
    └── ENCSR000COR1_1_trimmed_fastqc.zip
    ```

Dosyalar artık normal dosya sisteminizde erişilebilir.

#### 2.1.4. Konteynerden çıkın

Konteynerden çıkmak için `exit` yazın.

```bash
exit
```

İsteminiz normale dönmelidir.

### 2.2. Ek örnekleri hizalayın

Hizalama komutları, bölüm 1.2'de çalıştırdığımızla aynıdır.
Orijinal dizin dosyaları artık var olmayan bir konteyner içinde oluşturulduğundan, daha önce kaydettiğimiz tarball'dan genom dizinini çıkarmamız gerekir.

#### 2.2.1. Konteyneri başlatın

Bu konteyner imajını bölüm 1.2'de zaten çektik; bu nedenle doğrudan başlatabiliriz:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

İsteminiz, konteynerin içinde olduğunuzu belirtmek için değişir.

#### 2.2.2. Genom dizinini çıkarın

Bağlı dosya sistemine kaydettiğimiz tarball'dan genom dizin dosyalarını çıkarın:

```bash
tar -xzf /data/genome_index.tar.gz
```

Bu, çalışma dizininde `genome_index.*` dosyalarını geri yükler.

#### 2.2.3. Ek örnekler üzerinde hizalama çalıştırın

Yeni kırpılmış iki örnek üzerinde HISAT2 hizalamasını birbiri ardına çalıştırın.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ2_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ2_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ2_1_trimmed.bam
```

??? success "Komut çıktısı"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 18736
    		Aligned 0 time: 1531 (8.17%)
    		Aligned 1 time: 16726 (89.27%)
    		Aligned >1 times: 479 (2.56%)
    	Overall alignment rate: 91.83%
    ```

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COR1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COR1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COR1_1_trimmed.bam
```

??? success "Komut çıktısı"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 38056
    		Aligned 0 time: 2311 (6.07%)
    		Aligned 1 time: 33289 (87.47%)
    		Aligned >1 times: 2456 (6.45%)
    	Overall alignment rate: 93.93%
    ```

Bu tamamlandığında, çalışma dizininde her iki örnek için de BAM ve log dosyalarına sahip olmalısınız.

#### 2.2.4. Çıktı dosyalarını taşıyın

Hizalama çıktı dosyalarını bölüm 1'de kullandığımız aynı dizine taşıyın.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/aligned
```

??? abstract "Dizin içeriği"

    ```console
    /data/aligned
    ├── ENCSR000COQ1_1_trimmed.bam
    ├── ENCSR000COQ1_1_trimmed.hisat2.log
    ├── ENCSR000COQ2_1_trimmed.bam
    ├── ENCSR000COQ2_1_trimmed.hisat2.log
    ├── ENCSR000COR1_1_trimmed.bam
    └── ENCSR000COR1_1_trimmed.hisat2.log
    ```

Dosyalar artık normal dosya sisteminizde erişilebilir.

#### 2.2.5. Konteynerden çıkın

Konteynerden çıkmak için `exit` yazın.

```bash
exit
```

İsteminiz normale dönmelidir.

### 2.3. Kapsamlı bir QC raporu oluşturun

Artık üç örnek için QC, kırpma ve hizalama çıktısına sahip olduğumuza göre, bunları tek bir raporda toplamak için MultiQC kullanabiliriz.
MultiQC, uyumlu QC raporları için dizinlerde arama yapar ve bulduğu her şeyi toplar.

#### 2.3.1. Konteyneri çekin

`multiqc` yüklü bir konteyner imajını çekelim:

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Komut çıktısı"

    ```console
    a3c26f6199d64b7c: Pulling from library/pip_multiqc
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
    2ed162b168e8: Pull complete
    ca06fe148f21: Pull complete
    Digest: sha256:af0e9de56896805aa2a065f7650362956f4213d99e95314f6fec472c6a3bf091
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

Bazı katmanların `Already exists` gösterdiğini fark edeceksiniz; bunun nedeni, bu katmanların daha önce çektiğimiz konteyner imajlarıyla paylaşılmasıdır.
Sonuç olarak, bu çekme işlemi öncekilerden daha hızlı olmalıdır.

#### 2.3.2. Konteyneri etkileşimli olarak başlatın

Konteyneri, daha önce olduğu gibi veri dizini bağlı olarak etkileşimli olarak başlatın.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

İsteminiz, konteynerin içinde olduğunuzu belirtmek için değişecektir.

#### 2.3.3. MultiQC komutunu çalıştırın

`multiqc`'yi çalıştırın; üç örneğin tümü için QC ile ilgili çıktı dosyalarını sakladığımız dizinlere işaret edin.
`-n` bayrağı çıktı raporunun adını ayarlar.

```bash
multiqc /data/reads /data/trimmed /data/aligned -n all_samples_QC
```

??? success "Komut çıktısı"

    ```console hl_lines="8 9 10 11 12"

    /// MultiQC 🔍 v1.32

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 36/36
               hisat2 | Found 3 reports
             cutadapt | Found 3 reports
               fastqc | Found 3 reports
        write_results | Data        : all_samples_QC_data
        write_results | Report      : all_samples_QC.html
              multiqc | MultiQC complete
    ```

Burada aracın üç örneğin tümü için QC raporlarını bulduğunu görüyoruz: `fastqc` ile yaptığımız ilk QC, `cutadapt`'ten (`trim_galore` aracılığıyla) kırpma sonrası raporlar ve `hisat2` tarafından üretilen hizalama özetleri.

Çıktı dosyaları çalışma dizinindedir:

```bash
ls all_samples_QC*
```

??? abstract "Dizin içeriği"

    ```console
    all_samples_QC.html

    all_samples_QC_data:
    cutadapt_filtered_reads_plot.txt                     multiqc.log
    cutadapt_trimmed_sequences_plot_3_Counts.txt         multiqc.parquet
    cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc_citations.txt
    fastqc-status-check-heatmap.txt                      multiqc_cutadapt.txt
    fastqc_adapter_content_plot.txt                      multiqc_data.json
    fastqc_overrepresented_sequences_plot.txt            multiqc_fastqc.txt
    fastqc_per_base_n_content_plot.txt                   multiqc_general_stats.txt
    fastqc_per_base_sequence_quality_plot.txt            multiqc_hisat2.txt
    fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_software_versions.txt
    fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_sources.txt
    fastqc_per_sequence_quality_scores_plot.txt
    fastqc_sequence_counts_plot.txt
    fastqc_sequence_duplication_levels_plot.txt
    fastqc_top_overrepresented_sequences_table.txt
    hisat2_se_plot.txt
    llms-full.txt
    ```

Ana çıktı, temel metrikleri içeren bir veri dizini ile birlikte `all_samples_QC.html` raporudur.

#### 2.3.4. Çıktı dosyalarını taşıyın

Raporu ve veri dizinini bağlı dosya sistemine taşıyın.

```bash
mkdir /data/multiqc
mv all_samples_QC* /data/multiqc
```

Dosyalar artık normal dosya sisteminizde erişilebilir.

#### 2.3.5. Konteynerden çıkın

Konteynerden çıkmak için `exit` yazın.

```bash
exit
```

İsteminiz normale dönmelidir.
Bu, tüm RNAseq işleme komutlarının testini sonlandırır.

---

### Özet

FastQC, Trim Galore, HISAT2 ve MultiQC komutlarını ilgili konteynerlerinde nasıl çalıştıracağınızı biliyorsunuz; birden fazla örneği işlemeyi ve QC raporlarını toplamayı da içerir.

### Sırada ne var?

Bir mola verin, ardından aynı komutları çalışmayı yürütmek için konteynerler kullanan iş akışlarına nasıl saracağınızı öğrenmek için [Bölüm 2](./02_single-sample.md)'ye geçin.
