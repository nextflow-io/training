# Bölüm 1: Yöntemlere genel bakış ve manuel test

Varyant çağırma (variant calling), bir genom dizisindeki varyasyonları bir referans genoma göre tanımlamayı amaçlayan bir genomik analiz yöntemidir.
Burada, tüm genom dizileme verilerinde kısa germline varyantları, _yani_ SNP'leri ve indelleri çağırmak için tasarlanmış araçları ve yöntemleri kullanacağız.

![GATK boru hattı](img/gatk-pipeline.png)

Tam bir varyant çağırma boru hattı tipik olarak referansa haritalama (bazen genom hizalaması olarak da adlandırılır) ve varyant filtreleme ve önceliklendirme dahil olmak üzere birçok adım içerir.
Basitlik için, bu kursta yalnızca varyant çağırma kısmına odaklanacağız.

### Yöntemler

Size germline SNP'leri ve indelleri tanımlamak için tüm genom dizileme örneklerine varyant çağırmayı uygulamanın iki yolunu göstereceğiz.
İlk olarak, her örnekten bağımsız olarak varyantları çağıran basit bir **örnek başına yaklaşımla** başlayacağız.
Ardından, birden fazla örneği birlikte analiz ederek daha doğru ve bilgilendirici sonuçlar üreten daha sofistike bir **ortak çağırma yaklaşımını** göstereceğiz.

Her iki yaklaşım için herhangi bir iş akışı kodu yazmaya dalmadan önce, komutları bazı test verileri üzerinde manuel olarak deneyeceğiz.

### Veri seti

Aşağıdaki veri ve ilgili kaynakları sağlıyoruz:

- İnsan kromozom 20'sinin küçük bir bölgesinden (hg19/b37'den) oluşan **bir referans genom** ve yardımcı dosyaları (indeks ve dizin sözlüğü).
- Bir aile üçlüsüne (anne, baba ve oğul) karşılık gelen **üç tüm genom dizileme örneği**, bunlar dosya boyutlarını küçük tutmak için kromozom 20'deki küçük bir veri dilimine indirilmiştir.
  Bu, zaten referans genoma haritalanmış Illumina kısa okuma dizileme verisidir ve [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) formatında sağlanmıştır (Binary Alignment Map, SAM'ın sıkıştırılmış bir versiyonu, Sequence Alignment Map).
- Örneklerimizin varyantları çağırmak için uygun veri içerdiği genomdaki koordinatlar olan **bir genomik aralıklar listesi**, BED formatında sağlanmıştır.

### Yazılım

İlgili iki ana araç, dizilim hizalama dosyalarını işlemek için yaygın olarak kullanılan bir araç seti olan [Samtools](https://www.htslib.org/) ve Broad Institute'da geliştirilen varyant keşfi için bir dizi araç olan [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit)'dir.

Bu araçlar GitHub Codespaces ortamında yüklü değildir, bu yüzden onları konteynerlar aracılığıyla kullanacağız (bkz. [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Not"

    `nf4-science/genomics` dizininde olduğunuzdan emin olun, böylece `pwd` yazdığınızda gösterilen yolun son kısmı `genomics` olmalıdır.

---

## 1. Örnek başına varyant çağırma

Örnek başına varyant çağırma, her örneği bağımsız olarak işler: varyant çağırıcı bir seferde bir örnek için dizileme verilerini inceler ve örneğin referanstan farklı olduğu konumları tanımlar.

Bu bölümde, örnek başına varyant çağırma yaklaşımını oluşturan iki komutu test ediyoruz: Samtools ile bir BAM dosyasını indeksleme ve GATK HaplotypeCaller ile varyantları çağırma.
Bunlar, bu kursun 2. Bölümünde bir Nextflow iş akışına sarmalayacağımız komutlardır.

1. [Samtools](https://www.htslib.org/) kullanarak bir BAM girdi dosyası için bir indeks dosyası oluşturun
2. VCF (Variant Call Format) formatında örnek başına varyant çağrıları oluşturmak için indekslenmiş BAM dosyası üzerinde GATK HaplotypeCaller'ı çalıştırın

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

İki komutu yalnızca bir örnek üzerinde test ederek başlıyoruz.

### 1.1. Samtools ile bir BAM girdi dosyasını indeksleyin

İndeks dosyaları, biyoenformatik dosya formatlarının yaygın bir özelliğidir; GATK gibi araçların dosyanın tamamını okumak zorunda kalmadan verilerin bir alt kümesine erişmesini sağlayan ana dosyanın yapısı hakkında bilgi içerirler.
Bu, bu dosyaların ne kadar büyük olabileceği nedeniyle önemlidir.

BAM dosyaları genellikle bir indeks olmadan sağlanır, bu nedenle birçok analiz iş akışındaki ilk adım `samtools index` kullanarak bir tane oluşturmaktır.

Bir Samtools konteynırı çekeceğiz, etkileşimli olarak çalıştıracağız ve BAM dosyalarından biri üzerinde `samtools index` komutunu çalıştıracağız.

#### 1.1.1. Samtools konteynırını çekin

Samtools konteyner imajını indirmek için `docker pull` komutunu çalıştırın:

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

??? success "Komut çıktısı"

    ```console
    1.20--b5dfbd93de237464: Pulling from library/samtools
    6360b3717211: Pull complete
    2ec3f7ad9b3c: Pull complete
    7716ca300600: Pull complete
    4f4fb700ef54: Pull complete
    8c61d418774c: Pull complete
    03dae77ff45c: Pull complete
    aab7f787139d: Pull complete
    4f4fb700ef54: Pull complete
    837d55536720: Pull complete
    897362c12ca7: Pull complete
    3893cbe24e91: Pull complete
    d1b61e94977b: Pull complete
    c72ff66fb90f: Pull complete
    0e0388f29b6d: Pull complete
    Digest: sha256:bbfc45b4f228975bde86cba95e303dd94ecf2fdacea5bfb2e2f34b0d7b141e41
    Status: Downloaded newer image for community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    ```

Bu imajı daha önce indirmediyseniz, tamamlanması bir dakika sürebilir.
Tamamlandıktan sonra, konteyner imajının yerel bir kopyasına sahip olursunuz.

#### 1.1.2. Samtools konteynırını etkileşimli olarak çalıştırın

Konteynırı etkileşimli olarak çalıştırmak için `-it` bayraklarıyla `docker run` kullanın.
`-v ./data:/data` seçeneği, araçların girdi dosyalarına erişebilmesi için yerel `data` dizinini konteynıra bağlar.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

İsteğiniz `(base) root@a1b2c3d4e5f6:/tmp#` gibi bir şeye dönüşür, bu da artık konteyner içinde olduğunuzu gösterir.
Veri dosyalarına `/data` altından erişilebilir.

#### 1.1.3. İndeksleme komutunu çalıştırın

[Samtools belgeleri](https://www.htslib.org/doc/samtools-index.html) bize bir BAM dosyasını indekslemek için çalıştırmamız gereken komut satırını verir.

Yalnızca girdi dosyasını sağlamamız gerekir; araç, girdi dosya adına `.bai` ekleyerek çıktı için otomatik olarak bir ad oluşturacaktır.

```bash
samtools index /data/bam/reads_mother.bam
```

??? abstract "Dizin içerikleri"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

Şimdi orijinal BAM girdi dosyasıyla aynı dizinde `reads_mother.bam.bai` adlı bir dosya görmelisiniz.

#### 1.1.4. Samtools konteynırından çıkın

Konteynerdan çıkmak için `exit` yazın.

```bash
exit
```

İsteğiniz artık konteynırı başlatmadan önceki haline dönmüş olmalıdır.

### 1.2. GATK HaplotypeCaller ile varyantları çağırın

Bir GATK konteynırı çekeceğiz, etkileşimli olarak çalıştıracağız ve az önce indekslediğimiz BAM dosyası üzerinde `gatk HaplotypeCaller` komutunu çalıştıracağız.

#### 1.2.1. GATK konteynırını çekin

GATK konteyner imajını indirmek için `docker pull` komutunu çalıştırın:

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

??? success "Komut çıktısı"

    Bazı katmanlar `Already exists` gösterir çünkü daha önce çektiğimiz Samtools konteyner imajıyla paylaşılıyorlar.

    ```console
    4.5.0.0--730ee8817e436867: Pulling from library/gatk4
    6360b3717211: Already exists
    2ec3f7ad9b3c: Already exists
    7716ca300600: Already exists
    4f4fb700ef54: Already exists
    8c61d418774c: Already exists
    03dae77ff45c: Already exists
    aab7f787139d: Already exists
    4f4fb700ef54: Already exists
    837d55536720: Already exists
    897362c12ca7: Already exists
    3893cbe24e91: Already exists
    d1b61e94977b: Already exists
    e5c558f54708: Pull complete
    087cce32d294: Pull complete
    Digest: sha256:e33413b9100f834fcc62fd5bc9edc1e881e820aafa606e09301eac2303d8724b
    Status: Downloaded newer image for community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    ```

Bu, ilk çekişten daha hızlı olmalıdır çünkü iki konteyner imajı katmanlarının çoğunu paylaşır.

#### 1.2.2. GATK konteynırını etkileşimli olarak çalıştırın

GATK konteynırını, Samtools için yaptığımız gibi veri dizini bağlanmış şekilde etkileşimli olarak çalıştırın.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

İsteğiniz artık GATK konteynerinin içinde olduğunuzu gösterecek şekilde değişir.

#### 1.2.3. Varyant çağırma komutunu çalıştırın

[GATK belgeleri](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) bize bir BAM dosyası üzerinde varyant çağırma gerçekleştirmek için çalıştırmamız gereken komut satırını verir.

BAM girdi dosyasını (`-I`) ve referans genomu (`-R`), çıktı dosyası için bir ad (`-O`) ve analiz edilecek genomik aralıkların bir listesini (`-L`) sağlamamız gerekir.

Ancak, indeks dosyasının yolunu belirtmemize gerek yoktur; araç, yerleşik adlandırma ve birlikte bulundurma kuralına dayanarak aynı dizinde otomatik olarak arayacaktır.
Aynı şey referans genomun yardımcı dosyaları (indeks ve dizin sözlüğü dosyaları, `*.fai` ve `*.dict`) için de geçerlidir.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

??? success "Komut çıktısı"

    Araç ayrıntılı günlük çıktısı üretir. Vurgulanan satırlar başarılı tamamlanmayı onaylar.

    ```console hl_lines="37 51 56 57"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.vcf -L /data/ref/intervals.bed
    00:27:50.687 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:27:50.854 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.858 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:27:50.858 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:27:50.858 INFO  HaplotypeCaller - Executing as root@a1fe8ff42d07 on Linux v6.10.14-linuxkit amd64
    00:27:50.858 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:27:50.859 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:27:50 AM GMT
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.861 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:27:50.861 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:27:50.861 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:27:50.863 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:27:50.864 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:27:50.864 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:27:50.864 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:27:50.864 INFO  HaplotypeCaller - Requester pays: disabled
    00:27:50.865 INFO  HaplotypeCaller - Initializing engine
    00:27:50.991 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:27:51.016 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:27:51.029 INFO  HaplotypeCaller - Done initializing engine
    00:27:51.040 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:27:51.042 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:27:51.042 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:27:51.046 INFO  HaplotypeCallerEngine - Disabling physical phasing, which is supported only for reference-model confidence output
    00:27:51.063 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:27:51.085 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:27:51.086 INFO  IntelPairHmm - Available threads: 10
    00:27:51.086 INFO  IntelPairHmm - Requested threads: 4
    00:27:51.086 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:27:51.128 INFO  ProgressMeter - Starting traversal
    00:27:51.136 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:27:51.882 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:27:52.969 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:27:52.971 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35           1145.7
    00:27:52.971 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:27:52.976 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003346916
    00:27:52.976 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.045731709
    00:27:52.977 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:27:52.981 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:27:52 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.04 minutes.
    Runtime.totalMemory()=203423744
    ```

Çıktı dosyası `reads_mother.vcf`, konteyner içindeki çalışma dizininizde oluşturulur, bu nedenle çıktı dosyası yolunu değiştirmediğiniz sürece VS Code dosya gezgininde görmezsiniz.
Ancak, küçük bir test dosyasıdır, bu nedenle açmak ve içeriği görüntülemek için `cat` yapabilirsiniz.
Dosyanın başlangıcına kadar kaydırırsanız, birçok satır meta veriden oluşan bir başlık ve ardından satır başına bir varyant çağrısı listesi bulacaksınız.

??? abstract "Dosya içeriği"

    ```console title="reads_mother.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Her satır, örneğin dizileme verilerinde tanımlanan olası bir varyantı açıklar. VCF formatını yorumlamak için rehberlik için [bu faydalı makaleye](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/) bakın.

Çıktı VCF dosyasına, GATK tarafından otomatik olarak oluşturulan `reads_mother.vcf.idx` adlı bir indeks dosyası eşlik eder.
BAM indeks dosyasıyla aynı işleve sahiptir; araçların tüm dosyayı yüklemeden veri alt kümelerini aramasına ve almasına izin verir.

#### 1.2.4. GATK konteynırından çıkın

Konteynerdan çıkmak için `exit` yazın.

```bash
exit
```

İsteğiniz normale dönmüş olmalıdır.
Bu, örnek başına varyant çağırma testini sonlandırır.

---

## 2. Bir kohort üzerinde ortak çağırma

Az önce kullandığımız varyant çağırma yaklaşımı örnek başına varyant çağrıları üretir.
Bu, her örnekten izole edilmiş varyantlara bakmak için iyidir, ancak sınırlı bilgi verir.
Varyant çağrılarının birden fazla örnek arasında nasıl farklılık gösterdiğine bakmak genellikle daha ilginçtir.
GATK bu amaç için ortak varyant çağırma (joint variant calling) adı verilen alternatif bir yöntem sunar.

Ortak varyant çağırma, her örnek için GVCF (Genomic VCF) adı verilen özel bir varyant çıktısı oluşturmayı, ardından tüm örneklerden gelen GVCF verilerini birleştirmeyi ve bir 'ortak genotipleme' istatistiksel analizi çalıştırmayı içerir.

![Ortak analiz](img/joint-calling.png)

Bir örneğin GVCF'sinin özel yanı, yalnızca programın varyasyon kanıtı bulduğu konumları değil, genomun hedeflenen alanındaki tüm konumlar hakkındaki dizilim veri istatistiklerini özetleyen kayıtlar içermesidir.
Bu, ortak genotipleme hesaplaması için kritiktir ([daha fazla okuma](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

GVCF, ek bir parametreyle (`-ERC GVCF`) az önce test ettiğimiz aynı araç olan GATK HaplotypeCaller tarafından üretilir.
GVCF'leri birleştirmek, örnek başına çağrıları bir veri deposuna (bir veritabanına benzer) birleştiren GATK GenomicsDBImport ile yapılır.
Ardından gerçek 'ortak genotipleme' analizi GATK GenotypeGVCFs ile yapılır.

Burada GVCF'ler oluşturmak ve ortak genotipleme çalıştırmak için gereken komutları test ediyoruz.
Bunlar, bu kursun 3. Bölümünde bir Nextflow iş akışına sarmalayacağımız komutlardır.

1. Samtools kullanarak her BAM girdi dosyası için bir indeks dosyası oluşturun
2. Örnek başına genomik varyant çağrılarının bir GVCF'sini oluşturmak için her BAM girdi dosyası üzerinde GATK HaplotypeCaller'ı çalıştırın
3. Tüm GVCF'leri toplayın ve bunları bir GenomicsDB veri deposunda birleştirin
4. Kohort düzeyinde bir VCF üretmek için birleştirilmiş GVCF veri deposu üzerinde ortak genotipleme çalıştırın

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

Şimdi tüm üç BAM dosyasını indekslemeyle başlayarak bu komutların hepsini test etmemiz gerekiyor.

### 2.1. Üç örneğin tümü için BAM dosyalarını indeksleyin

Yukarıdaki ilk bölümde, yalnızca bir BAM dosyasını indeksledik.
Şimdi GATK HaplotypeCaller'ın bunları işleyebilmesi için üç örneğin tümünü indekslememiz gerekiyor.

#### 2.1.1. Samtools konteynırını etkileşimli olarak çalıştırın

Samtools konteyner imajını zaten çektik, bu yüzden doğrudan çalıştırabiliriz:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

İsteğiniz, daha önce olduğu gibi veri dizini bağlı şekilde konteyner içinde olduğunuzu gösterecek şekilde değişir.

#### 2.1.2. Üç örneğin tümünde indeksleme komutunu çalıştırın

Üç BAM dosyasının her birinde indeksleme komutunu çalıştırın:

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

??? abstract "Dizin içerikleri"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

Bu, indeks dosyalarını karşılık gelen BAM dosyalarıyla aynı dizinde üretmelidir.

#### 2.1.3. Samtools konteynırından çıkın

Konteynerdan çıkmak için `exit` yazın.

```bash
exit
```

İsteğiniz normale dönmüş olmalıdır.

### 2.2. Üç örneğin tümü için GVCF'ler oluşturun

Ortak genotipleme adımını çalıştırmak için üç örneğin tümü için GVCF'lere ihtiyacımız var.

#### 2.2.1. GATK konteynırını etkileşimli olarak çalıştırın

GATK konteyner imajını daha önce zaten çektik, bu yüzden doğrudan çalıştırabiliriz:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

İsteğiniz GATK konteynerinin içinde olduğunuzu gösterecek şekilde değişir.

#### 2.2.2. GVCF seçeneğiyle varyant çağırma komutunu çalıştırın

Bir genomik VCF (GVCF) üretmek için, HaplotypeCaller'ın GVCF modunu açan `-ERC GVCF` seçeneğini temel komuta ekliyoruz.

Ayrıca çıktı dosyası için dosya uzantısını `.vcf`'den `.g.vcf`'ye değiştiriyoruz.
Bu teknik olarak bir gereklilik değildir, ancak şiddetle önerilen bir kuraldır.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

??? success "Komut çıktısı"

    ```console hl_lines="39 53 58 59"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.g.vcf -L /data/ref/intervals.bed -ERC GVCF
    00:28:03.593 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:03.765 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.768 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:03.768 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:03.768 INFO  HaplotypeCaller - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:03.768 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:03.769 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:28:03 AM GMT
    00:28:03.769 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.770 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.772 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:28:03.773 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:28:03.773 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:03.774 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:03.774 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:28:03.774 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:28:03.775 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:28:03.775 INFO  HaplotypeCaller - Requester pays: disabled
    00:28:03.776 INFO  HaplotypeCaller - Initializing engine
    00:28:03.896 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:03.919 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:03.934 INFO  HaplotypeCaller - Done initializing engine
    00:28:03.935 INFO  HaplotypeCallerEngine - Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled
    00:28:03.943 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:28:03.945 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:28:03.946 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:28:03.955 INFO  HaplotypeCallerEngine - Standard Emitting and Calling confidence set to -0.0 for reference-model confidence output
    00:28:03.956 INFO  HaplotypeCallerEngine - All sites annotated with PLs forced to true for reference-model confidence output
    00:28:03.972 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:28:03.993 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:28:03.994 INFO  IntelPairHmm - Available threads: 10
    00:28:03.994 INFO  IntelPairHmm - Requested threads: 4
    00:28:03.994 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:28:04.044 INFO  ProgressMeter - Starting traversal
    00:28:04.070 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:28:04.874 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:28:06.535 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:28:06.537 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35            851.6
    00:28:06.538 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:28:06.543 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003648749
    00:28:06.544 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.031498916
    00:28:06.544 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:28:06.547 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:28:06 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.05 minutes.
    Runtime.totalMemory()=281018368
    ```

Bu, konteyner içindeki mevcut çalışma dizininde GVCF çıktı dosyası `reads_mother.g.vcf`'yi oluşturur.

İçeriği görüntülemek için `cat` yaparsanız, bölüm 1'de oluşturduğumuz eşdeğer VCF'den çok daha uzun olduğunu göreceksiniz. Dosyanın başlangıcına bile kaydıramazsınız ve satırların çoğu VCF'de gördüklerimizden oldukça farklı görünür.

??? abstract "Dosya içeriği"

    ```console title="reads_mother.g.vcf" linenums="1674"
    20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
    20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
    20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
    ```

Bunlar, varyant çağırıcının varyasyon kanıtı bulamadığı varyant olmayan bölgeleri temsil eder, bu nedenle varyasyonun yokluğundaki güven düzeyini tanımlayan bazı istatistikleri yakaladı.
Bu, iki çok farklı durum rakamı arasında ayrım yapmayı mümkün kılar: (1) örneğin homozigot-referans olduğunu gösteren kaliteli veri vardır ve (2) her iki şekilde de bir belirleme yapmak için yeterli kaliteli veri yoktur.

Bir GVCF'de, tipik olarak bunların arasına serpiştirilmiş daha az sayıda varyant kaydıyla birlikte bu tür varyant olmayan satırlar çoktur.
Gerçek bir varyant çağrısını bulmak için dosyanın yalnızca ilk 176 satırını yüklemek için GVCF üzerinde `head -176` çalıştırmayı deneyin.

??? abstract "Dosya içeriği"

    ```console title="reads_mother.g.vcf" linenums="174"
    20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
    20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
    20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
    ```

İkinci satır, dosyadaki ilk varyant kaydını gösterir, bu da daha önce baktığımız VCF dosyasındaki ilk varyanta karşılık gelir.

Tıpkı orijinal VCF gibi, çıktı GVCF dosyasına da `reads_mother.g.vcf.idx` adlı bir indeks dosyası eşlik eder.

#### 2.2.3. İşlemi diğer iki örnek üzerinde tekrarlayın

Kalan iki örnek için aşağıdaki komutları birer birer çalıştırarak GVCF'ler oluşturun.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Bu tamamlandığında, mevcut dizininizde `.g.vcf` ile biten üç dosyanız (örnek başına bir tane) ve `.g.vcf.idx` ile biten ilgili indeks dosyalarınız olmalıdır.

Ama konteynerdan çıkmayın!
Bir sonraki adımda aynı konteyneri kullanacağız.

### 2.3. Ortak genotipleme çalıştırın

Artık tüm GVCF'lere sahip olduğumuza göre, bir örnek kohortu için varyant çağrıları oluşturmak için ortak genotipleme yaklaşımını deneyebiliriz.
Tüm GVCF'lerden gelen verileri bir veri deposunda birleştirmeyi ve ardından ortak çağrılmış varyantların nihai VCF'sini oluşturmak için ortak genotipleme analizinin gerçekleştirilmesini içeren iki adımlı bir yöntemdir.

#### 2.3.1. Tüm örnek başına GVCF'leri birleştirin

Bu ilk adım, tüm GVCF'lerden gelen verileri bir GenomicsDB veri deposunda birleştirmek için GenomicsDBImport adlı başka bir GATK aracını kullanır.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

??? success "Komut çıktısı"

    ```console hl_lines="33 36 37 39 40"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenomicsDBImport -V reads_mother.g.vcf -V reads_father.g.vcf -V reads_son.g.vcf -L /data/ref/intervals.bed --genomicsdb-workspace-path family_trio_gdb
    00:28:20.772 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:20.914 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.917 INFO  GenomicsDBImport - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:20.917 INFO  GenomicsDBImport - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:20.917 INFO  GenomicsDBImport - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:20.917 INFO  GenomicsDBImport - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:20.918 INFO  GenomicsDBImport - Start Date/Time: February 8, 2026 at 12:28:20 AM GMT
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.920 INFO  GenomicsDBImport - HTSJDK Version: 4.1.0
    00:28:20.921 INFO  GenomicsDBImport - Picard Version: 3.1.1
    00:28:20.921 INFO  GenomicsDBImport - Built for Spark Version: 3.5.0
    00:28:20.922 INFO  GenomicsDBImport - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:20.923 INFO  GenomicsDBImport - Deflater: IntelDeflater
    00:28:20.924 INFO  GenomicsDBImport - Inflater: IntelInflater
    00:28:20.924 INFO  GenomicsDBImport - GCS max retries/reopens: 20
    00:28:20.924 INFO  GenomicsDBImport - Requester pays: disabled
    00:28:20.925 INFO  GenomicsDBImport - Initializing engine
    00:28:21.144 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:21.152 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:21.157 INFO  GenomicsDBImport - Done initializing engine
    00:28:21.287 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:21.290 INFO  GenomicsDBImport - Vid Map JSON file will be written to /tmp/family_trio_gdb/vidmap.json
    00:28:21.290 INFO  GenomicsDBImport - Callset Map JSON file will be written to /tmp/family_trio_gdb/callset.json
    00:28:21.291 INFO  GenomicsDBImport - Complete VCF Header will be written to /tmp/family_trio_gdb/vcfheader.vcf
    00:28:21.291 INFO  GenomicsDBImport - Importing to workspace - /tmp/family_trio_gdb
    00:28:21.453 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.757 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.859 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.979 INFO  GenomicsDBImport - Done importing batch 1/1
    00:28:21.988 INFO  GenomicsDBImport - Import completed!
    00:28:21.988 INFO  GenomicsDBImport - Shutting down engine
    [February 8, 2026 at 12:28:21 AM GMT] org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport done. Elapsed time: 0.02 minutes.
    Runtime.totalMemory()=305135616
    ```

Bu adımın çıktısı, birleştirilmiş varyant verilerini birden fazla farklı dosya biçiminde tutan daha fazla iç içe dizin içeren bir dizin kümesini içeren etkin bir dizindir.
Etrafında dolaşabilirsiniz ancak bu veri deposu formatının insanlar tarafından doğrudan okunması amaçlanmadığını hızlıca göreceksiniz.

!!! note "Not"

    GATK, gerektiğinde veri deposundan varyant çağrı verilerini incelemeyi ve çıkarmayı mümkün kılan araçlar içerir.

#### 2.3.2. Ortak genotipleme analizinin gerçekleştirilmesi

Bu ikinci adım, kohorttaki tüm örneklerde mevcut veriler ışığında varyant istatistiklerini ve bireysel genotipleri yeniden hesaplamak için GenotypeGVCFs adlı başka bir GATK aracını kullanır.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

??? success "Komut çıktısı"

    ```console hl_lines="30 35 37 38"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenotypeGVCFs -R /data/ref/ref.fasta -V gendb://family_trio_gdb -O family_trio.vcf
    00:28:24.625 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:24.798 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.801 INFO  GenotypeGVCFs - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:24.801 INFO  GenotypeGVCFs - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:24.801 INFO  GenotypeGVCFs - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:24.801 INFO  GenotypeGVCFs - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:24.802 INFO  GenotypeGVCFs - Start Date/Time: February 8, 2026 at 12:28:24 AM GMT
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.804 INFO  GenotypeGVCFs - HTSJDK Version: 4.1.0
    00:28:24.804 INFO  GenotypeGVCFs - Picard Version: 3.1.1
    00:28:24.804 INFO  GenotypeGVCFs - Built for Spark Version: 3.5.0
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:24.806 INFO  GenotypeGVCFs - Deflater: IntelDeflater
    00:28:24.806 INFO  GenotypeGVCFs - Inflater: IntelInflater
    00:28:24.807 INFO  GenotypeGVCFs - GCS max retries/reopens: 20
    00:28:24.807 INFO  GenotypeGVCFs - Requester pays: disabled
    00:28:24.808 INFO  GenotypeGVCFs - Initializing engine
    00:28:25.023 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:25.081 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field InbreedingCoeff  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAC  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAF  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.109 INFO  GenotypeGVCFs - Done initializing engine
    00:28:25.184 INFO  ProgressMeter - Starting traversal
    00:28:25.187 INFO  ProgressMeter -        Current Locus  Elapsed Minutes    Variants Processed  Variants/Minute
    00:28:25.446 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    GENOMICSDB_TIMER,GenomicsDB iterator next() timer,Wall-clock time(s),0.15034835899999904,Cpu time(s),0.1355218420000006
    00:28:26.189 INFO  ProgressMeter - 20_10037292_10066351:13953              0.0                  3390         202994.0
    00:28:26.190 INFO  ProgressMeter - Traversal complete. Processed 3390 total variants in 0.0 minutes.
    00:28:26.194 INFO  GenotypeGVCFs - Shutting down engine
    [February 8, 2026 at 12:28:26 AM GMT] org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs done. Elapsed time: 0.03 minutes.
    Runtime.totalMemory()=296747008
    ```

Bu, konteyner içindeki mevcut çalışma dizininde VCF çıktı dosyası `family_trio.vcf`'yi oluşturur.
Makul derecede küçük başka bir dosyadır, bu nedenle içeriğini görüntülemek için bu dosyayı `cat` yapabilir ve ilk birkaç varyant satırını bulmak için yukarı kaydırabilirsiniz.

??? abstract "Dosya içeriği"

    ```console title="family_trio.vcf" linenums="40"
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
    20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
    20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
    20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
    ```

Bu, daha önce oluşturduğumuz VCF'ye benzer görünüyor, bu sefer üç örneğin tümü için genotip düzeyinde bilgiye sahibiz.
Dosyadaki son üç sütun, alfabetik sırayla listelenen örnekler için genotip bloklarıdır.

Test aile üçlümüz için çok ilk varyant için çağrılan genotiplere bakarsak, babanın heterozigot-varyant (`0/1`), anne ve oğlun ise her ikisinin de homozigot-varyant (`1/1`) olduğunu görüyoruz.

Bu nihayetinde veri setinden çıkarmak istediğimiz bilgidir!

#### 2.3.3. GATK konteynırından çıkın

Konteynerdan çıkmak için `exit` yazın.

```bash
exit
```

İsteğiniz normale dönmüş olmalıdır.
Bu, varyant çağırma komutlarının manuel testini sonlandırır.

---

### Özet

Samtools indeksleme ve GATK varyant çağırma komutlarını ilgili konteynerlarında nasıl test edeceğinizi biliyorsunuz; buna GVCF'ler oluşturma ve birden fazla örnek üzerinde ortak genotipleme çalıştırma dahildir.

### Sırada ne var?

Aynı komutları, işi yürütmek için konteynerlar kullanan iş akışlarına nasıl saracağınızı öğrenin.
