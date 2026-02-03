# Bölüm 1: Örnek başına varyant çağırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu kursun ilk bölümünde, bireysel sekanslama örneklerine GATK varyant çağırma uygulayan basit bir varyant çağırma pipeline'ı oluşturmayı göstereceğiz.

### Yöntem genel bakışı

Varyant çağırma, bir referans genoma göre bir genom dizisindeki varyasyonları tanımlamayı amaçlayan bir genomik analiz yöntemidir.
Burada kısa varyantları, _yani_ SNP'leri ve indel'leri çağırmak için tasarlanmış araçları ve yöntemleri kullanacağız.

![GATK pipeline](img/gatk-pipeline.png)

Tam bir varyant çağırma pipeline'ı tipik olarak referansa haritalama (bazen genom hizalaması olarak da anılır) ve varyant filtreleme ve önceliklendirme dahil olmak üzere birçok adım içerir.
Kolaylık sağlamak için, bu kurs bölümünde sadece varyant çağırma kısmına odaklanacağız.

### Veri seti

Aşağıdaki veriyi ve ilgili kaynakları sağlıyoruz:

- İnsan kromozom 20'sinin (hg19/b37'den) küçük bir bölgesinden oluşan **bir referans genom** ve yardımcı dosyaları (dizin ve sekans sözlüğü).
- Bir aile üçlüsüne (anne, baba ve oğul) karşılık gelen **üç tam genom sekanslama örneği**, dosya boyutlarını küçük tutmak için kromozom 20 üzerindeki küçük bir veri dilimine alt kümelenmiştir.
  Bu, referans genoma zaten haritalanmış Illumina kısa okuma sekanslama verisidir ve [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) formatında (Binary Alignment Map, SAM'in sıkıştırılmış versiyonu, Sequence Alignment Map) sağlanmıştır.
- **Genomik aralıkların bir listesi**, yani örneklerimizin varyant çağırmaya uygun veriye sahip olduğu genomdaki koordinatlar, BED formatında sağlanmıştır.

### İş akışı

Bu kurs bölümünde, aşağıdakileri yapan bir iş akışı geliştireceğiz:

1. [Samtools](https://www.htslib.org/) kullanarak her BAM girdi dosyası için bir dizin dosyası oluşturun
2. Her BAM girdi dosyası üzerinde GATK HaplotypeCaller'ı çalıştırarak örnek başına varyant çağrılarını VCF (Variant Call Format) formatında oluşturun

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

!!! note

    Dizin dosyaları biyoinformatik dosya formatlarının yaygın bir özelliğidir; ana dosyanın yapısı hakkında bilgi içerirler ve GATK gibi araçların tüm dosyayı okumak zorunda kalmadan verinin bir alt kümesine erişmesine olanak tanır.
    Bu, bu dosyaların ne kadar büyük olabileceği düşünüldüğünde önemlidir.

---

## 0. Isınma: Samtools ve GATK komutlarını etkileşimli olarak test edin

Öncelikle komutları bir iş akışına sarmaya çalışmadan önce manuel olarak denemek istiyoruz.
İhtiyacımız olan araçlar (Samtools ve GATK) GitHub Codespaces ortamında yüklü değil, bu yüzden bunları konteynerler aracılığıyla kullanacağız (bkz. [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note

     `nf4-science/genomics` dizininde olduğunuzdan emin olun, böylece `pwd` yazdığınızda gösterilen yolun son kısmı `genomics` olsun.

### 0.1. Samtools ile bir BAM girdi dosyasını dizinleyin

Bir Samtools konteyneri çekip etkileşimli olarak başlatacağız ve BAM dosyalarından biri üzerinde `samtools index` komutunu çalıştıracağız.

#### 0.1.1. Samtools konteynerini çekin

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

#### 0.1.2. Samtools konteynerini etkileşimli olarak başlatın

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

#### 0.1.3. Dizinleme komutunu çalıştırın

[Samtools dokümantasyonu](https://www.htslib.org/doc/samtools-index.html) bize bir BAM dosyasını dizinlemek için çalıştırılacak komut satırını verir.

Sadece girdi dosyasını sağlamamız gerekir; araç, girdi dosya adına `.bai` ekleyerek çıktı için otomatik olarak bir isim oluşturacaktır.

```bash
samtools index /data/bam/reads_mother.bam
```

Bu hemen tamamlanmalıdır ve şimdi orijinal BAM girdi dosyasıyla aynı dizinde `reads_mother.bam.bai` adlı bir dosya görmelisiniz.

??? abstract "Dizin içeriği"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

#### 0.1.4. Samtools konteynerinden çıkın

```bash
exit
```

### 0.2. GATK HaplotypeCaller ile varyantları çağırın

Bir GATK konteyneri çekip etkileşimli olarak başlatacağız ve az önce dizinlediğimiz BAM dosyası üzerinde `gatk HaplotypeCaller` komutunu çalıştıracağız.

#### 0.2.1. GATK konteynerini çekin

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

#### 0.2.2. GATK konteynerini etkileşimli olarak başlatın

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

#### 0.2.3. Varyant çağırma komutunu çalıştırın

[GATK dokümantasyonu](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) bize bir BAM dosyası üzerinde varyant çağırma gerçekleştirmek için çalıştırılacak komut satırını verir.

BAM girdi dosyasını (`-I`) ve referans genomu (`-R`), çıktı dosyası için bir isim (`-O`) ve analiz edilecek genomik aralıkların bir listesini (`-L`) sağlamamız gerekir.

Ancak, dizin dosyasının yolunu belirtmemize gerek yok; araç, yerleşik adlandırma ve birlikte konumlandırma kuralına göre otomatik olarak aynı dizinde arayacaktır.
Aynısı referans genomun yardımcı dosyaları (dizin ve sekans sözlüğü dosyaları, `*.fai` ve `*.dict`) için de geçerlidir.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

Çıktı dosyası `reads_mother.vcf`, konteynerdeki çalışma dizininizin içinde oluşturulur, bu nedenle çıktı dosya yolunu değiştirmediğiniz sürece VS Code dosya gezgininde görmezsiniz.
Ancak küçük bir test dosyasıdır, bu nedenle içeriği açmak ve görüntülemek için `cat` komutunu kullanabilirsiniz.
Dosyanın başlangıcına kadar kaydırırsanız, birçok meta veri satırından oluşan bir başlık ve ardından satır başına bir olmak üzere varyant çağrılarının bir listesini bulacaksınız.

```console title="reads_mother.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

Her satır, örneğin sekanslama verisinde tanımlanan olası bir varyantı açıklar. VCF formatını yorumlama konusunda rehberlik için [bu faydalı makaleye](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/) bakın.

Çıktı VCF dosyasına, GATK tarafından otomatik olarak oluşturulan `reads_mother.vcf.idx` adlı bir dizin dosyası eşlik eder.
BAM dizin dosyasıyla aynı işleve sahiptir, araçların tüm dosyayı yüklemeden veri alt kümelerini aramasına ve almasına olanak tanır.

#### 0.2.4. GATK konteynerinden çıkın

```bash
exit
```

### Çıkarımlar

Samtools dizinleme ve GATK varyant çağırma komutlarını ilgili konteynerlerinde nasıl test edeceğinizi biliyorsunuz.

### Sırada ne var?

Aynı komutları, işi yürütmek için konteynerler kullanan iki adımlı bir iş akışına nasıl saracağınızı öğrenin.

---

## 1. Bir BAM dosyası üzerinde Samtools index çalıştıran tek aşamalı bir iş akışı yazın

Size iş akışının ana bölümlerini özetleyen `genomics-1.nf` adlı bir iş akışı dosyası sağlıyoruz.
İşlevsel değildir; amacı sadece gerçek iş akışını yazmak için kullanacağınız bir iskelet olarak hizmet etmektir.

### 1.1. Dizinleme işlemini tanımlayın

Dizinleme işlemini açıklayan `SAMTOOLS_INDEX` adını vereceğimiz bir işlem yazarak başlayalım.

```groovy title="genomics-1.nf" linenums="9"
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

Bu eğitim serisinin Bölüm 1 ve Bölüm 2'sinde öğrendiğiniz tüm parçaları tanımalısınız.

Bu işlem, `input_bam` girdisi aracılığıyla bir dosya yolu iletmemizi gerektirecek, o halde bunu bir sonraki adımda ayarlayalım.

### 1.2. Bir girdi parametre bildirimi ekleyin

Dosyanın üst kısmında, `Pipeline parameters` bölümü altında, `reads_bam` adlı bir CLI parametresi bildiriyor ve ona varsayılan bir değer veriyoruz.
Bu şekilde, tembel olabiliriz ve pipeline'ı başlatmak için komutu yazdığımızda girdiyi belirtmeyiz (geliştirme amaçları için).

```groovy title="genomics-1.nf" linenums="3"
/*
 * Pipeline parametreleri
 */
params {
    // Birincil girdi
    reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
}
```

Şimdi hazır bir işlemimiz ve ona çalışması için bir girdi verecek bir parametremiz var, o halde bu şeyleri birlikte bağlayalım.

!!! note

    `${projectDir}`, mevcut Nextflow iş akışı betiğinin (`genomics-1.nf`) bulunduğu dizini işaret eden yerleşik bir Nextflow değişkenidir.

    Bu, mutlak yolları sabit kodlamadan iş akışı deposuna dahil edilen dosyalara, veri dizinlerine ve diğer kaynaklara başvurmayı kolaylaştırır.

### 1.3. SAMTOOLS_INDEX'i çalıştırmak için workflow bloğu ekleyin

`workflow` bloğunda, `SAMTOOLS_INDEX` işlemine girdi beslemek için bir **kanal** kurmamız gerekir; ardından o kanalın içeriği üzerinde çalışmak için işlemin kendisini çağırabiliriz.

```groovy title="genomics-1.nf" linenums="24"
workflow {

    main:
    // Girdi kanalı oluştur (CLI parametresi aracılığıyla tek dosya)
    reads_ch = channel.fromPath(params.reads_bam)

    // Girdi BAM dosyası için dizin dosyası oluştur
    SAMTOOLS_INDEX(reads_ch)

    publish:
    bam_index = SAMTOOLS_INDEX.out
}
```

Workflow bloğunun iki bölümü vardır:

- `main:` kanal işlemlerini ve işlem çağrılarını içerir
- `publish:` hangi çıktıların yayınlanması gerektiğini bildirir, bunları adlandırılmış hedeflere atar

[Hello Channels](../../hello_nextflow/02_hello_channels.md) bölümünde kullandığımız `.fromPath` kanal fabrikasını kullandığımızı fark edeceksiniz.
Gerçekten de çok benzer bir şey yapıyoruz.
Fark, Nextflow'a dosya yolunun kendisini bir girdi öğesi olarak kanala yüklemesini söylememizdir, içeriğini okumak yerine.

### 1.4. Sonuçların nereye yayınlanacağını tanımlamak için bir output bloğu ekleyin

Workflow bloğundan sonra, iş akışı çıktılarının nereye yayınlanacağını belirten bir `output` bloğu ekliyoruz.

```groovy title="genomics-1.nf" linenums="37"
output {
    bam_index {
        path '.'
    }
}
```

`publish:` bölümündeki her adlandırılmış hedef (`bam_index` gibi), temel çıktı dizinine göre çıktı yolunu yapılandırabileceğiniz kendi bloğuna sahip olur.

!!! note

    Burada kullandığımız veri dosyaları çok küçük olsa da, genomik verilerde çok büyük olabilirler.
    Varsayılan olarak, Nextflow gereksiz dosya kopyalarından kaçınmak için yayın dizinindeki çıktı dosyalarına sembolik bağlantılar oluşturur.
    `mode` seçeneğini kullanarak bu davranışı değiştirebilirsiniz (örn. `mode 'copy'`) gerçek kopyalar oluşturmak için.
    Sembolik bağlantıların `work` dizininizi temizlediğinizde bozulacağını unutmayın, bu nedenle üretim iş akışları için `mode 'copy'` kullanmak isteyebilirsiniz.

### 1.5. Çıktı dizinini yapılandırın

Temel çıktı dizini `outputDir` yapılandırma seçeneği aracılığıyla ayarlanır. Bunu `nextflow.config` dosyasına ekleyin:

=== "Sonra"

    ```groovy title="nextflow.config" hl_lines="2"
    docker.enabled = true
    outputDir = 'results_genomics'
    ```

=== "Önce"

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

### 1.6. Dizinleme adımının çalıştığını doğrulamak için iş akışını çalıştırın

Hadi iş akışını çalıştıralım! Hatırlatma olarak, girdi parametresini bildirdiğimizde girdi için varsayılan bir değer ayarladığımız için komut satırında bir girdi belirtmemize gerek yok.

```bash
nextflow run genomics-1.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

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
    results_genomics/
    └── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ```

İşte burada!

### Çıkarımlar

Bir genomik aracı tek adımlı bir Nextflow iş akışına nasıl saracağınızı ve bir konteyner kullanarak çalıştırmayı biliyorsunuz.

### Sırada ne var?

İlkinin çıktısını tüketen ikinci bir adım ekleyin.

---

## 2. Dizinlenmiş BAM dosyası üzerinde GATK HaplotypeCaller çalıştırmak için ikinci bir işlem ekleyin

Artık girdi dosyamız için bir dizinimiz olduğuna göre, iş akışının ilginç kısmı olan varyant çağırma adımını kurmaya geçebiliriz.

### 2.1. Varyant çağırma işlemini tanımlayın

Varyant çağırma işlemini açıklayan `GATK_HAPLOTYPECALLER` adını vereceğimiz bir işlem yazalım.

```groovy title="genomics-1.nf" linenums="44"
/*
 * GATK HaplotypeCaller ile varyantları çağır
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

Burada her bir çıktı kanalımızı benzersiz şekilde adlandırmak için yeni bir sözdizimi (`emit:`) kullandığımızı fark edeceksiniz ve bunun nedenleri yakında netleşecektir.

Bu komut oldukça fazla girdi alır, çünkü GATK basit bir dizinleme işine kıyasla analizi gerçekleştirmek için daha fazla bilgiye ihtiyaç duyar.
Ancak girdi bloğunda tanımlanan girdilerden daha fazlasının GATK komutunda listelendiğini fark edeceksiniz. Bunun nedeni nedir?

!!! note

    GATK, BAM dizin dosyasını ve referans genomun yardımcı dosyalarını aramayı bilir çünkü bu dosyaları çevreleyen kurallara aşindir.
    Ancak Nextflow alan bağımsız olacak şekilde tasarlanmıştır ve biyoinformatik dosya formatı gereksinimleri hakkında hiçbir şey bilmez.

Nextflow'a bu dosyaları çalışma zamanında çalışma dizininde aşamalı olarak yerleştirmesi gerektiğini açıkça söylememiz gerekir; aksi takdirde yapmayacaktır ve GATK (haklı olarak) dizin dosyalarının eksik olduğu hakkında bir hata verecektir.

Benzer şekilde, çıktı VCF'nin dizin dosyasını (`"${input_bam}.vcf.idx"` dosyası) açıkça listelememiz gerekir, böylece Nextflow sonraki adımlarda gerekli olması durumunda bu dosyayı takip etmeyi bilir.

### 2.2. Yardımcı girdiler için tanımlar ekleyin

Yeni işlemimiz sağlanacak bir avuç ek dosya beklediğinden, `Pipeline parameters` bölümü altında bunlar için bazı CLI parametreleri ve bazı varsayılan değerler (daha öncekiyle aynı nedenlerle) ayarlıyoruz.

```groovy title="genomics-1.nf" linenums="8"
    // Yardımcı dosyalar
    reference: Path = "${projectDir}/data/ref/ref.fasta"
    reference_index: Path = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict: Path = "${projectDir}/data/ref/ref.dict"
    intervals: Path = "${projectDir}/data/ref/intervals.bed"
```

### 2.3. Yardımcı dosya yollarını tutacak değişkenler oluşturun

Ana veri girdileri kanallar aracılığıyla dinamik olarak aktarılırken, yardımcı dosyaları ele almanın iki yaklaşımı vardır. Önerilen yaklaşım, veri akışını daha net ve tutarlı hale getiren açık kanallar oluşturmaktır. Alternatif olarak, file() fonksiyonu daha basit durumlar için kullanılabilir, özellikle aynı dosyaya birden fazla işlemde başvurmanız gerektiğinde - ancak bunun hala dolaylı olarak kanallar oluşturduğuna dikkat edin. <!-- TODO: Açıklık: bu hala yazılan girdilerle gerekli mi? -->

Bunu workflow bloğuna ekleyin (`reads_ch` oluşturulduktan sonra, `main:` bölümünün içine):

```groovy title="genomics-1.nf" linenums="79"
    // Yardımcı dosyalar için dosya yollarını yükle (referans ve aralıklar)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)
```

Bu, yardımcı dosya yollarının bunlara ihtiyaç duyan herhangi bir işleme girdi olarak sağlanabilmesi için hazır hale getirir.

### 2.4. GATK_HAPLOTYPECALLER'ı çalıştırmak için workflow bloğuna bir çağrı ekleyin

Artık ikinci işlemimizi kurduk ve tüm girdiler ve yardımcı dosyalar hazır ve mevcut, workflow gövdesine `GATK_HAPLOTYPECALLER` işlemine bir çağrı ekleyebiliriz.

```groovy title="genomics-1.nf" linenums="88"
    // Dizinlenmiş BAM dosyasından varyantları çağır
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )
```

Bu eğitim serisinin Bölüm 1'inden `*.out` sözdizimini tanımalısınız; Nextflow'a `SAMTOOLS_INDEX` tarafından çıktılanan kanalı almasını ve bunu `GATK_HAPLOTYPECALLER` işlem çağrısına bağlamasını söylüyoruz.

!!! note

    Girdilerin işlem çağrısında işlemin girdi bloğunda listelendikleriyle tamamen aynı sırada sağlandığını fark edeceksiniz.
    Nextflow'da girdiler konumsaldır, yani aynı sıraya _uymalısınız_; ve elbette aynı sayıda öğe olmalıdır.

### 2.5. publish bölümünü ve output bloğunu güncelleyin

VCF çıktılarını içerecek şekilde `publish:` bölümünü güncellememiz ve `output` bloğuna karşılık gelen hedefler eklememiz gerekiyor.

```groovy title="genomics-1.nf" linenums="99"
    publish:
    bam_index = SAMTOOLS_INDEX.out
    vcf = GATK_HAPLOTYPECALLER.out.vcf
    vcf_idx = GATK_HAPLOTYPECALLER.out.idx
}

output {
    bam_index {
        path '.'
    }
    vcf {
        path '.'
    }
    vcf_idx {
        path '.'
    }
}
```

### 2.6. Varyant çağırma adımının çalıştığını doğrulamak için iş akışını çalıştırın

Genişletilmiş iş akışını `-resume` ile çalıştıralım, böylece dizinleme adımını tekrar çalıştırmamız gerekmez.

```bash
nextflow run genomics-1.nf -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Şimdi konsol çıktısına bakarsak, listelenen iki işlemi görüyoruz.

İlk işlem önbellekleme sayesinde beklendiği gibi atlandı, ikinci işlem ise yepyeni olduğu için çalıştırıldı.

Çıktı dosyalarını sonuçlar dizininde (çalışma dizinine sembolik bağlantılar olarak) bulacaksınız.

??? abstract "Dizin içeriği"

    ```console
    results_genomics/
    ├── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */cf/36f756*/reads_mother.bam.vcf
    └── reads_mother.bam.vcf.idx -> */cf/36f756*/reads_mother.bam.vcf.idx
    ```

VCF dosyasını açarsanız, GATK komutunu doğrudan konteynerde çalıştırarak oluşturduğunuz dosyayla aynı içeriği görmelisiniz.

```console title="reads_mother.bam.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

Bu, çalışmamızdaki her örnek için oluşturmayı önemsediğimiz çıktıdır.

### Çıkarımlar

Gerçek analiz işi yapan ve yardımcı dosyalar gibi genomik dosya formatı özellikleriyle başa çıkabilen çok basit iki adımlı bir iş akışı nasıl yapacağınızı biliyorsunuz.

### Sırada ne var?

İş akışının birden fazla örneği toplu olarak işlemesini sağlayın.

---

## 3. İş akışını bir örnek grubu üzerinde çalışacak şekilde uyarlayın

Tek bir örnek üzerinde işlemeyi otomatikleştirebilen bir iş akışına sahip olmak güzel, ama ya 1000 örneğiniz varsa?
Tüm örneklerinizde döngü yapan bir bash betiği yazmanız mı gerekiyor?

Hayır, şükürler olsun! Sadece kodda küçük bir değişiklik yapın ve Nextflow bunu da sizin için halledecektir.

### 3.1. Girdi parametre bildirimini üç örneği listeleyen bir diziye çevirin

Gelin girdi BAM dosyası bildirimindeki o varsayılan dosya yolunu, `Pipeline parameters` bölümünün altında üç test örneğimiz için dosya yollarını listeleyen bir diziye çevirelim.

=== "Sonra"

    ```groovy title="genomics-1.nf" linenums="7"
    // Birincil girdi (üç örnek dizisi)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

=== "Önce"

    ```groovy title="genomics-1.nf" linenums="7"
        // Birincil girdi
        reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
    ```

!!! note

    Yazılan parametre bildirimleri (`reads_bam: Path` gibi) kullanırken, bir dizi değeri atayamazsınız.
    Diziler için, tür açıklamasını atlayın.

Ve aslında yapmamız gereken tek şey bu, çünkü workflow gövdesinde kullandığımız kanal fabrikası (`.fromPath`) girdi kanalına yüklemek için birden fazla dosya yolunu kabul etmekten tek bir tane yüklemek kadar mutlu.

!!! note

    Normalde, örnek listesini iş akışı dosyanıza sabit kodlamak istemezsiniz, ancak işleri basit tutmak için burada bunu yapıyoruz.
    Bu eğitim serisinin ilerleyen bölümlerinde girdileri işlemek için daha zarif yollar sunacağız.

### 3.2. Üç örnek üzerinde çalıştığını doğrulamak için iş akışını çalıştırın

Artık boru tesisatı üç test örneğinin tamamında çalışacak şekilde ayarlandığına göre iş akışını çalıştırmayı deneyelim.

```bash
nextflow run genomics-1.nf -resume
```

Komik şey: bu _çalışabilir_, VEYA _başarısız olabilir_. Örneğin, işte başarılı bir çalıştırma:

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

İş akışı çalıştırmanız başarılı olduysa, şuna benzer bir hata alana kadar tekrar çalıştırın:

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

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

GATK komut hata çıktısına bakarsanız, şuna benzer bir satır olacaktır:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Bu garip, iş akışının ilk adımında BAM dosyalarını açıkça dizinlediğimizi düşünürsek. Boru tesisatında bir sorun mu olabilir?

#### 3.2.1. İlgili çağrılar için çalışma dizinlerini kontrol edin

Konsol çıktısında listelenen başarısız `GATK_HAPLOTYPECALLER` işlem çağrısı için çalışma dizininin içine bakalım.

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

Ne demek? Nextflow bu işlem çağrısının çalışma dizininde bir dizin dosyası aşamalı olarak yerleştirmiş, ancak yanlış olanı. Bu nasıl olabilir?

#### 3.2.2. Kanal içeriğini incelemek için [view() operatörünü](https://www.nextflow.io/docs/latest/reference/operator.html#view) kullanın

Workflow gövdesinde `GATK_HAPLOTYPER` işlem çağrısından önce bu iki satırı ekleyin:

```groovy title="genomics-1.nf" linenums="84"
    // geçici tanılama
    reads_ch.view()
    SAMTOOLS_INDEX.out.view()
```

Ardından iş akışı komutunu tekrar çalıştırın.

```bash
nextflow run genomics-1.nf
```

Bir kez daha, bu başarılı olabilir veya başarısız olabilir. İşte başarılı bir çalıştırma:

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [fervent_pasteur] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    [a2/dbd8d5] GATK_HAPLOTYPECALLER (3) | 3 of 3 ✔
    ```

Ve işte başarısız bir tane:

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [angry_hamilton] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    [a3/cf3a89] GATK_HAPLOTYPECALLER (3) | 1 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (3)'
    ...
    ```

Tekrar başarısız olması için birkaç kez çalıştırmanız gerekebilir.
Bu hata tutarlı bir şekilde yeniden üretilmeyecektir çünkü bireysel işlem çağrılarının yürütme sürelerindeki bazı değişkenliklere bağlıdır.

Eklediğimiz iki `.view()` çağrısının çıktısı başarısız bir çalıştırma için şöyle görünür:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

İlk üç satır girdi kanalına, ikincisi ise çıktı kanalına karşılık gelir.
Üç örnek için BAM dosyalarının ve dizin dosyalarının aynı sırada listelenmediğini görebilirsiniz!

!!! note

    Birden fazla öğe içeren bir kanal üzerinde bir Nextflow işlemi çağırdığınızda, Nextflow yürütmeyi mümkün olduğunca paralelleştirmeye çalışacak ve çıktıları kullanılabilir oldukları sırayla toplayacaktır.
    Sonuç olarak, karşılık gelen çıktılar orijinal girdilerin verildiğinden farklı bir sırada toplanabilir.

Şu anda yazıldığı şekliyle, iş akışı betiğimiz dizin dosyalarının dizinleme adımından girdilerin verildiği anne/baba/oğul sırasıyla aynı sırada listelenerek çıkacağını varsayıyor.
Ancak bunun böyle olacağının garantisi yoktur, bu yüzden bazen (her zaman olmasa da) yanlış dosyalar ikinci adımda eşleştirilir.

Bunu düzeltmek için, BAM dosyalarının ve dizin dosyalarının kanallar aracılığıyla birlikte seyahat etmesini sağlamamız gerekir.

!!! tip

    İş akışı kodundaki `view()` ifadeleri hiçbir şey yapmaz, bu nedenle bunları içeride bırakmak bir sorun değildir.
    Ancak konsol çıktınızı karıştıracaklar, bu nedenle sorunu gidermekle işiniz bittiğinde bunları kaldırmanızı öneririz.

### 3.3. SAMTOOLS_INDEX işleminin çıktısını, girdi dosyasını ve dizinini birlikte tutan bir demete değiştirin

Bir BAM dosyasının ve dizininin yakından ilişkili kalmasını sağlamanın en basit yolu, onları dizin görevinden çıkan bir demete paketlemektir.

!!! note

    Bir **demet**, bir fonksiyondan birden fazla değer döndürmek için yaygın olarak kullanılan sonlu, sıralı bir öğe listesidir. Demetler, birden fazla girdi veya çıktıyı işlemler arasında aktarırken ilişkilerini ve sıralarını korumak için özellikle kullanışlıdır.

İlk olarak, `SAMTOOLS_INDEX` işleminin çıktısını, çıktı bildiriminde BAM dosyasını içerecek şekilde değiştirelim.

=== "Sonra"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Önce"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        path "${input_bam}.bai"
    ```

Bu şekilde, her dizin dosyası orijinal BAM dosyasıyla sıkı bir şekilde eşleştirilecek ve dizinleme adımının genel çıktısı, dosya çiftlerini içeren tek bir kanal olacaktır.

### 3.4. GATK_HAPLOTYPECALLER işleminin girdisini bir demet olacak şekilde değiştirin

İş akışındaki ilk işlemin çıktısının 'şeklini' değiştirdiğimize göre, ikinci işlemin girdi tanımını eşleşecek şekilde güncellememiz gerekiyor.

Özellikle, daha önce `GATK_HAPLOTYPECALLER` işleminin girdi bloğunda iki ayrı girdi yolu bildirdiğimiz yerde, şimdi `SAMTOOLS_INDEX` tarafından yayılan demetin yapısıyla eşleşen tek bir girdi bildiriyoruz.

=== "Sonra"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Önce"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        path input_bam
        path input_bam_index
    ```

Elbette, `GATK_HAPLOTYPECALLER`ın beklediği girdilerin şeklini değiştirdiğimize göre, işlem çağrısını buna göre workflow gövdesinde güncellememiz gerekiyor.

### 3.5. Workflow bloğundaki GATK_HAPLOTYPECALLER çağrısını güncelleyin

Artık BAM dosyası `SAMTOOLS_INDEX` tarafından kanal çıktısına paketlendiğinden, orijinal `reads_ch`'yi `GATK_HAPLOTYPECALLER` işlemine sağlamamıza gerek yok.

Sonuç olarak, o satırı basitçe silebiliriz.

=== "Sonra"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
    ```

=== "Önce"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
    ```

Dizin uyumsuzluğu sorununu çözmek için gerekli olan tüm yeniden kablolama budur.

### 3.6. Demet için publish bölümünü ve output bloğunu güncelleyin

`SAMTOOLS_INDEX.out` artık hem BAM'i hem de dizinini içeren bir demet olduğundan, her iki dosya da birlikte yayınlanacaktır.
Artık her iki dosyayı da içerdiğini yansıtmak için hedefi `bam_index`'ten `indexed_bam`'e yeniden adlandırıyoruz.

=== "Sonra"

    ```groovy title="genomics-1.nf" hl_lines="2"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
    ```

=== "Önce"

    ```groovy title="genomics-1.nf"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    ```

Ayrıca yeni hedef adını kullanmak için output bloğunu güncellememiz gerekiyor:

=== "Sonra"

    ```groovy title="genomics-1.nf" hl_lines="2"
    output {
        indexed_bam {
            path '.'
        }
    ```

=== "Önce"

    ```groovy title="genomics-1.nf"
    output {
        bam_index {
            path '.'
        }
    ```

### 3.7. Her seferinde üç örnek üzerinde doğru çalıştığını doğrulamak için iş akışını çalıştırın

Elbette, kanıt pudingin içindedir, bu yüzden bunun ileriye dönük güvenilir bir şekilde çalışacağından emin olmak için iş akışını birkaç kez daha çalıştıralım.

```bash
nextflow run genomics-1.nf
```

Bu sefer (ve her seferinde) her şey doğru şekilde çalışmalıdır:

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Sonuçlar dizini artık her örnek için (demetten) hem BAM hem de BAI dosyalarını ve VCF çıktılarını içermektedir:

??? abstract "Sonuçlar dizini içeriği"

    ```console
    results_genomics/
    ├── reads_father.bam -> */60/e2614c*/reads_father.bam
    ├── reads_father.bam.bai -> */60/e2614c*/reads_father.bam.bai
    ├── reads_father.bam.vcf -> */b8/91b3c8*/reads_father.bam.vcf
    ├── reads_father.bam.vcf.idx -> */b8/91b3c8*/reads_father.bam.vcf.idx
    ├── reads_mother.bam -> */3e/fededc*/reads_mother.bam
    ├── reads_mother.bam.bai -> */3e/fededc*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */32/5ca037*/reads_mother.bam.vcf
    ├── reads_mother.bam.vcf.idx -> */32/5ca037*/reads_mother.bam.vcf.idx
    ├── reads_son.bam -> */3c/36d1c2*/reads_son.bam
    ├── reads_son.bam.bai -> */3c/36d1c2*/reads_son.bam.bai
    ├── reads_son.bam.vcf -> */d7/a6b046*/reads_son.bam.vcf
    └── reads_son.bam.vcf.idx -> */d7/a6b046*/reads_son.bam.vcf.idx
    ```

İsterseniz, `SAMTOOLS_INDEX` çıktı kanalının içeriğinin nasıl göründüğüne bakmak için `.view()`'i tekrar kullanabilirsiniz:

```groovy title="genomics-1.nf" linenums="92"
SAMTOOLS_INDEX.out.view()
```

Kanalın beklenen üç demeti içerdiğini göreceksiniz (okunabilirlik için dosya yolları kısaltılmıştır).

```console title="Çıktı"
[*/60/e2614c*/reads_father.bam, */60/e2614c*/reads_father.bam.bai]
[*/3e/fededc*/reads_mother.bam, */3e/fededc*/reads_mother.bam.bai]
[*/3c/36d1c2*/reads_son.bam, */3c/36d1c2*/reads_son.bam.bai]
```

Bu, ileriye dönük çok daha güvenli olacaktır.

### Çıkarımlar

İş akışınızı birden fazla örnek üzerinde (bağımsız olarak) çalıştırmayı biliyorsunuz.

### Sırada ne var?

Örnekleri toplu olarak işlemeyi kolaylaştırın.

---

## 4. İş akışının bir grup girdi dosyası içeren bir metin dosyasını kabul etmesini sağlayın

Bir iş akışına birden fazla veri girdi dosyası sağlamanın çok yaygın bir yolu, dosya yollarını içeren bir metin dosyasıyla yapmaktır.
Satır başına bir dosya yolu listeleyen ve başka hiçbir şey olmayan basit bir metin dosyası kadar basit olabilir veya dosya ek meta veriler içerebilir, bu durumda genellikle örnek listesi olarak adlandırılır.

Burada size basit durumu nasıl yapacağınızı göstereceğiz.

### 4.1. Girdi dosya yollarını listeleyen sağlanan metin dosyasını inceleyin

`data/` dizininde bulabileceğiniz `sample_bams.txt` adlı girdi dosya yollarını listeleyen bir metin dosyası zaten yaptık.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Gördüğünüz gibi, satır başına bir dosya yolu listeledik ve bunlar mutlak yollardır.

!!! note

    Burada kullandığımız dosyalar sadece GitHub Codespaces'inizin yerel dosya sistemindedir, ancak bulut depolamadaki dosyalara da işaret edebiliriz.

### 4.2. Parametre varsayılanını güncelleyin

`reads_bam` girdi parametremiz için varsayılan değeri `sample_bams.txt` dosyasına işaret edecek şekilde değiştirelim.

=== "Sonra"

    ```groovy title="genomics-1.nf" linenums="7"
        // Birincil girdi (girdi dosyalarının dosyası, satır başına bir tane)
        reads_bam: Path = "${projectDir}/data/sample_bams.txt"
    ```

=== "Önce"

    ```groovy title="genomics-1.nf" linenums="7"
    // Birincil girdi (üç örnek dizisi)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

Bu şekilde tembel olmaya devam edebiliriz, ancak dosya listesi artık iş akışı kodunun kendisinde yaşamıyor, bu da doğru yönde büyük bir adımdır.

### 4.3. Bir dosyadan satırları okumak için kanal fabrikasını güncelleyin

Şu anda, girdi kanal fabrikamız kendisine verdiğimiz herhangi bir dosyayı dizinleme işlemine beslemek istediğimiz veri girdileri olarak ele alıyor.
Şimdi ona girdi dosya yollarını listeleyen bir dosya verdiğimize göre, dosyayı ayrıştırmak ve içer
