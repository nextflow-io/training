# Bölüm 2: Bir kohort üzerinde birleşik çağrı yapma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu kursun ilk bölümünde, tamamen doğrusal olan ve her örneğin verilerini diğerlerinden bağımsız olarak işleyen bir varyant çağırma pipeline'ı oluşturdunuz.
Ancak gerçek bir genomik kullanım durumunda, tipik olarak birden fazla örneğin varyant çağrılarına birlikte bakmanız gerekecektir.

Bu ikinci bölümde, Bölüm 1'deki pipeline'ı temel alarak GATK ile birleşik varyant çağırma işlemini uygulamak için kanalları ve kanal operatörlerini nasıl kullanacağınızı gösteriyoruz.

### Yönteme genel bakış

Bu kursun ilk bölümünde kullandığımız GATK varyant çağırma yöntemi basitçe örnek başına varyant çağrıları üretti.
Her örneğin varyantlarına yalnızca izole olarak bakmak istiyorsanız bu iyidir, ancak bu sınırlı bilgi verir.
Varyant çağrılarının birden fazla örnek arasında nasıl farklılaştığına bakmak genellikle daha ilgi çekicidir ve bunu yapmak için GATK, burada gösterdiğimiz birleşik varyant çağırma adında alternatif bir yöntem sunar.

Birleşik varyant çağırma, her örnek için GVCF (Genomik VCF için) adı verilen özel bir varyant çıktısı türü oluşturmayı, ardından tüm örneklerden GVCF verilerini birleştirmeyi ve son olarak bir 'birleşik genotipleme' istatistiksel analizi çalıştırmayı içerir.

![Birleşik analiz](img/joint-calling.png)

Bir örneğin GVCF'sinin özel yanı, genomun hedeflenen alanındaki tüm pozisyonlar hakkında dizi veri istatistiklerini özetleyen kayıtlar içermesidir, yalnızca programın varyasyon kanıtı bulduğu pozisyonları değil.
Bu, birleşik genotipleme hesaplaması için kritiktir ([daha fazla okuma](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

GVCF, Bölüm 1'de kullandığımız aynı araç olan GATK HaplotypeCaller tarafından ek bir parametre (`-ERC GVCF`) ile üretilir.
GVCF'leri birleştirmek, örnek başına çağrıları bir veri deposuna (bir veritabanına benzer) birleştiren GATK GenomicsDBImport ile yapılır, ardından gerçek 'birleşik genotipleme' analizi GATK GenotypeGVCFs ile yapılır.

### İş akışı

Özetlemek gerekirse, kursun bu bölümünde aşağıdakileri yapan bir iş akışı geliştireceğiz:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Samtools kullanarak her BAM girdi dosyası için bir dizin dosyası oluşturma
2. Her BAM girdi dosyası üzerinde örnek başına genomik varyant çağrılarının bir GVCF'sini oluşturmak için GATK HaplotypeCaller'ı çalıştırma
3. Tüm GVCF'leri toplayıp bunları bir GenomicsDB veri deposunda birleştirme
4. Kohort düzeyinde bir VCF üretmek için birleştirilmiş GVCF veri deposu üzerinde birleşik genotipleme çalıştırma

Bunu Bölüm 1'dekiyle aynı veri setine uygulayacağız.

---

## 0. Isınma: Samtools ve GATK'yi doğrudan çalıştırın

Daha önce olduğu gibi, komutları bir iş akışına sarmaya çalışmadan önce manuel olarak denemek istiyoruz.

!!! note

     Doğru çalışma dizininde olduğunuzdan emin olun:
     `cd /workspaces/training/nf4-science/genomics`

### 0.1. Samtools ile bir BAM girdi dosyasını indeksleyin

Bu ilk adım Bölüm 1'dekiyle aynıdır, dolayısıyla çok tanıdık gelmeli, ancak bu sefer bunu üç örnek için de yapmamız gerekiyor.

!!! note

    Teknik olarak pipeline'ımız aracılığıyla üç örnek için zaten dizin dosyaları oluşturduk, bu yüzden gidip bunları sonuçlar dizininden çıkarabiliriz. Ancak bunu manuel olarak yeniden yapmak daha temizdir ve sadece bir dakika sürecektir.

#### 0.1.1. Samtools konteynerini etkileşimli olarak başlatın

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

#### 0.1.2. Üç örnek için indeksleme komutunu çalıştırın

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

Daha önce olduğu gibi, bu dizin dosyalarını ilgili BAM dosyalarıyla aynı dizinde üretmelidir.

??? abstract "Dizin içeriği"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

Artık üç örneğin tümü için dizin dosyalarımız olduğuna göre, her biri için GVCF'leri oluşturmaya geçebiliriz.

#### 0.1.3. Samtools konteynerinden çıkın

```bash
exit
```

### 0.2. GATK HaplotypeCaller ile GVCF modunda varyantları çağırın

Bu ikinci adım Bölüm 1: Merhaba Genomik'te yaptığımıza çok benzer, ancak şimdi GATK'yi 'GVCF modunda' çalıştıracağız.

#### 0.2.1. GATK konteynerini etkileşimli olarak başlatın

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

#### 0.2.2. Varyant çağırma komutunu GVCF seçeneğiyle çalıştırın

Genomik bir VCF (GVCF) üretmek için temel komuta `-ERC GVCF` seçeneğini ekliyoruz, bu da HaplotypeCaller'ın GVCF modunu açar.

Ayrıca çıktı dosyası için dosya uzantısını `.vcf`'den `.g.vcf`'ye değiştiriyoruz.
Bu teknik olarak bir gereklilik değil, ancak şiddetle tavsiye edilen bir kuraldır.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

Bu, konteynerdeki mevcut çalışma dizininde `reads_mother.g.vcf` GVCF çıktı dosyasını oluşturur.

İçeriğini görüntülemek için `cat` yaparsanız, Bölüm 1'de oluşturduğumuz eşdeğer VCF'den çok daha uzun olduğunu göreceksiniz. Dosyanın başlangıcına bile kaydıramazsınız ve satırların çoğu Bölüm 1'deki VCF'de gördüklerimizden oldukça farklı görünüyor.

```console title="Çıktı" linenums="1674"
20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
```

Bunlar, varyant çağırıcının varyasyon kanıtı bulmadığı varyant olmayan bölgeleri temsil eder, bu nedenle varyasyonun yokluğuna olan güven düzeyini açıklayan bazı istatistikler yakalamıştır. Bu, iki çok farklı durum rakamını ayırt etmeyi mümkün kılar: (1) örneğin homozigot-referans olduğunu gösteren kaliteli veriler vardır ve (2) her iki şekilde de bir belirleme yapmak için yeterli iyi veri mevcut değildir.

Bir GVCF'de, tipik olarak aralarına serpiştirilmiş daha az sayıda varyant kaydı ile birlikte bu tür çok sayıda varyant olmayan satır vardır. Gerçek bir varyant çağrısını bulmak için dosyanın ilk 176 satırını yüklemek için GVCF üzerinde `head -176` çalıştırmayı deneyin.

```console title="Çıktı" linenums="174"
20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
```

İkinci satır dosyadaki ilk varyant kaydını gösterir; bu, Bölüm 1'de baktığımız VCF dosyasındaki ilk varyanta karşılık gelir.

Orijinal VCF'nin olduğu gibi, çıktı GVCF dosyasına da `reads_mother.g.vcf.idx` adlı bir dizin dosyası eşlik eder.

#### 0.2.3. İşlemi diğer iki örnek üzerinde tekrarlayın

Birleşik genotipleme adımını test etmek için üç örneğin tümü için GVCF'lere ihtiyacımız var, bu yüzden şimdi bunları manuel olarak oluşturalım.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

Bu tamamlandığında, mevcut dizininizde `.g.vcf` ile biten üç dosyanız (örnek başına bir) ve `.g.vcf.idx` ile biten ilgili dizin dosyalarınız olmalıdır.

### 0.3. Birleşik genotipleme çalıştırın

Artık tüm GVCF'lere sahip olduğumuza göre, sonunda bir örnek kohortu için varyant çağrıları oluşturmak için birleşik genotipleme yaklaşımını deneyebiliriz.
Hatırlatma olarak, tüm GVCF'lerden verileri bir veri deposunda birleştirmeyi ve ardından birleşik çağrılmış varyantların nihai VCF'sini oluşturmak için asıl birleşik genotipleme analizini çalıştırmayı içeren iki adımlı bir yöntemdir.

#### 0.3.1. Tüm örnek başına GVCF'leri birleştirin

Bu ilk adım, tüm GVCF'lerden verileri bir GenomicsDB veri deposunda birleştirmek için GenomicsDBImport adlı başka bir GATK aracını kullanır.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

Bu adımın çıktısı, etkin bir şekilde, birleştirilmiş varyant verilerini birden fazla farklı dosya biçiminde tutan daha fazla iç içe dizin içeren bir dizindir.
Etrafına göz atabilirsiniz ancak bu veri deposu formatının doğrudan insanlar tarafından okunması amaçlanmadığını hızla göreceksiniz.

!!! note

    GATK, gerektiğinde veri deposundan varyant çağrı verilerini incelemeyi ve çıkarmayı mümkün kılan araçlar içerir.

#### 0.3.2. Asıl birleşik genotipleme analizini çalıştırın

Bu ikinci adım, kohorttaki tüm örnekler arasında mevcut veriler ışığında varyant istatistiklerini ve bireysel genotipleri yeniden hesaplamak için GenotypeGVCFs adlı başka bir GATK aracını kullanır.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

<!--
??? success "Komut çıktısı"

    ```console

    ```
-->

Bu, konteynerdeki mevcut çalışma dizininde `family_trio.vcf` VCF çıktı dosyasını oluşturur.
Oldukça küçük bir başka dosya, bu yüzden içeriğini görüntülemek için bu dosyayı `cat` yapabilir ve ilk birkaç varyant satırını bulmak için yukarı kaydırabilirsiniz.

```console title="family_trio.vcf" linenums="40"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
```

Bu, Bölüm 1'de oluşturduğumuz orijinal VCF'ye daha çok benziyor, ancak bu sefer üç örneğin tümü için genotip düzeyinde bilgiye sahibiz.
Dosyadaki son üç sütun, alfabetik sırada listelenen örnekler için genotip bloklarıdır.

İlk varyant için test aile üçlümüz için çağrılan genotiplere bakarsak, babanın heterozigot-varyant (`0/1`) ve anne ve oğulun her ikisinin de homozigot-varyant (`1/1`) olduğunu görüyoruz.

Sonuçta veri setinden çıkarmak istediğimiz bilgi budur! O halde gidelim tüm bunları bir Nextflow iş akışına saralım ki bunu ölçekte yapabilelim.

#### 0.3.3. GATK konteynerinden çıkın

```bash
exit
```

### Çıkarım

Örnek başına varyant çağırmada yer alan bireysel komutları terminalde çalıştırarak istediğiniz bilgiyi üreteceklerini doğrulamayı biliyorsunuz.

### Sırada ne var?

Bu komutları gerçek bir pipeline'a sarın.

---

## 1. Örnek başına varyant çağırma adımını GVCF üretecek şekilde değiştirin

İyi haber şu ki, Bölüm 1'de bu işin bir kısmını yapan bir iş akışı yazdığımız için baştan başlamamıza gerek yok.
Ancak bu pipeline VCF dosyaları üretirken, şimdi birleşik genotipleme yapmak için GVCF dosyaları istiyoruz.
Bu yüzden GVCF varyant çağırma modunu açarak ve çıktı dosyası uzantısını güncelleyerek başlamamız gerekiyor.

!!! note

    Kolaylık sağlamak için, Bölüm 1'in sonunda olduğu gibi GATK iş akışının yeni bir kopyasıyla çalışacağız, ancak farklı bir ad altında: `genomics-2.nf`.

### 1.1. HaplotypeCaller'a bir GVCF yaymasını söyleyin ve çıktı uzantısını güncelleyin

Kod düzenleyicide `genomics-2.nf` dosyasını açalım.
Çok tanıdık gelmeli, ancak beklendiği gibi çalıştığından emin olmak istiyorsanız çalıştırmaktan çekinmeyin.

İki değişiklik yaparak başlayacağız:

- GATK HaplotypeCaller komutuna `-ERC GVCF` parametresini ekleyin;
- GATK kuralına göre ilgili `.g.vcf` uzantısını kullanmak için çıktı dosyası yolunu güncelleyin.

`-ERC GVCF` eklediğinizde önceki satırın sonuna bir ters eğik çizgi (`\`) eklediğinizden emin olun.

=== "Sonra"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4 6"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "Önce"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

HaplotypeCaller'ı VCF'ler yerine GVCF'ler oluşturacak şekilde değiştirmek için gereken tek şey bu, değil mi?

### 1.2. GVCF'leri oluşturabildiğinizi doğrulamak için pipeline'ı çalıştırın

Nextflow yürütme komutu daha öncekiyle aynıdır, iş akışı dosyası adının kendisi hariç.
Bunu uygun şekilde güncellediğinizden emin olun.

```bash
nextflow run genomics-2.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_venter] DSL2 - revision: a2d6f6f09f

    executor >  local (6)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [72/3249ca] GATK_HAPLOTYPECALLER (3) | 0 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Missing output file(s) `reads_son.bam.vcf` expected by process `GATK_HAPLOTYPECALLER (2)`

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_son.bam         -O reads_son.bam.g.vcf         -L intervals.bed         -ERC GVCF
    ```

Ve çıktı... tamamen kırmızı! Ah hayır.

Yürütülen komut doğru, bu yüzden GATK aracının davranışını değiştirmek için bunun yeterli olduğu konusunda haklıydık.
Ancak eksik çıktı dosyasıyla ilgili o satıra bakın. Bir şey fark ettiniz mi?

Evet, Nextflow'a yeni bir dosya adı beklemesini söylemeyi unuttuk. Hata.

### 1.3. Süreç çıktıları bloğundaki çıktı dosyası uzantısını da güncelleyin

Çünkü araç komutunun kendisinde dosya uzantısını değiştirmek yeterli değildir, Nextflow'a beklenen çıktı dosya adının değiştiğini de söylemeniz gerekir.

=== "Sonra"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Önce"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

### 1.4. Yeni GVCF çıktıları için yayın hedeflerini güncelleyin

Artık VCF'ler yerine GVCF'ler ürettiğimize göre, iş akışının `publish:` bölümünü daha açıklayıcı isimler kullanacak şekilde güncellememiz gerekir.
Ayrıca netlik için GVCF dosyalarını kendi alt dizinlerine düzenleyeceğiz.

=== "Sonra"

    ```groovy title="genomics-2.nf" linenums="88" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Önce"

    ```groovy title="genomics-2.nf" linenums="88"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

### 1.5. Yeni dizin yapısı için çıktı bloğunu güncelleyin

GVCF dosyalarını bir `gvcf` alt dizinine koymak için `output` bloğunu da güncellememiz gerekiyor.

=== "Sonra"

    ```groovy title="genomics-2.nf" linenums="94" hl_lines="3 5 6 8 9"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

=== "Önce"

    ```groovy title="genomics-2.nf" linenums="94"
    output {
        indexed_bam {
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

### 1.6. Pipeline'ı tekrar çalıştırın

Bu sefer `-resume` ile çalıştıralım.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (3)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Bu sefer çalışıyor.

Nextflow çıktısının kendisi (normal VCF modunda başarılı bir çalıştırmaya kıyasla) farklı görünmüyor, ancak şimdi alt dizinlerde düzenlenmiş `.g.vcf` dosyalarını ve ilgili dizin dosyalarını, üç örneğin tümü için bulabiliriz.

??? abstract "Dizin içeriği (sembolik bağlantılar kısaltılmış)"

    ```console
    results_genomics/
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

GVCF dosyalarından birini açıp içinde gezinirseniz, GATK HaplotypeCaller'ın istendiği gibi GVCF dosyaları ürettiğini doğrulayabilirsiniz.

### Çıkarım

Tamam, bu Nextflow öğrenimi açısından minimaldi...
Ancak süreç çıktı bloğunun önemini yinelemek için güzel bir fırsattı!

### Sırada ne var?

Bir kanalın içeriğini toplamayı ve bunları tek bir girdi olarak bir sonraki sürece aktarmayı öğrenin.

---

## 2. Tüm örnekler arasında GVCF verilerini toplayın ve birleştirin

Şimdi tüm örnek başına GVCF'lerden verileri, yapmak istediğimiz birleşik genotipleme analizini destekleyen bir forma birleştirmemiz gerekiyor.

### 2.1. GVCF'leri birleştirecek süreci tanımlayın

Isınma bölümünde daha önce yaptığımızı hatırlayalım, GVCF'leri birleştirmek, sözde GenomicsDB formatında bir veri deposu üretecek olan GATK aracı GenomicsDBImport'un işidir.

Bunun nasıl çalışacağını tanımlamak için, ısınma bölümünde daha önce kullandığımız komuta dayalı olarak yeni bir süreç yazalım.

```groovy title="genomics-2.nf" linenums="66"
/*
 * GVCF'leri GenomicsDB veri deposunda birleştirin
 */
process GATK_GENOMICSDB {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    input:
    path all_gvcfs
    path all_idxs
    path interval_list
    val cohort_name

    output:
    path "${cohort_name}_gdb"

    script:
    """
    gatk GenomicsDBImport \
        -V ${all_gvcfs} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
}
```

Ne düşünüyorsunuz, makul görünüyor mu?

Hadi bağlayalım ve ne olacağını görelim.

### 2.2. Varsayılan değeri olan bir `cohort_name` parametresi ekleyin

Kohort için rastgele bir ad sağlamamız gerekiyor.
Eğitim serisinin ilerleyen bölümlerinde bu tür şeyler için örnek meta verilerinin nasıl kullanılacağını öğreneceksiniz, ancak şimdilik sadece `params` kullanarak bir CLI parametresi bildiriyoruz ve kolaylık için ona varsayılan bir değer veriyoruz.

```groovy title="genomics-2.nf" linenums="16"
    // Nihai çıktı dosyası için temel ad
    cohort_name: String = "family_trio"
```

### 2.3. GATK_HAPLOTYPECALLER çıktılarını örnekler arasında toplayın

`GATK_HAPLOTYPECALLER` sürecinin çıktı kanalını olduğu gibi takmış olsaydık, Nextflow süreci her örnek GVCF'si üzerinde ayrı ayrı çağırırdı.
Ancak, Nextflow'un üçünü de tek bir süreç çağrısına birlikte teslim edecek şekilde tüm üç GVCF'yi (ve dizin dosyalarını) paketlemek istiyoruz.

İyi haber: Bunu `collect()` kanal operatörünü kullanarak yapabiliriz. `workflow` gövdesine, GATK_HAPLOTYPECALLER çağrısından hemen sonra aşağıdaki satırları ekleyelim:

```groovy title="genomics-2.nf" linenums="118"
// Örnekler arasında varyant çağırma çıktılarını toplayın
all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
```

Bu biraz karmaşık görünüyor mu? Hadi bunu parçalara ayıralım ve sade bir dile çevirelim.

1. `.out` özelliğini kullanarak başvurulan `GATK_HAPLOTYPECALLER` sürecinin çıktı kanalını alıyoruz.
2. Kanaldan çıkan her 'öğe' bir çift dosyadır: GVCF ve dizin dosyası, bu sırayla çünkü süreç çıktı bloğunda bu sırada listeleniyorlar. Kolaylıkla, son oturumda bu sürecin çıktılarını adlandırdığımız için (`emit:` kullanarak), bir yandan `.vcf` ekleyerek GVCF'leri, diğer yandan `.out` özelliğinden sonra `.idx` ekleyerek dizin dosyalarını seçebiliriz. Bu çıktıları adlandırmamış olsaydık, sırasıyla `.out[0]` ve `.out[1]` ile başvurmamız gerekirdi.
3. Tüm GVCF dosyalarını `all_gvcfs_ch` adlı yeni bir kanalda tek bir öğe halinde paketlemek için `collect()` kanal operatörünü ekliyoruz ve dizin dosyalarıyla aynısını yaparak `all_idxs_ch` adlı yeni kanalı oluşturuyoruz.

!!! tip

    Burada tam olarak ne olduğunu gözünüzde canlandırmakta zorlanıyorsanız, kanal operatörlerini uygulamadan önce ve sonra kanalların içeriğini incelemek için `view()` operatörünü kullanabileceğinizi unutmayın.

Elde edilen `all_gvcfs_ch` ve `all_idxs_ch` kanalları, yeni yazdığımız `GATK_GENOMICSDB` sürecine takacağımız şeylerdir.

!!! note

    Merak ediyorsanız, GVCF'leri ve dizin dosyalarını ayrı ayrı topluyoruz çünkü GATK GenomicsDBImport komutu yalnızca GVCF dosya yollarını görmek istiyor. Neyse ki, Nextflow tüm dosyaları yürütme için birlikte sahnelediğinden, Bölüm 1'de BAM'ler ve dizinleri için yaptığımız gibi dosyaların sırası konusunda endişelenmemize gerek yok.

### 2.4. GATK_GENOMICSDB'yi çalıştırmak için iş akışı bloğuna bir çağrı ekleyin

Bir sürecimiz var ve girdi kanallarımız var. Sadece süreç çağrısını eklememiz gerekiyor.

```groovy title="genomics-2.nf" linenums="122"
    // GVCF'leri bir GenomicsDB veri deposunda birleştirin
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
```

Tamam, her şey bağlandı.

### 2.5. İş akışını çalıştırın

Bakalım bu çalışacak mı.

```bash
nextflow run genomics-2.nf -resume
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [disturbed_bell] DSL2 - revision: 57942246cc

    executor >  local (1)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [51/d350ea] GATK_GENOMICSDB          | 0 of 1
    ERROR ~ Error executing process > 'GATK_GENOMICSDB'

    Caused by:
      Process `GATK_GENOMICSDB` terminated with an error exit status (1)

    Command executed:

      gatk GenomicsDBImport         -V reads_son.bam.g.vcf reads_father.bam.g.vcf reads_mother.bam.g.vcf         -L intervals.bed         --genomicsdb-workspace-path family_trio_gdb
    ```

`-resume` ile çalıştırdığımız için oldukça hızlı çalışıyor, ancak başarısız oluyor!

Ah. İyi tarafından bakarsak, Nextflow'un `GATK_GENOMICSDB` sürecini aldığını ve özellikle onu yalnızca bir kez çağırdığını görüyoruz.
Bu, `collect()` yaklaşımının bir noktaya kadar çalıştığını gösteriyor.
Ama, ve büyük bir ama, süreç çağrısı başarısız oldu.

Yukarıdaki konsol çıktısını incelediğimizde, yürütülen komutun doğru olmadığını görebiliriz.

Hatayı görebiliyor musunuz?
Bu kısma bakın: `-V reads_father.bam.g.vcf reads_son.bam.g.vcf reads_mother.bam.g.vcf`

`gatk GenomicsDBImport`'a tek bir `-V` argümanı için birden fazla GVCF dosyası verdik, ancak araç her GVCF dosyası için ayrı bir `-V` argümanı bekliyor.

Hatırlatma olarak, konteynerde çalıştırdığımız komut buydu:

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

Bu demektir ki GVCF dosyaları paketimizi bir şekilde düzgün biçimlendirilmiş bir komut dizisine dönüştürmemiz gerekiyor.

### 2.6. Her girdi GVCF'si için ayrı bir `-V` argümanı olan bir komut satırı oluşturun

Nextflow'un Groovy'ye dayalı olmasının işe yaradığı yer burası, çünkü gerekli komut dizisini oluşturmak için oldukça basit bazı dizi manipülasyonları kullanmamıza izin verecek.

Özellikle, bu sözdizimini kullanarak: `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Bir kez daha, bunu bileşenlerine ayıralım.

1. İlk olarak, `all_gvcfs` girdi kanalının içeriğini alıyoruz ve üzerine `.collect()` uyguluyoruz (daha önceki gibi).
2. Bu, paketteki her bir GVCF dosya yolunu **closure**'a, `{ gvcf -> "-V ${gvcf}" }` iletmemize izin verir; burada `gvcf` o GVCF dosya yolunu ifade eder.
   Closure, dosya yolunun başına `-V ` eklemek için kullandığımız bir mini fonksiyondur, `"-V ${gvcf}"` biçiminde.
3. Ardından ayırıcı olarak tek bir boşlukla üç dizeyi birleştirmek için `.join(' ')` kullanırız.

Somut bir örnekle şöyle görünür:

1. Üç dosyamız var:

   `[A.ext, B.ext, C.ext]`

2. Closure, dizeleri oluşturmak için her birini değiştirir:

   `"-V A.ext", "-V B.ext", "-V C.ext"`

3. `.join(' ')` işlemi nihai dizeyi oluşturur:

   `"-V A.ext -V B.ext -V C.ext"`

Bu dizeye sahip olduğumuzda, onu `def` anahtar sözcüğüyle tanımlanan yerel bir değişkene, `gvcfs_line`'a atayabiliriz:

`def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Tamam, dizi manipülasyon şeyimiz var. Nereye koyacağız?

Bunu süreç tanımının içinde bir yere koymak istiyoruz, çünkü GVCF dosya yollarını sürece kanalize ettikten _sonra_ yapmak istiyoruz.
Bunun nedeni, Nextflow'un dosyaların kendilerini yürütme için doğru şekilde sahneleyebilmesi için onları dosya yolları olarak görmesi gerektiğidir.

Ama süreçte _tam olarak nereye_ ekleyebiliriz?

İlginç gerçek: `script:`'ten sonra ve `"""` öncesinde rastgele kod ekleyebilirsiniz!

Harika, o halde dizi manipülasyon satırımızı oraya ekleyelim ve oluşturduğu birleştirilmiş dizeyi kullanmak için `gatk GenomicsDBImport` komutunu güncelleyelim.

=== "Sonra"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="2 5"
        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

=== "Önce"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="4"
        script:
        """
        gatk GenomicsDBImport \
            -V ${all_gvcfs} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

Bu, girdileri `gatk GenomicsDBImport`'a doğru şekilde sağlamak için gereken her şey olmalı.

!!! tip

    `gatk GenomicsDBImport` komutunu güncellediğinizde, `${gvcfs_line}` değişkenini değiştirirken `-V ` önekini kaldırdığınızdan emin olun.

### 2.7. GenomicsDB çıktısını beklendiği gibi oluşturduğunu doğrulamak için iş akışını çalıştırın

Pekala, bakalım bu sorunu çözdü mü.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [peaceful_gates] DSL2 - revision: ca0bf847ed

    executor >  local (1)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [76/d13861] GATK_GENOMICSDB          | 1 of 1 ✔
    ```

Aha! Şimdi çalışıyor gibi görünüyor.

İlk iki adım başarıyla atlandı ve üçüncü adım bu sefer harika çalıştı.
GenomicsDB veri deposu çalışma dizininde oluşturuldu ancak sonuçlara yayınlanmadı, çünkü bu sadece birleşik genotipleme için kullanacağımız bir ara formattır.

Bu arada, çıktının tek bir dosya yerine bir dizin olmasını işlemek için özel bir şey yapmamız gerekmedi.

### Çıkarım

Artık bir kanaldan çıktıları nasıl toplayacağınızı ve bunları başka bir sürece tek bir girdi olarak nasıl paketleyeceğinizi biliyorsunuz.
Ayrıca girdileri uygun sözdizimi ile belirli bir araca sağlamak için bir komut satırının nasıl oluşturulacağını biliyorsunuz.

### Sırada ne var?

Aynı sürece ikinci bir komut eklemeyi öğrenin.

---

## 3. Birleşik genotipleme adımını aynı sürecin parçası olarak çalıştırın

Artık birleştirilmiş genomik varyant çağrılarına sahip olduğumuza göre, gerçekten önemsediğimiz nihai çıktıyı üretecek birleşik genotipleme aracını çalıştırabiliriz: kohort düzeyindeki varyant çağrılarının VCF'si.

Lojistik nedenlerle, birleşik genotiplemeyi aynı sürecin içine dahil etmeye karar veriyoruz.

### 3.1. Süreci GATK_GENOMICSDB'den GATK_JOINTGENOTYPING olarak yeniden adlandırın

Süreç birden fazla araç çalıştıracağından, adını tek bir araç adından ziyade genel işleme atıfta bulunacak şekilde değiştiriyoruz.

=== "Sonra"

    ```groovy title="genomics-2.nf"
    /*
     * GVCF'leri GenomicsDB veri deposunda birleştirin ve kohort düzeyinde çağrılar üretmek için birleşik genotipleme çalıştırın
     */
    process GATK_JOINTGENOTYPING {
    ```

=== "Önce"

    ```groovy title="genomics-2.nf"
    /*
     * GVCF'leri GenomicsDB veri deposunda birleştirin
     */
    process GATK_GENOMICSDB {
    ```

Süreç adlarınızı mümkün olduğunca açıklayıcı tutmayı unutmayın, meslektaşlarınız için —ve gelecekteki kendiniz için— okunabilirliği en üst düzeye çıkarmak için!

### 3.2. Birleşik genotipleme komutunu GATK_JOINTGENOTYPING sürecine ekleyin

Script bölümünün içinde birinci komuttan sonra ikinci komutu basitçe ekleyin.

=== "Sonra"

    ```groovy title="genomics-2.nf" linenums="89"  hl_lines="6-10"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb

        gatk GenotypeGVCFs \
            -R ${ref_fasta} \
            -V gendb://${cohort_name}_gdb \
            -L ${interval_list} \
            -O ${cohort_name}.joint.vcf
        """
    ```

=== "Önce"

    ```groovy title="genomics-2.nf" linenums="89"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

İki komut, terminalden manuel olarak çalıştırsaydık aynı şekilde seri olarak çalıştırılacaktır.

### 3.3. Referans genom dosyalarını GATK_JOINTGENOTYPING süreç girdi tanımlarına ekleyin

İkinci komut referans genom dosyalarını gerektirdiğinden, bunları süreç girdilerine eklememiz gerekiyor.

=== "Sonra"

    ```groovy title="genomics-2.nf" linenums="78"  hl_lines="6-8"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict
    ```

=== "Önce"

    ```groovy title="genomics-2.nf" linenums="78"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
    ```

Bunları yazmak can sıkıcı görünebilir, ancak unutmayın, bunları yalnızca bir kez yazarsınız ve ardından iş akışını bir milyon kez çalıştırabilirsiniz. Değer mi?

### 3.4. Kohort düzeyindeki varyant çağrılarının VCF'sini yayacak şekilde süreç çıktı tanımını güncelleyin

GenomicsDB veri deposunu kaydetmeyi gerçekten umursamıyoruz, bu sadece lojistik nedenlerle var olan bir ara formattır, bu yüzden istersek onu çıktı bloğundan kaldırabiliriz.

Gerçekten ilgilendiğimiz çıktı, birleşik genotipleme komutu tarafından üretilen VCF'dir.

=== "Sonra"

    ```groovy title="genomics-2.nf" linenums="87" hl_lines="2 3"
        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx
    ```

=== "Önce"

    ```groovy title="genomics-2.nf" linenums="87"
        output:
        path "${cohort_name}_gdb"
    ```

Neredeyse bitti!

### 3.5. Süreç çağrısını GATK_GENOMICSDB'den GATK_JOINTGENOTYPING'e güncelleyin

İş akışı gövdesindeki süreç çağrısını GATK_GENOMICSDB'den GATK_JOINTGENOTYPING'e yeniden adlandırmayı unutmayalım. Ve bunu yaparken, birleşik genotipleme aracına sağlamamız gerektiğinden referans genom dosyalarını da girdi olarak eklemeliyiz.

=== "Sonra"

    ```groovy title="genomics-2.nf" linenums="126"
    // GVCF'leri bir GenomicsDB veri deposunda birleştirin ve birleşik genotipleme uygulayın
    GATK_JOINTGENOTYPING(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file
    )
    ```

=== "Önce"

    ```groovy title="genomics-2.nf" linenums="126"
    // GVCF'leri bir GenomicsDB veri deposunda birleştirin
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
    ```

Artık süreç tamamen bağlandı.

### 3.6. Birleşik VCF'yi yayın bölümüne ekleyin

Yeni süreçten birleşik VCF çıktılarını yayınlamamız gerekiyor.
İş akışının `publish:` bölümüne bu satırları ekleyin:

```groovy title="genomics-2.nf" linenums="145"
    joint_vcf = GATK_JOINTGENOTYPING.out.vcf
    joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
```

### 3.7. Birleşik VCF hedeflerini çıktı bloğuna ekleyin

Son olarak, birleşik VCF dosyaları için çıktı hedefleri ekleyin.
Bu nihai çıktı olduğundan bunları sonuçlar dizininin köküne koyacağız.

```groovy title="genomics-2.nf" linenums="157"
    joint_vcf {
        path '.'
    }
    joint_vcf_idx {
        path '.'
    }
```

Artık her şey tamamen bağlanmış olmalı.

### 3.8. İş akışını çalıştırın

Sonunda, değiştirilmiş iş akışını çalıştırabiliriz...

```bash
nextflow run genomics-2.nf -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Ve çalışıyor!

Nihai çıktı dosyasını, `family_trio.joint.vcf`'i (ve dosya dizinini), sonuçlar dizininde bulacaksınız.

??? abstract "Dizin içeriği (sembolik bağlantılar kısaltılmış)"

    ```console
    results_genomics/
    ├── family_trio.joint.vcf -> */a6/7cc8ed*/family_trio.joint.vcf
    ├── family_trio.joint.vcf.idx -> */a6/7cc8ed*/family_trio.joint.vcf.idx
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Şüpheci tipteyseniz, birleşik VCF dosyasını açmak için tıklayabilir ve iş akışının bu bölümün başında araçları manuel olarak çalıştırarak elde ettiğiniz aynı varyant çağrılarını oluşturduğunu doğrulayabilirsiniz.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS
```
