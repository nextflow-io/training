# Bölüm 3: Kohort üzerinde birleşik çağırma

Bölüm 2'de, her bir örneğin verisini bağımsız olarak işleyen örnek başına varyant çağırma boru hattı oluşturdunuz.
Şimdi bunu genişleterek [Bölüm 1](01_method.md)'de ele alınan birleşik varyant çağırmayı uygulayacağız.

## Görev

Kursun bu bölümünde, iş akışını aşağıdakileri yapacak şekilde genişleteceğiz:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Samtools kullanarak her BAM girdi dosyası için bir indeks dosyası oluşturun
2. Her BAM girdi dosyası üzerinde GATK HaplotypeCaller'ı çalıştırarak örnek başına genomik varyant çağrılarının bir GVCF'sini oluşturun
3. Tüm GVCF'leri toplayın ve bunları bir GenomicsDB veri deposunda birleştirin
4. Kohort seviyesinde bir VCF üretmek için birleştirilmiş GVCF veri deposu üzerinde birleşik genotipleme çalıştırın

Bu bölüm doğrudan Bölüm 2 tarafından üretilen iş akışı üzerine inşa edilir.

??? info "Bu bölümden nasıl başlanır"

    Kursun bu bölümü [Bölüm 2: Örnek başına varyant çağırma](./02_per_sample_variant_calling.md) bölümünü tamamladığınızı ve çalışan bir `genomics.nf` boru hattınız olduğunu varsayar.

    Eğer Bölüm 2'yi tamamlamadıysanız veya bu bölüm için yeni başlamak istiyorsanız, Bölüm 2 çözümünü başlangıç noktanız olarak kullanabilirsiniz.
    `nf4-science/genomics/` dizini içinden şu komutları çalıştırın:

    ```bash
    cp solutions/part2/genomics-2.nf genomics.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Bu size eksiksiz bir örnek başına varyant çağırma iş akışı verir.
    Aşağıdaki komutu çalıştırarak başarıyla çalıştığını test edebilirsiniz:

    ```bash
    nextflow run genomics.nf -profile test
    ```

## Ders planı

Bunu iki adıma ayırdık:

1. **Örnek başına varyant çağırma adımını bir GVCF üretecek şekilde değiştirin.**
   Bu, süreç komutlarını ve çıktılarını güncellemeyi kapsar.
2. **Örnek başına GVCF'leri birleştiren ve genotipleyenbir birleşik genotipleme adımı ekleyin.**
   Bu, `collect()` operatörünü, komut satırı oluşturma için Groovy closure'larını ve çoklu komutlu süreçleri tanıtır.

!!! note

     Doğru çalışma dizininde olduğunuzdan emin olun:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. Örnek başına varyant çağırma adımını bir GVCF üretecek şekilde değiştirin

Bölüm 2'deki boru hattı VCF dosyaları üretir, ancak birleşik çağırma GVCF dosyaları gerektirir.
GVCF varyant çağırma modunu açmamız ve çıktı dosya uzantısını güncellememiz gerekiyor.

[Bölüm 1](01_method.md)'den GVCF varyant çağırma komutunu hatırlayın:

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Bölüm 2'de sardığımız temel HaplotypeCaller komutuna kıyasla, farklar `-ERC GVCF` parametresi ve `.g.vcf` çıktı uzantısıdır.

### 1.1. HaplotypeCaller'a bir GVCF yayınlamasını söyleyin ve çıktı uzantısını güncelleyin

İki değişiklik yapmak için `modules/gatk_haplotypecaller.nf` modül dosyasını açın:

- GATK HaplotypeCaller komutuna `-ERC GVCF` parametresini ekleyin;
- GATK konvansiyonuna uygun olarak çıktı dosya yolunu `.g.vcf` uzantısını kullanacak şekilde güncelleyin.

`-ERC GVCF` eklediğinizde önceki satırın sonuna bir ters eğik çizgi (`\`) eklediğinizden emin olun.

=== "Sonra"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5 7"
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

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

Ayrıca yeni dosya uzantısıyla eşleşmesi için output bloğunu güncellememiz gerekiyor.
Komut çıktısını `.vcf`'den `.g.vcf`'ye değiştirdiğimiz için, süreç `output:` bloğu aynı değişikliği yansıtmalıdır.

### 1.2. Süreç çıktıları bloğunda çıktı dosya uzantısını güncelleyin

=== "Sonra"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Önce"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

Ayrıca iş akışının yayınlama ve çıktı yapılandırmasını yeni GVCF çıktılarını yansıtacak şekilde güncellememiz gerekiyor.

### 1.3. Yeni GVCF çıktıları için yayınlama hedeflerini güncelleyin

Artık VCF'ler yerine GVCF'ler ürettiğimiz için, iş akışının `publish:` bölümünü daha açıklayıcı isimler kullanacak şekilde güncellememiz gerekir.
Ayrıca netlik için GVCF dosyalarını kendi alt dizinlerinde düzenleyeceğiz.

=== "Sonra"

    ```groovy title="genomics.nf" linenums="47" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="47"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Şimdi eşleşmesi için output bloğunu güncelleyin.

### 1.4. Yeni dizin yapısı için output bloğunu güncelleyin

GVCF dosyalarını bir `gvcf` alt dizinine koymak için `output` bloğunu da güncellememiz gerekiyor.

=== "Sonra"

    ```groovy title="genomics.nf" linenums="53" hl_lines="3 5 6 8 9"
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

    ```groovy title="genomics.nf" linenums="53"
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

Modül, yayınlama hedefleri ve output bloğunun tümü güncellendiğine göre, değişiklikleri test edebiliriz.

### 1.5. Boru hattını çalıştırın

Değişikliklerin çalıştığını doğrulamak için iş akışını çalıştırın.

```bash
nextflow run genomics.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (6)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Nextflow çıktısı daha öncekiyle aynı görünüyor, ancak `.g.vcf` dosyaları ve indeks dosyaları artık alt dizinlerde düzenlenmiş durumda.

??? abstract "Dizin içeriği (sembolik bağlantılar kısaltılmış)"

    ```console
    results/
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

### Özet

Bir araç komutunun çıktı dosya adını değiştirdiğinizde, süreç `output:` bloğu ve yayınlama/çıktı yapılandırması eşleşecek şekilde güncellenmelidir.

### Sırada ne var?

Bir kanalın içeriğini toplamayı ve bunları bir sonraki sürece tek bir girdi olarak aktarmayı öğrenin.

---

## 2. Birleşik genotipleme adımı ekleyin

Şimdi örnek başına GVCF'leri toplamamız, bunları bir GenomicsDB veri deposunda birleştirmemiz ve kohort seviyesinde bir VCF üretmek için birleşik genotipleme çalıştırmamız gerekiyor.
[Bölüm 1](01_method.md)'de ele alındığı gibi, bu iki araçlı bir işlemdir: GenomicsDBImport GVCF'leri birleştirir, ardından GenotypeGVCFs nihai varyant çağrılarını üretir.
Her iki aracı da `GATK_JOINTGENOTYPING` adlı tek bir süreçte saracağız.

[Bölüm 1](01_method.md)'den iki komutu hatırlayın:

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

İlk komut örnek başına GVCF'leri ve bir intervals dosyasını alır ve bir GenomicsDB veri deposu üretir.
İkincisi o veri deposunu, bir referans genomu alır ve nihai kohort seviyesi VCF'yi üretir.
Konteyner URI'si HaplotypeCaller ile aynıdır: `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

### 2.1. Girdileri ayarlayın

Birleşik genotipleme süreci henüz sahip olmadığımız iki tür girdi gerektirir: keyfi bir kohort adı ve tüm örneklerden toplanan GVCF çıktıları bir arada paketlenmiş.

#### 2.1.1. Bir `cohort_name` parametresi ekleyin

Kohort için keyfi bir ad sağlamamız gerekiyor.
Eğitim serisinin ilerleyen bölümlerinde bu tür şeyler için örnek meta verilerini nasıl kullanacağınızı öğreneceksiniz, ancak şimdilik sadece `params` kullanarak bir CLI parametresi bildiriyoruz ve kolaylık için ona varsayılan bir değer veriyoruz.

=== "Sonra"

    ```groovy title="genomics.nf" linenums="14" hl_lines="3-4"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"

        // Base name for final output file
        cohort_name: String = "family_trio"
    }
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="14"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"
    }
    ```

#### 2.1.2. HaplotypeCaller çıktılarını örnekler arasında toplayın

`GATK_HAPLOTYPECALLER` çıktı kanalını doğrudan yeni sürece bağlarsak, Nextflow süreci her bir örnek GVCF'si üzerinde ayrı ayrı çağırır.
Her üç GVCF'yi (ve indeks dosyalarını) paketlemek istiyoruz, böylece Nextflow hepsini birlikte tek bir süreç çağrısına verir.

Bunu `collect()` kanal operatörünü kullanarak yapabiliriz.
GATK_HAPLOTYPECALLER çağrısından hemen sonra, `workflow` gövdesine aşağıdaki satırları ekleyin:

=== "Sonra"

    ```groovy title="genomics.nf" hl_lines="4-6"
            intervals_file
        )

        // Collect variant calling outputs across samples
        all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

=== "Önce"

    ```groovy title="genomics.nf"
            intervals_file
        )
    ```

Bunu açıklayalım:

1. `.out` özelliğini kullanarak `GATK_HAPLOTYPECALLER`'dan çıktı kanalını alıyoruz.
2. Bölüm 1'de `emit:` kullanarak çıktıları isimlendirdiğimiz için, GVCF'leri `.vcf` ile ve indeks dosyalarını `.idx` ile seçebiliriz. İsimlendirilmiş çıktılar olmadan, `.out[0]` ve `.out[1]` kullanmak zorunda kalırdık.
3. `collect()` operatörü tüm dosyaları tek bir elemanda paketler, dolayısıyla `all_gvcfs_ch` üç GVCF'yi birlikte içerir ve `all_idxs_ch` üç indeks dosyasını birlikte içerir.

GVCF'leri ve indeks dosyalarını ayrı ayrı toplayabiliriz (bunları tuple'larda birlikte tutmanın aksine) çünkü Nextflow tüm girdi dosyalarını yürütme için bir arada hazırlayacak, dolayısıyla indeks dosyaları GVCF'lerin yanında mevcut olacak.

!!! tip

    Kanal operatörleri uygulamadan önce ve sonra kanalların içeriğini incelemek için `view()` operatörünü kullanabilirsiniz.

### 2.2. Birleşik genotipleme sürecini yazın ve iş akışında çağırın

Bölüm 2'de kullandığımız aynı modeli izleyerek, süreç tanımını bir modül dosyasına yazacağız, iş akışına içe aktaracağız ve hazırladığımız girdiler üzerinde çağıracağız.

#### 2.2.1. Her GVCF'ye bir `-V` argümanı vermek için bir string oluşturun

Süreç tanımını doldurmaya başlamadan önce, çözmemiz gereken bir şey var.
GenomicsDBImport komutu her GVCF dosyası için ayrı bir `-V` argümanı bekler, şöyle:

```bash
gatk GenomicsDBImport \
    -V reads_mother.bam.g.vcf \
    -V reads_father.bam.g.vcf \
    -V reads_son.bam.g.vcf \
    ...
```

Eğer `-V ${all_gvcfs_ch}` yazsaydık, Nextflow sadece dosya adlarını birleştirir ve komutun o kısmı şöyle görünürdü:

```groovy
-V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf
```

Ancak string'in şöyle görünmesi gerekiyor:

```groovy
-V reads_mother.bam.g.vcf -V reads_father.bam.g.vcf -V reads_son.bam.g.vcf
```

Önemlisi, bu string'i toplanan kanaldaki hangi dosyalar olursa olsun dinamik olarak oluşturmamız gerekiyor.
Nextflow (Groovy aracılığıyla) bunu yapmanın özet bir yolunu sağlar:

```groovy
def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
```

Bunu açıklayalım:

1. `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }` her dosya yolu üzerinde yinelenir ve önüne `-V ` ekler, `["-V A.g.vcf", "-V B.g.vcf", "-V C.g.vcf"]` üretir.
2. `.join(' ')` bunları boşluklarla birleştirir: `"-V A.g.vcf -V B.g.vcf -V C.g.vcf"`.
3. Sonuç, komut şablonuna interpolasyon yapabileceğimiz yerel bir değişkene (`def` ile tanımlanan) `gvcfs_line` atanır.

Bu satır, sürecin `script:` bloğunun içine, komut şablonundan önce gider.
`script:` ile komut şablonunun açılış `"""` arasına keyfi Groovy kodu yerleştirebilirsiniz.

Ardından sürecin `script:` bloğunda o string'in tamamına `gvcfs_line` olarak atıfta bulunabileceksiniz.

#### 2.2.2. Birleşik genotipleme süreci için modülü doldurun

Şimdi tam süreci yazmaya başlayabiliriz.

`modules/gatk_jointgenotyping.nf` dosyasını açın ve süreç tanımının ana hatlarını inceleyin.

Yukarıda verilen bilgileri kullanarak süreç tanımını doldurun, ardından çalışmanızı aşağıdaki "Sonra" sekmesindeki çözümle kontrol edin.

=== "Önce"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * GVCF'leri GenomicsDB veri deposunda birleştir ve kohort seviyesi çağrıları üretmek için birleşik genotipleme çalıştır
     */
    process GATK_JOINTGENOTYPING {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Sonra"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1" hl_lines="8 11-17 20-21 24-25 29-33"
    #!/usr/bin/env nextflow

    /*
     * GVCF'leri GenomicsDB veri deposunda birleştir ve kohort seviyesi çağrıları üretmek için birleşik genotipleme çalıştır
     */
    process GATK_JOINTGENOTYPING {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx

        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
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
    }
    ```

Burada dikkat çekmeye değer birkaç şey var.

Daha önce olduğu gibi, komutlar doğrudan referans vermese de birkaç girdi listelenmiştir: `all_idxs`, `ref_index` ve `ref_dict`.
Bunları listelemek, Nextflow'un bu dosyaları komutlarda görünen dosyaların yanında çalışma dizininde hazırlamasını sağlar, GATK bunları isimlendirme konvansiyonlarına göre bulmayı bekler.

`gvcfs_line` değişkeni, GenomicsDBImport için `-V` argümanlarını oluşturmak üzere yukarıda açıklanan Groovy closure'ı kullanır.

Bu süreç terminalde yapacağınız gibi iki komutu seri olarak çalıştırır.
GenomicsDBImport örnek başına GVCF'leri bir veri deposunda birleştirir, ardından GenotypeGVCFs o veri deposunu okur ve nihai kohort seviyesi VCF'yi üretir.
GenomicsDB veri deposu (`${cohort_name}_gdb`) yalnızca süreç içinde kullanılan bir ara yapıdır; output bloğunda görünmez.

Bunu tamamladıktan sonra, süreç kullanıma hazırdır.
Bunu iş akışında kullanmak için modülü içe aktarmanız ve bir süreç çağrısı eklemeniz gerekir.

#### 2.2.3. Modülü içe aktarın

`genomics.nf` dosyasına, mevcut import ifadelerinin altına import ifadesini ekleyin:

=== "Sonra"

    ```groovy title="genomics.nf" linenums="21" hl_lines="3"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    include { GATK_JOINTGENOTYPING } from './modules/gatk_jointgenotyping.nf'
    ```

=== "Önce"

    ```groovy title="genomics.nf" linenums="21"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

Süreç artık iş akışı kapsamında kullanılabilir.

#### 2.2.4. Süreç çağrısını ekleyin

`collect()` satırlarından sonra iş akışı gövdesine `GATK_JOINTGENOTYPING` çağrısını ekleyin:

=== "Sonra"

    ```groovy title="genomics.nf" hl_lines="3-12"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

        // Combine GVCFs into a GenomicsDB data store and apply joint genotyping
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

    ```groovy title="genomics.nf"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

Süreç artık tamamen bağlanmış durumda.
Şimdi çıktıların nasıl yayınlandığını yapılandırıyoruz.

### 2.3. Çıktı işlemeyi yapılandırın

Birleşik VCF çıktılarını yayınlamamız gerekiyor.
Birleşik genotipleme sonuçları için yayınlama hedefleri ve output bloğu girdileri ekleyin.

#### 2.3.1. Birleşik VCF için yayınlama hedefleri ekleyin

İş akışının `publish:` bölümüne birleşik VCF ve indeksini ekleyin:

=== "Sonra"

    ```groovy title="genomics.nf" hl_lines="5-6"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
        joint_vcf = GATK_JOINTGENOTYPING.out.vcf
        joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
    ```

=== "Önce"

    ```groovy title="genomics.nf"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Şimdi eşleşmesi için output bloğunu güncelleyin.

#### 2.3.2. Birleşik VCF için output bloğu girdileri ekleyin

Birleşik VCF dosyaları için girdiler ekleyin.
Bu nihai çıktı olduğu için bunları results dizininin kökünde koyacağız.

=== "Sonra"

    ```groovy title="genomics.nf" hl_lines="11-16"
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
        joint_vcf {
            path '.'
        }
        joint_vcf_idx {
            path '.'
        }
    }
    ```

=== "Önce"

    ```groovy title="genomics.nf"
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

Süreç, yayınlama hedefleri ve output bloğunun tümü yerinde olduğuna göre, eksiksiz iş akışını test edebiliriz.

### 2.4. İş akışını çalıştırın

Her şeyin çalıştığını doğrulamak için iş akışını çalıştırın.

```bash
nextflow run genomics.nf -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

İlk iki adım önceki çalıştırmadan önbelleğe alınmış durumda ve yeni `GATK_JOINTGENOTYPING` adımı her üç örnekten toplanan girdiler üzerinde bir kez çalışıyor.
Nihai çıktı dosyası, `family_trio.joint.vcf` (ve indeksi), results dizinindedir.

??? abstract "Dizin içeriği (sembolik bağlantılar kısaltılmış)"

    ```console
    results/
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

Birleşik VCF dosyasını açarsanız, iş akışının beklenen varyant çağrılarını ürettiğini doğrulayabilirsiniz.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

Artık otomatik, tamamen yeniden üretilebilir bir birleşik varyant çağırma iş akışına sahipsiniz!

!!! note

    Size verdiğimiz veri dosyalarının kromozom 20'nin yalnızca küçük bir bölümünü kapsadığını unutmayın.
    Gerçek bir varyant çağrı setinin boyutu milyonlarca varyantla sayılır.
    Bu yüzden eğitim amaçları için yalnızca verilerin küçük alt kümelerini kullanıyoruz!

### Özet

Bir kanaldan çıktıları toplamayı ve bunları başka bir sürece tek bir girdi olarak paketlemeyi biliyorsunuz.
Ayrıca Groovy closure'ları kullanarak komut satırı oluşturmayı ve tek bir süreçte birden fazla komut çalıştırmayı da biliyorsunuz.

### Sırada ne var?

Kendinizi tebrik edin! Nextflow for Genomics kursunu tamamladınız.

Öğrendiklerinizi gözden geçirmek ve sırada ne olduğunu öğrenmek için son [kurs özetine](./next_steps.md) gidin.
