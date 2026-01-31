# Bölüm 3: Kodu modüllere taşıma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu kursun ilk bölümünde, tamamen doğrusal olan ve her örneğin verisini diğerlerinden bağımsız olarak işleyen bir varyant çağırma pipeline'ı oluşturdunuz.

İkinci bölümde, Bölüm 1'deki pipeline üzerine inşa ederek GATK ile ortak varyant çağırmayı uygulamak için kanalları ve kanal operatörlerini nasıl kullanacağınızı gösterdik.

Bu bölümde, o workflow'daki kodu modüllere nasıl dönüştüreceğinizi göstereceğiz. Bu eğitim bölümünü takip etmek için Bölüm 1 ve Bölüm 2'yi ve modüllerin temellerini kapsayan [Merhaba Modüller](../../../hello_nextflow/hello_modules.md) bölümünü tamamlamış olmanız gerekir.

---

## 0. Isınma

Workflow'umuzu geliştirmeye başladığımızda, her şeyi tek bir kod dosyasına koyduk.
Şimdi kodumuzu **modülerleştirme** zamanı, _yani_ process tanımlarını modüllere çıkarma.

Bölüm 2'deki ile aynı workflow ile başlayacağız, bunu sizin için `genomics-3.nf` dosyasında sağladık.

!!! note "Not"

     Doğru çalışma dizininde olduğunuzdan emin olun:
     `cd /workspaces/training/nf4-science/genomics`

Başlangıç noktasını doğrulamak için workflow'u çalıştırın:

```bash
nextflow run genomics-3.nf -resume
```

```console title="Çıktı"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [serene_borg] DSL2 - revision: 0cbebb67a1

executor >  local (7)
[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (1) | 3 of 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

Proje dizininizin içinde artık bir `work` dizini ve bir `results_genomics` dizini olacak.

### Çıkarım

Workflow'unuzu modülerleştirmeye başlamaya hazırsınız.

### Sırada ne var?

Genomics workflow'unun process'lerini modüllere taşıyın.

---

## 1. Process'leri modüllere taşıma

[Merhaba Modüller](../../../hello_nextflow/hello_modules.md) bölümünde öğrendiğiniz gibi, process tanımını herhangi bir dizindeki kendi dosyasına kopyalayarak bir modül oluşturabilirsiniz ve bu dosyaya istediğiniz adı verebilirsiniz.

Daha sonra (özellikle teste geldiğimizde) netleşecek nedenlerle, bu eğitimde dosyayı `main.nf` olarak adlandırma ve araç kiti ile komuttan sonra adlandırılmış bir dizin yapısına yerleştirme konvansiyonunu takip edeceğiz.

### 1.1. `SAMTOOLS_INDEX` process'i için bir modül oluşturma

`SAMTOOLS_INDEX` process'i durumunda, 'samtools' araç kiti ve 'index' komuttur. Bu yüzden, `modules/samtools/index` dizin yapısını oluşturacağız ve `SAMTOOLS_INDEX` process tanımını bu dizinin içindeki `main.nf` dosyasına koyacağız.

```bash
mkdir -p modules/samtools/index
touch modules/samtools/index/main.nf
```

`main.nf` dosyasını açın ve `SAMTOOLS_INDEX` process tanımını içine kopyalayın.

```groovy title="modules/samtools/index/main.nf" linenums="1"
/*
 * BAM dizin dosyası oluştur
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    tuple path(input_bam), path("${input_bam}.bai")

    script:
    """
    samtools index '$input_bam'
    """
}
```

Ardından, `SAMTOOLS_INDEX` process tanımını `genomics-3.nf` dosyasından kaldırın ve bir sonraki process tanımından önce modül için bir import bildirimi ekleyin, şu şekilde:

=== "Sonra"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1 2"
    // Modülleri dahil et
    include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'

    /*
     * GATK HaplotypeCaller ile varyant çağır
     */
    process GATK_HAPLOTYPECALLER {
    ```

=== "Önce"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1"
    /*
     * GATK HaplotypeCaller ile varyant çağır
     */
    process GATK_HAPLOTYPECALLER {
    ```

Artık workflow'u tekrar çalıştırabilirsiniz ve öncekiyle aynı şekilde çalışması gerekir. `-resume` bayrağını sağlarsanız, hiçbir yeni görevin çalıştırılması bile gerekmez:

```bash
nextflow run genomics-3.nf -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-3.nf` [sleepy_snyder] DSL2 - revision: aa68d06c43

    [0f/71b55e] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [f1/18971b] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    [0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
    ```

### 1.2. `GATK_HAPLOTYPECALLER` ve `GATK_JOINTGENOTYPING` process'leri için modüller oluşturma

Kalan process'ler için aynı adımları tekrarlayın.
Her process için:

1. Dizin yapısını oluşturun (`modules/gatk/haplotypecaller/` ve `modules/gatk/jointgenotyping/`)
2. Process tanımını içeren bir `main.nf` dosyası oluşturun
3. Process tanımını `genomics-3.nf` dosyasından kaldırın
4. Modül için bir import bildirimi ekleyin

İşiniz bittiğinde, şunu çalıştırarak modüller dizin yapınızın doğru olduğunu kontrol edin:

```bash
tree modules/
```

??? abstract "Dizin içeriği"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf

    5 directories, 3 files
    ```

Ana workflow dosyasında, parametreler bölümünden sonra şuna benzer bir şey de olmalıdır:

```
include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'

workflow {
```

### Çıkarım

Genomics workflow'u örnek olarak kullanarak bir workflow'u modülerleştirme pratiği yaptınız.

### Sırada ne var?

Modülerleştirilmiş workflow'u test edin.

---

## 2. Modülerleştirilmiş workflow'u test etme

Her şeyin hala çalıştığını doğrulamak için modülerleştirilmiş workflow'u çalıştırın.

```bash
nextflow run genomics-3.nf -resume
```

```console title="Çıktı"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [astonishing_venter] DSL2 - revision: ca27264c13

[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
```

Her şey hala çalışıyor, pipeline'ın yeniden başlatılabilirliği dahil.
Sonuçlar `results_genomics` dizinine yayınlanmaya devam ediyor.

```console title="Dizin içeriği"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
```

### Çıkarım

Bir workflow'u modülerleştirdiniz ve öncekiyle aynı şekilde çalıştığını doğruladınız.

### Sırada ne var?

Öğrendiklerinizi gözden geçirin ve teste ilerleyin.

---

## 3. Özet

Workflow'u modülerleştirdiniz ve pipeline'ın çalışma şeklinde hiçbir şey değişmedi.
Bu kasıtlıdır: kodu işlevselliğini etkilemeden yeniden yapılandırdınız.

Modüller yalnızca process mantığını içerir, bu da onları temiz ve yeniden kullanılabilir yapar.
Ana betik neyin nereye yayınlanacağını kontrol ederken, modüller hesaplama görevi üzerine odaklanmış kalır.

Kodunuzu sürdürmeyi kolaylaştıracak şeyler için bir temel oluşturdunuz.
Örneğin, artık nf-test framework'ü kullanarak pipeline'ınıza testler ekleyebilirsiniz.
Bu kursun bir sonraki bölümünde bakacağımız şey budur.
