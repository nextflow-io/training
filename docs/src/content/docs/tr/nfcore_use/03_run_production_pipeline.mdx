# Bölüm 3: Bir Üretim Pipeline'ı Çalıştırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[Bölüm 2](./02_configure_execution.md)'de nf-core/demo için parametrelerin nasıl ayarlanacağını ve yapılandırmanın nasıl özelleştirileceğini öğrendiniz.
Şimdi öğrendiklerinizi gerçek bir üretim pipeline'ına, nf-core/rnaseq'e uygulayacağız.

---

## 1. nf-core/rnaseq'i İndirme ve Çalıştırma

Şimdiye kadar eğitim amacıyla tasarlanmış minimal bir pipeline olan `nf-core/demo`'yu kullandık.
Şimdi gerçek bir üretim pipeline'ı indirip test profiliyle çalıştıracağız.

`nf-core/rnaseq` pipeline'ı, toplu RNA dizileme analizinin temel adımlarını gerçekleştirir: kalite kontrolü, adaptör kırpma, okuma hizalama ve gen düzeyinde niceleme.
Bugüne kadar en yaygın kullanılan nf-core pipeline'ı olma özelliğini taşımaktadır.

### 1.1. Pipeline'ı İndirme

İndirmek için aşağıdaki komutu çalıştırın.

```bash
nextflow pull nf-core/rnaseq
```

??? success "Komut çıktısı"

    ```console
    Checking nf-core/rnaseq ...
     downloaded from https://github.com/nf-core/rnaseq.git - revision: 1f03b53ef7 [master]
    ```

Pipeline artık yerel olarak önbelleğe alınmış ve çalıştırılmaya hazırdır.

### 1.2. Test Profilini Çalıştırma

Test profili ve Docker ile çalıştırın:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir rnaseq-results
```

??? failure "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [suspicious_dijkstra] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    [-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT  -
    Plus 47 more processes waiting for tasks…

    Execution cancelled -- Finishing pending tasks before exit
    -[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (WT_REP2)'

    Caused by:
      Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB


    Command executed:

      fq lint \
          --disable-validator P001 \
          SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > WT_REP2.fq_lint.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/nfcore-use/work/xx/xxxxxxxxxxxxxxxxxxxxxx

    Container:
      quay.io/biocontainers/fq:0.12.0--h9ee0642_0

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting
    ```

Bu hatadaki kritik satır şudur:

```console
Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB
```

Varsayılan Codespaces makinesi 8 GB RAM'e sahiptir; bu, Docker Desktop için de tipik varsayılan değerdir.
Pipeline, `FQ_LINT` süreci için 12 GB talep etmektedir; bu, makinenin sağlayabileceğinden fazladır.

Bu 12 GB değeri, `conf/base.config` dosyasında tanımlanan `process_low` kaynak etiketinden gelmektedir:

```groovy title="conf/base.config"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

Daha büyük bir makine türü kullanmak bir seçenek olabilir; ancak test amacıyla mevcut donanım üzerinde çalışabilmek istiyoruz.
Daha iyi bir yaklaşım, özel bir yapılandırma dosyasıyla varsayılan kaynak değerlerini geçersiz kılmaktır.

### 1.3. Özel Yapılandırmayla Yeniden Çalıştırma

Etiket tabanlı kaynak varsayılanlarını geçersiz kılan özel bir yapılandırma dosyası sağlıyoruz.

??? full-code "laptop.config"

    ```groovy title="laptop.config"
    process {
        withLabel: 'process_low' {
            cpus   = 2
            memory = 6.GB
        }
        withLabel: 'process_medium' {
            cpus   = 4
            memory = 6.GB
        }
        withLabel: 'process_high' {
            cpus   = 6
            memory = 6.GB
        }
        withLabel: 'process_high_memory' {
            memory = 6.GB
        }
    }
    ```

[Bölüm 2](./02_configure_execution.md)'de tek bir süreci ada göre hedeflemek için `withName:` kullanımı tanıtılmıştı.
Burada ise aynı etiketi paylaşan tüm süreçleri aynı anda hedeflemek için `withLabel:` kullanıyoruz.

Bu dosya çalışma dizininizde zaten mevcuttur.
Geçersiz kılmaları uygulamak için `-c` ile birlikte kullanın:

```bash
nextflow run nf-core/rnaseq -profile test,docker -c laptop.config --outdir rnaseq-results
```

??? success "Komut çıktısı (pipeline başlatılıyor)"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [romantic_faraday] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    executor >  local (7)
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FQ_LINT (RAP1_IAA_30M_REP1)   | 3 of 5
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQC (RAP1_IAA_30M_REP1)    | 2 of 5
    ...
    ```

Pipeline artık çalışıyor ve görevlerin tek tek tamamlandığını izleyebilirsiniz.
Bu minimal test veri kümesinde, toplamda 200'den fazla görev yürütülerek 15–20 dakika içinde tamamlanacaktır.

Gerçek RNA-seq deneyleri genellikle onlarca örnek içerir ve saatler ya da günler sürer.
Nextflow; HPC zamanlayıcılarını (SLURM, PBS, LSF) ve bulut platformlarını (AWS, Google Cloud, Azure) destekler. Bu platformlar, işi birçok düğüme dağıtarak toplam çalışma süresini önemli ölçüde kısaltabilir.
Ancak bu ortamların kurulumu ciddi bir karmaşıklık gerektirir.

Seqera Platform (Nextflow'un geliştiricileri tarafından geliştirilen), Nextflow pipeline'larını HPC veya bulut altyapısında (kendi altyapınızda ya da sizin adınıza yönetilen bir altyapıda) başlatmak için web tabanlı bir arayüz sunar. Hesaplama ve veri yönetimi özellikleriyle büyük ölçekli pipeline çalıştırma sürecini kolaylaştırır.

!!! tip "İpucu"

    Akademik araştırmacılar, [Seqera akademik programı](https://seqera.io/academic-program/) aracılığıyla Seqera Platform'a ücretsiz erişebilir.

### Özetle

`nf-core/rnaseq`'i indirdiniz, nf-core kaynak etiketlerinin nasıl çalıştığını gördünüz ve bunları özel bir yapılandırma dosyasıyla nasıl geçersiz kılacağınızı öğrendiniz.
Daha da önemlisi, yerel çalıştırmanın gerçek ölçekli analizler için bir varış noktası değil, bir başlangıç noktası olduğunu gördünüz.

### Sırada ne var?

nf-core pipeline'larını çalıştırmanın temellerini öğrendiniz.
Buradan nereye gideceğinizi öğrenmek için [Sonraki adımlar](next_steps.md) bölümüne bakın.

---

## Özet

Bu bölümde şunları öğrendiniz:

- Üretim ölçeğinde bir pipeline'ı (nf-core/rnaseq) indirip çalıştırmak ve varsayılan kaynak etiketlerini geçersiz kılmak
