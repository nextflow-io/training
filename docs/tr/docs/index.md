---
title: Ana Sayfa
description: Nextflow topluluk eğitim portalına hoş geldiniz!
hide:
  - toc
  - footer
---

# Nextflow Eğitimi

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Kendi kendine öğrenme kursları__

    ---

    **Nextflow topluluk eğitim portalına hoş geldiniz!**

    Aşağıdaki kursları kendi hızınızda, web tabanlı ortamımızda veya kendi ortamınızda tamamlayın.
    Her kurs uygulamalıdır ve bağımsız olarak tamamlayabileceğiniz hedefe yönelik alıştırmalar içerir.

    [Kurslara göz atın :material-arrow-down:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-account-group-outline:{ .lg .middle } __Eğitim Etkinlikleri__

    ---

    **Kendi kendine öğrenmenin ötesinde bir şey mi arıyorsunuz?**

    Yapılandırılmış eğitim etkinliklerini, kendi eğitimlerinizi düzenlemek için rehberleri ve açık kaynak lisansımızla katkı politikamızı bulun.

    [Eğitim etkinliklerini görüntüleyin :material-arrow-right:](training_events.md){ .md-button .md-button--secondary .mt-1 }

</div>

!!! note "Yapay Zeka Destekli Çeviri"

    Bu çeviri yapay zeka kullanılarak oluşturulmuş ve insan çevirmenler tarafından gözden geçirilmiştir.
    Geri bildirimlerinizi ve iyileştirme önerilerinizi memnuniyetle karşılıyoruz.
    Daha fazla bilgi için [çeviri kılavuzumuza](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) bakın.

## Nextflow eğitim kursları kataloğu

<div class="grid cards" markdown>

-   :material-account:{ .lg .middle } __Kullanıcılar için__

    ---

    ### :material-play-circle:{.nextflow-primary} Pipeline'ları çalıştırın {.mt-1}

    Herhangi bir kod yazmadan mevcut pipeline'ları çalıştırmayı öğrenin.

    ??? courses "**Nextflow Run:** Nextflow ile pipeline çalıştırın"

        Kod anlamayı gerektirmeyen, Nextflow pipeline'larını çalıştırmaya hızlı bir giriş. Pipeline başlatmayı, çıktıları almayı, konteyner kullanmayı ve temel düzeyde yürütmeyi yapılandırmayı kapsar.

        [Eğitimi görüntüleyin :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ??? courses "**Use nf-core:** Topluluk tarafından denetlenmiş pipeline'ları bulun ve çalıştırın"

        nf-core topluluk projesinden pipeline'ları bulmaya, çalıştırmaya ve yapılandırmaya hızlı bir giriş. Minimal bir demo pipeline'dan başlayarak üretim ölçeğinde bir analiz pipeline'ına kadar ilerler.

        [Eğitimi görüntüleyin :material-arrow-right:](nfcore_use/index.md){ .md-button .md-button--secondary }

    ??? courses "**Scale with Seqera:** Pipeline'ları büyük ölçekte başlatın ve izleyin"

        Seqera Platform ile Nextflow pipeline'larını hem web arayüzünden hem de komut satırından başlatmaya ve izlemeye uygulamalı bir giriş.

        [Eğitimi görüntüleyin :material-arrow-right:](seqera_scale/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-tune:{.nextflow-primary} Yürütmeyi yönetin {.mt-1}

    Pipeline yürütmesini etkin biçimde yönetmeyi öğrenin.

    ??? courses "**Configure Execution:** Kaynakları, yeniden denemeleri ve yürütme profillerini yapılandırın"

        Nextflow pipeline yürütmesini yapılandırmaya uygulamalı bir giriş: farklı hesaplama ortamlarına uyum sağlama, kaynak tahsislerini ve yeniden denemeleri kontrol etme ve önceden tanımlanmış yapılandırma profilleri arasında geçiş yapma.

        [Eğitimi görüntüleyin :material-arrow-right:](config_exec/index.md){ .md-button .md-button--secondary }

    !!! info compact "Daha fazla konu geliyor"

        Performans ayarı, HPC/bulut yürütmesi ve daha fazlası bu bölüm için planlanmaktadır.
        Bir sonraki konuyu belirlemek için [kısa ilgi anketimize](https://seqera.typeform.com/to/JCs91e8v) oy verin.

-   :material-code-tags:{ .lg .middle } __Geliştiriciler için__

    ---

    ### :material-wrench:{.nextflow-primary} Pipeline yazın {.mt-1}

    Kendi Nextflow pipeline'larınızı geliştirmeyi öğrenin.

    ??? courses "**Hello Nextflow:** Sıfırdan kendi pipeline'larınızı geliştirin"

        Bu kurs, basit ama tam işlevli pipeline'lar geliştirmeye yetecek düzeyde Nextflow dilinin temel bileşenlerini; ayrıca pipeline tasarımı, geliştirme ve yapılandırma uygulamalarının kilit unsurlarını kapsar.

        [Eğitimi görüntüleyin :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Build with nf-core:** nf-core araçlarını ve kurallarını kullanın"

        [nf-core](https://nf-co.re/) uyumlu pipeline'lar geliştirmeyi öğrenmek isteyen Nextflow geliştiricileri için.
        Kurs; nf-core şablonunu ve geliştirme en iyi uygulamalarını kullanan, mevcut nf-core modüllerinden yararlanan, basit ama tam işlevli pipeline'lar geliştirmeye yetecek düzeyde nf-core pipeline'larının yapısını kapsar.

        [Eğitimi görüntüleyin :material-arrow-right:](nfcore_build/index.md){ .md-button .md-button--secondary }

    ??? catalog "**Side Quests:** Gelişmiş Nextflow konularına dalın"

        Belirli konularda bilgi ve becerilerini genişletmek ya da derinleştirmek isteyen Nextflow geliştiricileri için tasarlanmış bağımsız mini kurslardan oluşan bir koleksiyon.
        Doğrusal bir sırayla sunulmakla birlikte herhangi bir sırayla alınabilir (bağımlılıklar için her mini kursun giriş bölümüne bakın).

        [Side Quests'e göz atın :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Bilim için Nextflow {.mt-1}

    Belirli bilimsel uygulamalar için Nextflow pipeline'ları geliştirmeyi öğrenin.

    ??? courses "**Genomics:** Varyant çağırma pipeline'ı geliştirin"

        Kendi genomik pipeline'larını geliştirmeyi öğrenmek isteyen araştırmacılar için bir kurs. Temel Nextflow geliştirme kalıplarını göstermek amacıyla varyant çağırma kullanım senaryosunu ele alır.

        [Eğitimi görüntüleyin :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**RNAseq:** Toplu RNAseq işleme pipeline'ı geliştirin"

        Kendi RNAseq pipeline'larını geliştirmeyi öğrenmek isteyen araştırmacılar için bir kurs. Temel Nextflow geliştirme kalıplarını göstermek amacıyla toplu RNAseq işleme kullanım senaryosunu ele alır.

        [Eğitimi görüntüleyin :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Bioimaging:** Görüntüleme pipeline'larını çalıştırın ve yapılandırın"

        Biyogörüntüleme pipeline'larını çalıştırmayı ve yapılandırmayı öğrenmek isteyen araştırmacılar için bir kurs. Temel Nextflow kullanım kalıplarını göstermek amacıyla nf-core/molkart'ı ele alır.

        [Eğitimi görüntüleyin :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

</div>

## Kurulum ve Yardım

<div class="grid cards mb-4" markdown>

-   :material-cog-outline:{ .lg .middle } __Eğitim Ortamı__

    ---

    Nextflow eğitimleri için ortamınızı kurma seçenekleri.

    [Eğitim ortamlarını görüntüleyin :material-arrow-right:](envsetup/index.md){ .md-button .md-button--secondary }

-   :material-tag-outline:{ .lg .middle } __Nextflow sürümleri__

    ---

    Nextflow'un sözdizimi sürümlerinin gelişimini anlama ve yönetme.

    [Sürüm gereksinimlerini kontrol edin :material-arrow-right:](info/nxf_versions.md){ .md-button .md-button--secondary }

-   :material-file-code-outline:{ .lg .middle } __Hello pipeline__

    ---

    Hello pipeline'ının ne yaptığına ve nasıl yapılandırıldığına ilişkin özet.

    [Özeti okuyun :material-arrow-right:](info/hello_pipeline.md){ .md-button .md-button--secondary }

-   :material-lifebuoy:{ .lg .middle } __Yardım alma__

    ---

    Nextflow eğitiminde bir sorunla karşılaştığınızda başvurabileceğiniz kaynaklar.

    [Yardım bulun :material-arrow-right:](help.md){ .md-button .md-button--secondary }

</div>
