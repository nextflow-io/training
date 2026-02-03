---
title: Ana Sayfa
description: Nextflow topluluk eğitim portalına hoş geldiniz!
hide:
  - toc
  - footer
---

# Nextflow Eğitimi

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Bireysel çalışma kursları__

    ---

    **Nextflow topluluk eğitim portalına hoş geldiniz!**

    Aşağıda listelenen eğitim kursları, bireysel çalışma kaynağı olarak kullanılabilecek şekilde tasarlanmıştır.
    Bu kursları istediğiniz zaman, Github Codespaces aracılığıyla sunduğumuz web tabanlı ortamda veya kendi ortamınızda kendi başınıza çalışabilirsiniz.

    [Kursları keşfedin :material-arrow-right:](#nextflow-egitim-kurslari-katalogu){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Ek bilgiler__

    ---

    ??? warning "Sürüm uyumluluğu"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **Ocak 2026 itibarıyla, aksi belirtilmedikçe tüm Nextflow eğitim kurslarımız strict syntax etkinleştirilmiş Nextflow sürüm 25.10.2 veya üstünü gerektirir.**

        Sürüm gereksinimleri ve strict syntax hakkında daha fazla bilgi için lütfen [Nextflow dokümantasyonu geçiş rehberine](https://nextflow.io/docs/latest/strict-syntax.html) bakın.

        Önceki söz dizimine karşılık gelen eğitim materyallerinin eski sürümlerine bu web sayfasının menü çubuğundaki sürüm seçici aracılığıyla ulaşabilirsiniz.

    ??? terminal "Ortam seçenekleri"

        Eğitim için ihtiyacınız olan her şeyin önceden yüklendiği, Github Codespaces aracılığıyla erişilebilen web tabanlı bir eğitim ortamı sunuyoruz (ücretsiz bir GitHub hesabı gerektirir).

        [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Bu sizin ihtiyaçlarınızı karşılamıyorsa, lütfen diğer [Ortam seçeneklerine](./envsetup/index.md) bakın.

    ??? learning "Eğitim etkinlikleri"

        Nextflow eğitimini yapılandırılmış bir etkinlik kapsamında almayı tercih ediyorsanız, bunun için birçok fırsat mevcuttur. Aşağıdaki seçeneklere göz atmanızı öneririz:

        - Topluluk ekibi tarafından üç ayda bir düzenlenen **[Eğitim Haftaları]()**
        - **[Seqera Etkinlikleri](https://seqera.io/events/)** Seqera tarafından düzenlenen yüz yüze eğitim etkinliklerini içerir ('Seqera Sessions' ve 'Nextflow Summit' için arama yapın)
        - **[Nextflow Elçileri]()** yerel topluluklarınız için etkinlikler düzenler
        - **[nf-core etkinlikleri](https://nf-co.re/events)** topluluk hackathon'larını içerir

    ??? people "Eğitmenler için bilgiler"

        Kendi eğitimlerinizi yürüten bir eğitmen iseniz, uygun atıf yaptığınız sürece materyallerimizi doğrudan eğitim portalından kullanabilirsiniz. Ayrıntılar için aşağıdaki 'Katkılar ve atıf' bölümüne bakın.

        Ayrıca, eğitim çalışmalarınızı nasıl daha iyi destekleyebileceğimizi öğrenmek isteriz! Lütfen bize [community@seqera.io](mailto:community@seqera.io) adresinden veya topluluk forumundan (bkz. [Yardım](help.md) sayfası) ulaşın.

    ??? licensing "Açık kaynak lisansı ve katkı politikası"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Bu eğitim materyali [Seqera](https://seqera.io) tarafından geliştirilip sürdürülmektedir ve topluluk yararına açık kaynak lisansı ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) altında yayınlanmaktadır. Bu materyali lisans kapsamı dışında kullanmak istiyorsanız (ticari kullanım ve yeniden dağıtım kısıtlamalarına dikkat edin), lütfen talebinizi görüşmek için [community@seqera.io](mailto:community@seqera.io) adresinden bize ulaşın.

        Topluluktan gelen iyileştirmeleri, düzeltmeleri ve hata raporlarını memnuniyetle karşılıyoruz. Her sayfanın sağ üst köşesinde kod deposuna bağlanan bir :material-file-edit-outline: simgesi bulunur; burada sorunları bildirebilir veya pull request aracılığıyla eğitim kaynak materyalinde değişiklik önerebilirsiniz. Daha fazla ayrıntı için depodaki `README.md` dosyasına bakın.

</div>

!!! note "Yapay Zeka Destekli Çeviri"

    Bu çeviri yapay zeka kullanılarak oluşturulmuş ve insan çevirmenler tarafından gözden geçirilmiştir.
    Geri bildirimlerinizi ve iyileştirme önerilerinizi memnuniyetle karşılıyoruz.
    Daha fazla bilgi için [çeviri kılavuzumuza](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) bakın.

## Nextflow eğitim kursları kataloğu

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Giriş seviyesi__

    ---

    ### :material-compass:{.nextflow-primary} Yeni Başlayanlar için Nextflow {.mt-1}

    Nextflow'a tamamen yeni olanlar için tasarlanmış alana bağımlı olmayan kurslar. Her kurs, öğrencilerin becerilerini kademeli olarak geliştirmelerine yardımcı olmak için tasarlanmış bir dizi eğitim modülünden oluşur.

    ??? courses "**Hello Nextflow:** Kendi iş akışlarınızı geliştirmeyi öğrenin"

        Bu kurs, basit ama tam işlevsel iş akışları geliştirmeye yetecek kadar ayrıntıyla Nextflow dilinin temel bileşenlerini, ayrıca iş akışı tasarımı, geliştirme ve yapılandırma uygulamalarının temel öğelerini kapsar.

        [Hello Nextflow eğitimine başlayın :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Mevcut iş akışlarını çalıştırmayı öğrenin"

        Hello Nextflow geliştirici kursuna dayanan, ancak kod üzerine daha az odaklanan, Nextflow iş akışlarını çalıştırma ve yapılandırmaya kısa bir giriş. Çalıştırma, çıktılar, temel kod yapısı ve farklı hesaplama ortamları için yapılandırmayı kapsar.

        [Nextflow Run eğitimine başlayın :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Bilim için Nextflow {.mt-1}

    'Hello Nextflow'da sunulan kavram ve bileşenleri belirli bilimsel kullanım durumlarına uygulamayı öğrenin.

    ??? courses "**Genomik için Nextflow** (varyant çağırma)"

        Kendi genomik iş akışlarını geliştirmek isteyen araştırmacılar için. Bu kurs, basit ama işlevsel bir genomik iş akışının nasıl geliştirileceğini göstermek için varyant çağırma kullanım örneğini kullanır.

        [Genomik için Nextflow eğitimine başlayın :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**RNAseq için Nextflow** (toplu RNAseq)"

        Kendi RNAseq iş akışlarını geliştirmek isteyen araştırmacılar için. Bu kurs, basit ama işlevsel bir RNAseq iş akışının nasıl geliştirileceğini göstermek için toplu RNAseq işleme kullanım örneğini kullanır.

        [RNAseq için Nextflow eğitimine başlayın :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Görüntüleme için Nextflow** (mekansal omikler)"

        Analiz iş akışlarını çalıştırmayı ve özelleştirmeyi öğrenmek isteyen görüntüleme ve mekansal omikler alanındaki araştırmacılar için. Bu kurs, Nextflow iş akışlarının nasıl çalıştırılacağını, yapılandırılacağını ve girdilerinin nasıl yönetileceğini göstermek için biyolojik olarak ilgili bir iş akışı sağlamak üzere nf-core/molkart iş akışını kullanır.

        [Görüntüleme için Nextflow eğitimine başlayın :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __İleri seviye__

    ---

    ### :material-bridge:{.nextflow-primary} Nextflow'dan nf-core'a {.mt-1}

    [nf-core](https://nf-co.re/) topluluk projesinden kod ve en iyi uygulamaları kullanmayı öğrenin.

    Bu kurslar, Nextflow temellerinden nf-core en iyi uygulamalarına geçmenize yardımcı olur.
    nf-core topluluğunun iş akışlarını nasıl ve neden oluşturduğunu ve bu tekniklere nasıl katkıda bulunabileceğinizi ve yeniden kullanabileceğinizi anlayın.

    ??? courses "**Hello nf-core:** nf-core ile başlayın"

        [nf-core](https://nf-co.re/) uyumlu iş akışlarını çalıştırmayı ve geliştirmeyi öğrenmek isteyen geliştiriciler için. Bu kurs, nf-core şablonunu ve geliştirme en iyi uygulamalarını izleyen basit ama tam işlevsel iş akışları geliştirmeyi sağlayacak kadar ayrıntıyla nf-core iş akışlarının yapısını ve mevcut nf-core modüllerinin kullanımını kapsar.

        [Hello nf-core eğitimine başlayın :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} İleri Nextflow Eğitimi {.mt-1}

    Gerçek dünya kullanım durumlarını ele almak için Nextflow iş akışlarını geliştirme ve dağıtma konusunda ileri kavramları ve mekanizmaları öğrenin.

    ??? courses "**Side Quest'ler:** Bağımsız konulara derinlemesine dalışlar"

        Belirli konularda beceri yelpazelerini genişletmek ve/veya derinleştirmek isteyen Nextflow geliştiricileri için bağımsız mini kurslar. Doğrusal olarak sunulurlar ancak herhangi bir sırayla alınabilirler (her mini kurs genel bakışındaki bağımlılıklara bakın).

        [Side Quest'leri keşfedin :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Eğitim Koleksiyonları:** Side Quest'ler aracılığıyla önerilen öğrenme yolları"

        Eğitim Koleksiyonları, belirli bir tema veya kullanım durumu etrafında kapsamlı bir öğrenme deneyimi sağlamak için birden fazla Side Quest'i birleştirir.

        [Eğitim Koleksiyonlarını keşfedin :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "Arşivlenmiş eğitim materyallerini mi arıyorsunuz?"

    Eski eğitim materyalleri (Temel Eğitim, İleri Eğitim ve diğer deneysel kurslar) Nextflow 3.0 strict syntax ile uyumsuz oldukları için eğitim portalından kaldırılmıştır.
    Bu materyallere erişmeniz gerekiyorsa, Ocak 2026 öncesi [git geçmişinde](https://github.com/nextflow-io/training) mevcutturlar.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
