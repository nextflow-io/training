---
title: Ana Sayfa
description: Nextflow topluluk eğitim portalına hoş geldiniz!
hide:
  - toc
  - footer
---

# Nextflow Eğitimi

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Self-servis kurslar__

    ---

    **Nextflow topluluk eğitim portalına hoş geldiniz!**

    Aşağıda listelenen eğitim kursları, self-servis bir kaynak olarak kullanılabilecek şekilde tasarlanmıştır.
    İstediğiniz zaman, Github Codespaces aracılığıyla sağladığımız web tabanlı ortamda veya kendi ortamınızda bu kursları kendi başınıza tamamlayabilirsiniz.

    [Kursları keşfedin :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Ek bilgiler__

    ---

    ??? warning "Sürüm uyumluluğu"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **Ocak 2026 itibarıyla, aksi belirtilmedikçe tüm Nextflow eğitim kurslarımız Nextflow sürüm 25.10.2 veya daha yenisini ve strict syntax'ın etkinleştirilmesini gerektirir.**

        Sürüm gereksinimleri ve strict syntax hakkında daha fazla bilgi için lütfen [Nextflow dokümanları migrasyon kılavuzuna](https://nextflow.io/docs/latest/strict-syntax.html) bakın.

        Eski söz dizimine karşılık gelen eğitim materyalinin eski sürümlerine, bu web sayfasının menü çubuğundaki sürüm seçici aracılığıyla erişilebilir.

    ??? terminal "Ortam seçenekleri"

        Eğitim için ihtiyacınız olan her şeyin önceden yüklenmiş olduğu, Github Codespaces üzerinden erişilebilen (ücretsiz bir GitHub hesabı gerektirir) web tabanlı bir eğitim ortamı sağlıyoruz.

        [![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Bu sizin ihtiyaçlarınıza uygun değilse, lütfen diğer [Ortam seçeneklerine](./envsetup/index.md) bakın.

    ??? learning "Eğitim etkinlikleri"

        Nextflow eğitimini yapılandırılmış bir etkinliğin parçası olarak almayı tercih ediyorsanız, bunu yapmanız için birçok fırsat bulunmaktadır. Aşağıdaki seçeneklere göz atmanızı öneririz:

        - Topluluk ekibi tarafından üç ayda bir düzenlenen **[Training Weeks]()**
        - **[Seqera Events](https://seqera.io/events/)**, Seqera tarafından düzenlenen yüz yüze eğitim etkinliklerini içerir ('Seqera Sessions' ve 'Nextflow Summit' için arama yapın)
        - **[Nextflow Ambassadors]()**, yerel toplulukları için etkinlikler düzenler
        - **[nf-core events](https://nf-co.re/events)**, topluluk hackathon'larını içerir

    ??? people "Eğitmenler için bilgiler"

        Kendi eğitimlerinizi yürüten bir eğitmen iseniz, uygun şekilde atıf yaptığınız sürece materyallerimizi doğrudan eğitim portalından kullanabilirsiniz. Ayrıntılar için aşağıdaki 'Katkılar ve teşekkürler' bölümüne bakın.

        Ayrıca, eğitim çabalarınızı nasıl daha iyi destekleyebileceğimiz konusunda sizden haber almak isteriz! Lütfen [community@seqera.io](mailto:community@seqera.io) adresinden veya topluluk forumundan bizimle iletişime geçin ([Yardım](help.md) sayfasına bakın).

    ??? licensing "Açık kaynak lisansı ve katkı politikası"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Bu eğitim materyali [Seqera](https://seqera.io) tarafından geliştirilip sürdürülmekte ve topluluğun yararına açık kaynak lisansı ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) altında yayınlanmaktadır. Bu materyali lisansın kapsamı dışında bir şekilde kullanmak istiyorsanız (ticari kullanım ve yeniden dağıtım konusundaki sınırlamalara dikkat edin), lütfen talebinizi görüşmek için [community@seqera.io](mailto:community@seqera.io) adresinden bizimle iletişime geçin.

        Topluluktan gelen iyileştirmeleri, düzeltmeleri ve hata raporlarını memnuniyetle karşılıyoruz. Her sayfanın sağ üst köşesinde, eğitim kaynak materyaline ilişkin sorunları bildirebileceğiniz veya pull request aracılığıyla değişiklik önerebileceğiniz kod deposuna bağlanan bir :material-file-edit-outline: simgesi bulunmaktadır. Daha fazla ayrıntı için depodaki `README.md` dosyasına bakın.

</div>

!!! note "Yapay Zeka Destekli Çeviri"

    Bu çeviri yapay zeka kullanılarak oluşturulmuş ve insan çevirmenler tarafından gözden geçirilmiştir.
    Geri bildirimlerinizi ve iyileştirme önerilerinizi memnuniyetle karşılıyoruz.
    Daha fazla bilgi için [çeviri kılavuzumuza](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) bakın.

## Nextflow eğitim kursları kataloğu

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Başlangıç seviyesi__

    ---

    ### :material-compass:{.nextflow-primary} Yeni Başlayanlar için Nextflow {.mt-1}

    Nextflow'a tamamen yeni olanlar için tasarlanmış, alana özgü olmayan kurslar. Her kurs, öğrencilerin becerilerini aşamalı olarak geliştirmelerine yardımcı olmak için tasarlanmış bir dizi eğitim modülünden oluşur.

    ??? courses "**Hello Nextflow:** Kendi pipeline'larınızı geliştirmeyi öğrenin"

        Bu kurs, basit ama tamamen işlevsel pipeline'lar geliştirmeyi sağlayacak kadar ayrıntılı bir şekilde Nextflow dilinin temel bileşenlerini ve ayrıca pipeline tasarımı, geliştirme ve yapılandırma uygulamalarının temel unsurlarını kapsar.

        [Hello Nextflow eğitimine başlayın :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Mevcut pipeline'ları çalıştırmayı öğrenin"

        Nextflow pipeline'larını çalıştırma ve yapılandırmaya yönelik kısa bir giriş, Hello Nextflow geliştirici kursuna dayanır ancak kod üzerinde daha az odaklanır. Yürütme, çıktılar, temel kod yapısı ve farklı hesaplama ortamları için yapılandırmayı kapsar.

        [Nextflow Run eğitimine başlayın :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Bilim için Nextflow {.mt-1}

    'Hello Nextflow'da sunulan kavram ve bileşenleri belirli bilimsel kullanım durumlarına uygulamayı öğrenin.

    ??? courses "**Nextflow for Genomics** (varyant çağırma)"

        Kendi genomik pipeline'larını geliştirmeyi öğrenmek isteyen araştırmacılar için. Kurs, basit ama işlevsel bir genomik pipeline'ın nasıl geliştirileceğini göstermek için bir varyant çağırma kullanım durumu kullanır.

        [Nextflow for Genomics eğitimine başlayın :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for RNAseq** (bulk RNAseq)"

        Kendi RNAseq pipeline'larını geliştirmeyi öğrenmek isteyen araştırmacılar için. Kurs, basit ama işlevsel bir RNAseq pipeline'ının nasıl geliştirileceğini göstermek için bir bulk RNAseq işleme kullanım durumu kullanır.

        [Nextflow for RNAseq eğitimine başlayın :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for Imaging** (uzamsal omikler)"

        Analiz pipeline'larını çalıştırmayı ve özelleştirmeyi öğrenmek isteyen görüntüleme ve uzamsal omikler alanındaki araştırmacılar için. Kurs, Nextflow pipeline'larını nasıl çalıştıracağınızı, yapılandıracağınızı ve girdileri nasıl yöneteceğinizi göstermek için biyolojik olarak ilgili bir pipeline olan nf-core/molkart pipeline'ını kullanır.

        [Nextflow for Imaging eğitimine başlayın :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __İleri seviye__

    ---

    ### :material-bridge:{.nextflow-primary} Nextflow'dan nf-core'a {.mt-1}

    [nf-core](https://nf-co.re/) topluluk projesinden kod ve en iyi uygulamaları kullanmayı öğrenin.

    Bu kurslar, Nextflow temellerinden nf-core en iyi uygulamalarına geçmenize yardımcı olur.
    nf-core topluluğunun pipeline'ları nasıl ve neden oluşturduğunu ve bu tekniklere nasıl katkıda bulunabileceğinizi ve bunları nasıl yeniden kullanabileceğinizi anlayın.

    ??? courses "**Hello nf-core:** nf-core ile başlayın"

        [nf-core](https://nf-co.re/) uyumlu pipeline'ları çalıştırmayı ve geliştirmeyi öğrenmek isteyen geliştiriciler için. Kurs, nf-core şablonunu ve geliştirme en iyi uygulamalarını takip eden basit ama tamamen işlevsel pipeline'lar geliştirmeyi ve mevcut nf-core modüllerini kullanmayı sağlayacak kadar ayrıntılı bir şekilde nf-core pipeline'larının yapısını kapsar.

        [Hello nf-core eğitimine başlayın :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} İleri Düzey Nextflow Eğitimi {.mt-1}

    Gerçek dünya kullanım durumlarını ele almak için Nextflow pipeline'larını geliştirme ve dağıtma konusunda ileri düzey kavramları ve mekanizmaları öğrenin.

    ??? courses "**Side Quests:** Bağımsız konulara derinlemesine dalışlar"

        Beceri aralıklarını genişletmek ve/veya belirli konularda becerilerini derinleştirmek isteyen Nextflow geliştiricileri için tasarlanmış bağımsız mini-kurslar. Doğrusal olarak sunulurlar ancak herhangi bir sırada alınabilirler (her mini-kurs genel bakışındaki bağımlılıklara bakın).

        [Side Quests'e göz atın :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Side Quests aracılığıyla önerilen öğrenme yolları"

        Training Collections, belirli bir tema veya kullanım durumu etrafında kapsamlı bir öğrenme deneyimi sağlamak için birden fazla Side Quest'i birleştirir.

        [Training Collections'a göz atın :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "Arşivlenmiş eğitim materyallerini mi arıyorsunuz?"

    Eski eğitim materyalleri (Fundamentals Training, Advanced Training ve diğer deneysel kurslar), Nextflow 3.0 strict syntax ile uyumsuz oldukları için eğitim portalından kaldırılmıştır.
    Bu materyallere erişmeniz gerekiyorsa, Ocak 2026'dan önceki [git geçmişinde](https://github.com/nextflow-io/training) mevcutturlar.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
