# Nextflow Eğitimi

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Kendi kendine öğrenme kursları__

    ---

    **Nextflow topluluk eğitim portalına hoş geldiniz!**

    Aşağıda listelenen eğitim kursları, kendi kendine öğrenme kaynağı olarak kullanılabilecek şekilde tasarlanmıştır.
    İstediğiniz zaman, Github Codespaces aracılığıyla sağladığımız web tabanlı ortamda veya kendi ortamınızda bu kursları kendi başınıza tamamlayabilirsiniz.

    [Kursları keşfedin :material-arrow-right:](#nextflow-egitim-kurslari-katalogu){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Ek bilgiler__

    ---

    ??? warning "Sürüm uyumluluğu"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **Ocak 2026 itibarıyla, aksi belirtilmedikçe tüm Nextflow eğitim kurslarımız Nextflow sürüm 25.10.2 veya üzeri ve katı sözdizimi etkinleştirilmiş olarak gerektirir.**

        Sürüm gereksinimleri ve katı sözdizimi hakkında daha fazla bilgi için lütfen [Nextflow dokümanları geçiş kılavuzuna](https://nextflow.io/docs/latest/strict-syntax.html) bakınız.

        Önceki sözdizimine karşılık gelen eğitim materyalinin eski sürümleri, bu web sayfasının menü çubuğundaki sürüm seçici aracılığıyla mevcuttur.

    ??? terminal "Ortam seçenekleri"

        Eğitim için ihtiyacınız olan her şeyin önceden yüklenmiş olduğu, Github Codespaces aracılığıyla erişilebilen web tabanlı bir eğitim ortamı sağlıyoruz (ücretsiz bir GitHub hesabı gerektirir).

        [![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Bu ihtiyaçlarınıza uygun değilse, lütfen diğer [Ortam seçeneklerine](./envsetup/index.md) bakınız.

    ??? learning "Eğitim etkinlikleri"

        Nextflow eğitimini yapılandırılmış bir etkinliğin parçası olarak almayı tercih ediyorsanız, bunu yapmanız için birçok fırsat bulunmaktadır. Aşağıdaki seçeneklere göz atmanızı öneririz:

        - Topluluk ekibi tarafından üç ayda bir düzenlenen **[Eğitim Haftaları]()**
        - **[Seqera Etkinlikleri](https://seqera.io/events/)** Seqera tarafından düzenlenen yüz yüze eğitim etkinliklerini içerir ('Seqera Sessions' ve 'Nextflow Summit' için arama yapın)
        - **[Nextflow Elçileri]()** yerel toplulukları için etkinlikler düzenler
        - **[nf-core etkinlikleri](https://nf-co.re/events)** topluluk hackathon'larını içerir

    ??? people "Eğitmenler için bilgiler"

        Kendi eğitimlerinizi yürüten bir eğitmen iseniz, uygun şekilde atıfta bulunduğunuz sürece materyallerimizi doğrudan eğitim portalından kullanabilirsiniz. Ayrıntılar için aşağıdaki 'Katkılar ve teşekkürler' bölümüne bakınız.

        Ayrıca, eğitim çabalarınızı nasıl daha iyi destekleyebileceğimiz konusunda sizden haber almak isteriz! Lütfen [community@seqera.io](mailto:community@seqera.io) adresinden veya topluluk forumundan bizimle iletişime geçin ([Yardım](help.md) sayfasına bakınız).

    ??? licensing "Açık kaynak lisansı ve katkı politikası"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Bu eğitim materyali [Seqera](https://seqera.io) tarafından geliştirilip sürdürülmekte ve topluluğun yararına açık kaynak lisansı ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) altında yayınlanmaktadır. Bu materyali lisans kapsamı dışında bir şekilde kullanmak istiyorsanız (ticari kullanım ve yeniden dağıtım konusundaki sınırlamalara dikkat ediniz), lütfen talebinizi görüşmek için [community@seqera.io](mailto:community@seqera.io) adresinden bizimle iletişime geçiniz.

        Topluluktan gelen iyileştirmeleri, düzeltmeleri ve hata raporlarını memnuniyetle karşılıyoruz. Her sayfanın sağ üst köşesinde, eğitim kaynak materyaline sorun bildirebileceğiniz veya pull request aracılığıyla değişiklik önerebileceğiniz kod deposuna bağlanan bir :material-file-edit-outline: simgesi bulunmaktadır. Daha fazla ayrıntı için depodaki `README.md` dosyasına bakınız.

</div>

!!! note "Yapay Zeka Destekli Çeviri"

    Bu çeviri, insan gözetiminde yapay zeka kullanılarak oluşturulmuştur.
    Geri bildirimlerinizi ve iyileştirme önerilerinizi memnuniyetle karşılıyoruz.
    [Çeviri kılavuzuna](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) bakınız.

## Nextflow eğitim kursları kataloğu

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Başlangıç seviyesi__

    ---

    ### :material-compass:{.nextflow-primary} Yeni Başlayanlar için Nextflow {.mt-1}

    Nextflow'a tamamen yeni olanlar için tasarlanmış, alana özgü olmayan kurslar. Her kurs, öğrencilerin becerilerini aşamalı olarak geliştirmelerine yardımcı olmak için tasarlanmış bir dizi eğitim modülünden oluşur.

    ??? courses "**Hello Nextflow:** Kendi boru hatlarınızı geliştirmeyi öğrenin"

        Bu kurs, basit ama tamamen işlevsel boru hatları geliştirmeyi sağlayacak kadar ayrıntılı olarak Nextflow dilinin temel bileşenlerini ve ayrıca boru hattı tasarımı, geliştirme ve yapılandırma uygulamalarının temel unsurlarını kapsar.

        [Hello Nextflow eğitimine başlayın :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Mevcut boru hatlarını çalıştırmayı öğrenin"

        Nextflow boru hatlarını çalıştırma ve yapılandırmaya yönelik kısa bir giriş, Hello Nextflow geliştirici kursuna dayanır ancak kod üzerinde daha az odaklanır. Yürütme, çıktılar, temel kod yapısı ve farklı hesaplama ortamları için yapılandırmayı kapsar.

        [Nextflow Run eğitimine başlayın :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Bilim için Nextflow {.mt-1}

    'Hello Nextflow'da sunulan kavram ve bileşenleri belirli bilimsel kullanım durumlarına uygulamayı öğrenin.

    ??? courses "**Genomik için Nextflow** (varyant çağırma)"

        Kendi genomik boru hatlarını geliştirmeyi öğrenmek isteyen araştırmacılar için. Kurs, basit ama işlevsel bir genomik boru hattının nasıl geliştirileceğini göstermek için bir varyant çağırma kullanım durumu kullanır.

        [Genomik için Nextflow eğitimine başlayın :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**RNAseq için Nextflow** (bulk RNAseq)"

        Kendi RNAseq boru hatlarını geliştirmeyi öğrenmek isteyen araştırmacılar için. Kurs, basit ama işlevsel bir RNAseq boru hattının nasıl geliştirileceğini göstermek için bir bulk RNAseq işleme kullanım durumu kullanır.

        [RNAseq için Nextflow eğitimine başlayın :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Görüntüleme için Nextflow** (uzamsal omikler)"

        Analiz boru hatlarını çalıştırmayı ve özelleştirmeyi öğrenmek isteyen görüntüleme ve uzamsal omikler alanındaki araştırmacılar için. Kurs, Nextflow boru hatlarını nasıl çalıştıracağınızı, yapılandıracağınızı ve girdileri nasıl yöneteceğinizi göstermek için biyolojik olarak ilgili bir boru hattı sağlamak amacıyla nf-core/molkart boru hattını kullanır.

        [Görüntüleme için Nextflow eğitimine başlayın :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __İleri seviye__

    ---

    ### :material-bridge:{.nextflow-primary} Nextflow'dan nf-core'a {.mt-1}

    [nf-core](https://nf-co.re/) topluluk projesinden kod ve en iyi uygulamaları kullanmayı öğrenin.

    Bu kurslar, Nextflow temellerinden nf-core en iyi uygulamalarına geçmenize yardımcı olur.
    nf-core topluluğunun boru hatlarını nasıl ve neden oluşturduğunu ve bu tekniklere nasıl katkıda bulunabileceğinizi ve bunları nasıl yeniden kullanabileceğinizi anlayın.

    ??? courses "**Hello nf-core:** nf-core ile başlayın"

        [nf-core](https://nf-co.re/) uyumlu boru hatlarını çalıştırmayı ve geliştirmeyi öğrenmek isteyen geliştiriciler için. Kurs, nf-core şablonunu ve geliştirme en iyi uygulamalarını takip eden basit ama tamamen işlevsel boru hatları geliştirmeyi ve mevcut nf-core modüllerini kullanmayı sağlayacak kadar ayrıntılı olarak nf-core boru hatlarının yapısını kapsar.

        [Hello nf-core eğitimine başlayın :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} İleri Düzey Nextflow Eğitimi {.mt-1}

    Gerçek dünya kullanım durumlarını ele almak için Nextflow boru hatlarını geliştirme ve dağıtma konusunda ileri düzey kavram ve mekanizmaları öğrenin.

    ??? courses "**Yan Görevler:** Bağımsız konulara derinlemesine dalışlar"

        Belirli konularda beceri aralıklarını genişletmek ve/veya derinleştirmek isteyen Nextflow geliştiricileri için tasarlanmış bağımsız mini kurslar. Doğrusal olarak sunulurlar ancak herhangi bir sırada alınabilirler (her mini kurs genel bakışındaki bağımlılıklara bakınız).

        [Yan Görevlere göz atın :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Eğitim Koleksiyonları:** Yan Görevler arasında önerilen öğrenme yolları"

        Eğitim Koleksiyonları, belirli bir tema veya kullanım durumu etrafında kapsamlı bir öğrenme deneyimi sağlamak için birden fazla Yan Görevi birleştirir.

        [Eğitim Koleksiyonlarına göz atın :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "Arşivlenmiş eğitim materyallerini mi arıyorsunuz?"

    Eski eğitim materyalleri (Temel Eğitim, İleri Düzey Eğitim ve diğer deneysel kurslar), Nextflow 3.0 katı sözdizimi ile uyumsuz oldukları için eğitim portalından kaldırılmıştır.
    Bu materyallere erişmeniz gerekiyorsa, Ocak 2026 öncesindeki [git geçmişinde](https://github.com/nextflow-io/training) mevcutturlar.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
