---
title: Hello nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - nf-core pipeline'larını almak, çalıştırmak ve yürütülmesini yönetmek
    - nf-core pipeline'larının kod yapısını ve proje organizasyonunu açıklamak
    - Bir şablondan temel nf-core uyumlu pipeline oluşturmak
    - Düz bir Nextflow iş akışını nf-core standartlarına uyacak şekilde yükseltmek
    - nf-core uyumlu bir pipeline'a nf-core modülleri eklemek
    - Kendi modüllerinizi nf-core'a katkıda bulunmak
    - nf-core araçlarını kullanarak girdileri ve parametreleri doğrulamak
  audience_prerequisites:
    - "**Hedef Kitle:** Bu kurs, temel Nextflow'a zaten aşina olan ve nf-core kaynaklarını ve en iyi uygulamalarını kullanmayı öğrenmek isteyen katılımcılar için tasarlanmıştır."
    - "**Beceriler:** Komut satırı, temel betik yazma kavramları ve yaygın dosya formatları hakkında bilgi sahibi olunduğu varsayılmaktadır."
    - "**Kurslar:** [Hello Nextflow](../hello_nextflow/index.md) kursunu veya eşdeğerini tamamlamış olmalıdır."
    - "**Alan:** Alıştırmaların tümü alandan bağımsızdır, dolayısıyla önceden bilimsel bilgi gerekmemektedir."
---

# Hello nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello nf-core, nf-core kaynaklarını ve en iyi uygulamalarını kullanmaya yönelik uygulamalı bir giriştir.**

![nf-core logo](./img/nf-core-logo.png#only-light)
![nf-core logo](./img/nf-core-logo-darkbg.png#only-dark)

Pratik örnekler ve rehberli alıştırmalar üzerinden çalışarak, nf-core uyumlu modüller ve pipeline'lar kullanmayı ve geliştirmeyi, ayrıca nf-core araçlarını etkili bir şekilde kullanmayı öğreneceksiniz.

nf-core en iyi uygulamalarına göre pipeline geliştirmeye başlamak için gereken becerileri ve özgüveni kazanacaksınız.

<!-- additional_information -->

## Kursa genel bakış

Bu kurs, bilgileri kademeli olarak tanıtmak üzere yapılandırılmış, hedefe yönelik alıştırmalar içeren, uygulamalı bir şekilde tasarlanmıştır.

Nextflow kullanılarak oluşturulmuş, küratörlüğü yapılmış bir bilimsel pipeline seti geliştirmek ve sürdürmek için bir topluluk çabası olan [**nf-core**](https://nf-co.re/)'un yanı sıra, açık geliştirme, test etme ve akran değerlendirmesini teşvik eden ilgili araçlar ve yönergeler ile tanıştırılacaksınız ([Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x), [Genome Biol 26, 228 (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9)).

nf-core topluluğu tarafından geliştirilen pipeline'lar modüler, ölçeklenebilir ve taşınabilir olacak şekilde tasarlanmıştır, bu da araştırmacıların kendi verilerini ve hesaplama kaynaklarını kullanarak bunları kolayca uyarlayıp çalıştırmalarına olanak tanır.
Proje tarafından uygulanan en iyi uygulama yönergeleri, pipeline'ların sağlam, iyi belgelenmiş ve gerçek dünya veri setlerine karşı doğrulanmış olmasını sağlar.
Bu, bilimsel analizlerin güvenilirliğini ve tekrarlanabilirliğini artırmaya yardımcı olur ve sonuçta araştırmacıların bilimsel keşiflerini hızlandırmasını sağlar.

Bu kursta nf-core pipeline'ları hakkında bilmeniz gereken her şeyi ele almayacağız, çünkü nf-core, topluluk tarafından yıllar içinde geliştirilen birçok özellik ve konvansiyonu kapsar.
Bunun yerine, başlamanıza ve nf-core'un nasıl çalıştığını anlamanıza yardımcı olacak temel kavramlara odaklanacağız.

### Ders planı

Bunu, nf-core kaynaklarını kullanmanın belirli yönlerine odaklanan beş bölüme ayırdık.

| Kurs bölümü                                                           | Özet                                                                                                                                                         | Tahmini süre |
| --------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------ |
| [Bölüm 1: Demo pipeline çalıştırma](./01_run_demo.md)                 | Mevcut bir nf-core pipeline'ı çalıştırın ve bu pipeline'ları temel Nextflow iş akışlarından farklı kılan şeyleri anlamak için kod yapısını inceleyin         | 30 dakika    |
| [Bölüm 2: Hello'yu nf-core için yeniden yazma](./02_rewrite_hello.md) | Mevcut bir iş akışını nf-core şablon iskeletine uyarlayın, [Hello Nextflow](../hello_nextflow/index.md) kursunda üretilen basit iş akışından başlayarak      | 60 dakika    |
| [Bölüm 3: nf-core modülü kullanma](./03_use_module.md)                | Topluluk modülleri kütüphanesini keşfedin ve yaygın biyoinformatik araçlarını sarmalayan önceden oluşturulmuş, test edilmiş modülleri entegre etmeyi öğrenin | 30 dakika    |
| [Bölüm 4: nf-core modülü oluşturma](./04_make_module.md)              | nf-core tarafından belirlenen belirli yapı, adlandırma konvansiyonları ve metadata gereksinimlerini kullanarak kendi nf-core tarzı modülünüzü oluşturun      | 30 dakika    |
| [Bölüm 5: Girdi doğrulama ekleme](./05_input_validation.md)           | nf-schema kullanarak hem komut satırı parametreleri hem de girdi veri dosyaları için girdi doğrulama uygulayın                                               | 30 dakika    |

Bu kursun sonunda, nf-core projesi tarafından sunulan muazzam kaynak zenginliğinden yararlanabileceksiniz.

Kursa başlamaya hazır mısınız?

[Öğrenmeye başlayın :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
