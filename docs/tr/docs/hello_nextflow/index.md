---
title: Hello Nextflow
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Nextflow iş akışlarının çalıştırılmasını başlatma ve yönetme
    - Nextflow tarafından oluşturulan çıktıları (sonuçlar) ve günlük dosyalarını bulma ve yorumlama
    - Temel sorunları giderme
    - Temel Nextflow bileşenlerinden basit bir çok adımlı iş akışı oluşturma
    - Temel kanal fabrikaları ve operatör türlerini ayırt etme ve bunları basit bir iş akışında etkili bir şekilde kullanma
    - HPC ve bulut dahil olmak üzere yaygın hesaplama platformlarında çalışacak şekilde iş akışı yürütmeyi yapılandırma
    - Kod modülerliği ve yazılım konteynerleri dahil olmak üzere iş akışlarını FAIR yapan tekrarlanabilirlik, taşınabilirlik ve kod yeniden kullanımı için en iyi uygulamaları uygulama
  audience_prerequisites:
    - "**Hedef kitle:** Bu kurs, Nextflow'a tamamen yeni olan ve kendi iş akışlarını geliştirmek isteyen öğrenciler için tasarlanmıştır."
    - "**Beceriler:** Komut satırı, temel betik kavramları ve yaygın dosya formatları ile biraz aşinalık varsayılmaktadır."
    - "**Alan:** Alıştırmaların tümü alana bağımlı değildir, bu nedenle önceden bilimsel bilgi gerekmez."
  videos_playlist: https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n
---

# Hello Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello Nextflow, tekrarlanabilir ve ölçeklenebilir veri analizi iş akışları oluşturmaya pratik bir giriştir.**

Pratik örnekler ve rehberli alıştırmalar üzerinde çalışarak, süreçleri tanımlama, bunları iş akışlarına bağlama, dosya ve yazılım bağımlılıklarını yönetme, yürütmeyi zahmetsizce paralelleştirme ve farklı hesaplama ortamlarında iş akışlarını çalıştırma dahil olmak üzere Nextflow ile iş akışları geliştirmenin temellerini öğreneceksiniz.

Nextflow ile kendi iş akışlarınızı geliştirmeye ve çalıştırmaya başlamak için gereken beceri ve özgüveni kazanacaksınız.

<!-- additional_information -->

## Kurs genel bakışı

Bu kurs, bilgileri kademeli olarak tanıtan hedef odaklı alıştırmalarla uygulamalı olacak şekilde tasarlanmıştır.

Bazı metin girdilerini alan, birkaç dönüştürme adımı çalıştıran ve dönüştürülmüş metni söyleyen bir karakterin ASCII resmini içeren tek bir metin dosyası çıktısı üreten basit bir Nextflow iş akışı geliştireceksiniz.

### Ders planı

Sizi kavramlar ve kodla bunaltmamak için, bunu her biri Nextflow ile iş akışları geliştirmenin belirli yönlerine odaklanan altı bölüme ayırdık.

| Kurs bölümü                                           | Özet                                                                                                               | Tahmini süre |
| ----------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ | ------------ |
| [Bölüm 1: Hello World](./01_hello_world.md)           | Bir Nextflow iş akışını birleştirme ve çalıştırmaya dahil olan temel bileşenler ve ilkeler                         | 30 dk        |
| [Bölüm 2: Hello Channels](./02_hello_channels.md)     | Girdileri işlemek ve yürütmeyi zahmetsizce paralelleştirmek için kanalları ve operatörleri kullanma                | 45 dk        |
| [Bölüm 3: Hello Workflow](./03_hello_workflow.md)     | Birden fazla adımı birbirine bağlamak ve adımlar arasında veri aktarımını yönetmek için kanalları kullanma         | 60 dk        |
| [Bölüm 4: Hello Modules](./04_hello_modules.md)       | Yeniden kullanılabilirliği artırmak ve bakım yükünü azaltmak için kod modülerliği ilkelerini uygulama              | 20 dk        |
| [Bölüm 5: Hello Containers](./05_hello_containers.md) | Yazılım bağımlılıklarını yönetmek ve tekrarlanabilirliği artırmak için konteynerleri bir mekanizma olarak kullanma | 60 dk        |
| [Bölüm 6: Hello Config](./06_hello_config.md)         | Farklı hesaplama ortamlarında iş akışı davranışını özelleştirme ve kullanımı optimize etme                         | 60 dk        |

Bu kursun sonunda, bilimsel hesaplama ihtiyaçlarınız için tekrarlanabilir iş akışları geliştirme yolculuğunuzdaki sonraki adımları atmaya iyi hazırlanmış olacaksınız.

Kursa başlamaya hazır mısınız?

[Başlayın :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
