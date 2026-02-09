---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Nextflow iş akışlarının yürütülmesini başlatmak ve yönetmek
    - Çıktıları (sonuçları) ve günlük dosyalarını bulmak ve yorumlamak
    - Basit çok adımlı bir iş akışında temel Nextflow bileşenlerini tanımak
    - Pipeline yürütmesini HPC ve bulut dahil yaygın bilgi işlem platformlarında çalışacak şekilde yapılandırmak
    - Kod modülerliği ve yazılım konteynırları dahil olmak üzere, pipeline'ları FAIR yapan tekrarlanabilirlik, taşınabilirlik ve kod yeniden kullanımı için en iyi uygulamaları özetlemek
  audience_prerequisites:
    - "**Hedef Kitle:** Bu kurs, Nextflow'a tamamen yeni olan ve mevcut pipeline'ları çalıştırmak isteyen öğrenciler için tasarlanmıştır."
    - "**Beceriler:** Komut satırı, temel betik kavramları ve yaygın dosya formatları hakkında bir miktar aşinalık varsayılmaktadır."
    - "**Alan:** Tüm alıştırmalar alan-bağımsızdır, bu nedenle önceden bilimsel bilgi gerekmemektedir."
---

# Nextflow Run

**Nextflow Run, tekrarlanabilir ve ölçeklenebilir veri analizi iş akışlarını çalıştırmaya yönelik uygulamalı bir giriştir.**

Pratik örnekler ve rehberli alıştırmalar üzerinden çalışarak, pipeline'ları nasıl yürüteceğiniz, dosyaları ve yazılım bağımlılıklarını nasıl yöneteceğiniz, yürütmeyi zahmetsizce nasıl paralel hale getireceğiniz ve iş akışlarını farklı bilgi işlem ortamlarında nasıl çalıştıracağınız dahil olmak üzere Nextflow kullanımının temellerini öğreneceksiniz.

Nextflow ile iş akışlarını çalıştırmaya başlamak için gereken becerileri ve güveni kazanacaksınız.

<!-- additional_information -->

## Kursa genel bakış

### Neler yapacaksınız

Bu kurs uygulamalıdır ve bilgiyi kademeli olarak tanıtmak üzere yapılandırılmış hedefe yönelik alıştırmalar içerir.

Metin girdilerini işleyen bir Nextflow pipeline'ının birkaç sürümünü çalıştıracaksınız.
Tek bir adımdan oluşan basit bir sürümle başlayacak ve sonunda tablo halinde metin girdileri içeren bir CSV dosyasını alan, birkaç dönüştürme adımı çalıştıran ve dönüştürülmüş metni söyleyen bir karakterin ASCII resmini içeren tek bir metin dosyası çıktısı veren çok adımlı bir sürüme ilerleyeceksiniz.

Bu kurs, pipeline'ları çalıştırmaya odaklanır (temel `nextflow run` komutundan sonra adlandırılmıştır).
Nextflow pipeline'ları geliştirmeye yönelik bir giriş arıyorsanız, [Hello Nextflow](../hello_nextflow/index.md) bölümüne bakın.

### Ders planı

Bunu, Nextflow ile yazılmış pipeline'ları çalıştırma ve yönetmenin belirli yönlerine odaklanacak üç bölüme ayırdık.

| Kurs bölümü                                       | Özet                                                                                                                                | Tahmini süre |
| ------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- | ------------ |
| [Bölüm 1: Temel işlemleri çalıştırma](./01_basics.md) | Basit bir iş akışının yürütülmesini başlatma ve yönetme                                                                             | 30 dakika    |
| [Bölüm 2: Gerçek pipeline'ları çalıştırma](./02_pipeline.md) | Karmaşık girdileri işleme, çok adımlı iş akışlarını çalıştırma, konteyner kullanma ve yürütmeyi zahmetsizce paralel hale getirme | 60 dakika    |
| [Bölüm 3: Yapılandırmayı çalıştırma](./03_config.md) | Pipeline davranışını özelleştirme ve farklı hesaplama ortamlarında kullanımı optimize etme                                          | 60 dakika    |

Bu kursun sonunda, bilimsel hesaplama ihtiyaçlarınız için tekrarlanabilir iş akışlarını çalıştırma yolculuğunuzda sonraki adımlarla başa çıkmak için iyi hazırlanmış olacaksınız.

Kursa başlamaya hazır mısınız?

[Öğrenmeye başlayın :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
