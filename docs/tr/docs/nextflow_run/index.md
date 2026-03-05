---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Nextflow workflow'larının çalıştırılması ve yönetimi
    - Çıktıların (sonuçlar) ve log dosyalarının bulunması ve yorumlanması
    - Basit çok adımlı bir workflow'da temel Nextflow bileşenlerinin tanınması
    - Pipeline çalıştırmasının HPC ve bulut dahil yaygın hesaplama platformlarında yapılandırılması
    - Kod modülerliği ve yazılım konteynerları dahil, pipeline'ları FAIR yapan tekrar üretilebilirlik, taşınabilirlik ve kod yeniden kullanımı için en iyi uygulamaların özetlenmesi
  audience_prerequisites:
    - "**Hedef Kitle:** Bu kurs, Nextflow'a tamamen yeni başlayan ve mevcut pipeline'ları çalıştırmak isteyen öğrenciler için tasarlanmıştır."
    - "**Beceriler:** Komut satırı, temel betik yazma kavramları ve yaygın dosya formatları hakkında bir miktar bilgi sahibi olunduğu varsayılmaktadır."
    - "**Alan:** Alıştırmaların tümü alandan bağımsızdır, bu nedenle önceden bilimsel bilgi gerekmez."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run, tekrar üretilebilir ve ölçeklenebilir veri analizi workflow'larını çalıştırmaya yönelik uygulamalı bir giriştir.**

Pratik örnekler ve rehberli alıştırmalar üzerinde çalışarak, pipeline'ları çalıştırma, dosyaları ve yazılım bağımlılıklarını yönetme, çalıştırmayı zahmetsizce paralelleştirme ve farklı hesaplama ortamlarında workflow'ları çalıştırma dahil olmak üzere Nextflow kullanmanın temellerini öğreneceksiniz.

Workflow'ları Nextflow ile çalıştırmaya başlamak için gerekli becerileri ve özgüveni kazanacaksınız.

<!-- additional_information -->

## Kurs genel bakış

### Ne yapacaksınız

Bu kurs, bilgileri kademeli olarak tanıtmak üzere yapılandırılmış hedefe yönelik alıştırmalarla uygulamalıdır.

Metin girdilerini işleyen bir Nextflow pipeline'ının birkaç versiyonunu çalıştıracaksınız.
Tek bir adımdan oluşan basit bir versiyonla başlayacak ve sonunda tablo şeklinde metin girdileri içeren bir CSV dosyası alan, birkaç dönüşüm adımı çalıştıran ve dönüştürülmüş metni söyleyen bir karakterin ASCII resmini içeren tek bir metin dosyası çıktı veren çok adımlı bir versiyona ilerleyeceksiniz.

Bu kurs, pipeline'ları çalıştırmaya odaklanır (temel `nextflow run` komutunun adını almıştır).
Nextflow pipeline'ları geliştirmeye giriş arıyorsanız, [Hello Nextflow](../hello_nextflow/index.md) bölümüne bakın.

### Ders planı

Bunu, Nextflow'da yazılmış pipeline'ları çalıştırmanın ve yönetmenin belirli yönlerine odaklanacak üç bölüme ayırdık.

| Kurs bölümü                                                  | Özet                                                                                                                              | Tahmini süre |
| ------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------- | ------------ |
| [Bölüm 1: Temel işlemleri çalıştırma](./01_basics.md)        | Basit bir workflow'un başlatılması ve çalıştırılmasının yönetilmesi                                                               | 30 dk        |
| [Bölüm 2: Gerçek pipeline'ları çalıştırma](./02_pipeline.md) | Karmaşık girdilerin işlenmesi, çok adımlı workflow'ların çalıştırılması, konteynerların kullanılması ve zahmetsiz paralelleştirme | 60 dk        |
| [Bölüm 3: Çalıştırma yapılandırması](./03_config.md)         | Pipeline davranışının özelleştirilmesi ve farklı hesaplama ortamlarında kullanımın optimize edilmesi                              | 60 dk        |

Bu kursun sonunda, bilimsel hesaplama ihtiyaçlarınız için tekrar üretilebilir workflow'ları çalıştırma yolculuğunuzdaki sonraki adımları atmaya hazır olacaksınız.

Kursu almaya hazır mısınız?

[Öğrenmeye başla :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
