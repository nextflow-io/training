---
title: Genomik için Nextflow
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Tek bir örneğe varyant çağırma uygulamak için doğrusal bir iş akışı yazın
    - İndeks dosyaları ve referans genom kaynakları gibi yardımcı dosyaları uygun şekilde yönetin
    - Örnek başına varyant çağırmayı paralelleştirmek için Nextflow'un veri akışı paradigmasından yararlanın
    - İlgili kanal operatörlerini kullanarak çok örnekli ortak çağırma uygulayın
  audience_prerequisites:
    - "**Hedef Kitle:** Bu kurs, veri analizi boru hatları geliştirmek veya özelleştirmek isteyen genomik ve ilgili alanlardaki araştırmacılar için tasarlanmıştır."
    - "**Beceriler:** Komut satırı, temel betik kavramları ve yaygın genomik dosya formatları hakkında bir miktar bilgi sahibi olunduğu varsayılmaktadır."
    - "**Ön Koşullar:** [Hello Nextflow](../../hello_nextflow/) kursunda ele alınan temel Nextflow kavramları ve araçları."
---

# Genomik için Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Gerçek dünya genomik kullanım senaryosuna Nextflow uygulaması: GATK ile varyant çağırma.**

Bu kurs, [Hello Nextflow](../../hello_nextflow/) başlangıç eğitimi üzerine inşa edilmiştir ve Nextflow'un genomik alanının özel bağlamında nasıl kullanılacağını gösterir.
Yüksek verimli dizileme verilerini analiz etmek için yaygın olarak kullanılan bir yazılım paketi olan [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit) ile bir varyant çağırma boru hattı uygulayacaksınız.

<!-- additional_information -->

## Kursa genel bakış

Bu kurs uygulamalıdır ve bilgiyi kademeli olarak tanıtmak üzere yapılandırılmış hedefe yönelik alıştırmalar içerir.

Metodolojiyi anlamak için terminalde varyant çağırma araçlarını manuel olarak çalıştırarak başlayacak, ardından analizi otomatikleştiren ve ölçeklendiren bir Nextflow boru hattını aşamalı olarak oluşturacaksınız.

### Ders planı

Bunu, her biri Nextflow'u bir genomik kullanım senaryosuna uygulamanın belirli yönlerine odaklanan üç bölüme ayırdık.

| Kurs bölümü                                                                 | Özet                                                                                                                    | Tahmini süre |
| --------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- | ------------ |
| [Bölüm 1: Metoda genel bakış](./01_method.md)                               | Varyant çağırma metodolojisini anlama ve araçları manuel olarak çalıştırma                                              | 30 dakika    |
| [Bölüm 2: Örnek başına varyant çağırma](./02_per_sample_variant_calling.md) | BAM dosyalarını indeksleyen ve varyantları çağıran bir boru hattı oluşturma, ardından birden fazla örneğe ölçeklendirme | 60 dakika    |
| [Bölüm 3: Kohort üzerinde ortak çağırma](./03_joint_calling.md)             | Örnek başına çıktıları toplamak için kanal operatörlerini kullanarak çok örnekli ortak genotipleme ekleme               | 45 dakika    |

Bu kursun sonunda, temel Nextflow kavramlarını ve araçlarını tipik bir genomik kullanım senaryosuna uygulayabileceksiniz.

Kursa başlamaya hazır mısınız?

[Başlayın :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
