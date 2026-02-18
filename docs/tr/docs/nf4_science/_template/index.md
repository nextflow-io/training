---
title: {DOMAIN} için Nextflow
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Tek bir örneğe {METHOD} uygulamak için doğrusal bir iş akışı yazın
    - {ACCESSORY_FILES} gibi yardımcı dosyaları uygun şekilde yönetin
    - Örnek başına işlemeyi paralelleştirmek için Nextflow'un veri akışı paradigmasından yararlanın
    - İlgili kanal operatörlerini kullanarak çok örnekli toplama işlemini uygulayın
  audience_prerequisites:
    - "**Hedef Kitle:** Bu kurs, veri analizi boru hatları geliştirmek veya özelleştirmek isteyen {DOMAIN} ve ilgili alanlardaki araştırmacılar için tasarlanmıştır."
    - "**Beceriler:** Komut satırı, temel betik kavramları ve yaygın {DOMAIN} dosya formatları hakkında bir miktar bilgi sahibi olunduğu varsayılmaktadır."
    - "**Ön Koşullar:** [Hello Nextflow](../../hello_nextflow/) kursunda ele alınan temel Nextflow kavramları ve araçları."
---

# {DOMAIN} için Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow'u gerçek dünya {DOMAIN} kullanım senaryosuna uygulayan uygulamalı bir kurs: {METHOD_SHORT_DESCRIPTION}.**

Bu kurs, [Hello Nextflow](../../hello_nextflow/) başlangıç eğitimi üzerine inşa edilmiştir ve Nextflow'un {DOMAIN} alanının özel bağlamında nasıl kullanılacağını gösterir.
[{TOOL_A}]({TOOL_A_URL}) ve [{TOOL_B}]({TOOL_B_URL}) ile bir {METHOD} boru hattı uygulayacaksınız.

<!-- additional_information -->

## Kursa genel bakış

Bu kurs uygulamalıdır ve bilgiyi kademeli olarak tanıtmak üzere yapılandırılmış hedef odaklı alıştırmalar içerir.

Metodolojiyi anlamak için analiz araçlarını terminalde manuel olarak çalıştırarak başlayacak, ardından analizi otomatikleştiren ve ölçeklendiren bir Nextflow boru hattını aşamalı olarak oluşturacaksınız.

### Ders planı

Bunu, her biri Nextflow'u bir {DOMAIN} kullanım senaryosuna uygulamanın belirli yönlerine odaklanan üç bölüme ayırdık.

| Kurs bölümü                                          | Özet                                                                                                          | Tahmini süre |
| ---------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ------------ |
| [Bölüm 1: Metoda genel bakış](./01_method.md)        | {METHOD} metodolojisini anlama ve araçları manuel olarak çalıştırma                                           | 30 dakika    |
| [Bölüm 2: Tek örnekli işleme](./02_single_sample.md) | {PART2_SUMMARY} yapan bir boru hattı oluşturma, ardından birden fazla örneğe ölçeklendirme                    | 60 dakika    |
| [Bölüm 3: Çok örnekli toplama](./03_multi_sample.md) | Örnek başına çıktıları toplamak için kanal operatörlerini kullanarak çok örnekli {AGGREGATION_SUMMARY} ekleme | 45 dakika    |

Bu kursun sonunda, temel Nextflow kavramlarını ve araçlarını tipik bir {DOMAIN} kullanım senaryosuna uygulayabileceksiniz.

Kursa başlamaya hazır mısınız?

[Başlayın :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
