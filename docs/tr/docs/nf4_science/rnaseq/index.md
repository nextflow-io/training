---
title: RNAseq
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Temel RNAseq işleme ve QC yöntemlerini uygulamak için doğrusal bir iş akışı yazmak
    - FASTQ gibi alana özgü dosyaları ve genom referans kaynaklarını uygun şekilde ele almak
    - Tek uçlu ve çift uçlu dizileme verilerini işlemek
    - Nextflow'un veri akışı paradigmasından yararlanarak örnek başına RNAseq işlemesini paralelleştirmek
    - İlgili kanal operatörleri kullanarak birden fazla adım ve örnek genelinde QC raporlarını birleştirmek
  audience_prerequisites:
    - "**Hedef Kitle:** Bu kurs, veri analizi pipeline'ları geliştirmek veya özelleştirmek isteyen transkriptomik ve ilgili alanlardaki araştırmacılar için tasarlanmıştır."
    - "**Beceriler:** Komut satırı, temel betik kavramları ve yaygın RNAseq dosya formatları ile bir miktar aşinalık varsayılmaktadır."
    - "**Ön Koşullar:** [Hello Nextflow](../../hello_nextflow/) eğitiminde ele alınan temel Nextflow kavramları ve araçları."
---

# RNAseq için Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow'u gerçek dünya transkriptomik kullanım durumuna uygulayan uygulamalı bir kurs: Trim Galore, HISAT2 ve FastQC ile toplu RNAseq işleme.**

Bu kurs, [Hello Nextflow](../../hello_nextflow/) başlangıç eğitimi üzerine inşa edilmiş olup, Nextflow'un toplu RNAseq analizi bağlamında nasıl kullanılacağını göstermektedir.
Adaptör dizilerini kırpan, okumaları genom referansına hizalayan ve birkaç aşamada kalite kontrol (QC) gerçekleştiren bir işleme pipeline'ı uygulayacaksınız.

<!-- additional_information -->

## Kursa genel bakış

Bu kurs uygulamalı olup, bilgiyi kademeli olarak tanıtmak üzere yapılandırılmış hedef odaklı alıştırmalar içermektedir.

İşleme araçlarını terminalde manuel olarak çalıştırarak metodolojiyi anlamakla başlayacak, ardından analizi otomatikleştiren ve ölçeklendiren bir Nextflow pipeline'ını aşamalı olarak oluşturacaksınız.

### Ders planı

Bunu, her biri Nextflow'u bir RNAseq kullanım durumuna uygulamanın belirli yönlerine odaklanan üç bölüme ayırdık.

| Kurs bölümü                                                     | Özet                                                                                                                           | Tahmini süre |
| --------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ | ------------ |
| [Bölüm 1: Yönteme genel bakış](./01_method.md)                  | RNAseq işleme metodolojisini anlamak ve araçları manuel olarak çalıştırmak                                                     | 30 dakika    |
| [Bölüm 2: Tek örnekli uygulama](./02_single-sample.md)          | Tek bir örneği kırpan, hizalayan ve QC yapan bir pipeline oluşturmak, ardından birden fazla örneği işlemek için ölçeklendirmek | 60 dakika    |
| [Bölüm 3: Çok örnekli çift uçlu uygulama](./03_multi-sample.md) | Pipeline'ı çift uçlu veriyi işlemek ve örnekler genelinde QC raporlarını birleştirmek için genişletmek                         | 45 dakika    |

Bu kursun sonunda, temel Nextflow kavramlarını ve araçlarını tipik bir RNAseq kullanım durumuna uygulayabileceksiniz.

Kursa başlamaya hazır mısınız?

[Başlayın :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
