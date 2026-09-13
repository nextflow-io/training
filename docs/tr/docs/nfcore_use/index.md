---
title: nf-core Kullanımı
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - nf-core topluluk pipeline'larını bulun, indirin ve çalıştırın
    - Parametreler ve yapılandırma dosyaları kullanarak pipeline çalıştırmasını yapılandırın
    - nf-core pipeline'larının parametreleri ve girdi verilerini nasıl doğruladığını anlayın
    - Üretim ölçeğinde bir pipeline (nf-core/rnaseq) çalıştırın ve varsayılan kaynak tahsislerini geçersiz kılın
  audience_prerequisites:
    - "**Hedef Kitle:** Bu kurs, yerel Nextflow pipeline'larını nasıl çalıştıracağını bilen ve nf-core'a yeni olan, mevcut topluluk pipeline'larını çalıştırmak isteyen öğrenciler için tasarlanmıştır."
    - "**Beceriler:** Komut satırına, temel betik yazma kavramlarına ve yaygın dosya formatlarına belirli düzeyde aşinalık varsayılmaktadır."
    - "**Kurslar:** [Nextflow Run](../nextflow_run/index.md) kursunu tamamlamış olmanız ya da `nextflow run` ile yerel bir pipeline çalıştırma konusunda rahat olmanız gerekmektedir."
    - "**Alan:** Alıştırmalar biyoinformatik pipeline'ları kullanmaktadır; ancak önceden bilimsel alan bilgisi gerekmemektedir."
---

# nf-core Kullanımı

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**nf-core Kullanımı, nf-core topluluk pipeline'larını bulmaya, çalıştırmaya ve yapılandırmaya yönelik uygulamalı bir giriş kursudur.**

Pratik örnekler ve rehberli alıştırmalar üzerinden çalışarak nf-core pipeline'larını bulmayı ve indirmeyi, yerleşik test profilleri kullanarak çalıştırmayı ve parametreler ile yapılandırma dosyaları aracılığıyla çalıştırmayı özelleştirmeyi öğreneceksiniz.

Bu kursu tamamladığınızda, kendi analizleriniz için nf-core pipeline'larını çalıştırmaya başlayacak beceri ve özgüvene sahip olacaksınız.

<!-- additional_information -->

## Kursa genel bakış

Bu kurs uygulamalıdır; bilgiyi kademeli olarak tanıtmak üzere yapılandırılmış, hedefe yönelik alıştırmalar içermektedir.

Eğitim amaçlı nf-core projesi tarafından sürdürülen minimal bir pipeline olan `nf-core/demo` ile başlayacak, ardından öğrendiklerinizi toplu RNA dizileme analizi için yaygın olarak kullanılan bir üretim pipeline'ı olan `nf-core/rnaseq`'e uygulayacaksınız.

Bu kurs, pipeline çalıştırmaya odaklanmaktadır.
nf-core uyumlu pipeline geliştirmeye giriş arıyorsanız [Build with nf-core](../nfcore_build/index.md) sayfasına bakınız.

### Ders planı

| Kurs bölümü                                                          | Özet                                                                                                          | Tahmini süre |
| -------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ------------ |
| [Bölüm 1: Demo pipeline çalıştırma](./01_run_demo.md)                | Bir nf-core pipeline'ı bulun, indirin ve test profili kullanarak çalıştırın                                   | 20 dakika    |
| [Bölüm 2: Pipeline çalıştırmasını yapılandırma](./02_configure_execution.md) | Parametreleri ayarlayın, doğrulamayı anlayın ve kaynak tahsisi ile araç argümanlarını özelleştirin    | 20 dakika    |
| [Bölüm 3: Üretim pipeline'ı çalıştırma](./03_run_production_pipeline.md) | nf-core/rnaseq'i indirin ve çalıştırın; varsayılan kaynak tahsislerini geçersiz kılın                    | 20 dakika    |

Bu kursun sonunda nf-core projesinin sunduğu zengin topluluk pipeline'larından yararlanabileceksiniz.

Kursa başlamaya hazır mısınız?

[Öğrenmeye başlayın :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
