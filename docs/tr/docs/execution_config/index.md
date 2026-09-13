---
title: Yürütme Yapılandırması
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Docker ve Conda arasında yazılım paketleme teknolojisini değiştirme
    - Bir yürütme platformu seçme ve Nextflow'un görev yürütmesini buna nasıl uyarladığını anlama
    - İşlem kaynağı tahsislerini kontrol etme ve başarısız olan görevleri otomatik olarak yeniden deneme
    - Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profiller tanımlama ve birleştirme
  audience_prerequisites:
    - "**Hedef Kitle:** Bu kurs, yerel Nextflow pipeline'larını nasıl başlatacağını bilen ve yürütmeyi daha ayrıntılı yapılandırmak isteyen öğrenenler için tasarlanmıştır."
    - "**Beceriler:** Komut satırına belirli düzeyde aşinalık varsayılmaktadır."
    - "**Kurslar:** [Nextflow Run](../nextflow_run/index.md) kursunu tamamlamış olmak ya da `nextflow run` ile yerel bir pipeline çalıştırma konusunda rahat olmak gerekmektedir."
---

# Yürütme Yapılandırması

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Yürütme Yapılandırması, Nextflow pipeline yürütmesini farklı hesaplama ortamlarına uyarlamaya yönelik uygulamalı bir girişdir.**

Hedefe yönelik alıştırmalar üzerinden çalışarak yazılım paketleme teknolojisini nasıl değiştireceğinizi, bir yürütme platformu nasıl seçeceğinizi, işlem kaynağı tahsislerini ve yeniden denemeleri nasıl kontrol edeceğinizi ve yapılandırmayı çalışma zamanında geçiş yapılabilir profiller halinde nasıl paketleyeceğinizi öğreneceksiniz.

Nextflow pipeline yürütmesini profesyonelce yapılandırmak için gereken beceri ve özgüveni kazanacaksınız.

<!-- additional_information -->

## Kursa genel bakış

Bu kurs uygulamalıdır ve [Nextflow Run](../nextflow_run/index.md) kursunda ele alınan beceriler üzerine inşa edilmektedir.

O kurstan alınan çok adımlı pipeline'ı ele alacak ve yapılandırmasını farklı hesaplama ortamlarına uyacak şekilde aşamalı olarak uyarlayacaksınız; ardından her şeyi çalışma zamanında geçiş yapabileceğiniz profiller halinde paketleyeceksiniz.

### Ders planı

| Kurs bölümü                                                                              | Özet                                                                                    | Tahmini süre |
| ---------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- | ------------ |
| [Bölüm 1: Hesaplama ortamınıza uyarlama](./01_packaging_and_execution.md)                | Yazılım paketleme teknolojisini değiştirme ve bir yürütme platformu seçme               | 20 dakika    |
| [Bölüm 2: İşlem kaynaklarını ve hataları yönetme](./02_resources_and_retries.md)         | Kaynak tahsislerini kontrol etme ve başarısız olan görevleri otomatik olarak yeniden deneme | 15 dakika    |
| [Bölüm 3: Yapılandırmaları değiştirmek için profil kullanma](./03_profiles.md)           | Profiller tanımlama ve birleştirme, tam olarak çözümlenmiş yapılandırmayı inceleme      | 15 dakika    |

Bu kursun sonunda, çeşitli hesaplama ortamları için Nextflow pipeline'larını yapılandırma ve aralarında minimum güçlükle geçiş yapma konusunda kendinize güveneceksiniz.

Kursa başlamaya hazır mısınız?

[Öğrenmeye başlayın :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
