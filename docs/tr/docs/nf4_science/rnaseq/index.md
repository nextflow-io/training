---
title: RNAseq için Nextflow
hide:
  - toc
---

# RNAseq için Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu eğitim kursu, transkriptomik ve ilgili alanlarda veri analizi pipeline'ları geliştirmek veya özelleştirmek isteyen araştırmacılar için tasarlanmıştır.
[Hello Nextflow](../../hello_nextflow/) başlangıç eğitimi üzerine inşa edilmiş olup, Nextflow'un toplu RNAseq analizi bağlamında nasıl kullanılacağını göstermektedir.

Özellikle, bu kurs adaptör dizilerini kırpmak, okumaları genom referansına hizalamak ve birkaç aşamada kalite kontrol (QC) gerçekleştirmek için basit bir toplu RNAseq işleme pipeline'ının nasıl uygulanacağını göstermektedir.

Hadi başlayalım! Eğitim ortamını başlatmak için aşağıdaki "Open in GitHub Codespaces" düğmesine tıklayın (tercihen ayrı bir sekmede), ardından yüklenirken okumaya devam edin.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Öğrenim hedefleri

Bu kursu tamamlayarak, temel Nextflow kavramlarını ve araçlarını tipik bir RNAseq kullanım durumuna nasıl uygulayacağınızı öğreneceksiniz.

Bu atölye sonunda şunları yapabileceksiniz:

- Temel RNAseq işleme ve QC yöntemlerini uygulamak için doğrusal bir workflow yazmak
- FASTQ gibi alana özgü dosyaları ve genom referans kaynaklarını uygun şekilde ele almak
- Tek uçlu ve çift uçlu dizileme verilerini işlemek
- Nextflow'un dataflow paradigmasından yararlanarak örnek başına RNAseq işlemesini paralelleştirmek
- İlgili channel operatörlerini kullanarak birden fazla adım ve örnek genelinde QC raporlarını birleştirmek

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Ön koşullar

Kurs, aşağıdakilerle ilgili minimal bir aşinalık varsaymaktadır:

- Bu bilimsel alanda yaygın olarak kullanılan araçlar ve dosya formatları
- Komut satırı deneyimi
- [Hello Nextflow](../../hello_nextflow/) başlangıç eğitiminde ele alınan temel Nextflow kavramları ve araçları.

Teknik gereksinimler ve ortam kurulumu için [Ortam Kurulumu](../../envsetup/) mini kursuna bakın.
