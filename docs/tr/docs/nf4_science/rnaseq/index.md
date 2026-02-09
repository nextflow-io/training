---
title: RNAseq için Nextflow
hide:
  - toc
---

# RNAseq için Nextflow

Bu eğitim kursu, transkriptomik ve ilgili alanlarda veri analiz boru hatlarını geliştirmek veya özelleştirmek isteyen araştırmacılar için tasarlanmıştır.
[Hello Nextflow](../../hello_nextflow/) başlangıç eğitimi üzerine inşa edilmiştir ve Nextflow'un bulk RNAseq analizi bağlamında nasıl kullanılacağını gösterir.

Özellikle, bu kurs adaptör dizilerini kırpmak, okumaları bir genom referansına hizalamak ve çeşitli aşamalarda kalite kontrolü (QC) gerçekleştirmek için basit bir bulk RNAseq işleme boru hattının nasıl uygulanacağını gösterir.

Hadi başlayalım! Eğitim ortamını başlatmak için aşağıdaki "Open in GitHub Codespaces" düğmesine tıklayın (tercihen ayrı bir sekmede), ardından yüklenirken okumaya devam edin.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Öğrenme hedefleri

Bu kursu tamamlayarak, temel Nextflow kavramlarını ve araçlarını tipik bir RNAseq kullanım senaryosuna nasıl uygulayacağınızı öğreneceksiniz.

Bu atölye çalışmasının sonunda şunları yapabileceksiniz:

- Temel RNAseq işleme ve QC yöntemlerini uygulamak için doğrusal bir iş akışı yazmak
- FASTQ ve referans genom kaynakları gibi alana özgü dosyaları uygun şekilde işlemek
- Tek uçlu ve çift uçlu dizileme verilerini işlemek
- Örnek başına RNAseq işlemeyi paralelleştirmek için Nextflow'un veri akışı paradigmasından yararlanmak
- İlgili kanal operatörlerini kullanarak birden fazla adım ve örnek genelinde QC raporlarını toplamak

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Ön koşullar

Kurs, aşağıdakilerle ilgili minimal bir aşinalık varsayar:

- Bu bilimsel alanda yaygın olarak kullanılan araçlar ve dosya formatları
- Komut satırı deneyimi
- [Hello Nextflow](../../hello_nextflow/) başlangıç eğitiminde ele alınan temel Nextflow kavramları ve araçları

Teknik gereksinimler ve ortam kurulumu için [Ortam Kurulumu](../../envsetup/) mini kursuna bakın.
