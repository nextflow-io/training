```yaml
---
title: Görüntüleme için Nextflow Çalıştırma
hide:
  - toc
---

# Görüntüleme için Nextflow Çalıştırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu eğitim kursu, veri analizi pipeline'larını çalıştırmak ve özelleştirmekle ilgilenen görüntüleme ve uzamsal biyoloji alanındaki araştırmacılar için tasarlanmıştır.
Molecular Cartography uzamsal transkriptomik verilerini işlemek için kullanılan bir pipeline olan [nf-core/molkart](https://nf-co.re/molkart) kullanarak workflow'ları çalıştırma, düzenleme ve yapılandırma ile ilgili temel Nextflow kavramlarını öğretir.
Burada öğrendiğiniz beceriler herhangi bir Nextflow veya nf-core pipeline'ına aktarılabilir.

Hadi başlayalım! Eğitim ortamını başlatmak için aşağıdaki "Open in GitHub Codespaces" düğmesine tıklayın (tercihen ayrı bir sekmede), ardından yüklenirken okumaya devam edin.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Öğrenme hedefleri

Bu kursu tamamlayarak, görüntüleme analizi pipeline'larını çalıştırmak için temel Nextflow kavramlarını ve araçlarını nasıl uygulayacağınızı öğreneceksiniz.

Bu workshop'un sonunda şunları yapabileceksiniz:

- Bir Nextflow workflow'unu yerel olarak başlatmak ve yürütmeyi izlemek
- Nextflow tarafından oluşturulan çıktıları (sonuçları) ve log dosyalarını bulmak ve yorumlamak
- Bir nf-core pipeline'ını test verileri ve özel girdilerle çalıştırmak
- Profiller ve parametre dosyalarını kullanarak pipeline yürütmesini yapılandırmak
- Örnek sayfaları (samplesheet) ve komut satırı parametrelerini kullanarak girdileri yönetmek

## Hedef kitle ve ön koşullar

Bu kurs aşağıdakiler konusunda minimal bir aşinalık varsayar:

- Komut satırı deneyimi
- Görüntüleme dosya formatları hakkında temel bilgi (TIFF görüntüleri, tablo verileri)

Teknik gereksinimler ve ortam kurulumu için [Ortam Kurulumu](../../envsetup/) mini kursuna bakın.
```
