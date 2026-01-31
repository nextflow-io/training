# Genomik için Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu eğitim kursu, veri analizi pipeline'larını geliştirmek veya özelleştirmek isteyen genomik ve ilgili alanlardaki araştırmacılar için tasarlanmıştır.
[Hello Nextflow](../../hello_nextflow/) başlangıç eğitimi üzerine inşa edilmiştir ve Nextflow'un genomik alanının özel bağlamında nasıl kullanılacağını gösterir.

Özellikle bu kurs, yüksek verimli dizileme verilerini analiz etmek için yaygın olarak kullanılan bir yazılım paketi olan [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit) ile basit bir varyant çağırma pipeline'ının nasıl uygulanacağını göstermektedir.

Haydi başlayalım! Eğitim ortamını başlatmak için aşağıdaki "Open in GitHub Codespaces" düğmesine tıklayın (tercihen ayrı bir sekmede), ardından yüklenirken okumaya devam edin.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Öğrenme hedefleri

Bu kursta çalışarak, temel Nextflow kavramlarını ve araçlarını tipik bir genomik kullanım senaryosuna nasıl uygulayacağınızı öğreneceksiniz.

Bu atölye çalışmasının sonunda şunları yapabileceksiniz:

- Tek bir örneğe varyant çağırma uygulamak için doğrusal bir workflow yazma
- Dizin dosyaları ve referans genom kaynakları gibi yardımcı dosyaları uygun şekilde işleme
- Örnek başına varyant çağırmayı paralelleştirmek için Nextflow'un veri akışı paradigmasından yararlanma
- İlgili channel operatörlerini kullanarak çok örnekli varyant çağırma uygulama
- Genomik-spesifik özgünlükleri uygun şekilde işleyen adım başına ve uçtan uca pipeline testleri uygulama

<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Ön koşullar

Kurs aşağıdakilerle ilgili minimal bir aşinalık varsayar:

- Bu bilimsel alanda yaygın olarak kullanılan araçlar ve dosya formatları
- Komut satırı deneyimi
- [Hello Nextflow](../../hello_nextflow/) başlangıç eğitiminde kapsanan temel Nextflow kavramları ve araçları

Teknik gereksinimler ve ortam kurulumu için [Ortam Kurulumu](../../envsetup/) mini-kursuna bakın.
