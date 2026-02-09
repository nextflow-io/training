# Genomik için Nextflow

**Gerçek dünya genomik kullanım senaryosuna Nextflow uygulaması: GATK ile varyant çağırma üzerine uygulamalı bir kurs.**

Bu kurs, [Hello Nextflow](../../hello_nextflow/) başlangıç eğitimi üzerine inşa edilmiştir ve Nextflow'un genomik alanının özel bağlamında nasıl kullanılacağını gösterir.
Yüksek verimli dizileme verilerini analiz etmek için yaygın olarak kullanılan bir yazılım paketi olan [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit) ile bir varyant çağırma pipeline'ı uygulayacaksınız.

<!-- additional_information -->

## Kurs genel bakış

Bu kurs, bilgiyi kademeli olarak tanıtmak üzere yapılandırılmış hedef odaklı alıştırmalarla uygulamalıdır.

Metodolojini anlamak için önce terminalde varyant çağırma araçlarını manuel olarak çalıştırarak başlayacak, ardından analizi otomatikleştiren ve ölçeklendiren bir Nextflow pipeline'ını aşamalı olarak oluşturacaksınız.

### Ders planı

Bunu, her biri Nextflow'u bir genomik kullanım senaryosuna uygulamanın belirli yönlerine odaklanan üç bölüme ayırdık.

| Kurs bölümü                                                                 | Özet                                                                                                              | Tahmini süre |
| --------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- | ------------ |
| [Bölüm 1: Metod genel bakış](./01_method.md)                                | Varyant çağırma metodolojisini anlama ve araçları manuel olarak çalıştırma                                        | 30 dakika    |
| [Bölüm 2: Örnek başına varyant çağırma](./02_per_sample_variant_calling.md) | BAM dosyalarını indeksleyen ve varyant çağıran bir pipeline oluşturma, ardından birden fazla örneğe ölçeklendirme | 60 dakika    |
| [Bölüm 3: Kohort üzerinde joint calling](./03_joint_calling.md)             | Örnek başına çıktıları toplamak için kanal operatörlerini kullanarak çok örnekli joint genotyping ekleme          | 45 dakika    |

Bu kursun sonunda, temel Nextflow kavramlarını ve araçlarını tipik bir genomik kullanım senaryosuna uygulayabileceksiniz.

Kursa başlamaya hazır mısınız?

[Başlayın :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
