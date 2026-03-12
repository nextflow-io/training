---
title: Nextflow for {DOMAIN}
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Write a linear workflow to apply {METHOD} to a single sample
    - Handle accessory files such as {ACCESSORY_FILES} appropriately
    - Leverage Nextflow's dataflow paradigm to parallelize per-sample processing
    - Implement multi-sample aggregation using relevant channel operators
  audience_prerequisites:
    - "**Audience:** This course is designed for researchers in {DOMAIN} and related fields who want to develop or customize data analysis pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts, and common {DOMAIN} file formats is assumed."
    - "**Prerequisites:** Foundational Nextflow concepts and tooling covered in [Hello Nextflow](../../hello_nextflow/)."
---

---

title: {DOMAIN} için Nextflow
hide:

- toc
  page_type: index_page
  index_type: course
  additional_information:
  technical_requirements: true
  learning_objectives: - Tek bir örneğe {METHOD} uygulamak için doğrusal bir iş akışı yazın - {ACCESSORY_FILES} gibi yardımcı dosyaları uygun şekilde ele alın - Örnek başına işlemeyi paralel hale getirmek için Nextflow'un veri akışı paradigmasından yararlanın - İlgili kanal operatörlerini kullanarak çok örnekli toplama işlemi uygulayın
  audience_prerequisites: - "**Hedef Kitle:** Bu kurs, veri analizi boru hatlarını geliştirmek veya özelleştirmek isteyen {DOMAIN} ve ilgili alanlardaki araştırmacılar için tasarlanmıştır." - "**Beceriler:** Komut satırına, temel betik yazma kavramlarına ve yaygın {DOMAIN} dosya formatlarına belirli düzeyde aşinalık varsayılmaktadır." - "**Ön Koşullar:** [Hello Nextflow](../../hello_nextflow/) kursunda ele alınan temel Nextflow kavramları ve araçları."

---

# {DOMAIN} için Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow'u gerçek dünya {DOMAIN} kullanım senaryosuna uygulayan uygulamalı bir kurs: {METHOD_SHORT_DESCRIPTION}.**

Bu kurs, [Hello Nextflow](../../hello_nextflow/) başlangıç eğitimi üzerine inşa edilmiştir ve Nextflow'un {DOMAIN} alanına özgü bağlamda nasıl kullanılacağını göstermektedir.
[{TOOL_A}]({TOOL_A_URL}) ve [{TOOL_B}]({TOOL_B_URL}) araçlarıyla bir {METHOD} boru hattı uygulayacaksınız.

<!-- additional_information -->

## Kursa genel bakış

Bu kurs uygulamalıdır; bilgiyi kademeli olarak sunmak üzere yapılandırılmış, hedefe yönelik alıştırmalar içermektedir.

Metodolojiye hakim olmak için analiz araçlarını terminalde manuel olarak çalıştırarak başlayacak, ardından analizi otomatikleştiren ve ölçeklendiren bir Nextflow boru hattını aşamalı olarak oluşturacaksınız.

### Ders planı

Bu kursu, Nextflow'un {DOMAIN} kullanım senaryosuna uygulanmasının belirli yönlerine odaklanan üç bölüme ayırdık.

| Kurs bölümü                                          | Özet                                                                                                          | Tahmini süre |
| ---------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ------------ |
| [Bölüm 1: Yönteme genel bakış](./01_method.md)       | {METHOD} metodolojisini anlama ve araçları manuel olarak çalıştırma                                           | 30 dakika    |
| [Bölüm 2: Tek örnekli işleme](./02_single_sample.md) | {PART2_SUMMARY} yapan bir boru hattı oluşturma ve ardından birden fazla örneğe ölçeklendirme                  | 60 dakika    |
| [Bölüm 3: Çok örnekli toplama](./03_multi_sample.md) | Örnek başına çıktıları toplamak için kanal operatörlerini kullanarak çok örnekli {AGGREGATION_SUMMARY} ekleme | 45 dakika    |

Bu kursun sonunda, temel Nextflow kavramlarını ve araçlarını tipik bir {DOMAIN} kullanım senaryosuna uygulayabileceksiniz.

Kursa başlamaya hazır mısınız?

[Başlayın :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
