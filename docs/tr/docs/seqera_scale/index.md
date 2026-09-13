---
title: Scale with Seqera
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Sign up for Seqera Platform and explore the Community Showcase
    - Add a pipeline to a workspace and launch it from the web interface
    - Authenticate and launch pipelines from the command line with the `tw` CLI
    - Register a GitHub-hosted pipeline and launch it both ways
  audience_prerequisites:
    - "**Audience:** This course is designed for learners who want to run Nextflow pipelines at scale using Seqera Platform."
    - "**Skills:** Familiarity with running nf-core pipelines from the command line is assumed."
    - "**Courses:** Must have completed [Nextflow Run](../nextflow_run/index.md) and [Use nf-core](../nfcore_use/index.md), or otherwise be comfortable running local and `nf-core/rnaseq` pipelines."
---

---
title: Seqera ile Ölçeklendirme
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Seqera Platform'a kaydolun ve Community Showcase'i keşfedin
    - Bir çalışma alanına pipeline ekleyin ve web arayüzünden başlatın
    - `tw` CLI ile kimlik doğrulaması yapın ve pipeline'ları komut satırından başlatın
    - GitHub'da barındırılan bir pipeline'ı kaydedin ve her iki yöntemle de başlatın
  audience_prerequisites:
    - "**Hedef Kitle:** Bu kurs, Seqera Platform kullanarak Nextflow pipeline'larını büyük ölçekte çalıştırmak isteyen öğrenenler için tasarlanmıştır."
    - "**Beceriler:** nf-core pipeline'larını komut satırından çalıştırma konusunda temel bilgiye sahip olunduğu varsayılmaktadır."
    - "**Kurslar:** [Nextflow Run](../nextflow_run/index.md) ve [Use nf-core](../nfcore_use/index.md) kurslarını tamamlamış ya da yerel ve `nf-core/rnaseq` pipeline'larını çalıştırma konusunda deneyimli olunması gerekmektedir."
---

# Seqera ile Ölçeklendirme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Seqera ile Ölçeklendirme, Seqera Platform kullanarak Nextflow pipeline'larını başlatmaya ve izlemeye yönelik uygulamalı bir giriş kursudur.**

Pratik örnekler üzerinden çalışarak Seqera Platform'a erişim kurabilecek, hem web arayüzünden hem de komut satırından üretim ölçeğinde bir pipeline başlatabilecek ve çalışma alanınıza yeni bir pipeline ekleyebileceksiniz.

Bu kursun sonunda Seqera Platform'da kendi pipeline'larınızı çalıştırma ve izleme konusunda gerekli becerileri ve özgüveni kazanmış olacaksınız.

<!-- additional_information -->

## Kursa genel bakış

Bu kurs uygulamalı bir yapıya sahiptir ve [Use nf-core](../nfcore_use/index.md) kursunda çalıştırdığınız pipeline'lar üzerine inşa edilmektedir.

Seqera Platform'a kaydolarak ve üretim ölçeğinde bir pipeline olan `nf-core/rnaseq`'i web arayüzünden başlatarak işe başlayacaksınız.
Ardından aynı işlemi bir terminalden gerçekleştirmek için `tw` komut satırı aracına geçecek ve son olarak yeni bir pipeline olan `nf-core/demo`'yu kaydedip her iki yöntemle de başlatacaksınız.

### Ders planı

| Kurs bölümü                                                                | Özet                                                                                                          | Tahmini süre |
| -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ------------ |
| [Bölüm 1: Pipeline'ları web arayüzünden başlatma](./01_run_with_seqera.md) | Seqera Platform erişimini kurun ve üretim ölçeğinde bir pipeline'ı web arayüzünden başlatın                   | 20 dakika    |
| [Bölüm 2: Pipeline'ları komut satırından başlatma](./02_launch_from_cli.md) | `tw` CLI ile kimlik doğrulaması yapın, kaydedilmiş bir pipeline başlatın ve CLI'dan yeni bir pipeline kaydedin | 25 dakika    |

Bu kursun sonunda, web arayüzünü mi yoksa komut satırını mı tercih ettiğinizden bağımsız olarak Seqera Platform'da Nextflow pipeline'larını başlatma ve izleme konusunda kendinize güveneceksiniz.

Kursa başlamaya hazır mısınız?

[Öğrenmeye başlayın :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
