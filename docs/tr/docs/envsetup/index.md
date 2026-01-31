---
title: Ortam seçenekleri
description: Nextflow eğitimleri için ortamınızı ayarlama seçenekleri
hide:
  - toc
  - footer
---

# Ortam seçenekleri

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Öğrencilerin yazılım yönetimi ile zaman ve çaba harcamak yerine Nextflow öğrenmeye odaklanmalarını sağlayan tutarlı ve kapsamlı olarak test edilmiş bir ortam sağlamayı hedefliyoruz.
Bu amaçla, tüm kurslarımızda çalışmak için gerekli tüm yazılımları, kod dosyalarını ve örnek verileri içeren konteynerize bir ortam geliştirdik.

Bu konteynerize ortam, Github Codespaces üzerinde hazır olarak veya Devcontainers eklentisi ile VS Code'da yerel olarak çalıştırılabilir.

<div class="grid cards" markdown>

- :material-cloud-outline:{ .lg .middle } **Github Codespaces**

  ***

  GitHub Codespaces, buluttaki sanal makineler tarafından desteklenen, tüm araçlar ve veriler dahil önceden oluşturulmuş bir eğitim ortamı sağlamamıza olanak tanıyan web tabanlı bir hizmettir. Github hesabı olan herkes için ücretsiz olarak erişilebilir.

  [Github Codespaces kullanın :material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

- :material-laptop:{ .lg .middle } **Yerel Devcontainer'lar**

  ***

  Devcontainer'lı VS Code, tüm eğitim araçları önceden yapılandırılmış, yerel olarak çalışan konteynerize bir geliştirme ortamı sağlar. Codespaces ile aynı önceden oluşturulmuş ortamı sunar ancak tamamen yerel donanımınızda çalışır.

  [Devcontainer'ları yerel olarak kullanın :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Manuel kurulum talimatları

Yukarıdaki seçeneklerin hiçbiri ihtiyaçlarınıza uymuyorsa, yazılım bağımlılıklarını manuel olarak yükleyerek ve eğitim deposunu klonlayarak bu ortamı kendi yerel sisteminizde çoğaltabilirsiniz.

[Manuel kurulum :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Gitpod'un kullanımdan kaldırılması"

    Nextflow Eğitimi, Şubat 2025'e kadar [Gitpod](https://gitpod.io) kullanıyordu.
    Ancak Gitpod yapımcıları, ücretsiz işlevselliği [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex) sistemi lehine emekliye ayırmaya karar verdi.
    Bu nedenle, önceden kurulum gerektirmeden tek tıklamayla geliştirici ortamı sunan GitHub Codespaces'e geçtik.

    Gitpod'a ne zaman kaydolduğunuza ve hizmeti tam olarak ne zaman emekliye ayıracaklarına bağlı olarak, eski bulut IDE'lerinde eğitimi başlatabilirsiniz, ancak ileriye dönük güvenilir erişimi garanti edemeyiz:
    [Gitpod'da aç](https://gitpod.io/#https://github.com/nextflow-io/training).
