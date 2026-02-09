# Ortam seçenekleri

Nextflow eğitimleri için ortamınızı kurma seçenekleri

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Öğrencilerin yazılım yönetimine zaman ve çaba harcamak zorunda kalmadan Nextflow öğrenmeye odaklanmalarını sağlayan tutarlı ve kapsamlı test edilmiş bir ortam sunmayı hedefliyoruz.
Bu amaçla, tüm kurslarımızda çalışmak için gerekli tüm yazılımları, kod dosyalarını ve örnek verileri içeren konteynerleştirilmiş bir ortam geliştirdik.

Bu konteynerleştirilmiş ortam, Github Codespaces üzerinde veya Devcontainers uzantısı ile VS Code'da yerel olarak kullanıma hazır şekilde çalıştırılabilir.

<div class="grid cards" markdown>

-   :material-cloud-outline:{ .lg .middle } __Github Codespaces__

    ---

    GitHub Codespaces, buluttaki sanal makineler tarafından desteklenen, tüm araçlar ve veriler dahil olmak üzere eğitim için önceden hazırlanmış bir ortam sunmamızı sağlayan web tabanlı bir hizmettir. Github hesabı olan herkes tarafından ücretsiz olarak erişilebilir.

    [Github Codespaces kullanın:material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

-   :material-laptop:{ .lg .middle } __Yerel Devcontainers__

    ---

    Devcontainers ile VS Code, tüm eğitim araçlarının önceden yapılandırıldığı, yerel olarak çalışan konteynerleştirilmiş bir geliştirme ortamı sağlar. Codespaces ile aynı önceden hazırlanmış ortamı sunar ancak tamamen yerel donanımınızda çalışır.

    [Devcontainers'ı yerel olarak kullanın :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Manuel kurulum talimatları

Yukarıdaki seçeneklerin hiçbiri ihtiyaçlarınıza uygun değilse, yazılım bağımlılıklarını manuel olarak yükleyerek ve eğitim deposunu klonlayarak bu ortamı kendi yerel sisteminizde çoğaltabilirsiniz.

[Manuel kurulum :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Gitpod'un Kullanımdan Kaldırılması"

    Nextflow Training, Şubat 2025'e kadar [Gitpod](https://gitpod.io) kullanıyordu.
    Ancak Gitpod'un yapımcıları, [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex) sistemini tercih ederek ücretsiz işlevselliği sonlandırmaya karar verdi.
    Bu nedenle, önceden kurulum gerektirmeyen tek tıklamayla geliştirici ortamı da sunan GitHub Codespaces kullanmaya geçtik.

    Gitpod'a ne zaman kaydolduğunuza ve hizmeti tam olarak ne zaman sonlandırdıklarına bağlı olarak, eğitimi eski bulut IDE'lerinde başlatabilirsiniz, ancak ileriye dönük güvenilir erişimi garanti edemeyiz:
    [Gitpod'da Aç](https://gitpod.io/#https://github.com/nextflow-io/training).
