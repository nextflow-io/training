# Manuel kurulum

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Eğitimi çalıştırmak için ihtiyacınız olan her şeyi kendi yerel ortamınızda manuel olarak kurmak mümkündür.

Burada bunu standart POSIX uyumlu sistemlerde (dizüstü bilgisayar gibi kişisel bir makine varsayarak) nasıl yapacağınızı belgeledik.
Bazı ayrıntıların özel sisteminize bağlı olarak farklı olabileceğini unutmayın.

!!! tip "İpucu"

    Devam etmeden önce, [Devcontainer'lar yaklaşımını](03_devcontainer.md) düşündünüz mü?
    Manuel kurulum gerektirmeden gerekli tüm araçları ve bağımlılıkları sağlar.

## Genel yazılım gereksinimleri

Nextflow, Java yüklü herhangi bir POSIX uyumlu sistemde (Linux, macOS, Linux için Windows Alt Sistemi vb.) kullanılabilir.
Eğitim kurslarımızın birkaç ek gereksinimi vardır.

Toplamda, aşağıdaki yazılımların yüklü olması gerekir:

- Bash veya eşdeğer kabuk
- [Java 11 (veya sonrası, 21'e kadar)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (veya sonrası)
- [Nextflow eklentisi](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow) ile [VSCode](https://code.visualstudio.com)

VSCode uygulaması teknik olarak isteğe bağlıdır, ancak kurslar üzerinde çalışırken ve genel olarak Nextflow geliştirme çalışmalarınız için kullanmanızı şiddetle öneririz.

Nextflow dokümantasyon kılavuzu, bu bağımlılıkları yüklemek için [Ortam kurulumu](https://www.nextflow.io/docs/latest/developer-env.html) altında talimatlar sağlar.

## Nextflow ve nf-core araçları

Aşağıda bağlantısı verilen makalelerde ayrıntılı olarak açıklandığı gibi Nextflow'un kendisini ve nf-core araçlarını yüklemeniz gerekecektir:

- [Nextflow kurulumu](https://www.nextflow.io/docs/latest/install.html)
- [nf-core araçları](https://nf-co.re/docs/nf-core-tools/installation)

Nextflow için kendi kendine kurulum seçeneğini ve nf-core araçları için PyPI seçeneğini kullanmanızı öneririz.

!!! warning "Sürüm uyumluluğu"

    <!-- Any update to this content needs to be copied to the home page -->
    **Ocak 2026 itibarıyla, aksi belirtilmedikçe tüm Nextflow eğitim kurslarımız strict v2 syntax etkinleştirilmiş Nextflow sürüm 25.10.2 veya üstünü gerektirir.**

    Sürüm gereksinimleri ve strict v2 syntax hakkında daha fazla bilgi için lütfen [Nextflow sürümleri](../info/nxf_versions.md) rehberine bakın.

    Önceki söz dizimine karşılık gelen eğitim materyallerinin eski sürümlerine bu web sayfasının menü çubuğundaki sürüm seçici aracılığıyla ulaşabilirsiniz.

## Eğitim materyalleri

Eğitim materyallerini indirmenin en kolay yolu, bu komutu kullanarak tüm depoyu klonlamaktır:

```bash
git clone https://github.com/nextflow-io/training.git
```

Her kursun kendi dizini vardır.
Bir kurs üzerinde çalışmak için, bir terminal penceresi açın (ideal olarak VSCode uygulamasının içinden) ve ilgili dizine `cd` ile geçin.

Daha sonra web sitesinde verilen kurs talimatlarını takip edebilirsiniz.
