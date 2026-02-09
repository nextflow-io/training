# Manuel kurulum

Eğitim için ihtiyacınız olan her şeyi kendi yerel ortamınıza manuel olarak kurmak mümkündür.

Burada, standart POSIX uyumlu sistemlerde (dizüstü bilgisayar gibi kişisel bir makine varsayarak) bunun nasıl yapılacağını belgeledik.
Bazı ayrıntıların özel sisteminize bağlı olarak farklı olabileceğini unutmayın.

!!! tip "İpucu"

    Devam etmeden önce, [Devcontainers yaklaşımını](03_devcontainer.md) kullanmayı düşündünüz mü?
    Manuel kurulum gerektirmeden gerekli tüm araçları ve bağımlılıkları sağlar.

## Genel yazılım gereksinimleri

Nextflow, Java yüklü herhangi bir POSIX uyumlu sistemde (Linux, macOS, Windows Subsystem for Linux, vb.) kullanılabilir.
Eğitim kurslarımızın birkaç ek gereksinimi vardır.

Toplamda, aşağıdaki yazılımların kurulu olması gerekir:

- Bash veya eşdeğer kabuk
- [Java 11 (veya daha yeni, 21'e kadar)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (veya daha yeni)
- [VSCode](https://code.visualstudio.com) ve [Nextflow eklentisi](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

VSCode uygulaması teknik olarak isteğe bağlıdır ancak kurslar üzerinde çalışırken ve genel olarak Nextflow geliştirme çalışmalarınız için kullanmanızı şiddetle öneririz.

Nextflow dokümantasyon kılavuzu, bu bağımlılıkların kurulumu için [Ortam kurulumu](https://www.nextflow.io/docs/latest/developer-env.html) altında talimatlar sağlar.

## Nextflow ve nf-core araçları

Aşağıda bağlantısı verilen makalelerde ayrıntılı olarak açıklandığı gibi, Nextflow'un kendisini ve nf-core araçlarını kurmanız gerekecektir:

- [Nextflow kurulumu](https://www.nextflow.io/docs/latest/install.html)
- [nf-core araçları](https://nf-co.re/docs/nf-core-tools/installation)

Nextflow için kendi kendine kurulum seçeneğini ve nf-core araçları için PyPI seçeneğini kullanmanızı öneririz.

!!! warning "Uyarı: Sürüm uyumluluğu"

    <!-- Bu içeriğe yapılan herhangi bir güncellemenin ana sayfaya kopyalanması gerekir -->
    **Ocak 2026 itibarıyla, tüm Nextflow eğitim kurslarımız, aksi belirtilmedikçe, katı v2 sözdizimi etkinleştirilmiş Nextflow sürüm 25.10.2 veya daha yenisini gerektirir.**

    Sürüm gereksinimleri ve katı v2 sözdizimi hakkında daha fazla bilgi için lütfen [Nextflow sürümleri](../info/nxf_versions.md) kılavuzuna bakın.

    Önceki sözdizimine karşılık gelen eğitim materyalinin eski sürümlerine bu web sayfasının menü çubuğundaki sürüm seçici aracılığıyla erişilebilir.

## Eğitim materyalleri

Eğitim materyallerini indirmenin en kolay yolu, tüm depoyu bu komutu kullanarak klonlamaktır:

```bash
git clone https://github.com/nextflow-io/training.git
```

Her kursun kendi dizini vardır.
Bir kurs üzerinde çalışmak için bir terminal penceresi açın (tercihen VSCode uygulamasının içinden) ve ilgili dizine `cd` ile gidin.

Daha sonra web sitesinde sağlanan kurs talimatlarını takip edebilirsiniz.
