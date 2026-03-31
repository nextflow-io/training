---
title: Nextflow sürümleri
description: Nextflow'un sözdizimi sürümlerinin evrimini anlama ve yönetme
hide:
  - toc
  - footer
---

## Desteklenen güncel Nextflow sözdizimi sürümü ve gereksinimleri

Eğitim portalının 3.0 sürümü itibarıyla, kurs indeks sayfasında aksi belirtilmedikçe tüm eğitim kurslarımız Nextflow'un 25.10.2 sürümünü temel almaktadır (sürüm bildirimi içermeyebilecek kullanımdan kaldırılmış veya arşivlenmiş materyaller bu kapsam dışındadır).

Kurslar artık iş akışı düzeyinde tiplendirilmiş girdiler ve iş akışı düzeyinde çıktı yönergeleri kullandığından, V2 sözdizimi ayrıştırıcısının kullanılması gerekmektedir.
[Github Codespaces](../envsetup/01_setup.md) veya [yerel devcontainer'lar](../envsetup/03_devcontainer.md) aracılığıyla sağladığımız ortamı kullanmayı planlıyorsanız, kurs talimatlarında özellikle belirtilmedikçe herhangi bir işlem yapmanıza gerek yoktur.
Ancak eğitimleri kendi ortamınızda ([Manuel kurulum](../envsetup/02_local.md)) tamamlamayı planlıyorsanız, v2 sözdizimi ayrıştırıcısı etkinleştirilmiş şekilde Nextflow 25.10.2 veya daha yeni bir sürümü kullandığınızdan emin olmanız gerekmektedir.

## Eğitim materyallerinin eski sürümleri

Eğitim materyallerimiz Şubat 2025'ten itibaren sürümlendirilmektedir.

**25.10.2 öncesi** Nextflow sürümleriyle çalışan eğitim materyallerinin eski sürümlerine, her sayfanın üst kısmındaki eğitim materyallerinin numaralı sürümünü gösteren açılır menü öğesi aracılığıyla erişebilirsiniz.
Eğitim materyallerinin eski bir sürümünü seçtiğinizde, eğitim ortamına yönelik bağlantılar otomatik olarak ortamın ilgili sürümünü belirtecektir.

## Nextflow sözdizimi sürümleri hakkında diğer bilgiler

Nextflow'un zaman zaman birbiriyle karıştırılan iki farklı sürümlendirme kavramı bulunmaktadır: **DSL sürümleri** ve **sözdizimi ayrıştırıcı sürümleri**.

**DSL1 ile DSL2**, Nextflow pipeline'ları yazmanın temelden farklı iki yolunu ifade eder.
DSL1, süreçlerin kanallar aracılığıyla örtük olarak birbirine bağlandığı özgün sözdizimidir.
Nextflow 20.07'de tanıtılan DSL2 ise modülerlik özellikleri ekledi: diğer dosyalardan süreç ve iş akışı içe aktarabilme, açık `workflow` blokları ve adlandırılmış süreç çıktıları.
DSL1, Nextflow 22.03'te kullanımdan kaldırıldı ve 22.12'de tamamen kaldırıldı.
Tüm modern Nextflow kodları DSL2 kullanmaktadır.

**Sözdizimi ayrıştırıcı v1 ile v2**, her ikisi de DSL2 koduyla çalışan farklı ayrıştırıcıları ifade eder.
v1 ayrıştırıcısı, özgün ve daha esnek ayrıştırıcıdır.
v2 ayrıştırıcısı daha katıdır ve statik tipleme (tiplendirilmiş girdiler ve çıktılar) ile iş akışı düzeyinde çıktı yönergeleri gibi yeni dil özelliklerini etkinleştirir.
v2 ayrıştırıcısı aynı zamanda daha iyi hata mesajları sağlar ve daha fazla hatayı çalışma zamanı yerine ayrıştırma aşamasında yakalar.
v2 ayrıştırıcısı, Nextflow 26.04'te varsayılan hale gelecektir.

Özetle: DSL2 yazdığınız dildir; sözdizimi ayrıştırıcı sürümü ise bu dilin ne kadar katı yorumlandığını ve hangi gelişmiş özelliklerin kullanılabilir olduğunu belirler.

### Nextflow sürümünü kontrol etme ve ayarlama

Sisteminizde yüklü olan Nextflow sürümünü `nextflow --version` komutunu kullanarak kontrol edebilirsiniz.

Nextflow sürümünüzü nasıl güncelleyeceğiniz hakkında daha fazla bilgi için lütfen [Nextflow'u Güncelleme](https://www.nextflow.io/docs/latest/updating-nextflow.html) başlıklı referans belgelerine bakınız.

### v2 sözdizimi ayrıştırıcısını etkinleştirme

Mevcut oturumunuz için v2 sözdizimi ayrıştırıcısını **etkinleştirmek** üzere terminalinizde aşağıdaki komutu çalıştırın:

```bash
export NXF_SYNTAX_PARSER=v2
```

Bunu kalıcı hale getirmek için (v2'nin Nextflow 26.04'te varsayılan olmasını beklerken), export komutunu kabuk profilinize ekleyin (`~/.bashrc`, `~/.zshrc` vb.):

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

`NXF_SYNTAX_PARSER=v2` ortam değişkeninin geçici bir gereksinim olduğunu unutmayın.
Nextflow 26.04 ve sonrasında v2 ayrıştırıcısı varsayılan hale gelecek ve bu ayar artık gerekli olmayacaktır.

### v2 sözdizimi ayrıştırıcısını devre dışı bırakma

Mevcut oturumunuz için v2 sözdizimi ayrıştırıcısını **devre dışı bırakmak** üzere terminalinizde aşağıdaki komutu çalıştırın:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- 26.04 sonrasındaki sürümlerde devre dışı bırakmak mümkün olacak mı? -->

### Mevcut kodun taşınması

Mevcut kodunuzu Nextflow'un daha yeni sürümleriyle uyumlu hale getirmeye ilişkin rehberlik için lütfen referans belgelerindeki [Taşıma Notları](https://www.nextflow.io/docs/latest/migrations/index.html)'na bakınız.

Aşağıdaki iki makale, en son sürüme geçiş için özellikle yararlıdır:

- [İş akışı çıktılarına geçiş](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Statik tiplere geçiş](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Bu özelliklerin her ikisi de eğitim materyallerinin 3.0 sürümünden itibaren başlangıç düzeyi eğitiminin bir parçası olarak ele alınmaktadır.

Taşımayı planladığınız Nextflow kodunun nesline bağlı olarak, işlemlerin büyük bölümünü `nextflow lint -format` komutunu kullanarak Nextflow linter aracılığıyla gerçekleştirebilirsiniz.
Daha fazla ayrıntı için [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) CLI referansına bakınız.

Bu bilgilerin işinize yaramasını umuyoruz.
Yardıma ihtiyaç duyarsanız Slack'te veya forumda bize ulaşın.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
