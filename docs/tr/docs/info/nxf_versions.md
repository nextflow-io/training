---
title: Nextflow sürümleri
description: Nextflow'un sözdizimi sürümlerinin evrimini anlama ve yönetme
hide:
  - toc
  - footer
---

## Mevcut desteklenen Nextflow sözdizimi sürümü ve gereksinimleri

Eğitim portalının 3.0 sürümü itibarıyla, kurs dizin sayfasında aksi belirtilmedikçe (kullanımdan kaldırılmış veya arşivlenmiş materyaller hariç, bunlar sürüm bildirimi içermeyebilir) tüm eğitim kurslarımız Nextflow'un 25.10.2 sürüm yayınına dayanmaktadır.

Kurslar artık iş akışı seviyesinde tiplendirilmiş girdilerin yanı sıra iş akışı seviyesinde çıktı yönergelerini kullandığından, V2 sözdizimi ayrıştırıcısının kullanılmasını gerektirmektedir.
[Github Codespaces](../envsetup/01_setup.md) veya [yerel devcontainer'lar](../envsetup/03_devcontainer.md) aracılığıyla sağladığımız ortamı kullanmayı planlıyorsanız, kurs talimatlarında özellikle belirtilmedikçe herhangi bir şey yapmanıza gerek yoktur.
Ancak, eğitimleri kendi ortamınızda ([Manuel kurulum](../envsetup/02_local.md)) çalışmayı planlıyorsanız, v2 sözdizimi ayrıştırıcısı etkinleştirilmiş Nextflow sürüm 25.10.2 veya daha yenisini kullandığınızdan emin olmanız gerekecektir.

## Eğitim materyallerinin eski sürümleri

Eğitim materyallerimiz Şubat 2025'ten beri sürümlendirilmektedir.

**25.10.2'den önceki** Nextflow sürümleriyle çalışan eğitim materyallerinin eski sürümlerine, her sayfanın üst kısmında eğitim materyallerinin numaralı sürümünü gösteren açılır menü öğesi aracılığıyla erişebilirsiniz.
Eğitim materyallerinin eski bir sürümünü seçtiğinizde, eğitim ortamına olan bağlantılar otomatik olarak ortamın ilgili sürümünü belirtecektir.

## Nextflow sözdizimi sürümleri hakkında diğer bilgiler

Nextflow'un bazen karıştırılan iki farklı sürümleme kavramı vardır: **DSL sürümleri** ve **sözdizimi ayrıştırıcısı sürümleri**.

**DSL1 vs DSL2**, Nextflow pipeline'larını yazmanın temelde farklı yollarını ifade eder.
DSL1, süreçlerin kanallar aracılığıyla örtük olarak bağlandığı orijinal sözdizimiydi.
Nextflow 20.07'de tanıtılan DSL2, modülerlik özellikleri ekledi: süreçleri ve iş akışlarını diğer dosyalardan içe aktarma yeteneği, açık `workflow` blokları ve adlandırılmış süreç çıktıları.
DSL1, Nextflow 22.03'te kullanımdan kaldırıldı ve 22.12'de kaldırıldı.
Tüm modern Nextflow kodu DSL2 kullanır.

**Sözdizimi ayrıştırıcısı v1 vs v2**, her ikisi de DSL2 koduyla çalışan farklı ayrıştırıcıları ifade eder.
v1 ayrıştırıcısı, orijinal, daha esnek ayrıştırıcıdır.
v2 ayrıştırıcısı daha katıdır ve statik tipleme (tiplendirilmiş girdiler ve çıktılar) ve iş akışı seviyesinde çıktı yönergeleri gibi yeni dil özelliklerini etkinleştirir.
v2 ayrıştırıcısı ayrıca daha iyi hata mesajları sağlar ve daha fazla hatayı çalışma zamanında değil ayrıştırma zamanında yakalar.
v2 ayrıştırıcısı, Nextflow 26.04'te varsayılan hale gelecektir.

Özetle: DSL2 yazdığınız dildir; sözdizimi ayrıştırıcısı sürümü, bu dilin ne kadar katı yorumlandığını ve hangi gelişmiş özelliklerin mevcut olduğunu belirler.

### Nextflow sürümünü kontrol etme ve ayarlama

Sisteminizde hangi Nextflow sürümünün kurulu olduğunu `nextflow --version` komutunu kullanarak kontrol edebilirsiniz.

Nextflow sürümünüzü nasıl güncelleyeceğiniz hakkında daha fazla bilgi için lütfen [Nextflow'u Güncelleme](https://www.nextflow.io/docs/latest/updating-nextflow.html) referans belgelerine bakın.

### v2 sözdizimi ayrıştırıcısını etkinleştirme

Mevcut oturumunuz için v2 sözdizimi ayrıştırıcısını **etkinleştirmek** için terminalinizde aşağıdaki komutu çalıştırın:

```bash
export NXF_SYNTAX_PARSER=v2
```

Bunu kalıcı hale getirmek için (v2'nin Nextflow 26.04'te varsayılan hale gelmesine kadar), export komutunu kabuk profilinize (`~/.bashrc`, `~/.zshrc`, vb.) ekleyin:

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

`NXF_SYNTAX_PARSER=v2` ortam değişkeninin geçici bir gereklilik olduğunu unutmayın.
Nextflow 26.04'ten itibaren, v2 ayrıştırıcısı varsayılan hale gelecek ve bu ayara artık ihtiyaç duyulmayacaktır.

### v2 sözdizimi ayrıştırıcısını devre dışı bırakma

Mevcut oturumunuz için v2 sözdizimi ayrıştırıcısını **devre dışı bırakmak** için terminalinizde aşağıdaki komutu çalıştırın:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- 26.04'ten sonraki sürümlerde devre dışı bırakmak mümkün olacak mı? -->

### Mevcut kodu taşıma

Mevcut kodun Nextflow'un daha yeni sürümlerine uygun hale getirilmesi için taşıma konusunda rehberlik için lütfen referans belgelerindeki [Taşıma Notları](https://www.nextflow.io/docs/latest/migrations/index.html)'na bakın.

Bu iki makale, en son sürüme taşıma için özellikle yararlıdır:

- [İş akışı çıktılarına taşıma](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Statik tiplere taşıma](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Bu özelliklerin her ikisi de eğitim materyallerinin 3.0 sürümünden başlayarak başlangıç eğitiminin bir parçası olarak ele alınmaktadır.

Taşımayı planladığınız Nextflow kodunun nesline bağlı olarak, çoğunu `nextflow lint -format` komutunu kullanarak Nextflow linter'ı ile yapabilirsiniz.
Daha fazla ayrıntı için [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) CLI referansına bakın.

Bunun yararlı olacağını umuyoruz.
Yardıma ihtiyacınız varsa, Slack veya forum üzerinden ulaşın.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
