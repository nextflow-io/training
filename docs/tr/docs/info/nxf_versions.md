---
title: Nextflow sürümleri
description: Nextflow sözdizimi sürümlerinin evrimini anlama ve yönetme
hide:
  - toc
  - footer
---

## Şu anda desteklenen Nextflow sözdizimi sürümü ve gereksinimleri

Eğitim portalının 3.0 sürümünden itibaren, tüm eğitim kurslarımız, kurs indeks sayfasında aksi belirtilmedikçe Nextflow'un 25.10.2 sürümüne dayanmaktadır (sürüm bildirimi içermeyen kullanımdan kaldırılmış veya arşivlenmiş materyaller hariç).

Kurslar artık workflow seviyesinde tipli girdiler ve workflow seviyesinde çıktı yönergeleri kullandığından, V2 sözdizimi ayrıştırıcısının kullanılması gerekmektedir.
[Github Codespaces](../envsetup/01_setup.md) veya [yerel devcontainer'lar](../envsetup/03_devcontainer.md) aracılığıyla sağladığımız ortamı kullanmayı planlıyorsanız, kurs talimatlarında özellikle belirtilmedikçe bir şey yapmanıza gerek yoktur.
Ancak, eğitimleri kendi ortamınızda çalışmayı planlıyorsanız ([Manuel kurulum](../envsetup/02_local.md)), v2 sözdizimi ayrıştırıcısı etkinleştirilmiş olarak Nextflow 25.10.2 veya sonraki bir sürümü kullandığınızdan emin olmanız gerekecektir.

## Eğitim materyallerinin eski sürümleri

Eğitim materyallerimiz Şubat 2025'ten beri sürümlendirilmektedir.

**25.10.2 öncesi** Nextflow sürümleriyle çalışan eğitim materyallerinin eski sürümlerine, her sayfanın üstündeki eğitim materyallerinin numaralı sürümünü gösteren açılır menü öğesi aracılığıyla erişebilirsiniz.
Eğitim materyallerinin daha eski bir sürümünü seçtiğinizde, eğitim ortamına bağlantılar otomatik olarak ortamın ilgili sürümünü belirtecektir.

## Nextflow sözdizimi sürümleri hakkında diğer bilgiler

Nextflow'un bazen karıştırılan iki farklı sürümleme kavramı vardır: **DSL sürümleri** ve **sözdizimi ayrıştırıcı sürümleri**.

**DSL1 ve DSL2**, Nextflow pipeline'ları yazmanın temelden farklı yollarını ifade eder.
DSL1, process'lerin channel'lar aracılığıyla örtük olarak bağlandığı orijinal sözdizimiydi.
Nextflow 20.07'de tanıtılan DSL2, modülerlik özellikleri ekledi: diğer dosyalardan process'leri ve workflow'ları içe aktarma yeteneği, açık `workflow` blokları ve adlandırılmış process çıktıları.
DSL1, Nextflow 22.03'te kullanımdan kaldırıldı ve 22.12'de kaldırıldı.
Tüm modern Nextflow kodu DSL2 kullanır.

**Sözdizimi ayrıştırıcı v1 ve v2**, her ikisi de DSL2 koduyla çalışan farklı ayrıştırıcıları ifade eder.
v1 ayrıştırıcısı, daha toleranslı olan orijinal ayrıştırıcıdır.
v2 ayrıştırıcısı daha katıdır ve statik tipleme (tipli girdiler ve çıktılar) ile workflow seviyesinde çıktı yönergeleri gibi yeni dil özelliklerini etkinleştirir.
v2 ayrıştırıcısı ayrıca daha iyi hata mesajları sağlar ve çalışma zamanı yerine ayrıştırma zamanında daha fazla hata yakalar.
v2 ayrıştırıcısı Nextflow 26.04'te varsayılan olacaktır.

Özetle: DSL2 yazdığınız dildir; sözdizimi ayrıştırıcı sürümü, bu dilin ne kadar katı yorumlandığını ve hangi gelişmiş özelliklerin mevcut olduğunu belirler.

### Nextflow sürümünü kontrol etme ve ayarlama

Sisteminizde hangi Nextflow sürümünün yüklü olduğunu `nextflow --version` komutuyla kontrol edebilirsiniz.

Nextflow sürümünüzü nasıl güncelleyeceğiniz hakkında daha fazla bilgi için lütfen referans belgelerindeki [Nextflow'u Güncelleme](https://www.nextflow.io/docs/latest/updating-nextflow.html) bölümüne bakın.

### v2 sözdizimi ayrıştırıcısını etkinleştirme

Mevcut oturumunuz için v2 sözdizimi ayrıştırıcısını **etkinleştirmek** için terminalinizde aşağıdaki komutu çalıştırın:

```bash
export NXF_SYNTAX_PARSER=v2
```

Bunu kalıcı hale getirmek için (v2'nin Nextflow 26.04'te varsayılan olmasını beklerken), export komutunu shell profilinize (`~/.bashrc`, `~/.zshrc`, vb.) ekleyin:

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

`NXF_SYNTAX_PARSER=v2` ortam değişkeninin geçici bir gereklilik olduğunu unutmayın.
Nextflow 26.04'ten itibaren v2 ayrıştırıcısı varsayılan olacak ve bu ayar artık gerekli olmayacaktır.

### v2 sözdizimi ayrıştırıcısını devre dışı bırakma

Mevcut oturumunuz için v2 sözdizimi ayrıştırıcısını **devre dışı bırakmak** için terminalinizde aşağıdaki komutu çalıştırın:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### Mevcut kodun taşınması

Mevcut kodun daha yeni Nextflow sürümlerine uyacak şekilde taşınmasına ilişkin rehberlik için lütfen referans belgelerindeki [Taşıma Notları](https://www.nextflow.io/docs/latest/migrations/index.html)'na bakın.

Bu iki makale, en son sürüme taşıma için özellikle faydalıdır:

- [Workflow çıktılarına taşıma](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Statik tiplere taşıma](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Bu özelliklerin her ikisi de eğitim materyallerinin 3.0 sürümünden başlayarak başlangıç eğitiminin bir parçası olarak ele alınmaktadır.

Taşımayı planladığınız Nextflow kodunun nesline bağlı olarak, `nextflow lint -format` komutunu kullanarak Nextflow linter'ı ile çoğunu yapabilirsiniz.
Daha fazla ayrıntı için [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) CLI referansına bakın.

Bunun faydalı olacağını umuyoruz.
Yardıma ihtiyacınız olursa Slack'te veya forumda bize ulaşın.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
