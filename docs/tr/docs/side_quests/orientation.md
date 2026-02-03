# Yönlendirme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

GitHub Codespaces ortamı, bu eğitim kursunu tamamlamak için gerekli tüm yazılımları, kodları ve verileri içerir, bu nedenle kendiniz herhangi bir şey yüklemenize gerek yoktur.
Ancak, oturum açmak için (ücretsiz) bir hesaba ihtiyacınız var ve arayüze aşina olmak için birkaç dakikanızı ayırmalısınız.

Henüz yapmadıysanız, lütfen daha ileri gitmeden önce [bu bağlantıyı](../../envsetup/) takip edin.

## Sağlanan materyaller

Bu eğitim kursu boyunca, `side-quests/` dizininde çalışacağız.
Bu dizin, ihtiyacınız olacak tüm kod dosyalarını, test verilerini ve ek dosyaları içerir.

Bu dizinin içeriğini keşfetmekten çekinmeyin; bunu yapmanın en kolay yolu, GitHub Codespaces çalışma alanının sol tarafındaki dosya gezginini kullanmaktır.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.
Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde temsil etmek için `tree` çıktısını kullanıyoruz, bazen netlik için küçük değişikliklerle.

Burada ikinci seviyeye kadar bir içindekiler tablosu oluşturuyoruz:

```bash
tree . -L 2
```

Bunu `side-quests` içinde çalıştırırsanız, aşağıdaki çıktıyı görmelisiniz:

```console title="Dizin içeriği"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Başlamak için bilmeniz gerekenler şunlardır:**

- **Her dizin ayrı bir yan göreve karşılık gelir.**
  İçerikleri ilgili yan görevin sayfasında ayrıntılı olarak açıklanmıştır.

- **`solutions` dizini**, her yan görevin çeşitli adımlarını tamamlamaktan kaynaklanan tamamlanmış workflow ve/veya modül scriptlerini içerir.
  Çalışmanızı kontrol etmek ve sorunları gidermek için referans olarak kullanılmaları amaçlanmıştır.

!!!tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız, buraya geri dönmek için her zaman bu komutu çalıştırabilirsiniz:

    ```bash
    cd /workspaces/training/side-quests
    ```

Şimdi, kursa başlamak için bu sayfanın sağ alt köşesindeki oka tıklayın.
