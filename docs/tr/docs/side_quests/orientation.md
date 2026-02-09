# Yönlendirme

GitHub Codespaces ortamı, bu eğitim kursunu tamamlamak için gerekli tüm yazılımı, kodu ve veriyi içerir, bu nedenle kendiniz herhangi bir şey yüklemenize gerek yoktur.
Ancak, giriş yapmak için (ücretsiz) bir hesaba ihtiyacınız vardır ve arayüze aşina olmak için birkaç dakikanızı ayırmalısınız.

Henüz yapmadıysanız, lütfen daha ileri gitmeden önce [bu bağlantıyı](../../envsetup/) takip edin.

## Sağlanan materyaller

Bu eğitim kursu boyunca, `side-quests/` dizininde çalışacağız.
Bu dizin, ihtiyaç duyacağınız tüm kod dosyalarını, test verilerini ve yardımcı dosyaları içerir.

Bu dizinin içeriğini keşfetmekten çekinmeyin; bunu yapmanın en kolay yolu GitHub Codespaces çalışma alanının sol tarafındaki dosya gezginini kullanmaktır.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.
Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde temsil etmek için `tree` çıktısını kullanıyoruz, bazen netlik için küçük değişikliklerle.

Burada ikinci seviyeye kadar bir içindekiler tablosu oluşturuyoruz:

```bash
tree . -L 2
```

Bunu `side-quests` içinde çalıştırırsanız, aşağıdaki çıktıyı görmelisiniz:

```console title="Directory contents"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Başlamak için bilmeniz gerekenler şunlardır:**

- **Her dizin bireysel bir yan göreve karşılık gelir.**
  İçerikleri ilgili yan görevin sayfasında detaylandırılmıştır.

- **`solutions` dizini**, her yan görevin çeşitli adımlarını tamamlayarak elde edilen tamamlanmış iş akışı ve/veya modül betiklerini içerir.
  Çalışmanızı kontrol etmek ve sorunları gidermek için referans olarak kullanılmaları amaçlanmıştır.

!!!tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız, buraya geri dönmek için her zaman şu komutu çalıştırabilirsiniz:

    ```bash
    cd /workspaces/training/side-quests
    ```

Şimdi, kursa başlamak için bu sayfanın sağ alt köşesindeki oka tıklayın.
