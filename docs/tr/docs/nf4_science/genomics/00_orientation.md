# Yönlendirme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Eğitim ortamı, bu eğitim kursunda çalışmak için gerekli tüm yazılımı, kodu ve veriyi içerir, bu nedenle kendiniz herhangi bir şey yüklemenize gerek yoktur.
Ancak, giriş yapmak için (ücretsiz) bir hesaba ihtiyacınız var ve arayüzü tanımak için birkaç dakikanızı ayırmalısınız.

Henüz yapmadıysanız, lütfen daha ileri gitmeden önce [bu bağlantıyı](../../../envsetup/) takip edin.

## Sağlanan Materyaller

Bu eğitim kursu boyunca, `nf4-science/genomics/` dizininde çalışacağız; eğitim çalışma alanını açtığınızda bu dizine geçmeniz gerekir.
Bu dizin, ihtiyaç duyacağınız tüm kod dosyalarını, test verilerini ve yardımcı dosyaları içerir.

Bu dizinin içeriğini keşfetmekten çekinmeyin; bunu yapmanın en kolay yolu, VSCode arayüzündeki eğitim çalışma alanının sol tarafındaki dosya gezginini kullanmaktır.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.
Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde temsil etmek için `tree` çıktısını kullanıyoruz, bazen netlik için küçük değişikliklerle.

Burada ikinci seviyeye kadar bir içindekiler tablosu oluşturuyoruz:

```bash
tree . -L 2
```

Bunu `nf4-science/genomics` içinde çalıştırırsanız, aşağıdaki çıktıyı görmelisiniz:

```console title="Dizin içeriği"

.
├── data
│   ├── bam
│   ├── ref
│   ├── sample_bams.txt
│   └── samplesheet.csv
├── genomics-1.nf
├── genomics-2.nf
├── genomics-3.nf
├── genomics-4.nf
├── nextflow.config
└── solutions
    ├── modules
    ├── nf-test.config
    └── tests

6 directories, 8 files

```

!!!note "Not"

    Bu çok fazla gibi görünse de endişelenmeyin; kursun her adımında ilgili parçaları ele alacağız.
    Bu sadece size bir genel bakış sunmak için.

**Başlamak için bilmeniz gerekenler burada özetlenmiştir:**

- **`.nf` dosyaları**, kursun hangi bölümünde kullanıldıklarına göre adlandırılmış iş akışı betikleridir.

- **`nextflow.config` dosyası**, minimal ortam özelliklerini ayarlayan bir yapılandırma dosyasıdır.
  Şimdilik bunu göz ardı edebilirsiniz.

- **`data` dizini**, girdi verilerini ve ilgili kaynakları içerir; kursta daha sonra açıklanacaktır.

- **`solutions` dizini**, kursun 3. ve 4. Bölümlerinden kaynaklanan modül dosyalarını ve test yapılandırmalarını içerir.
  Çalışmanızı kontrol etmek ve sorunları gidermek için referans olarak kullanılması amaçlanmıştır.

!!!tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız, her zaman buraya dönmek için bu komutu çalıştırabilirsiniz:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Şimdi, kursa başlamak için bu sayfanın sağ alt köşesindeki oka tıklayın.
