# Yönlendirme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu yönlendirme, "Open in GitHub Codespaces" düğmesine tıklayarak eğitim ortamını zaten açtığınızı varsayar.
Henüz açmadıysanız, lütfen şimdi açın; ideal olarak ikinci bir tarayıcı penceresi veya sekmesinde açarak bu talimatlara geri dönebilirsiniz.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "Makine boyutu gereksinimi"

    Bu eğitim kursu için Codespace oluştururken mutlaka **8 çekirdekli makine** seçtiğinizden emin olun. Biyogörüntüleme iş akışları ek hesaplama kaynakları gerektirir.

## GitHub Codespaces

GitHub Codespaces ortamı, bu eğitim kursunda çalışmak için gerekli tüm yazılımı, kodu ve veriyi içerir, bu nedenle kendiniz herhangi bir şey yüklemenize gerek yoktur.
Ancak, giriş yapmak için (ücretsiz) bir GitHub hesabına ihtiyacınız var ve arayüze aşina değilseniz, [GitHub Codespaces Yönlendirmesi](../../envsetup/index.md) mini kursunu tamamlayarak birkaç dakikanızı ayırıp arayüze alışmalısınız.

## Docker imajlarını önceden indirin

Codespace'inizi açtıktan sonra, bu eğitim kursu için ihtiyaç duyacağımız tüm Docker imajlarını önceden indirelim.
Bu, daha sonra zaman kazandıracak ve iş akışlarının sorunsuz çalışmasını sağlayacaktır.

Yeni bir terminal sekmesi açın ve aşağıdaki komutu çalıştırın:

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

Bu komut, gerekli tüm Docker imajlarını arka planda indirecektir.
Bu çalışırken yönlendirmenin geri kalanıyla devam edebilirsiniz.

!!!tip "İpucu"

    `-stub` bayrağı, pipeline'ın gerçek veri işlemeden hızlıca çalışmasını sağlar, bu da imajları indirmek için mükemmeldir. Terminal sekmesinde ilerlemeyi izleyebilirsiniz.

## Çalışma dizini

Bu eğitim kursu boyunca, `nf4-science/imaging/` dizininde çalışacağız.

Şimdi terminalde şu komutu çalıştırarak dizini değiştirin:

```bash
cd nf4-science/imaging/
```

!!!tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız, GitHub Codespaces eğitim ortamında çalıştığınızı varsayarak, her zaman tam yolu kullanarak geri dönebilirsiniz:

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**Şimdi, kursa başlamak için bu sayfanın sağ alt köşesindeki oka tıklayın.**
