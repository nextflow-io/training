# Başlarken

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=gZxlXgkVxuLEzOsC" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube kanalında [tüm oynatma listesini](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) görün.

:green_book: Video metni [burada](./transcripts/00_orientation.md) mevcuttur.
///

!!! tip "İpucu"

    YouTube videolarının süper güçleri var!

    - :fontawesome-solid-closed-captioning: Yüksek kaliteli (manuel olarak düzenlenmiş) altyazılar. :material-subtitles: simgesiyle açabilirsiniz
    - :material-bookmark: Sayfa başlıklarına karşılık gelen zaman çizelgesinde video bölümleri.

## Eğitim ortamını başlatın

GitHub Codespaces üzerinde sağladığımız önceden hazırlanmış ortamı kullanmak için aşağıdaki "Open in GitHub Codespaces" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../envsetup/index.md) sayfasına bakın.

Eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (ekipmanınıza bağlı olarak sağ tıklama, ctrl-tıklama veya cmd-tıklama kullanın), böylece ortam yüklenirken okumaya devam edebilirsiniz.
Kursu tamamlamak için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, eğitim kursunda çalışmak için gerekli tüm yazılımları, kodları ve verileri içerir, bu nedenle kendiniz herhangi bir şey yüklemenize gerek yoktur.

Codespace, bir dosya sistemi gezgini, bir kod düzenleyici ve bir terminal kabuğu içeren bir VSCode arayüzü ile kurulmuştur.
Kurs boyunca verilen tüm talimatlar (örneğin 'dosyayı açın', 'kodu düzenleyin' veya 'bu komutu çalıştırın') aksi belirtilmedikçe VSCode arayüzünün bu üç bölümüne atıfta bulunur.

Bu kursu kendi başınıza çalışıyorsanız, lütfen daha fazla ayrıntı için [ortam temelleri](../envsetup/01_setup.md) ile tanışın.

### Sürüm gereksinimleri

Bu eğitim, Nextflow 25.10.2 veya sonrası için **v2 sözdizimi ayrıştırıcısı ETKİNLEŞTİRİLMİŞ** olarak tasarlanmıştır.
Yerel veya özel bir ortam kullanıyorsanız, lütfen [burada](../info/nxf_versions.md) belgelenen doğru ayarları kullandığınızdan emin olun.

## Çalışmaya hazırlanın

Codespace'iniz çalışmaya başladığında, eğitime dalmadan önce yapmanız gereken iki şey var: bu özel kurs için çalışma dizininizi ayarlayın ve sağlanan materyallere bir göz atın.

### Çalışma dizinini ayarlayın

Varsayılan olarak, codespace tüm eğitim kurslarının kök dizininde ayarlanmış çalışma dizini ile açılır, ancak bu kurs için `hello-nextflow/` dizininde çalışacağız.

Terminalde şu komutu çalıştırarak şimdi dizini değiştirin:

```bash
cd hello-nextflow/
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz, böylece dosya gezgini kenar çubuğunda yalnızca ilgili dosyalar görünür:

```bash
code .
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örneğin codespace'iniz uyku moduna geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak, her zaman tam yolu kullanarak geri dönebilirsiniz:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Şimdi içeriğe bir göz atalım.

### Sağlanan materyalleri keşfedin

Bu dizinin içeriğini, eğitim çalışma alanının sol tarafındaki dosya gezginini kullanarak keşfedebilirsiniz.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.

Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde temsil etmek için `tree` çıktısını kullanırız, bazen netlik için küçük değişikliklerle.

Burada ikinci seviyeye kadar bir içindekiler tablosu oluşturuyoruz:

```bash
tree . -L 2
```

??? abstract "Dizin içeriği"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── test-params.yaml
    ```

Bölümü genişletmek ve içeriğini görüntülemek için renkli kutuya tıklayın.
Beklenen komut çıktısını özlü bir şekilde dahil etmek için bunun gibi daraltılabilir bölümler kullanırız.

- **`.nf` dosyaları**, kursun hangi bölümünde kullanıldıklarına göre adlandırılmış iş akışı betikleridir.

- **`nextflow.config` dosyası**, minimal ortam özelliklerini ayarlayan bir yapılandırma dosyasıdır.
  Şimdilik görmezden gelebilirsiniz.

- **`data/` altındaki `greetings.csv` dosyası**, kursun çoğunda kullanacağımız girdi verilerini içerir. İlk kez tanıttığımızda Bölüm 2'de (Kanallar) açıklanmaktadır.

- **`test-params.*` dosyaları**, Bölüm 6'da (Yapılandırma) kullanacağımız yapılandırma dosyalarıdır. Şimdilik görmezden gelebilirsiniz.

- **`solutions` dizini**, kursun her adımından kaynaklanan tamamlanmış iş akışı betiklerini içerir.
  Çalışmanızı kontrol etmek ve sorunları gidermek için referans olarak kullanılmaları amaçlanmıştır.

## Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Ortamım çalışır durumda
- [ ] Çalışma dizinini uygun şekilde ayarladım

Tüm kutuları işaretleyebiliyorsanız, hazırsınız demektir.

**[Bölüm 1: Hello World](./01_hello_world.md) bölümüne devam etmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
