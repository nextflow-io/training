# Başlarken

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=gZxlXgkVxuLEzOsC" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube kanalında [tüm oynatma listesine](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) bakın.

:green_book: Video transkripti [burada](./transcripts/00_orientation.md) mevcuttur.
///

!!! tip "İpucu"

    YouTube videolarının bazı süper güçleri var!

    - :fontawesome-solid-closed-captioning: Yüksek kaliteli (manuel olarak düzenlenmiş) altyazılar. :material-subtitles: simgesi ile açabilirsiniz
    - :material-bookmark: Sayfa başlıklarına karşılık gelen zaman çizelgesinde video bölümleri.

## Eğitim ortamı başlatma

GitHub Codespaces'te sunduğumuz önceden oluşturulmuş ortamı kullanmak için aşağıdaki "Open in GitHub Codespaces" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../envsetup/index.md) bölümüne bakın.

Ortam yüklenirken okumaya devam edebilmeniz için eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (ekipmanınıza bağlı olarak sağ tıklama, ctrl-tıklama veya cmd-tıklama kullanın).
Kurs boyunca çalışmak için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, eğitim kursu boyunca çalışmak için gerekli tüm yazılım, kod ve verileri içerir, böylece kendiniz bir şey yüklemeniz gerekmez.

Codespace, bir dosya sistemi gezgini, kod düzenleyici ve terminal kabuğu içeren VSCode arayüzü ile kurulmuştur.
Kurs boyunca verilen tüm talimatlar (örn. 'dosyayı açın', 'kodu düzenleyin' veya 'bu komutu çalıştırın') aksi belirtilmedikçe VScode arayüzünün bu üç bölümüne atıfta bulunur.

Bu kursu kendiniz çalışıyorsanız, daha fazla ayrıntı için lütfen [ortam temelleri](../envsetup/01_setup.md) ile tanışın.

### Sürüm gereksinimleri

Bu eğitim, **v2 syntax parser ETKİNLEŞTİRİLMİŞ** Nextflow 25.10.2 veya üstü için tasarlanmıştır.
Yerel veya özel bir ortam kullanıyorsanız, lütfen [burada](../info/nxf_versions.md) belgelenen doğru ayarları kullandığınızdan emin olun.

## Çalışmaya hazırlanın

Codespace'iniz çalıştıktan sonra, eğitime dalmadan önce yapmanız gereken iki şey var: bu belirli kurs için çalışma dizininizi ayarlayın ve sağlanan materyallere bir göz atın.

### Çalışma dizinini ayarlayın

Varsayılan olarak, codespace tüm eğitim kurslarının kökünde çalışma dizini ayarlanmış olarak açılır, ancak bu kurs için `hello-nextflow/` dizininde çalışacağız.

Terminalde bu komutu çalıştırarak şimdi dizini değiştirin:

```bash
cd hello-nextflow/
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz, böylece dosya gezgini kenar çubuğunda yalnızca ilgili dosyalar görünür:

```bash
code .
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örn. codespace'iniz uykuya geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak her zaman tam yolu kullanarak geri dönebilirsiniz:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Şimdi içeriğe bir göz atalım.

### Sağlanan materyalleri keşfedin

Bu dizinin içeriğini eğitim çalışma alanının sol tarafındaki dosya gezginini kullanarak keşfedebilirsiniz.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.

Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde temsil etmek için `tree` çıktısını kullanıyoruz, bazen netlik için küçük değişikliklerle.

Burada ikinci seviyeye kadar içindekiler tablosu oluşturuyoruz:

```bash
tree . -L 2
```

??? abstract "Dizin içerikleri"

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
Bu tür daraltılabilir bölümleri, beklenen komut çıktısını özlü bir şekilde dahil etmek için kullanıyoruz.

- **`.nf` dosyaları**, kursun hangi bölümünde kullanıldıklarına göre adlandırılmış iş akışı betikleridir.

- **`nextflow.config` dosyası**, minimal ortam özelliklerini ayarlayan bir yapılandırma dosyasıdır.
  Şimdilik görmezden gelebilirsiniz.

- **`data/` altındaki `greetings.csv` dosyası**, kursun çoğunda kullanacağımız girdi verilerini içerir. Bölüm 2'de (Channels) ilk kez tanıtıldığında açıklanmaktadır.

- **`test-params.*` dosyaları**, Bölüm 6'da (Configuration) kullanacağımız yapılandırma dosyalarıdır. Şimdilik görmezden gelebilirsiniz.

- **`solutions` dizini**, kursun her adımından elde edilen tamamlanmış iş akışı betiklerini içerir.
  Bunlar, çalışmanızı kontrol etmek ve herhangi bir sorunu gidermek için referans olarak kullanılmak üzere tasarlanmıştır.

## Hazırlık kontrol listesi

Başlamaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Ortamım hazır ve çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım

Tüm kutuları işaretleyebiliyorsanız, hazırsınız.

**[Bölüm 1: Hello World](./01_hello_world.md)'e devam etmek için, bu sayfanın sağ alt köşesindeki oka tıklayın.**
