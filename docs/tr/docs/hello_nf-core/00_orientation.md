# Başlarken

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eğitim ortamını başlatma

GitHub Codespaces üzerinde sağladığımız hazır ortamı kullanmak için aşağıdaki "GitHub Codespaces'te Aç" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../envsetup/index.md) sayfasına bakın.

Ortam yüklenirken okumaya devam edebilmeniz için eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (ekipmanınıza bağlı olarak sağ tıklama, ctrl+tıklama veya cmd+tıklama kullanın).
Kursu takip edebilmek için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, eğitim kursu boyunca çalışmak için gerekli tüm yazılım, kod ve veriyi içerir, dolayısıyla kendiniz herhangi bir şey yüklemenize gerek yoktur.

Codespace, bir dosya sistemi gezgini, bir kod düzenleyici ve bir terminal kabuğu içeren bir VSCode arayüzü ile kurulmuştur.
Kurs boyunca verilen tüm talimatlar (örn. 'dosyayı açın', 'kodu düzenleyin' veya 'bu komutu çalıştırın') aksi belirtilmedikçe VSCode arayüzünün bu üç bölümünü ifade eder.

Bu kursu kendi başınıza takip ediyorsanız, lütfen daha fazla ayrıntı için [ortam temellerine](../envsetup/01_setup.md) göz atın.

### Sürüm gereksinimleri

Bu eğitim, **v2 sözdizimi ayrıştırıcısı DEVRE DIŞI** olan **Nextflow 25.10.2** veya sonrası için tasarlanmıştır.

#### Eğitim ortamımızı kullanıyorsanız:

Daha ileriye gitmeden önce aşağıdaki komutu çalıştırmanız GEREKMEKTEDİR:

```bash
export NXF_SYNTAX_PARSER=v1
```

#### Yerel veya özel bir ortam kullanıyorsanız:

Lütfen [burada](../info/nxf_versions.md) belgelenen doğru ayarları kullandığınızdan emin olun.

Eğitim ayrıca **nf-core tools 3.4.1** gerektirir.
nf-core araçlarının farklı bir sürümünü kullanırsanız, takip etmekte zorluk yaşayabilirsiniz.

Ortamınızda hangi sürümün yüklü olduğunu `nf-core --version` komutunu kullanarak kontrol edebilirsiniz.

## Çalışmaya hazır olun

Codespace'iniz çalışmaya başladığında, eğitime dalmadan önce yapmanız gereken iki şey vardır: bu özel kurs için çalışma dizininizi ayarlamak ve sağlanan materyallere göz atmak.

### Çalışma dizinini ayarlama

Varsayılan olarak, codespace tüm eğitim kurslarının kök dizininde ayarlanmış çalışma dizini ile açılır, ancak bu kurs için `hello-nf-core/` dizininde çalışacağız.

Terminalde şu komutu çalıştırarak şimdi dizini değiştirin:

```bash
cd hello-nf-core/
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örn. codespace'iniz uykuya geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak, geri dönmek için her zaman tam yolu kullanabilirsiniz:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

Şimdi bu dizinin içeriğine göz atalım.

### Sağlanan materyalleri keşfetme

Bu dizinin içeriğini, eğitim çalışma alanının sol tarafındaki dosya gezginini kullanarak keşfedebilirsiniz.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.

Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde göstermek için `tree` çıktısını kullanıyoruz, bazen netlik için küçük değişikliklerle.

Burada ikinci seviyeye kadar bir içindekiler tablosu oluşturuyoruz:

```bash
tree . -L 2
```

??? abstract "Dizin içeriği"

    ```console
    .
    ├── greetings.csv
    ├── original-hello
    │   ├── hello.nf
    │   ├── modules
    │   └── nextflow.config
    └── solutions
        ├── composable-hello
        ├── core-hello-part2
        ├── core-hello-part3
        ├── core-hello-part4
        ├── core-hello-part5
        └── core-hello-start
    ```

Bölümü genişletmek ve içeriğini görüntülemek için renkli kutuya tıklayın.
Beklenen komut çıktısını özlü bir şekilde dahil etmek için bunun gibi daraltılabilir bölümler kullanıyoruz.

- **`greetings.csv` dosyası**, test amaçlı kullandığımız bazı minimal sütunsal verileri içeren bir CSV'dir.

- **`original-hello` dizini**, Hello Nextflow eğitim serisinin tamamı boyunca çalışarak üretilen kaynak kodun bir kopyasını içerir (Docker etkinleştirilmiş olarak).

- **`solutions` dizini**, kursun her adımından kaynaklanan tamamlanmış workflow betiklerini içerir.
  Bunlar, çalışmanızı kontrol etmek ve herhangi bir sorunu gidermek için referans olarak kullanılmak üzere tasarlanmıştır.

## Hazırlık kontrol listesi

Başlamaya hazır olduğunuzu mu düşünüyorsunuz?

- [ ] Bu kursun hedefini ve ön koşullarını anlıyorum
- [ ] Ortamım çalışır durumda
- [ ] Sözdizimi ayrıştırıcısının **v1** olarak ayarlandığından emin oldum
- [ ] Çalışma dizinini uygun şekilde ayarladım

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

**Bölüm 1'e devam etmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
