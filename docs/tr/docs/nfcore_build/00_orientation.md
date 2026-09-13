# Başlarken

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eğitim ortamını başlatma

GitHub Codespaces üzerinde sağladığımız hazır ortamı kullanmak için aşağıdaki "GitHub Codespaces'te Aç" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../envsetup/index.md) sayfasına bakın.

Ortam yüklenirken okumaya devam edebilmeniz için eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (ekipmanınıza bağlı olarak sağ tıklama, ctrl+tıklama veya cmd+tıklama kullanın).
Kursu takip edebilmek için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, eğitim kursu boyunca çalışmak için gerekli tüm yazılım, kod ve veriyi içerir; dolayısıyla kendiniz herhangi bir şey yüklemenize gerek yoktur.

Codespace, bir dosya sistemi gezgini, bir kod düzenleyici ve bir terminal kabuğu içeren bir VSCode arayüzü ile kurulmuştur.
Kurs boyunca verilen tüm talimatlar (örn. 'dosyayı açın', 'kodu düzenleyin' veya 'bu komutu çalıştırın') aksi belirtilmedikçe VSCode arayüzünün bu üç bölümünü ifade eder.

Bu kursu kendi başınıza takip ediyorsanız, lütfen daha fazla ayrıntı için [ortam temellerine](../envsetup/01_setup.md) göz atın.

### Sürüm gereksinimleri

Bu eğitim, varsayılan olarak Nextflow 26.04 sürümünden itibaren etkin olan **v2 sözdizimi ayrıştırıcısıyla** birlikte **Nextflow 25.10.2 veya sonrası** ile çalışmaktadır.
Eğitim ortamımızda herhangi bir şey yapmanıza gerek yoktur: v2 ayrıştırıcısıyla Nextflow 26.04.4 çalıştırılmaktadır. Yerel veya özel bir ortam kullanıyorsanız [sürüm notlarına](../info/nxf_versions.md) bakın.

Eğitim ayrıca **nf-core tools 4.0.2** gerektirir.
nf-core araçlarının farklı bir sürümünü kullanırsanız, takip etmekte zorluk yaşayabilirsiniz.

Ortamınızda hangi sürümün yüklü olduğunu `nf-core --version` komutunu kullanarak kontrol edebilirsiniz.

!!! warning "v2 ayrıştırıcı uyumluluğu"

    Birçok nf-core pipeline'ı henüz v2 sözdizimi ayrıştırıcısını desteklememektedir.
    Bu kursta kullanılanlar dışında bir nf-core pipeline'ı çalıştırır ve hatalarla karşılaşırsanız, `export NXF_SYNTAX_PARSER=v1` ayarını yaparak v1 ayrıştırıcısına geçmeniz gerekebilir.
    Ayrıntılar için [sürüm notlarına](../info/nxf_versions.md) bakın.

## Çalışmaya hazır olun

Codespace'iniz çalışmaya başladığında, eğitime dalmadan önce yapmanız gereken iki şey vardır: bu özel kurs için çalışma dizininizi ayarlamak ve sağlanan materyallere göz atmak.

### Çalışma dizinini ayarlama

Varsayılan olarak, codespace tüm eğitim kurslarının kök dizininde ayarlanmış çalışma dizini ile açılır; ancak bu kurs için `nfcore-build/` dizininde çalışacağız.

Terminalde şu komutu çalıştırarak şimdi dizini değiştirin:

```bash
cd nfcore-build/
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örn. codespace'iniz uykuya geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak, geri dönmek için her zaman tam yolu kullanabilirsiniz:

    ```bash
    cd /workspaces/training/nfcore-build
    ```

Şimdi bu dizinin içeriğini keşfedelim.

### Sağlanan materyalleri keşfetme

Bu dizinin içeriğini, eğitim çalışma alanının sol tarafındaki dosya gezginini kullanarak keşfedebilirsiniz.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.

Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde göstermek için `tree` çıktısını kullanıyoruz; bazen netlik için küçük değişikliklerle.

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
        ├── core-hello-part1
        ├── core-hello-part2
        ├── core-hello-part3
        ├── core-hello-part4
        └── core-hello-start
    ```

Bölümü genişletmek ve içeriğini görüntülemek için renkli kutuya tıklayın.
Beklenen komut çıktısını özlü bir şekilde dahil etmek için bunun gibi daraltılabilir bölümler kullanıyoruz.

- **`greetings.csv` dosyası**, test amaçlı kullandığımız bazı minimal sütunsal verileri içeren bir CSV'dir.

- **`original-hello` dizini**, Hello Nextflow eğitim serisinin tamamı boyunca çalışarak üretilen kaynak kodun bir kopyasını içerir (Docker etkinleştirilmiş olarak).

- **`solutions` dizini**, kursun her adımından kaynaklanan tamamlanmış iş akışı betiklerini içerir.
  Bunlar, çalışmanızı kontrol etmek ve herhangi bir sorunu gidermek için referans olarak kullanılmak üzere tasarlanmıştır.

### Pipeline referansını ayarlama

Bu kursun 1. Bölümü, referans olarak `nf-core/demo` pipeline'ının kaynak kodunu kullanır.
Pipeline'ı indirin ve kolayca göz atabilmek için bir kısayol oluşturun:

```bash
nextflow pull nf-core/demo
ln -s $NXF_HOME/assets pipelines
```

[nf-core'u Kullanma](../nfcore_use/index.md) kursunu daha önce tamamladıysanız, `nf-core/demo` yerel olarak önbellekte mevcuttur; ancak yukarıdaki sembolik bağlantının bu dizinde yeniden oluşturulması gerekir.

## Hazırlık kontrol listesi

Başlamaya hazır olduğunuzu mu düşünüyorsunuz?

- [ ] Bu kursun hedefini ve ön koşullarını anlıyorum
- [ ] Ortamım çalışır durumda
- [ ] nf-core tools 4.0.2 kullandığımı doğruladım (`nf-core --version` ile kontrol edin)
- [ ] Çalışma dizinini uygun şekilde ayarladım

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

**Bölüm 1'e devam etmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
