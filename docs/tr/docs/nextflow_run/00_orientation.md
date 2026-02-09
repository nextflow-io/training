# Başlarken

## Eğitim ortamını başlatın

GitHub Codespaces üzerinde sağladığımız önceden hazırlanmış ortamı kullanmak için aşağıdaki "GitHub Codespaces'te Aç" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../envsetup/index.md) sayfasına bakın.

Eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (ekipmanınıza bağlı olarak sağ tıklama, ctrl-tıklama veya cmd-tıklama kullanın), böylece ortam yüklenirken okumaya devam edebilirsiniz.
Kursu tamamlamak için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, eğitim kursunu tamamlamak için gerekli tüm yazılımları, kodları ve verileri içerir, bu nedenle kendiniz herhangi bir şey yüklemenize gerek yoktur.

Codespace, bir dosya sistemi gezgini, bir kod düzenleyici ve bir terminal kabuğu içeren VSCode arayüzü ile kurulmuştur.
Kurs boyunca verilen tüm talimatlar (örneğin 'dosyayı açın', 'kodu düzenleyin' veya 'bu komutu çalıştırın'), aksi belirtilmedikçe VSCode arayüzünün bu üç bölümüne atıfta bulunur.

Bu kursu kendi başınıza tamamlıyorsanız, lütfen daha fazla ayrıntı için [ortam temelleri](../envsetup/01_setup.md) ile tanışın.

### Sürüm gereksinimleri

Bu eğitim, Nextflow 25.10.2 veya sonrası için **v2 sözdizimi ayrıştırıcısı ETKİNLEŞTİRİLMİŞ olarak** tasarlanmıştır.
Yerel veya özel bir ortam kullanıyorsanız, lütfen [burada](../info/nxf_versions.md) belgelenen doğru ayarları kullandığınızdan emin olun.

## Çalışmaya hazırlanın

Codespace'iniz çalışmaya başladığında, eğitime dalmadan önce yapmanız gereken iki şey var: bu özel kurs için çalışma dizininizi ayarlayın ve sağlanan materyallere göz atın.

### Çalışma dizinini ayarlayın

Varsayılan olarak, codespace tüm eğitim kurslarının kök dizininde çalışma dizini ayarlanmış olarak açılır, ancak bu kurs için `nextflow-run/` dizininde çalışacağız.

Terminalde şu komutu çalıştırarak şimdi dizini değiştirin:

```bash
cd nextflow-run/
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz, böylece dosya gezgini kenar çubuğunda yalnızca ilgili dosyalar görünür:

```bash
code .
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örneğin codespace'iniz uyku moduna geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak, her zaman tam yolu kullanarak geri dönebilirsiniz:

    ```bash
    cd /workspaces/training/nextflow-run
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
    ├── 1-hello.nf
    ├── 2a-inputs.nf
    ├── 2b-multistep.nf
    ├── 2c-modules.nf
    ├── 2d-container.nf
    ├── 3-main.nf
    ├── data
    │   └── greetings.csv
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 3-main.nf
    │   ├── modules
    │   └── nextflow.config
    ├── test-params.json
    └── test-params.yaml
    ```

Bölümü genişletmek ve içeriğini görüntülemek için renkli kutuya tıklayın.
Beklenen komut çıktısını ve ayrıca dizin ve dosya içeriklerini özlü bir şekilde görüntülemek için bunun gibi daraltılabilir bölümler kullanırız.

- **`.nf` dosyaları**, kursun hangi bölümünde kullanıldıklarına göre numaralandırılmış iş akışı betikleridir.

- **`nextflow.config` dosyası**, minimal ortam özelliklerini ayarlayan bir yapılandırma dosyasıdır.
  Şimdilik görmezden gelebilirsiniz.

- **`data/` altındaki `greetings.csv` dosyası**, kursun çoğunda kullanacağımız girdi verilerini içerir. Bölüm 2'de (Pipeline'ları çalıştırma) ilk kez tanıttığımızda açıklanmaktadır.

- **`test-params.*` dosyaları**, Bölüm 3'te (Yapılandırma) kullanacağımız yapılandırma dosyalarıdır. Şimdilik görmezden gelebilirsiniz.

- **`solutions` dizini**, kursu tamamlamanın sonucunda ortaya çıkan iş akışının ve yardımcı dosyalarının (config ve modüller) son halini içerir.
  Çalışmanızı kontrol etmek ve sorunları gidermek için referans olarak kullanılmaları amaçlanmıştır.

## Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Ortamım çalışır durumda
- [ ] Çalışma dizinini uygun şekilde ayarladım

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

**[Bölüm 1: Temel işlemleri çalıştırma](./01_basics.md) bölümüne devam etmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
