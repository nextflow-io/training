# Başlangıç

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eğitim ortamını başlatın

GitHub Codespaces'te sağladığımız önceden oluşturulmuş ortamı kullanmak için aşağıdaki "Open in GitHub Codespaces" butonuna tıklayın. Diğer seçenekler için [Ortam seçenekleri](../envsetup/index.md) bölümüne bakın.

Ortam yüklenirken okumaya devam edebilmeniz için eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (cihazınıza bağlı olarak sağ tıklama, ctrl-tıklama veya cmd-tıklama kullanın).
Kurs boyunca çalışmak için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, eğitim kursunda çalışmak için gerekli tüm yazılım, kod ve verileri içerir, böylece kendiniz hiçbir şey yüklemeniz gerekmez.

Codespace, dosya sistemi gezgini, kod editörü ve terminal shell içeren bir VSCode arayüzüyle kurulmuştur.
Kurs sırasında verilen tüm talimatlar (örn. 'dosyayı açın', 'kodu düzenleyin' veya 'bu komutu çalıştırın') aksi belirtilmedikçe VScode arayüzünün bu üç bölümüne atıfta bulunur.

Bu kursu kendi başınıza çalışıyorsanız, daha fazla ayrıntı için [ortam temelleri](../envsetup/01_setup.md) ile tanışın.

### Versiyon gereksinimleri

Bu eğitim, v2 sözdizimi ayrıştırıcısı ETKİNLEŞTİRİLMİŞ olarak Nextflow 25.10.2 veya üstü için tasarlanmıştır.
Yerel veya özel bir ortam kullanıyorsanız, [burada](../info/nxf_versions.md) belgelendiği gibi doğru ayarları kullandığınızdan emin olun.

## Çalışmaya hazırlanın

Codespace'iniz çalışmaya başladığında, eğitime dalmadan önce yapmanız gereken iki şey var: bu kursa özel çalışma dizininizi ayarlayın ve sağlanan materyallere göz atın.

### Çalışma dizinini ayarlayın

Varsayılan olarak, codespace tüm eğitim kurslarının kökünde çalışma dizini ayarlanmış şekilde açılır, ancak bu kurs için `nextflow-run/` dizininde çalışacağız.

Şimdi terminalde bu komutu çalıştırarak dizini değiştirin:

```bash
cd nextflow-run/
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz, böylece dosya gezgini kenar çubuğunda yalnızca ilgili dosyalar görünür:

```bash
code .
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örn. codespace'iniz uyursa), Github Codespaces eğitim ortamında çalıştığınızı varsayarak, ona geri dönmek için her zaman tam yolu kullanabilirsiniz:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

Şimdi içeriğine bir göz atalım.

### Sağlanan materyalleri keşfedin

Bu dizinin içeriğini eğitim çalışma alanının sol tarafındaki dosya gezginini kullanarak keşfedebilirsiniz.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.

Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde temsil etmek için `tree` çıktısını kullanıyoruz, bazen netlik için küçük değişikliklerle.

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
Beklenen komut çıktısının yanı sıra dizin ve dosya içeriğini özlü bir şekilde görüntülemek için bu tür daraltılabilir bölümler kullanıyoruz.

- **`.nf` dosyaları**, kursun hangi bölümünde kullanıldıklarına göre numaralandırılmış workflow betikleridir.

- **`nextflow.config` dosyası**, minimum ortam özelliklerini ayarlayan bir yapılandırma dosyasıdır.
  Şimdilik görmezden gelebilirsiniz.

- **`data/` altındaki `greetings.csv` dosyası**, kursun büyük bölümünde kullanacağımız girdi verilerini içerir. İlk kez tanıttığımız Bölüm 2'de (Pipeline'ları çalıştırma) açıklanmıştır.

- **`test-params.*`** dosyaları, Bölüm 3'te (Yapılandırma) kullanacağımız yapılandırma dosyalarıdır. Şimdilik görmezden gelebilirsiniz.

- **`solutions` dizini**, kursu tamamlamanın sonucunda oluşan workflow'un ve yardımcı dosyalarının (yapılandırma ve modüller) son halini içerir.
  Çalışmanızı kontrol etmek ve herhangi bir sorunu gidermek için referans olarak kullanılması amaçlanmıştır.

## Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu mu düşünüyorsunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Ortamım çalışıyor
- [ ] Çalışma dizinimi uygun şekilde ayarladım

Tüm kutuları işaretleyebiliyorsanız, hazırsınız.

**[Bölüm 1: Temel İşlemleri Çalıştırma](./01_basics.md) bölümüne devam etmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
