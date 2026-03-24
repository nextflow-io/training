# Başlarken

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eğitim ortamını başlatın

GitHub Codespaces üzerinde sağladığımız önceden hazırlanmış ortamı kullanmak için aşağıdaki "GitHub Codespaces'te Aç" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../../envsetup/index.md) sayfasına bakın.

Eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (ekipmanınıza bağlı olarak sağ tıklama, ctrl-tıklama veya cmd-tıklama kullanın), böylece ortam yüklenirken okumaya devam edebilirsiniz.
Kurs boyunca çalışmak için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, eğitim kursu boyunca çalışmak için gerekli tüm yazılımı, kodu ve veriyi içerir; bu nedenle kendiniz herhangi bir şey yüklemenize gerek yoktur.

Codespace, bir dosya sistemi gezgini, bir kod düzenleyici ve bir terminal kabuğu içeren bir VSCode arayüzü ile kurulmuştur.
Kurs boyunca verilen tüm talimatlar (örneğin 'dosyayı açın', 'kodu düzenleyin' veya 'bu komutu çalıştırın'), aksi belirtilmedikçe VSCode arayüzünün bu üç bölümüne atıfta bulunur.

Bu kursu kendi başınıza çalışıyorsanız, lütfen daha fazla ayrıntı için [ortam temellerini](../../envsetup/01_setup.md) inceleyin.

### Sürüm gereksinimleri

Bu eğitim, Nextflow 25.10.2 veya sonrası için **v2 sözdizimi ayrıştırıcısı ETKİNLEŞTİRİLMİŞ** olarak tasarlanmıştır.
Yerel veya özel bir ortam kullanıyorsanız, lütfen [burada](../../info/nxf_versions.md) belgelenen doğru ayarları kullandığınızdan emin olun.

## Çalışmaya hazırlanın

Codespace'iniz çalışmaya başladığında, eğitime dalmadan önce yapmanız gereken iki şey var: bu özel kurs için çalışma dizininizi ayarlayın ve sağlanan materyallere bir göz atın.

### Çalışma dizinini ayarlayın

Varsayılan olarak, codespace tüm eğitim kurslarının kök dizininde ayarlanmış çalışma dizini ile açılır; ancak bu kurs için `nf4-science/{DOMAIN_DIR}/` dizininde çalışacağız.

Terminalde şu komutu çalıştırarak şimdi dizini değiştirin:

```bash
cd nf4-science/{DOMAIN_DIR}/
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz, böylece dosya gezgini kenar çubuğunda yalnızca ilgili dosyalar görünür:

```bash
code .
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örneğin codespace'iniz uyku moduna geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak, ona dönmek için her zaman tam yolu kullanabilirsiniz:

    ```bash
    cd /workspaces/training/nf4-science/{DOMAIN_DIR}
    ```

Şimdi içeriğe bir göz atalım.

### Sağlanan materyalleri keşfedin

Bu dizinin içeriğini, eğitim çalışma alanının sol tarafındaki dosya gezginini kullanarak keşfedebilirsiniz.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.

Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde temsil etmek için `tree` çıktısını kullanırız; bazen netlik için küçük değişikliklerle.

Burada ikinci seviyeye kadar bir içindekiler tablosu oluşturuyoruz:

```bash
tree . -L 2
```

??? abstract "Dizin içeriği"

    ```console
    .
    ├── data
    │   ├── {DATA_SUBDIRS}
    │   └── samplesheet.csv
    ├── {DOMAIN_DIR}.nf
    ├── modules
    │   ├── {TOOL_A_MODULE}.nf
    │   └── {TOOL_B_MODULE}.nf
    ├── nextflow.config
    └── solutions
        ├── modules
        ├── part2
        └── part3

    N directories, N files
    ```

Bölümü genişletmek ve içeriğini görüntülemek için renkli kutuya tıklayın.
Beklenen komut çıktısını ve ayrıca dizin ve dosya içeriklerini özlü bir şekilde görüntülemek için bunun gibi daraltılabilir bölümler kullanırız.

- **`{DOMAIN_DIR}.nf` dosyası**, kurs boyunca oluşturacağınız bir iş akışı betiğidir.

- **`modules` dizini**, kurs boyunca dolduracağınız iskelet modül dosyalarını içerir.

- **`nextflow.config` dosyası**, minimum ortam özelliklerini ayarlayan bir yapılandırma dosyasıdır.
  Şimdilik göz ardı edebilirsiniz.

- **`data` dizini**, girdi verilerini ve ilgili kaynakları içerir; kursta daha sonra açıklanacaktır.

- **`solutions` dizini**, tamamlanmış modül dosyalarını ve bir sonraki bölüm için başlangıç noktası olarak kullanılabilecek bölüme özel çözümleri içerir.
  Çalışmanızı kontrol etmek ve sorunları gidermek için referans olarak kullanılmaları amaçlanmıştır.

## Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu mu düşünüyorsunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Ortamım çalışır durumda
- [ ] Çalışma dizinini uygun şekilde ayarladım

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

**[Bölüm 1: Yönteme genel bakış](./01_method.md) sayfasına devam etmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
