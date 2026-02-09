# Başlarken

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eğitim ortamını başlatın

GitHub Codespaces üzerinde sağladığımız önceden hazırlanmış ortamı kullanmak için aşağıdaki "GitHub Codespaces'te Aç" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../../envsetup/index.md) sayfasına bakın.

Eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (ekipmanınıza bağlı olarak sağ tıklama, ctrl-tıklama veya cmd-tıklama kullanın), böylece ortam yüklenirken okumaya devam edebilirsiniz.
Kurs boyunca çalışmak için bu talimatları paralel olarak açık tutmanız gerekecek.

[![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, eğitim kursunda çalışmak için gerekli tüm yazılımı, kodu ve veriyi içerir, bu nedenle kendiniz herhangi bir şey yüklemenize gerek yoktur.

Codespace bir VSCode arayüzü ile kurulmuştur; bu arayüz bir dosya sistemi gezgini, bir kod düzenleyici ve bir terminal kabuğu içerir.
Kurs boyunca verilen tüm talimatlar (örn. 'dosyayı açın', 'kodu düzenleyin' veya 'bu komutu çalıştırın') aksi belirtilmedikçe VSCode arayüzünün bu üç bölümüne atıfta bulunur.

Bu kursu kendi başınıza çalışıyorsanız, lütfen daha fazla ayrıntı için [ortam temellerine](../../envsetup/01_setup.md) göz atın.

### Sürüm gereksinimleri

Bu eğitim, Nextflow 25.10.2 veya sonrası için **v2 sözdizimi ayrıştırıcısının ETKİN olduğu** sürüm için tasarlanmıştır.
Yerel veya özel bir ortam kullanıyorsanız, lütfen [burada](../../info/nxf_versions.md) belgelenen doğru ayarları kullandığınızdan emin olun.

## Çalışmaya hazır olun

Codespace'iniz çalışmaya başladığında, eğitime dalmadan önce yapmanız gereken iki şey var: bu özel kurs için çalışma dizininizi ayarlayın ve sağlanan materyallere göz atın.

### Çalışma dizinini ayarlayın

Varsayılan olarak, codespace tüm eğitim kurslarının kökünde ayarlanmış çalışma dizini ile açılır, ancak bu kurs için `nf4-science/genomics/` dizininde çalışacağız.

Terminalde şu komutu çalıştırarak şimdi dizini değiştirin:

```bash
cd nf4-science/genomics/
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz, böylece dosya gezgini kenar çubuğunda yalnızca ilgili dosyalar gösterilir:

```bash
code .
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örn. codespace'iniz uykuya geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak, ona dönmek için her zaman tam yolu kullanabilirsiniz:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Şimdi içeriğe bir göz atalım.

### Sağlanan materyalleri keşfedin

Bu dizinin içeriğini, eğitim çalışma alanının sol tarafındaki dosya gezginini kullanarak keşfedebilirsiniz.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.

Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde temsil etmek için `tree` çıktısını kullanıyoruz, bazen netlik için küçük değişikliklerle.

Burada ikinci seviyeye kadar bir içindekiler tablosu oluşturuyoruz:

```bash
tree . -L 2
```

??? abstract "Dizin içeriği"

    ```console
    .
    ├── data
    │   ├── bam
    │   ├── ref
    │   ├── sample_bams.txt
    │   └── samplesheet.csv
    ├── genomics.nf
    ├── modules
    │   ├── gatk_haplotypecaller.nf
    │   └── samtools_index.nf
    ├── nextflow.config
    └── solutions
        ├── modules
        ├── nf-test.config
        ├── part2
        └── tests

    8 directories, 8 files
    ```

Bölümü genişletmek ve içeriğini görüntülemek için renkli kutuya tıklayın.
Beklenen komut çıktısını ve ayrıca dizin ve dosya içeriklerini özlü bir şekilde görüntülemek için bu şekilde katlanabilir bölümler kullanıyoruz.

- **`genomics.nf` dosyası**, kurs boyunca oluşturacağınız bir iş akışı betiğidir.

- **`modules` dizini**, kurs boyunca dolduracağınız iskelet modül dosyalarını içerir.

- **`nextflow.config` dosyası**, minimal ortam özelliklerini ayarlayan bir yapılandırma dosyasıdır.
  Şimdilik bunu göz ardı edebilirsiniz.

- **`data` dizini**, girdi verilerini ve ilgili kaynakları içerir; kursta daha sonra açıklanacaktır.

- **`solutions` dizini**, tamamlanmış modül dosyalarını ve Bölüm 3 için başlangıç noktası olarak hizmet edebilecek bir Bölüm 2 çözümünü içerir.
  Bunlar, çalışmanızı kontrol etmek ve sorunları gidermek için referans olarak kullanılması amaçlanmıştır.

## Hazırlık kontrol listesi

Dalışa hazır olduğunuzu mu düşünüyorsunuz?

- [ ] Bu kursun hedefini ve ön koşullarını anlıyorum
- [ ] Ortamım çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım

Tüm kutuları işaretleyebiliyorsanız, hazırsınız.

**[Bölüm 1: Metoda genel bakış ve manuel test](./01_method.md) bölümüne devam etmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
