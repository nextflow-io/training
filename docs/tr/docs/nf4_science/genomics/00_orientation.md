# Başlarken

## Eğitim ortamını başlatın

GitHub Codespaces üzerinde sağladığımız önceden hazırlanmış ortamı kullanmak için aşağıdaki "GitHub Codespaces'te Aç" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../../envsetup/index.md) bölümüne bakın.

Eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (ekipmanınıza bağlı olarak sağ tıklama, ctrl-tıklama veya cmd-tıklama kullanın), böylece ortam yüklenirken okumaya devam edebilirsiniz.
Kursu tamamlamak için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, eğitim kursunda çalışmak için gerekli tüm yazılımları, kodları ve verileri içerir, bu nedenle kendiniz herhangi bir şey yüklemenize gerek yoktur.

Codespace, bir dosya sistemi gezgini, bir kod düzenleyici ve bir terminal kabuğu içeren bir VSCode arayüzü ile kurulmuştur.
Kurs boyunca verilen tüm talimatlar (örneğin 'dosyayı açın', 'kodu düzenleyin' veya 'bu komutu çalıştırın'), aksi belirtilmedikçe VSCode arayüzünün bu üç bölümüne atıfta bulunur.

Bu kursu kendi başınıza tamamlıyorsanız, lütfen daha fazla ayrıntı için [ortam temelleri](../../envsetup/01_setup.md) ile kendinizi tanıştırın.

### Sürüm gereksinimleri

Bu eğitim, Nextflow 25.10.2 veya sonrası için **v2 sözdizimi ayrıştırıcısı ETKİNLEŞTİRİLMİŞ** olarak tasarlanmıştır.
Yerel veya özel bir ortam kullanıyorsanız, lütfen [burada](../../info/nxf_versions.md) belgelenen doğru ayarları kullandığınızdan emin olun.

## Çalışmaya hazırlanın

Codespace'iniz çalışmaya başladığında, eğitime dalmadan önce yapmanız gereken iki şey var: bu özel kurs için çalışma dizininizi ayarlayın ve sağlanan materyallere bir göz atın.

### Çalışma dizinini ayarlayın

Varsayılan olarak, codespace tüm eğitim kurslarının kök dizininde ayarlanmış çalışma dizini ile açılır, ancak bu kurs için `nf4-science/genomics/` dizininde çalışacağız.

Terminalde şu komutu çalıştırarak şimdi dizini değiştirin:

```bash
cd nf4-science/genomics/
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz, böylece dosya gezgini kenar çubuğunda yalnızca ilgili dosyalar görünür:

```bash
code .
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örneğin codespace'iniz uyku moduna geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak, ona dönmek için her zaman tam yolu kullanabilirsiniz:

    ```bash
    cd /workspaces/training/nf4-science/genomics
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
Beklenen komut çıktısını ve ayrıca dizin ve dosya içeriklerini özlü bir şekilde görüntülemek için bunun gibi daraltılabilir bölümler kullanırız.

- **`genomics.nf` dosyası**, kurs boyunca oluşturacağınız bir iş akışı betiğidir.

- **`modules` dizini**, kurs sırasında dolduracağınız iskelet modül dosyalarını içerir.

- **`nextflow.config` dosyası**, minimum ortam özelliklerini ayarlayan bir yapılandırma dosyasıdır.
  Şimdilik göz ardı edebilirsiniz.

- **`data` dizini**, kursun ilerleyen bölümlerinde açıklanan girdi verilerini ve ilgili kaynakları içerir.

- **`solutions` dizini**, tamamlanmış modül dosyalarını ve Bölüm 3 için başlangıç noktası olarak kullanılabilecek bir Bölüm 2 çözümünü içerir.
  Çalışmanızı kontrol etmek ve sorunları gidermek için bir referans olarak kullanılmaları amaçlanmıştır.

## Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu mu düşünüyorsunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Ortamım çalışır durumda
- [ ] Çalışma dizinini uygun şekilde ayarladım

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

**[Bölüm 1: Yönteme genel bakış ve manuel test](./01_method.md) bölümüne devam etmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
