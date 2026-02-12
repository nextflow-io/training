# Başlarken

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eğitim ortamını başlatın

GitHub Codespaces üzerinde sağladığımız önceden hazırlanmış ortamı kullanmak için aşağıdaki "GitHub Codespaces'te Aç" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../../envsetup/index.md) sayfasına bakın.

Eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (ekipmanınıza bağlı olarak sağ tıklama, ctrl-tıklama veya cmd-tıklama kullanın), böylece ortam yüklenirken okumaya devam edebilirsiniz.
Kursu tamamlamak için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam temelleri

Bu eğitim ortamı, eğitim kursunu tamamlamak için gerekli tüm yazılımı, kodu ve veriyi içermektedir, dolayısıyla kendiniz herhangi bir şey yüklemenize gerek yoktur.

Codespace, bir dosya sistemi gezgini, bir kod düzenleyici ve bir terminal kabuğu içeren bir VSCode arayüzü ile kurulmuştur.
Kurs boyunca verilen tüm talimatlar (örneğin 'dosyayı açın', 'kodu düzenleyin' veya 'bu komutu çalıştırın'), aksi belirtilmedikçe VSCode arayüzünün bu üç bölümüne atıfta bulunur.

Bu kursu kendi başınıza tamamlıyorsanız, lütfen daha fazla ayrıntı için [ortam temellerine](../../envsetup/01_setup.md) göz atın.

### Sürüm gereksinimleri

Bu eğitim, Nextflow 25.10.2 veya sonrası için **v2 sözdizimi ayrıştırıcısı ETKİNLEŞTİRİLMİŞ** olarak tasarlanmıştır.
Yerel veya özel bir ortam kullanıyorsanız, lütfen [burada](../../info/nxf_versions.md) belgelenen doğru ayarları kullandığınızdan emin olun.

## Çalışmaya hazırlanın

Codespace'iniz çalışmaya başladığında, eğitime dalmadan önce yapmanız gereken iki şey var: bu özel kurs için çalışma dizininizi ayarlayın ve sağlanan materyallere bir göz atın.

### Çalışma dizinini ayarlayın

Varsayılan olarak, codespace tüm eğitim kurslarının kök dizininde açılır, ancak bu kurs için `nf4-science/rnaseq/` dizininde çalışacağız.

Terminalde şu komutu çalıştırarak şimdi dizini değiştirin:

```bash
cd nf4-science/rnaseq/
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz, böylece dosya gezgini kenar çubuğunda yalnızca ilgili dosyalar görünür:

```bash
code .
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örneğin codespace'iniz uyku moduna geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak, oraya geri dönmek için her zaman tam yolu kullanabilirsiniz:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Şimdi içeriğe bir göz atalım.

### Sağlanan materyalleri keşfedin

Bu dizinin içeriğini, eğitim çalışma alanının sol tarafındaki dosya gezginini kullanarak keşfedebilirsiniz.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.

Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde göstermek için `tree` çıktısını kullanıyoruz, bazen netlik için küçük değişikliklerle.

Burada üçüncü seviyeye kadar bir içindekiler tablosu oluşturuyoruz:

```bash
tree . -L 3
```

??? abstract "Dizin içeriği"

    ```console
    .
    ├── data
    │   ├── genome.fa
    │   ├── paired-end.csv
    │   ├── reads
    │   │   ├── ENCSR000COQ1_1.fastq.gz
    │   │   ├── ENCSR000COQ1_2.fastq.gz
    │   │   ├── ENCSR000COQ2_1.fastq.gz
    │   │   ├── ENCSR000COQ2_2.fastq.gz
    │   │   ├── ENCSR000COR1_1.fastq.gz
    │   │   ├── ENCSR000COR1_2.fastq.gz
    │   │   ├── ENCSR000COR2_1.fastq.gz
    │   │   ├── ENCSR000COR2_2.fastq.gz
    │   │   ├── ENCSR000CPO1_1.fastq.gz
    │   │   ├── ENCSR000CPO1_2.fastq.gz
    │   │   ├── ENCSR000CPO2_1.fastq.gz
    │   │   └── ENCSR000CPO2_2.fastq.gz
    │   └── single-end.csv
    ├── nextflow.config
    ├── rnaseq.nf
    └── solutions
        ├── modules
        │   ├── fastqc.nf
        │   ├── fastqc_pe.nf
        │   ├── hisat2_align.nf
        │   ├── hisat2_align_pe.nf
        │   ├── multiqc.nf
        │   ├── trim_galore.nf
        │   └── trim_galore_pe.nf
        ├── rnaseq-2.1.nf
        ├── rnaseq-2.2.nf
        ├── rnaseq-2.3.nf
        ├── rnaseq-3.1.nf
        ├── rnaseq-3.2.nf
        └── rnaseq_pe-3.3.nf
    ```

Bölümü genişletmek ve içeriğini görüntülemek için renkli kutuya tıklayın.
Beklenen komut çıktısını ve dizin ve dosya içeriklerini özlü bir şekilde görüntülemek için bu gibi daraltılabilir bölümler kullanıyoruz.

- **`rnaseq.nf` dosyası**, kurs boyunca oluşturacağınız bir iş akışı betiği için bir iskelet içerir.

- **`modules` dizini**, kurs boyunca dolduracağınız süreç modülleri için iskeletler içerir.

- **`nextflow.config` dosyası**, minimum ortam özelliklerini ayarlayan bir yapılandırma dosyasıdır.
  Şimdilik bunu göz ardı edebilirsiniz.

- **`data` dizini**, girdi verilerini ve ilgili kaynakları içerir, kursun ilerleyen bölümlerinde açıklanacaktır.

- **`solutions` dizini**, kursun her adımından kaynaklanan tamamlanmış iş akışı betiklerini ve modüllerini içerir.
  Bunlar, çalışmanızı kontrol etmek ve sorunları gidermek için referans olarak kullanılmak üzere tasarlanmıştır.
  Bölüm 2 çözümü, Bölüm 3 için bir başlangıç noktası olarak kullanılabilir.

## Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Ortamım çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım

Tüm kutuları işaretleyebiliyorsanız, hazırsınız demektir.

**[Bölüm 1: Yönteme genel bakış](./01_method.md) sayfasına devam etmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
