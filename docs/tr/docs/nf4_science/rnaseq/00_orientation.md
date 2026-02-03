# Oryantasyon

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Eğitim ortamı, bu eğitim kursunu tamamlamak için gerekli tüm yazılımı, kodu ve veriyi içermektedir, dolayısıyla kendiniz herhangi bir şey yüklemenize gerek yoktur.
Ancak, giriş yapmak için bir (ücretsiz) hesaba ihtiyacınız var ve arayüze aşina olmak için birkaç dakikanızı ayırmalısınız.

Henüz yapmadıysanız, lütfen daha ileri gitmeden önce [Ortam Kurulumu](../../envsetup/) mini kursunu tamamlayın.

## Sağlanan materyaller

Bu eğitim kursu boyunca, eğitim çalışma alanını açtığınızda içine girmeniz gereken `nf4-science/rnaseq/` dizininde çalışacağız.
Bu dizin, ihtiyaç duyacağınız tüm kod dosyalarını, test verilerini ve yardımcı dosyaları içermektedir.

Bu dizinin içeriğini keşfetmekten çekinmeyin; bunu yapmanın en kolay yolu VSCode arayüzündeki eğitim çalışma alanının sol tarafındaki dosya gezginini kullanmaktır.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.
Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde göstermek için `tree` çıktısını kullanıyoruz, bazen netlik için küçük değişikliklerle.

Burada ikinci seviyeye kadar bir içindekiler tablosu oluşturuyoruz:

```bash
tree . -L 3
```

??? success "Dizin içeriği"

    ```console
    rnaseq
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

!!!note

    Bu çok fazla gibi görünüyorsa endişelenmeyin; kursun her adımında ilgili parçaları ele alacağız.
    Bu sadece size bir genel bakış sağlamak için.

**İşte başlamak için bilmeniz gerekenlerın bir özeti:**

- **`rnaseq.nf` dosyası**, geliştireceğimiz iş akışı betiğinin ana hatlarıdır.

- **`nextflow.config` dosyası**, minimum ortam özelliklerini ayarlayan bir yapılandırma dosyasıdır. Şimdilik bunu göz ardı edebilirsiniz.

- **`data` dizini**, girdi verilerini ve ilgili kaynakları içerir:

  - İnsan kromozomu 20'nin küçük bir bölgesinden oluşan (hg19/b37'den) `genome.fa` adlı _bir referans genom_.
  - Dosya boyutlarını küçük tutmak için küçük bir bölgeye indirilmiş, `reads/` dizinindeki _RNAseq verileri_.
  - Örnek veri dosyalarının kimliklerini ve yollarını listeleyen, toplu işleme için _CSV dosyaları_.

- **`solutions` dizini**, kursun her adımından kaynaklanan tamamlanmış iş akışı betiklerini ve modüllerini içerir.
  Bunlar, çalışmanızı kontrol etmek ve sorunları gidermek için referans olarak kullanılmak üzere tasarlanmıştır.
  Dosya adındaki numara, kursun ilgili bölümünün adımına karşılık gelir.

!!!tip

    Herhangi bir nedenle bu dizinden çıkarsanız, oraya geri dönmek için her zaman şu komutu çalıştırabilirsiniz:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Şimdi, kursa başlamak için bu sayfanın sağ alt köşesindeki oka tıklayın.
