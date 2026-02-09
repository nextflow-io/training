# Yönlendirme

Eğitim ortamı, bu eğitim kursunu tamamlamak için gerekli tüm yazılımı, kodu ve veriyi içerir, bu nedenle kendiniz herhangi bir şey yüklemenize gerek yoktur.
Ancak, giriş yapmak için (ücretsiz) bir hesaba ihtiyacınız vardır ve arayüzü tanımak için birkaç dakikanızı ayırmalısınız.

Henüz yapmadıysanız, lütfen daha ileri gitmeden önce [Ortam Kurulumu](../../envsetup/) mini-kursunu tamamlayın.

## Sağlanan materyaller

Bu eğitim kursu boyunca, `nf4-science/rnaseq/` dizininde çalışacağız. Eğitim çalışma alanını açtığınızda bu dizine geçmeniz gerekir.
Bu dizin, ihtiyaç duyacağınız tüm kod dosyalarını, test verilerini ve yardımcı dosyaları içerir.

Bu dizinin içeriğini keşfetmekten çekinmeyin; bunu yapmanın en kolay yolu, VSCode arayüzündeki eğitim çalışma alanının sol tarafındaki dosya gezginini kullanmaktır.
Alternatif olarak, `tree` komutunu kullanabilirsiniz.
Kurs boyunca, dizin yapısını ve içeriğini okunabilir bir biçimde temsil etmek için `tree` çıktısını kullanıyoruz, bazen netlik için küçük değişikliklerle.

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

    Bu çok fazla görünüyorsa endişelenmeyin; kursun her adımında ilgili parçaları inceleyeceğiz.
    Bu sadece size genel bir bakış sağlamak içindir.

**Başlamak için bilmeniz gerekenler şunlardır:**

- **`rnaseq.nf` dosyası**, geliştireceğimiz iş akışı betiğinin taslağıdır.

- **`nextflow.config` dosyası**, minimum ortam özelliklerini ayarlayan bir yapılandırma dosyasıdır. Şimdilik bunu görmezden gelebilirsiniz.

- **`data` dizini**, girdi verilerini ve ilgili kaynakları içerir:

  - İnsan kromozom 20'sinin küçük bir bölgesinden (hg19/b37'den) oluşan `genome.fa` adlı _bir referans genom_.
  - Dosya boyutlarını küçük tutmak için küçük bir bölgeye indirgenmiş _RNAseq verileri_, `reads/` dizininde.
  - Toplu işleme için örnek veri dosyalarının kimliklerini ve yollarını listeleyen _CSV dosyaları_.

- **`solutions` dizini**, kursun her adımından kaynaklanan tamamlanmış iş akışı betiklerini ve modüllerini içerir.
  Bunlar, çalışmanızı kontrol etmek ve sorunları gidermek için bir referans olarak kullanılmak üzere tasarlanmıştır.
  Dosya adındaki numara, kursun ilgili bölümünün adımına karşılık gelir.

!!!tip

    Herhangi bir nedenle bu dizinden çıkarsanız, her zaman bu komutu çalıştırarak geri dönebilirsiniz:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Şimdi, kursa başlamak için bu sayfanın sağ alt köşesindeki oka tıklayın.
