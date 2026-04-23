# Kurs Özeti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow for RNAseq eğitim kursunu tamamladığınız için tebrikler!

## Yolculuğunuz

Metodolojisini anlamak için terminalde RNAseq işleme araçlarını manuel olarak çalıştırarak başladınız.
Ardından süreci otomatikleştirmek için tek örnekli bir Nextflow pipeline'ı oluşturdunuz, bunu paralel olarak birden fazla örneği işleyecek şekilde ölçeklendirdiniz ve paired-end veriyi işlemek ve örnekler arasında QC raporlarını toplamak için genişlettiniz.

### Oluşturduklarınız

- Girdi olarak FASTQ dosyalarını alan ve çıktı olarak kırpılmış okumalar, hizalamalar ve toplanmış QC raporları üreten bir RNAseq işleme pipeline'ı.
- Ayrı modül dosyalarında saklanan kırpma (Trim Galore), hizalama (HISAT2), kalite kontrol (FastQC) ve rapor toplama (MultiQC) için süreçler.
- Pipeline, Nextflow'un veri akışı paradigmasını kullanarak girdi örneklerinin işlenmesini otomatik olarak paralelleştirir.
- Son pipeline, paired-end dizileme verisini işler.

### Kazanılan Beceriler

Bu uygulamalı kurs boyunca şunları öğrendiniz:

- Temel RNAseq işleme ve QC yöntemlerini uygulamak için doğrusal bir iş akışı yazmak
- FASTQ ve referans genom kaynakları gibi alana özgü dosyaları uygun şekilde işlemek
- Single-end ve paired-end dizileme verisini işlemek
- Örnek başına RNAseq işlemeyi paralelleştirmek için Nextflow'un veri akışı paradigmasından yararlanmak
- İlgili kanal operatörlerini kullanarak birden fazla adım ve örnek arasında QC raporlarını toplamak

Artık kendi çalışmanızda RNAseq analiz iş akışlarına Nextflow'u uygulamaya başlamak için donanımlısınız.

## Becerilerinizi Geliştirmek İçin Sonraki Adımlar

Sırada ne yapacağınıza dair en iyi önerilerimiz:

- [Nextflow for Science](../index.md) ile Nextflow'u diğer bilimsel analiz kullanım durumlarına uygulayın
- [Hello nf-core](../../hello_nf-core/index.md) ile nf-core'a başlayın
- [Side Quests](../../side_quests/index.md) ile daha gelişmiş Nextflow özelliklerini keşfedin

Son olarak, Nextflow'un yaratıcıları tarafından geliştirilen ve iş akışlarınızı başlatmayı ve yönetmeyi, ayrıca verilerinizi yönetmeyi ve herhangi bir ortamda etkileşimli olarak analiz çalıştırmayı daha da kolay hale getiren bulut tabanlı bir platform olan [**Seqera Platform**](https://seqera.io/)'a göz atmanızı öneririz.

## Yardım Alma

Yardım kaynakları ve topluluk desteği için [Yardım sayfası](../../help.md)'na bakın.

## Geri Bildirim Anketi

Devam etmeden önce, lütfen kurs anketini tamamlamak için bir dakikanızı ayırın! Geri bildiriminiz, eğitim materyallerimizi herkes için geliştirmemize yardımcı olur.

[Ankete katılın :material-arrow-right:](survey.md){ .md-button .md-button--primary }
