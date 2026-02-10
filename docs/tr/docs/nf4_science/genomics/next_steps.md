# Kurs özeti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow for Genomics eğitim kursunu tamamladığınız için tebrikler! 🎉

## Yolculuğunuz

Metodolojiyi anlamak için terminalde varyant çağırma araçlarını manuel olarak çalıştırarak başladınız.
Ardından süreci otomatikleştirmek için tek örnekli bir Nextflow pipeline'ı oluşturdunuz, bunu paralel olarak birden fazla örneği işleyecek şekilde ölçeklendirdiniz ve kanal operatörleri kullanarak çok örnekli ortak genotipleme eklediniz.

### Oluşturduklarınız

- Girdi olarak BAM dosyalarını alan ve çıktı olarak ortak çağrılmış VCF'ler üreten bir varyant çağırma pipeline'ı.
- Ayrı modül dosyalarında saklanan üç süreç (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER` ve `GATK_JOINTGENOTYPING`).
- Nextflow'un veri akışı paradigmasını kullanarak herhangi bir sayıda girdi örneğine otomatik olarak ölçeklenen pipeline.
- Sonuçlar `results/` adlı bir dizine yayınlanır.

### Kazanılan beceriler

Bu uygulamalı kurs boyunca şunları öğrendiniz:

- Tek bir örneğe varyant çağırma uygulamak için doğrusal bir iş akışı yazmak
- İndeks dosyaları ve referans genom kaynakları gibi yardımcı dosyaları uygun şekilde işlemek
- Örnek başına varyant çağırmayı paralelleştirmek için Nextflow'un veri akışı paradigmasından yararlanmak
- İlgili kanal operatörlerini kullanarak çok örnekli ortak çağırma uygulamak

Artık kendi çalışmalarınızda Nextflow'u genomik analiz iş akışlarına uygulamaya başlamak için donanımlısınız.

## Becerilerinizi geliştirmek için sonraki adımlar

Sırada ne yapmanızı önerdiğimize dair en iyi önerilerimiz:

- [Nextflow for Science](../index.md) ile Nextflow'u diğer bilimsel analiz kullanım durumlarına uygulayın
- [Hello nf-core](../../hello_nf-core/index.md) ile nf-core'a başlayın
- [Side Quests](../../side_quests/index.md) ile daha gelişmiş Nextflow özelliklerini keşfedin

Son olarak, Nextflow'un yaratıcıları tarafından geliştirilen ve iş akışlarınızı başlatmayı ve yönetmeyi, ayrıca verilerinizi yönetmeyi ve analizleri herhangi bir ortamda etkileşimli olarak çalıştırmayı daha da kolaylaştıran bulut tabanlı bir platform olan [**Seqera Platform**](https://seqera.io/)'a göz atmanızı öneririz.

## Yardım almak

Yardım kaynakları ve topluluk desteği için [Yardım sayfasına](../../help.md) bakın.

## Geri bildirim anketi

Devam etmeden önce, lütfen kurs anketini tamamlamak için bir dakikanızı ayırın! Geri bildiriminiz, eğitim materyallerimizi herkes için geliştirmemize yardımcı olur.

[Ankete katılın :material-arrow-right:](survey.md){ .md-button .md-button--primary }
