# Kurs özeti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello Nextflow eğitim kursunu tamamladığınız için tebrikler! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=Lp8hS8RdaMwbp5j5&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [Tüm oynatma listesini Nextflow YouTube kanalında](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) izleyin.

:green_book: Video ile birlikte [video transkriptini](./transcripts/07_next_steps.md) okuyabilirsiniz.
///

## Yolculuğunuz

Sabit kodlanmış bir komut çalıştıran çok basit bir iş akışı ile başladınız.
Altı bölüm boyunca, bu basit iş akışını kanallar, operatörler, konteynerler için yerleşik destek ve yapılandırma seçenekleri dahil olmak üzere Nextflow'un temel özelliklerini kullanan modüler çok adımlı bir boru hattına dönüştürdünüz.

### Ne inşa ettiniz

- Hello iş akışının son hali, metin selamlamaları içeren bir CSV dosyasını girdi olarak alır.
- Dört adım, ayrı modül dosyalarında saklanan Nextflow süreçleri (`sayHello`, `convertToUpper`, `collectGreetings` ve `cowpy`) olarak uygulanmıştır.
- Sonuçlar `results/` adlı bir dizine yayınlanır.
- Boru hattının son çıktısı, büyük harfli selamlamaları söyleyen bir karakterin ASCII sanatını içeren düz metin dosyasıdır.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Her selamlamayı kendi çıktı dosyasına yazar (örn. "Hello-output.txt")
2. **`convertToUpper`:** Her selamlamayı büyük harfe dönüştürür (örn. "HELLO")
3. **`collectGreetings`:** Tüm büyük harfli selamlamaları tek bir toplu dosyada toplar
4. **`cowpy`:** `cowpy` aracını kullanarak ASCII sanatı oluşturur

İş akışı yapılandırması, girdileri ve parametreleri esnek ve tekrarlanabilir bir şekilde sağlamayı destekler.

### Edinilen beceriler

Bu uygulamalı kurs boyunca şunları öğrendiniz:

- Basit bir çok adımlı iş akışı oluşturmaya yetecek temel Nextflow bileşenlerini tanımlama ve kullanma
- Operatörler ve kanal fabrikaları gibi sonraki adım kavramlarını açıklama
- Nextflow iş akışını yerel olarak başlatma
- Nextflow tarafından oluşturulan çıktıları (sonuçlar) ve günlük dosyalarını bulma ve yorumlama
- Temel sorunları giderme

Artık Nextflow'da kendi boru hatlarınızı geliştirmeye başlamak için temel bilgilerle donatıldınız.

## Becerilerinizi geliştirmek için sonraki adımlar

Bundan sonra ne yapılacağına dair en iyi 3 önerimiz:

- [Nextflow for Science](../nf4_science/index.md) ile Nextflow'u bilimsel bir analiz kullanım durumuna uygulayın
- [Hello nf-core](../hello_nf-core/index.md) ile nf-core'a başlayın
- [Side Quests](../side_quests/index.md) ile daha gelişmiş Nextflow özelliklerini keşfedin

Son olarak, Nextflow'un yaratıcıları tarafından geliştirilen, iş akışlarınızı başlatmayı ve yönetmeyi, verilerinizi yönetmeyi ve herhangi bir ortamda etkileşimli analizler çalıştırmayı daha da kolaylaştıran bulut tabanlı bir platform olan [**Seqera Platform**](https://seqera.io/)'a göz atmanızı öneririz.

## Geri bildirim anketi

Devam etmeden önce, lütfen kurs anketini tamamlamak için bir dakikanızı ayırın! Geri bildirimleriniz, eğitim materyallerimizi herkes için geliştirmemize yardımcı olur.

[Ankete katılın :material-arrow-right:](survey.md){ .md-button .md-button--primary }
