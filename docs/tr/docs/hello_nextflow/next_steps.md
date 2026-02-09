# Kurs özeti

Hello Nextflow eğitim kursunu tamamladığınız için tebrikler! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [Nextflow YouTube kanalındaki tüm oynatma listesini](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) izleyin.

:green_book: Videoyla birlikte [video transkriptini](./transcripts/07_next_steps.md) okuyabilirsiniz.
///

## Yolculuğunuz

Sabit kodlanmış bir komut çalıştıran çok basit bir iş akışı ile başladınız.
Altı bölüm boyunca, bu basit iş akışını kanallar, operatörler, konteynerler için yerleşik destek ve yapılandırma seçenekleri gibi Nextflow'un temel özelliklerini kullanan modüler çok adımlı bir boru hattına dönüştürdünüz.

### Ne oluşturdunuz

- Hello iş akışının son hali, girdi olarak metin selamlamaları içeren bir CSV dosyası alır.
- Dört adım, ayrı modül dosyalarında saklanan Nextflow süreçleri (`sayHello`, `convertToUpper`, `collectGreetings` ve `cowpy`) olarak uygulanmıştır.
- Sonuçlar `results/` adlı bir dizine yayınlanır.
- Boru hattının nihai çıktısı, büyük harfli selamlamaları söyleyen bir karakterin ASCII sanatını içeren düz metin dosyasıdır.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Her selamlamayı kendi çıktı dosyasına yazar (_örn._ "Hello-output.txt")
2. **`convertToUpper`:** Her selamlamayı büyük harfe dönüştürür (_örn._ "HELLO")
3. **`collectGreetings`:** Tüm büyük harfli selamlamaları tek bir toplu dosyada toplar
4. **`cowpy`:** `cowpy` aracını kullanarak ASCII sanatı oluşturur

İş akışı yapılandırması, girdileri ve parametreleri esnek ve tekrarlanabilir bir şekilde sağlamayı destekler.

### Kazanılan beceriler

Bu uygulamalı kurs boyunca şunları öğrendiniz:

- Basit bir çok adımlı iş akışı oluşturmak için yeterli temel Nextflow bileşenlerini tanımlamak ve kullanmak
- Operatörler ve kanal fabrikaları gibi sonraki adım kavramlarını tanımlamak
- Bir Nextflow iş akışını yerel olarak başlatmak
- Nextflow tarafından oluşturulan çıktıları (sonuçları) ve günlük dosyalarını bulmak ve yorumlamak
- Temel sorunları gidermek

Artık Nextflow'da kendi boru hatlarınızı geliştirmeye başlamak için temel bilgiye sahipsiniz.

## Becerilerinizi geliştirmek için sonraki adımlar

Sırada ne yapacağınıza dair en iyi 3 önerimiz:

- [Nextflow for Science](../nf4_science/index.md) ile Nextflow'u bilimsel bir analiz kullanım senaryosuna uygulayın
- [Hello nf-core](../hello_nf-core/index.md) ile nf-core'a başlayın
- [Side Quests](../side_quests/index.md) ile daha gelişmiş Nextflow özelliklerini keşfedin

Son olarak, Nextflow'un yaratıcıları tarafından geliştirilen ve iş akışlarınızı başlatmayı ve yönetmeyi, ayrıca verilerinizi yönetmeyi ve analizleri herhangi bir ortamda etkileşimli olarak çalıştırmayı daha da kolaylaştıran bulut tabanlı bir platform olan [**Seqera Platform**](https://seqera.io/)'a göz atmanızı öneririz.

## Geri bildirim anketi

Devam etmeden önce, lütfen kurs anketini tamamlamak için bir dakikanızı ayırın! Geri bildiriminiz, eğitim materyallerimizi herkes için geliştirmemize yardımcı olur.

[Ankete katılın :material-arrow-right:](survey.md){ .md-button .md-button--primary }
