# Kurs özeti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow Run eğitim kursunu tamamladığınız için tebrikler!

<!-- placeholder for video -->

## Yolculuğunuz

Çok temel bir workflow ile başladınız ve onu çalıştırmayı, çıktıları bulmayı ve çalışmasını yönetmeyi öğrendiniz.
Ardından, bu workflow'un giderek daha karmaşık versiyonları üzerinde çalıştınız ve channel'lar ile operatörler, kod modülerleştirme ve konteynerlar dahil Nextflow pipeline'larına güç veren temel kavramları ve mekanizmaları tanımayı öğrendiniz.
Son olarak, tercihlerinize ve hesaplama altyapınıza uyması için bir pipeline'ın yapılandırmasını nasıl özelleştireceğinizi öğrendiniz.

### Ne öğrendiniz

Artık Hello pipeline'ının çalışmasını yönetebilir, nasıl yapılandırıldığını tanımlayabilir ve dahil olan ana kod parçalarını tanımlayabilirsiniz.

- Hello workflow'unun son hali, girdi olarak metin selamlamaları içeren bir CSV dosyası alır.
- Dört adım, ayrı modül dosyalarında saklanan Nextflow process'leri olarak uygulanır (`sayHello`, `convertToUpper`, `collectGreetings` ve `cowpy`).
- Sonuçlar `results/` adlı bir dizine yayınlanır.
- Pipeline'ın nihai çıktısı, büyük harfli selamlamaları söyleyen bir karakterin ASCII art'ını içeren düz metin dosyasıdır.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Her selamlamayı kendi çıktı dosyasına yazar (örn. "Hello-output.txt")
2. **`convertToUpper`:** Her selamlamayı büyük harfe dönüştürür (örn. "HELLO")
3. **`collectGreetings`:** Tüm büyük harfli selamlamaları tek bir batch dosyasına toplar
4. **`cowpy`:** `cowpy` aracını kullanarak ASCII art oluşturur

Workflow yapılandırması, girdilerin ve parametrelerin esnek, tekrar üretilebilir bir şekilde sağlanmasını destekler.

### Kazanılan beceriler

Bu uygulamalı kurs sayesinde şunları öğrendiniz:

- Nextflow workflow'unu yerel olarak başlatma
- Nextflow tarafından oluşturulan çıktıları (sonuçlar) ve log dosyalarını bulma ve yorumlama
- Basit çok adımlı bir workflow'u oluşturan temel Nextflow bileşenlerini tanıma
- Operatörler ve channel factory'leri gibi ileri adım kavramlarını tanımlama
- Pipeline'ları farklı hesaplama ortamları için yapılandırma

Artık mevcut Nextflow pipeline'larını kendi çalışmanıza entegre etmeye başlamak için temel bilgilerle donatıldınız.

## Becerilerinizi geliştirmek için sonraki adımlar

Sırada ne yapacağınız için en iyi önerilerimiz:

- Sadece Nextflow çalıştırmayın, yazın! [Hello Nextflow](../hello_nextflow/index.md) ile Nextflow geliştiricisi olun
- Bilimsel analiz kullanım durumlarına Nextflow uygulayın: [Nextflow for Science](../nf4_science/index.md)
- [Hello nf-core](../hello_nf-core/index.md) ile nf-core'a başlayın
- [Debugging Side Quest](../side_quests/debugging.md) ile hata ayıklama tekniklerini öğrenin

Son olarak, Nextflow'un yaratıcıları tarafından geliştirilen bulut tabanlı bir platform olan [**Seqera Platform**](https://seqera.io/)'a bakmanızı öneririz; workflow'larınızı başlatmayı ve yönetmeyi, verilerinizi yönetmeyi ve herhangi bir ortamda etkileşimli analizler çalıştırmayı çok daha kolay hale getirir.

## Yardım alma

Yardım kaynakları ve topluluk desteği için [Yardım sayfasına](../help.md) bakın.

## Geri bildirim anketi

Devam etmeden önce, lütfen kurs anketini doldurmak için bir dakikanızı ayırın! Geri bildirimleriniz, eğitim materyallerimizi herkes için iyileştirmemize yardımcı olur.

[Ankete katıl :material-arrow-right:](survey.md){ .md-button .md-button--primary }
