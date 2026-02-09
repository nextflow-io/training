# Kurs özeti

Nextflow Run eğitim kursunu tamamladığınız için tebrikler! 🎉

<!-- placeholder for video -->

## Yolculuğunuz

Çok basit bir iş akışı ile başladınız ve onu çalıştırmayı, çıktıları bulmayı ve yürütmesini yönetmeyi öğrendiniz.
Ardından, bu iş akışının giderek daha karmaşık versiyonları üzerinde çalıştınız ve kanallar ve operatörler, kod modülerleştirme ve konteynerler dahil olmak üzere Nextflow pipeline'larını güçlendiren temel kavramları ve mekanizmaları tanımayı öğrendiniz.
Son olarak, bir pipeline'ın yapılandırmasını tercihlerinize ve hesaplama altyapınıza uyacak şekilde nasıl özelleştireceğinizi öğrendiniz.

### Öğrendikleriniz

Artık Hello pipeline'ının yürütmesini yönetebilir, nasıl yapılandırıldığını açıklayabilir ve ilgili ana kod parçalarını tanımlayabilirsiniz.

- Hello iş akışının son hali, girdi olarak metin selamlamaları içeren bir CSV dosyası alır.
- Dört adım, ayrı modül dosyalarında saklanan Nextflow süreçleri (`sayHello`, `convertToUpper`, `collectGreetings` ve `cowpy`) olarak uygulanmıştır.
- Sonuçlar `results/` adlı bir dizinde yayınlanır.
- Pipeline'ın nihai çıktısı, büyük harfli selamlamaları söyleyen bir karakterin ASCII sanatını içeren düz metin dosyasıdır.

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

- Bir Nextflow iş akışını yerel olarak başlatmak
- Nextflow tarafından oluşturulan çıktıları (sonuçları) ve günlük dosyalarını bulmak ve yorumlamak
- Basit çok adımlı bir iş akışını oluşturan temel Nextflow bileşenlerini tanımak
- Operatörler ve kanal fabrikaları gibi ileri düzey kavramları açıklamak
- Pipeline'ları farklı hesaplama ortamları için yapılandırmak

Artık mevcut Nextflow pipeline'larını kendi çalışmanıza entegre etmeye başlamak için temel bilgiye sahipsiniz.

## Becerilerinizi geliştirmek için sonraki adımlar

Sırada ne yapmanızı önerdiğimize dair en iyi önerilerimiz:

- Sadece Nextflow çalıştırmayın, yazın! [Hello Nextflow](../hello_nextflow/index.md) ile Nextflow geliştiricisi olun
- [Nextflow for Science](../nf4_science/index.md) ile Nextflow'u bilimsel analiz kullanım senaryosuna uygulayın
- [Hello nf-core](../hello_nf-core/index.md) ile nf-core'a başlayın
- [Debugging Side Quest](../side_quests/debugging.md) ile hata ayıklama tekniklerini öğrenin

Son olarak, Nextflow'un yaratıcıları tarafından geliştirilen ve iş akışlarınızı başlatmayı ve yönetmeyi, ayrıca verilerinizi yönetmeyi ve analizleri herhangi bir ortamda etkileşimli olarak çalıştırmayı daha da kolaylaştıran bulut tabanlı bir platform olan [**Seqera Platform**](https://seqera.io/)'a göz atmanızı öneririz.

## Yardım almak

Yardım kaynakları ve topluluk desteği için [Yardım sayfasına](../help.md) bakın.

## Geri bildirim anketi

Devam etmeden önce, lütfen kurs anketini tamamlamak için bir dakikanızı ayırın! Geri bildiriminiz, eğitim materyallerimizi herkes için geliştirmemize yardımcı olur.

[Ankete katılın :material-arrow-right:](survey.md){ .md-button .md-button--primary }
