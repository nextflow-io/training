# Bölüm 3: Kaynak profilleme ve optimizasyon

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

BU BİR YER TUTUCUDUR

!!!note "Not"

    Bu eğitim modülü yeniden geliştirilme aşamasındadır.

---

YAPILACAK

### 1.1. Kaynak kullanım raporu oluşturmak için iş akışını çalıştırın

Nextflow'un raporu otomatik olarak oluşturmasını sağlamak için komut satırınıza `-with-report <dosyaadi>.html` ekleyin.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-1.html
```

Rapor bir html dosyasıdır, tarayıcınızda indirebilir ve açabilirsiniz. Ayrıca soldaki dosya gezgininde sağ tıklayıp `Show preview`'e tıklayarak VS Code içinde görüntüleyebilirsiniz.

Raporu incelemek için birkaç dakika ayırın ve kaynak ayarlamalarını yapabileceğiniz bazı fırsatları belirlemeye çalışın.
Kullanım sonuçlarını ayrılan kaynakların yüzdesi olarak gösteren sekmelere tıkladığınızdan emin olun.
Tüm mevcut özellikleri açıklayan [dokümantasyon](https://www.nextflow.io/docs/latest/reports.html) mevcuttur.

<!-- TODO: insert images -->

Bir gözlem, `GATK_JOINTGENOTYPING`'in CPU açısından oldukça açgözlü göründüğüdür; bu mantıklıdır çünkü çok sayıda karmaşık hesaplama gerçekleştirir.
Bu nedenle bunu artırmayı deneyebilir ve çalışma süresini kısaltıp kısaltmadığına bakabiliriz.

Ancak, bellek tahsislerinde hedefi aşmış görünüyoruz; tüm görevler yalnızca onlara verdiğimizin bir kısmını kullanıyor.
Bunu geri çekmeli ve bazı kaynakları tasarruf etmeliyiz.

### 1.2. Belirli bir görev için kaynak tahsislerini ayarlayın

Belirli bir görev için kaynak tahsisleri belirtmek üzere `withName` process seçicisini kullanabiliriz.
Process bloğunda tek başına olduğunda sözdizimi şu şekilde görünür:

```groovy title="Sözdizimi"
process {
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Bunu `nextflow.config` dosyasındaki mevcut process bloğuna ekleyelim.

```groovy title="nextflow.config" linenums="11"
process {
    // tüm görevler için varsayılanlar
    cpus = 2
    memory = 2.GB
    // belirli bir görev için tahsisler
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Bu şekilde belirtildiğinde, varsayılan ayarlar tüm görevlere uygulanacaktır **hariç** `GATK_JOINTGENOTYPING` görevi, bu özel bir durumdur ve çok daha fazla CPU alır.
Umarım bunun bir etkisi olur.

### 1.3. Değiştirilmiş yapılandırma ile tekrar çalıştırın

İş akışını değiştirilmiş yapılandırma ile ve raporlama bayrağı açık şekilde tekrar çalıştıralım, ancak onları ayırt edebilmek için rapora farklı bir ad verdiğimize dikkat edin.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-2.html
```

Bir kez daha, çalışma süresinde önemli bir fark fark etmeyebilirsiniz, çünkü bu çok küçük bir iş yükü ve araçlar 'gerçek' işi yapmaktan çok yardımcı görevlerde daha fazla zaman harcıyor.

Ancak, ikinci rapor kaynak kullanımımızın artık daha dengeli olduğunu gösteriyor.

<!-- **TODO: screenshots?** -->

Gördüğünüz gibi, bu yaklaşım görevleriniz farklı kaynak gereksinimleri olduğunda faydalıdır. Tahminlere değil, gerçek verilere dayanarak her görev için ayarladığınız kaynak tahsislerini doğru boyutlandırmanızı sağlar.

!!!note "Not"

    Bu, kaynaklarınızı optimize etmek için yapabileceklerinizin sadece küçük bir tadımlığıdır.
    Nextflow'un kendisinde kaynak sınırlamaları nedeniyle başarısız olan işleri yeniden denemek için gerçekten güzel bir [dinamik yeniden deneme mantığı](https://www.nextflow.io/docs/latest/process.html#dynamic-task-resources) bulunmaktadır.
    Ek olarak, Seqera Platform kaynak tahsislerinizi otomatik olarak optimize etmek için yapay zeka destekli araçlar sunar.

    Bu eğitim kursunun yaklaşan bir bölümünde bu yaklaşımların her ikisini de ele alacağız.

Bununla birlikte, hangi computing executor ve bilgi işlem altyapısını kullandığınıza bağlı olarak tahsis edebileceğiniz (veya etmeniz gereken) konusunda bazı kısıtlamalar olabilir. Örneğin, kümeniz başka yerde çalıştırırken geçerli olmayan belirli sınırlar içinde kalmanızı gerektirebilir.
