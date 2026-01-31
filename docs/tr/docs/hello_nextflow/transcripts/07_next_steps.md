# Sonraki Adımlar - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/xHOcx_4Ancg?si=Lp8hS8RdaMwbp5j5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa sadece transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../05_hello_containers.md) geri dönün.

## Hoş Geldiniz

Tebrikler, başardınız. Hello Nextflow adlı ilk Nextflow eğitim kursunu tamamladınız.

Aferin. Sonuna kadar devam ettiğiniz için teşekkür ederiz ve Nextflow öğrenmek için harcadığınız zaman ve çabayı gerçekten takdir ediyoruz. Bunun çalışmalarınız için faydalı olacağını umuyoruz.

## Sonraki Adımlar

Sonraki adımlar için eğitim portalını takip edin: training.nextflow.io. Oraya sürekli olarak yeni kurs materyalleri ekliyoruz ve bunları yeniliyoruz. Böylece daha ileri düzey eğitimler veya ilgilendiğiniz araştırma alanına özel eğitimler bulabilirsiniz.

Özellikle Nextflow for Science sayfasına göz atın. Bu sayfada, Hello Nextflow'da öğrendiklerinizi belirli kullanım senaryolarına genişleten, bağımsız kısa kurslar serisi bulunmaktadır.

Genomik için bir tane ve ayrıca RNA-seq için bir tane var. Yakında daha fazlasını getirmeyi umuyoruz.

## Yan Görevler

Hello Nextflow'da bahsedebileceğimiz ama çok fazla detay olacak birçok şey var. Bu konuların bazılarını, belirli konularda kısa kurslar olan Yan Görevler içine koyuyoruz.

İlginizi çekebilecek içerikler içeren daha büyük temel eğitim ve ileri düzey eğitim kursları da bulunmaktadır.

## nf-core

Bu kursta bir iki kez bahsedildi, ancak kesinlikle nf-core projesine göz atın. Farklı veri türleri için yüzün üzerinde pipeline var, bu nedenle kendi pipeline'ınızı oluşturmanıza gerek kalmayabilir.

Ayrıca nf-core geliştirici araçlarını kullandığınızda, bir pipeline oluşturabileceğiniz ve bu modülleri saniyeler içinde içe aktarabileceğiniz bine yakın süreç modülü bulunmaktadır.

## Seqera Platform

Son olarak, Seqera Platform için kısa bir tanıtım. Bu, hiç şüphesiz Nextflow'u pratikte çalıştırmanın en iyi yoludur. Bulut tabanlı bir platformdur, ancak kendi hesaplama altyapınızı eklersiniz; ister AWS, Google Batch veya Azure'da kendi bulut hesabınız olsun, ister kendi HPC'niz olsun. Ücretsiz katman herkes için mevcuttur ve akademisyenseniz, ücretsiz profesyonel düzey erişim için akademik programımıza başvurabilirsiniz.

Seqera Platform, iş akışlarını başlatmak ve izlemek için sadece bir grafik arayüzün ötesine geçer. İnteraktif oturumlar çalıştırmak için Data Studios gibi ek araçlar ve bulut üzerinden verilere daha hızlı ve daha ucuz erişim sağlayan Fusion gibi temel araçlar da vardır.

## Destek ve etkinlikler

Unutmayın, herhangi bir sorunla karşılaşırsanız, sadece community.seqera.io adresine gidin. Oradaki forumumuz gerçekten aktiftir, Nextflow topluluğu çok güçlüdür ve yardıma hazır insanlar neredeyse her zaman mevcuttur; ister eğitim, ister Nextflow'un günlük kullanımıyla ilgili herhangi bir şey olsun.

Ve tabii ki, harika bir sonraki adım, nf-core hackathon'larından biri veya Nextflow Summit etkinliklerinden birine katılmaktır. Çok eğlenceliler ve orada sizinle tanışıp Nextflow'u ne için kullandığınız hakkında konuşmak gerçekten güzel olurdu.

## Teşekkürler

Bu eğitim materyalinin yazılmasında yer alan herkese büyük bir teşekkür etmek istiyorum. Ben sundum, ama gerçekten tüm zor işler Seqera'daki eğitim ekibi tarafından yapıldı. Özellikle, bu materyali yeniden yazmak için büyük miktarda çalışma yapan Geraldine'e.

Ayrıca Marcel, Ken, Adam, John, bilimsel geliştirme ekibinden diğerleri ve topluluktaki diğerlerine.

## Geri Bildirim Anketi

Kursu tamamladığınıza göre, ne düşündüğünüzü bilmek isteriz. training.nextflow.io'da, Hello Nextflow bölümünün altında bir geri bildirim anketi bulacaksınız.

Sadece dört soru var ama bizim için gerçekten önemli. Hiçbir şey yoksa, bize kabaca kaç kişinin eğitim aldığını söyler. Ayrıca beğenip beğenmediğinizi de söyler ve herhangi bir öneriniz varsa, lütfen sonunda bunları bırakın. Her bir gönderimi okuyoruz.

Herhangi bir hata görürseniz, her şey GitHub'da açık kaynaklıdır, bu nedenle bir issue oluşturabilir veya bir pull request yapabilir veya forumda bize bir mesaj bırakabilirsiniz. Ne düşündüğünüzü ve nasıl geliştirebileceğimizi duymayı çok isteriz. Tekrar teşekkürler. Yakında görüşmek üzere.
