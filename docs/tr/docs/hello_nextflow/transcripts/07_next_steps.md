# Sonraki Adımlar - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=Lp8hS8RdaMwbp5j5&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../next_steps.md) geri dönün.

## Hoş Geldiniz

​

Tebrikler, başardınız!

Sona ulaştınız ve Hello Nextflow eğitim kursunu tamamladınız. Umarız keyif almışsınızdır. Sonuna kadar bizimle kaldığınız için çok teşekkür ederiz ve Nextflow öğrenmek için harcadığınız zaman ve çabayı gerçekten takdir ediyoruz. Umarız çalışmalarınız için faydalı olacaktır.

## training.nextflow.io'daki Diğer Kurslar

training.nextflow.io'ya tekrar gelmeyi unutmayın. Sürekli olarak yeni kısa kurslar ekliyoruz ve ayrıca burada zaten mevcut olan materyallerin çoğunu yeniliyoruz. Dolayısıyla bu Hello Nextflow eğitim kursu zaman içinde güncellenecektir.

Bu özellikle önemlidir çünkü Nextflow'daki sözdizimini güncelliyoruz ve 2026'da oldukça fazla yeni özellik gelecek, bu nedenle bu kurs 2027'de bir sonraki sefer yaptığımızda biraz farklı görünecek ve hissettirecektir.

Özellikle "Nextflow for Science" sayfasına dikkat çekmek istiyorum. Bunlar kısa kurslardır ve bu Hello Nextflow kursunu takip edecek şekilde tasarlanmıştır. Ve Nextflow'un genomik veya RNAseq gibi farklı kullanım durumlarıyla nasıl kullanılacağını gösterirler. Sürekli olarak daha fazla bilimsel kullanım durumu eklemeye çalışıyoruz.

Ayrıca Yan Görevler de var. Hello Nextflow gibi bir kurs geliştirdiğimizde, kapsayabileceğimiz çok şey var ve her şeyi kapsam içinde tutmak zor. Bu nedenle, insanlar için ilginç olduğunu düşündüğümüz ve daha derinlemesine ele alınmayı hak eden belirli bir konu varsa, bunu bir Yan Göreve koyuyoruz.

Gidin ve inceleyin ve çalışmanızla ilgili olabilecek farklı şeyler varsa, nf-test veya meta verilerle farklı şeyler yapmak ve yaygın betik kalıpları gibi, Yan Görevlere göz atın ve daha fazla öğrenmek için faydalı olup olmayacağını görün.

Ayrıca nf-core hakkında bir kurs var. Umarız bu noktada projeye aşinasınızdır, ancak değilseniz, gidin ve inceleyin. Farklı analiz türleri ve farklı veri türleri için yaklaşık 150 farklı boru hattı var, bu nedenle ihtiyaç duyduğunuz veri analizi türü için kullanıma hazır bir boru hattının olması tamamen mümkün.

Önemli olarak, nf-core'da ayrıca bileşenler var, yaklaşık 1700 farklı modül, farklı süreçler ve araçlar için sarmalayıcılar. Ve nf-core ile birlikte gelen araçlarla, bunları karıştırıp eşleştirebilir ve Lego parçaları gibi kendi boru hattınızı oluşturabilirsiniz. Çok daha hızlı ve daha tekrarlanabilir.

## Seqera Platform

Nextflow kullanımınızı ölçeklendirdikçe, Seqera Platform'u inceleyin, Nextflow'u çalıştırmanın en iyi yolu. Kendi altyapınızda çalıştırabilirsiniz, yani HPC veya AWS, Azure, Google Cloud, Oracle ve daha fazlası. Ayrıca herhangi bir bilgi işlem altyapısını yönetmek istemiyorsanız kendi Seqera Compute'umuzu da kullanabilirsiniz.

Seqera Platform, sizin için ortamı oluşturan Batch Forge gibi özelliklerle bu karmaşık bulut altyapılarının kurulumunu gerçekten basitleştirir. Ve ayrıca gözlemlenebilirlik, denetim günlüğü ve uyumluluk konusunda gerçekten yardımcı olur.

Disk erişimini ve veri aktarımlarını optimize eden Fusion gibi teknolojilerle boru hatlarının nesnel olarak daha ucuz ve daha hızlı çalışmasını sağlar. Ve ayrıca boru hatlarınızın yapılandırmasının mümkün olduğunca sıkı ayarlandığından emin olmak için boru hattı optimizasyonu vardır.

Boru hatlarını çalıştırmanın yanı sıra tamamen farklı özellikler de var. Etkileşimli analizler çalıştırabileceğiniz ve oluşturduğunuz herhangi bir özel docker imajından ortamlar oluşturabileceğiniz Studios'umuz var. Ve farklı dosya sistemlerinizi nerede olurlarsa olsunlar keşfetmenize yardımcı olan Data Explorer.

Seqera Platform'un ücretsiz bir katmanı var, bu nedenle bu özelliklerin hemen hemen tamamını şu anda ücretsiz olarak kullanabilirsiniz. Ve kurumsal e-posta adresinizle kaydolursanız size Seqera Compute ile yüz dolarlık ücretsiz hesaplama kredisi bile vereceğiz. Son olarak, bir akademik program var, bu nedenle bir üniversitede çalışıyorsanız, fiyatlandırma sayfasına göz atın, oradaki formu bulun ve bize bildirin, sizi ücretsiz olarak Cloud Pro'ya yükselteceğiz.

## Topluluk Yardımı ve Etkinlikler

Tamam. İleriye dönük olarak. Nextflow ile ilgili herhangi bir desteğe ihtiyacınız olursa, community.seqera.io'ya göz atın. Gerçekten aktif ve sizi orada görmeyi ve farklı sorunlarınızı ve kullanım durumlarınızı tartışmayı umuyoruz ve belki şimdi diğer bazı insanlara bile yardımcı olabilirsiniz.

Ayrıca devam eden birçok etkinliğimiz var. nf-core ve Nextflow'dan gelen topluluk etkinliklerimiz var. Mart ayında çevrimiçi ve dağıtılmış bir nf-core hackathon'umuz var, geçen yıl dünyanın her yerinden sitelerle bin kişinin üzerinde katılımcı oldu. Bu nedenle lütfen yapabiliyorsanız bize katılın.

Ayrıca Nextflow Summit etkinliklerimiz var, biri Boston'da ve sonra Barcelona'da ve çevrimiçi bir etkinliğimiz var. İnsanların Nextflow'u gerçekten büyük, çılgın ve heyecan verici farklı şekillerde kullanmasını duyabileceğiniz harika konuşmalar. Ve ayrıca bunlarla ilişkili hackathonlar ve yüz yüze eğitimler var.

## Nextflow Podcast ve Blog

Nextflow ekosisteminde olup bitenlerden haberdar olmak istiyorsanız, seqera.io/blog'u mutlaka inceleyin.

Orada toplulukta çalışan insanlardan topluluk blog yazılarını ve ayrıca Seqera'dan Nextflow ve oluşturduğumuz diğer araçlara yönelik güncellemeler hakkında blog yazılarını duyabileceğiniz bir Nextflow bölümü var.

Ayrıca benim kişisel projem olan Nextflow Podcast'e bir tanıtım yapmak istiyorum. Spotify'da veya Apple Music'te veya YouTube'da inceleyin. Nextflow ile veya ilişkili teknolojilerle çalışan veya topluluktaki diğer insanlarla sohbet ettiğim yeni bölümler periyodik olarak yayınlıyoruz. Ve işlerin nasıl çalıştığına ve insanların ne yaptığına dair gerçek teknik, derinlemesine incelemeler yapıyoruz. Bu nedenle ilgileniyorsanız, onları inceleyin. Gerçekten eğlenceliler.

## Teşekkürler

Tamam, bir dizi teşekkür etmek istiyorum. Seqera'daki eğitim ekibi bu materyalden sorumludur. Bir kameranın önünde oturuyorum, ancak gerçekten tüm zor işler diğer insanlar tarafından yapıldı. Özellikle Hello Nextflow ve diğerleri için bu eğitim materyali kursunu yazan ve yenileyen Geraldine'e özel bir teşekkür. Ve ayrıca özellikle yeni Nextflow sözdizimi için sözdizimini güncellemede ve ayrıca kursların çoğunu kendisi yazmada gerçekten yardımcı olan Jon'a. Rike, Rob, Florian ve diğerleri gibi bilimsel geliştirme ekibindeki diğerlerinin üzerinde çalıştığımız materyale büyük katkıları oldu.

Ayrıca topluluktaki insanlara teşekkür etmek istiyorum. Örneğin, çok yeni olan yeni çeviriler, büyükelçi programındaki ve başka yerlerdeki insanlar tarafından büyük ölçüde etkilenmiştir. Ve gerçekten, eğitim materyalinin açık kaynak doğası, oldukça sık gelen ve bize gerçekten yardımcı olan pull request'lerin ve issue'ların olduğu anlamına gelir.

## Anket

Şimdi bitirdiğinize göre, henüz yapmadıysanız, lütfen hızlıca geri bildirim anketini yapın. training.nextflow.io web sitesinde Hello Nextflow bölümünün hemen altında.

Sadece beş soru. Gerçekten, gerçekten hızlı, ancak bu kabaca kaç kişinin eğitimi yaptığını takip etmemizi sağlıyor ve ayrıca bize eğitim materyalini nasıl geliştirebileceğimizi söyleyebilirsiniz. Tüm yanıtları gerçekten kontrol ediyoruz, bu nedenle oradaki girdilerinize gerçekten değer veriyoruz.

## Hoşça Kalın

Bir kez daha, bu kursa ve bu yolculuğa katıldığınız için çok çok teşekkürler. Eğitim materyalinde geliştirilmesi gerektiğini düşündüğünüz bir şey fark ettiyseniz bir GitHub issue veya Pull Request bırakın. Ve sizi başka bir Nextflow eğitim kursunda veya belki bir hackathon'da veya bir etkinlikte görmeyi gerçekten umuyorum. Tekrar teşekkürler.​
