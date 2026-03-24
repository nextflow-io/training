# Sonraki Adımlar - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=Lp8hS8RdaMwbp5j5&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım tam talimatlar için [kurs materyaline](../next_steps.md) dönün.

## Hoş Geldiniz

​

Tebrikler, başardınız!

Sona ulaştınız ve Hello Nextflow eğitim kursunu tamamladınız. Umarım keyif almışsınızdır. Sonuna kadar bizimle kaldığınız için çok teşekkür ederiz ve Nextflow'u öğrenmek için harcadığınız zaman ve çabayı gerçekten takdir ediyoruz. İşiniz için faydalı olacağını umuyoruz.

## training.nextflow.io'daki Diğer Kurslar

training.nextflow.io'yu ziyaret etmeye devam etmeyi unutmayın. Sürekli olarak yeni kısa kurslar ekliyoruz ve ayrıca zaten burada olan materyallerin çoğunu yeniliyoruz. Dolayısıyla bu Hello Nextflow eğitim kursu zamanla güncellenecektir.

Bu özellikle önemlidir çünkü Nextflow'da sözdizimini güncelliyoruz ve 2026'da oldukça fazla yeni özellik gelecek, dolayısıyla bu kurs 2027'de bir sonraki yaptığımızda biraz farklı görünecek ve hissettirecek.

Özellikle "Nextflow for Science" sayfasına dikkat çekmek istiyorum. Bunlar kısa kurslardır ve bu Hello Nextflow kursunu takip edecek şekilde tasarlanmıştır. Genomik, RNAseq veya her türlü farklı şey olsun, belirli farklı kullanım durumlarıyla Nextflow'un nasıl kullanılacağını gösterirler. Sürekli olarak daha fazla bilimsel kullanım durumu eklemeye çalışıyoruz.

Ayrıca Yan Görevler de var. Hello Nextflow gibi bir kurs geliştirdiğimizde, kapsayabileceğimiz çok şey var ve her şeyi kapsam içinde tutmak zor. Dolayısıyla insanlar için ilginç olduğunu düşündüğümüz, daha derinlemesine ele alınmayı hak eden belirli bir konu varsa, bunu bir Yan Görev'e koyuyoruz.

Gidin ve inceleyin; nf-test veya metadata ile farklı şeyler yapmak ve yaygın betik kalıpları gibi işinizle alakalı olabilecek farklı şeyler varsa, Yan Görevlere göz atın ve daha fazla bilgi edinmek için faydalı olup olmayacağına bakın.

Ayrıca nf-core hakkında bir kurs var. Umarım bu noktada projeye aşinasınızdır, ama değilseniz, gidin ve inceleyin. Farklı analiz türleri ve farklı veri türleri için neredeyse 150 farklı pipeline var; dolayısıyla ihtiyacınız olan veri analizi türü için kullanıma hazır bir pipeline olması tamamen mümkün.

Önemlisi, nf-core'da ayrıca bileşenler var; neredeyse 1700 farklı modül, farklı süreçler ve araçlar için sarmalayıcılar. nf-core ile gelen araçlarla bunları karıştırabilir ve Lego parçaları gibi kendi pipeline'ınızı oluşturabilirsiniz. Çok daha hızlı ve daha tekrarlanabilir.

## Seqera Platform

Nextflow ile kullanımınızı ölçeklendirirken Seqera Platform'u inceleyin; Nextflow'u çalıştırmanın en iyi yolu budur. Kendi altyapınızda, yani HPC veya AWS, Azure, Google Cloud, Oracle ve daha fazlasında çalıştırabilirsiniz. Hiçbir bilişim altyapısını yönetmek istemiyorsanız kendi Seqera Compute'umuzu da kullanabilirsiniz.

Seqera Platform, ortamı sizin için oluşturan Batch Forge gibi özelliklerle bu karmaşık bulut altyapılarının kurulumunu gerçekten basitleştirir. Ayrıca gözlemlenebilirlik, denetim günlüğü ve uyumluluk konusunda da gerçekten yardımcı olur.

Disk erişimini ve veri aktarımlarını optimize eden Fusion gibi teknolojilerle pipeline'ların nesnel olarak daha ucuz ve daha hızlı çalışmasını sağlar. Pipeline'larınızın yapılandırmasının mümkün olduğunca sıkı ayarlandığından emin olmak için pipeline optimizasyonu da mevcuttur.

Pipeline'ları çalıştırmanın dışında tamamen farklı özellikler de var. Etkileşimli analizler çalıştırabileceğiniz ve oluşturduğunuz herhangi bir özel Docker imajından ortamlar oluşturabileceğiniz Studios'umuz var. Ayrıca nerede olurlarsa olsunlar farklı dosya sistemlerinizi keşfetmenize yardımcı olan Data Explorer da mevcut.

Seqera Platform'un ücretsiz bir katmanı var; dolayısıyla bu özelliklerin neredeyse tamamını şu anda ücretsiz olarak kullanabilirsiniz. Kurumsal e-posta adresinizle kaydolursanız Seqera Compute ile yüz dolarlık ücretsiz işlem kredisi de veriyoruz. Son olarak, bir akademik program var; bir üniversitede çalışıyorsanız fiyatlandırma sayfasına göz atın, oradaki formu bulun ve bize bildirin, sizi Cloud Pro'ya ücretsiz olarak yükselteceğiz.

## Topluluk Yardımı ve Etkinlikler

Tamam. İleriye dönük olarak. Nextflow ile ilgili herhangi bir desteğe ihtiyacınız olursa community.seqera.io'ya göz atın. Gerçekten aktif bir topluluk; sizi orada görmeyi ve farklı sorunlarınızı ve kullanım durumlarınızı tartışmayı umuyoruz. Belki artık başkalarına bile yardımcı olabilirsiniz.

Ayrıca devam eden birçok etkinliğimiz var. nf-core ve Nextflow'dan gelen topluluk etkinliklerimiz var. Mart ayında çevrimiçi ve dağıtılmış bir nf-core hackathon'umuz var; geçen yıl dünya çapında sitelerden bine yakın kişi katıldı. Yapabilirseniz lütfen bize katılın.

Ayrıca Nextflow Summit etkinliklerimiz var; biri Boston'da, ardından Barcelona'da ve çevrimiçi bir etkinliğimiz var. İnsanların Nextflow'u gerçekten büyük, vahşi ve heyecan verici farklı şekillerde kullandığını duyabileceğiniz harika konuşmalar. Bunlarla ilişkili hackathon'lar ve yüz yüze eğitimler de var.

## Nextflow Podcast ve Blog

Nextflow ekosisteminde olup bitenlerden haberdar olmak istiyorsanız seqera.io/blog'u mutlaka kontrol edin.

Orada Nextflow için bir bölüm var; topluluktaki çalışan insanlardan topluluk blog yazılarını ve Nextflow ile ürettiğimiz diğer araçlara yönelik güncellemeler hakkında Seqera'dan blog yazılarını okuyabilirsiniz.

Ayrıca benim tutkum olan projem için de bir çağrı yapmak istiyorum: Nextflow Podcast. Bunu Spotify'da, Apple Music'te veya YouTube'da kontrol edin. Nextflow veya ilişkili teknolojilerle çalışan ya da topluluktaki insanlarla sohbet ettiğim yeni bölümler periyodik olarak yayınlıyoruz. Şeylerin nasıl çalıştığı ve insanların ne yaptığı hakkında gerçek teknik, derin dalışlar yapıyoruz. İlgileniyorsanız, kontrol edin. Gerçekten eğlenceli.

## Teşekkürler

Tamam, bir dizi teşekkür etmek istiyorum. Seqera'daki eğitim ekibi bu materyalden sorumlu. Ben bir kameranın önünde oturuyorum, ama gerçekten tüm zor iş o diğer insanlar tarafından yapıldı. Özellikle Hello Nextflow ve diğerleri için bu eğitim materyali kursunu yazan ve yenileyen Geraldine için özel bir çağrı. Ayrıca özellikle yeni Nextflow sözdizimi için sözdizimini güncellemede yardımcı olan ve kursların çoğunu kendisi yazan Jon için de. Rike, Rob, Florian ve diğerleri gibi bilimsel geliştirme ekibindeki kişiler de üzerinde çalıştığımız materyale büyük katkıda bulundu.

Ayrıca topluluktaki insanlara da teşekkür etmek istiyorum. Örneğin çok yeni olan yeni çeviriler, büyükelçi programındaki ve başka yerlerdeki insanlar tarafından büyük ölçüde etkilendi. Eğitim materyalinin açık kaynak doğası, oldukça sık gelen pull request'leri ve issue'ları beraberinde getiriyor; bu da gerçekten bize yardımcı oluyor.

## Anket

Şimdi bitirdiğinize göre, henüz yapmadıysanız lütfen hızlıca geri bildirim anketini doldurun. training.nextflow.io web sitesinde Hello Nextflow bölümünün hemen altında yer alıyor.

Sadece beş soru. Gerçekten çok hızlı; ama bu, yaklaşık olarak kaç kişinin eğitimi yaptığını takip etmemizi sağlıyor ve ayrıca eğitim materyalini nasıl geliştirebileceğimizi bize söyleyebiliyorsunuz. Tüm yanıtları gerçekten kontrol ediyoruz; dolayısıyla oradaki girdilerinizi gerçekten değerlendiriyoruz.

## Hoşça Kalın

Bir kez daha, bu kursa ve bu yolculuğa katıldığınız için çok çok teşekkürler. Eğitim materyalinde geliştirilmesi gerektiğini düşündüğünüz bir şey gördüyseniz bir GitHub issue veya Pull Request bırakın. Sizi başka bir Nextflow eğitim kursunda ya da belki bir hackathon'da veya bir etkinlikte görmeyi gerçekten umuyorum. Tekrar teşekkürler.​
