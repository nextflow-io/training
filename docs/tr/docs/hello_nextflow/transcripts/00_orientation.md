# Oryantasyon - Video Metni

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli not"

    Bu sayfa yalnızca metni göstermektedir. Tam adım adım talimatlar için [kurs materyaline](../00_orientation.md) geri dönün.

## Hoş Geldiniz

Merhaba, Hello Nextflow'a hoş geldiniz. Benim adım Phil Ewels. Seqera'da Açık Kaynak Ürün Müdürüyüm ve bugün bu ilk Nextflow eğitim kursunda size rehberlik etmekten büyük mutluluk duyuyorum.

Nextflow'un temellerini inceleyeceğiz, pipeline'ları nasıl yazacağınızı, çalıştıracağınızı ve yapılandıracağınızı açıklayacağız.

Ve kendi basit çok adımlı pipeline'ınızı oluşturacaksınız. Operatörler ve channel fabrikaları gibi terimleri ele alacağız ve kursun sonunda kendi biyoinformatik pipeline'larınızı oluşturmaya başlamak için hazır olacaksınız.

Herhangi bir sorunuz varsa, lütfen community.seqera.io üzerinden bize ulaşın. Gerçekten aktif bir Nextflow topluluğumuz var, eğitime ayrılmış bir bölüm var, bu yüzden nerede takıldığınızı bize bildirin ve birisi yardımcı olabilecektir.

Tamam. Hadi başlayalım.

## Eğitim Web Sitesi

Nextflow kursları için tüm eğitim materyalleri training.nextflow.io adresinde bulunmaktadır. Web tarayıcınızda oraya gidebilirsiniz. Şimdi onu açın ve birlikte göz atalım.

Bunu 2.1.1 sürümüyle çalıştıracağım. Burada küçük güncellemeler ve düzeltmeler yayınlıyoruz, bu yüzden biraz farklıysa endişelenmeyin, ancak materyal çok fazla değiştiyse, konuşacağım materyalin tam sürümünü seçmek için en üstteki bu sürüm seçiciyi her zaman kullanabilirsiniz.

Daha çok açık mod insanıysanız, web sitesinin temasını buradan değiştirebilirsiniz.

Çevirileri burada görün, ancak kayıt sırasında bu yeni materyali kapsayan gerçekten yalnızca İngilizce var.

Ve GitHub'da eğitim web sitesi ve üzerinde çalışacağımız her şey için tüm kaynak kodunu görün.

Buradaki ana sayfa, sahip olduğumuz tüm farklı eğitim materyali kurslarını listeler. Aşağı kaydırdığımda, burada yapacağımız Hello Nextflow kursuyla birlikte Yeni Başlayanlar için Nextflow'u göreceğiz. Benzer şekilde çalışan sahip olduğumuz diğer tüm kursları da görebilirsiniz.

## Ortam Kurulumu

Aslında en üstteki bu ilkini kullanarak başlayacağım, bu tüm eğitim kursları için ortak ve özellikle ortamımızı kurmakla ilgili.

Tıkladığımda, beni bu bölüme götürüyor ve yerel olarak geliştirme için talimatları görebiliyoruz. Kendi dizüstü bilgisayarınızı kendi VS Code kopyanız ve kendi yazılım kurulumlarınızla kullanmak istiyorsanız veya çoğu insanın yapmasını beklediğimiz şey, GitHub Codespaces adlı bir şey kullanmaktır.

Codespaces, GitHub tarafından sağlanan, bulutta bir web sunucusu çalıştıran bir hizmettir ve ona bağlanabilirsiniz. Bu sunucuda VS code yüklü, bunu web tarayıcınızda çalıştırabilirsiniz veya isterseniz yerel VS code kurulumunuza bağlayabilirsiniz. Tüm hesaplama, tüm dosyalar, tüm düzenleme uzaktan gerçekleşir, bu da ihtiyacınız olan tüm yazılımın önceden yüklenmiş geldiği ve herkes için aynı olduğu anlamına gelir.

## GitHub Codespace Oluşturma

İhtiyacımız olan her şeyle codespace'i oluşturmak için dokümantasyon materyalinde "GitHub Codespaces'te Aç" yazan düğmeleri arayın. Şimdi ona tıklayacağım, yeni bir sekmede açacağım. Ve bu web sayfasıyla karşılaşıyorum. Şimdi bunun nextflow-io training ile ayarlanmış olduğunu görebilirsiniz.

Sadece yeni codespace oluştur'a tıklayabilirim. Ancak aslında Nextflow eğitimi için iki yerine dört CPU'lu biraz daha büyük bir makine kullanmayı öneriyoruz. Materyalin hangi sürümünü kullandığını değiştirebilirsiniz. Yani bu, bağlantıyı takip ettiğim dokümanların sürümü olduğu için 2.1.1'e varsayılan oluyor. Ancak istersem deponun belirli bir dalına da ayarlayabilirim.

Şimdi codespace oluştur'a tıklayacağım. Ve benim için ortamı kurmaya başlayacak.

## Codespace oluşturma

Şimdi, bunu ilk kez yaptığınızda çok uzun sürecek, bu yüzden şimdi gidip bir fincan çay almanın tam zamanı. Kendinizi rahat ettirin, yanınızda oturan kişiyle sohbet edin.

İlgileniyorsanız, kurulum günlüklerini görmek için buraya codespace oluşturuluyor'a tıklayabilirsiniz. Ve burada ihtiyacım olan her şeyi içeren bir Docker imajını çektiğini ve ortamı yapılandırdığını görebilirsiniz.

Şimdi, bu şekilde sadece ilk kez bir codespace oluşturduğunuzda beklemeniz gerekir. github.com/codespaces adresine giderseniz, açık olan tüm farklı Codespace'leri göreceksiniz. İşte yeni oluşturduğum. Bir dahaki sefere bunu yaptığınızda, buraya gidebilir ve önceki codespace'i seçebilir ve doğrudan içine atlayabilirsiniz. Ve mevcut ortamı ısıtmak çok, çok daha hızlı bir işlemdir. Bu ayrıca VS Code'a ve dosyalara yaptığınız tüm değişiklikleri de tutacak, böylece ayrılıp geri dönerseniz ilerlemenizi kaybetmezsiniz.

Diğer işlemleri yapmak için buradaki üç noktaya tıklayabilirsiniz. Örneğin, iki CPU ile yapılandırdıysanız ve şimdi dört istiyorsanız, makine türünü değiştirebilirsiniz. Veya sıfırdan ve temiz başlamak istiyorsanız, codespace'i silebilirsiniz.

## VS Code'a Giriş

Tamam, Codespaces ortamımı kurmayı bitirdi ve şimdi web tarayıcısında VS Code ile karşılaşıyorum.

VS code'a alışkınsanız. Bu çok tanıdık gelecek, daha önce kullanmadıysanız, oldukça basittir. Sayfanın farkında olmanız gereken birkaç farklı bölümü var.

Burada solda yan çubuğumuz var. Eğitim deposundan GitHub deposundaki tüm farklı dosyalarla Explorer'ın kurulduğunu görebilirsiniz.

Soldaki bu düğmelerde farklı araçlar olabilir. Yan çubukta. Projedeki tüm dosyaları arayabilirim. Git ile çalışabilirim, GitHub ile çalışabilirim, bunun gibi tüm farklı şeyler.

En üstte burada ana menü var. Dosya gezgini çoğunlukla burada sahip olacağımız şey ve bu dosyalardan herhangi birine sağ tıklayabilir ve beklediğiniz normal şeyleri yapabilirsiniz. Kes kopyala gibi bazı uyarıları tıklamanız gerekebilir ve yerel makinenize de indirebilirsiniz.

Codespace yüklendiğinde, bize bu ana alanda markdown dosyasının bir önizlemesini verir. Bu, github.com'da oluşturulanla aynı. Bunu kapatabilir ve Readme dosyasına çift tıklarsam, kod düzenleyicide kod olarak açıldığını göreceksiniz ve diğer herhangi bir dosyada olduğu gibi, bu kodu doğrudan düzenleyebiliriz.

Son olarak burada altta terminal penceremiz var. Oluşturulurken günlüklere bakıyordum, bu yüzden gösterdiği mevcut şey bu. Ayrıca yeni bir terminal oturumu başlatmak için bu artı düğmesine basabilirim. Bu benim makinemde çalışmıyor. Unutmayın, bu bulutta çalışıyor ve ikinin derinliğinde tree üç yaparsam, solda olan tüm aynı dosyaları burada göreceksiniz.

## Sadece "hello-nextflow" dosyalarını gösterme

Bu GitHub deposu, sadece yaptığımız değil, tüm farklı eğitim setlerini içerir. İsterseniz, sadece Hello Nextflow klasörüne odaklanabilirsiniz. Bunu biraz temizlemenin bir yolu, menü dosyasına gitmek ve ardından çalışma alanına klasör eklemektir.

Bunu tıklıyoruz, training'e gidiyoruz. Hello nextflow, ve ekle'ye tıklayın. Ekranınızı yenileyecek. Ve sonra Explorer'da, şimdi iki farklı çalışma alanımız var, daha önce training için sahip olduğumuz ve sadece Hello Nextflow ile bir tane.

İsterseniz, training'e sağ tıklayabilir ve çalışma alanından klasörü kaldır'a tıklayarak yan çubuktan tamamen kaldırabilirsiniz.

Şimdi yan çubukta sadece bu belirli eğitim kursu için dosyalarımız var. O uyarıyı gizleyebilir ve şimdi burada terminalde de aynı şeyi yapabilir ve dizin değiştirmek için CD yapabilirim. Hello, Nextflow. Ve yine, yan çubukta olan aynı dosyalara burada sahibiz.

## Hello Nextflow: dosyalar

Hello Nextflow kursu için bu dosyalara bakıyoruz.

Nextflow için olan bir dizi .nf dosyamız var ve eğitim kursunun her bölümü için bu dosyalardan biri var. Bu dosyalar üzerinde çalışacağız ve alıştırmalarda onları değiştireceğiz.

Ayrıca sadece bu ortamda Nextflow'u çalıştırmak için temel yapılandırma ayarlarına sahip bir nextflow.config dosyamız var, bu noktada gerçekten endişelenmenize gerek yok. Veri işlemek için kullanacağımız bir greetings.csv dosyası, bu kursun bir sonraki bölümünde tanıtılacak ve altıncı bölümde kullanılacak bir test-params.json dosyası, şimdilik göz ardı edebilirsiniz.

Bu Nextflow dosyaları sadece her alıştırmanın başlangıcıdır. Bitirdiklerinde nasıl görünmeleri gerektiğini görmek isterseniz, bir solutions dizinine gidebilirsiniz ve eğitim kursunun her bölümü için cevaplar vardır, böylece hedeflediğiniz şeyin çalışan bir sürümünü görebilirsiniz.

## Terminal açma

Herhangi bir noktada terminali kapatırsanız ve nasıl geri döneceğinizi hatırlayamazsanız, endişelenmeyin. En üstteki sağdaki bu düğmeler, çalışma alanında farklı panelleri açar ve kapatır. Bu yüzden alt panel için buna tıklayın ve yeniden görünecektir. Ve burada terminal seçildiğinden emin olun. Terminali tam ekran yapmak için buradaki bu düğmeye, sağ taraftaki terminaldeki oka da tıklayabilirsiniz.

Bunu çok sık yaptığımı göreceksiniz çünkü metni okuyabilmeniz için VS Code'u yakınlaştırdım. Ekran boyutunuza bağlı olarak, bunu yapmanız gerekebilir veya gerekmeyebilir. Aynı şey yan paneli küçültmek için de geçerlidir.

Tamam. Ortam için bu kadar yeterli sanırım. Başlamaya hazırız. Bir sonraki videoda birinci bölüm için benimle buluşun.

[Sonraki video metni :octicons-arrow-right-24:](01_hello_world.md)
