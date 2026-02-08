# Yönelim - Video Metni

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli not"

    Bu sayfa yalnızca metni göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../00_orientation.md) dönün.

## Hoş Geldiniz

Merhaba ve Hello Nextflow'a hoş geldiniz. Benim adım Phil Ewels. Nextflow'un arkasındaki şirket olan Seqera'da Açık Kaynak Yazılım Ürün Yöneticisiyim.

Bu kurs, Nextflow ile iş akışları oluşturmaya yönelik uygulamalı bir giriştir. Nextflow'a tamamen yeni başlayanlar ve kendi pipeline'larını geliştirmek isteyenler için tasarlanmıştır.

Örneklerin hepsi basit metin işlemedir, böylece alan uzmanlığına ihtiyaç duymadan Nextflow kavramlarına odaklanabilirsiniz, sadece biraz komut satırı aşinalığı yeterlidir.

Nextflow'un temellerini inceleyeceğiz: süreçler yazmak, bunları çok adımlı iş akışlarına bağlamak, yazılım bağımlılıklarını konteynerlerle yönetmek ve pipeline'ları farklı bilgi işlem ortamları için yapılandırmak. Sonunda, sıfırdan çalışan bir pipeline oluşturmuş olacaksınız.

Bu kurs pipeline _geliştirmeye_ odaklanır. Eğer sadece koda çok fazla dalmadan mevcut pipeline'ları _çalıştırmak_ istiyorsanız, size daha uygun olabilecek daha kısa bir "Nextflow Run" kursumuz var.

Burada temelleri kavradıktan sonra, bu kavramları gerçek bilimsel analizlere uygulayan devam kurslarımız da var. Size nf-core topluluğunun pipeline'larını ve en iyi uygulamalarını kullanmayı öğreteceğiz.

Sıkışıp kalırsanız community.seqera.io'ya gidin. Orada yalnızca eğitim sorularına ayrılmış bir bölümü olan aktif bir topluluk forumu var. İstediğiniz zaman kullanabilirsiniz, ancak ayrıca yardım etmek için özel olarak görevli kişilerle üç ayda bir eğitim haftaları düzenliyoruz. Dolayısıyla bunlardan biri sırasında eğitim alıyorsanız, kesinlikle çekinmeyin ve yardım isteyin.

Ayrıca yardım için Seqera AI'ya sormayı deneyebilirsiniz. Nextflow kodunu açıklamada ve hata ayıklamada size yardımcı olmada harikadır.

Nextflow'u ölçekte çalıştırmaya hazır olduğunuzda, bunu yapmanın en iyi yeri Seqera Platform'dur. Herhangi bir satıcı kilitlenmesi olmadan altyapınızda çalışır; pipeline başlatmadan gerçek zamanlı izlemeye, etkileşimli analiz ortamlarına kadar her şeyle. Ama şimdilik sadece temellere odaklanalım.

Pekala, hadi başlayalım.

## training.nextflow.io

Tamam. İlk belirtilmesi gereken şey, training.nextflow.io'daki tüm eğitim kurslarının çok etkileşimli olmasıdır. Fikir şu ki, eğitim materyalini ve talimatlarımı takip ediyorsunuz ve eğitim materyalini birlikte inceliyoruz. Bu yüzden iki şeye ihtiyacınız olacak: dizüstü bilgisayarınıza ve bu web sitesinin açık olmasına ihtiyacınız olacak. Ve hemen hemen bu kadar.

Bu, bunu kaydettiğimde bugün göründüğü gibi ana sayfadır. Farklı şeylere, arka plana ve sahip olduğumuz farklı kurslara bir genel bakış görebilirsiniz, bu liste her zaman büyüyor.

Yeni başlayanlar için Nextflow bizim bulunduğumuz yer. Burada iki kurs var, farklı bir kurs olan Nextflow Run ve bizim önemsediğimiz Hello Nextflow.

Ve ayrıca kenar çubuğundaki tüm farklı kursları görebilirsiniz. Hello Nextflow'a geçebilirim ve birlikte üzerinde çalışacağımız tüm farklı bölümleri görebilirim.

Burada dikkat edilmesi gereken birkaç önemli şey daha var. Öncelikle, eğitim materyali sürümlüdür, yani burada görebilirsiniz. 3.0 latest (en son) diyor, bu kaydettiğim sırada en son kararlı sürüm. Bu zaman içinde değişecek. Zaman içinde yeni kurslar yayınlıyoruz ve materyali güncelliyoruz. Yani 3.1 veya 3.2 ise çok fazla endişelenmeyin. 4.0 ise muhtemelen yeni bir video var ve belki de onu bulmalısınız çünkü muhtemelen önemli güncellemeler olacaktır.

Üstteki başka bir açılır menü de bu dil menüsü. Şimdi bu sürüm 3.0 için yepyeni. Daha önce İnsanlar tarafından elle yapılan çevrilmiş materyali aldık ve bunu bir LLM'ye aktardık ve LLM çevirisi kullanarak eğitim materyalinin farklı çevirilerini sürdürmek için bu tamamen yeni altyapıyı kurduk.

Yani artık burada tüm bu fantastik çevirilere sahibiz. Yani Korece dinlemek istiyorsanız, tüm web sitesini Korece yükleyebilirsiniz. Ve orada takip edebilirsiniz. Aynısı tüm bu diğer diller, Hintçe ve Almanca vb. için de geçerlidir. Ben İngilizce takip edeceğim. Bu, materyali yazdığımız birincil dil gibi.

Açık mod tercih ederseniz birkaç düğme daha. Karanlık mod yerine, burada üstte açık modda web sitesini takip edebilirsiniz.

Ve sonra da baktığımız her şey, açık kaynaklı olan nextflow-io/training adlı tek bir GitHub deposundadır. Ve bu düğmeye herhangi bir noktada tıklarsanız, GitHub deposuna gidecektir. Buna bir dakika içinde geri döneceğiz.

## GitHub Codespaces'i Kurmak

Tamam, artık bunu tarayıcı sekmesinde açtınız. Hadi Hello Nextflow'a gidelim ve tıklayalım. Giriş sayfasında bazı gereksinimleri, genel bakışı ve kabaca ele alacağımız şeylerin ders planını görebilirsiniz ve sonra başlangıca dalacağız.

Bu etkileşimli öğreticiyi yapmanın farklı yolları vardır. İsterseniz, bunu kendi bilgisayarınızda kendi Nextflow kurulumunuzla yerel olarak yapmaktan memnunsunuz. Ortam Seçenekleri'ne tıklarsak, bunu yerel Devcontainer'lar kullanarak veya ayrıca tüm yazılımı manuel kurulumla yerel olarak kurarak nasıl yapacağınıza dair daha fazla ayrıntı görebilirsiniz.

Bunu Seqera Studios ile düzgün çalışır hale getirmek için çalışıyoruz, bu da başka bir seçenek. Ancak şu anda en yaygın olan GitHub Codespaces kullanmaktır.

Codespaces, GitHub tarafından çalıştırılan uzak bir sunucuda bir sandbox ortamı kurar. Ve belirli bir kullanım miktarı için ücretsizdir, bu genellikle eğitim için iyidir. Ve size bir VS Code örneği, depodaki tüm dosyalara erişebileceğiniz, Nextflow'u ve her şeyi çalıştırabileceğiniz bir IDE kuracaktır. Ve Codespaces'i sizin için önceden yapılandırdık. Yani ihtiyacınız olan her şeye sahip.

Bunun güzelliği, bir Codespace kurmak için sadece bir tıklama olması. Herkes için aynı ve tüm ön koşulların zaten yüklendiğini biliyoruz, bu yüzden güzel ve hızlı.

Yani yapılacak ilk şey "Başlarken"e gitmek. _Codespaces'te Aç_ yazan bu düğmeyi arayın. Yeni bir sekmede açmak için komut + tıklama yapacağım ve bizi GitHub'a götürüyor.

İşte böyle görünüyor. Tüm seçenekleri burada sizin için ayarladık. İsterseniz, seçenekleri değiştir'e tıklayabilirsiniz. Burada yapabileceğiniz bazı şeyler var. Örneğin, bellek tükenmesi nedeniyle çökerse veya buna benzer bir şey varsa daha büyük bir örnek makine verebilirsiniz. Veya eğitim materyalinin belirli sürümlerini ayarlayın. Ama genellikle burada kurduğumuzla devam edebilirsiniz ve görebilirsiniz. Bu durumda, 3.0 sürümünü kullanıyor.

Bu yüzden yeni Codespace oluştur'a tıklayacağım. Ve bu beni içeri alıyor.

Ayrıca, orada devam edilecek Codespace yok diyor. Daha önce bir Codespace oluşturduysam, eğitim materyalindeki o düğmeye tekrar tıklamak beni aynı sayfaya götürür ve zaten çalışan tüm Codespace'leri listeler. Sonra hemen onlara geri dönüp kaldığınız yerden devam edebilirsiniz. Yani dizüstü bilgisayarınızı kapatsanız bile sorun değil.

Birkaç dakikalık hareketsizlikten sonra kendilerini otomatik olarak kapatırlar, ama sorun değil. Sadece yeniden başlatabilirsiniz.

Yeni bir Codespace başlattığınızda, bu sayfada böyle oturacak ve oldukça uzun süre yüklenecek. Yani şimdi kısa bir mola vermek için iyi bir zaman. Belki tuvalete gitmeyi unuttunuz veya başlamadan önce bir fincan çay istiyorsunuz? Bunu beklerken şimdi gidin, çünkü bir süre orada dönecek.

Hızlıca yüklenmesini beklerken, github.com/codespaces'e de gideceğim ve sadece şu anda çalışan tüm farklı Codespace'leri görebileceğiniz genel bakış sayfasını göstereceğim.

Burada nextflow-io/training için bir tane olduğunu görebilirsiniz. Değişiklik yok, çünkü henüz içinde bir şey yapmadım. Kullandığı kaynak miktarı ve şu anda kurulum yaptığını görebilirsiniz. Buraya gidebilir, bu küçük açılır menüye tıklayabilir ve sil'e tıklayabilirsiniz. Yani yanlışlıkla birden fazla Codespace kurarsanız ve bazılarını kullanmıyorsanız, eskileri silebilir ve temizleyebilirsiniz.

Son olarak, buna girmenin bir yolu daha. GitHub deposuna gidersek. Ve bu herhangi bir GitHub deposu için çalışır. Kod'a tıklayın. Depoyu yerel olarak klonlamak için komutlarınız var. Ve Codespaces adlı bir sekme var. Ve yine, yeni bir tane oluşturabilirsiniz ve zaten çalışan olanları görebilirsiniz.

Yani yine, Codespace'inizi nasıl oluşturduğunuzu unutursanız, her zaman bu şekilde geri dönebilirsiniz.

## VS Code Arayüzü

Tamam, oluşturucular bitti ve şimdi GitHub Codespaces'i yüklemeye başlıyor. Her zaman bu kadar uzun sürmez, bu yüzden endişelenmeyin. Sadece Codespace'i ilk kez oluşturduğunuzda. Zaten var olan bir Codespace'e geri atladığınızda çok daha hızlı.

Bu ilk sefereyse çok sabırsız olmayın, bize bir arayüz vermeye başlamasına rağmen henüz bitmedi.

Ama beklerken, size arayüzden VS Code'a biraz yabancıysanız bahsedeceğim.

İlk olarak, yapay zeka işleri için sohbet kenar çubuğu var, buna ihtiyacımız yok. Bu yüzden onu kapatacağım, ondan kurtulacağım ve biraz yer açacağım.

Solda, oluşturduğumuz çalışma alanı olan Git deposundaki tüm dosyaları bize gösteren dosya gezgini var. Not edin, bunlar yerel dosyalar değil. Bunların hepsi çalıştığımız uzak sunucuda. Yerel dosyaları sürükleyip bırakabilirsiniz ve benzeri şeyler yapabilirsiniz, ancak çoğunlukla bugün bunu düşünmeyeceğiz. Sadece tamamen uzaktan çalışıyoruz.

Bu kenar çubuğunda başka araçlar da var, örneğin arama. Yani bir depodaki tüm dosyaları tek seferde arayabilirsiniz. Ve eğitim deposunda geliştirme çalışması yapıyor olsaydık, Git ile kaynak kontrolü entegrasyonu ve hata ayıklama ve diğer şeyler yapabilirdik.

Diğer şeyler, burada yukarıda, sadece eğitim materyali için olan readme'nin bir önizlemesini yükleyen ana kod düzenleme penceresi var. Yani bu durumda markdown görüntülüyor, ama normalde bu bir kod düzenleyici olacak.

Ve bunun altında, tüm komutlarımızı çalıştıracağımız ve doğrudan Nextflow ile etkileşime gireceğimiz terminal var.

Codespace'teki her şey önceden yüklenmiş, bu nedenle Nextflow komutu zaten orada ve benzeri.

Tamam. Buraya geldiğinizde, yaklaşık olarak yapılmış olmalı. Şimdi Nextflow dil sunucusunu indirdiğini ve Nextflow uzantısı da dahil olmak üzere VS code'da bizim için bazı uzantılar kurduğunu görebilirsiniz, bu yararlı olacak. Bu yüzden onu kapatabilir ve README.md'yi kapatabilir.

Ve şimdi sol tarafta biraz daha fazla şey olduğunu görebilirsiniz. Burada biraz yakınlaştırdım, ama yakınlaştırmayı kaldırırsam, düğmelerden birinin Nextflow simgesiyle Nextflow yazdığını görebilirsiniz. ve bunun içinde projeyi keşfetmek için bazı güzel şeyler var ve şeyler, bunlara daha sonra geri döneceğiz.

Tamam. Bu panellerden herhangi birini kaybederseniz, sağ üstteki bu düğmeler gerçekten kullanışlıdır ve bunlar sadece şeyleri gösterir ve gizler. Bu, Explorer'ı gösterir ve gizler, alttaki terminali gösterir ve gizler. Ve benzeri.

Bunları oldukça fazla kullanacağım çünkü çok fazla yakınlaştırdım, bu yüzden ekranimdaki tüm metni görmenize yardımcı olmaya çalışıyorum ve bu yüzden terminalle tam ekrana gidebilmek ve sonra koda bakarken gizleyebilmek yararlı. Ama çoğu zaman tüm bunları aynı anda açık tutabilirsiniz.

Tamam, başka ne bakmak var? Çok fazla değil. Nextflow'un, dediğim gibi, kurulu olduğunu unutmayın. Bu yüzden "nextflow -version" yazabilirim ve kurulu sürümümüzün hangisi olduğunu söyleyerek gelmelidir.

Burada başka şeyler de yüklenmiş. Örneğin, her bölümün sonunda, web sitesinde bir dizi sınav sorusu var. Ve isterseniz quiz yazarak terminalde de bunları yapabilirsiniz.

Kullanacağım başka bazı klavye kısayolları da var, merak ederseniz diye. Örneğin, hemen önce Mac'imde cmd+K tuşuna bastım ve bu, önceki tüm çıktıdan kurtulmak için terminali temizledi. Bu yüzden bu, şeyleri temiz tutmak için güzel. Eğer bunu yaptığımı görürseniz, işte böyle yapıyorum.

Ve ayrıca terminale yeniyseniz, otomatik tamamlamak için tab kullanabileceğinizi unutmayın, bunu yolları otomatik tamamlamak için çok yapacağım.

Yani solda burada Hello Nextflow adlı bir klasör görebiliyorum, üzerinde çalışacağımız bu. "ls" yaparsam dosyaları listelemek için, "hel" yapabilirim, tab'a basarım, otomatik tamamlar. Ve bu, yolları tamamlamanın çok hızlı bir yolu.

## Sadece Hello Nextflow Klasörünü Açmak

Tamam. Bu harika. Yine de bu depoda çok şey var.

Web sitesi oluşturmak için tüm dosyalar var ve burada birden fazla farklı kurs var ve bunu bu kökten yapabilir ve sadece "Hello Nextflow" klasörüne tıklayabilirsiniz. Ama sadece buna odaklanmak güzel.

Bunu burada bir sürü tıklama ve proje dizini ayarlama ve benzeri şeylerle çalışma alanınız olarak ayarlayabilirsiniz. Ama en kolay yolu, VS Code'u başlatmak için CLI komutu olan code yazmak ve sonra "hello-nextflow".

Bu yeni bir tarayıcı sekmesi açacak ve eskisini kapatabilirsiniz. Ve tamamen aynı görünüyor. Ama şimdi bu alt dizinde olduğumuzu görebilirsiniz ve diğer tüm dosyalar görünmez ve daha temiz bir kurulumumuz var.

Burada ayrıca geçerli çalışma dizininin şimdi Hello Nextflow klasörü içinde olduğunu görebilirsiniz. Yani güzel ve temiz. Yanlış yerde olmaktan endişelenmemize gerek yok. Tamam.

## 2026 için Yeni Nextflow Sözdizimi

Bu noktada bahsetmem gereken özel bir şey var. Şu anda, 2026'nın başında, Nextflow'a farklı özellikler getirmeye başlıyoruz ve büyük yenilerden biri Nextflow içinde yeni bir dil sözdizimi ayrıştırıcısı.

Temel olarak Nextflow dosyalarınızı okuyan ve bunu çalışma zamanı için anlayan motor. Sözdiziminde bazı değişiklikler var ve doğru sözdizimi ayrıştırıcısı etkinleştirilmiş Nextflow kullanmanız gerçekten önemli.

Bunun için iki şeye ihtiyacınız var. Güncel bir Nextflow sürümüne ihtiyacınız var ve etkinleştirildiğinden emin olmanız gerekiyor.

Tekrar "nextflow -version" yaparsam, Codespaces'in 25.10.2 ile çalıştığını göreceksiniz ve 25.10, bu şeyleri kullanabilmek için minimum sürümdür.

Benim için henüz çıkmamış ama yakında çıkacak olan 26.04 kullanıyorsanız. O zaman bu, varsayılan olarak yeni sözdizimi ayrıştırıcısını çalıştıracak ve başka bir şey yapmanıza gerek yok.

Ama 25.10 çalıştırıyorsanız, katı sözdizimi ayrıştırıcısını veya v2 sözdizimi ayrıştırıcısını etkinleştirmeniz gerekir.

Bu bir ortam değişkeniyle yapılır. Codespaces'te zaten ayarlanmış, bu yüzden bir şey yapmanıza gerek yok. Ama yerel olarak çalıştırıyorsanız, bunu ayarlamanız gerekir ve "echo $NXF_SYNTAX_PARSER" yaparak bunu doğrulayabilirim ve v2 olarak ayarlanmalı.

Yani yerel olarak çalıştırıyorsanız, sadece "export NXF_SYNTAX_PARSER=v2" yapın. Basit olarak. Ama bunu yapmayı unutmayın. çünkü aksi takdirde ilerledikçe bazı garip tutarsızlıklar ve hatalar göreceksiniz.

Nextflow sürümü ve sözdizimi ayrıştırıcısı etrafındaki bu şeylerden herhangi biri hakkında hiç emin değilseniz, öncelikle, Codespaces'teyseniz endişelenmenize gerek olmadığını unutmayın. Her şey düzgün şekilde kurulmalı. Ama ikinci olarak, Nextflow eğitim materyaline giderseniz, aşağı inerseniz, sürüm gereksinimleri hakkında konuşun, burada keşif sürümlerinin etrafındaki yardım sayfasına götüren bir bağlantı var ve bu, her şeyi ayrıntılı olarak açıklıyor.

Bir dakikanız varsa bunu okumaya değer. çünkü Nextflow kullanmaya başladığınızda duyabileceğiniz bazı farklı terimlerin ne olduğunu açıklamaya yardımcı olur. DSL1, DSL2, sözdizimi ayrıştırıcı bir, sözdizimi ayrıştırıcı iki ve benzeri şeyler. Bu yüzden sadece bunun üzerinden bir kontrol yapmaya değer ve bu, az önce söylediklerimin bir kısmını tekrarlıyor.

Ayrıca daha önce Nextflow kodu yazdıysanız ve tazeleme için geri geliyorsanız gerçekten yararlı. Size değişen bazı şeyler söyler ve sizi Nextflow kodunuzu nasıl güncelleyeceğinizi anlatan Nextflow belgelerinin bölümlerine bağlar.

## Kurs Dosyaları

Tamam. Tanımamız gereken son şey sadece bu dizinde olan dosyaları görmek. Kenar çubuğuna ya bakabilirsiniz ya da genellikle eğitim materyalinde, tree komutunu kullanırız, -L, bakmak için seviye sayısıdır. İki diyeceğiz ve bunu tam ekran yaparsam, bu temelde orada kenar çubuğunda gördüğümüzün aynısını yansıtıyor, ancak noktayla başlayan gizli dosyaları hariç tutuyor.

Yani \*.nf dosyaları, Nextflow anlamına gelir. Yani bunlar Nextflow betik dosyalarıdır ve burada eğitim materyalinin her farklı bölümü için bir başlangıç dosyası var, açacağımız ve keşfedeceğimiz ve sonra düzenleyeceğimiz.

Bu dosyaları ilerledikçe değiştireceğiz ve böylece her bölümün sonunda, dosyalar bir sonraki bölümün başlangıcıyla hemen hemen aynı görünmelidir. Ama size bu farklı dosyaları veriyoruz, böylece her zaman taze başlayabilir ve sözdizimini karıştırmak konusunda çok fazla endişelenmezsiniz.

Kesinlikle işe yaraması gereken bir şeyle karşılaştırmanız gerekiyorsa. Çözümler klasörüne bakabilirsiniz ve bu, bölümlerin her biri için nihai durum gibidir, böylece yazdıklarınızı orada olanla karşılaştırabilirsiniz.

Bir veri dizini var. Bu, kursta örnek, girdi verisi olarak kullanacağımız sadece bir greetings.csv dosyasına sahip ve kursun ilerleyen bölümlerinde açıklayacağımız bir yapılandırma dosyası ve bazı parametreler gibi şeyler.

## Bitirme

Tamam, umarım artık her şey çalışıyor. Ekranınız benimki gibi görünüyor ve her şeye nasıl ulaşacağınızı ve tüm farklı dosyaların ne olduğunu anlıyorsunuz.

Başlarken sayfasının en altına kaydırırsanız, ne yaptığımı anladığımı söylemesi gereken küçük onay kutusu. Ortamım çalışıyor ve ayarlandınız, çalışma dizininizi "Hello Nextflow" klasörüne düzgün şekilde ayarladınız.

Bunların hepsini işaretlediyseniz ve yeşil görünüyorlarsa. Bir sonraki videoya ve bir sonraki bölüme geçebiliriz, bu birinci bölüm. Hello World. Birazdan görüşürüz.
