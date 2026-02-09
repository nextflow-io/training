# Yönlendirme - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli not"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../00_orientation.md) geri dönün.

## Hoş Geldiniz

Merhaba ve Hello Nextflow kursuna hoş geldiniz. Benim adım Phil Ewels. Nextflow'un arkasındaki şirket olan Seqera'da Açık Kaynak Yazılım Ürün Müdürüyüm.

Bu kurs, Nextflow ile iş akışları oluşturmaya yönelik uygulamalı bir giriş niteliğindedir. Nextflow'a tamamen yeni olan ve kendi boru hatlarını geliştirmek isteyen kişiler için tasarlanmıştır.

Örneklerin tümü basit metin işleme içerir, böylece alan uzmanlığına ihtiyaç duymadan Nextflow kavramlarına odaklanabilirsiniz, sadece biraz komut satırı aşinalığı yeterlidir.

Nextflow'un temellerini ele alacağız: süreçler yazmak, bunları çok adımlı iş akışlarına bağlamak, konteynerlerle yazılım bağımlılıklarını yönetmek ve farklı bilgi işlem ortamları için boru hatlarını yapılandırmak. Sonunda, sıfırdan çalışan bir boru hattı oluşturmuş olacaksınız.

Bu kurs boru hatları _geliştirmeye_ odaklanır. Koda fazla dalmadan mevcut boru hatlarını _çalıştırmak_ istiyorsanız, size daha uygun olabilecek daha kısa bir "Nextflow Run" kursumuz var.

Burada temelleri öğrendikten sonra, bu kavramları gerçek bilimsel analize uygulayan devam kurslarımız da var. Size nf-core topluluğunun boru hatlarını ve en iyi uygulamalarını nasıl kullanacağınızı öğreteceğiz.

Takılırsanız community.seqera.io adresine gidin. Orada sadece eğitim sorularına ayrılmış bir bölümü olan aktif bir topluluk forumu var. İstediğiniz zaman kullanabilirsiniz, ancak ayrıca özellikle yardım etmek için hazır kişilerle üç ayda bir eğitim haftaları düzenliyoruz. Dolayısıyla eğitimi bunlardan biri sırasında yapıyorsanız, kesinlikle çekinmeyin ve yardım isteyin.

Ayrıca yardım için Seqera AI'ya sormayı deneyebilirsiniz. Nextflow kodunu açıklamada ve hata ayıklamada size yardımcı olmada harikadır.

Nextflow'u ölçekte çalıştırmaya hazır olduğunuzda, bunu yapmanın en iyi yeri Seqera Platform'dur. Herhangi bir satıcı kilitlenmesi olmadan altyapınızda çalışır, boru hattı başlatmadan gerçek zamanlı izlemeye, etkileşimli analiz ortamlarına kadar her şeyle. Ama şimdilik sadece temellere odaklanalım.

Pekala, hadi başlayalım.

## training.nextflow.io

Tamam. İlk dikkat edilmesi gereken nokta, training.nextflow.io'daki tüm eğitim kurslarının çok etkileşimli olmasıdır. Fikir şu ki, eğitim materyalini ve talimatlarımı takip ediyorsunuz ve eğitim materyalini birlikte inceliyoruz. Yani iki şeye ihtiyacınız olacak: dizüstü bilgisayarınıza ve bu web sitesinin açık olmasına ihtiyacınız olacak. Ve hemen hemen bu kadar.

Yani bu, bunu kaydettiğimde göründüğü şekliyle ana sayfa. Farklı şeylere, arka plana ve sahip olduğumuz farklı kurslara genel bir bakış görebilirsiniz, liste her zaman büyüyor.

Yeni başlayanlar için Nextflow, bulunduğumuz yer. Burada iki kurs var, farklı bir kurs olan Nextflow Run ve önemsediğimiz Hello Nextflow.

Ayrıca kenar çubuğunda tüm farklı kursları görebilirsiniz. Hello Nextflow'a atlayabilirim ve birlikte üzerinde çalışacağımız tüm farklı bölümleri görebiliriz.

Burada dikkat edilmesi gereken birkaç önemli şey daha var. İlk olarak, eğitim materyali sürümlendirilmiştir, yani burada görebilirsiniz. 3.0 en son diyor, ki bu benim kaydettiğim sırada en son kararlı sürüm. Bu zamanla değişecek. Yeni kurslar yayınlıyoruz ve materyali zamanla güncelliyoruz. Yani 3.1 veya 3.2 ise, fazla endişelenmeyin. 4.0 ise, muhtemelen yeni bir video var ve belki onu bulmalısınız çünkü muhtemelen önemli güncellemeler olacaktır.

Üstteki bir başka açılır menü de bu dil menüsü. Şimdi bu sürüm 3.0 için yepyeni. Daha önce İnsanlar tarafından elle yapılan çevrilmiş materyali aldık ve bunu bir LLM'ye geçirdik ve LLM çevirisi kullanarak eğitim materyalinin farklı çevirilerini sürdürmek için bu yeni altyapıyı kurduk.

Yani şimdi burada tüm bu harika çeviriler var. Yani Korece dinlemek isterseniz, tüm web sitesini Korece yükleyebilirsiniz. Ve orada takip edebilirsiniz. Tüm bu diğer diller için de aynı, Hintçe ve Almanca ve benzeri. Ben İngilizce takip edeceğim. Bu, materyali yazdığımız birincil dil gibi.

Açık mod tercih ederseniz birkaç başka düğme daha var. Karanlık mod yerine, burada üstte açık modda web sitesini takip edebilirsiniz.

Ve sonra baktığımız her şey, açık kaynak olan nextflow-io/training adlı tek bir GitHub deposundadır. Ve herhangi bir noktada bu düğmeye tıklarsanız, GitHub deposuna gidecektir. Buna bir dakika içinde geri döneceğiz.

## GitHub Codespaces Kurulumu

Tamam, şimdi bunu tarayıcı sekmesinde açtınız. Hadi Hello Nextflow'a gidelim ve tıklayalım. Giriş sayfasında, gereksinimlerden bazılarını, genel bakışı ve kabaca neyi kapsayacağımızın ders planını görebilirsiniz, sonra başlangıca dalacağız.

Bu etkileşimli öğreticiyi yapmanın farklı yolları var. İsterseniz, bunu kendi bilgisayarınızda kendi Nextflow kurulumunuzla yerel olarak yapmaktan memnunsunuz. Ortam Seçeneklerine tıklarsak, bunu yerel Devcontainer'lar kullanarak veya manuel kurulumla tüm yazılımı yerel olarak kurarak nasıl yapacağınıza dair daha fazla ayrıntı görebilirsiniz.

Bunu Seqera Studios ile güzel bir şekilde çalıştırmak için çalışıyoruz, yani bu başka bir seçenek. Ama şu anda en yaygın olanı GitHub Codespaces kullanmak.

Codespaces, GitHub tarafından çalıştırılan uzak bir sunucuda bir sandbox ortamı kurar. Ve belirli bir kullanım miktarı için ücretsizdir, bu genellikle eğitim için yeterlidir. Ve size bir VS Code örneği, depodan tüm dosyalara erişebileceğiniz, Nextflow'u ve her şeyi çalıştırabileceğiniz bir IDE kuracaktır. Ve Codespaces'i sizin için önceden yapılandırdık. Yani ihtiyacınız olan her şeye sahip.

Bunun güzelliği, bir Codespace kurmak için sadece bir tıklama olması. Herkes için aynı ve tüm ön koşulların zaten yüklü olduğunu biliyoruz, yani güzel ve hızlı.

Yani yapılacak ilk şey "Başlarken"e gitmek. _Codespaces'te Aç_ yazan bu düğmeyi arayın. Yeni bir sekmede açmak için komut + tıklayacağım ve bizi GitHub'a götürüyor.

İşte böyle görünüyor. Tüm seçenekleri burada sizin için ayarladığımızı görebilirsiniz. İsterseniz, seçenekleri değiştir'e tıklayabilirsiniz. Burada yapabileceğiniz bazı şeyler var. Örneğin, bellek bitmesi nedeniyle çökerse veya buna benzer bir şey olursa daha büyük bir örnek makine verebilirsiniz. Veya eğitim materyalinin belirli sürümlerini ayarlayın. Ama genellikle burada kurduğumuz şeyle gidebilirsiniz ve görebilirsiniz. Bu durumda 3.0 sürümünü kullanıyor.

Yani yeni Codespace oluştur'a tıklayacağım. Ve bu beni içeri alıyor.

Ayrıca dikkat edin, orada devam edecek Codespace yok diyor. Daha önce bir Codespace oluşturduysam, eğitim materyalindeki o düğmeye tekrar tıklamak beni aynı sayfaya götürecek ve zaten çalışan tüm Codespace'leri listeleyecek. Sonra doğrudan onlara geri atlayabilir ve kaldığınız yerden devam edebilirsiniz. Yani dizüstü bilgisayarınızı kapattıysanız sorun değil.

Birkaç dakika hareketsizlikten sonra kendilerini otomatik olarak kapatırlar, ama sorun değil. Sadece yeniden başlatabilirsiniz.

Yeni bir Codespace başlattığınızda, bu sayfada böyle oturacak ve oldukça uzun süre yüklenecek. Yani şimdi kısa bir mola vermek için iyi bir zaman. Belki tuvalete gitmeyi unuttunuz veya başlamadan önce bir fincan çay istiyorsunuz? Bunu beklerken şimdi gidin, çünkü bir süre orada dönecek.

Yüklenmeyi beklerken hızlıca, github.com/codespaces'e de gideceğim ve şu anda çalışan tüm farklı Codespace'leri görebileceğiniz genel bakış sayfasını göstereceğim.

Burada nextflow-io/training için bir tane olduğunu görebilirsiniz. Değişiklik yok, çünkü henüz içinde hiçbir şey yapmadım. Kullandığı kaynak miktarı ve şu anda kurulum yaptığını görebilirsiniz. Buraya gidebilir, bu küçük açılır menüye tıklayabilir ve sil'e tıklayabilirsiniz. Yani yanlışlıkla birden fazla Codespace kurduysanız ve bazılarını kullanmıyorsanız, eskileri silebilir ve temizleyebilirsiniz.

Son olarak, buna girmenin bir yolu daha. GitHub deposuna gidersek. Ve bu herhangi bir GitHub deposu için çalışır. Kod'a tıklayın. Depoyu yerel olarak klonlamak için komutlarınız olabilir. Ve Codespaces adlı bir sekme var. Ve yine, yeni bir tane oluşturabilirsiniz ve zaten çalışan olanları görebilirsiniz.

Yani yine, Codespace'inizi nasıl oluşturduğunuzu unutursanız, her zaman bu şekilde geri dönebilirsiniz.

## VS Code Arayüzü

Tamam, oluşturucular bitti ve şimdi GitHub Codespaces'i yüklemeye başlıyor. Her zaman bu kadar uzun sürmez, yani endişelenmeyin. Sadece ilk kez Codespace oluşturduğunuzda. Zaten var olan birine geri dönerseniz, çok daha hızlıdır.

Bu ilk seferse çok sabırsız olmayın, henüz bitmedi, bize bir arayüz vermeye başlamasına rağmen.

Ama son şeylerin kurulmasını beklerken, VS Code'a biraz aşina değilseniz size arayüzü anlatacağım.

İlk olarak, ihtiyacımız olmayan AI şeyleri için sohbet kenar çubuğu var. Yani onu kapatacağım, ondan kurtulacağım ve biraz yer açacağım.

Solda, oluşturduğumuz çalışma alanı olan Git deposundaki tüm dosyaları bize gösteren dosya gezgini var. Bunların yerel dosyalar olmadığını unutmayın. Bunların hepsi çalıştığımız uzak sunucuda. Yerel dosyaları sürükleyip bırakabilirsiniz ve benzeri şeyler, ama çoğunlukla bugün bunu düşünmeyeceğiz. Sadece tamamen uzaktan çalışıyoruz.

Bu kenar çubuğunda başka araçlar da var, örneğin arama. Yani bir depodaki tüm dosyaları tek seferde arayabilirsiniz. Ve eğitim deposunda geliştirme çalışması yapıyor olsaydık, Git ile kaynak kontrolü entegrasyonu ve hata ayıklama ve diğer şeyler yapabilirdik.

Diğer şeyler, burada yukarıda sadece readme'nin önizlemesini yükleyen ana tür kod düzenleme penceresi var, bu eğitim materyali için. Yani bu durumda markdown görüntülüyor, ama normalde bu bir kod düzenleyici olacak.

Ve bunun altında tüm komutlarımızı çalıştıracağımız ve doğrudan Nextflow ile etkileşime gireceğimiz terminal var.

Codespace'teki her şey önceden yüklenmiş, yani Nextflow komutu zaten orada ve benzeri.

Tamam. Buraya geldiğinizde, hemen hemen bitmiş olmalı. Şimdi Nextflow dil sunucusunu indirdiğini ve VS code'da bizim için Nextflow uzantısı da dahil olmak üzere bazı uzantılar kurduğunu görebilirsiniz, bu yararlı olacak. Yani onu kapatabilir ve README.md'yi kapatabilir.

Ve şimdi sol tarafta biraz daha fazla şey olduğunu görebilirsiniz. Burada biraz yakınlaştırdım, ama uzaklaştırırsam düğmelerden birinin Nextflow simgesiyle Nextflow yazdığını görebilirsiniz. ve bunun içinde projeyi keşfetmek ve daha sonra geri döneceğimiz şeyler için bazı güzel şeyler var.

Tamam. Bu panellerden herhangi birini kaybederseniz, sağ üstteki bu düğmeler gerçekten yararlıdır ve bunlar sadece şeyleri gösterir ve gizler. Yani bu Gezgin'i gösterir ve gizler, alttaki terminali gösterir ve gizler. Ve benzeri.

Bunları oldukça fazla kullanacağım çünkü çok yakınlaştırdım, bu yüzden ekranımdaki tüm metni görmenize yardımcı olmaya çalışıyorum ve bu yüzden terminalle tam ekrana geçebilmek ve sonra koda bakarken gizleyebilmek yararlı. Ama çoğu zaman tüm bu şeyleri aynı anda açık tutabilirsiniz.

Tamam, başka ne var? Fazla bir şey yok. Nextflow'un, dediğim gibi, yüklü olduğunu unutmayın. Yani "nextflow -version" yazabilirim ve hangi sürümü yüklediğimizi söyleyerek gelmelidir.

Burada yüklü başka şeyler de var. Örneğin, her bölümün sonunda, web sitesinde bir dizi sınav sorusu var. Ve isterseniz bunları terminalde de quiz yazarak yapabilirsiniz.

Kullanacağım başka bazı klavye kısayolları var, merak ediyorsanız diye. Örneğin, az önce Mac'imde cmd+K'ye bastım ve bu, önceki tüm çıktıdan kurtulmak için terminali temizledi. Yani bu, işleri temiz tutmak için güzel. Bunu yaptığımı görürseniz işte böyle yapıyorum.

Ayrıca terminalde yeniyseniz, otomatik tamamlamak için sekme kullanabileceğinizi unutmayın, bunu yolları otomatik tamamlamak için çok yapacağım.

Yani burada solda Hello Nextflow adlı bir klasör görebiliyorum, üzerinde çalışacağımız şey bu. "ls" yaparsam dosyaları listelemek için, "hel" yapabilirim, sekmeye basarım, otomatik tamamlar. Ve bu, yolları tamamlamanın çok hızlı bir yolu.

## Sadece Hello Nextflow Klasörünü Açma

Tamam. Bu harika. Yine de bu depoda çok şey var.

Web sitesini oluşturmak için tüm dosyalar var ve burada birden fazla farklı kurs var ve bunu bu rotadan yapabilir ve sadece "Hello Nextflow" klasörüne tıklayabilirsiniz. Ama aslında sadece buna odaklanmak güzel.

Bunu burada bir sürü tıklamayla ve bir proje dizini ayarlayarak ve benzeri şeylerle çalışma alanınız olarak ayarlayabilirsiniz. Ama en kolay yol, VS Code'u başlatmak için CLI komutu olan code yazmak ve sonra "hello-nextflow".

Bu yeni bir tarayıcı sekmesi açacak ve eskisini kapatabilirsiniz. Ve tamamen aynı görünüyor. Ama şimdi bu alt dizinde olduğumuzu görebilirsiniz ve diğer tüm dosyalar görünmez ve daha temiz bir kurulumumuz var.

Burada ayrıca mevcut çalışma dizininin şimdi Hello Nextflow klasörü içinde olduğunu görebilirsiniz. Yani güzel ve temiz. Yanlış yerde olmak konusunda endişelenmemize gerek yok. Tamam.

## 2026 için Yeni Nextflow Sözdizimi

Bu noktada bahsetmem gereken özel bir şey var. Şu anda, 2026'nın başında, Nextflow'a farklı özellikler getirmeye başlıyoruz ve büyük yenilerden biri Nextflow içindeki yeni bir dil sözdizimi ayrıştırıcısı.

Temel olarak Nextflow dosyalarınızı okuyan ve bunu çalışma zamanı için anlayan motor. Sözdiziminde bazı değişiklikler var ve Nextflow'u doğru sözdizimi ayrıştırıcısı etkinleştirilmiş olarak kullanmanız gerçekten önemli.

Bunun için iki şeye ihtiyacınız var. Güncel bir Nextflow sürümüne ihtiyacınız var ve etkinleştirildiğinden emin olmanız gerekiyor.

"nextflow -version" tekrar yaparsam, Codespaces'in 25.10.2 ile çalıştığını göreceksiniz ve 25.10, bu şeyleri kullanabilmek için minimum sürümdür.

Benim için henüz çıkmamış ama yakında çıkacak olan 26.04 kullanıyorsanız. O zaman bu varsayılan olarak yeni sözdizimi ayrıştırıcısını çalıştıracak ve başka bir şey yapmanıza gerek yok.

Ama 25.10 çalıştırıyorsanız, katı sözdizimi ayrıştırıcısını veya v2 sözdizimi ayrıştırıcısını etkinleştirmeniz gerekir.

Bu bir ortam değişkeniyle yapılır. Codespaces'te zaten ayarlanmış, yani bir şey yapmanıza gerek yok. Ama yerel olarak çalıştırıyorsanız, bunu ayarlamanız gerekir ve bunu "echo $NXF_SYNTAX_PARSER" yaparak doğrulayabilirim ve v2 olarak ayarlanmalıdır.

Yani yerel olarak çalıştırıyorsanız, sadece "export NXF_SYNTAX_PARSER=v2" yapın. Bu kadar basit. Ama bunu yapmayı unutmayın. çünkü aksi takdirde ilerledikçe bazı garip tutarsızlıklar ve hatalar göreceksiniz.

Nextflow sürümü ve sözdizimi ayrıştırıcısı etrafındaki bu şeylerden herhangi biri hakkında hiç emin değilseniz, ilk olarak, Codespaces'teyseniz endişelenmenize gerek olmadığını unutmayın. Her şey düzgün bir şekilde kurulmuş olmalı. Ama ikinci olarak, Nextflow eğitim materyaline giderseniz, aşağı inerseniz, sürüm gereksinimleri hakkında konuşun, burada sürümleri keşfet etrafındaki yardım sayfasına götüren bir bağlantı var ve bu tüm bunları ayrıntılı olarak ele alıyor.

Bir dakikanız varsa bunu okumaya değer. çünkü Nextflow kullanmaya başladığınızda duyabileceğiniz bazı farklı terimlerin ne olduğunu netleştirmeye yardımcı olur. DSL1, DSL2, sözdizimi ayrıştırıcı bir, sözdizimi ayrıştırıcı iki ve benzeri gibi şeyler. Yani sadece bunun üzerinden bir kontrol yapmaya değer ve bu az önce söylediklerimin bir kısmını tekrarlıyor.

Ayrıca daha önce Nextflow kodu yazdıysanız ve bir tazeleme için geri geliyorsanız gerçekten yararlıdır. Size değişen bazı şeyleri söyler ve Nextflow kodunuzu nasıl güncelleyeceğinizi söyleyen Nextflow belgelerinin bölümlerine bağlantı verir.

## Kurs Dosyaları

Tamam. Kendimizi tanıştırmak için son şey, bu dizindeki dosyaları görmek. Kenar çubuğuna bakabilir veya genellikle eğitim materyalinde, tree komutunu kullanırız, -L, bakmak için seviye sayısıdır. İki diyeceğiz ve bunu tam ekran yaparsam, bunun temelde kenar çubuğunda gördüğümüz şeyi tam olarak yansıttığını göreceksiniz, ancak bir noktayla başlayan gizli dosyaları hariç tutar.

Yani \*.nf dosyaları, Nextflow anlamına gelir. Yani bunlar Nextflow betik dosyalarıdır ve burada eğitim materyalinin farklı bölümlerinin her biri için bir başlangıç dosyası var, açacağımız ve keşfedeceğimiz ve sonra düzenleyeceğimiz.

İlerledikçe bu dosyaları değiştireceğiz ve bu yüzden her bölümün sonunda, dosyalar bir sonraki bölümün başlangıcıyla hemen hemen aynı görünmelidir. Ama size bu farklı dosyaları veriyoruz, böylece her zaman taze başlayabilir ve sözdizimini karıştırmak konusunda fazla endişelenmezsiniz.

Kesinlikle çalışması gereken bir şeyle karşılaştırmanız gerekiyorsa. Çözümler klasörüne bakabilirsiniz ve bu, bölümlerin her biri için son durum gibidir, böylece yazdığınızı oradakilerle karşılaştırabilirsiniz.

Bir veri dizini var. Bu, kursun bir kısmında örnek girdi verisi olarak kullanacağımız sadece bir greetings.csv dosyasına sahip ve bir yapılandırma dosyası ve daha sonra kursta açıklayacağımız bazı parametreler gibi şeyler.

## Tamamlama

Tamam, yani şimdi umarım her şey çalışıyor. Ekranınız benimkiyle aynı görünüyor ve her şeye nasıl ulaşacağınızı ve tüm farklı dosyaların ne olduğunu anlıyorsunuz.

Başlarken sayfasının en altına kaydırırsanız, ne yaptığımı anladığımı söylemesi gereken küçük onay kutusu. Ortamım çalışıyor ve ayarlandınız, çalışma dizininiz "Hello Nextflow" klasörüne düzgün bir şekilde ayarlandı.

Tüm bunları işaretlediyseniz ve yeşil görünüyorlarsa. Bir sonraki videoya ve bir sonraki bölüme devam edebiliriz, bu birinci bölüm. Hello World. Birazdan görüşürüz.
