# Bölüm 2: Merhaba Kanallar - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!! note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım tam talimatlar için [kurs materyaline](../02_hello_channels.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca gösterge amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba ve Hello Nextflow'un 2. Bölümüne tekrar hoş geldiniz. Bu bölümün adı Merhaba Kanallar.

Kanallar, Nextflow iş akışınızdaki yapıştırıcı gibidirler. Nextflow'un tüm bilgileri aktarmak ve iş akışınızı yönetmek için kullandığı, farklı süreçleri bir arada tutan parçalardır.

Kanalların bir başka parçası da operatörlerdir. Bunlar temelde kanallarda içerikleri değiştirmek için kullanabileceğimiz fonksiyonlardır. Hadi VS Code'a dalıp nerede olduğumuza bakalım.

Bu VS Code'da çok yakınlaştım, bu yüzden işleri temiz ve düzenli tutmak için tüm _.nextflow\*_ dosyalarını ve _work/_ dizinini ve results/ dizinini ve Birinci Bölüm'deki her şeyi kaldırdım. Ve burada sıfırdan başlıyorum. Ama bunun için endişelenmeyin. İstemiyorsanız, o dosyaları orada bırakabilirsiniz. Herhangi bir soruna neden olmazlar.

Bu bölüm için _hello-channels.nf_ üzerinde çalışmaya başlayacağız ve bunu açarsam, daha önce üzerinde çalıştığımız dosyaya çok benzer görünmeli. Betiğin farklı kısımlarının farklı yerlerde olması mümkün, ama her şey temelde aynı olmalı.

Farklı olan bir şey, buradaki output bloğundaki yolun şimdi bu bölüm için _hello_channels_ olması, bu da sonuç dosyalarının eğer hâlâ oradaysa sonuçlarınızın farklı bir alt dizininde saklanacağı anlamına gelir. Bu yüzden çıktılar hakkında kafanız karışmadan başlamak için güzel ve temiz bir yer olmalı.

Tamam, bu iş akışını çalıştırdığımızda bu betiğin ne yaptığını hızlıca hatırlayalım. _"nextflow run hello-channels.nf"_ yapıyoruz. _"--input myinput"_ yapabiliriz ve bunu çalıştırdığımızda, yukarıdaki sayHello sürecinde greeting'e giden ve output.txt'ye kaydedilen bu params.input parametresini kullanacak. Ve bunu sonuç dosyasında görebiliriz. Harika.

## 1. Değişken girdileri açıkça bir kanal aracılığıyla sağlayın

Bu güzel. Ama oldukça basit. Bu parametrede bir değişkenimiz var, bu değişken bir kez çalışan bir sürece gidiyor ve gerçekten ölçeklenmiyor. Ve ona birçok farklı dosya veremiyoruz. Ona birçok farklı selamlama veremiyoruz. Sadece bir tanemiz var.

Gerçekte, Nextflow tamamen analizinizi ölçeklendirmekle ilgilidir. Yani muhtemelen birden fazla şey yapmasını istiyorsunuz. Ve bunu _kanallarla_ yapıyoruz.

Kanallar, Nextflow'u öğrenen birçok kişi için biraz benzersiz bir kavramdır. Fonksiyonel programlama kavramlarından gelir ve kafanıza yerleşmesi biraz zaman alabilir, ama bir kez tıkladığında, gerçekten Nextflow'un gücünü açarlar ve iş akışlarınızı nasıl yazacağınızın anahtarıdır.

## 1.1. Bir girdi kanalı oluşturun

Bu betiği alıp sadece bir _param_ yerine bir _kanal_ kullanmasını sağlayalım.

İş akışına gidiyoruz, bu bizim iş akışı mantığımızın işleri birbirine bağlama hakkında olduğu yer. Ve buraya gireceğim ve yeni bir kanal oluşturacağım.

Yeni bir kanal oluştur.

Ve ona "_greeting_ch"_ adını vereceğim. Bu değişkenin bir kanal olduğunu hatırlamanız için böyle "_\_ch"_ yapmak bir gelenektir. Ama istediğiniz adı verebilirsiniz.

Ve sonra eşittir diyeceğim ve _"channel.of"_ yapacağım.

Channel, kanallarla ilgili her şey için isim alanı gibidir. Küçük harf "c", eğer daha önce Nextflow kullanıyorsanız. Ve _".of"_, kanal oluşturmanın bir yolu olan Channel factory adı verilen bir şeydir.

Birçok farklı channel factory vardır. Burada sadece "." yaparsam, VS Code'un bir sürüsünü önerdiğini görebilirsiniz, ama _".of"_ en basit olanıdır ve sadece buraya bir girdi alır.

Yani bazı parantezler yapabilirim ve _"Hello Channels!"_ diyeceğim.

Harika. Bir kanalım var. Fantastik. Kaydet'e basabilirim, tekrar çalıştırabilirim, ama ilginç bir şey olmayacak. VS Code bana burada turuncu bir uyarı çizgisi verdi ve bunu kurduğunuzu söyledi: bunu oluşturdunuz, ama aslında hiçbir şey için kullanmadınız. Bu kanal tüketilmiyor.

Tamam, peki nasıl kullanırız? Çok basit. Bunu alacağım, kopyalayacağım ve _params.input_'u sileceğim ve bunun yerine _"greeting_ch"_ koyacağım. Yani bu kanalı sayHello'ya girdi olarak geçireceğiz.

Şimdilik bu dizgiyi sabit kodladığımı unutmayın. Bu, son bölümün sonunda kullandığımız güzel parametreden sonra biraz geri bir adım, ama mantığı görebilmeniz için işleri basit tutmak için.

Tamam, terminalime gireceğim ve iş akışını tekrar çalıştıracağım. Bu sefer herhangi bir _"--input"_ olmadan, ve çalışacak ve oluşturduğumuz kanalı kullanacak ve umarım _results/hello_channels/_ içinde bir dosyamız olmalı ve şimdi "Hello Channels!" diyor. Fantastik. Yani kanalımızdan beklediğimiz bu. Harika.

## 1.4. Kanal içeriğini incelemek için view() kullanın

Buraya eklenecek bir şey daha, kanallarda kullanabileceğimiz "_.view"_ adlı başka bir fonksiyona hızlı bir giriş.

Bu, Python veya alışkın olabileceğiniz diğer dillerdeki _print_ komutuna benzer ve çalıştırdığımızda bu kanalın içeriğini terminale döker.

Yani "_.view"_ yapın, ve sonra iş akışını tekrar çalıştırırsam, bu kanalın içeriğinin ne olduğunu, oluşturduğumuz zamanda, terminale yazdırmalı.

Nitekim, terminale buraya yazdırdığını görebilirsiniz. _"Hello Channels!"_.

İsterseniz bunları satırlar arasında bölebileceğinizi unutmayın ve aslında, Nextflow otomatik biçimlendiricisi bunu sizin için yapmaya çalışacaktır. Boşluk burada gerçekten önemli değil, bu yüzden bunları birbiri ardına zincirleyebilirsiniz.

## 2. İş akışını birden fazla girdi değeri üzerinde çalışacak şekilde değiştirin

Tamam, kanalımızda bir şey var ki bu güzel, ama temelde daha öncekiyle aynı. O halde biraz daha karmaşık hale getirelim. Kanalımıza birkaç şey daha ekleyelim.

"_.of()"_ channel factory birden fazla öğe alabilir, o halde birkaç tane daha yazalım. _Hello, Bonjour, Hej_ yapacağız. Ve sonra bu iş akışını tekrar çalıştırabiliriz ve ne olacağını göreceğiz.

Tekrar çalışmalı. Ve şimdi yazdırdık. _"Hello", "Bonjour"_ ve _"Hej"_ terminale view ifademizle. Fantastik.

## 2.1.2. Komutu çalıştırın ve log çıktısına bakın

Bu noktada işimizin bittiğini düşünebilirsiniz. Ama aslında burada bizi yakalayacak bir tuzak var. Çıktı dosyamıza buraya bakarsak. _"Hello"_ içinde olduğunu görebilirsiniz, ama diğer çıktılardan hiçbiri yok. Aslında, sadece bu.

Bu iş akışını birden çok kez çalıştırırsak, bazen _"Bonjour"_ olduğunu, bazen _"Hej"_ olduğunu bile görebiliriz. Biraz rastgele.

Terminale bakarsak, üç kez çalıştığını görebiliriz ve farklı view çıktılarını görebiliriz. Ama work dizinine gidersem, _"cat work"_ yapabilirim. Bu hash'i koy ve genişlet ve _output.txt_. Work dizinindeki bu dosyanın results dizininden farklı olduğunu görebilirsiniz ve bu _"Hej"_. Yani burada tam olarak çalışmayan bir şey var.

Ve anahtar şu ki, üç görev çalıştı. Nextflow çıktısı, işlem devam ederken bunu özetlemeye çalışır, böylece tüm terminalinizi tamamen ele geçirmez ve bu ANSI Günlüğü, temelde diğer görevlerin üzerine yazan ANSI kaçış kodlarını kullanır. Yani size sadece güncellenen son olanı gösterir.

## 2.1.3. Komutu -ansi-log false seçeneği ile tekrar çalıştırın

Bunu gerçekten daha iyi anlamak için yapabileceğimiz birkaç şey var. Work dizininin içine bakabiliriz ve oradaki tüm farklı work dizinlerini görebilirsiniz, ama bu biraz kafa karıştırıcı çünkü farklı Nextflow yürütme çalıştırmalarıyla karışacak.

Ya da Nextflow'a ANSI kaçış kodlarını kullanmamasını söyleyebiliriz.

Yani komutu tekrar çalıştırırsam, ama bu sefer _"-ansi-log false"_ diyerek kapatıyorum, ayrıca _$NO_COLOR_ veya _"$NXF_ANSI_LOG=false"_ ortam değişkenlerini de kullanabilirim. O zaman bu kaçış kodlarından herhangi biri olmadan Nextflow günlüğünün daha eski tarzını kullanır. Akıllı güncellemeler olmadan doğrudan terminale yazdırır.

Ve şimdi çalışan bu üç sürecin tümünü görebiliriz. Ve her birinin kendi görev hash'i var. Ve bu work dizinlerine girersek, belirttiğimiz üç farklı selamlamayı göreceğiz.

Yani bu artık biraz daha mantıklı. Umarım Nextflow'un bunu yaptığını anlarsınız, sadece bu work dizinleriyle size terminalde gösterdiği şeyde biraz akıllıydı.

Ancak, bu work dizinleriyle ilgili bir sorunu düzeltti, ancak çıktı dosyasıyla ilgili bir sorunu düzeltmedi. Hâlâ sadece _"Hello"_ yazan bir çıktı dosyamız var.

## 2.2. Çıktı dosyası adlarının benzersiz olacağından emin olun

Şimdi bunu anlamak için iş akışı betiğimize geri dönmemiz gerekiyor. Burada kanalımızı oluşturuyoruz, onu sürecimize geçiriyoruz ve sürece bakarsak, selamlamayı _"output.txt"_ adlı bir dosyaya yazıyoruz ve o çıktı dosyasını aşağıdaki output bloğuna geri geçirerek yayınlıyoruz.

Ancak, bu süreç her üç kez çalıştığında bu üç farklı görev. Hepsi _"output.txt"_ adlı bir dosya oluşturur, tüm bu çıktı dosyaları results dizinine yayınlanır ve hepsi birbirinin üzerine yazar. Yani orada aldığınız sonuç dosyası her ne ise, sadece oluşturulan son olandır, ancak diğerlerinin hepsini sildi. Bu gerçekten istediğimiz şey değil.

## 2.2.1. Dinamik bir çıktı dosyası adı oluşturun

Bunu ele almanın farklı yolları var, ama şimdilik en basiti sadece farklı benzersiz dosya adları oluşturmak. Böylece görev her farklı selamlama ile çalıştığında, farklı bir çıktı dosyası oluşturacak, bu da yayınlandığında artık çakışmayacak. Ve sonra üç benzersiz çıktı dosyası alacağız.

Bunu tamamen aynı şekilde yapıyoruz. Bu değişkeni script bloğu içinde herhangi bir yerde kullanabiliriz ve birden çok kez kullanabiliriz.

Yani bunu buraya yapıştırabilirim, _"$\{greeting\}\_output.txt"_, ve sonra bunu buraya da yapıştırmam gerekiyor çünkü artık _output.txt_ adlı bir dosya oluşturmuyoruz. Yani bunu güncellemezsen, Nextflow hiç oluşturulmamış bir dosya beklediğini söyleyen bir hatayla çökecektir.

Yani aynısını orada yapmam gerekiyor ve bu değişkenin anlaşılması için tek tırnak değil çift tırnak kullanmam gerekiyor.

Tamam, hadi deneyelim ve çalışıp çalışmadığına bakalım. İş akışını tekrar çalıştıracağız. Umarım bize üç farklı work dizinindeki üç farklı görevi gösterecektir. Ve nitekim, sol tarafta burada results klasörünün içinde artık üç farklı dosyayla üç farklı dosya adı ve her biri beklediğimiz farklı içeriklerle var. Yani dosyalar artık birbirinin üzerine yazmıyor ve her şey beklediğimiz gibi orada. 

Burada yaşadığımız bu biraz önemsiz bir kurulum, ama dosya yayınlamanın nasıl çalıştığı hakkında anlamanız gereken bazı temel kavramları vurguluyor ve tuzak olarak düşebileceğiniz bazı şeyler. Umarım kendi iş akışlarınızda bundan kaçınabilirsiniz.

Ayrıca burada yaptığımız şeyin gerçek hayat durumlarında biraz pratik olmadığını belirtmekte fayda var. Bazı girdi verilerini aldık ve bu verileri kullanıyoruz, ama aynı zamanda dosyayı bu verilerden sonra adlandırıyoruz, ki bunu genellikle yapamazsınız.

Yani gerçek daha olgun Nextflow iş hatlarında, genellikle belirli bir örnekle ilişkili tüm meta verileri içeren bir meta nesnesini etrafta taşırsınız. Daha sonra buna dayalı olarak dinamik dosya adları oluşturabilirsiniz, ki bu çok daha pratiktir.

Bunu en iyi uygulamalarla nasıl yapacağınızla ilgileniyorsanız, _training.nextflow.io_'da özellikle meta veriler ve meta haritalar hakkında bir yan görev var, daha fazla detay için oraya girebilirsiniz.

## 3. Bir dizi aracılığıyla birden fazla girdi sağlayın

Tamam. Şimdi kanalların nasıl yapılandırıldığı ve kodlama dilindeki diğer veri yapıları türlerinden nasıl farklı olduğu hakkında biraz keşfedeceğiz. Ve diğer dillerden geldiyseniz tanıdık bir kavram olabilecek bir diziyi nasıl potansiyel olarak kullanabileceğimi düşüneceğim.

Bir kanalda bir dizi kullanabilir miyim? Hadi deneyelim. Bir dizi oluşturacağım ve bunu dokümantasyondan kopyaladım, _"greetings_array"_ ve _"Hello", "Bonjour"_ ve _"Holà"_. Ve sonra bunu sabit kodlanmış dizgilerim yerine buraya koyacağım. Yani "channel.of" _"greetings_array"_ diyeceğim, bu diziyi bir kanala geçirerek. Hadi deneyelim.

Terminali açın ve iş hattını çalıştırın.

Tamam. View ifadesinin beklendiği gibi dizimizi yazdırdığını görebilirsiniz, ama sonra tüm bu kırmızı metin, veya hâlâ _"-ansi-log"_'u kapatmışsanız kırmızı olmayacak, ama tüm bu kırmızı metin bize bir şeylerin ters gittiğini söylüyor.

Artık burada güzel bir yeşil onay işaretimiz yok. Kırmızı bir çarpımız var ve bunu biraz daha genişletirsem okunması daha kolay olsun, Nextflow bize neyin yanlış gittiğini söylüyor.

Hadi bunu bölüm bölüm inceleyelim. Hatanın nedeni diyor ve sonra hatanın nedeni, eksik çıktı dosyaları. Yani temelde o output bloğu bu dosyanın oluşturulması gerektiğini söyledi ama oluşturulmadı. Sonra bu yürütülen komut diyor. Yani bu temelde o _.command.sh_ dosyasının içeriği. Tüm bu değişkenler konulduktan sonra böyle görünüyordu.

Ve burada echo komutumuzu görebilirsiniz aslında sadece bir kez çalıştırıldı ve tüm diziyi kullandı, ama bir dizi temsili olarak, ki bu gerçekten istediğimiz şey değildi.

Ve sonra komut öyle çıktı ve bu bizim gidip dosyaları görebileceğimiz, biraz daha anlayabileceğimiz work diziniydi.

Tamam. O zaman olan şey. Nextflow sadece bu tüm diziyi sürece tek bir kanal elemanı olarak geçirdi, bu da sürecin sadece bir kez çalıştığı anlamına geliyordu. Bir görevi vardı ve verileri beklediğimiz yapıda kullanmadı.

## 3.2. Kanal içeriğini dönüştürmek için bir operatör kullanın

Yani kullanılabilmeden önce bu kanala bir şey yapmamız gerekiyor. Ve bu, kanal içeriklerini manipüle etmek için kanallarda kullanabileceğimiz özel fonksiyonlar olan operatörleri kullanmanın aşamasını hazırlıyor.

Bu durumda, _flatten_ adı verilen bir şey kullanacağız. Bunu burada kanalın sonuna geçiriyoruz. Yani kanalı oluşturuyoruz ve sonra _flatten_ çalıştırıyoruz. Ve yine üzerine gelirseniz, bu komutun dokümantasyonunu VS Code'da hemen gösterir, bu çok yararlı. Tüm bu dokümantasyonları Nextflow web sitesinde, dokümantasyonda da bulabilirsiniz.

Şimdi bu kodu çalıştırabilirdim ve çalışıp çalışmadığını görebilirdim, ama bu aynı zamanda operatörler içinde ve Nextflow kodu içinde closure adı verilen dinamik kod yapmayı tanıtmak için güzel bir fırsat.

Yani _flatten_ çalıştırmadan önce burada bir view komutu ekleyeceğim. Ve burada dinamik closure olan bu kıvrımlı parantezler var. Ve burada sadece view operatörü bağlamında yürütülecek bazı keyfi kodlar var.

Burada, bu greeting al diyor, ki bu view operatörünün girdileri, ve bu burada. Bunu istediğim gibi çağırabilirim, buna _"foo"_ diyebilirim ve sadece daha sonra _"foo"_ olarak ona atıfta bulunmam gerekir. Ve sonra bununla, bunu döndür diyorum.

Ve sonra bir değişken için flatten'den önce diyen bir dizi döndürüyor. çok basit.

Şimdi bunlardan bir tane daha ekleyeceğim tamamen aynı şekilde, ama _flatten_'den sonra diyeceğim.

Yani bunun yaptığı, bu sırayla çalıştığı için, _flatten_ çalıştırmadan önce kanalın nasıl göründüğünü göreceksiniz, ve sonra _flatten_ çalıştırdıktan sonra tekrar.

Ve sonra bu greeting kanalı hâlâ oluşturuldu, yani hâlâ sürece geçirilecek. Ve umarım artık iş akışı çalışacaktır. Hadi deneyelim.

Harika. Öncelikle, iş hattı bu sefer çökmedi. Düzgün çalışan üç sürecimiz vardı ve küçük bir onay işaretimiz var. Ve sonra view ifadelerimizin çalıştığını görebiliriz.

_flatten_'den önce var, ki bu daha önce hatadan gördüğümüz dizi, ve sonra _flatten_ çağrıldıktan sonra üç kez var, _"Hello", "Bonjour"_ ve dizideki diğer üç ayrı eleman var, ki bunlar artık umduğumuz gibi kanaldaki üç ayrı eleman.

Ve _view_ operatörünün üç kez çalıştırıldığını görebilirsiniz. Ve bu, _flatten_'den sonra bu kanalın artık üç elemanı olduğu için. Yani operatör üç kez çağrılır.

Çok hızlı, channel factory'leri oluşturmadan önce, _"."_ yaptım dediğimi belirtmek isterim, ve kanal oluşturmanın birçok farklı yolunun olduğunu gördük ve bunlardan biri "_fromList"_ adında. Ve bu aslında özellikle bu aynı operasyonu yapmak için tasarlandı. Yani sadece from list greetings away yapabilirdik, ve bu çalışacaktır. Biraz daha temiz ve güzel bir sözdizimi. Ama bu gösterimin amaçları için, kanalın nasıl manipüle edildiğini ve farklı operatörlerin bir kanalın içeriğinde içeriği nasıl değiştirebileceğini görebilmeniz için adım adım yapmak istedik.

## 4. Girdi değerlerini bir CSV dosyasından okuyun

Tamam, bunu biraz daha gerçekçi nasıl yapabiliriz? Muhtemelen Nextflow iş hattınızda sabit kodlanmış dizilerle çok fazla kod oluşturmak istemeyeceksiniz. Muhtemelen başlattığınızda verileri dışarıdan almak isteyeceksiniz ve bu veriler neredeyse kesinlikle dosyalarda olacak.

Yani yapacağımız bir sonraki şey, bunu tekrarlayacağız, ama verileri tek bir CLI parametresinden veya sabit kodlanmış bir dizgiden veya diziden almak yerine, bir dosyadan alacağız.

Yani greetings away'den kurtulalım. Ve şimdi bu channel factory'yi tekrar değiştireceğiz. Az önce aralarından seçim yapabileceğiniz bir grup olduğunu söyledim ve _".fromPath"_ adlı bir tane var. Ve ona, bu durumda, daha önce kullandığımız girdimiz olan _params.input_ al diyeceğim.

Şimdi o parametre henüz kullanılmaya hazır değil. Hâlâ bunun bir dizi olduğunu söylüyoruz ve varsayılan olarak burada sabit kodlanmış, ama o dizgiyi üzerine yazabiliriz. Şimdi bunun bir dosya olmasını istiyoruz. Yani tür farklı. Artık bir _String_ değil. Bir _Path_.

Ve sonra isterseniz varsayılanı yine bir Path'e ayarlayabiliriz. Ve soldaki keşfet'e bakarsam, bu depoda, bu çalışma dizininde, data adında bir dizinin olduğunu görebilirsiniz. Orada _"greetings.csv"_ adlı bir dosyam var.

Yani varsayılanı burada _"data/greetings.csv"_'ye ayarlayabilirim. Şimdi, bu iş hattını herhangi bir komut satırı seçeneği olmadan tekrar çalıştırdığımda, bu varsayılan değeri kullanacak. Bunun bir path olduğunu biliyor, bu yüzden bunu bir dizgi değil bir path olarak işlemesi gerektiğini biliyor.

Ve sonra bunu bu _params.input_'tan bir channel factory'ye geçirecek ve kanalımızı oluşturacak, bu da sonra _sayHello_ adlı bu süreçte kullanılacak. Hadi deneyelim.

Tamam. Başarısız. Endişelenmeyin. Bu bekleniyordu. Ve eğitim materyalini takip ediyorsanız, orada da beklendiğini göreceksiniz. Burada ne oluyor görelim.

İş hattını çalıştırmaya çalıştı. Süreci yürütmeye çalıştı ve daha önce gördüğümüze oldukça benzer bir hata aldı.

Burada şöyle diyor: _echo_ çalıştırmaya çalıştık, ama bu CSV dosyasının içeriğini yankılamak yerine, sadece yolu yankıladı. Ve buradaki bu CSV dosyasına giden tam mutlak yolu görebilirsiniz.

Ve sonra nitekim, bunu bu gerçekten karmaşık yola yazmaya çalıştığı için, ne yapacağını gerçekten bilmiyordu. Ve süreç work dizininin kapsamı dışındaydı.

Başlangıçta Nextflow'un her yürütülen görevi özel bir work dizini içinde kapsüllediğini söyledim. Ve o work dizininin dışında olan verileri yazmaya çalışırsanız, Nextflow güvenlik önlemi olarak sizi durduracak. Ve burada olan bu. Mutlak bir yola yazmaya çalıştık ve Nextflow başarısız oldu ve bizi engelledi.

## 4.2. Dosyayı ayrıştırmak için splitCsv() operatörünü kullanın

Tamam, bu kanala bir bakalım ve neye benzediğini görelim. _".view"_ yapabiliriz, ve bunu web sitesinden kopyaladım. Yani _.view_, ve burada dinamik bir closure var ve girdi olarak "_csv"_ değişken adını söylüyoruz. Yani bu kanal içeriği, ve splitCsv'den önce diyoruz, ve böyle görünüyor.

Tekrar çalıştırırsam, yine başarısız olacak, ama bu kanalın içinde ne olduğunu bize gösterecek. Özellikle heyecan verici değil. Bu _path_ değişkeni. Yani bir terminale yazdırıldığı için burada sadece bir dizi olduğunu görebilirsiniz, ama bu dosya hakkındaki bilgileri ve meta verileri içeren bir _path_ nesnesi.

Dosyanın meta verilerini girdiye geçirmek istemiyoruz. O dosyanın içeriğini geçirmek istiyoruz. _greetings.csv_ dosyasına bakarsak, burada bu farklı değişkenlerin olduğunu görebilirsiniz. _Hello, Bonjour, Holà_ tekrar. Ve bunlar gerçekten sürece geçirmek istediğimiz şeyler, sadece dosyanın kendisi tek bir nesne olarak değil.

Yani bu CSV dosyasını ayrıştırmamız gerekiyor. Onu açmamız, CSV dosyasının içeriğine ulaşmamız ve sonra içerikleri kanal içinde sürece geçirmemiz gerekiyor.

Günlük mesajından muhtemelen anlayabileceğiniz gibi, başka bir operatör, başka bir kanal operatörü olan _splitCsv_'yi kullanmak istiyoruz. Yani "_dot" "s"_ yaparsam, ve sonra otomatik önerildiğini görebilirsiniz. Hata, _splitCsv_ ve bazı parantezler.

Ve sonra _splitCsv_'den sonra, nasıl göründüğünü görebilmemiz için başka bir _view_ ifadesi koyacağım. Hadi iş hattını çalıştıralım ve ne aldığımızı görelim.

Tamam. Hâlâ başarısız oldu, ama yeni ve heyecan verici bir şekilde, bu ilerleme.

Bu sefer yine betiğimizle ilgili bazı sorunlarımız var, ki bu render edildi. Şimdi. Artık son yolu almadık, ama bir değişken dizisi aldık, ki bu daha önce bir diziyi sabit girdi olarak geçirdiğimizde sahip olduğumuz hataya çok benziyor.

View operatöründen gelen günlüğümüzle, _splitCsv_'den öncenin path olduğunu görebiliriz. Ve nitekim, _splitCsv_'den sonra, üç farklı çıktımız var ve bu çıktıların her biri _greetings.csv_ dosyasındaki satırların her birine çok benziyor, ki bu mantıklı.

Yani burada olan, Nextflow'un bu CSV dosyasını ayrıştırdığı, bize CSV dosyasının her satırı için bir dizi olmak üzere üç nesne verdi. Yani sonra üç kez tek bir dizi değeri yerine bir değişken dizisini kanala geçirdik.

Tamam, yani geçen sefer bu sorunu yaşadığımızda, _flatten_ kullandık. Hadi çok hızlıca. Flatten'ı deneyelim ve ne olacağını görelim.

Bu değişkenlere ne dersem diyebilirim. Yani ona _myarray_ diyeceğim çünkü artık gerçekten bir CSV değil. Hadi tekrar çalıştırmayı deneyelim ve _flatten_ ile ne olacağını görelim.

Yani bu sefer çalışacağız, CSV'yi üç dizi nesnesine ayrıştırdık, ve sonra düzleştirdik. Ve bu sefer geçti. Ve Nextflow iş hattı çalıştı. Ancak _flatten_'in gerçekten işi abartıp her şeyi düzleştirdiğini görebilirsiniz. Ve her satır için üç bağımsız dizi girdisi alırız. Ve bu yüzden süreci CSV'nin her satırı için üç kez çalıştırdı. Ve şimdi bir sürü sonuç dosyamız var, ve 123, 456, ve her türlü şeyler, sadece CSV'nin ilk sütunu değil, ki bu gerçekten istediğimiz şeydi.

## 4.3. Selamlamaları çıkarmak için map() operatörünü kullanın

Peki sadece ilk sütuna nasıl ulaşırız? Eğer flatten burada çok basitse, aslında özelleştirebileceğimiz ve CSV'den ne istediğimizi söyleyebileceğimiz daha karmaşık bir operatöre ihtiyacımız var.

Bunu yapmak için _map_ kullanacağız. Temelde _map_, bana verilen her eleman üzerinde bazı kodları, bazı fonksiyonları çalıştır ve üzerinde bir tür dönüşüm yap der. Ve çok esnek olduğu için, Nextflow kodunda her zaman ortaya çıktığını göreceksiniz.

Tek başına hiçbir şey yapmaz. Yani normal parantez istemiyoruz, burada bir closure istiyoruz ve ne yapacağımızı söylememiz gerekiyor. Yani _"row"_ diyeceğim, çünkü CSV'den satırlar veriliyor, yani mantıklı bir değişken adı. Girdi. Ve sadece o dizinin ilk elemanını döndürmek istiyorum.

Nextflow'daki diziler sıfır tabanlıdır, yani sadece ilk eleman diyeceğiz, ki bu sıfırıncı satır. İkinci sütunu isteseydik, bir veya üçüncü sütun iki olabilirdi, ve böyle devam eder. Burada istediğimizi döndürebiliriz, ama sadece ilk değeri döndüreceğim.

Ve şimdi, iş hattını tekrar çalıştırabiliriz ve beklediğimizi yapıp yapmadığını görebiliriz.

Nitekim, _splitCsv_'den sonra dizilerimiz var, ve sonra _map_'ten sonra sadece _"Hello", "Bonjour"_ ve _"Holà"_ olmak üzere güzel temiz dizgilerimiz var. Ve iş hattı şimdi istediğimizi yapıyor. Fantastik.

Yani tüm bu view komutlarından artık kurtulabiliriz. Artık onlara ihtiyacımız yok.

## Özet

Hata ayıklamamızı bitirdik ve sonuçta bu kodla karşılaştık. CLI parametremiz olan _input_'u alarak, ki bu bir _Path_ olarak sınıflandırıldı. Nextflow yolu bulur, yükler ve CSV dosyasını anlar. Tüm farklı satırları döndürür. Ve sonra o satırın sadece ilk elemanını, kanal içeriklerini veren, sürece geçirilen kanala haritalıyoruz.

Ve süreç kanaldaki her eleman üzerinde çalışır, ki bu üç. Ve süreci üç kez çalıştırır, ona üç görev verir. Ve o sonuçlar daha sonra iş akışından yayınlanır, süreç çıktısı tarafından alınır. Bir iş akışından yayınlanır ve output bloğunda _"hello_channels"_ adlı bir alt dizine kaydedilir.

Oldukça havalı. Artık gerçek bir analiz için çalıştırabileceğiniz gerçek bir Nextflow iş hattına daha yakın bir şeye geliyoruz.

## Özet

Tamam. Umarım artık Nextflow kanallarının ve operatörlerin ne olduğunu ve operatörlerin kanallarda nasıl çalıştığını ve bunları nasıl oluşturabileceğinizi hissediyorsunuzdur.

Bu videonun başında söylediğim gibi, kanallar Nextflow'un yapıştırıcısıdır. Ve burada farklı girdileri alabileceğimizi, onları manipüle edebileceğimizi ve bu verileri alabileceğimizi ve sonra onları aşağı akış iş akışı mantığına geçirebileceğimizi görebilirsiniz.

Ve buradaki bu workflow bloğu gerçekten tüm o paralelleştirmeyi ve tüm akıllı mantığı oluşturduğunuz yerdir ve Nextflow'a iş akışı DAG'inizi nasıl oluşturacağını ve iş hattınızı nasıl yöneteceğini açıklarsınız.

Kanallar kafanıza yerleştirilmesi en kolay kavram değildir. Bu yüzden bir mola verin, bunun hakkında biraz düşünün, belki materyali tekrar okuyun ve bu kavramları gerçekten anladığınızdan emin olun çünkü bu Nextflow'u anlamanız için anahtardır ve kanalları, farklı kanal operatörlerini ve farklı channel factory'leri ne kadar iyi anlarsanız. Nextflow yazarken o kadar eğleneceksiniz ve iş hatlarınız o kadar güçlü olacak.

Bu Python veya diğer dillerdeki normal programlama ile aynı değil. Burada _if_ ifadeleri kullanmıyoruz, bu kanallar ve operatörler kullanan fonksiyonel akış programlamadır. Yani biraz farklı, ama aynı zamanda süper güçlü.

Bu bölümün sonu. Gidin hızlı bir mola verin ve bir sonraki videoda sizi Part three'de göreceğim, Merhaba İş Akışı'ndan geçeceğiz ve iş akışları hakkında biraz daha konuşacağız.

Önceki bölüm gibi, burada web sayfasının altında birkaç test sorusu var, bu yüzden bunlara hızlıca göz atabilirsiniz ve az önce yaptığımız materyalin tüm farklı bölümlerini anladığınızdan emin olabilirsiniz. Ve bunun dışında, bir sonraki videoda sizinle görüşeceğim. Çok teşekkürler.

Tamam.
