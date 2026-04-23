# Bölüm 2: Hello Channels - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [eğitim materyaline](../02_hello_channels.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba ve Hello Nextflow'un 2. Bölümüne tekrar hoş geldiniz. Bu bölümün adı Hello Channels.

Kanallar, Nextflow pipeline'ınızdaki yapıştırıcı gibidir. Nextflow'un tüm bilgileri aktarmak ve iş akışınızı düzenlemek için kullandığı, tüm farklı süreçleri bir arada tutan parçalardır.

Kanalların bir başka bileşeni de operatörlerdir. Bunlar temelde kanallardaki içerikleri değiştirmek için kullanabileceğimiz fonksiyonlardır. Hadi VS Code'a geçip nerede olduğumuza bakalım.

VS Code'da çok yakınlaştırdım; bu yüzden işleri temiz ve düzenli tutmak için tüm _.nextflow\*_ dosyalarını, _work/_ dizinini, results/ dizinini ve Birinci Bölüm'den kalan her şeyi kaldırdım. Burada sıfırdan başlıyorum. Ama buna çok takılmayın. İstemiyorsanız bu dosyaları bırakabilirsiniz; herhangi bir soruna yol açmazlar.

Bu bölüm için _hello-channels.nf_ üzerinde çalışmaya başlayacağız. Bu dosyayı açarsanız, daha önce üzerinde çalıştığımız dosyaya çok benzer görünmeli. Farklı bölümlerin betiğin farklı yerlerinde olması mümkün; ancak her şey temelde aynı olmalı.

Farklı olan bir şey, output bloğundaki yolun artık bu bölüm için _hello_channels_ olması. Bu, hâlâ orada varsa sonuç dosyalarının results dizininde farklı bir alt dizinde saklanacağı anlamına geliyor. Çıktılar konusunda kafanız karışmadan başlamak için güzel ve temiz bir yer olmalı.

Tamam, bu iş akışını çalıştırdığımızda bu betiğin ne yaptığını hızlıca hatırlayalım. _"nextflow run hello-channels.nf"_ komutunu çalıştırıyoruz. _"--input myinput"_ ekleyebiliriz; çalıştırdığımızda, yukarıdaki sayHello sürecine girdi olarak geçirilen params.input parametresini kullanacak. Bu parametre greeting'e gidiyor ve output.txt'ye kaydediliyor. Bunu results dosyasında görebiliriz. Harika.

## 1. Değişken girdileri açıkça bir kanal aracılığıyla sağlayın

Bu güzel. Ama oldukça basit. Bu parametrede bir kez çalışan bir sürece giden tek bir değişkenimiz var ve bu gerçekten ölçeklenmiyor. Ona birçok farklı dosya veremiyoruz, birçok farklı selamlama veremiyoruz. Sadece bir tane var.

Gerçekte, Nextflow tamamen analizinizi ölçeklendirmekle ilgilidir. Bu yüzden muhtemelen birden fazla şey yapmasını istersiniz. Bunu da _kanallarla_ yapıyoruz.

Kanallar, Nextflow öğrenen birçok kişi için biraz benzersiz bir kavramdır. Fonksiyonel programlama kavramlarından gelir ve kafanızda oturması biraz zaman alabilir; ancak bir kez anladığınızda, Nextflow'un gücünü gerçekten açığa çıkarırlar ve iş akışlarınızı nasıl yazacağınızın anahtarıdır.

## 1.1. Bir girdi kanalı oluşturun

Bu betiği alıp yalnızca bir _param_ yerine bir _kanal_ kullanacak şekilde değiştirerek başlayalım.

Tüm iş akışı mantığımızın bulunduğu workflow'a gidiyoruz. Buraya gireceğim ve yeni bir kanal oluşturacağım.

Yeni bir kanal oluşturun.

Buna "_greeting_ch"_ diyeceğim. "_\_ch"_ eklemek, bu değişkenin bir kanal olduğunu hatırlamanızı sağlayan bir gelenektir. Ama istediğiniz gibi adlandırabilirsiniz.

Sonra eşittir diyeceğim ve _"channel.of"_ yapacağım.

Channel, kanallarla ilgili her şey için ad alanı gibidir. Daha önce Nextflow kullandıysanız küçük "c" harfi. _".of"_ ise Channel factory denen bir şeydir; bu temelde kanal oluşturmanın bir yoludur.

Birçok farklı kanal fabrikası var. Burada sadece "." yaparsam, VS Code'un bir sürü öneri yaptığını görebilirsiniz; ancak _".of"_ en basit olanıdır ve burada sadece bir girdi alır.

Parantez açabilirim ve _"Hello Channels!"_ diyeceğim.

Harika. Bir kanalım var. Fantastik. Kaydet'e basabilirim, tekrar çalıştırabilirim; ama ilginç bir şey olmayacak. VS Code burada bana turuncu bir uyarı çizgisi verdi ve şunu söyledi: bunu oluşturdunuz, ama aslında hiçbir şey için kullanmadınız. Bu kanal tüketilmiyor.

Tamam, peki bunu nasıl kullanırız? Çok basit. Bunu alacağım, kopyalayacağım ve _params.input_'u sileceğim; yerine _"greeting_ch"_ koyacağım. Yani bu kanalı sayHello'ya girdi olarak geçireceğiz.

Şimdilik bu string'i sabit kodladığımı unutmayın. Bu, son bölümün sonunda kullandığımız güzel parametreden sonra biraz geri adım atmak gibi; ama mantığı görebilmeniz için işleri basit tutuyor.

Tamam, terminale gideceğim ve iş akışını tekrar çalıştıracağım. Bu sefer herhangi bir _"--input"_ olmadan çalışacak ve oluşturduğumuz bu kanalı kullanacak. Umarım _results/hello_channels/_ içinde bir dosyamız olmalı ve şimdi "Hello Channels!" diyor. Fantastik. Kanalımızdan beklediğimiz de buydu. Harika.

## 1.4. Kanal içeriğini incelemek için view() kullanın

Buraya eklenecek bir şey daha var: kanallarda kullanabileceğimiz "_.view"_ adlı başka bir fonksiyona hızlı bir giriş.

Bu, Python'daki _print_ komutuna veya alışık olabileceğiniz diğer dillerdeki benzer komutlara karşılık gelir ve çalıştırdığımızda bu kanalın içeriğini terminale döker.

"_.view"_ yapın; sonra iş akışını tekrar çalıştırırsam, oluşturduğumuz anda bu kanalın içeriğinin ne olduğunu terminale yazdırmalı.

Nitekim, terminale burada yazdırıldığını görebilirsiniz. _"Hello Channels!"_.

İsterseniz bunları satırlar arasında bölebileceğinizi unutmayın; aslında Nextflow otomatik biçimlendiricisi bunu sizin için yapmaya çalışacaktır. Boşluk burada gerçekten önemli değil, bu yüzden bunları birbiri ardına zincirleyebilirsiniz.

## 2. İş akışını birden fazla girdi değeri üzerinde çalışacak şekilde değiştirin

Tamam, kanalımızda bir şey var; bu güzel, ama temelde daha öncekiyle aynı. O halde biraz daha karmaşık hale getirelim. Kanalımıza birkaç şey daha ekleyelim.

"_.of()"_ kanal fabrikası birden fazla öğe alabilir; o halde birkaç tane daha yazalım. _Hello, Bonjour, Hej_ yapacağız. Sonra bu iş akışını tekrar çalıştırabiliriz ve ne olacağını göreceğiz.

Tekrar çalışmalı. View ifademizle terminale _"Hello"_, _"Bonjour"_ ve _"Hej"_ yazdırdık. Fantastik.

## 2.1.2. Komutu çalıştırın ve log çıktısına bakın

Bu noktada işimizin bittiğini düşünebilirsiniz. Ama aslında burada bizi tökezletecek küçük bir tuzak var. Çıktı dosyamıza bakarsak, içinde _"Hello"_ var; ama diğer çıktılardan hiçbiri yok. Aslında, sadece bu.

Bu iş akışını birden çok kez çalıştırırsak, bazen _"Bonjour"_ olduğunu, bazen _"Hej"_ olduğunu bile görebiliriz. Biraz rastgele.

Terminale bakarsak, üç kez çalıştığını ve farklı view çıktılarını görebiliriz. Ama work dizinine gidersem, _"cat work"_ yapabilirim. Bu hash'i koyun ve genişletin ve _output.txt_. Work dizinindeki bu dosyanın results dizininden farklı olduğunu görebilirsiniz; bu dosyada _"Hej"_ yazıyor. Yani burada tam olarak çalışmayan bir şey var.

Anahtar şu ki, çalışan üç görevimiz var. Nextflow çıktısı, tüm terminalinizi tamamen ele geçirmemesi için işlem devam ederken bunu özetlemeye çalışır; bu ANSI Logging, ANSI kaçış kodlarını kullanır ve temelde diğer görevlerin üzerine yazmıştır. Bu yüzden size güncellenen son olanı gösterir.

## 2.1.3. Komutu -ansi-log false seçeneğiyle tekrar çalıştırın

Bunu gerçekten daha iyi anlamak için yapabileceğimiz birkaç şey var. Work dizininin kendisine bakabiliriz ve oradaki tüm farklı work dizinlerini görebilirsiniz; ama bu biraz kafa karıştırıcı çünkü farklı Nextflow yürütme çalıştırmalarıyla karışacak.

Ya da Nextflow'a ANSI kaçış kodlarını kullanmamasını söyleyebiliriz.

Komutu tekrar çalıştırırsam, ama bu sefer kapatmak için _"-ansi-log false"_ dersem, _$NO_COLOR_ veya _"$NXF_ANSI_LOG=false"_ ortam değişkenlerini de kullanabilirim. O zaman bu kaçış kodları olmadan Nextflow'un daha eski tarz günlüğünü kullanır; hiçbir akıllı güncelleme olmadan doğrudan terminale yazdırır.

Ve şimdi çalışan bu üç sürecin tamamını görebiliriz. Her birinin kendi görev hash'i var. Bu work dizinlerine girersek, belirttiğimiz üç farklı selamlamayı göreceğiz.

Bu artık biraz daha mantıklı. Umarım Nextflow'un bunu yaptığını anlıyorsunuzdur; sadece bu work dizinleriyle terminalde size gösterdiği konuda biraz akıllıydı.

Ancak bu, work dizinleriyle ilgili bir sorunu çözdü; ama çıktı dosyasıyla ilgili sorunu çözmedi. Hâlâ sadece _"Hello"_ yazan tek bir çıktı dosyamız var.

## 2.2. Çıktı dosya adlarının benzersiz olacağından emin olun

Bunu anlamak için iş akışı betiğimize geri dönmemiz gerekiyor. Kanalımızı burada oluşturuyoruz, sürecimize geçiriyoruz ve sürece bakarsak, selamlamayı _"output.txt"_ adlı bir dosyaya yazıyoruz; bu çıktı dosyasını çıktı bloğuna geri geçirip yayımlıyoruz.

Ancak bu süreç, bu üç farklı görevi üç kez çalıştırır. Hepsi _"output.txt"_ adlı bir dosya oluşturur; tüm bu çıktı dosyaları results dizinine yayımlanır ve hepsi birbirinin üzerine yazar. Yani orada aldığınız sonuç dosyası, oluşturulan son dosyadır; diğerlerinin hepsini silmiştir. Bu gerçekten istediğimiz şey değil.

## 2.2.1. Dinamik bir çıktı dosya adı oluşturun

Bunu ele almanın farklı yolları var; ama şimdilik en basiti sadece farklı benzersiz dosya adları oluşturmak. Böylece görev her farklı selamlama ile çalıştığında farklı bir çıktı dosyası oluşturacak; bu da yayımlandığında artık çakışmayacak. Ve sonra üç benzersiz çıktı dosyası alacağız.

Bunu tamamen aynı şekilde yapıyoruz. Bu değişkeni script bloğu içinde herhangi bir yerde kullanabiliriz ve birden çok kez kullanabiliriz.

Buraya yapıştırabilirim, _"$\{greeting\}\_output.txt"_; sonra bunu buraya da yapıştırmam gerekiyor çünkü artık _output.txt_ adlı bir dosya oluşturmuyoruz. Bunu güncellemezseniz, Nextflow hiç oluşturulmamış bir dosya beklediğini söyleyen bir hatayla çökecek.

Yani orada da aynısını yapmam gerekiyor ve bu değişkenin anlaşılması için tek tırnak değil çift tırnak kullanmam gerekiyor.

Tamam, deneyelim ve işe yarayıp yaramadığına bakalım. İş akışını tekrar çalıştıracağız. Umarım bize üç farklı work dizinindeki üç farklı görevi gösterecek. Nitekim, solda results klasöründe görebilirsiniz. Şimdi üç farklı dosya adıyla üç farklı dosyamız var ve her biri beklediğimiz farklı içeriklere sahip. Yani dosyalar artık birbirinin üzerine yazmıyor ve her şey beklediğimiz gibi orada.

Burada geçtiğimiz bu biraz önemsiz bir kurulum; ama dosya yayımlamanın nasıl çalıştığı hakkında anlamanız gereken bazı temel kavramların altını çiziyor ve tuzak olarak düşebileceğiniz bazı şeyleri gösteriyor. Umarım bunu kendi iş akışlarınızda önleyebilirsiniz.

Burada yaptığımızın gerçek hayat durumlarında biraz pratik olmadığını da belirtmekte fayda var. Bazı girdi verilerini aldık ve bu veriyi kullanıyoruz; ama aynı zamanda dosyayı bu veriden sonra adlandırıyoruz ki bunu genellikle yapamazsınız.

Yani gerçek, daha olgun Nextflow pipeline'larında, genellikle belirli bir örnekle ilişkili tüm meta verileri içeren bir meta nesnesi geçirirsiniz. Daha sonra buna dayalı olarak dinamik dosya adları oluşturabilirsiniz; bu çok daha pratiktir.

Bunu en iyi uygulamalarla nasıl yapacağınızla ilgileniyorsanız, _training.nextflow.io_'da özellikle meta veriler ve meta haritalar hakkında bir yan görev var; daha fazla ayrıntı için oraya bakabilirsiniz.

## 3. Bir dizi aracılığıyla birden fazla girdi sağlayın

Tamam. Şimdi kanalların nasıl yapılandırıldığı ve kodlama dilindeki diğer veri yapılarından nasıl farklı olduğu hakkında biraz keşfedeceğiz. Diğer dillerden geldiyseniz tanıdık bir kavram olabilecek bir diziyi potansiyel olarak nasıl kullanabileceğimi düşüneceğim.

Bir kanalda dizi kullanabilir miyim? Deneyelim. Bir dizi oluşturacağım; bunu dokümanlardan kopyaladım: _"greetings_array"_ ve _"Hello"_, _"Bonjour"_ ve _"Holà"_. Sonra bunu sabit kodlanmış string'lerimin yerine buraya koyacağım. Yani "channel.of" _"greetings_array"_ diyeceğim; bu diziyi bir kanala geçireceğim. Deneyelim.

Terminali açın ve pipeline'ı çalıştırın.

Tamam. View ifadesinin burada beklediğimiz gibi dizimizi yazdırdığını görebilirsiniz; ama sonra tüm bu kırmızı metin — ya da hâlâ _"-ansi-log"_ kapalıysa kırmızı olmayacak — bize bir şeylerin yanlış gittiğini söylüyor.

Artık burada güzel bir yeşil onay işaretimiz yok. Kırmızı bir çarpımız var ve bunu biraz daha geniş yaparsam okumak daha kolay olur; Nextflow bize neyin yanlış gittiğini söylüyor.

Bunu bölüm bölüm inceleyelim. Hatanın nedeni olduğunu söylüyor ve sonra hatanın nedeni: eksik çıktı dosyaları. Yani temelde o çıktı bloğu bu dosyanın oluşturulması gerektiğini söyledi ve oluşturulmadı. Sonra bu yürütülen komut olduğunu söylüyor. Yani bu temelde o _.command.sh_ dosyasının içeriği; tüm bu değişkenler konulduktan sonra nasıl göründüğü bu.

Ve burada echo komutumuzun aslında sadece bir kez çalıştırıldığını ve tüm diziyi kullandığını görebilirsiniz; ama bir string gösteriminde, ki bu gerçekten istediğimiz şey değildi.

Sonra komut böyle çıktı ve bu, dosyaları görmek ve biraz daha anlamak için gidebileceğimiz work diziniydi.

Tamam. O zaman olan şey şuydu: Nextflow bu tüm diziyi sürece tek bir kanal öğesi olarak geçirdi; bu da sürecin sadece bir kez çalıştığı anlamına geliyordu. Bir görevi vardı ve veriyi beklediğimiz yapıda kullanmadı.

## 3.2. Kanal içeriğini dönüştürmek için bir operatör kullanın

Bu yüzden kullanılmadan önce bu kanala bir şey yapmamız gerekiyor. Ve bu, kanal içeriğini manipüle etmek için kanallarda kullanabileceğimiz özel fonksiyonlar olan operatörleri kullanmak için sahneyi hazırlıyor.

Bu durumda, _flatten_ adlı bir şey kullanacağız. Bunu burada kanalın sonuna geçiriyoruz. Yani kanalı oluşturuyoruz ve sonra _flatten_ çalıştırıyoruz. Yine, üzerine gelirseniz, bu komut için belgeleri doğrudan VS Code'da gösterir; bu çok yardımcı. Tüm bu belgeleri Nextflow web sitesinde, dokümantasyonda da bulabilirsiniz.

Bu kodu şimdi çalıştırabilir ve işe yarayıp yaramadığını görebilirim; ama aynı zamanda operatörler içinde ve Nextflow kodu içinde dinamik kod yapmanın nasıl yapılacağını tanıtmak için güzel bir fırsat; bunlara closure denir.

Bu yüzden _flatten_ çalıştırmadan önce buraya bir view komutu ekleyeceğim. Ve burada bu süslü parantezler var; bu dinamik closure. Ve burada sadece view operatörü bağlamında yürütülecek bazı rastgele kodlar var.

Burada, view operatörünün girdileri olan greeting'i al diyor ve bu burada. Buna istediğim gibi isim verebilirdim; buna _"foo"_ diyebilirdim ve sadece daha sonra _"foo"_ olarak ona atıfta bulunmam gerekir. Sonra bununla, bunu döndür diyorum.

Ve sonra bir değişken için flatten'dan önce yazan bir string döndürüyorum. Çok basit.

Şimdi bunun tamamen aynısını ekleyeceğim; ama _flatten_'dan sonra diyeceğim.

Yani bunun yaptığı şey şu: bu sırayla çalıştığı için, _flatten_ çalıştırmadan önce kanalın nasıl göründüğünü göreceksiniz ve sonra _flatten_ çalıştırdıktan sonra tekrar.

Ve sonra bu greeting kanalı hâlâ oluşturuldu, bu yüzden hâlâ sürece geçirilecek. Ve umarım şimdi iş akışı çalışacak. Deneyelim.

Harika. İlk olarak, pipeline bu sefer çökmedi. Düzgün çalışan üç sürecimiz vardı ve küçük bir onay işaretimiz var. Ve sonra view ifadelerimizin işe yaradığını görebiliriz.

_flatten_'dan önce var; bu daha önce başarısızlıktan gördüğümüz dizi. Sonra _flatten_ çağrıldıktan sonra üç kez var; _"Hello"_, _"Bonjour"_ ve dizideki diğer üç ayrı öğe var; bunlar şimdi umduğumuz gibi kanaldaki üç ayrı öğe.

Ve _view_ operatörünün üç kez çalıştırıldığını görebilirsiniz. Bunun nedeni, _flatten_'dan sonra bu kanalın artık üç öğesi olması. Bu yüzden operatör üç kez çağrılır.

Çok hızlı bir şekilde şunu belirtmek isterim: daha önce kanal fabrikaları oluştururken _"."_ yaptığımı söylemiştim ve kanallar oluşturmanın birçok farklı yolu olduğunu gördük; bunlardan biri "_fromList"_ olarak adlandırılır. Bu aslında özellikle bu aynı işlemi yapmak için tasarlanmıştır. Yani sadece from list greetings away yapabilirdik ve bu işe yarayacaktır. Biraz daha temiz ve güzel bir sözdizimi. Ama bu gösterimin amaçları için, kanalın nasıl manipüle edildiğini ve farklı operatörlerin bir kanalın içeriğini nasıl değiştirebileceğini görebilmeniz için adım adım yapmak istedik.

## 4. Girdi değerlerini bir CSV dosyasından okuyun

Tamam, bunu biraz daha gerçekçi nasıl yapabiliriz? Muhtemelen Nextflow pipeline'ınızda sabit kodlanmış dizilerle çok fazla kod oluşturmak istemeyeceksiniz. Muhtemelen başlattığınızda veriyi dışarıdan almak isteyeceksiniz ve bu veri neredeyse kesinlikle dosyalarda olacak.

Yani yapacağımız bir sonraki şey, bunu çoğaltacağız; ama veriyi tek bir CLI parametresinden veya sabit kodlanmış bir string ya da diziden almak yerine, bir dosyadan alacağız.

Greetings array'imizden kurtulalım. Ve şimdi bu kanal fabrikasını tekrar değiştireceğiz. Seçebileceğiniz bir sürü olduğunu söyledim ve _".fromPath"_ adlı bir tane var. Ona söyleyeceğim ki bu durumda, daha önce kullandığımız girdimiz olan _params.input_'u alsın.

Şimdi bu parametre henüz kullanılmaya hazır değil. Hâlâ bunun bir string olduğunu söylüyoruz ve burada bir varsayılanla sabit kodlanmış; ama bu string'i üzerine yazabiliriz. Şimdi bunun yerine bir dosya olmasını istiyoruz. Yani tür farklı. Artık bir _String_ değil. Bir _Path_.

Ve sonra istediğimiz takdirde varsayılanı yine bir Path'e ayarlayabiliriz. Soldaki explore'a bakarsam, bu depoda, bu çalışma dizininde, data adlı bir dizinimin olduğunu görebilirsiniz. Orada _"greetings.csv"_ adlı bir dosyam var.

Yani varsayılanı burada _"data/greetings.csv"_ olarak ayarlayabilirim. Şimdi, bu pipeline'ı herhangi bir komut satırı seçeneği olmadan tekrar çalıştırdığımda, bu varsayılan değeri kullanacak. Bunun bir path olduğunu biliyor; bu yüzden bunu bir string değil bir path olarak ele alması gerektiğini biliyor.

Ve sonra bunu bu _params.input_'tan bir kanal fabrikasına geçirecek ve kanalımızı oluşturacak; bu daha sonra _sayHello_ adlı bu süreçte kullanılacak. Deneyelim.

Tamam. Başarısız oldu. Endişelenmeyin. Bu bekleniyordu. Eğitim materyalini takip ediyorsanız, orada da beklendiğini göreceksiniz. Burada ne oluyor bakalım.

Pipeline'ı çalıştırmaya çalıştı. Süreci yürütmeye çalıştı ve daha önce gördüğümüze oldukça benzer bir hata aldı.

Burada diyor ki: _echo_ çalıştırmaya çalıştık; ama bu CSV dosyasının içeriğini yankılamak yerine, sadece yolu yankıladı. Ve görebilirsiniz ki bu CSV dosyasına tam mutlak yol burada.

Ve nitekim, bunu bu gerçekten karmaşık yola yazmaya çalıştığı için, gerçekten ne yapacağını bilmedi. Ve süreç work dizininin kapsamı dışındaydı.

Başlangıçta Nextflow'un her yürütülen görevi özel bir work dizini içinde kapsüllediğini söylemiştim. Ve bu work dizininin dışında olan verilere yazmaya çalışırsanız, Nextflow sizi bir güvenlik önlemi olarak durduracaktır. Ve burada olan da bu. Mutlak bir yola yazmaya çalıştık ve Nextflow başarısız oldu; bizi engelledi.

## 4.2. Dosyayı ayrıştırmak için splitCsv() operatörünü kullanın

Tamam, bu kanala bir bakalım ve nasıl göründüğünü görelim. _".view"_ yapabiliriz; bunu web sitesinden kopyaladım. Yani _.view_, ve burada dinamik bir closure'ımız var; girdi olarak "_csv"_ değişken adını söylüyoruz. Yani bu kanal içeriği ve splitCsv'den önce diyoruz; nasıl göründüğü bu.

Tekrar çalıştırırsam, yine başarısız olacak; ama bu kanalın içinde ne olduğunu bize gösterecek. Özellikle heyecan verici değil. O _path_ değişkeni. Yani görebilirsiniz ki terminale yazdırıldığı için burada sadece bir string; ama bu dosya hakkındaki bilgileri ve meta verileri içeren bir _path_ nesnesi.

Dosyanın meta verilerini girdiye geçirmek istemiyoruz. O dosyanın içeriğini geçirmek istiyoruz. _greetings.csv_ dosyasına bakarsak, burada bu farklı değişkenlerin olduğunu görebilirsiniz. Yine _Hello, Bonjour, Holà_. Ve bunlar gerçekten sürece geçirmek istediğimiz şeyler; sadece dosyanın kendisi tek bir nesne olarak değil.

Bu yüzden bu CSV dosyasını ayrıştırmamız gerekiyor. Onu açmamız, CSV dosyasının içeriğine ulaşmamız ve sonra içerikleri kanal içinde sürece geçirmemiz gerekiyor.

Log mesajından muhtemelen anlayabilirsiniz; başka bir kanal operatörü olan _splitCsv_'yi kullanmak istiyoruz. Yani "_dot" "s"_ yaparsam, otomatik önerildiğini görebilirsiniz. Oops, _splitCsv_ ve bazı parantezler.

Ve sonra _splitCsv_'den sonra, nasıl göründüğünü görebilmemiz için başka bir _view_ ifadesi koyacağım. Pipeline'ı çalıştıralım ve ne elde ettiğimize bakalım.

Tamam. Hâlâ başarısız oldu; ama yeni ve heyecan verici bir şekilde, bu ilerleme.

Bu sefer yine betiğimizle ilgili bir sorunumuz var; bu render edildi. Şimdi. Artık son yolu almadık; ama bir değişkenler dizisi aldık; bu da daha önce bir diziyi sabit girdi olarak geçirdiğimizde aldığımız hataya çok benziyor.

View operatöründen gelen günlüğümüzle, _splitCsv_'den önce yolun olduğunu görebiliriz. Ve nitekim, _splitCsv_'den sonra, üç farklı çıktımız var ve bu çıktıların her biri _greetings.csv_ dosyasındaki satırların her birine çok benziyor; bu mantıklı.

Yani burada olan şey, Nextflow'un bu CSV dosyasını ayrıştırması ve bize CSV dosyasının her satırı için bir dizi olmak üzere üç nesne vermesi. Yani sonra üç kez tek bir string değeri yerine kanala bir değişkenler dizisi geçirdik.

Tamam, geçen sefer bu sorunumuz olduğunda, _flatten_ kullandık. Hadi çok hızlı bir şekilde. Flatten'ı deneyelim ve ne olacağını görelim.

Bu değişkenleri ne olursa olsun çağırabilirim. Yani buna _myarray_ diyeceğim çünkü artık gerçekten bir CSV değil. Tekrar çalıştırmayı deneyelim ve _flatten_ ile ne olacağını görelim.

Yani bu sefer çalışacağız; CSV'yi üç dizi nesnesine ayrıştırdık ve sonra düzleştirdik. Ve bu sefer geçti. Nextflow pipeline'ı çalıştı. Ancak _flatten_'ın gerçekten işe koyulduğunu ve her şeyi düzleştirdiğini görebilirsiniz. Ve böylece her satır için üç bağımsız dizi girişi alıyoruz. Bu yüzden süreci her CSV satırı için üç kez çalıştırdı. Ve şimdi bir sürü sonuç dosyamız var; 123, 456 ve her türlü şey; sadece gerçekten istediğimiz CSV'nin ilk sütunu değil.

## 4.3. Selamlamaları çıkarmak için map() operatörünü kullanın

Peki sadece ilk sütuna nasıl ulaşırız? Flatten burada çok basitse, aslında özelleştirebileceğimiz ve CSV'den ne istediğimizi söyleyebileceğimiz daha karmaşık bir operatöre ihtiyacımız var.

Bunu yapmak için _map_ kullanacağız. Temelde _map_ sadece diyor ki: bana verilen her öğe üzerinde bir kod, bir fonksiyon çalıştır ve üzerinde bir tür dönüşüm yap. Ve çok esnek olduğu için, Nextflow kodunda her zaman karşınıza çıktığını göreceksiniz.

Kendi başına hiçbir şey yapmaz. Yani normal parantez istemiyoruz; burada bir closure istiyoruz ve ne yapacağını söylememiz gerekiyor. Yani _"row"_ diyeceğim; çünkü CSV'den satırlar veriliyor, bu yüzden mantıklı bir değişken adı. Girdi. Ve sadece o dizinin ilk öğesini döndürmek istiyorum.

Nextflow'daki diziler sıfır tabanlıdır; bu yüzden sadece ilk öğe olduğunu söyleyeceğiz; bu sıfır satırı. İkinci sütunu isteseydik, bir; üçüncü sütun iki olabilir ve böyle devam eder. Burada istediğimizi döndürebiliriz; ama sadece ilk değeri döndüreceğim.

Ve şimdi, pipeline'ı tekrar çalıştırabiliriz ve beklediğimizi yapıp yapmadığını görebiliriz.

Nitekim, _splitCsv_'den sonra dizilerimiz var ve sonra _map_'ten sonra güzel temiz string'lerimiz var; sadece _"Hello"_, _"Bonjour"_ ve _"Holà"_. Ve pipeline şimdi istediğimizi yapıyor. Fantastik.

Artık tüm bu view komutlarından kurtulabiliriz. Artık onlara ihtiyacımız yok.

## Özet

Hata ayıklamamızı bitirdik ve sonunda elde ettiğimiz kod bu. _input_ adlı CLI parametremizi alıyoruz; bu bir _Path_ olarak sınıflandırılmış. Nextflow yolu bulur, yükler ve CSV dosyasını anlar. Tüm farklı satırları döndürür. Ve sonra o satırın sadece ilk öğesini kanala eşleriz; bu bize kanal içeriğini veren bir tür; bu sürece geçirilir.

Ve süreç kanaldaki her öğe üzerinde çalışır; bu üç. Süreci üç kez çalıştırır; ona üç görev verir. Ve bu sonuçlar daha sonra iş akışından yayımlanır; süreç çıktısı tarafından alınır. Bir iş akışından yayımlanır ve çıktı bloğunda _"hello_channels"_ adlı bir alt dizine kaydedilir.

Oldukça havalı. Şimdi gerçek bir analiz için çalıştırabileceğiniz gerçek hayattaki bir Nextflow pipeline'ına daha yakın bir şeye ulaşıyoruz.

## Özetle

Tamam. Umarım şimdi Nextflow kanallarının ve operatörlerin ne olduğu, operatörlerin kanallar üzerinde nasıl çalıştığı ve bunları nasıl oluşturabileceğiniz hakkında bir fikir ediniyorsunuzdur.

Kanallar, bu videonun başında söylediğim gibi, Nextflow'un yapıştırıcısıdır. Ve burada farklı girdiler alabileceğimizi, bunları manipüle edebileceğimizi ve bu veriyi alıp aşağı akış iş akışı mantığına geçirebileceğimizi görebilirsiniz.

Ve buradaki bu workflow bloğu gerçekten tüm o paralelleştirmeyi ve tüm akıllı mantığı oluşturduğunuz ve Nextflow'a iş akışı DAG'ınızı nasıl oluşturacağını ve pipeline'ınızı nasıl düzenleyeceğini açıkladığınız yerdir.

Kanallar kafanızda oturması en kolay kavram değil. Bu yüzden bir mola verin, bunun hakkında biraz düşünün, belki materyali tekrar okuyun ve bu kavramları gerçekten anladığınızdan emin olun; çünkü bu Nextflow anlayışınızın anahtarıdır. Kanalları ve farklı kanal operatörlerini ve farklı kanal fabrikalarını ne kadar iyi anlarsanız, Nextflow yazarken o kadar çok eğleneceksiniz ve pipeline'larınız o kadar güçlü olacak.

Bu, Python'daki veya diğer dillerdeki normal programlama ile aynı değil. Burada _if_ ifadeleri kullanmıyoruz; bu kanallar ve operatörler kullanarak fonksiyonel akış programlamadır. Yani biraz farklı; ama aynı zamanda süper güçlü.

Bu bölümün sonu. Gidin hızlı bir mola verin ve bir sonraki videoda sizi göreceğim; üçüncü bölüm için Hello Workflow'dan geçeceğiz ve iş akışları hakkında biraz daha konuşacağız.

Önceki bölüm gibi, burada web sayfasının altında birkaç sınav sorusu var; bu yüzden bunları hızlıca gözden geçirebilir ve az önce yaptığımız materyalin tüm farklı bölümlerini anladığınızdan emin olabilirsiniz. Bunun dışında, sizi bir sonraki videoda göreceğim. Çok teşekkür ederim.

Tamam.

​
