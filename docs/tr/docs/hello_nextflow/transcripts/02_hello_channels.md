# Bölüm 2: Hello Channels - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!! note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../02_hello_channels.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba ve Hello Nextflow'un 2. Bölümüne tekrar hoş geldiniz. Bu bölümün adı Hello Channels.

Kanallar, Nextflow boru hattınızdaki yapıştırıcı gibidir. Nextflow'un tüm bilgileri aktarmak ve iş akışınızı düzenlemek için kullandığı, tüm farklı süreçleri bir arada tutan parçalardır.

Kanalların bir başka bölümü de operatörlerdir. Bunlar temel olarak içerikleri değiştirmek için kanallar üzerinde kullanabileceğimiz fonksiyonlardır. Hadi VS Code'a dalıp nerede olduğumuza bakalım.

Bu VS Code'da çok yakınlaştırdım, bu yüzden işleri temiz ve düzenli tutmak için tüm _.nextflow\*_ dosyalarını ve _work/_ dizinini, _results/_ dizinini ve Birinci Bölüm'den her şeyi kaldırdım. Ve burada yeni baştan başlıyorum. Ama bunun için çok endişelenmeyin. İstemiyorsanız, bu dosyaları orada bırakabilirsiniz. Herhangi bir soruna neden olmazlar.

Bu bölüm için _hello-channels.nf_ üzerinde çalışmaya başlayacağız ve bunu açarsam, daha önce üzerinde çalıştığımız dosyaya çok benzer görünmeli. Farklı bölümlerin betiğin farklı yerlerinde olması mümkün, ancak her şey temelde aynı olmalı.

Farklı olan bir şey, buradaki çıktı bloğundaki yolun artık bu bölüm için _hello_channels_ olması, bu da sonuç dosyalarının hala orada varsa sonuçlarınızda farklı bir alt dizinde saklanacağı anlamına gelir. Bu yüzden çıktılar hakkında kafanız karışmadan başlamak için güzel ve temiz bir yer olmalı.

Tamam, bu betiği çalıştırdığımızda ne yaptığını hızlıca hatırlayalım. _"nextflow run hello-channels.nf"_ yapıyoruz. _"--input myinput"_ yapabiliriz ve bunu çalıştırdığımızda, yukarıdaki sayHello sürecinin değişkeni olan bu params.input parametresini kullanacak, bu da greeting'e gider ve output.txt'ye kaydedilir. Ve bunu sonuç dosyasında görebiliriz. Harika.

## 1. Değişken girdileri açıkça bir kanal aracılığıyla sağlayın

Bu güzel. Ama oldukça basit. Bu parametrede bir kez çalışan bir sürece giden bir değişkenimiz var ve gerçekten ölçeklenmiyor. Ve burada oluşturmak için ona birçok farklı dosya veremiyoruz. Ona birçok farklı selamlama veremiyoruz. Sadece bir tane var.

Gerçekte, Nextflow tamamen analizinizi ölçeklendirmekle ilgilidir. Bu yüzden muhtemelen birden fazla şey yapmasını istersiniz. Ve bunu _kanallarla_ yapıyoruz.

Kanallar, Nextflow'u öğrenen birçok kişi için biraz benzersiz bir kavramdır. Fonksiyonel programlama kavramlarından gelir ve kafanızda oturması biraz zaman alabilir, ancak bir kez anladığınızda, Nextflow'un gücünü gerçekten açığa çıkarırlar ve iş akışlarınızı nasıl yazacağınızın anahtarıdır.

## 1.1. Bir girdi kanalı oluşturun

Bu betiği alıp sadece bir _param_ yerine bir _kanal_ kullanmasını sağlayarak başlayalım.

İşleri bir araya getirme konusundaki tüm iş akışı mantığımızın bulunduğu workflow'a gidiyoruz. Ve buraya gireceğim ve yeni bir kanal oluşturacağım.

Yeni bir kanal oluştur.

Ve buna "_greeting_ch"_ diyeceğim. Bu değişkenin bir kanal olduğunu hatırlayabilmeniz için "_\_ch"_ yapmak bir gelenektir. Ama istediğiniz gibi adlandırabilirsiniz.

Ve sonra eşittir diyeceğim ve _"Channel.of"_ yapacağım.

Channel, kanallarla ilgili her şey için ad alanı gibidir. Daha önce Nextflow kullanıyorsanız küçük harf "c". Ve _".of"_, kanal oluşturmanın bir yolu olan Channel factory denen bir şeydir.

Birçok farklı kanal fabrikası var. Burada sadece "." yaparsam, VS Code'un bir sürü öneri yaptığını görebilirsiniz, ancak _".of"_ en basit olanıdır ve sadece burada bir girdi alır.

Yani bazı parantezler yapabilirim ve _"Hello Channels!"_ diyeceğim.

Harika. Bir kanalım var. Fantastik. Kaydet'e basabilirim, tekrar çalıştırabilirim, ama ilginç bir şey olmayacak. VS Code burada bana turuncu bir uyarı çizgisi verdi ve bunun ayarlandığını söyledi: bunu oluşturdunuz, ama aslında hiçbir şey için kullanmadınız. Bu kanal tüketilmiyor.

Tamam, peki bunu nasıl kullanırız? Çok basit. Bunu alacağım, kopyalayacağım ve _params.input_'u sileceğim ve bunun yerine _"greeting_ch"_ koyacağım. Yani bu kanalı sayHello'ya girdi olarak geçireceğiz.

Şimdilik bu dizeyi sabit kodladığımı unutmayın. Bu, son bölümün sonunda kullandığımız güzel parametreden sonra biraz geri adım, ama mantığı görebilmeniz için işleri basit tutuyor.

Tamam, terminalime gideceğim ve iş akışını tekrar çalıştıracağım. Bu sefer herhangi bir _"--input"_ olmadan, ve çalışacak ve oluşturduğumuz bu kanalı kullanacak ve umarım _results/hello_channels/_ içinde bir dosyamız olmalı ve şimdi "Hello Channels!" diyor. Fantastik. Yani kanalımızdan umduğumuz şey bu. Harika.

## 1.4. Kanal içeriğini incelemek için view() kullanın

Burada eklenecek bir şey daha var, kanallarda kullanabileceğimiz "_.view"_ adlı başka bir fonksiyona hızlı bir giriş.

Bu, alışık olabileceğiniz Python veya diğer dillerdeki _print_ komutuna benzer ve çalıştırdığımızda bu kanalın içeriğini terminale döker.

Yani "_.view"_ yapın ve sonra iş akışını tekrar çalıştırırsam, bu kanalın içeriğinin ne olduğunu, oluşturduğumuz zamandaki halini terminale yazdırmalı.

Nitekim, terminale buraya yazdırıldığını görebilirsiniz. _"Hello Channels!"_.

İsterseniz bunları satırlar arasında bölebileceğinizi unutmayın ve aslında Nextflow otomatik biçimlendiricisi bunu sizin için yapmaya çalışacaktır. Boşluk burada gerçekten önemli değil, bu yüzden bunları birbiri ardına zincirleyebilirsiniz.

## 2. İş akışını birden fazla girdi değeri üzerinde çalışacak şekilde değiştirin

Tamam, yani kanalımızda güzel olan bir şey var, ama temelde daha öncekiyle aynı. O halde biraz daha karmaşık hale getirelim. Kanalımıza birkaç şey daha ekleyelim.

"_.of()"_ kanal fabrikası birden fazla öğe alabilir, o halde birkaç tane daha yazalım. _Hello, Bonjour, Hej_ yapacağız. Ve sonra bu iş akışını tekrar çalıştırabiliriz ve ne olacağını göreceğiz.

Tekrar çalışmalı. Ve şimdi yazdırdık. View ifademizle terminale _"Hello", "Bonjour"_ ve _"Hej"_. Fantastik.

## 2.1.2. Komutu çalıştırın ve log çıktısına bakın

Bu noktada işimizin bittiğini düşünebilirsiniz. Ama aslında burada bizi yakalayacak küçük bir tuzak var. Çıktı dosyamıza bakarsak. İçinde _"Hello"_ olduğunu görebilirsiniz, ama diğer çıktılardan hiçbiri yok. Aslında, sadece bu.

Bu iş akışını birden çok kez çalıştırırsak, bazen _"Bonjour"_ olduğunu, bazen _"Hej"_ olduğunu bile görebiliriz. Biraz rastgele.

Terminale bakarsak, üç kez çalıştığını ve farklı view çıktılarını görebiliriz. Ama work dizinine gidersem, _"cat work"_ yapabilirim. Bu hash'i koyun ve genişletin ve _output.txt_. Bu dosyanın work dizinindeki dosyanın results dizininden farklı olduğunu görebilirsiniz ve bu _"Hej"_. Yani burada tam olarak çalışmayan bir şey var.

Ve anahtar şu ki, üç görev çalıştı. Nextflow çıktısı, işlem devam ederken bunu özetlemeye çalışır, böylece terminalinizin tamamını tamamen ele geçirmez ve bu ANSI Logging, ANSI kaçış kodlarını kullanır, temelde diğer görevlerin üzerine yazmıştır. Bu yüzden size güncellenen son olanı gösterir.

## 2.1.3. Komutu -ansi-log false seçeneğiyle tekrar çalıştırın

Bunu gerçekten daha iyi anlamak için yapabileceğimiz birkaç şey var. Work dizininin kendisine bakabiliriz ve oradaki tüm farklı work dizinlerini görebilirsiniz, ama bu biraz kafa karıştırıcı çünkü farklı Nextflow yürütme çalıştırmalarıyla karışacak.

Ya da Nextflow'a ANSI kaçış kodlarını kullanmamasını söyleyebiliriz.

Yani komutu tekrar çalıştırırsam, ama bu sefer _"-ansi-log false"_ diyerek kapatırsam, _$NO_COLOR_ veya _"$NXF_ANSI_LOG=false"_ ortam değişkenlerini de kullanabilirim. O zaman bu kaçış kodları olmadan Nextflow günlüğünün daha eski tarzını kullanır. Akıllı güncellemeler olmadan doğrudan terminale yazdırır.

Ve şimdi çalışan bu üç sürecin tümünü görebiliriz. Ve her birinin kendi görev hash'i var. Ve bu work dizinlerine girersek, belirttiğimiz üç farklı selamlamayı göreceğiz.

Yani bu artık biraz daha mantıklı. Umarım Nextflow'un bunu yaptığını anlarsınız, sadece bu work dizinleriyle terminalde size gösterdiği konuda biraz akıllıydı.

Ancak, bu work dizinleriyle ilgili bir sorunu çözdü, ama çıktı dosyasıyla ilgili bir sorunu çözmedi. Hala sadece _"Hello"_ yazan bir çıktı dosyamız var.

## 2.2. Çıktı dosya adlarının benzersiz olacağından emin olun

Şimdi bunu anlamak için iş akışı betiğimize geri dönmemiz gerekiyor. Kanalımızı burada oluşturuyoruz, sürecimize geçiriyoruz ve sürece bakarsak, selamlamayı _"output.txt"_ adlı bir dosyaya yazıyoruz ve bu çıktı dosyasını buradaki çıktı bloğuna geri geçiriyoruz, yayınlıyoruz.

Ancak, bu süreç her üç kez çalıştığında bu üç farklı görev. Hepsi _"output.txt"_ adlı bir dosya oluşturur, tüm bu çıktı dosyaları results dizinine yayınlanır ve hepsi birbirinin üzerine yazar. Yani orada aldığınız sonuç dosyası her neyse, sadece oluşturulan son dosyadır, ama diğerlerinin hepsini silmiştir. Bu gerçekten istediğimiz şey değil.

## 2.2.1. Dinamik bir çıktı dosya adı oluşturun

Bunu ele almanın farklı yolları var, ama şimdilik en basiti sadece farklı benzersiz dosya adları oluşturmak. Böylece görev her farklı selamlama ile çalıştığında, farklı bir çıktı dosyası oluşturacak, bu da yayınlandığında artık çakışmayacak. Ve sonra üç benzersiz çıktı dosyası alacağız.

Bunu tamamen aynı şekilde yapıyoruz. Bu değişkeni script bloğu içinde herhangi bir yerde kullanabiliriz ve birden çok kez kullanabiliriz.

Yani bunu buraya yapıştırabilirim, _"$\{greeting\}\_output.txt"_, ve sonra bunu buraya da yapıştırmam gerekiyor çünkü artık _output.txt_ adlı bir dosya oluşturmuyoruz. Yani bunu güncellemezsen, Nextflow hiç oluşturulmamış bir dosya beklediğini söyleyen bir hatayla çökecek.

Yani orada da aynısını yapmam gerekiyor ve bu değişkenin anlaşılması için tek tırnak değil çift tırnak kullanmam gerekiyor.

Tamam, deneyelim ve işe yarayıp yaramadığına bakalım. İş akışını tekrar çalıştıracağız. Umarım bize üç farklı work dizinindeki üç farklı görevi gösterecek. Ve nitekim, burada solda results klasöründe görebilirsiniz. Şimdi üç farklı dosya adıyla üç farklı dosyamız var ve her biri beklediğimiz farklı içeriklerle. Yani dosyalar artık birbirinin üzerine yazmıyor ve her şey beklediğimiz gibi orada.

Burada geçtiğimiz bu biraz önemsiz bir kurulum, ama dosya yayınlamanın nasıl çalıştığı hakkında anlamanız gereken bazı temel kavramların altını çiziyor ve tuzak olarak düşebileceğiniz bazı şeyler. Umarım bunu kendi iş akışlarınızda önleyebilirsiniz.

Burada yaptığımız şeyin gerçek hayat durumlarında biraz pratik olmadığını da belirtmekte fayda var. Bazı girdi verilerini aldık ve bu veriyi kullanıyoruz, ama aynı zamanda dosyayı bu veriden sonra adlandırıyoruz, ki bunu genellikle yapamazsınız.

Yani gerçek daha olgun Nextflow boru hatlarında, genellikle belirli bir örnekle ilişkili tüm meta verileri içeren bir meta nesnesi geçirirsiniz. Daha sonra buna dayalı olarak dinamik dosya adları oluşturabilirsiniz, bu çok daha pratiktir.

Bunu en iyi uygulamalarla nasıl yapacağınızla ilgileniyorsanız, _training.nextflow.io_'da özellikle meta veriler ve meta haritalar hakkında bir yan görev var, bu yüzden daha fazla ayrıntı için oraya bakabilirsiniz.

## 3. Bir dizi aracılığıyla birden fazla girdi sağlayın

Tamam. Şimdi kanalların nasıl yapılandırıldığı ve kodlama dilindeki diğer veri yapılarından nasıl farklı olduğu hakkında biraz keşfedeceğiz. Ve diğer dillerden geldiyseniz tanıdık bir kavram olabilecek bir diziyi nasıl kullanabileceğimi düşüneceğim.

Bir kanalda bir dizi kullanabilir miyim? Deneyelim. Bir dizi oluşturacağım ve bunu dokümanlardan kopyaladım, _"greetings_array"_ ve _"Hello", "Bonjour"_ ve _"Holà"_. Ve sonra bunu sabit kodlanmış dizelerim yerine buraya koyacağım. Yani "Channel.of" _"greetings_array"_ diyeceğim, bu diziyi bir kanala geçireceğim. Deneyelim.

Terminali açın ve boru hattını çalıştırın.

Tamam. View ifadesinin dizimizi beklendiği gibi yazdırdığını görebilirsiniz, ama sonra tüm bu kırmızı metin, ya da hala _"-ansi-log"_ kapalıysa kırmızı olmayacak, ama tüm bu kırmızı metin bize bir şeylerin yanlış gittiğini söylüyor.

Artık burada güzel bir yeşil onay işaretimiz yok. Kırmızı bir çarpımız var ve bunu biraz daha geniş yaparsam okumak daha kolay olur, Nextflow bize neyin yanlış gittiğini söylüyor.

Hadi bunu bölüm bölüm inceleyelim. Hatanın nedeni olduğunu söylüyor ve sonra hatanın nedeni, eksik çıktı dosyaları. Yani temelde o çıktı bloğu bu dosyanın oluşturulması gerektiğini söyledi ve oluşturulmadı. Sonra bu yürütülen komut olduğunu söylüyor. Yani bu temelde o _.command.sh_ dosyasının içeriği. Tüm bu değişkenler konulduktan sonra böyle görünüyordu.

Ve burada echo komutumuzun aslında sadece bir kez çalıştırıldığını ve tüm diziyi kullandığını görebilirsiniz, ama bir dize gösteriminde, ki bu gerçekten istediğimiz şey değildi.

Ve sonra komut böyle çıktı ve bu dosyaları görmek için gidebileceğimiz work diziniydi, biraz daha anlamak için.

Tamam. O zaman olan şey şuydu. Nextflow bu tüm diziyi sürece tek bir kanal öğesi olarak geçirdi, bu da sürecin sadece bir kez çalıştığı anlamına geliyordu. Bir görevi vardı ve veriyi beklediğimiz yapıda kullanmadı.

## 3.2. Kanal içeriğini dönüştürmek için bir operatör kullanın

Yani kullanılmadan önce bu kanala bir şey yapmamız gerekiyor. Ve bu, kanal içeriklerini manipüle etmek için kullanabileceğimiz özel fonksiyonlar olan operatörleri kullanmak için sahneyi hazırlıyor.

Bu durumda, _flatten_ adında bir şey kullanacağız. Bunu burada kanalın sonuna geçiriyoruz. Yani kanalı oluşturuyoruz ve sonra _flatten_ çalıştırıyoruz. Ve yine, üzerine gelirseniz, VS Code'da bu komutun belgelerini hemen gösterir, bu çok yardımcı. Tüm bu belgeleri Nextflow web sitesinde, belgelerde de bulabilirsiniz.

Şimdi bu kodu çalıştırabilir ve işe yarayıp yaramadığını görebilirim, ama aynı zamanda operatörler içinde ve Nextflow kodu içinde dinamik kod yapmanın nasıl yapılacağını tanıtmak için güzel bir fırsat, bunlara closure denir.

Yani _flatten_ çalıştırmadan önce buraya bir view komutu ekleyeceğim. Ve burada bu süslü parantezler var, bu dinamik closure. Ve burada sadece view operatörü bağlamında yürütülecek bazı rastgele kodlar var.

Burada, view operatörünün girdileri olan greeting'i al diyor ve bu burada. Buna istediğim gibi isim verebilirdim, buna _"foo"_ diyebilirdim ve sadece daha sonra _"foo"_ olarak ona atıfta bulunmam gerekir. Ve sonra bununla, bunu döndür diyorum.

Ve sonra bir değişken için flatten'dan önce yazan bir dize döndürüyorum. çok basit.

Şimdi bunun tamamen aynısını ekleyeceğim, ama _flatten_'dan sonra diyeceğim.

Yani bunun yaptığı şey, bu sırayla çalıştığı için, _flatten_ çalıştırmadan önce kanalın nasıl göründüğünü göreceksiniz ve sonra _flatten_ çalıştırdıktan sonra tekrar.

Ve sonra bu greeting kanalı hala oluşturuldu, bu yüzden hala sürece geçirilecek. Ve umarım şimdi iş akışı çalışacak. Deneyelim.

Harika. İlk olarak, bu sefer boru hattı çökmedi. Düzgün çalışan üç sürecimiz vardı ve küçük bir onay işaretimiz var. Ve sonra view ifadelerimizin çalıştığını görebiliriz.

_flatten_'dan önce var, bu daha önce başarısızlıktan gördüğümüz dizi, ve sonra _flatten_ çağrıldıktan sonra üç kez var, _"Hello", "Bonjour"_ ve dizideki diğer üç ayrı öğe var, bunlar şimdi umduğumuz gibi kanaldaki üç ayrı öğe.

Ve _view_ operatörünün üç kez çalıştırıldığını görebilirsiniz. Ve bunun nedeni _flatten_'dan sonra bu kanalın artık üç öğesi olması. Ve bu yüzden operatör üç kez çağrılır.

Çok hızlı bir şekilde, daha önce kanal fabrikaları oluştururken _"."_ yaptığımı söylemek isterim ve kanallar oluşturmanın birçok farklı yolu olduğunu gördük ve bunlardan birine "_fromList"_ denir. Ve bu aslında özellikle bu aynı işlemi yapmak için tasarlanmıştır. Yani sadece from list greetings away yapabilirdik ve bu işe yarayacaktır. Biraz daha temiz ve güzel bir sözdizimi. Ama bu gösterimin amaçları için, kanalın nasıl manipüle edildiğini ve farklı operatörlerin bir kanalın içeriğindeki içeriği nasıl değiştirebileceğini görebilmeniz için adım adım yapmak istedik.

## 4. Girdi değerlerini bir CSV dosyasından okuyun

Tamam, bunu biraz daha gerçekçi nasıl yapabiliriz? Muhtemelen Nextflow boru hattınızda sabit kodlanmış dizilerle bir sürü kod oluşturmak istemeyeceksiniz. Muhtemelen başlattığınızda veriyi dışarıdan almak isteyeceksiniz ve bu veri neredeyse kesinlikle dosyalarda olacak.

Yani yapacağımız bir sonraki şey, bunu çoğaltacağız, ama veriyi tek bir CLI parametresinden veya sabit kodlanmış bir dizeden veya diziden almak yerine, bir dosyadan alacağız.

Hadi greetings away'den kurtulalım. Ve şimdi bu kanal fabrikasını tekrar değiştireceğiz. Seçim yapabileceğiniz bir sürü olduğunu söyledim ve _".fromPath"_ adında bir tane var. Ve ona söyleyeceğim, bu durumda, daha önce kullandığımız girdimiz olan _params.input_'u alsın.

Şimdi bu parametre henüz kullanılmaya hazır değil. Hala bunun bir dize olduğunu söylüyoruz ve burada bir varsayılanla sabit kodlanmış, ama bu dizeyi geçersiz kılabiliriz. Şimdi bunun bir dosya olmasını istiyoruz. Yani tür farklı. Artık bir _String_ değil. Bir _Path_.

Ve sonra istediğimiz gibi varsayılanı tekrar bir Path'e ayarlayabiliriz. Ve soldaki keşfet'e bakarsam, bu depoda, bu çalışma dizininde, data adında bir dizinimin olduğunu görebilirsiniz. Orada _"greetings.csv"_ adında bir dosyam var.

Yani varsayılanı burada _"data/greetings.csv"_ olarak ayarlayabilirim. Şimdi, bu boru hattını herhangi bir komut satırı seçeneği olmadan tekrar çalıştırdığımda, bu varsayılan değeri kullanacak. Bunun bir yol olduğunu biliyor, bu yüzden bunu bir dize olarak değil bir yol olarak ele alması gerektiğini biliyor.

Ve sonra bunu bu _params.input_'tan bir kanal fabrikasına geçirecek ve kanalımızı oluşturacak, bu da _sayHello_ adlı bu süreçte kullanılacak. Deneyelim.

Tamam. Başarısız oldu. Endişelenmeyin. Bu bekleniyordu. Ve eğitim materyalini takip ediyorsanız, orada da beklendiğini göreceksiniz. Burada ne oluyor bakalım.

Boru hattını çalıştırmaya çalıştı. Süreci yürütmeye çalıştı ve daha önce gördüğümüze oldukça benzer bir hata aldı.

Burada diyor ki: \_echo çalıştırmaya çalıştık, \_ama bu CSV dosyasının içeriğini yankılamak yerine, sadece yolu yankıladı. Ve burada bu CSV dosyasına giden tam mutlak yolu görebilirsiniz.

Ve sonra elbette, bunu bu gerçekten karmaşık yola yazmaya çalıştığı için, ne yapacağını gerçekten bilmiyordu. Ve süreç work dizininin kapsamı dışındaydı.

Başlangıçta Nextflow'un her yürütülen görevi özel bir work dizini içinde kapsüllediğini söylemiştim. Ve bu work dizininin dışında olan verilere yazmaya çalışırsanız, Nextflow sizi bir güvenlik önlemi olarak durduracaktır. Ve burada olan da bu. Mutlak bir yola yazmaya çalıştık ve Nextflow başarısız oldu ve bizi engelledi.

## 4.2. Dosyayı ayrıştırmak için splitCsv() operatörünü kullanın

Tamam, bu kanala bir bakalım ve nasıl göründüğünü görelim. _".view"_ yapabiliriz ve bunu web sitesinden kopyaladım. Yani _.view_ ve burada dinamik bir closure'ımız var ve girdi olarak "_csv"_ değişken adını söylüyoruz. Yani bu kanal içerikleri ve splitCsv'den önce diyoruz ve böyle görünüyor.

Tekrar çalıştırırsam, yine başarısız olacak, ama bu kanalın içinde ne olduğunu gösterecek. Özellikle heyecan verici değil. Bu _path_ değişkeni. Yani burada bir terminale yazdırıldığı için sadece bir dize olduğunu görebilirsiniz, ama bu dosya hakkındaki bilgileri ve meta verileri içeren bir _path_ nesnesi.

Dosyanın meta verilerini girdiye geçirmek istemiyoruz. Bu dosyanın içeriğini geçirmek istiyoruz. _greetings.csv_ dosyasına bakarsak, burada bu farklı değişkenlerin olduğunu görebilirsiniz. Yine _Hello, Bonjour, Holà_. Ve bunlar gerçekten sürecimize geçirmek istediğimiz şeyler, sadece dosyanın kendisi tek bir nesne olarak değil.

Yani bu CSV dosyasını ayrıştırmamız gerekiyor. Onu açmamız, CSV dosyasının içeriğine ulaşmamız ve sonra içerikleri kanal içinde sürece geçirmemiz gerekiyor.

Log mesajından muhtemelen anlayabileceğiniz gibi, başka bir operatör, başka bir kanal operatörü olan _splitCsv_'yi kullanmak istiyoruz. Yani "_dot" "s"_ yaparsam ve sonra otomatik önerildiğini görebilirsiniz. Hoop, _splitCsv_ ve bazı parantezler.

Ve sonra _splitCsv_'den sonra, başka bir _view_ ifadesi koyacağım, böylece sonrasında nasıl göründüğünü görebiliriz. Boru hattını çalıştıralım ve ne aldığımızı görelim.

Tamam. Hala başarısız oldu, ama yeni ve heyecan verici bir şekilde, bu ilerleme.

Bu sefer yine betiğimizle ilgili bir sorunumuz var, bu render edildi. Şimdi. Artık son yolu almadık, ama daha önce bir diziyi sabit girdi olarak geçirdiğimizde yaşadığımız hataya çok benzeyen bir değişkenler dizisi aldık.

View operatöründen gelen günlüğümüzle, _splitCsv_'den önce yolun olduğunu görebiliriz. Ve nitekim, _splitCsv_'den sonra, üç farklı çıktımız var ve bu çıktıların her biri _greetings.csv_ dosyasındaki satırların her birine çok benziyor, bu mantıklı.

Yani burada olan şey, Nextflow'un bu CSV dosyasını ayrıştırması, bize CSV dosyasının her satırı için bir dizi olan üç nesne vermesi. Yani sonra üç kez bir dize değeri yerine bir değişkenler dizisini kanala geçirdik.

Tamam, geçen sefer bu sorunumuz olduğunda, _flatten_ kullandık. Hadi çok hızlı bir şekilde. Flatten'ı deneyelim ve ne olacağını görelim.

Bu değişkenleri ne olursa olsun çağırabilirim. Yani buna _myarray_ diyeceğim çünkü artık gerçekten bir CSV değil. Tekrar çalıştırmayı deneyelim ve _flatten_ ile ne olacağını görelim.

Yani bu sefer çalışacağız, CSV'yi üç dizi nesnesine ayrıştırdık ve sonra düzleştirdik. Ve bu sefer geçti. Ve Nextflow boru hattı çalıştı. Ancak _flatten_'ın gerçekten işe koyulduğunu ve her şeyi düzleştirdiğini görebilirsiniz. Ve böylece her satır için üç bağımsız dizi girişi alıyoruz. Ve böylece süreci her CSV satırı için üç kez çalıştırdı. Ve şimdi bir sürü sonuç dosyamız var ve 123, 456 ve her türlü şey, sadece gerçekten istediğimiz CSV'nin ilk sütunu değil.

## 4.3. Selamlamaları çıkarmak için map() operatörünü kullanın

Peki sadece ilk sütuna nasıl ulaşırız? Flatten burada çok basitse, aslında özelleştirebileceğimiz ve CSV'den ne istediğimizi söyleyebileceğimiz daha karmaşık bir operatöre ihtiyacımız var.

Bunu yapmak için _map_ kullanacağız. Temel olarak _map_ sadece diyor ki, bana verilen her öğe üzerinde bir kod, bir fonksiyon çalıştır ve üzerinde bir tür dönüşüm yap. Ve çok esnek olduğu için, Nextflow kodunda her zaman karşınıza çıktığını göreceksiniz.

Kendi başına hiçbir şey yapmaz. Yani normal parantez istemiyoruz, burada bir closure istiyoruz ve ne yapacağını söylememiz gerekiyor. Yani _"row"_ diyeceğim, çünkü CSV'den satırlar veriliyor, bu yüzden mantıklı bir değişken adı. Girdi. Ve sadece o dizinin ilk öğesini döndürmek istiyorum.

Nextflow'daki diziler sıfır tabanlıdır, bu yüzden sadece ilk öğeyi söyleyeceğiz, bu satır sıfır. İkinci sütunu isteseydik, bir olabilir veya üçüncü sütun iki olabilir, ve böyle devam eder. Burada istediğimizi döndürebiliriz, ama sadece ilk değeri döndüreceğim.

Ve şimdi, boru hattını tekrar çalıştırabiliriz ve beklediğimizi yapıp yapmadığını görebiliriz.

Nitekim, _splitCsv_'den sonra dizilerimizi aldık ve sonra _map_'ten sonra güzel temiz dizelerimiz var, sadece _"Hello", "Bonjour"_ ve _"Holà"_. Ve boru hattı şimdi istediğimizi yapıyor. Fantastik.

Yani artık tüm bu view komutlarından kurtulabiliriz. Artık onlara ihtiyacımız yok.

## Özet

Hata ayıklamamızı bitirdik ve sonunda elde ettiğimiz kod bu. _input_ adlı CLI parametremizi alıyoruz, bu bir _Path_ olarak sınıflandırılmış. Nextflow yolu bulur, yükler ve CSV dosyasını anlar. Tüm farklı satırları döndürür. Ve sonra o satırın sadece ilk öğesini kanala eşleriz, bu da bize kanal içeriklerini veren, sürece geçirilen şey.

Ve süreç kanaldaki her öğe üzerinde çalışır, bu üç. Ve süreci üç kez çalıştırır, ona üç görev verir. Ve bu sonuçlar daha sonra iş akışından yayınlanır, süreç çıktısı tarafından alınır. Bir iş akışından yayınlanır ve çıktı bloğunda _"hello_channels"_ adlı bir alt dizine kaydedilir.

Oldukça havalı. Şimdi gerçek bir analiz için çalıştırabileceğiniz gerçek hayattaki bir Nextflow boru hattına daha yakın bir şeye ulaşıyoruz.

## Özet

Tamam. Umarım şimdi Nextflow kanallarının ve operatörlerin ne olduğu ve operatörlerin kanallar üzerinde nasıl çalıştığı ve bunları nasıl oluşturabileceğiniz hakkında bir fikir ediniyorsunuzdur.

Bu videonun başında söylediğim gibi kanallar, Nextflow'un yapıştırıcısıdır. Ve burada farklı girdiler alabileceğimizi ve bunları manipüle edebileceğimizi ve bu veriyi alıp sonra aşağı akış iş akışı mantığına geçirebileceğimizi görebilirsiniz.

Ve buradaki bu workflow bloğu gerçekten tüm o paralelleştirmeyi ve tüm akıllı mantığı oluşturduğunuz yerdir ve Nextflow'a iş akışı DAG'ınızı nasıl oluşturacağınızı ve boru hattınızı nasıl düzenleyeceğinizi açıklarsınız.

Kanallar kafanızda oturması en kolay kavram değil. Bu yüzden bir mola verin, bunun hakkında biraz düşünün, belki materyali tekrar okuyun ve bu kavramları gerçekten anladığınızdan emin olun çünkü bu Nextflow anlayışınızın anahtarıdır ve kanalları ve farklı kanal operatörlerini ve farklı kanal fabrikalarını ne kadar iyi anlarsanız. Nextflow yazarken o kadar çok eğleneceksiniz ve boru hatlarınız o kadar güçlü olacak.

Bu Python veya diğer dillerde normal programlama ile aynı değil. Burada _if_ ifadeleri kullanmıyoruz, bu kanallar ve operatörler kullanarak fonksiyonel akış programlamasıdır. Yani biraz farklı, ama aynı zamanda süper güçlü.

Bu bölümün sonu. Gidin hızlı bir mola verin ve bir sonraki videoda sizi göreceğim, üçüncü bölüm için Hello Workflow'dan geçeceğiz ve iş akışları hakkında biraz daha konuşacağız.

Önceki bölüm gibi, burada web sayfasının altında birkaç sınav sorusu var, bu yüzden bunları hızlıca gözden geçirebilir ve az önce yaptığımız materyalin tüm farklı bölümlerini anladığınızdan emin olabilirsiniz. Ve bunun dışında, bir sonraki videoda görüşürüz. Çok teşekkür ederim.

Tamam.
