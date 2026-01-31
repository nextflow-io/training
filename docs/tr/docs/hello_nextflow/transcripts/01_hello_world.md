# Bölüm 1: Merhaba Dünya - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../01_hello_world.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş geldiniz

Merhaba, Merhaba Nextflow'un Birinci Bölümüne hoş geldiniz.

Altı bölümden oluşan bu kursun ilk bölümünde, Nextflow'un en temel konularına gireceğiz. Terminalde bazı komutlar çalıştırarak başlayacağız ve ardından bu Bash komutlarını alıp bunları bir Nextflow betiğine nasıl dönüştüreceğimizi göreceğiz.

İlk Nextflow işlem hattımızı çalıştırmayı deneyeceğiz, Nextflow'un ne yaptığını, nerede çalıştığını, hangi dosyaları oluşturduğunu ve bu dosyaların amacının ne olduğunu göreceğiz.

Pekala, hadi başlayalım.

## training.nextflow.io

Öncelikle, training.nextflow.io adresine gidin. Daha önce olduğu gibi, tüm materyaller burada yazılı ve ben bunları adım adım çalışacağım. Eğitimin adımlarını yaparken ekranımı göstereceğim, ancak söylediğim her şey eğitim materyalinde yazılı, böylece kendi hızınızda takip edebilir ve hepsini orada yazılı bulabilirsiniz.

Bu videonun ayrıca video altyazıları etkinleştirilmiş, bu yüzden bunları açıp söylediklerimi tam olarak takip edebilirsiniz.

Tamam, hadi Merhaba Nextflow'a gidelim. Bugün yapacağımız kurs bu ve oryantasyonu ilk videoda zaten yaptık, bu yüzden doğrudan birinci bölüme geçeceğiz. Merhaba Dünya.

Tamam, şimdi bu eğitim materyalinden ayrılıp Code Spaces ortamıma geçeceğim. Bu, ilk videoda kurduğumuz şey. Umarım sizin sisteminizde de buna çok benzer bir şey vardır. VS Code kullanıyorum ve eğitim materyaline bakıyorum ve hello Nextflow dizinine geçiş yaptım.

## 0. Isınma: Merhaba Dünya'yı doğrudan çalıştırın

Tamam. Herkese tanıdık gelecek birkaç temel şeyle başlayalım. Terminalde çok basit bir komut yazarak başlayacağım. Burada 'echo Merhaba Dünya!' yazacağım, enter'a basacağım ve sürpriz yok, terminal istediğimi yapıyor ve o dizgiyi döndürüyor. Merhaba dünya.

Tamam, sonra o komutu geri almak için yukarı tuşuna basacağım ve biraz daha düzenleyeceğim. Bu sefer çıktıyı bir dosyaya yönlendirelim. Bunun yerine output.txt dosyasına yazacağım ve enter'a basacağım, bu sefer terminalde hiçbir şey yok çünkü çıktı terminale gelmedi. Dosyaya gitti.

Sonra o dosyayı 'cat output.txt' yaparak okuyabilirim, dosya adını otomatik genişletmek için tab tuşuna basıyorum ve işte. Dosya orada.

Bu dosyayı VS Code'daki kenar çubuğunda, dosya gezgininde de görebilirim. Çift tıklayıp burada açabilirim. Hiçbir şeye tıklamadan VS Code'da açmak isterseniz, "code" ve ardından "output.txt" de yapabilirsiniz ve aynı şeyi yapar.

Harika. Bu ilk adım. Çok basit.

## 1. Merhaba Dünya iş akışı başlangıç betiğini inceleyin

Tamam. Şimdi aynı şeyi yapacağız, ancak doğrudan terminalde değil, Nextflow'da.

Başlangıç için ilk örnek betiği kullanacağız, bu dosyanın adı Merhaba Dünya. Terminalde görmek için "ls" yapabilirim ve Mac'teyim, bu yüzden o dosyayı açmak için command tuşuyla tıklayabilirim veya kenar çubuğunda burada çift tıklayabilirdim.

Bu dosyada görebileceğimiz birkaç şey var. En üstte, bunun bir Nextflow dosyası olduğunu ve bu şekilde çalıştırılabileceğini söyleyen bir hash ifadesi var. Burada açık gri renkte, yürütmeyi etkilemeyen ve sadece betiği okumamıza yardımcı olan bazı yorumlar var.

Ve sonra iki ana yapı var. Burada bir process ve bir workflow var.

Nextflow'daki process'ler, işlem hattının adımlarıdır. Aslında mantığı uygulayan ve işlemi yapan parçalardır.

Alttaki workflow ise bu process'leri birbirine bağlar ve iş akışının mantığını, her şeyin birbirine nasıl bağlandığını yönetir.

Bir process'e bakarak başlayacağız. Birazdan workflow'a geri döneceğiz.

## 1.2 Process tanımı

Her process, process anahtar kelimesiyle başlar. Bir adı vardır ve sonra bazı süslü parantezler vardır ve bu süslü parantezlerin içindeki her şey o tek process'tir.

Bir process'in script bölümü olmalıdır ve burada, hesaplama ortamında gerçekten yürütülen kodun parçası olan çok satırlı bir dize içinde bir bash kod parçası bulunur.

Burada ayrıca Nextflow'a betik tarafından hangi dosyaların oluşturulmasının beklendiğini söyleyen bir output ifadesi var. Buradaki çıktının, Nextflow'a bunun bir değer veya dize değil, bir dosya olduğunu söyleyen path anahtar kelimesine sahip olduğunu unutmayın.

Script bloğu içinde, bu sadece normal bir bash ifadesidir ve terminalde yazdığımızla tamamen aynıdır. output.txt adlı bir dosyaya merhaba dünya yazdırıyoruz. Bu output.txt daha sonra output tanımı tarafından yakalanır. Output tanımı aslında hiçbir şey yapmıyor. Sadece Nextflow'a neyi bekleyeceğini söylüyor ve bu dosya oluşturulmazsa, Nextflow bir hata verecektir.

Bu örneğin harika olmadığını unutmayın çünkü dosya adını burada sabit kodladık, output.txt ve output.txt. Bunlardan biri değiştirilirse, bu iş akışımızda bir hataya neden olur.

Bunu değişkenlerle yapmanın daha iyi bir yolu var, bunu birazdan ele alacağız.

## 1.3 Workflow tanımı

Tamam. Workflow'a geçersek, bir yorumum olduğunu ve sonra sayHello adlı process'i çalıştırdığımızı görebiliriz. Bu, yukarıdaki ile aynı anahtar kelimedir. Bu, bir iş akışının olabileceği kadar basittir. Değişken girdi olmadan tek bir process çağırıyoruz, bu yüzden onu başka bir şeye bağlamıyoruz. Bu kursun ilerleyen bölümlerinde, değişken girdiler kullanarak ve kanallarla bağlantı kurarak bunu nasıl daha güçlü hale getireceğimizden bahsedeceğiz.

## 2. İş akışını çalıştırın

Tamam, ihtiyacımız olan tek şey bu. Hadi çalıştırıp ne olacağını görelim. Terminali temizleyeceğim ve sonra "nextflow run" yapacağım ve dosya adını çağıracağım, bu da hello-world.nf. Bir Nextflow işlem hattını çalıştırmak için ihtiyacımız olan tek şey bu. Bu işlem hattı herhangi bir girdi almıyor, bu yüzden başka argümanlara ihtiyacımız yok.

Hadi enter'a basalım ve ne olacağını görelim.

Tamam. Umarım buna benzer bir çıktınız olmuştur. Nextflow'un çalıştığını ve hangi sürümü kullandığını söyleyen birkaç bilgi parçamız var. Hangi betiğin başlatıldığını söylüyor ve bu belirli iş akışı çalıştırması için rastgele oluşturulmuş bir ad veriyor. Bu durumda, benimki "gloomy_crick" olarak adlandırıldı.

Bunun en önemli kısmı, işlem hattında hangi adımların çalıştığını söylemesidir. sayHello adlı process'imizin çalıştığını görebilirsiniz ve bir kez çalıştı ve yüzde yüz tamamlandı.

Buradaki bu kısım, o belirli iş akışı görevi için hash'tir. Her process bir veya daha fazla kez çalışır ve bu çalıştırmaların her birine görev denir.

## 2.2. Çıktıyı ve kayıtları work dizininde bulun

Her görev, çalıştığı kendi izole dizinine sahiptir, bu nedenle iş akışının geri kalanının yürütülmesinden ayrıdır. Bu hash, work dizinindeki dosya yapısına karşılık gelir. "tree work" yaparsam, a0 ve sonra kısa hash'in daha uzun bir versiyonunu ve sonra output.txt dosyamızı görebiliriz. Bunu kenar çubuğunda da görebilirsiniz.

Kenar çubuğunda burada bazı ek dosyalar olduğunu görebilirsiniz. Bunların terminalde görünmemesinin nedeni, bunların gizli dosyalar olmasıdır, nokta ile başlarlar. Ve gerçekten de, "tree -a" (tümü için) ve "work" yaparsam, onları burada görebiliriz.

Bu nokta dosyaları, Nextflow'un oluşturduğu her work dizininde bulunur ve her birinin biraz farklı bir görevi vardır. İlk olarak .command.begin, Nextflow için çalıştırılmadan önce görevi kuran bazı talimatlar içerir. .command.run, Nextflow'un kendisi tarafından yürütülen gerçek talimatlardır. Sonra .command.sh muhtemelen en ilginç olanıdır. Bu, process bloğu betiğimizden çözülen betiktir.

Açarsam, output.txt dosyasına "echo Merhaba Dünya" yaptığımızı görebilirsiniz. Bu, bu durumda process'imizle tamamen aynıdır, ancak Nextflow kodumuzdaki herhangi bir değişkenimiz varsa, her görevin farklı bir .command.sh'si olacak ve bu değişkenlerin nasıl çözüldüğünü görebilirsiniz.

Diğer dosyalar görevin nasıl yürütüldüğüyle ilgilidir. Yani .command.err, .log ve .out standart hata, standart çıktı ve ikisinin birleşimidir. Ve .exitcode, Nextflow'a bu görevin hangi çıkış koduyla yürütüldüğünü, başarılı olup olmadığını söyler.

Son olarak, output.txt dosyamız var ve evet, "Merhaba Dünya" beklediğimiz şey bu ve oluşturulan bu.

Tamam, harika. Bu sizin ilk Nextflow çalıştırmanızdı. Tebrikler. Gerçekten bu kadar basit.

Sırada, her seferinde işlem hattının nasıl çalışacağı konusunda bir değişiklik yapmak istediğimizde kodu düzenlemek zorunda kalmamak için bunu biraz daha rahat nasıl yapacağımıza geçeceğiz.

## 3. İş akışı çalıştırmalarını yönetin

Bu dizin yapısı, tüm görevleri ayrı tutmak ve her şeyi düzenli tutmak için harikadır, ancak tabii ki çıktı dosyalarınızı bulmak çok kullanışlı değildir. İşlem hattınızın sonuçlarını bulmaya çalışarak bir sürü iç içe dizin arasında dolaşmak istemezsiniz.

## 3.1. Çıktıları yayınlayın

İyi haber şu ki yapmanız gerekmiyor. Work dizinleri gerçekten sadece Nextflow'un kendisinin kullanması içindir. Bu yüzden yapacağımız şey, Nextflow için "publishDir" adlı bir fonksiyon kullanmak.

İş akışımıza geri dönüyoruz, process'e gidiyoruz. Burada directive adı verilen yeni bir ifade ekleyebiliriz. Nextflow'un bu işlevselliğin nasıl çalıştığını artıran process'lerin üstünde bulunan bu şeylere böyle der ve kullanacağımız publishDir adındadır.

Burada yazmaya başladığımı görebilirsiniz ve VS Code için Nextflow uzantısı bana directive'i önerdi, bu yüzden sadece enter'a basabilirim.

Tamam, bunu "results" adlı bir dizinle takip edeceğim ve çıktı dosyalarını oraya kopyalamasını söyleyeceğiz. Bu yüzden mode copy diyeceğim. Harika. Kaydet'e basacağım ve iş akışını tekrar çalıştıralım.

nextflow run hello-world.nf

Tamamen aynı şekilde çalışıyor. Bu sefer biraz farklı bir hash'imiz olduğunu unutmayın. Nextflow, iş akışını her çalıştırdığınızda farklı bir hash kullanır. Ve sonuç olarak farklı bir work dizinleri setimiz var. Alanlar, EB adında olanlar, ancak tüm dosyaların aynı olduğunu görebilirsiniz. Ancak, bu sefer yeni olan şey, "results" adlı bir dizinimizin de olması.

Burada "results" içinde çıktı dosyamız var. Nextflow'a bunu söyledik. Sonuç dosyalarını "results" adlı bir dizine kaydet ve oraya kopyala dedik. Ve bu şimdi bulmak çok daha kolay. İş akışını başlattığımız yerin yanında orada ve Nextflow'un gerçek yürütmeyi nerede veya nasıl çalıştırdığından bağımsız olarak, tüm farklı dosyalar orada istediğiniz gibi organize edilebilir.

publishDir'in sembolik bağlantıları işleyebildiğini unutmayın, bu paylaşılan bir dosya sisteminde çalışıyorsanız ve alanda tasarruf etmek istiyorsanız iyidir. Ayrıca, bir process tarafından oluşturulan tüm dosyaları çıktı olarak tanımlamanız gerekmez.

Nextflow yalnızca bu output bloğunda tanımlanan şeyleri kopyalayacaktır. Bu nedenle, bu process'in aşağı akışında gerekmeyen adım tarafından oluşturulan ara dosyalarınız varsa, bunları çıktıda tanımlamazsınız ve publishDir'de görünmezler. Bu nedenle, bu bir işlem hattından çıktı dosyalarınızı temiz tutmanın ve iş yerinin bittiği anda ara dosyaları kolayca silmenin bir yoludur.

Hızlı bir not. Workflow output definitions adında, sonunda publishDir'in yerini alacak yeni Nextflow sözdizimi geliyor. Bu bize iş akışından tüm çıktıları workflow bloğunun içinde işlem hattı düzeyinde tanımlamanın bir yolunu veriyor. Denemek isterseniz bu Nextflow dokümanlarında açıklanmıştır. Ancak şimdilik, publishDir bir süre daha kullanılacak, bu yüzden 2025 eğitiminde hala bu var.

## 3.2. İş akışını -resume ile yeniden başlatın

Tamam. Work dizininin burada şimdi iş akışını her çalıştırdığımızda farklı hash'lerle iki sonuç setine sahip olduğunu belirttim. Bu iyi. Ancak bazen gerekli olmadığında her seferinde adımları yeniden hesaplamak istemeyiz.

Belki iş akışınızı yinelemeli olarak oluşturuyorsunuz ve adımlar ekliyorsunuz ve ilk adımların önbelleğe alınmış sürümleri yeniden kullanmasını istiyorsunuz. Veya belki hesaplama sisteminizde iş akışınızın ortasında bir şeyler ters gitti ve kaldığı yerden devam etmesini istiyorsunuz, ancak zaten tamamladığı adımları atlamasını istiyorsunuz.

Nextflow'un bunun için resume adında yerleşik işlevselliği var. Hadi deneyelim. Öncelikle, ne olduğunu hatırlayabilmemiz için work dizinine bir göz atacağım.

Ve sonra "nextflow run hello-world.nf" yapacağım ve burada tek bir komut ekleyeceğim, "-resume".

Tek çizgi olduğuna dikkat edin, bu gerçekten önemli. Çalıştıracağım ve çıktı temelde tamamen aynı görünecek, birkaç küçük farkla.

Burada gri renkte "cached" yazıyor. Bu, Nextflow'un bu sefer görevi çalıştırmadığı anlamına gelir. Gereksinimlerimizle eşleşen bir şey buldu ve adımı yeniden çalıştırmak yerine bu çıktıları doğrudan yeniden kullandı.

Ve evet, buradaki hash'e bakarsanız, bunun önceki bir çalıştırmadan sahip olduğumuz mevcut hash'e karşılık geldiğini görebilirsiniz.

## 3.3. Eski work dizinlerini silin

Tamam. Ancak yinelemeli olarak geliştirme yapıyorsanız, bu iş akışı dosyalarının çoğunu oluşturacaksınız. Alanda sıkıntı çekiyorsanız bu bir sorun olabilir.

Nextflow, birkaç yardımcı komutla bu work dizinlerini temizlememize yardımcı olabilir. "nextflow log" yaparsam. Bu bana bu dizinde yaptığım tüm farklı iş akışı çalıştırmalarının bir listesini verir ve burada çalıştırma adları var. Çalıştırdığımız ilk olan gloomy quick'i görebilirsiniz ve sonra bu iki yeni.

Şimdi o adı alıp "nextflow clean" komutuyla kullanabiliriz. Tek bir çalıştırma adı belirtebilirim. Ya da daha da iyisi, Nextflow'a "-before" ile tek bir iş akışı adından önceki her şeyi silmesini söyleyebilirim ve "stupefied_shaw" yazacağım. Bu benim en son çalıştırmam, "-n".

"-n" komutu Nextflow'a bunu gerçekten hiçbir şeyi silmeden kuru bir çalıştırma olarak yapmasını söyledi ve hangi hash dizinlerinin kaldırılmış olacağını söylüyor. Evet, bu sadece ilk çalıştırmadan olan. İkinci çalıştırmaların her ikisi de aynı hash dizinini kullanıyor.

Yine çalıştıracağım, ama şimdi kuru çalıştırma için "-n" yerine, "-f" yapacağım, force için ve o hash dizinini kaldırdı. Şimdi "tree work" yaparsam, sadece bu çıktı dosyasının kaldığını görebiliriz.

Harika. Böylece orada bir sürü disk alanı temizlemeyi başardık.

Work dizinlerini silerken dikkat edilmesi gereken birkaç şey, sonuçlar dizininize sembolik bağlantı şeyleri yaparsanız, o sembolik bağlantı kaynakları şimdi silinecek ve sonuçlarınız sonsuza kadar gidecek. Bu yüzden copy modunu kullanmak daha güvenli bir şeydir ve genel olarak önerdiğimiz şeydir.

İkincisi, Nextflow'un resume işlevselliği bu work dizinlerine dayanır. Bu nedenle bunları silerseniz ve Nextflow'u tekrar çalıştırırsanız, resume işlevselliği artık çalışmayacaktır. Bu nedenle, hangi şeylere ihtiyacınız olabileceğini veya olmayabileceğini takip etmek size kalmış ve yalnızca bunun güvenli olduğundan emin olduğunuzda şeyleri silin.

Yapabileceğimiz diğer şey, iş akışı çalıştırmamızı bitirdikten ve artık ihtiyacımız olmadığından emin olduğumuzda tüm work dizinini silmektir.

Yani "rm -r work" yapabilirim. Orada önemli bir şey olmadığını biliyorum. Results dizininde kopyaladığımız yerde önem verdiğim sonuçlarım var. Ve bu yüzden work dizinini silmek güvenliydi. Bu yaklaşımlardan hangisini kullanacağınız size kalmış.

## 4. Komut satırında iletilen değişken girdiyi kullanın

Tamam, sırada ne var? İş akışı betiğimizde burada output.txt dosyasında bazı değerleri sabit kodladığımızı ve bunun için daha iyi bir yol olabileceğini belirtmiştim.

Hadi buna başlayalım. Yapacağımız üç şey var. Process'e yeni bir girdi ekleyeceğiz. Process betiğine bu girdiyi nasıl kullanacağını söyleyeceğiz ve sonra Nextflow'u çalıştırırken komut satırında bir flag ile dinamik olarak kullanabilmemiz için iş akışında bunu bağlayacağız.

İlk önce öncelikler. Hadi buraya bir input bloğu ekleyelim. Çıktı ile aynı. Bu, process için yeni bir bölüm ve "val greeting" diyeceğim.

Burada "val" dediğime dikkat edin, bu bunun bir path değil, bir değişken olduğunu söylüyor.

Sonra betiğe inebilirim ve sonra bu sabit kodlanmış metni buradan çıkarıp $greeting yapabilirim. Bu diğer programlama dilleri gibi çalışır. Burada bir değişken tanımlıyoruz ve bu script bloğu içinde ona başvuruyoruz. Nextflow bu process'i çalıştırdığında, değişken enterpolasyon yapılacak. Ve o .command.sh dosyasına gittiğimizde, bunun yerine gerçek sabit kodlanmış dizgiyi göreceğiz.

## 4.1.3. Bir CLI parametresi ayarlayın ve bunu process çağrısına girdi olarak sağlayın

Tamam, ama değişkeni nerede sağlıyoruz? Sonra workflow bölümüne gidiyoruz ve buradaki uzantının artık bir girdi beklediğini söylediğini ve bana bir uyarı verdiğini görebilirsiniz.

Şimdi, yapabileceğimiz en basit şey sadece sabit kodlamak. "Merhaba Dünya" yazabilir ve o dize girdisini process'e sağlayabilirim. Ama yine, bu gerçekten hiçbir sorunu çözmezdi. Her seferinde bir şeyi değiştirmek istediğimizde yine geri dönüp işlem hattı kodunu düzenlemek zorunda kalırdık, bu da iyi değil.

İyi haber şu ki Nextflow'un parametreler adı verilen komut satırı argümanlarını işlemek için yerleşik bir sistemi var. Bu yüzden bunun yerine, params adında bu özel değişkenlerden birini kullanabilirim ve istediğim gibi adlandırabilirim, ama iş akışı mantığıyla eşleşmesi için greeting diyeceğim.

Kaydet'e bas ve bununla ne yapabileceğimizi görelim.

Terminale geri dönersem. Yani "nextflow run hello-world.nf" yapıyoruz. Daha önce olduğu gibi, ancak temel fark --greeting yapmamız

Burada iki çizgi olduğuna dikkat edin çünkü bu bir parametredir. Daha önce iş akışını devam ettirdiğimizde, bu tek bir çizgiydi. Bunun nedeni resume'un temel bir Nextflow seçeneği olması ve bunun işlem hattımıza özgü bir parametre olmasıdır.

İkisini karıştırmayın. Bunu yapmak kolay. --resume yerine --resume yapsaydınız, o zaman bu "params.resume" olurdu, bu da hiçbir şey yapmazdı. Aynı şekilde, burada tek bir çizgi yapsaydınız, Nextflow bunu anahtar bir argüman olarak tanımazdı.

Yani params.greeting'e karşılık gelen --greeting.

Şimdi bunu istediğim herhangi bir metinle takip edebilirim. Yani şu anda İsveç'teyim, bu yüzden "Hej världen" diyeceğim.

Hadi çalıştıralım, ne olacağını görelim, hakikatin anı.

Tamam, daha önce olduğu gibi, process'in tekrar çalıştığını görebilirsiniz, tek bir yürütmeyle sayHello.

Bu, publishDir "results" dizinindeki dosyanın üzerine yazmış olacak. Bu yüzden dosyaları yeniden çalıştırırken dikkatli olun çünkü yayınlanan havadaki şeyler üzerine yazılacak.

Şimdi "code results/output.txt" yapabilirim ve evet, çıktımız güncellendi ve şimdi "Hej världen" yazıyor.

## 4.2. Komut satırı parametreleri için varsayılan değerleri kullanın

Tamam, bu harika. Ancak şimdi sorun, iş akışımızın her zaman bu parametreyi tanımlamamıza bağımlı olması ve iş akışınız için şeylerin mantıklı bir şekilde çalışması için varsayılanları geçersiz kılmadığınız sürece mantıklı varsayılanlara sahip olmanın güzel olmasıdır.

Bunu yapma şeklimiz, iş akışı betiğimizde parametre için varsayılan bir değer ayarlamaktır.

Bu yüzden hello-world.nf dosyama geri dönersem, workflow'un hemen üstündeki betiğe girebilirim, "params.greeting" yazabilir ve bunu diğer değişkenler gibi tanımlayabilirim. Hadi buraya bir dize koyalım ve "Holà mundo!" diyelim.

Şimdi bu parametrenin tanımlanmış varsayılanı var, burada kullanılacak, veya hala --greeting ile komut satırında geçersiz kılabiliriz, daha önce yaptığımız gibi.

Hadi çalışıp çalışmadığını kontrol edelim. "nextflow run hello-world.nf"

Bu sefer komut satırı argümanı yok ve doğru şeyi yapıp yapmadığını kontrol edin.

"code results/output.txt". Ve işte orada. Varsayılanımızı aldık.

Tamam, tekrar deneyelim, size yalan söylemediğimi kontrol edelim. Tekrar çalıştıralım, ama --greeting yapalım ve eğitim materyalinden örneği kullanalım, "Konnichiwa!" diyelim.

İş akışını yeniden çalıştırır ve evet, üstteki çıktı dosyamız komut satırında sağladığımız yeni değerle güncellendi.

Harika. Bu, herhangi bir Nextflow iş akışı yazmak için gerçek merkezi bir yönüdür. İşlem hattı kodunuzda mantıklı varsayılanlar tanımlamak, ancak terminalde komut satırı argümanlarına sahip olarak son kullanıcı için yapılandırmayı çok kolay hale getirmek.

Son kullanıcının yapılandırmayı birden fazla farklı yerde geçersiz kılabileceğini unutmayın. Yaptığınız her Nextflow çalıştırmasına uygulanan ana dizininizde bir yapılandırma dosyanız olabilir. Başlatma dizininde bir yapılandırma dosyanız olabilir. İşlem hattı dizininde bir yapılandırma dosyanız olabilir. Tüm bu farklı yapılandırma konumları, Nextflow dokümanlarında açıklanan belirli bir sırayla yüklenir.

Tamam, birinci bölümün sonu bu. Nextflow'da bir process ve workflow ile ilk iş akışı betiğimize sahip olduk. Girdilere, çıktılara, betiklere ve yayınlamaya ve parametreleri ve bir girdi kanalını process'imize nasıl bağlayacağımıza baktık.

Tebrikler, Nextflow kodu yazmaya doğru ilk adımınız tamamlandı.

Biraz mola verin ve birkaç dakika sonra ikinci bölümde görüşürüz.

[Sonraki video transkripti :octicons-arrow-right-24:](02_hello_channels.md)
