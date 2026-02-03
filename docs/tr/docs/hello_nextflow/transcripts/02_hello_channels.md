# Bölüm 2: Merhaba Kanallar - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [eğitim materyaline](../02_hello_channels.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca yönlendirici amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba, Merhaba Nextflow'un ikinci bölümüne hoş geldiniz.

Bu bölümün adı Merhaba Kanallar. Nextflow'un bu temel parçası hakkında konuşacağız.

Kanallar, pipeline'ınızdaki farklı adımları birbirine bağlayan şeylerdir, verilerinizin ve mantığınızın iş akışınız boyunca akmasının yoludur.

Pekala, hadi başlayalım.

training.nextflow.io adresine giderek başlayalım

Yan çubukta Merhaba Nextflow'u bulun ve ikinci bölüme tıklayın. Merhaba Kanallar.

Tüm materyal burada yazılı, böylece kendi hızınızda ilerleyebilir ve kaçırmış olabileceğiniz herhangi bir şeyi yakalayabilirsiniz.

Web sitesini açtıktan sonra, Codespaces'i yükleyebilirsiniz ve son bölümün sonunda kaldığımız yerden devam edeceğiz.

## 0. Isınma: hello-channels.nf dosyasını çalıştırın

Bu bölüm için farklı bir dosyayı düzenleyeceğiz. Bu dosyanın adı Merhaba Kanallar, böylece yan çubukta bulabilirsiniz, açmak için çift tıklayın.

Şimdi, eğer birinci bölümden yeni geldiyseniz, bu dosya size çok tanıdık gelecektir. Buradaki başlangıç noktası temelde birinci bölümü bitirdiğimiz yer, sayHello adlı işlememiz, girdimiz, çıktımız, publishDir'imiz ve params.greeting'imiz ve basit iş akışımız ile.

Yeni bir dosya ile başlıyoruz, bu yüzden herkes için eşit bir zemin, ama isterseniz önceki dosyanızla devam edebilirsiniz.

Not edin, temiz bir başlangıç noktası olması için buradaki tüm .nextflow\* dosyalarını ve work dizinlerini de sildim. Bunu yapıp yapmadığınız önemli değil, size kalmış.

Pekala. Bu pipeline'ın beklediğimiz gibi çalıştığını kontrol ederek başlayalım. Burada terminali açacağım.

"nextflow run hello-channels.nf" yazıp enter'a basıyorum.

O küçük iş akışını çalıştıracak, sayHello adımımızı çalıştıracak, o hash ile bir work dizini oluşturacak ve işte results klasörümüz ve varsayılan params.greeting'imizden çıktı dosyamız.

Bu harika. Birinci bölümle tamamen aynı, beklediğimiz gibi çalışıyor.

## 1. Değişken girdileri bir kanal aracılığıyla açıkça sağlayın

Birinci bölümde, aslında zaten kanalları kullanıyordunuz, sadece farkında değildiniz. Burada bir string belirttiğimizde, Nextflow otomatik olarak bizim için o string'in etrafında bir kanal oluşturdu, sadece bir process çağırdığımızı bildiği için, bu yüzden bir girdi kanalına ihtiyacımız vardı.

Yapacağımız ilk şey, kanalın kendisini yazarak bunu açık hale getirmek.

## 1.1. Bir girdi kanalı oluşturun

Bu yüzden betiğin altındaki workflow'a gideceğim ve greeting_ch diyeceğim. Bu, bir değişken adının sonuna \_ch koymak için Nextflow kodunda sık kullandığımız bir gelenektir, böylece bir kanal olduğunu belirlemek kolaydır, ama bunu yapmak zorunda değilsiniz. Eşittir channel of Merhaba Kanallar.

Az önce kullandığımız şey, Nextflow dilinde "Channel Factory" denen bir şeydir. İşte bu şey, bu değişkeni yeni bir kanala ayarlıyoruz ve buradaki bu channel factory bize belirli bir şekilde bir kanal oluşturuyor.

Nextflow'un farklı türde girdilerden kanallar oluşturmak için bir avuç farklı channel factory'si vardır. Dot of en basit olanıdır ve ona verdiğimiz herhangi bir string'i alır.

VS Code'da bu kelimelerin üzerine geldiğimde, Nextflow uzantısının bana bu sözdiziminin ne yaptığını açıklayan bir açılır pencere verdiğini fark edin ve o açılır pencerenin altında da daha fazla oku metni var.

Buna tıklarsam, Nextflow dokümanlarını açacak. Yeni bir sekmede ve beni doğrudan bu spesifik şeyin dokümantasyonuna götürecek. Bu durumda channel.of için.

## 1.2. Kanalı process çağrısına girdi olarak ekleyin

Uzantının bize ayrıca bir uyarı verdiğini, burada yeni bir kanal oluşturduğumuzu ama herhangi bir şey tarafından kullanılmadığını söylediğini unutmayın.

Öyleyse, hadi bunu düzeltelim. Yeni kanal adını alacağım ve bu params.greeting'i yeni kanalımızla değiştireceğim.

Artık --greeting komut satırı bayrağını kullanmıyoruz, params.greeting kullanılmıyor, bu string'i tekrar sabit kodlamaya geri dönüyoruz. Bu tamam. Sadece işleri basit tutmaya çalışıyorum. Daha sonra geri gelip params'ı tekrar kullanacağız.

## 1.3. Workflow komutunu tekrar çalıştırın

Pekala, bunun çalıştığını bir kontrol edelim. Terminali açıyorum ve tekrar not ediyorum. Nextflow run hello channels. output.txt'yi kontrol edin, ve işte orada.

Harika, biraz sıkıcı bir örnek, daha önce yaptığımızla tamamen aynı şeyi yapıyor, ama şimdi en azından mantık biraz daha net. Yeni bir kanal yazmak konusunda açık oluyoruz.

Aynı şeyi yapmak için etkili bir şekilde daha fazla kod yazdık. Ancak kanallarımızı nasıl oluşturduğumuz konusunda biraz daha karmaşık hale geldikçe bu daha mantıklı olmaya başlayacak.

## 2. İş akışını birden fazla girdi değerinde çalışacak şekilde değiştirin

Pekala, bunu biraz daha ilginç hale getirelim. Bir Nextflow pipeline'ını tek bir girdi üzerinde çalıştırmak istemeniz çok nadirdir, o yüzden ona birkaç girdi verelim.

## 2.1. Girdi kanalına birden fazla karşılama yükleyin

Buradan dokümanlardan. Bu farklı string'leri, üç tanesini kopyalayacağım. Hello, Bonjour, Olà. Oh, umut. Copilot birkaç tane daha öneriyor. Öyleyse bunları tab enter ile girelim.

Buradaki Nextflow dokümanları bize bu operatöre birden fazla değer verebileceğimizi söylüyor, bu yüzden çalışması gerekir, ama deneyelim ve ne olduğunu görelim.

## 2.1.2. Komutu çalıştırın ve günlük çıktısına bakın

Evet ve hayır. Bakalım. Burada beş görevden beşinin çalıştığını söylüyor, ama bize sadece bir hash gösteriyor, bu biraz garip. Sorun değil. Her şey burada beklendiği gibi. Varsayılan olarak. Nextflow, terminale ANSI kontrol kodları adı verilen özel bir çıktı türü kullanır, bu da çalıştırılan tüm farklı işlemlerin güzel bir sıkıştırılmış görünümünü vermek için belirli satırların üzerine yazar.

Bu, daha büyük iş akışlarınız olduğunda ve yüzlerce veya binlerce farklı örnek çalıştırdığınızda çok daha mantıklıdır. Terminalde o kadar çok çıktı oluşturabilirsiniz ki, bakmak imkansızdır, oysa bu güncellenen görünüm size gerçek zamanlı bir ilerleme sağlar.

## 2.1.3. Komutu tekrar -ansi-log false seçeneğiyle çalıştırın

İsterseniz, tekrar çalıştırabilirsiniz ve bu sefer "-ansi-log false" diyen tek bir tire ile ek bir Nextflow çekirdek argümanı kullanacağım. Bu, Nextflow günlük çıktısının önceki sürümünü kullanır. Ve burada başlatılan tüm bireysel işlemleri görebilirsiniz.

Bunu yapıp yapmamak size kalmış. Nextflow'dan gelen çıktı her iki durumda da tamamen aynıdır.

## 2.2. Çıktı dosya adlarının benzersiz olacağından emin olun

Pekala, çıktı dosyalarına bakalım o zaman, results'a gideceğiz. Ama sadece tek bir çıktı dosyası var. Ne oldu? İşlemin birçok kez çalıştığını gördük. Work dizinine gidebilir ve tüm farklı hash'leri görebiliriz, tüm görevler düzgün bir şekilde yürütüldü. Ama burada işlememizde hatırlarsanız, her şeyi bir output.txt dosyasına kaydediyor ve sonra bunu bu dizine yayınlıyoruz.

Yani aynı dosya beş kez oluşturuldu ve sonra beş kez üzerine yazıldı. Ve sadece hangi görevin en son yürütülmüş olduğuna sahibiz.

## 2.2.1. Dinamik bir çıktı dosya adı oluşturun

Bunu düzeltme yolumuz dinamik bir çıktı dosya adı kullanmaktır. Burada işlem içinde zaten greeting adında bir değişkenimiz var, bu yüzden bunu çıktı dosya adında kullanabiliriz. Bunu kopyalıyorum ve $greeting-output.txt yapıyorum.

Bunu tırnak içine alacağım, böylece bash burada sızabilecek herhangi bir boşluktan kafası karışmasın. Ve sonra aynı dosya adını alacağım ve buradaki çıktıyı güncelleyeceğim.

Çıktının bununla eşleşmesi gerçekten önemlidir, çünkü aksi takdirde bu dosya bulunamayacak ve Nextflow çökecektir.

Gerçekten önemli bir düzenleme daha yapacağım, bu tek tırnakları çift tırnaklara çevireceğim. Bunu yaptığımda kodun renginin değiştiğini fark edin. Bu değişken yalnızca çift tırnak kullanırsak genişletilir. Burada tek tırnak kullanırsam, gerçek bir değer olarak kullanılır ve $greeting-output adında tek bir dosya elde ederim, bu da istediğim şey değil.

## 2.2.2. İş akışını çalıştırın

Öyleyse çift tırnakları geri koyalım ve bir deneyelim.

Başlamadan önce dizinimı temizleyeceğim, böylece yeni dosyaları görmek kolay olur. .nextflow, work ve results adlı her şeyi sileceğim.

Ve o Nextflow komutunu tekrar çalıştıracağım ve hangi dosyaların oluşturulduğuna bakalım. Yani orada beş işlemi çalıştırıyor. Çok yakından izliyorsanız, çalışırken o satırın güncellendiğini görmüş olabilirsiniz.

Ve şimdi results dizinine gidebiliriz ve elbette beş farklı çıktımız var ve hepsi farklı karşılama ile öneklenmiş.

Bunların her birini açarsam, her birinin karşılık gelen karşılamayı içerdiğini göreceğiz. Harika. İstediğimiz bu.

## 3. Bir kanalın içeriğini dönüştürmek için bir operatör kullanın

Pekala, şimdi kanalların ne olduğunu ve channel factory'lerin ne olduğunu biliyoruz. Peki ya operatörler? Bu, Nextflow dilinin bir başka terimi olup, bize belirli şeyler yapmak için kanallar üzerinde işlem yapmamıza izin veren bir dizi fonksiyondur. Nextflow, kanalları çeşitli şekillerde manipüle etmemize izin veren bir operatör paketi ile birlikte gelir.

## 3.1. Kanala girdi olarak bir değer dizisi sağlayın

Bunu bir örnekle inceleyelim. Diyelim ki bu girdi string'lerini almak istiyoruz, ama bunları doğrudan bir channel factory'ye koymak yerine, onları bir dizi olarak tanımlamak istiyoruz.

## 3.1.1. Girdi değişkenini ayarlayın

Öyleyse bunları alacağım ve yukarıda yeni bir satır olarak yapacağım ve greetings, dizi diyeceğim.

İşte oldu. O dizi değişkenini alacağım ve onu channel.of'a koyacağım ve kaydedeceğim.

## 3.1.3. İş akışını çalıştırın

Şimdi, bakalım ne olacak. Terminalime geri dönüyorum. Sadece o geçici dosyaların hepsini tekrar temizleyeceğim. Ve iş akışını çalıştıralım.

İyi değil. Pekala. Bozuldu. Sorun değil. Bu sefer bozulmasını bekliyordum. Bir Nextflow iş akışı başarısız olduğunda neyin yanlış gittiğini hata ayıklamak, bir Nextflow geliştiricisi olmanın önemli bir parçasıdır. Bu çok olacak ve hata mesajının ne söylediğini ve bununla nasıl başa çıkılacağını anlamak önemlidir.

Nextflow hata mesajları aslında oldukça yapılandırılmıştır. Bize hangi işlemin yanlış gittiğini söyler. Bir sebep için bize bir hata mesajı verir. Çalıştırmaya çalıştığı komutun o belirli görev içinde ne olduğunu, çıkış durumunun ne olduğunu, çıktının o görev work dizininin nerede olduğunu söyler.

VS Code'da buna option, tıklayabiliceğimi ve yan çubukta açtığını unutmayın, böylece doğrudan oraya gidebilir ve .command.sh dosyası dahil olmak üzere önceki bölümde bahsettiğimiz tüm bu gizli dosyaları görüntüleyebilirim. Bunun burada yürütülen komutlarla aynı olduğunu görebilirsiniz.

Bu dosyaya bakarak, burada neyin yanlış gitmiş olabileceğine dair bir fikir edinebiliriz, dizideki her öğe için tek bir görev çalıştırmak yerine, sadece tüm diziyi bir string olarak bir kerede sağladı. Bu yüzden o diziyi kanala geçirmeden önce bireysel değerlere açmamız gerekiyor. Geri dönelim ve bir operatör kullanarak bunu yapıp yapamayacağımızı görelim.

## 3.2. Kanal içeriğini dönüştürmek için bir operatör kullanın

Bu durumda, diziyi kanala geçirmeden önce değiştirmeyeceğiz. Kanalı beklediğimiz şekilde davranacak şekilde ayarlayacağız. Bunu flatten operatörünü kullanarak yapacağız, dot yazmaya başlayabilirim ve VS Code uzantısının elimizde olan tüm farklı operatörleri önermeye başladığını görebiliriz.

## 3.2.1. flatten() operatörünü ekleyin

Ve flatten'ı seçeceğim. Bu bağlamda Nextflow için beyaz boşluğun önemli olmadığını unutmayın. Yani isterseniz bu operatörleri yeni bir satıra koyabilirsiniz. Böylece bunu buraya bırakabilirim ve ".of" altında oturacak şekilde girintileyebilirim ve insanların genellikle bu şekilde bir kanala çok sayıda operatör zincirlediğini ve okumayı kolaylaştırmak için bu şekilde girintilediğini göreceksiniz.

Ayrıca görebilirsiniz, daha önce olduğu gibi bunun üzerine gelebilir ve flatten operatörünün ne yaptığını okuyabilir ve ayrıca istersem dokümantasyona bir bağlantı izleyebilirim.

Bu operatör, içinde tek bir dizi olan bu kanalı alıyor ve dizi değerlerini ayırıyor.

## 3.2.2. Kanal içeriğini incelemek için view() ekleyin

view operatörünü kullanarak kanallara göz atabiliriz ve burada birkaç tane ekleyeceğim. Bu, diğer dillerde print ifadeleri kullanmak gibidir. Bu yüzden dot view yapacağım ve sonra bu kıvrımlı parantezleri kullanacağım.

Buna closure denir. Bu temelde view operatörüne, kanal içindeki her öğe üzerinde yürüteceği ek kod verir. Bu durumda, flatten'dan önce karşılama diyeceğim. Karşılama.

Burada sadece bu closure'ın kapsamı içinde olan bir değişken tanımlıyorum. Yani bu değişken sadece burada kullanılıyor ve istediğimi çağırabilirim. Önemli değil. Okumayı kolaylaştırmak için sadece greeting kullanıyorum.

Bazı Nextflow pipeline'larında, insanların "$it" adlı özel bir örtük değişken kullandığını görebilirsiniz. Bunun gibi. Bu, Nextflow kodu içinde özel bir değişkendir, bir değişken tanımı yapmak zorunda kalmamanız için bir kısayoldur. Ancak zamanla, bunun Nextflow'a yeni başlayanlar için çok net olmadığını düşünüyoruz ve artık "$it" kullanımını caydırıyoruz.

Bu yüzden greeting ile önceki davranışa sadık kalacağım ve bunu bu şekilde kullanacağım çünkü bu daha açık ve neler olup bittiği konusunda daha net.

Sonra bu satırı kopyalayacağım ve flatten argümanlarından sonra tamamen aynı şeyi tekrar yapacağım. view operatörü biraz özeldir çünkü öğeler üzerinde bir şey yapar, ancak aynı zamanda bunları bir sonraki operatöre geçirmeye devam eder, böylece bunu bu şekilde bir operasyon zincirinin ortasına zincirleyebiliriz ve oradaki durumu yazdıracak ve devam edecektir. Umarım bu bize flatten operatöründen önce ve sonra kanalın nasıl göründüğünü gösterecektir.

## 3.2.3. İş akışını çalıştırın

Hadi deneyelim. Temizle. Çalışma alanındaki her şeyi temizle. Pipeline'ı tekrar çalıştır.

Tamam, yani beş işlemimizi tekrar çalıştırdığını görebiliriz. Bir hata ile çökmedi, bu kesinlikle iyi. Ve şimdi flatten'dan öncesi var ve elbette dizimiz var ve flatten'dan sonra, dizinin her öğesi için bir kez olmak üzere beş kez yazdırılmış. Tam olarak umduğumuz şey bu. Yani bu gerçekten iyi bir haber. Ve bu tam olarak koddan beklediğimizle uyuyor.

Artık bu hata ayıklama ifadelerine ihtiyacımız yok, bu yüzden ya bunları yorumlayabilirim ya da silebilirim. Kodumu güzel ve temiz tutmak için onları sileceğim. Pekala, harika. Bu örnek artık güzel bir şekilde çalışıyor ve kanalların biraz daha karmaşık mantık yapabildiğini görmeye başlayabiliriz.

## 4. Bir CSV dosyasından girdi değerlerini ayrıştırmak için bir operatör kullanın

Şimdi bunu bir dizi girdi içeren bir dosya kullanarak yapmayı deneyeceğiz. Bu, bir örnek sayfası veya bir metadata CSV'si kullanarak Nextflow pipeline'ları yazmanın çok yaygın bir yoludur.

## 4.1. Betiği, karşılamaların kaynağı olarak bir CSV dosyası bekleyecek şekilde değiştirin

Yan çubuğa gidersem, örnek deposunda greetings.csv'yi görebilirsiniz ve bu, üç farklı karşılama ile sadece üç satır içeren çok, çok basit bir CSV dosyasıdır. Bakalım bu CSV dosyasını iş akışımız içinde kullanabilecek miyiz.

Şimdi birinci bölümde yaptığımız gibi params kullanmaya geri döneceğim, böylece bir komut satırı girdisi alabiliriz.

Bu greetings dizisini sileceğim.

## 4.1.1. Girdi parametresini CSV dosyasına işaret edecek şekilde değiştirin

params greeting'i dosya adına ayarlayacağım, bu da greetings.csv ve bu özel değişkeni kanalı oluşturmak için kullanacağım. Bunu oraya koyacağım ve hatalar kaybolacak. Bunun bu değişkeni varsayılan olarak şimdi ayarladığını unutmayın. Bu yüzden pipeline'ı herhangi bir argüman olmadan çalıştırırsam, greetings.csv kullanacak, ama istersem bu değişkeni üzerine yazmak için --greeting yapabilirim.

## 4.1.2. Bir dosyayı işlemek için tasarlanmış bir channel factory'ye geçin

Pekala, şimdi bir string veya bir string dizisi yerine bir dosya geçiriyoruz, bu yüzden muhtemelen farklı bir channel factory'ye ihtiyacımız var.

Şimdiye kadar kullandığımız "of"'tan kurtulacağız ve bunun yerine .fromPath kullanacağız. Bu tam olarak kulağa geldiği gibi bir şey yapar. Bir string dosya adı veya glob kullanarak değerler yerine yollarla bir kanal oluşturur. Ayrıca flatten operatörünü kaldıracağım çünkü artık bir dosya geçirdiğimize göre buna ihtiyacımız yok.

## 4.1.3. İş akışını çalıştırın

Kaydet'e basacağım, terminali açacağım, iş akışını çalıştıracağım ve sonra ne olduğunu göreceğim.

Pekala. Yine çöktü. Endişelenmeyin. Bunu da bekliyordum. Hata mesajına bakalım ve neyin yanlış gittiğini anlayıp anlayamayacağımızı görelim. Burada yürütülen komutu görebiliriz ve daha önce tüm dizinin yazdırıldığı gibi, şimdi dosyanın içeriğinden geçmek yerine dosya yolunun komuta yankılandığını görüyoruz.

## 4.2. Dosyayı ayrıştırmak için splitCsv() operatörünü kullanın

Yani bunun yerine dosyanın içeriğini kullanmak için başka bir operatöre ihtiyacımız var. Bunun için kullanacağımız operatör splitCsv olarak adlandırılır. Mantıklı, çünkü yüklediğimiz bir CSV dosyası.

## 4.2.1. Kanala splitCsv() uygulayın

Tamam, yani splitCsv. Parantezi kapat. Burada herhangi bir argümana ihtiyacımız yok. Ve yine, burada neler olup bittiğine dair bir içgörü vermek için bazı view operatörleri kullanacağım.

.view csv splitCsv'den sonra. splitCsv'den önce.

## 4.2.2. İş akışını tekrar çalıştırın

Pekala, bunu çalıştırmayı deneyelim ve ne olduğunu görelim.

Tamam, bu sefer biraz daha fazla çıktımız var, ama yine başarısız oldu. view ifadelerine bakabiliriz ve burada split CSV'den önce görüyorsunuz ve önceki hata mesajında gördüğümüz gibi bir dosya yolumuz var. split CSV'den sonra, şimdi CSV dosyasındaki üç satıra karşılık gelen üç değerimiz var.

Ancak, bu değerlerin her birinin köşeli parantezlerle çevrelendiğini görebilirsiniz. Yani bunların her biri kendi başına bir diziydi ve bu bize daha önce sahip olduğumuz aynı alanı verdi, burada sadece tek bir string yerine bir diziyi yankılamaya çalışıyor.

Bir CSV dosyası hakkında düşünürsek, bu bir anlamda mantıklıdır. Tipik olarak, bir CSV dosyasının satırları ve sütunları olacaktır, bu yüzden split CSV iki boyutlu dizi yapar. Dizinin ilk boyutu her satırdır ve sonra her satır için her sütun olan ikinci bir boyut vardır.

Yani burada her satırda sadece tek bir değerimiz var, bu yüzden tek bir sütunumuz var, bu yüzden dosyanın her satırı için tek öğeli bir dizimiz var.

Bu iyi. Ayrıştırılmış CSV dosyasının her satırı için o diziyi daraltmak için sadece başka bir operatöre ihtiyacımız var. Bunu temizleyelim. Terminalden kurtulalım ve ne yapabileceğimize bakalım.

## 4.3. Karşılamaları çıkarmak için map() operatörünü kullanın

Şimdi daha önce kullandığımız flatten operatörünü tekrar kullanabiliriz. Bunun bir diziyi bir dizi değere nasıl daraltabileceğini gördük, bu burada çok iyi çalışırdı. Ancak, iş akışları içinde çok yaygın olan map operatörü adlı başka bir operatörü göstermek için fırsatı kullanacağım.

## 4.3.1. Kanala map() uygulayın

Dot map yapacağım ve item item[0] yapacağım.

Eğer çok sayıda başka kod dili yazıyorsanız, map operatörüne zaten aşina olabilirsiniz. Bir dizi veya kanal gibi yinelenebilir bir şey alır ve bunun her değeri üzerinde bir işlem yapar.

Burada, bu closure'ın kapsamı içinde item adında bir değişken tanımlamamız gerektiğini ve sonra o dizinin sadece ilk değerini döndürmek istediğimizi söylüyoruz. Yani item index sıfır.

Bu etkili bir şekilde diziyi düzleştiriyor. Bunun nasıl daha karmaşık olacak şekilde genişletebileceğimizi görebilirsiniz: CSV dosyamızın altı sütunu olsaydı, ancak sadece dördüncü sütunla ilgileniyorsak, burada belirli bir indekse erişebiliriz. Veya değeri aşağı akış işlemeye geçirmeden önce üzerinde başka herhangi bir işlem türü yapabiliriz.

Yani map operatörü son derece esnektir ve uçuş halindeki kanalları değiştirmek için çok güçlüdür. Yürütmemizde ne yaptığını görebilmemiz için başka bir view ifadesi koyalım. O satırı tanımlayabilir ve aşağı taşıyabiliriz. Ve map'ten sonra.

## 4.3.2. İş akışını bir kez daha çalıştırın

Terminali açalım ve iş akışını çalıştırmayı deneyelim.

Tamam, bu sefer hata yok. Bu iyi bir işaret. Şimdi view ifadelerinden gelen tüm bu farklı çıktıları gözden geçirebiliriz. split CSV'den önce, tek bir yolumuz vardı. split CSV'den sonra, tek değerli dizilerimiz vardı ve sonra map'ten sonra, herhangi bir dizi sözdizimi olmadan sadece değerlerimiz var. results dizinine gidelim ve işte tam olarak istediğimiz gibi davranan çıktı dosyalarımız.

Burada küçük bir bonus var. view operatörlerinin yaptıkları çıktının sırasında biraz karışık olduğunu görebilirsiniz. Bunun nedeni Nextflow'un bu farklı görevlerin paralelleştirilmesini yapıyor olmasıdır. Yani CSV'yi böldükten sonra, bu kanalda üç öğe var ve bu üç öğenin işlenmesini otomatik olarak paralel olarak gerçekleştiriyor. Bu, çıktıların sırasının stokastik olduğu ve değişebileceği anlamına gelir. Bu durumda, sadece bazı view operatörlerinin sonraki adım tamamlandıktan sonra geri döndüğü oldu ve bu yüzden bu sırayla geldi.

Aynı iş akışını tekrar çalıştırırsam. O zaman elbette, farklı bir sırayla geldi ve bu sefer split CSV'lerimizi ve map'lerimizi beklediğimiz sırada aldık.

Bu yüzden sadece unutmayın, bir process görevinden gelen çıktıların sırasına güvenemezsiniz çünkü Nextflow bu paralelleştirmeyi sizin için otomatik olarak gerçekleştiriyor. Nextflow bunu sizin için veri akışı mantığıyla yapıyor ve bu Nextflow'un gerçek gücüdür.

Tamam, bu muhtemelen tüm eğitimin en önemli bölümlerinden biridir. Kanalları, channel factory'leri ve operatörleri anladığınızda, Nextflow'un gücünü ve onu bir programlama dili olarak benzersiz yapan şeyi kavramaya başlarsınız. Bu işlevsellik, Nextflow'un tüm iş akışlarınızı sizin için paralelleştirmesine ve çok temiz bir sözdizimi ve bir itme veri akışı modeliyle son derece karmaşık iş akışı mantığı oluşturmasına izin verir. İlk başta biraz garip bir kavram olabilir, ama böyle kod yazmaya alıştıktan sonra, hızla doğal hissedecek ve farkına varmadan harika iş akışları yazıyor olacaksınız.

Bir mola verin, bir fincan çay, etrafta dolaşın ve üçüncü bölüme geçelim, burada bu kavramları daha karmaşık iş akışlarına genişletmeye başlıyoruz. Bir sonraki videoda görüşürüz.

[Sonraki video transkripti :octicons-arrow-right-24:](03_hello_workflow.md)
