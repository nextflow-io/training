# Bölüm 5: Hello Containers - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım eksiksiz talimatlar için [kurs materyaline](../05_hello_containers.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş geldiniz ve arka plan

Merhaba ve Hello Nextflow'a tekrar hoş geldiniz. Bu, Hello Containers adlı beşinci bölüm. Ve bu kursum bölümünde, bir boru hattı için yazılım gereksinimlerini nasıl kapsülleyeceğimiz hakkında konuşacağız; böylece boru hattını çalıştıran kişilerin yazılımı yüklemeyi düşünmelerine gerek kalmayacak.

Benim kadar uzun süredir biyoinformatikte çalışıyorsanız, sık sık kötü eski günler dediğim zamanları hatırlayabilirsiniz; başkasının boru hattını çalıştırmak veya çalışmalarını tekrarlamak istediğinizde, kullandıkları tüm farklı yazılım araçlarını, aynı sürümlerde, makinenizde derlemek için saatler veya günler harcardınız ve bu bir kabusту. Gerçekten zordu.

Bir HPC üzerinde çalışıyorsanız, sistem yöneticilerinin sizin için yazılım yüklemeye çalıştığı ortam modüllerini kullanmış olabilirsiniz, ki bu iyi ama yine de kusurlu bir yöntemdi.

Ama şimdi bunu yapmanın daha iyi yolları var. Nextflow'un farklı yazılım konteyner teknolojileri için yerleşik desteği var. Docker en yaygın olanıdır. Bugün kullanacağımız da bu. Codespaces'te iyi çalışır. Yerel bilgisayarınızda iyi çalışır ve bulutta iyi çalışır.

Ama aynı zamanda Singularity veya Apptainer, HPC sistemlerinde çok yaygındır ve etkili bir şekilde tamamen aynı şekilde çalışır. Ya da Podman, Shifter, çok benzer olan bir dizi başka araç var.

Biraz benzer ama tam olarak aynı olmayan ve Nextflow'un desteklediği ekstra bir tanesi Conda. Ve Nextflow sizin için süreç bazında Conda ortamlarını yönetebilir, bu kendi Conda ortamlarınızı yapmaktan çok daha iyidir. Ve yine bir boru hattı ile birlikte gönderilebilir.

Bu bölüme konteyner teknolojileri, Docker ve nasıl çalıştıkları hakkında biraz konuşarak başlayacağız. Ve ilk yarıyı sadece Docker'da manuel olarak yapacağız, böylece perde arkasında neler olduğunu ve bunun nasıl çalıştığını anlayacaksınız. Çünkü Nextflow'un ne yaptığını ve iş akışınızın çalıştırıldığında ne yaptığını anlamak için bu gerçekten önemli.

Öyleyse, hadi Codespaces'imize geçelim. Şimdi her şeyi tekrar temizledim, ama Hello Containers'a atlarsanız, tüm scriptlerimizin ve her şeyin modüller bölümünün sonundaki gibi orada olduğunu görmelisiniz. Yani burada modüller dizininde oluşturduğum farklı modüllerimiz var.

Onlar hala orada. Çalışabilmesi için orada olmaları gerekiyor. Ve iş akışı ve çıktı, dosyalarınızın bu dizinde sonlanması için çıktı yayınlama yolunu Hello Containers olarak değiştirmemiz dışında tamamen aynı.

İsterseniz çalıştığını kontrol etmek için bunu şimdi çalıştırabiliriz veya terminale devam edebiliriz.

## 1. Bir konteyneri 'manuel olarak' kullanma

Konteynerlerimizi yönetmek için Docker kullanacağız ve Codespaces'imde kurulu olup olmadığını "docker -v" yaparak kontrol edebilirim, bu bana kurulu olan sürümü gösterir ve her şeyin düzgün çalıştığını gösterir.

Şimdi konteynerler ve Docker'ın gerçekten önemli olan iki kavramı var. Birine image denir, birine de container. Image, kullanacağınız tüm dosya sisteminin anlık görüntüsü gibidir ve container çalışan ortamdır. Bir image kullanarak bir container oluşturursunuz.

Bir kez o containerin içinde olduğunuzda, tipik olarak tam bir işletim sistemi gibi çalışır. Dış dünyadan kesilmiştir. Her şeyden ayrıdır ve bu iyi bir şeydir. Nextflow ile böylesine iyi bir tekrarlanabilirlik elde etmemizin nedeni budur.

Çünkü bir konteyner içinde çalışan görevler için, yerel sisteminizdeki herhangi bir yapılandırma dosyasından etkilenmezler. Diğer dış etkilerden, kendi küçük korumalı alanlarında çalışırlar. Dosyalar daha sonra çok, çok tekrarlanabilir bir şekilde üretilir çünkü her farklı bilgi işlem ortamında çalışan her kişi için tamamen aynı temel kütüphaneleri, tüm aynı bağımlılıkları, tamamen aynı yazılımı kullanıyorsunuz. Açıkçası bunun harika ve muhteşem olduğunu düşünüyorum ve hala bunun mümkün olduğu gerçeği aklımı başımdan alıyor.

## 1.1. Konteyner imajını çekme

Bazı Docker imajlarını ve Docker'ı deneyeceğiz, sisteminizde çalıştırdığınızda Docker, bilgisayarınızda veya bu durumda kod alanında, geçmişte indirilen ve kullanılan tüm farklı imajları ve üzerine inşa edildikleri farklı katmanları takip eden bir docker kayıt defteri tutar.

Docker ile yerel olarak hangi imajlara sahip olduğumuzu "docker image ls" yaparak görebiliriz. Ve bu durumda burada bir sürü Docker imajı olduğunu görebilirsiniz, bunların hepsi bu Codespaces'i kurmakla ilgili. Hepsi dev containerlar ve şeylerle ilgili. Yani bunlar hakkında çok fazla endişelenmenize gerek yok, ama bu kurs ilerledikçe daha fazla imaj ekleyip indirdikçe, o listeyi kontrol edebilirsiniz ve yerel kayıt defterinin çektiğimiz tüm bu şeyleri takip ettiğini göreceksiniz.

Ama yeni bir tane alacağız "docker pull" yaparak. Ve bu Docker'a web'den yeni bir imaj getirmesini söyler.

Ardından o konteyner için URI'yi koyuyoruz. Şimdi bu, yerel olarak oluşturduğunuz ve ardından internete gönderdiğiniz bir docker imajı olabilir. Başka birinin yaptığı bir imaj olabilir. Docker imajları yapmanın çok, çok, çok farklı yolları var, ama muhtemelen en basit yollardan biri bunu dış kaynak kullanmak ve bunu başka birinin sizin için yapmasını sağlamak.

Ve bu eğitimde kullanacağımız şey, Seqera'dan Seqera Containers adlı bir hizmettir.

Şimdi, Seqera Containers tamamen ücretsizdir ve Wave adını verdiğimiz, Nextflow'a tamamlayıcı bir şekilde konteynerleri yönetmek için oluşturulmuş açık kaynaklı bir yazılım parçası kullanır. Ve Nextflow ile uğraştığımız birçok yaygın kullanım durumunu ele alır.

İhtiyacımız olan yazılımın Conda'da, Bioconda'da veya conda-forge kanallarında veya diğer daha alana özel kanallarda paketlenmiş olması çok yaygındır. Ve Wave ve Seqera Containers bundan imaj oluşturmada gerçekten iyidir.

Yani bu web arayüzüne gidebilirim ve "cowpy" adlı bir paketle oynayacağız. Yani istediğim paketin adını yazıyorum. Arama yapıyor, Python paket dizininde buldu, bu yüzden bunu kullanabilirim. Ya da biraz daha beklersem, bioconda ve conda-forge'u arıyor. Ve görebilirsiniz, burada herhangi bir conda kanalı belirleyebilirim. Yani bir Nvidia kanalı veya başka bir şey bulmak isterseniz, bu da çalışmalı.

Ve ardından benim için bir docker imajı mı yoksa bir singularity imajı mı oluşturmasını istediğimi ve ayrıca hangi CPU mimarisini kullanmak istediğimi belirleyebilirim. Yani amd64 veya arm64.

Ve bioconda sonuçları listelendikten sonra, artık mevcut olan tüm farklı sürümleri de görebiliyorum. Bunu koyacağım. Ve şimdi aramaya devam edip Conda'dan daha fazla paket alabilirdim ve bu konteyneri istediğim gibi oluşturabilirdim, ama sadece bunu istiyorum. Yani Get Container'a tıklayacağım.

Şimdi, başka biri daha önce aynı konteyneri istemiş ve bir kayıt defterinden döndürüldü, bu yüzden hemen alıyoruz. Ama başka hiç kimse bu yazılım paketini veya bu yazılım paketi kombinasyonunu hiç istememişse, Wave ve Seqera Containers bunu bizim için anında oluşturur.

Bu URL'yi kopyalayabiliriz ve ayrıca yapı detaylarını görebiliriz. Ve bu bize hizmetin arka planda ne yaptığını gösterir. Bir conda ortam dosyası oluşturdu. Bir docker dosyası ve sonra işte bu, docker yapı sürecini çalıştırıyor. Ayrıca bir tarama, bir güvenlik taraması çalıştırdı, böylece herhangi bir CVE'yi görebilirsiniz. Ve bunun ne zaman oluşturulduğunu söyler.

Wave ve Seqera Containers bundan çok daha fazlasını yapabilir, ama bu en yaygın olan basit bir kullanım durumu türü. Ve bu imajların en az beş yıl boyunca barındırıldığını söylemeliyim. Yani bu URL'leri boru hatlarınıza yerleştirebilir ve yakın zamanda kaybolmayacaklarını bilebilirsiniz.

Yani cowpy için docker imajım için URL'mi aldım.

Şimdi "docker pull" o URL'yi yapabilirim ve tüm farklı katmanları getirip bu imajı yerel olarak kullanılabilir hale indirecek.

## 1.2. Konteyneri cowpy'yi tek seferlik komut olarak çalıştırmak için kullanma

Tamam, şimdi onu gerçekten kullanmayı deneyelim. Yani şimdi "docker pull" komutu yerine "docker run" komutu kullanacağım ve "--rm" bayrağını kullanacağım, bu sadece Docker'a benden yapmasını istediğim şeyi bitirdiğinde bu konteyneri kapatmasını söyler. Ve sonra konteyner için tanımlayıcıyı koyuyorum, bu sadece bir URI.

Ve sonra sonunda, Docker'ın konteyner içinde çalıştırmasını istediğim komutu belirtiyorum. Sadece cowpy diyeceğim, bu Conda Forge'dan yüklenen, imaj içinde mevcut olan aracın adı.

Enter'a basacağım ve işte buyurun. Bir sistemde cowpy'yi çalıştırdık. Bize biraz bilgi veren küçük bir ineğimiz var.

Şimdi dikkat edin cowpy yerel sistemimde kurulu değil. Yani eğer sadece tüm Docker şeyleri olmadan çalıştırırsam, command not found diyor. Yani bu bir imaj çekti. Docker kullanarak bir konteyner oluşturdu ve sonra o konteynerin içine gitti ve bu komutu bizim için çalıştırdı ve çıktıyı terminalimize geri verdi. Çok, çok havalı.

## 1.3. Konteyneri cowpy'yi interaktif olarak çalıştırmak için kullanma

Tamam, şimdi bir adım daha ileri gideceğiz ve bu konteyneri interaktif olarak çalıştıracağız ve etrafa bir göz atacağız, böylece konteyner içinde ne olduğunu görebiliriz.

Yani geri gidip çalıştırma komutumu alırsam ve sondaki cowpy'yi çıkaracağım, çünkü aslında cowpy çalıştırmak istemiyorum. Bir Bash terminali çalıştırmak istiyorum.

Ve sonra buraya geri gideceğim ve "-it" yapacağım, bu Interactive ve Terminal veya TTY anlamına gelir ve enter'a basacağım.

Ve şimdi istem görebilirsiniz, yazdığım kısımdan önceki kısım değişti. Bu, dizini söyleyen Codespaces istemiydi ve şimdi base ve root ve tmp diyor. Yani şimdi konteyner içindeyim ve eğer "ls" yaparsam, bu dizinde gördüğüm dosyaların çalışma alanımdaki dosyalardan farklı olduğunu göreceksiniz.

Ve aslında, yerel codespaces çalışma alanımdan veya yerel sürücümden hiçbir dosyayı Docker konteyneri içinde göremiyorum. Docker konteyner çalışma zamanı tamamen izole edilmiştir ve dışarıdaki bir host dosya sisteminden hiçbir dosyayı yazamaz veya okuyamaz.

Ancak, konteyner içinde yüklü olan yazılımı görebilir ve çalıştırabilirim. Yani cowpy'yi çalıştırabilirim ve cowpy'nin nasıl kullanılacağı hakkında biraz daha görebiliriz. Burada "cowpy 'Hello World'" yapabilirim ve bu, alıntımı aslında küçük bir konuşma balonunun içine koymasını söyler. Ve farklı inek türlerini de çalıştırabilirsiniz, yani inek olmak zorunda değil. "-c" yapabilirsiniz. Ve ben İsveç'teyim, bu yüzden bir geyik seçeceğim. Çok güzel. Ona boynuz verdi.

Ve eğitim dokümanlarında açıklandığını görebileceğiniz oynayabileceğiniz bir sürü farklı tane var.

## 1.3.4. Konteynere veri bağlama

Tamam. Dosya sistemimizdeki dosyalar üzerinde cowpy çalıştırabilseydik güzel olurdu.

Tabii ki, sadece konteynere sahip olmak ve hiçbir şeye erişim olmamak çok kullanışlı değil. Güvenli ve tekrarlanabilir olabilir, ama çok kullanışlı değil.

Peki bunu nasıl yaparız? Bu Docker konteynerinden exit yazarak çıkacağım ve istem bize şimdi tekrar normal Codespaces'imizde olduğumuzu söylüyor.

Ve aynı komutu tekrar çalıştıracağım. Ama bu sefer buraya bazı ek bayraklar ekleyeceğim. Ve önemli olanı "-v"dir, bu bir volume bağlamak anlamına gelir, bu temelde bir disk alanının bir parçası gibi bir şey.

"-v" iki kısım alır: bir string gibi bir şey var ve sonra iki nokta üst üste ve bir string. Ve ilk kısım konteynere bağlanması gereken yerel dosya sistemidir. Ve sonra ikinci kısım bunun konteyner içinde nereye gitmesi gerektiğidir.

Şimdi sadece tüm yerel dosya sistemimi buraya yüklemek istiyorum. Yani "." mevcut çalışma dizinidir. Yani sadece "." yapacağım ve sonra ":", ve sonra bunu konteyner içinde "my_project" adlı yeni bir dizine koyacağız. Bu gerçekten herhangi bir şey olarak adlandırılabilirdi.

Ve sonra tekrar çalıştıracağım.

Bırakıldığım çalışma dizininde, /tmp, dosyalar orada değil. Ama eğer "ls my_project" yaparsam, işte buyurun: Codespaces'te yerel olarak sahip olduğumuz tüm aynı dosyalar şimdi o yoldaki konteyner içinde mevcut.

Bu okuma ve yazma erişimidir, bu yüzden bu dizinde yeni dosyalar oluşturabilirim ve bunlar host dosya sistemimde görünecek. Bu özel dizin, o zaman konteyner dışındaymışım gibi tam olarak davranır, böylece şimdi okuyabilir ve yazabilir ve şeyler yapabilirim.

## 1.3.5. Bağlanan veriyi kullanma

Tamam, bunu yapabileceğimizi kanıtlayalım. "cat /my_project/data/greetings.csv" yapıyorum. Eğer hatırlarsanız bu dosya içerikleri böyle görünür. Şimdi bunu cowpy'ye aktarabilirim ve inek o dosyanın farklı çıktılarını küçük konuşma balonunda yazdıracak, bu biraz eğlenceli.

Yani görebilirsiniz, şimdi konteynerdeki yazılımı host sistemimizdeki dosyalarla etkileşim kurmak için kullanabiliriz.

Tamam, tekrar dışarı çıkalım ve eğitim materyalinin geri kalanıyla devam edelim.

## 2. Nextflow'da konteyner kullanma

Yani konteyner kullanmak gerçekten havalı. Umarım mantıklıdır. Ve bu konteynerlerin değerini ve analiz yazılımı çalıştırmak için neden kullanışlı olduğunu görebilirsiniz.

Ama tüm bu aynı süreci Nextflow içinde nasıl yaparız? Kendimiz bir sürü Docker komutu çalıştırmak istemiyoruz. Sadece Nextflow'un tüm bunları bizim için halletmesini istiyoruz.

Öyleyse bunu çözelim. Boru hattımıza cowpy çalıştırmak için yeni bir süreç ekleyeceğiz. Tamam, yeni sürecimiz için yeni bir modül oluşturalım. Yani modüllere gidin, buna cowPy.nf diyelim ve sonra eğitim materyalinden kodu kopyalayacağım.

Ama görebilirsiniz süreç çok basit. Şimdiye kadar yaptıklarımıza çok benziyor, girdi dosyamız olan bir yol ile bir girdi bloğumuz var ve ayrıca burada bir değer var, böylece bu bir karakter olacak, böylece isterseniz tekrar bir geyik kullanabiliriz.

Ve sonra bir çıktı, bu burada tek bir dosya, bir yol ve sonra bir script. Ve konteynerin içinde interaktif olarak yaptığımız şeyin aynısını yapıyoruz: girdi dosyasını okumak için "cat" yapıyoruz. O içerikleri cowpy'ye aktarıyoruz. O girdiye göre belirli bir karakter seçiyoruz, cowpy girdi dosyası adlı bir çıktı dosyasına yazıyoruz, bu daha sonra çıktıya yansıtılıyor.

Harika. Bunu dahil edelim. Yani include \{ COWPY \} from "./modules/cowpy.nf", buna cowpy mi dedim? Evet.

Ve sonra yeni sürecimizi iş akışının ana bloğunda çağıralım. Yani cowpy'yi çalıştıralım. Ve yeni cowpy sürecimizi alacağız ve collectGreetings.out diyeceğiz.

Ve eğer hatırlarsanız, bu modül için iki çıktı vardı. Biri outfile, diğeri report. VS Code uzantısı bunları bizim için otomatik öneriyor ve .outfile istiyoruz.

Her zaman buradan bu sürece atlayabilirsiniz. Ya üzerine gelirsiniz ve çıktıların ne olduğunu hızlıca gösterir. Ve daha fazla detay görmek isterseniz içine komut tıklayıp modül dosyasını açabilirsiniz.

Yani işte buyurun. Outfile orada ve bu yol. Yani bu şimdi cowpy sürecimiz için girdi dosyası olacak. Harika.

Şimdi eğer hatırlarsanız, bir cowpy sürecinin iki girdisi var. Ayrıca karakter için değer kanalı vardı. Yani buraya "params.character" ekleyebiliriz. İstersenm bunu sabit kodlayabilirdim, ama onu bir CLI seçeneği yapalım, böylece tire, tire character yapabiliriz.

Doğru. Şimdi az önce çağırdığımız girdi parametresini tanımlamam ve ona bir varsayılan vermem gerekiyor. Yani character, String. Ve geyiği seviyorum, bu yüzden varsayılan olarak moose'a ayarlayacağım.

Doğru, çalıştırmayı deneyelim. Yani eğer Nextflow run hello containers yaparsam, ne olduğunu göreceğiz.

Eski çalışma dizinleri etrafta dolaşsaydı tire resume kullanabilirdim. Ve yine, bu ilk süreçler önbelleğe alınmış olurdu ve biraz daha hızlı olurdu, ama temelde aynı olmalı.

Şimdi hemen görebiliriz ki yeni sürecimize geldiğinde bir hata fırlattı, bize burada cowpy sürecini çalıştırırken bir hata olduğunu söylüyor ve 127 çıkış durumu ile çıktı. Çalıştırmaya çalıştığı komut bu. Doğru görünüyor, beklediğimiz gibi görünüyor. O çıktı dosya adını alıyor, bu doğru görünüyor, bir moose karakteriyle çalıştırıyor ve kaydetmeye çalışıyor.

Ama burada komut hatasının cowpy komutunun bulunamadığını söylediğini görebilirsiniz. Ve bu mantıklı çünkü henüz Nextflow'a bir konteyner kullanmasını söylemedik. Sadece ona cowpy komutunu verdik. Ve daha önce söylediğim gibi, cowpy yerel sistemimizde kurulu değil. Yani çalıştırmayı denediğinde, başarısız oldu.

## 2.3.1. cowpy için bir konteyner belirtme

Nextflow'a mevcut bir konteyner olduğunu ve kullanabileceğini söylememiz gerekiyor. Peki bunu nasıl yaparız?

Modülümüze girelim, üste "container" adlı yeni bir bildirim ekleyeceğiz. Ve bunu bir stringe ayarlayacağız.

Şimdi, eğer hatırlarsanız, Seqera Containers'da, o URL'yi kopyalayabilirim ve bunu buradaki tırnak içine bırakıyorum.

Şimdi geri dönün ve tekrar çalıştırmayı deneyin.

Bakalım bu sefer çalışacak mı.

Ne yazık ki, süreç için bir konteyner tanımlamış olmamıza rağmen tamamen aynı şekilde başarısız oluyor. Yani docker imajımızı kullanabilmek için, iş akışını çalıştırdığımızda Nextflow'a Docker kullanımını etkinleştirmesini söylememiz gerekiyor.

Ve bunu yeni bir yapılandırma dosyası oluşturarak yapacağız. Yani touch nextflow.config diyeceğim.

Bu, boru hattını başlatırken çalışma dizinindeyse otomatik olarak yüklenecek özel bir dosya adıdır. Yani bu Nextflow dot config dosyasına girersem, aslında zaten var olduğunu görebilirsiniz, bunu unutmuştum. Ve burada zaten docker.enabled var, ama false'a ayarlı, bu varsayılan.

Yani eğer bunu bunun yerine True'ya eşit olacak şekilde değiştirirsem, docker.enabled. Ve bunların tümü için Nextflow dokümanlarında referans dokümanları var. Ve ayrıca VS Code uzantısı ile üzerine geldiğimde, bununla ilgili dokümanları çeker ve bunun ne anlama geldiğini ve nasıl ayarlanacağını söyler.

Yani şimdi true'ya ayarladık ve eğer Nextflow'u tekrar çalıştırırsam, Nextflow şimdi o docker imajını yerel olarak henüz yoksa bizim için çekmeyi ve ardından o sürecі o konteyner ortamıyla çalıştırmayı bilecek.

Ve başarıyla çalıştığını ve cowpy'nin yanında küçük bir tik işareti olduğunu görebiliriz. Harika. Eğer yukarı çıkıp sonuçlar dizinine bakarsam, dosya henüz orada değil. Ve bunun nedeni, diğerlerinin hepsi gibi bu çıktı dosyasını yayınlamamız gerekiyor.

Yani iş akışı içindeki yayınlama bloğuna gidiyoruz, mycowpy equals cowpy.out diyoruz.

Ve sonra burada çıktı bloğunda, mycowpy, küme parantezler path. Hoop. Hello containers. Mode, copy.

Şimdi tekrar çalıştırırsam, tamamen aynı şekilde çalışmalı. Tire resume kullanabilirdim ve her seferinde unutuyorum. Ve sonra yukarı çıkıyorum ve şimdi cowpy-COLLECTED adlı yeni bir dosya oluşturuldu ve işte benim geyiğim BONJOUR, HELLO, HOLA diyor. Harika.

Şimdi tabii ki şimdi ayrıca "--character" geçirebilirim. Farklı seçenekler neler? Sanırım bir Turkey var? Yani character Turkey kullanabilirim. Tamamen aynı şekilde çalışacak. Tire resume kullanmak için başka bir fırsatı kaçırdım ve şimdi dosyamızı yüklersek ve şimdi bir Turkey'imiz var. Harika.

## 2.3.4. Nextflow'un konteynerize edilmiş görevi nasıl başlattığını inceleme

Tamam. Son küçük şey. Hadi bu komutu tekrar hızlıca çalıştıralım, bu sefer resume, ve Nextflow'un tüm bunların bizim için çalışmasını sağlamak için perde arkasında ne yaptığını görmek için çalışma dizinine hızlıca bir göz atalım.

Bu sefer süper hızlı, hadi bu çalışma dizinine gidelim, cd work/. Şimdi eğer hatırlarsanız burada bir sürü nokta dosyamız var ve bu durumda ilgilendiğimiz, neredeyse hiç bakmamız gerekmediğini söylediğim .command.run adlı olan.

Eğer code dot command run yaparsam, editörde açacak. Ve bu dosyada arama yapabilirim ve aşağı kaydırırsam Docker run'ı görmeliyim. Ve Nextflow'un Docker etkinleştirildiğinde bizim için docker run komutunu yaptığını görebilirsiniz. Burada bir sürü farklı bayrak ve şey var, ama kendimiz çalıştırırken kullandığımız "-v" bayrağını görebilirsiniz. Ve yerel çalışma alanı dizinini konteynere bağladığını görebilirsiniz, böylece konteyner girdi dosyalarımıza erişebilir ve çıktıları kaydedebilir. Ve sonunda, ayrıca .command.sh'yi çalıştırıyor, bu içinde cowpy komutu olan üretilmiş scripttir.

Ve böylece Nextflow'un iş akışı mantığını aldığını görebilirsiniz, bu gerçekten önemsediğimiz şeydir, analizimize özgü olandır, ve sistemimizde Docker'ın çalışmasını sağlamak için tüm akıllı perde arkası işleri yapıyor.

Ve bunu gerçekten taşınabilir bir şekilde yapıyor, böylece boru hattının son kullanıcısı kullandıkları teknolojiyi değiştirebilir: Docker, Singularity, Apptainer, Conda. Bu gerçekten boru hattı mantığı için önemli değil, ama Nextflow tüm temel altyapı ihtiyaçlarını halledecek, böylece her yerde çalışır.

Ve bu gerçekten Nextflow'un süper gücüdür. Tekrarlanabilirlik ve taşınabilirlik. Ve Nextflow ile iş akışınızı gerçekten paylaşabilir ve diğer insanlar sistemlerinde çalıştırabilir ve sadece çalışacaktır.

Bu yapmak gerçekten, gerçekten zor bir şeydir ve şimdi siz de iş akışlarınızla bunu nasıl yapacağınızı biliyorsunuz.

Tamam, bu bölüm için bu kadar. Eğer bir kursun sonuna giderseniz, konteynerler hakkında tekrar bir sınav bulacaksınız. Umarım hepsi mantıklıydı. Analizle çalışmanın gerçekten havalı bir yolu. Ve konteynerlere yeniyseniz, umarım bunun gidilecek yol olduğuna sizi ikna etmişimdir ve asla geriye bakmayacaksınız.

Ama bununla, belki biraz mola verin ve bir kaç dakika içinde Hello Nextflow'un tamamı yapılandırma hakkında olan son altıncı bölümü geçmek için bana katılın.

Çok teşekkür ederim.
