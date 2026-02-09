# Bölüm 5: Hello Containers - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../05_hello_containers.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş geldiniz ve arka plan

Merhaba ve Hello Nextflow'a tekrar hoş geldiniz. Bu, Hello Containers adlı beşinci bölüm. Ve bu kurs bölümünde, bir boru hattının yazılım gereksinimlerini nasıl kapsülleyeceğimizden bahsedeceğiz, böylece boru hattını çalıştıran kişilerin yazılımı yükleme konusunda düşünmelerine gerek kalmayacak.

Biyoinformatikte benim kadar uzun süredir çalışıyorsanız, genellikle eski kötü günler dediğim zamanları hatırlayabilirsiniz; başka birinin boru hattını çalıştırmak veya çalışmalarını tekrarlamak istediğinizde, kullandıkları tüm farklı yazılım araçlarını aynı sürümlerde yüklemeye çalışarak, bunları makinenizde derlemeye çalışarak saatler veya günler harcardınız ve bu bir kabustu. Gerçekten zordu.

Bir HPC üzerinde çalışıyorsanız, sistem yöneticilerinin sizin için yazılım yüklemeye çalıştığı ortam modüllerini kullanmış olabilirsiniz, bu iyiydi ama yine de kusurlu.

Ama şimdi bunu yapmanın daha iyi yolları var. Nextflow, farklı yazılım konteyner teknolojileri için yerleşik desteğe sahip. Docker en yaygın olanı. Bugün kullanacağımız da bu. Codespaces'te iyi çalışıyor. Yerel bilgisayarınızda iyi çalışıyor ve bulutta iyi çalışıyor.

Ayrıca HPC sistemlerinde çok yaygın olan ve aslında tamamen aynı şekilde çalışan Singularity veya Apptainer var. Ya da Podman, Shifter, bunların hepsi çok benzer olan bir sürü başka teknoloji var.

Biraz benzer ama tam olarak aynı olmayan ve Nextflow'un desteklediği bir ekstra teknoloji Conda. Ve Nextflow, süreç bazında sizin için Conda ortamlarını yönetebilir, bu da kendi Conda ortamlarınızı yapmaktan çok daha iyi. Ve yine bir boru hattıyla birlikte gönderilebilir.

Bu bölüme konteyner teknolojileri, Docker ve bunların nasıl çalıştığı hakkında biraz konuşarak başlayacağız. Ve ilk yarısını sadece Docker'da manuel olarak yapacağız, böylece perde arkasında neler olduğunu ve bunun nasıl çalıştığını anlayacaksınız. Çünkü Nextflow'un ne yaptığını ve iş akışınızın yürütülürken ne yaptığını anlamak için bu gerçekten önemli.

Şimdi Codespaces'imize geçelim. Her şeyi tekrar temizledim, ama Hello Containers'a geçerseniz, tüm betiklerimizin ve her şeyin modüller bölümünün sonundaki gibi orada olduğunu görmelisiniz. Yani modüller dizininde oluşturduğum farklı modüllerimiz burada.

Hala oradalar. Çalışabilmesi için orada olmaları gerekiyor. ve iş akışı ve çıktı, çıktı yayınlama yolunu Hello Containers olarak değiştirmemiz dışında hepsi aynı, böylece dosyalarınız o dizinde oluşuyor.

İsterseniz çalıştığını kontrol etmek için bunu şimdi çalıştırabiliriz veya terminale devam edebiliriz.

## 1. Bir konteyneri 'manuel olarak' kullanma

Konteynerlerimizi yönetmek için Docker kullanacağız ve "docker -v" yaparak Codespaces'imde kurulu olup olmadığını kontrol edebilirim, bu bana kurulu sürümü ve her şeyi gösterir ve düzgün çalıştığını gösterir.

Şimdi konteynerler ve Docker'ın gerçekten önemli olan iki kavramı var. Biri image (imaj) olarak adlandırılır ve biri container (konteyner). İmaj, kullanacağınız tüm dosya sisteminin anlık görüntüsüdür ve konteyner çalışan ortamdır. Yani bir imaj kullanarak bir konteyner oluşturursunuz.

Bir konteynerin içine girdiğinizde, genellikle tam bir işletim sistemi gibi çalışır. Dış dünyadan kesilmiştir. Her şeyden ayrıdır ve bu iyi bir şeydir. Nextflow ile bu kadar iyi tekrarlanabilirlik elde etmemizin yolu budur.

Çünkü bir konteyner içinde çalışan görevler için, yerel sisteminizdeki herhangi bir yapılandırma dosyası tarafından kirletilmezler. Herhangi bir dış etki, kendi küçük korumalı alanlarında çalışırlar. Dosyalar daha sonra çok, çok tekrarlanabilir bir şekilde üretilir çünkü her farklı bilgi işlem ortamında çalışan her kişi için aynı temel kütüphaneleri, aynı bağımlılıkları, tamamen aynı yazılımı kullanıyorsunuz. Açıkçası bence bu harika ve muhteşem bir şey ve bunun mümkün olması. Ve hatta bugün bile hala aklımı başımdan alıyor.

## 1.1. Konteyner imajını çekme

Bazı Docker imajlarını ve Docker'ı deneyeceğiz, sisteminizde çalıştırdığınızda, bilgisayarınızda veya bu durumda kod alanında, geçmişte indirilmiş ve kullanılmış tüm farklı imajları ve bunların üzerine inşa edildikleri farklı katmanları takip eden bir docker kayıt defteri vardır.

Docker ile yerel olarak hangi imajlara sahip olduğumuzu "docker image ls" yaparak görebiliriz. Ve bu durumda burada bir sürü Docker imajı olduğunu görebilirsiniz, bunların hepsi bu Codespaces'i kurmakla ilgili. Hepsi dev konteynerler ve şeylerle ilgili. Yani onlar hakkında çok fazla endişelenmenize gerek yok, ama bu kurs ilerledikçe daha fazla imaj eklediğimizde ve indirdiğimizde, bu listeyi kontrol edebilirsiniz ve yerel kayıt defterinin çektiğimiz tüm bu şeyleri takip ettiğini göreceksiniz.

Ama "docker pull" yaparak yeni bir tane alacağız. Ve bu Docker'a web'den yeni bir imaj getirmesini söyler.

Daha sonra o konteyner için URI'yi koyuyoruz. Şimdi bu, yerel olarak oluşturduğunuz ve ardından internete gönderdiğiniz bir docker imajı olabilir. başka birinin yaptığı bir imaj olabilir. Docker imajları yapmanın birçok, birçok, birçok farklı yolu var, ama tartışmasız en basit yollardan biri bunu dış kaynak kullanmak ve başka birinin bunu sizin için yapmasını sağlamak.

Ve bu eğitimde kullanacağımız şey, Seqera'dan Seqera Containers adlı bir hizmet.

Şimdi, Seqera Containers tamamen ücretsiz ve geliştirdiğimiz Wave adlı açık kaynaklı bir yazılım parçası kullanıyor, bu Nextflow'a tamamlayıcı bir şekilde konteynerleri yönetmek için oluşturuldu. Ve Nextflow ile uğraştığımız yaygın kullanım durumlarının çoğunu ele alıyor.

İhtiyacımız olan yazılımın Conda'da, Bioconda'da veya conda-forge kanallarında veya diğer daha alana özgü kanallarda paketlenmiş olması çok yaygındır. Ve Wave ve Seqera Containers bundan imaj oluşturmada gerçekten iyidir.

Bu web arayüzüne gidebilirim ve "cowpy" adlı bir paketle oynayacağız. Yani istediğim paketin adını yazıyorum. Arama yapıyor, Python paket dizininde buldu, böylece onu kullanabilirim. Ya da biraz daha beklersem, bioconda ve conda-forge'u arıyor. Ve görebilirsiniz, burada herhangi bir conda kanalı belirtebilirim. Yani bir Nvidia kanalı veya başka bir şey bulmak istiyorsanız, bu da çalışmalı.

Ve sonra benim için bir docker imajı mı yoksa bir singularity imajı mı oluşturmasını istediğimi ve ayrıca hangi CPU mimarisini kullanmak istediğimi belirtebilirim. Yani amd64 veya arm64.

Ve bioconda sonuçları listelendikten sonra, artık mevcut olan tüm farklı sürümleri de görebiliyorum. Onu ekleyeceğim. Ve şimdi aramaya devam edebilir ve istediğim şekilde bu konteyneri oluşturmak için Conda'dan daha fazla paket alabilirim, ama sadece bunu istiyorum. Yani Get Container'a tıklayacağım.

Şimdi, başka biri daha önce aynı konteyneri istedi ve bir kayıt defterinden döndürüldü, bu yüzden hemen alıyoruz. Ama başka hiç kimse bu yazılım paketini veya bu yazılım paketleri kombinasyonunu istememiş olsaydı, Wave ve Seqera Containers bunu anında bizim için oluştururdu.

Bu URL'yi kopyalayabiliriz ve ayrıca yapı ayrıntılarını görüntüle'yi görebiliriz. Ve bu bize hizmetin arka planda ne yaptığını gösterir. Bir conda ortam dosyası oluşturdu. Bir docker dosyası ve sonra bu, docker yapı sürecini çalıştırıyor. Ayrıca bir tarama, bir güvenlik taraması çalıştırdı, böylece herhangi bir CVE'yi görebilirsiniz. Ve bunun ne zaman oluşturulduğunu söyler.

Wave ve Seqera Containers bundan çok daha fazlasını yapabilir, ama bu en yaygın olan basit bir kullanım durumu. Ve bu imajların en az beş yıl boyunca barındırıldığını söylemeliyim. Yani bu URL'leri boru hatlarınıza yerleştirebilir ve yakın zamanda kaybolmayacaklarını bilebilirsiniz.

Yani cowpy için docker imajım için URL'mi aldım.

Şimdi "docker pull" o URL'yi yapabilirim ve tüm farklı katmanları getirecek ve bu imajı yerel olarak kullanılabilir olması için indirecek.

## 1.2. Konteyneri tek seferlik bir komut olarak cowpy çalıştırmak için kullanma

Tamam, şimdi onu gerçekten kullanmayı deneyelim. Şimdi "docker pull" yerine bir "docker run" komutu kullanacağım ve "--rm" bayrağını kullanacağım, bu sadece Docker'a istediğim şeyi bitirdikten sonra bu konteyneri kapatmasını söyler. Ve sonra konteyner için tanımlayıcıyı koyuyorum, bu sadece bir URI.

Ve sonra sonunda, Docker'ın bu imajdan oluşturulan konteyner içinde çalıştırmasını istediğim komutu belirtiyorum. Sadece cowpy diyeceğim, bu Conda Forge'dan yüklenen aracın adı, imajın içinde mevcut.

Enter'a basacağım ve işte. Bir sistemde cowpy çalıştırdık. Bize bazı bilgiler veren küçük bir ineğimiz var.

Şimdi cowpy'nin yerel sistemimde kurulu olmadığını unutmayın. Yani sadece tüm Docker şeyleri olmadan çalıştırırsam, komut bulunamadı diyor. Yani bu bir imaj çekti. Docker kullanarak bir konteyner oluşturdu ve sonra o konteynerin içine girdi ve bu komutu bizim için çalıştırdı ve çıktıyı terminalimize geri verdi. Çok, çok havalı.

## 1.3. Konteyneri etkileşimli olarak cowpy çalıştırmak için kullanma

Tamam, şimdi bir adım daha ileri gideceğiz ve bu konteyneri etkileşimli olarak çalıştıracağız ve etrafı kurcalayacağız, böylece konteyner içinde neler olduğunu görebiliriz.

Yani yukarı gidip çalıştırma komutumu alırsam ve sonundaki cowpy'yi kaldıracağım, çünkü aslında cowpy çalıştırmak istemiyorum. Bir Bash terminali çalıştırmak istiyorum.

Ve sonra buraya geri döneceğim ve "-it" yapacağım, bu Interactive ve Terminal veya TTY anlamına gelir ve enter'a basacağım.

Ve şimdi istem, yazdığım yerden önceki kısmın değiştiğini görebilirsiniz. Bu, dizini söyleyen Codespaces istemiydi ve şimdi base ve root ve tmp diyor. Yani şimdi konteynerin içindeyim ve "ls" yaparsam, bu dizinde gördüğüm dosyaların çalışma alanımdaki dosyalardan farklı olduğunu göreceksiniz.

Ve aslında, yerel codespaces çalışma alanımdan veya Docker konteyneri içindeki yerel sürücümden hiçbir dosyayı göremiyorum. Docker konteyner çalışma zamanı tamamen izole edilmiştir ve dışarıdaki bir ana bilgisayar dosya sisteminden herhangi bir dosya yazamaz veya okuyamaz.

Ancak, konteyner içinde kurulu olan yazılımı görebilir ve çalıştırabilirim. Yani cowpy çalıştırabilirim ve cowpy'yi nasıl kullanacağımız hakkında biraz daha fazla şey görebiliriz. Burada "cowpy 'Hello World'" yapabilirim ve bu, alıntımı aslında küçük bir konuşma balonunun içine koymasını söyler. Ve ayrıca farklı türde inekler çalıştırabilirsiniz, yani bir inek olmak zorunda değil. Bir "-c" yapabilirsiniz. Ve ben İsveç'teyim, bu yüzden bir geyik seçeceğim. Çok güzel. Ona boynuzlar verdi.

Ve oynayabileceğiniz bir sürü farklı tane var, bunları eğitim dokümanlarında açıklandığını görebilirsiniz.

## 1.3.4. Verileri konteynere bağlama

Tamam. Dosya sistemimizdeki dosyalar üzerinde cowpy çalıştırabilseydik güzel olurdu.

Tabii ki, sadece konteynere sahip olmak ve hiçbir şeye erişim olmamak çok kullanışlı değil. Güvenli ve tekrarlanabilir olabilir, ama çok kullanışlı değil.

Peki bunu nasıl yaparız? Bu Docker konteynerinden exit yazarak çıkacağım ve istemin bize şimdi normal Codespaces'imize geri döndüğümüzü söylediğini görebilirsiniz.

Ve aynı komutu tekrar çalıştıracağım. Ama bu sefer buraya bazı ek bayraklar ekleyeceğim. Ve önemlisi "-v", bu bir volume bağlamak anlamına gelir, bu temelde bir disk alanının bir parçası gibi.

"-v" iki bölüm alır: bir string var ve sonra iki nokta üst üste ve bir string. Ve ilk kısım, konteynere bağlanması gereken yerel dosya sistemidir. Ve sonra ikinci kısım, bunun konteyner içinde nerede bitmesi gerektiğidir.

Şimdi sadece tüm yerel dosya sistemimi buraya yüklemek istiyorum. Yani "." mevcut çalışma dizinidir. Yani sadece "." yapacağım ve sonra ":", ve sonra bunu konteyner içinde "my_project" adlı yeni bir dizine koyacağız. Bu gerçekten herhangi bir şey olarak adlandırılabilir.

Ve sonra tekrar çalıştıracağım.

Bırakıldığım çalışma dizininde, /tmp, dosyalar orada değil. Ama "ls my_project" yaparsam, işte burada: yerel olarak Codespaces'imizde sahip olduğumuz tüm aynı dosyalar artık o yolda konteyner içinde mevcut.

Bu okuma ve yazma erişimidir, bu yüzden bu dizinde yeni dosyalar oluşturabilirim ve bunlar ana bilgisayar dosya sistemimde görünecektir. Bu özel dizin, o zaman konteynerin dışındaymışım gibi tam olarak davranır, böylece şimdi okuyabilir, yazabilir ve şeyler yapabilirim.

## 1.3.5. Bağlı verileri kullanma

Tamam, bunu yapabileceğimizi kanıtlayalım. "cat /my_project/data/greetings.csv" yapıyorum. Hatırlarsanız bu dosya içeriği böyle görünüyor. Şimdi bunu cowpy'ye aktarabilirim ve inek o dosyanın farklı çıktılarını küçük konuşma balonunda yazdıracak, bu biraz eğlenceli.

Yani görebilirsiniz, artık konteynerdeki yazılımı ana bilgisayar sistemimizdeki dosyalarla etkileşim kurmak için kullanabiliriz.

Tamam, geri çıkalım ve eğitim materyalinin geri kalanına devam edelim.

## 2. Nextflow'da konteyner kullanma

Yani konteyner kullanmak gerçekten havalı. Umarım bu mantıklıdır. Ve bu konteynerlerin değerini ve bunun analiz yazılımı çalıştırmak için neden yararlı olduğunu görebilirsiniz.

Ama tüm bu aynı süreci Nextflow içinde nasıl yaparız? Kendimiz bir sürü Docker komutu çalıştırmak istemiyoruz. Sadece Nextflow'un tüm bunları bizim için halletmesini istiyoruz.

Hadi bunu çözelim. Boru hattımıza cowpy çalıştırmak için yeni bir süreç ekleyeceğiz. Tamam, yeni sürecimiz için yeni bir modül oluşturalım. Yani modüllere gidin, buna cowPy.nf diyelim ve sonra kodu buradan eğitim materyalinden kopyalayacağım.

Ama görebilirsiniz süreç çok basit. Şimdiye kadar yaptıklarımıza çok benziyor, bir path ile bir input bloğumuz var, bu bizim girdi dosyamız ve ayrıca burada bir value var, böylece bu bir karakter olacak, istediğimizde tekrar bir geyik kullanabiliriz.

Ve sonra bir output, burada tek bir dosya, bir path ve sonra bir script. Ve konteynerin içinde etkileşimli olarak yaptığımız aynı şeyi yapıyoruz: girdi dosyasını okumak için "cat" yapıyoruz. O içeriği cowpy'ye aktarıyoruz. O girdiye göre belirli bir karakter seçiyoruz, cowpy girdi dosyası adlı bir çıktı dosyasına yazıyoruz, bu daha sonra çıktıya yansıtılıyor.

Harika. Bunu dahil edelim. Yani include \{ cowpy \} from "./modules/cowpy.nf", buna cowpy mi dedim? Evet.

Ve sonra yeni sürecimizi iş akışının ana bloğunda aşağıda çağıralım. Yani cowpy çalıştıralım. Ve yeni cowpy sürecimizi alacağız ve collectGreetings.out diyeceğiz.

Ve hatırlarsanız, bu modül için iki çıktı vardı. Biri outfile ve biri report adında. VS Code uzantısı bunları bizim için otomatik öneriyor ve .outfile istiyoruz.

Her zaman buraya bu sürece atlayabilirsiniz. Ya üzerine gelirsiniz ve çıktıların ne olduğunu hızlıca göstermelidir. Ve ayrıca içine komut tıklayabiliriz ve daha ayrıntılı görmek istiyorsanız modül dosyasını açacaktır.

İşte buradayız. Oradaki outfile bu ve bu path. Yani bu artık cowpy sürecimiz için girdi dosyası olacak. Harika.

Hatırlarsanız, bir cowpy sürecinin iki girdisi var. Ayrıca karakter için value kanalımız vardı. Yani buraya "params.character" ekleyebiliriz. Bunu sabit kodlayabilirdim, ama bunu bir CLI seçeneği yapalım, böylece tire, tire character yapabiliriz.

Doğru. Şimdi az önce çağırdığımız girdi parametresini tanımlamam ve ona bir varsayılan vermem gerekiyor. Yani character, String. Ve geyiği seviyorum, bu yüzden varsayılan olarak geyik olarak ayarlayacağım.

Doğru, çalıştırmayı deneyelim. Yani Nextflow run hello containers yaparsam, ne olacağını göreceğiz.

Eski çalışma dizinlerim etrafta dolaşıyorsa tire resume kullanabilirdim. Ve yine, bu ilk süreçler önbelleğe alınmış olurdu ve biraz daha hızlı olurdu, ama temelde aynı olmalı.

Şimdi hemen görebiliriz, yeni sürecimize geldiğinde bir hata fırlattı, burada cowpy sürecini yürütürken bir hata olduğunu söylüyor ve 127 çıkış durumuyla çıktı. Çalıştırmaya çalıştığı komut bu. Doğru görünüyor, beklediğimiz gibi görünüyor. O çıktı dosya adını alıyor, doğru görünüyor, bir geyik karakteriyle çalıştırıyor ve kaydetmeye çalışıyor.

Ama buradaki komut hatasının cowpy komutunun bulunamadığını söylediğini görebilirsiniz. Ve bu mantıklı çünkü aslında Nextflow'a henüz bir konteyner kullanmasını söylemedik. Sadece ona cowpy komutunu verdik. Ve daha önce söylediğim gibi, cowpy yerel sistemimizde kurulu değil. Yani çalıştırmaya çalıştığında başarısız oldu.

## 2.3.1. cowpy için bir konteyner belirtme

Nextflow'a kullanılabilir bir konteyner olduğunu ve onu kullanabileceğini söylememiz gerekiyor. Peki bunu nasıl yaparız?

Modülümüze girersek, en üste "container" adlı yeni bir bildirim ekleyeceğiz. Ve sonra bunu bir string'e ayarlayacağız.

Şimdi, hatırlarsanız, Seqera Containers'da, o URL'yi kopyalayabilirim ve onu buraya tırnak içine bırakıyorum.

Geri dön ve tekrar çalıştırmayı dene.

Bu sefer çalışıp çalışmadığını görelim.

Ne yazık ki, süreç için bir konteyner tanımlamış olmamıza rağmen tamamen aynı şekilde başarısız oluyor. Yani docker imajımızı kullanmak için, iş akışını çalıştırdığımızda Nextflow'a Docker kullanımını etkinleştirmesini söylememiz gerekiyor.

Ve bunu yeni bir yapılandırma dosyası oluşturarak yapacağız. Yani touch nextflow.config diyeceğim.

Bu, boru hattını başlatırken çalışma dizinindeyse otomatik olarak yüklenecek özel bir dosya adıdır. Yani bu Nextflow dot config dosyasına girersem, aslında zaten var olduğunu görebilirsiniz, bunu unutmuştum. Ve burada zaten docker.enabled var, ama false olarak ayarlanmış, bu varsayılan.

Yani bunu yerine True'ya eşit olarak değiştirirsem, docker.enabled. Ve bunların hepsi için Nextflow dokümanlarında referans dokümanları var. Ve ayrıca VS Code uzantısıyla üzerine geldiğimde görebilirsiniz, buna özgü dokümanları çeker ve bunun ne anlama geldiğini ve nasıl ayarlanacağını söyler.

Yani şimdi true olarak ayarladık ve Nextflow'u tekrar çalıştırırsam, Nextflow artık yerel olarak yoksa o docker imajını bizim için getirmeyi bilecek ve sonra o süreç o konteyner ortamıyla yürütecek.

Ve başarıyla çalıştığını görebiliriz ve cowpy'nin yanında küçük bir tik var. Harika. Yukarı gidip sonuçlar dizinine bakarsam, dosya henüz orada değil. Ve bunun nedeni, diğerlerinin hepsinde olduğu gibi bu çıktı dosyasını hala yayınlamamız gerekiyor.

Yani iş akışı içindeki yayınlama bloğuna gidiyoruz, mycowpy equals cowpy.out diyoruz.

Ve sonra burada çıktı bloğunda, mycowpy, süslü parantezler path. Hata. Hello containers. Mode, copy.

Şimdi tekrar çalıştırırsam, tamamen aynı şekilde çalışmalı. Tire resume kullanabilirdim ve her seferinde unutuyorum. Ve sonra yukarı gidiyorum ve şimdi cowpy-COLLECTED adlı yeni bir dosya oluşturuldu ve işte BONJOUR, HELLO, HOLA diyen geyiğim. Harika.

Şimdi tabii ki şimdi "--character" de geçebilirim. Farklı seçenekler neler? Sanırım bir Hindi var? Yani character Turkey kullanabilirim. Tamamen aynı şekilde çalışacak. Tire resume kullanmak için başka bir fırsatı kaçırdım ve şimdi dosyamızı yüklersek ve şimdi bir Hindimiz var. Harika.

## 2.3.4. Nextflow'un konteynerleştirilmiş görevi nasıl başlattığını inceleme

Tamam. Son küçük şey. Hadi bu komutu tekrar hızlıca çalıştıralım, bu sefer resume, ve tüm bunların bizim için çalışmasını sağlamak için Nextflow'un perde arkasında ne yaptığına hızlıca bir göz atalım.

Bu sefer süper hızlı, hadi bu çalışma dizinine gidelim, cd work/. Şimdi hatırlarsanız burada bir sürü nokta dosyamız var ve bu durumda ilgilendiğimiz, neredeyse hiç bakmamız gerekmediğini söylediğim .command.run adlı dosya.

Code dot command run yaparsam, editörde açılacak. Ve bu dosyada arama yapabilirim ve aşağı kaydırırsam Docker run'ı görmeliyim. Ve Nextflow'un docker run komutunu bizim için yaptığını görebilirsiniz, yapılandırmada Docker etkinleştirildiğinde. Burada bir sürü farklı bayrak ve şey var, ama kendimiz kullandığımızda kullandığımız "-v" bayrağını görebilirsiniz. Ve yerel çalışma alanı dizinini konteynere bağladığını görebilirsiniz, böylece konteyner girdi dosyalarımıza erişebilir ve çıktıları kaydedebilir. Ve sonunda, ayrıca .command.sh çalıştırıyor, bu içinde cowpy komutu olan oluşturulan betiktir.

Ve böylece Nextflow'un iş akışı mantığını aldığını görebilirsiniz, bu gerçekten önemsediğimiz şey, analizimize özgü olan ve Docker'ın sistemimizde çalışmasını sağlamak için tüm akıllı perde arkası işlerini yapıyor.

Ve bunu gerçekten taşınabilir bir şekilde yapıyor, böylece boru hattının son kullanıcısı kullandıkları teknolojiyi değiştirebilir: Docker, Singularity, Apptainer, Conda. Bu gerçekten boru hattı mantığı için önemli değil, ama Nextflow tüm temel altyapı ihtiyaçlarını halledecek, böylece her yerde çalışır.

Ve bu gerçekten Nextflow'un süper gücüdür. Tekrarlanabilirlik ve taşınabilirlik. Ve Nextflow ile iş akışınızı gerçekten paylaşabilirsiniz ve diğer insanlar onu kendi sistemlerinde çalıştırabilir ve sadece çalışacaktır.

Bu gerçekten, gerçekten zor bir şey ve şimdi bunu iş akışlarınızla nasıl yapacağınızı da biliyorsunuz.

Tamam, bu bölüm için bu kadar. Kursun sonuna giderseniz, konteynerler hakkında tekrar bir sınav bulacaksınız. Umarım hepsi mantıklıydı. Analizle çalışmanın gerçekten havalı bir yolu. Ve konteynerlerde yeniyseniz, umarım bunun gidilecek yol olduğuna sizi ikna etmişimdir ve asla geriye bakmayacaksınız.

Ama bununla, belki biraz ara verin ve birkaç dakika içinde Hello Nextflow'un son altıncı bölümünü geçmek için bana katılın, bu tamamen yapılandırma hakkında.

Çok teşekkür ederim.
