# Bölüm 5: Hello Containers - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım eksiksiz talimatlar için [kurs materyaline](../05_hello_containers.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş geldiniz ve arka plan

Merhaba ve Hello Nextflow'a tekrar hoş geldiniz. Bu, Hello Containers adlı beşinci bölüm. Kursun bu bölümünde, bir boru hattı için yazılım gereksinimlerini nasıl kapsülleyeceğimizi konuşacağız; böylece boru hattını çalıştıran kişilerin yazılımı yüklemeyi düşünmelerine gerek kalmayacak.

Benim kadar uzun süredir biyoinformatikte çalışıyorsanız, sık sık "kötü eski günler" dediğim zamanları hatırlayabilirsiniz. Başkasının boru hattını çalıştırmak veya çalışmalarını tekrarlamak istediğinizde, kullandıkları tüm farklı yazılım araçlarını aynı sürümlerde makinenizde derlemek için saatler veya günler harcardınız; bu gerçek bir kabusту. Gerçekten zordu.

Bir HPC üzerinde çalışıyorsanız, sistem yöneticilerinin sizin için yazılım yüklemeye çalıştığı ortam modüllerini kullanmış olabilirsiniz. Bu yöntem işe yarıyordu, ama yine de kusurlu bir çözümdü.

Ama şimdi bunu yapmanın daha iyi yolları var. Nextflow'un farklı yazılım konteyner teknolojileri için yerleşik desteği bulunuyor. Docker en yaygın olanıdır; bugün kullanacağımız da bu. Codespaces'te iyi çalışır, yerel bilgisayarınızda iyi çalışır ve bulutta da iyi çalışır.

Bunların yanı sıra, HPC sistemlerinde çok yaygın olan ve temelde tamamen aynı şekilde çalışan Singularity veya Apptainer da mevcuttur. Podman, Shifter ve bunlara çok benzeyen başka araçlar da var.

Nextflow'un desteklediği, benzer ama tam olarak aynı olmayan bir diğer seçenek ise Conda'dır. Nextflow, süreç bazında Conda ortamlarını sizin için yönetebilir; bu, kendi Conda ortamlarınızı kendiniz oluşturmaktan çok daha iyidir. Üstelik bir boru hattıyla birlikte gönderilebilir.

Bu bölüme konteyner teknolojileri, Docker ve bunların nasıl çalıştığı hakkında biraz konuşarak başlayacağız. İlk yarıyı yalnızca Docker'da manuel olarak yapacağız; böylece perde arkasında neler olduğunu ve bunun nasıl çalıştığını anlayacaksınız. Nextflow'un ne yaptığını ve iş akışınızın çalıştırıldığında ne yaptığını anlamak için bu gerçekten önemlidir.

Hadi Codespaces'imize geçelim. Her şeyi tekrar temizledim, ama Hello Containers'a atlarsanız tüm betiklerimizin ve her şeyin modüller bölümünün sonundaki gibi orada olduğunu görmelisiniz. Yani modüller dizininde oluşturduğum farklı modüllerimiz var.

Onlar hâlâ orada. Çalışabilmesi için orada olmaları gerekiyor. İş akışı ve çıktı da tamamen aynı; tek fark, dosyalarınızın o dizinde sonlanması için çıktı yayımlama yolunu Hello Containers olarak değiştirmiş olmamız.

İsterseniz çalıştığını kontrol etmek için bunu şimdi çalıştırabilirsiniz ya da terminale devam edebiliriz.

## 1. Bir konteyneri 'manuel olarak' kullanma

Konteynerlerimizi yönetmek için Docker kullanacağız. Codespaces'imde kurulu olup olmadığını `docker -v` komutuyla kontrol edebilirim; bu komut kurulu olan sürümü gösterir ve her şeyin düzgün çalıştığını doğrular.

Konteynerler ve Docker'da gerçekten önemli olan iki kavram var. Birine image, diğerine container denir. Image, kullanacağınız tüm dosya sisteminin anlık görüntüsü gibidir; container ise çalışan ortamdır. Bir image kullanarak container oluşturursunuz.

O container'ın içine girdiğinizde, tipik olarak tam bir işletim sistemi gibi çalışır. Dış dünyadan kesilmiştir; her şeyden ayrıdır ve bu iyi bir şeydir. Nextflow ile bu kadar iyi bir tekrarlanabilirlik elde etmemizin nedeni de budur.

Bir container içinde çalışan görevler, yerel sisteminizdeki herhangi bir yapılandırma dosyasından etkilenmez. Diğer dış etkilerden bağımsız olarak kendi küçük korumalı alanlarında çalışırlar. Dosyalar daha sonra son derece tekrarlanabilir bir şekilde üretilir; çünkü farklı bilgi işlem ortamlarında çalışan herkes için tamamen aynı temel kütüphaneleri, aynı bağımlılıkları ve tamamen aynı yazılımı kullanıyorsunuz. Açıkçası bunun harika ve muhteşem olduğunu düşünüyorum; hâlâ bunun mümkün olduğu gerçeği aklımı başımdan alıyor.

## 1.1. Konteyner imajını çekme

Docker imajlarını ve Docker'ı deneyeceğiz. Sisteminizde çalıştırdığınızda Docker, bilgisayarınızda veya bu durumda kod alanında, geçmişte indirilen ve kullanılan tüm farklı imajları ve üzerine inşa edildikleri farklı katmanları takip eden bir docker kayıt defteri tutar.

Docker ile yerel olarak hangi imajlara sahip olduğumuzu `docker image ls` komutuyla görebiliriz. Bu durumda burada bir sürü Docker imajı olduğunu görebilirsiniz; bunların hepsi bu Codespaces'i kurmakla ilgilidir. Hepsi dev container'lar ve benzeri şeylerle ilgili. Bu yüzden onlar hakkında çok fazla endişelenmenize gerek yok. Ama kurs ilerledikçe daha fazla imaj ekleyip indirdikçe o listeyi kontrol edebilirsiniz; yerel kayıt defterinin çektiğimiz tüm bu şeyleri takip ettiğini göreceksiniz.

Yeni bir tane almak için `docker pull` komutunu kullanacağız. Bu komut Docker'a web'den yeni bir imaj getirmesini söyler.

Ardından o container için URI'yi giriyoruz. Bu, yerel olarak oluşturup internete gönderdiğiniz bir docker imajı olabilir ya da başka birinin yaptığı bir imaj olabilir. Docker imajı yapmanın pek çok farklı yolu var, ama muhtemelen en basit yollardan biri bunu dış kaynak kullanarak başka birine yaptırmaktır.

Bu eğitimde kullanacağımız şey, Seqera'dan Seqera Containers adlı bir hizmettir.

Seqera Containers tamamen ücretsizdir ve Wave adını verdiğimiz, Nextflow'a tamamlayıcı bir şekilde konteynerleri yönetmek için geliştirilmiş açık kaynaklı bir yazılım parçası kullanır. Nextflow ile sık karşılaştığımız pek çok yaygın kullanım durumunu ele alır.

İhtiyacımız olan yazılımın Conda'da, Bioconda'da veya conda-forge kanallarında ya da daha alana özel kanallarda paketlenmiş olması çok yaygındır. Wave ve Seqera Containers, bu kaynaklardan imaj oluşturmada gerçekten başarılıdır.

Bu web arayüzüne gidebilirim; "cowpy" adlı bir paketle oynayacağız. İstediğim paketin adını yazıyorum. Arama yapıyor; Python paket dizininde buldu, bu yüzden bunu kullanabilirim. Biraz daha beklersem bioconda ve conda-forge'u da arıyor. Gördüğünüz gibi, burada herhangi bir conda kanalı belirleyebilirim. Nvidia kanalı veya başka bir şey bulmak isterseniz, bu da çalışmalı.

Ardından benim için bir docker imajı mı yoksa singularity imajı mı oluşturmasını istediğimi ve hangi CPU mimarisini kullanmak istediğimi belirleyebilirim: amd64 veya arm64.

Bioconda sonuçları listelendikten sonra mevcut olan tüm farklı sürümleri de görebiliyorum. Bunu seçeceğim. İstersem aramaya devam edip Conda'dan daha fazla paket alabilir ve bu container'ı istediğim gibi oluşturabilirdim, ama yalnızca bunu istiyorum. Bu yüzden Get Container'a tıklayacağım.

Başka biri daha önce aynı container'ı istemiş ve bir kayıt defterinden döndürülmüş, bu yüzden hemen alıyoruz. Ama başka hiç kimse bu yazılım paketini veya bu yazılım paketi kombinasyonunu hiç istememişse, Wave ve Seqera Containers bunu bizim için anında oluşturur.

Bu URL'yi kopyalayabiliriz ve ayrıca yapı ayrıntılarını görebiliriz. Bu, hizmetin arka planda ne yaptığını gösterir: bir conda ortam dosyası oluşturdu, bir docker dosyası oluşturdu ve ardından docker yapı sürecini çalıştırdı. Ayrıca bir güvenlik taraması çalıştırdı; böylece herhangi bir CVE'yi görebilirsiniz. Ne zaman oluşturulduğunu da söyler.

Wave ve Seqera Containers bundan çok daha fazlasını yapabilir, ama bu en yaygın olan basit bir kullanım durumudur. Bu imajların en az beş yıl boyunca barındırıldığını da belirtmeliyim. Yani bu URL'leri boru hatlarınıza yerleştirebilir ve yakın zamanda kaybolmayacaklarını bilebilirsiniz.

Böylece cowpy için docker imajımın URL'sini aldım.

Şimdi `docker pull` komutuyla o URL'yi çekebilirim; tüm farklı katmanları getirip bu imajı yerel olarak kullanılabilir hale indirecek.

## 1.2. Konteyneri cowpy'yi tek seferlik komut olarak çalıştırmak için kullanma

Tamam, şimdi onu gerçekten kullanmayı deneyelim. `docker pull` yerine `docker run` komutunu kullanacağım ve `--rm` bayrağını ekleyeceğim; bu yalnızca Docker'a istediğim işlemi bitirdiğinde bu container'ı kapatmasını söyler. Ardından container için tanımlayıcıyı giriyorum; bu yalnızca bir URI.

Sonunda, Docker'ın bu imajdan oluşturulan container içinde çalıştırmasını istediğim komutu belirtiyorum. Yalnızca `cowpy` diyeceğim; bu, conda-forge'dan yüklenen ve imaj içinde mevcut olan aracın adıdır.

Enter'a basıyorum ve işte buyurun. Bir sistemde cowpy'yi çalıştırdık. Bize biraz bilgi veren küçük bir ineğimiz var.

Dikkat edin: cowpy yerel sistemimde kurulu değil. Yani tüm Docker şeyleri olmadan çalıştırırsam "command not found" diyor. Demek ki bu bir imaj çekti, Docker kullanarak bir container oluşturdu, ardından o container'ın içine girdi ve bu komutu bizim için çalıştırdı; çıktıyı terminalimize geri verdi. Çok, çok havalı.

## 1.3. Konteyneri cowpy'yi interaktif olarak çalıştırmak için kullanma

Tamam, şimdi bir adım daha ileri gideceğiz. Bu container'ı interaktif olarak çalıştıracak ve etrafa bir göz atacağız; böylece container içinde ne olduğunu görebileceğiz.

Çalıştırma komutuma geri dönüp sondaki `cowpy`'yi çıkaracağım; çünkü aslında cowpy çalıştırmak istemiyorum. Bir Bash terminali çalıştırmak istiyorum.

Buraya geri gelip `-it` ekleyeceğim; bu Interactive ve Terminal veya TTY anlamına gelir. Enter'a basıyorum.

Şimdi istem değişti; yazdığım kısımdan önceki bölüm artık farklı görünüyor. Bu, dizini gösteren Codespaces istemiydi; şimdi ise base, root ve tmp yazıyor. Yani artık container'ın içindeyim. `ls` yaparsam, bu dizinde gördüğüm dosyaların çalışma alanımdaki dosyalardan farklı olduğunu göreceksiniz.

Aslında, yerel Codespaces çalışma alanımdan veya yerel sürücümden hiçbir dosyayı Docker container'ı içinde göremiyorum. Docker container çalışma zamanı tamamen izole edilmiştir; dışarıdaki bir host dosya sisteminden hiçbir dosyayı okuyamaz veya yazamaz.

Ancak container içinde yüklü olan yazılımı görebilir ve çalıştırabilirim. Yani cowpy'yi çalıştırabilirim ve cowpy'nin nasıl kullanılacağı hakkında biraz daha şey görebiliriz. Burada `cowpy 'Hello World'` yapabilirim; bu, alıntımı küçük bir konuşma balonunun içine koymasını söyler. Farklı türde inekler de çalıştırabilirsiniz; inek olmak zorunda değil. `-c` yapabilirsiniz. Ben İsveç'teyim, bu yüzden bir geyik seçeceğim. Çok güzel; ona boynuz verdi.

Eğitim dokümanlarında açıklandığını görebileceğiniz oynayabileceğiniz bir sürü farklı seçenek var.

## 1.3.4. Konteynere veri bağlama

Tamam. Dosya sistemimizdeki dosyalar üzerinde cowpy çalıştırabilseydik güzel olurdu.

Tabii ki, yalnızca container'a sahip olmak ve hiçbir şeye erişim olmaması çok kullanışlı değil. Güvenli ve tekrarlanabilir olabilir, ama çok kullanışlı değil.

Peki bunu nasıl yaparız? `exit` yazarak bu Docker container'ından çıkacağım; istem bize artık tekrar normal Codespaces'imizde olduğumuzu söylüyor.

Aynı komutu tekrar çalıştıracağım. Ama bu sefer buraya bazı ek bayraklar ekleyeceğim. Önemli olan `-v` bayrağıdır; bu bir volume bağlamak anlamına gelir; temelde bir disk alanının bir parçası gibidir.

`-v` iki kısım alır: bir string, ardından iki nokta üst üste ve bir string daha. İlk kısım container'a bağlanması gereken yerel dosya sistemidir. İkinci kısım ise bunun container içinde nereye gitmesi gerektiğidir.

Tüm yerel dosya sistemimi buraya yüklemek istiyorum. `.` mevcut çalışma dizinidir. Yani yalnızca `.` yapacağım, ardından `:` ve bunu container içinde `my_project` adlı yeni bir dizine koyacağız. Bu gerçekten herhangi bir şey olarak adlandırılabilirdi.

Sonra tekrar çalıştıracağım.

Bırakıldığım çalışma dizininde, yani `/tmp`'de, dosyalar orada değil. Ama `ls my_project` yaparsam, işte buyurun: Codespaces'te yerel olarak sahip olduğumuz tüm aynı dosyalar artık o yoldaki container içinde mevcut.

Bu okuma ve yazma erişimidir; bu dizinde yeni dosyalar oluşturabilirim ve bunlar host dosya sistemimde görünecek. Bu özel dizin, container dışındaymışım gibi tam olarak davranır; artık okuyabilir, yazabilir ve işlemler yapabilirim.

## 1.3.5. Bağlanan veriyi kullanma

Tamam, bunu yapabileceğimizi kanıtlayalım. `cat /my_project/data/greetings.csv` komutunu çalıştırıyorum. Hatırlarsanız bu dosyanın içeriği böyle görünür. Bunu cowpy'ye aktarabilirim; inek o dosyanın farklı çıktılarını küçük konuşma balonunda yazdıracak; bu biraz eğlenceli.

Gördüğünüz gibi, artık container'daki yazılımı host sistemimizin dosyalarıyla etkileşim kurmak için kullanabiliyoruz.

Tamam, tekrar dışarı çıkalım ve eğitim materyalinin geri kalanıyla devam edelim.

## 2. Nextflow'da konteyner kullanma

Konteyner kullanmak gerçekten havalı. Umarım mantıklıdır. Bu container'ların değerini ve analiz yazılımı çalıştırmak için neden kullanışlı olduklarını görebilirsiniz.

Ama tüm bu aynı süreci Nextflow içinde nasıl yaparız? Kendimiz bir sürü Docker komutu çalıştırmak istemiyoruz. Nextflow'un tüm bunları bizim için halletmesini istiyoruz.

Öyleyse bunu çözelim. Boru hattımıza cowpy çalıştırmak için yeni bir süreç ekleyeceğiz. Yeni sürecimiz için yeni bir modül oluşturalım. Modüller dizinine gidin, buna `cowPy.nf` diyelim ve ardından eğitim materyalinden kodu kopyalayacağım.

Sürecin çok basit olduğunu görebilirsiniz. Şimdiye kadar yaptıklarımıza çok benziyor: girdi dosyamız olan bir `path` ile bir girdi bloğumuz var ve ayrıca burada bir `value` var; bu bir karakter olacak, yani isterseniz yine bir geyik kullanabiliriz.

Ardından bir çıktı var; bu burada tek bir dosya, bir `path` ve sonra bir `script`. Container içinde interaktif olarak yaptığımız şeyin aynısını yapıyoruz: girdi dosyasını okumak için `cat` kullanıyoruz, o içerikleri cowpy'ye aktarıyoruz, o girdiye göre belirli bir karakter seçiyoruz ve cowpy girdi dosyası adlı bir çıktı dosyasına yazıyoruz; bu daha sonra çıktıya yansıtılıyor.

Harika. Bunu dahil edelim. Yani `include { cowpy } from "./modules/cowpy.nf"` — buna cowpy mi dedim? Evet.

Ardından yeni sürecimizi iş akışının ana bloğunda çağıralım. Cowpy'yi çalıştıralım. Yeni cowpy sürecimizi alacağız ve `collectGreetings.out` diyeceğiz.

Hatırlarsanız, bu modül için iki çıktı vardı: biri `outfile`, diğeri `report`. VS Code uzantısı bunları bizim için otomatik öneriyor ve `.outfile`'ı istiyoruz.

Her zaman bu sürece atlayabilirsiniz. Üzerine gelirseniz çıktıların ne olduğunu hızlıca gösterir. Daha fazla ayrıntı görmek isterseniz komut tıklayıp modül dosyasını açabilirsiniz.

İşte buyurun. `outfile` orada ve bu bir `path`. Yani bu artık cowpy sürecimiz için girdi dosyası olacak. Harika.

Hatırlarsanız, bir cowpy sürecinin iki girdisi var. Karakter için `value` kanalı da vardı. Yani buraya `params.character` ekleyebiliriz. İstersem bunu sabit kodlayabilirdim, ama bir CLI seçeneği yapalım; böylece `--character` kullanabiliriz.

Doğru. Şimdi az önce çağırdığımız girdi parametresini tanımlamam ve ona bir varsayılan vermem gerekiyor. Yani `character`, `String`. Geyiği seviyorum, bu yüzden varsayılan olarak `moose`'a ayarlayacağım.

Doğru, çalıştırmayı deneyelim. `nextflow run hello containers` yaparsam ne olduğunu göreceğiz.

Eski çalışma dizinleri etrafta dolaşsaydı `--resume` kullanabilirdim. Yine, bu ilk süreçler önbelleğe alınmış olurdu ve biraz daha hızlı olurdu, ama temelde aynı olmalı.

Şimdi hemen görebiliriz ki yeni sürecimize geldiğinde bir hata fırlattı; bize cowpy sürecini çalıştırırken bir hata olduğunu ve 127 çıkış durumuyla çıktığını söylüyor. Çalıştırmaya çalıştığı komut bu. Doğru görünüyor, beklediğimiz gibi: çıktı dosya adını alıyor, moose karakteriyle çalıştırıyor ve kaydetmeye çalışıyor.

Ama komut hatasının cowpy komutunun bulunamadığını söylediğini görebilirsiniz. Bu mantıklı; çünkü henüz Nextflow'a bir container kullanmasını söylemedik. Yalnızca ona cowpy komutunu verdik. Daha önce söylediğim gibi, cowpy yerel sistemimizde kurulu değil. Yani çalıştırmayı denediğinde başarısız oldu.

## 2.3.1. cowpy için bir konteyner belirtme

Nextflow'a mevcut bir container olduğunu ve kullanabileceğini söylememiz gerekiyor. Peki bunu nasıl yaparız?

Modülümüze girelim; üste `container` adlı yeni bir yönerge ekleyeceğiz ve bunu bir string'e ayarlayacağız.

Hatırlarsanız, Seqera Containers'da o URL'yi kopyalayabilirim ve buradaki tırnak işaretlerinin içine bırakıyorum.

Şimdi geri dönün ve tekrar çalıştırmayı deneyin.

Bakalım bu sefer çalışacak mı.

Ne yazık ki, süreç için bir container tanımlamış olmamıza rağmen tamamen aynı şekilde başarısız oluyor. Docker imajımızı kullanabilmek için, iş akışını çalıştırırken Nextflow'a Docker kullanımını etkinleştirmesini söylememiz gerekiyor.

Bunu yeni bir yapılandırma dosyası oluşturarak yapacağız. `touch nextflow.config` diyeceğim.

Bu, boru hattını başlatırken çalışma dizininde bulunursa otomatik olarak yüklenecek özel bir dosya adıdır. Bu `nextflow.config` dosyasına girersem, aslında zaten var olduğunu görebilirsiniz; bunu unutmuştum. Burada `docker.enabled` zaten var, ama `false`'a ayarlı; bu varsayılan değerdir.

Bunu `true`'ya eşit olacak şekilde değiştirirsem, `docker.enabled`. Bu yapılandırma kapsamlarının tümü için Nextflow dokümanlarında referans dokümanları var. Ayrıca VS Code uzantısıyla üzerine geldiğimde, bununla ilgili dokümanları çekiyor ve bunun ne anlama geldiğini ve nasıl ayarlanacağını söylüyor.

Şimdi `true`'ya ayarladık. Nextflow'u tekrar çalıştırırsam, Nextflow artık o docker imajını yerel olarak henüz yoksa bizim için çekmeyi ve ardından o süreci o container ortamıyla çalıştırmayı bilecek.

Başarıyla çalıştığını ve cowpy'nin yanında küçük bir tik işareti olduğunu görebiliriz. Harika. Yukarı çıkıp sonuçlar dizinine bakarsam, dosya henüz orada değil. Bunun nedeni, diğerleri gibi bu çıktı dosyasını da yayımlamamız gerekiyor olması.

Yani iş akışı içindeki yayımlama bloğuna gidiyoruz, `mycowpy equals cowpy.out` diyoruz.

Ardından aşağıdaki çıktı bloğunda, `mycowpy`, küme parantezler `path`. Hoop. Hello containers. Mode, copy.

Şimdi tekrar çalıştırırsam, tamamen aynı şekilde çalışmalı. `--resume` kullanabilirdim ve her seferinde unutuyorum. Yukarı çıkıyorum ve şimdi `cowpy-COLLECTED` adlı yeni bir dosya oluşturuldu; geyiğim BONJOUR, HELLO, HOLA diyor. Harika.

Tabii ki şimdi ayrıca `--character` geçirebilirim. Farklı seçenekler neler? Sanırım bir Turkey var? Yani `character Turkey` kullanabilirim. Tamamen aynı şekilde çalışacak. `--resume` kullanmak için başka bir fırsatı kaçırdım ve şimdi dosyamızı yüklersek bir Turkey'imiz var. Harika.

## 2.3.4. Nextflow'un konteynerize edilmiş görevi nasıl başlattığını inceleme

Tamam. Son küçük şey. Bu komutu tekrar hızlıca çalıştıralım, bu sefer `--resume` ile, ve Nextflow'un tüm bunların bizim için çalışmasını sağlamak için perde arkasında ne yaptığını görmek amacıyla çalışma dizinine hızlıca bir göz atalım.

Bu sefer süper hızlı. Hadi bu çalışma dizinine gidelim: `cd work/`. Hatırlarsanız burada bir sürü nokta dosyamız var; bu durumda ilgilendiğimiz, neredeyse hiç bakmamız gerekmediğini söylediğim `.command.run` adlı olandır.

`code .command.run` yaparsam editörde açacak. Bu dosyada arama yapabilirim; aşağı kaydırırsam `docker run`'ı görmeliyim. Nextflow'un Docker etkinleştirildiğinde yapılandırmada bizim için `docker run` komutunu çalıştırdığını görebilirsiniz. Burada bir sürü farklı bayrak ve şey var, ama kendimiz çalıştırırken kullandığımız `-v` bayrağını görebilirsiniz. Yerel çalışma alanı dizinini container'a bağladığını görebilirsiniz; böylece container girdi dosyalarımıza erişebilir ve çıktıları kaydedebilir. Sonunda, içinde cowpy komutu olan üretilmiş betik olan `.command.sh`'yi de çalıştırıyor.

Böylece Nextflow'un iş akışı mantığını aldığını görebilirsiniz; bu gerçekten önemsediğimiz şeydir, analizimize özgü olandır. Nextflow, sistemimizde Docker'ın çalışmasını sağlamak için tüm akıllı perde arkası işleri yapıyor.

Bunu gerçekten taşınabilir bir şekilde yapıyor; böylece boru hattının son kullanıcısı kullandıkları teknolojiyi değiştirebilir: Docker, Singularity, Apptainer, Conda. Bu, boru hattı mantığı açısından gerçekten önemli değil; ama Nextflow tüm temel altyapı ihtiyaçlarını halledecek, böylece her yerde çalışır.

Ve bu gerçekten Nextflow'un süper gücüdür: tekrarlanabilirlik ve taşınabilirlik. Nextflow ile iş akışınızı gerçekten paylaşabilirsiniz; diğer insanlar kendi sistemlerinde çalıştırabilir ve sadece çalışır.

Bu gerçekten, gerçekten zor bir şeydir ve artık siz de iş akışlarınızla bunu nasıl yapacağınızı biliyorsunuz.

Tamam, bu bölüm için bu kadar. Kursun sonuna giderseniz, konteynerler hakkında yine bir sınav bulacaksınız. Umarım hepsi mantıklıydı. Analizle çalışmanın gerçekten havalı bir yolu. Konteynerlere yeniyseniz, umarım bunun gidilecek yol olduğuna sizi ikna etmişimdir; asla geriye bakmayacaksınız.

Ama bununla birlikte, belki biraz mola verin ve birkaç dakika içinde Hello Nextflow'un tamamen yapılandırma hakkında olan son altıncı bölümünü geçmek için bana katılın.

Çok teşekkür ederim.
