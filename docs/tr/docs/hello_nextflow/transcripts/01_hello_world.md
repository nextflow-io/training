# Bölüm 1: Hello World - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım tam talimatlar için [ders materyaline](../01_hello_world.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgi amaçlıdır ve materyallerdeki tüm bölüm numaralarını kapsamayabilir.

## Hoş Geldiniz

Merhaba ve tekrar hoş geldiniz.

Şu anda "Hello Nextflow" kursunun "Hello World" adlı Birinci Bölümündesiniz. Bu bölümde Nextflow'un temel kavramlarını anlamaya başlayacağız.

Umarım artık Codespaces'te ya da eşdeğer bir ortamda VS Code çalışır durumdadır ve Explorer'da çalışma alanınızda Hello Nextflow klasörü ile tüm bu farklı dosyalar görünüyordur.

Terminalde yalnızca Bash kullanarak çok basit şeyler yaparak başlayacağız; ardından aynı şeyleri Nextflow içinde yapıp yapamayacağımıza bakacağız. Böylece sözdiziminin nasıl göründüğü hakkında bir fikir edineceksiniz.

## 0. Isınma

O hâlde gerçekten basit bir yerden başlayalım. Terminale bir şey yazdırmak için "echo" ile başlayalım: "Hello World". Enter'a basıyorum ve bu terminale gidiyor. Hello World. Umarım bu kursu izleyen herkes için sürpriz olmamıştır.

Tamam, bununla bir şeyler yapalım. Sadece terminale yazdırmak yerine bir dosyaya yazalım. Klavyemde yukarı ok tuşuna basıyorum; bu Bash geçmişini döngüsel olarak gösteriyor ve son komutumu veriyor. Sonuna küçük bir büyüktür işareti ekleyeceğim; bu, komutun çıktısını bir dosyaya yönlendiriyor. Dosyayı output.txt olarak adlandıracağım.

Komutu çalıştırmak için tekrar Enter'a basıyorum; bu sefer terminalde hiçbir şey görünmüyor, ancak sol tarafta output.txt adında yeni bir dosyanın belirdiğini görebiliyoruz.

Bunu terminalde cat gibi bir komutla görüntüleyebiliriz. cat output.txt yazıyorum ve gerçekten "Hello World" yazıyor. Ayrıca üzerine çift tıklayarak VS Code'daki kod düzenleyicide açabiliriz.

## 1.1. Kodu İnceleme

Pekâlâ. Basit olduğunu söylemiştim. Sırada ne var? Bu süreci alıp tekrar deneyelim; ancak bu sefer Nextflow içinde yapalım.

Dediğim gibi, bu kursun tüm farklı bölümleri bir betikle başlıyor ve bu bölümün betiği Hello World olarak adlandırılmış. O hâlde Hello World'ü bulacağım. Tek tıkladığımda önizleme gösteriyor; düzenleyicide açmak için çift tıklıyorum. Terminali de hızlıca kapatıyorum.

Bu çok basit bir betik; olabildiğince sade. Yalnızca 22 satır uzunluğunda ve temelde aynı şeyi yapıyor. Hatta bir kısmı tanıdık gelecektir. Az önce yazdığımız şey bu. Bash komutunun bir dosyaya yönlendirildiğini görebiliyoruz.

Tamam. Başka ne var? Ayrıca bu dosyada Nextflow'un bazı temel kavramlarını görmeye başlayabiliriz. Burada kırmızıyla bir process ve bir workflow görüyoruz. Bunlar Nextflow'daki özel anahtar kelimeler ve özel terminoloji.

## 1.1.1. Process Tanımı

Bir iş akışı içindeki farklı süreçler, iş akışınızın farklı mantıksal birimlerini kapsar. Her süreç tek bir şey yapar.

Çalıştırdığımızda bir görev ya da birden fazla görev oluşturur; bunlar bir boru hattının gerçek yürütme adımlarıdır. Tüm süreçler daha sonra altta gördüğümüz bir workflow bloğu içinde düzenlenir; bu örnekte yalnızca o tek süreci çalıştırır.

Süreç adı bu anahtar kelimeyi takip eder ve bu ad neredeyse her şey olabilir. Sürecin içeriği ise bu süslü parantezlerin içindedir.

Bir process için gerçek anlamda tek gereksinim, bir tür script veya exec bloğu içermesidir. Bu, üçlü tırnak işaretleri içindedir ve boru hattını çalıştırdığımızda çalışma dizinine yazılan Bash betiğidir; bilgisayarınızda veya sunucunuzda gerçekten çalışan şey budur.

Bu genellikle Bash'tir; ancak en üste farklı bir shebang ekleyebilirsiniz ve bu bir Python betiği ya da R betiği olabilir. Fark etmez. Bu betikteki her şey yürütülecektir.

Bu sürece eklediğimiz bir diğer şey ise çıktı bildirimidir. Bu, Nextflow'a bu sürecin output.txt adlı bir çıktı dosyası beklediğini söyler. Bunun bir path olduğunu belirtir; dolayısıyla bir dosya gibi ele alınmalıdır. Örneğin bu val olsaydı, bir değişken veya değer gibi ele alınacağını söylerdi.

Bunun bu dosyayı oluşturmadığını unutmayın. Gerçekte oluşturan, aşağıdaki betiktir. Bu yalnızca Nextflow'a bu dosya adıyla bir çıktı dosyası beklemesini söyler.

## 1.1.2. Workflow Tanımı

Tamam. Altta bir workflow var ve yine bir bildirim görüyoruz. Bu bildirim Main olarak adlandırılmış. Bu, bir betik bloğunun iş akışı karşılığıdır diyebiliriz. İş akışının bir şey yaptığı kısımdır. Bu örnekte sayHello adlı süreci çağırıyoruz.

Elbette boru hattınız normalde bundan çok daha karmaşık görünecektir. Muhtemelen birden fazla süreciniz olacak ve aralarındaki veri akışını düzenlemek için kanallar kullanacaksınız. Bunu bu kursun sonraki bölümlerinde ele alacağız; ancak şimdilik bu yeterli. Bu geçerli bir boru hattı ve çalışması gerekiyor.

VS Code'da buradaki DAG önizlemesine bile tıklayabilirim. DAG veya DAG, boru hattındaki veri akışı yapısının bir temsilidir; yan tarafta bir mermaid diyagramı olarak işlendiğini görebiliriz. Bu örnekte çok basit: iş akışını temsil eden bir kutu ve sayHello adlı bir süreç var; ancak ilerledikçe daha ilginç görünebilir.

## 1.2. İş Akışını Çalıştırma

Tamam, bu iş akışını çalıştırmayı deneyelim ve ne olduğunu görelim.

Altta terminali tekrar açacağım, çıktıyı temizleyeceğim ve Nextflow Run yazacağım. Ardından betik adını yazacağım: hello-world.nf. Enter'a basıyorum.

Tamam, en üstte Nextflow'un çalıştığını, hangi sürümün çalıştığını, betik adını ve diğer bilgileri gösteren standart bilgiler var.

Gerçekten aradığımız önemli şey _şurada_; yürütülen farklı görevlerin özeti.

Sizinkinde küçük yeşil bir tik işareti varsa, tebrikler. İlk boru hattınızı çalıştırdınız. Harika.

Burada çalışan sürecin adı yazıyor: Say Hello. Bir kez çalıştığını ve başarılı olduğunu söylüyor. Bu, ilerledikçe güncellenir; dolayısıyla daha büyük bir boru hattı çalıştırırken ilerleme burada gösterilir. Ancak bu çok küçük olduğu için neredeyse anında çalışıyor.

## 1.2.2. Çıktıyı ve Günlükleri Work Dizininde Bulma

Bir Nextflow boru hattı çalıştırdığınızda, bu süreçlerin her biri birbirine bağlanır ve her süreç, daha önce söylediğim gibi, bir veya birden fazla görev oluşturabilir. Bu örnekte bu süreçten tek bir görev elde ettik. Yalnızca bir kez çalıştı ve bu görev _hash_'i altında tamamlandı.

Nextflow, çalışma dizininizdeki dosyalarla doğrudan ilgilenmez; work adlı özel bir klasör oluşturur. "ls" yaparsam burada belirdiğini göreceğiz: _work_. Bu klasörün içinde çalışan her görev için alt dizinler bulunur. Bu, söz konusu hash ile eşleşiyor. "ls work/c4" yaparsam ve ardından kısaltılmış ama 203 ile başlayan kısmı görürsem, bu boru hattını çalıştırdığımızda bu süreç tarafından oluşturulan çalışma dizinidir. Bunu yan tarafta da görebilirsiniz.

Bu dosyaları listelediğimde output.txt dosyasının oluşturulduğunu görebilirsiniz. Burada da görebilirsiniz. Ayrıca normal "ls" ile görünmeyen bir sürü gizli dosya var.

output.txt'ye tıklarsam, çıktımız orada. Harika. Boru hattı çalıştı.

Özünde tek satırlık bir Bash betiği çalıştırmak için oldukça fazla şablon kod gibi görünebilir; ancak süreçlerimiz daha karmaşık hâle geldikçe daha anlamlı olacak. Nextflow'daki bu work dizini ve oluşturulan bu dosyalar, Nextflow'u bu kadar güçlü kılan şeyin gerçek omurgasıdır.

Her görev, boru hattının her öğesi diğer tüm görevlerden yalıtılmıştır. Yeniden üretilebilir. Birbirleriyle çakışmıyorlar ve her şey paralel çalışabiliyor. Alıştığınızda gerçekten güzel bir yöntem; çünkü bu yalıtım sayesinde tek bir göreve girip tam olarak ne olduğunu görebilir ve hata ayıklayabilirsiniz.

Work dizinindeki diğer dosyalara hızlıca bakalım. Yukarıdan aşağıya, _.command.begin_ adlı bir dosya var. Bu boş. Nextflow'un "tamam, görevi başlatıyorum" dediği, sentinel dosya olarak adlandırılan bir dosya. Burada ilginç bir şey yok.

Ardından _.command.error_, _.command.log_ ve _.command.out_ var. Bunların hepsi çalışan Bash komutunun veya betiğin çıktılarıdır. Bu standart hata. Bu standart çıktı. Bu ise ikisinin çıktıkları sırayla birleşimi. Böylece mantıksal sırayı elde edersiniz.

Tamam, bunlar da bu örnek için boştu; pek ilginç değil. Ancak _.command.run_'a gelince işler daha ilginç hâle geliyor.

Bu genellikle çok uzun bir betiktir. Nextflow'un gerçekte yürüttüğü şey budur. Buraya girerseniz Nextflow'un tüm iç mantığını görmeye başlarsınız; ne yaptığını ve sürecinizi nasıl yürüttüğünü anlarsınız. Bu, nerede çalıştığınıza bağlıdır: yerel olarak mı çalışıyoruz yoksa SLURM'a bir iş olarak mı gönderiyoruz? İkinci durumda en üstte SLURM başlıkları olacaktır. Tüm bu farklı kurulumlar.

Genel olarak bu dosyaya hiçbir zaman bakmanız gerekmez. Nextflow tarafından otomatik olarak oluşturulur ve boru hattınıza özgü gerçekten özel bir şey içermez. Ama çalışan şeyin özü budur.

Bir sonraki dosya çok daha ilginç. _.command.sh_, sürecinizden gelen oluşturulmuş betiktir; burada Nextflow'un Bash başlığını eklediğini ve ardından betik bloğumuzda yer alan komutumuzu yürüttüğünü görebilirsiniz.

_.command.run_ dosyasının yaptığı tek şey bu _.command.sh_ dosyasını çalıştırmaktır.

Bu gerçekten kullanışlı bir dosyadır; bir şeyi hata ayıklamaya çalışırken ve Nextflow boru hattınızın mantığının beklediğiniz şeyi yapıp yapmadığını kontrol ederken en çok baktığınız dosya genellikle budur.

Son olarak _.exitcode_ adlı bir dosya var; bu, görevden çıkış kodunu yakalar. Bu örnekte başarılıydı, dolayısıyla çıkış kodu sıfırdı.

Bir şeyler ters giderse, bellek yetersizliği veya başka bir sorun nedeniyle başarısız olursa, ne yanlış gittiğini anlamak için bu dosya çok kullanışlıdır.

## 1.3. İş Akışını Tekrar Çalıştırma

Work dizinleri hakkında anlaşılması gereken bir şey daha var: Bu boru hattını tekrar tekrar çalıştırmaya devam edersem, yani _"nextflow run hello-world.nf"_ yaparsam, tam olarak aynı şeyi yapacak; ancak bu sefer yeni bir görev kimliği olacak. Bu hash'in farklı olduğunu görebilirsiniz ve şimdi work'e bakarsam iki hash dizini var. Bunlar yine birbirinden ayrı.

Dolayısıyla bir Nextflow iş akışını her çalıştırdığınızda, cache kullanan resume'u kullanmadığınız sürece (buna daha sonra değineceğiz), bu süreçleri birbirinden ayrı yeni work dizinlerinde yeniden çalıştıracaktır. Dosya adı çakışması yaşamazsınız, bu tür sorunlarla karşılaşmazsınız. Her şey yalıtılmış ve temiz.

Bu dizine girersek tüm aynı dosyaları ve sıfırdan yeniden oluşturulmuş aynı _output.txt_'yi görebiliriz.

## 2. Çıktıları Yayımlama

Tamam, bu Nextflow'un kendi içinde harika; boru hattınızı çalıştırırken her şeyin birbirinden ayrı, temiz ve yönetilebilir olmasını sağlıyor.

Ancak sonuçlarınızı incelemeye çalışan biri için pek kullanışlı değil. Sonuç dosyalarınızı bulmak için binlerce farklı work dizinini karıştırmak istemezsiniz. Zaten bunu yapmanız da gerekmiyor. Work dizinleri, dosyalarınızın oluşturulduğu nihai konum olarak tasarlanmamıştır.

Bunu dosyalarımızı yayımlayarak yapıyoruz.

## 2.1.1. sayHello Sürecinin Çıktısını Bildirme

Betiğimize geri dönersek, iş akışı bloğumuzda çalışacağız. Hangi dosyaları beklediğimizi, hangi dosyalarla ilgilendiğimizi söyleyeceğiz; ardından output bloğu adlı yeni bir blok oluşturacağız.

Bu, Nextflow'un 26.04 sürümünde varsayılan olarak gelen sözdizimi ayrıştırıcısıyla birlikte gelen yeni sözdizimidir. Dolayısıyla daha önce biraz Nextflow kullandıysanız, bu yeni olan şeylerden biridir.

Main bloğumuz var ve şimdi publish yazacağım; Nextflow'a yayımlama hakkında ne bekleyeceğini söyleyeceğim. Buna _first_output_ diyeceğiz ve _sayHello.out_ olarak adlandıracağız.

Orada yanlışlıkla bir yazım hatası yaptım; ancak bu, Nextflow VS Code uzantısının bazı özelliklerini de vurgulamak için iyi bir fırsat. Hemen altında kırmızı dalgalı bir çizgi belirdi ve bir şeylerin yanlış olduğunu söylüyor. Üzerine geldiğimde bu değişkenin tanımlı olmadığını söylüyor.

Bu örnekte oldukça açık; bir yazım hatası yaptım. sayHello yazmak istedim ve dalgalı çizgi kayboluyor.

Şimdi mor renkte. Nextflow sözdizimi ayrıştırıcısı bunun bir süreç olduğunu biliyor ve üzerine geldiğimde bu sürecin nasıl göründüğünün kısaltılmış bir temsilini veriyor. Böylece bir bakışta herhangi bir girdi almadığını ve bu çıktıyı verdiğini hızlıca görebiliyorum. VS Code'da bu uzantıyla çalışmak, kod yazarken size çok sayıda bağlamsal bilgi sağlıyor.

Bu sürecin çıktısına _.out_ sözdizimi ile başvurabileceğimizi unutmayın. Şu an buna istediğimiz adı verebiliriz; bu yalnızca rastgele bir değişken adıdır.

## 2.1.2. Betiğe Bir output: Bloğu Ekleme

Önemli hâle geldiği yer, buradaki yeni bloğumuzdur; bu artık workflow bloğunun altında, workflow'un içinde değil. Yine süslü parantezler. Burada Nextflow'a iş akışı tarafından oluşturulan tüm dosyaları nereye koyacağını söylüyoruz.

Şimdi burada oluşturduğum değişken adını alacağım, oraya koyacağım ve bunun için süslü parantezler ekleyeceğim. Nextflow'a bir path kullanmasını söyleyeceğim. Hata. Tırnak içinde path. Ve nokta kullanacağım. Bu yalnızca Nextflow'a dosyayı sonuçlar dizininin köküne koymasını söylüyor. Herhangi bir alt dizin veya benzeri bir şey yok.

İş akışımızı tekrar çalıştırmayı deneyelim. _"nextflow run hello-world.nf"_ yaparsam, temelde tamamen aynı görünmesi gerekiyor. Nextflow'da gerçekten hiçbir şey değişmedi. Aynı şeyleri çalıştırıyor. Sadece bunları yine work dizinlerinde yapıyor.

Ancak şimdi _"ls results/"_ yaparsam, burada results adında yeni bir dizin oluşturulduğunu göreceksiniz; bu, iş akışı yayımlaması için varsayılan temel dizindir. İçinde output.txt adlı bir dosya var.

_"ls -l results"_ yaparsam, bunun aslında work dizinine yumuşak bağlantılı olduğunu göreceksiniz. Yani bu gerçek bir dosya değil; work dizinine bağlantılı ve tüm dosyaları orada toplamış.

## 2.2. Özel Bir Konum Belirleme

"Results", bu path için varsayılan addır. İş akışını tekrar çalıştırırsam ve bu sefer tek kısa çizgi yaparsam (bu bir çekirdek Nextflow seçeneği olduğu için), _"-output-dir **my**results"_ yazarsam. Kısaca _"-o"_ da yapılabilir. Bu, dosyaların depolandığı temel dizin için farklı bir konum belirleyecek ve bir kez daha _myresults/_'da bir _output.txt_ olacak.

Bu harika; ancak muhtemelen tüm dosyaların kökte olmasını istemiyoruz. Biraz düzen istiyoruz; dolayısıyla burada istediğimiz adda bir alt dizin de oluşturabiliriz. Diyelim ki _"path 'hello_world'"_ yazalım ve bunu tekrar çalıştırayım: _"nextflow run hello-world.nf"_. Sonuçlar dizinine ve bir alt dizine gitmesi gerekiyor; gerçekten de şimdi en üstte results altında _hello_world/_ ve _output.txt_ var.

Dikkat edilmesi gereken önemli bir nokta: eski _output.txt_ dosyası hâlâ orada. Bunu yaptığınızda sonuçlar dizini silinmiyor. Yalnızca yeni dosyalar oraya kopyalanıyor. Aynı dosya adına sahip dosyaların üzerine yazılır; ancak eski dosyalar temizlenmez. Dolayısıyla boru hatlarını yeniden çalıştırırken biraz dikkatli olmanız gerekiyor. Zaten orada olan dosyaların üzerine yazmak istemiyorsanız boş bir dizin kullandığınızdan emin olun.

## 2.3. Yayımlama Modunu copy Olarak Ayarlama

Tamam, bu dosyaların yumuşak bağlantılar olduğundan bahsetmiştim. _"ls -l results/hello_world/"_ yaparsam, work dizinine yumuşak bağlantı verdiğini görebilirsiniz. HPC gibi bir ortamda çalışıyorsanız ve bunlar gerçekten büyük dosyalardır ve bunları kopyalamak istemiyorsanız bu genellikle iyi bir şeydir; çünkü dosyaların dosya sisteminde yalnızca bir kez depolandığı anlamına gelir.

Ancak work dizinini silerseniz: _"rm -r work"_ yapıp oluşturulan tüm ara dosyaları temizlersem, şimdi bu dosyayı okumaya çalışırsam (_"results/hello_world/"_), artık var olmayan bir dosyaya yumuşak bağlantı gösterecek ve veriler sonsuza dek kaybolacak; bu pek de iyi olmayabilir.

Bu nedenle genellikle yapabiliyorsanız yumuşak bağlantı yerine dosyaları kopyalamanın iyi bir uygulama olduğunu söylüyorum; çünkü daha güvenli. Yalnızca bu work dizinlerini silmediğiniz sürece iki kat daha fazla disk alanı kullanacağının farkında olun.

Bunu output bloğuyla yapmak için buradaki ilk çıktıya gideceğim. Daha önce path'i ayarladım ve şimdi mode'u ayarlayacağım; yazarken VS Code uzantısının burada bir çıktı yönergesi olduğunu bilerek öneriler sunduğunu görebilirsiniz. copy diyeceğim. Kaydediyorum.

İş akışını yeniden çalıştıralım. Dosyaları yeniden oluşturacak, yeni work dizini.

Şimdi _"ls -l results/hello_world/"_ yaparsam, bunun gerçek bir dosya olduğunu ve artık yumuşak bağlantı olmadığını görebilirsiniz; Nextflow bunu kopyaladı. Bunu bilmek güzel. Dolayısıyla path ve mode, kendinizi oldukça sık yazarken bulacağınız şeyler.

Elbette bu çok basit. Bunu ilerledikçe daha karmaşık ve güçlü hâle getireceğiz; bu şeyleri nasıl dinamik ve çok ayrıntılı olmayan bir şekilde yapacağınızı göreceksiniz.

## 2.4. Süreç Düzeyindeki publishDir Yönergeleri Hakkında Not

Şimdi, buna başlarken söylediğim gibi, bu oldukça yeni bir sözdizimi biçimidir. Bunu kaydederken yalnızca Nextflow'un en son sürümlerinde mevcuttur ve Workflow Outputs olarak adlandırılmaktadır.

Bunu kullanırsanız harika. Nextflow içindeki pek çok başka harika özelliğin kilidini açar; örneğin bu dosyaların oluşturulurken kökenini izlemeye yardımcı olan Nextflow Lineage gibi. Yakında 26.04'te varsayılan olacak. Gelecekte ise iş akışlarınızı yazmanın tek yolu bu olacak.

Ancak şu anda bu geçiş aşamasındayken, doğada publishDir adlı bir şey kullanan boru hatlarıyla karşılaşabilirsiniz; bu eski yöntemdir ve iş akışı ile çıktı düzeyinde değil, süreç düzeyinde tanımlanır.

Bu bildirim temelde aynı şeyi söylüyor. Sonuç dosyalarını results adlı bir dizine yayımla ve copy modunu kullan diyor. Sözdiziminin çok benzer olduğunu görebilirsiniz. Ancak artık yeni boru hatları yazarken, AI sonuçlarında, belgelerde veya diğer boru hatlarında göreseniz bile bu publishDir yönergesini kullanmamaya çalışın; çünkü bu eski yöntemdir.

2026'da hepimiz workflow outputs kullanıyor olmalıyız.

Bunların tümü belgelenmiştir; bunu yapıyorsanız ve daha önce Nextflow kullandıysanız, buradaki Nextflow belgelerine gidebilirsiniz: nextflow.io/docs/. Eğiticiler bölümüne kaydırırsam, _Migrating to Workflow Outputs_ adlı bir eğitici var.

Gerçekten iyi. Tüm sözdizimini, eski sözdizimle nasıl eşdeğer olduğunu, neden değiştirdiğimizi ele alıyor; bir zaman çizelgesi ve her şey var. Çok sayıda örnekle tüm farklı senaryoları inceliyor. Böylece mevcut Nextflow kodunu yeni sözdizimine kolayca dönüştürebilirsiniz.

## 3.1. sayHello Sürecini Değişken Girdi Bekleyecek Şekilde Değiştirme

Tamam, basit betiğimiz var; bir süreç çalıştırıyor, bir dosya oluşturuyor, Nextflow'a bunun bir çıktı olduğunu söylüyor ve ardından Nextflow'a bu dosyayı nereye kaydedeceğini söylüyoruz. Bu iyi bir başlangıç.

Ancak her şey sabit kodlanmış olmasaydı daha ilginç olurdu. O hâlde, bu sürecin bir iş akışı başlatırken çalışma zamanında kontrol edebileceğimiz değişken bir girdi alabileceğini Nextflow'a nasıl söyleyeceğimizi düşünelim.

Bunu gerçekleştirmek için birkaç farklı şey yapmamız gerekiyor.

Öncelikle, bu sürece bir girdi değişkeni kabul edebileceğini söylememiz gerekiyor; yeni bir bildirim bloğu olarak buraya _input_ yazıyoruz. Buna _"val greeting"_ diyeceğiz.

val kısmı, aşağıdaki path'in karşılığıdır. Nextflow'a bunun bir değişken olduğunu söylüyor; bu örnekte bir string gibi. Üzerine tekrar gelirseniz, uzantı bunun ne anlama geldiğini söylüyor.

Ardından Nextflow'a bununla ne yapacağını söyleyeceğiz. Bir değişken olduğunu söylemek yeterli değil. Betik içinde bu değişkeni nasıl kullanacağınızı belirtmeniz gerekiyor. Bu nedenle buradaki sabit kodlanmış string'i kaldıracağım ve bir değişken koyacağım.

Bunu süslü parantezler olmadan hızlıca yapacağım; bunun izin verilen ve eski stil bir yöntem olduğunu göstermek için. Ancak yeni sözdizimi ile bunu gerçekten süslü parantezlerin içine koymayı öneriyoruz; bu şekilde bunun Nextflow tarafından burada enterpolasyon yapıldığı çok açık oluyor.

Harika. Yani _"input greeting"_, _$\{greeting\}_ içine gidiyor. Son olarak, iş akışı düzeyinde Nextflow'a bu sürecin artık bir girdi aldığını söylememiz gerekiyor. Bunu yapmak için temelde ona bir değişken vereceğiz.

## 3.2. Kullanıcı Girdisini Yakalamak İçin Bir Komut Satırı Parametresi Ayarlama

Bunu yeniden Hello World gibi sabit kodlayabiliriz ve bu gayet iyi çalışır; ancak açıkçası bize herhangi bir avantaj sağlamaz. Bunu çalışma zamanında yapılandırabilmek istedik; dolayısıyla Nextflow'u başlatırken CLI'da yapabilmek istiyoruz.

Bunu yapmanın yolu, _params_ adlı özel bir Nextflow kavramıdır. Buna _params.input_ diyeceğiz.

Bu, bu girdi değişkenini CLI'da açığa çıkarır ve Nextflow'u başlatırken çift kısa çizgi kullandığımız yer burasıdır.

Buna istediğim adı verebilirim; _hello_, _greeting_ diyebilirim. Fark etmez. Orada ne yaparsam, bir boru hattı başlatırken CLI seçeneği olarak açığa çıkacaktır. Bu Nextflow'un gerçek bir sihir numarası; çünkü bu parametrelerle iş akışı betiğinizi çok hızlı bir şekilde oluşturabileceğiniz anlamına geliyor ve esasen boru hattınız için özel bir CLI oluşturuyorsunuz; bu da başlatırken farklı seçenekleri anında özelleştirmeyi gerçekten kolaylaştırıyor.

O hâlde deneyelim. Terminale geri dönelim. _"nextflow run"_ komutumuz burada. Şimdi daha önce gördüğümüz _"params.input"_ ile eşleşen _"--input"_ yapacağım. Sanırım belgelerde Fransızca. Geraldine Fransızca konuşmayı seviyor. Ben İsveç'te yaşadığım için İsveççe yapacağım. "_Hej Världen_" diyeceğim ve Enter'a basacağım.

Tek tırnak veya çift tırnak kullanılabilir; bu yalnızca Bash'in bunu nasıl yorumladığını etkiler.

Nextflow boru hattını tamamen aynı şekilde çalıştırıyor. Çalışma dizinini ve her şeyin aynı olduğunu görebilirsiniz. Ancak şimdi _"results/hello_world/output"_'a gidersem, bunun yerine güzel İsveçcemizi görebiliriz.

Yani CLI'dan bir parametreye dinamik olarak bir girdi ilettik. Bunu sürece girdi olarak ilettik ve süreç bunu yorumlayarak bir betik bloğuna koydu; bu da o betik sonucunun çıktısını dinamik olarak değiştirdi. Oldukça harika.

Burada çok az sözdizimi ile oldukça karmaşık bir mantık. Umarım bunun artık nasıl ölçeklendiğini görebiliyorsunuzdur. Boru hatlarımızın mantığını ve özelleştirilebilirliğini Nextflow betiğine gerçekten bu şekilde yerleştiriyoruz.

## 3.4. Komut Satırı Parametreleri İçin Varsayılan Değerleri Kullanma

Tamam, bu harika. Ancak şu anda sorun şu: Bu boru hattını her çalıştırdığımda, çalışması için --input yapmam gerekiyor.

Bu parametresiz çalıştırmayı denersem, Nextflow bu parametreye ihtiyaç duyduğunu ancak ayarlanmadığını söyleyen bir hata verecek; dolayısıyla ne yapacağını bilmiyordu.

Bu arada bu harika yeni bir şey. Geçmişte Nextflow boş bir string ile çalışırdı ve anlaması güç tüm türden garip hatalar alırdınız. Ancak yeni Nextflow sözdizimi ayrıştırıcısında biraz daha dikkatli davranıyor ve hemen söylüyor.

Dolayısıyla her zaman her seçeneği belirtmek istemeyiz. Makul varsayılanlar belirtmek iyi bir uygulamadır. Peki bunu betiğimizde nasıl yapıyoruz?

Bunu yazarken, _params.input_'i kullandığımız yere doğrudan koyduğumuzu fark edeceksiniz. Dolayısıyla açık çözüm bir varsayılan tanımlamaktır; bunu iş akışındaki özel bir params bloğunda betiğin en üstünde yapıyoruz. Bu iş akışı betiğinde.

Yine burada bazı yeni sözdizimi var; dikkat edin. Bu gerçekten harika şeyler. Burada beklenen parametrenin adı var.

Ardından bu iki nokta üst üste karakterinden sonra değişkenin türünü tanımlıyoruz. Bunu yapmak zorunda değilsiniz; boş bırakabilirsiniz; ancak gerçekten güzel. Nextflow'a bir string beklediğimizi ve buna göre davranmasını söylüyor.

Bunun yerine bir sayı istiyorsak, örneğin float yazabiliriz; bu kayan noktalı bir sayı istediğimizi söyler. Bunu bir string ile çalıştırmayı denersek, float olmayan bir string verdiğimiz için hata verecektir. String olarak geçirirsek, bunun bir string olduğunu bilir. Baştaki sıfırları olsa ve tamamen sayısal olsa bile, yine de gerçek bir string olarak geçirecektir.

Dolayısıyla bu tür güvenliği Nextflow'un çok yeni bir özelliğidir; ancak kodunuzu yazmayı ve çalıştırmayı daha güvenli hâle getirmek için gerçekten güçlüdür.

Ardından bir eşittir işareti ve burada varsayılan değer var. Nextflow başlangıçta Barselona'da yazılmıştır; dolayısıyla varsayılan olarak biraz İspanyolca olması uygun görünüyor: _"Holà mundo!"_.

Tamam, o betiği kaydediyorum, geri dönüyorum, betiği _--input_ olmadan tekrar çalıştırıyorum. Bu sefer çalışması ve _results_'ta yeni dosyamızı oluşturması gerekiyor. Bu dosyada şimdi _"Holà mundo!"_ yazıyor.

Ancak bu yalnızca bir varsayılan; dolayısıyla daha önce yaptığımız şeyi hâlâ yapamayacağımız anlamına gelmiyor. Geri dönüp eski betiğimi bulursam, _"Hej Världen"_, komut satırında _--input_ yaptığım için bu varsayılanın üzerine yazacak ve output.txt dosyasında yeniden kullanacaktır.

Yani betikteki bu yalnızca ayarladığım varsayılan değerdir.

İş akışımızı daha karmaşık hâle getirip daha fazla parametre ekledikçe, betiğin en üstündeki bu params bloğu hepsini tek bir yerde toplamaya başlayacak.

Betiğinizde oldukça güzel bir simetri elde ediyorsunuz; burada tüm iş akışı girdileriniz en üstte ve iş akışı çıktılarınız en altta yer alıyor. İş akışınızın dış dünyayla arayüzünün ne olduğu çok açık. Dolayısıyla yeni sözdizimi ile yeni bir boru hattını çok hızlı bir şekilde alıp nasıl kullanacağınızı anlayabilirsiniz.

Son bir harika şey daha. Bununla bir varsayılan değer ayarlamak zorunda değiliz. params input yapıp varsayılan değer ayarlamazsak, Nextflow'a bu parametrenin zorunlu olduğunu söyler; yine boru hattı onsuz çalışmayı reddeder; ancak null hakkında bir şey yerine daha kullanışlı bir hata mesajı verir.

Yani girdinin beklendiğini, zorunlu olduğunu ancak komut satırında belirtilmediğini söylüyor. Çok güzel.

Tamam, umarım artık Nextflow boru hattınızı değişken girdiler ve parametrelerle nasıl kuracağınız, varsayılanı nasıl ayarlayacağınız, türleri nasıl ayarlayacağınız (Boolean doğru/yanlış bayrağı veya integer ya da farklı türler olabilir), bunları iş akışınıza nasıl ileteceğiniz, nereye gittiği ve ardından sürecinize nasıl enterpolasyon yapıldığı konusunda netlik kazandı. Ayrıca Nextflow'u başlatırken bunları komut satırında nasıl özelleştireceğinizi de biliyorsunuz. Bu, basit Bash komutumuza kıyasla daha ilginç görünmeye başlıyor.

## 4. İş Akışı Yürütmelerini Yönetme

Tamam. Sırada ne var? Bu bölümün son kısmında, tüm farklı iş akışı yürütmelerini nasıl yöneteceğimiz hakkında biraz konuşacağız. Kenar çubuğuma ve Explorer'da work'ün altına bakarsanız, bir sürü farklı boru hattı çalıştırdığımı ve bu work dizinlerinin oldukça uzadığını, çok sayıda olduğunu göreceksiniz.

Diğer bir şey de, daha önce söylediğim gibi, bu boru hattını her yeniden çalıştırdığımda yeni bir work dizinleri seti oluşturuyor ve tüm süreçleri sıfırdan yeniden çalıştırıyor; bu iyi bir şey. Bu amaçlanan davranış. Yeniden üretilebilir ve her şeyi taze bir şekilde yeniden oluşturuyor. Ancak açıkçası, çok uzun süren süreçler çalıştırıyorsanız, boru hattı yarı yolda çöktüyse veya boru hattının sonunda bir şeyi değiştirdiyseniz, her zaman baştan başlamak zorunda kalmak can sıkıcı.

## 4.1. Bir İş Akışını -resume ile Yeniden Başlatma

Neyse ki Nextflow, daha önce neyin çalıştırıldığını ve neyin mevcut olduğunu bilmekte gerçekten iyidir; eski sonuçları yeniden kullanmak çok basittir. Komutun sonuna yeni bir bayrak ekliyoruz: _"-resume"_.

Şimdi, input'ta iki kısa çizgi olduğuna dikkat edin; çünkü bu bir parametredir. resume'da yalnızca bir kısa çizgi var; çünkü bu bir çekirdek Nextflow seçeneğidir.

Uzun süredir Nextflow kullansanız bile insanları sürekli şaşırtıyor. Dolayısıyla her zaman bir veya iki kısa çizgiyi hatırlayın. Çekirdek bir Nextflow seçeneği olup olmadığına bağlıdır.

Tamam, şimdi _-resume_ yapıyorum ve tam olarak aynı iş akışını tekrar çalıştırıyorum. Bu sefer bir temel farkla neredeyse tamamen aynı görünmesi gerekiyor.

Buradaki çıktıda sonuçların önbelleğe alındığını görebilirsiniz. Aslında bu görev hash'i önceki çalıştırmayla tamamen aynı ve o work dizinini bütünüyle yeniden kullandı. Girdiler, çıktılar ve betik değiştirilmemişti. Dolayısıyla o dosyayı oradan alıyor ve boru hattında aşağı akış adımları varsa, bunları boru hattının bir sonraki adımına iletirdi.

Yani boru hattının tamamını baştan sona çalıştırıyor; ancak yapabildiği her görev için önbelleğe alınmış sonuçları kullanıyor.

_-resume_ yaptığınızda, çalışma dizininizdeki son boru hattı çalıştırmasını devam ettirir; ne olursa olsun. Ancak aslında orada yaptığınız herhangi bir önceki çalıştırmadan devam edebilirsiniz. Şimdiye kadar oldukça fazla yaptık.

## 4.2. Geçmiş Yürütmelerin Günlüğünü İnceleme

Hepsine bakmak için _"nextflow run"_ yerine _"nextflow log"_ yapabiliriz; bu bize tüm bu farklı çalıştırmaları gösteren güzel bir çıktı verecek. Ekranımı biraz küçültmem gerekiyor ki görebilelim; tüm bu farklı çalıştırmalar, ne zaman yaptığımız, oturum kimliği, komut ve her şey.

Buraya bakabilir ve bunlardan herhangi birinin çalıştırma adını alıp ardından o belirli olanlardan birini devam ettirebiliriz. Geri gidip _hungry_ekeblad_ adlı olanı devam ettirebilirim. Bunu _resume_'dan sonra koyuyorum.

Bu arada merak ediyorsanız, tüm bu sıfatlar ve bilim insanı adları Nextflow kaynak kodundadır. Gidip bularak favori bilim insanınızı ekleyerek Nextflow'a ilk pull request'inizi yapmak için gerçekten iyi bir yol.

Her neyse, bunu yaptım ve geri gidip bu iş akışı çalıştırmasından önbelleğe alınmış sonuçlara baktı, hâlâ yeniden kullanabileceğini fark etti ve kullandı. Yani önbelleğe alınmış sonuçları tekrar aldım.

## 4.3. Eski Work Dizinlerini Silme

Harika. Peki bu work dizinlerini temizlemek istersem ne olur? Burada çok sayıda var. Çok sayıda dosya var. Belki son birkaç boru hattı çalıştırmasından devam etmek istediğimi kesinlikle biliyorum; ancak öncekilerle ilgilenmiyorum.

O zaman buradan birini seçebilir ve başka bir Nextflow komutu kullanabilirim: _"nextflow clean"_. _"nextflow clean"_ yapacağım, _"-before"_ yapacağım ve bu örnekte _reverent_pike_ olan belirli çalıştırma adını belirteceğim; ardından _"-n"_ yapacağım; bu Nextflow'a yalnızca kuru çalıştırma yapmasını söylüyor. Yani gerçekte hiçbir şey yapmadan neyi sileceğini söylüyor. Bu work dizinlerini kaldıracaktı.

Bu mantıklı görünüyor. O hâlde aynı komutu tekrar yapacağım; ancak _"-n"_ yerine gerçekten temizleme yapmak için _"-f"_ kullanacağım. Bu sefer tüm bu dizinleri gerçekten kaldırdı. Work dizinlerine bakarsam, artık çok daha hafif görünüyor. Harika.

İşte yerel work dizinlerinizin tümünü önbelleği tamamen yok etmeden oldukça güvenli bir şekilde temizlemenin yolu. Yani isterseniz hâlâ devam edebilirsiniz.

Bu bayrakların ne işe yaradığını unutursanız, her Nextflow komutu için _"nextflow help"_ ve ardından komutun adını yapabilirsiniz. _"nextflow help clean"_ yaparsam, tüm farklı seçenekleri görebilirsiniz: _-after, -before, -but_; bu temizleme davranışını yapılandırmanın tüm farklı yolları. Oldukça harika.

## Özet

Tamam, Hello Nextflow'un birinci bölümünün sonu bu. Kursun oldukça yoğun bir başlangıcı; ancak umarım artık bir Nextflow betiğinin nasıl göründüğü hakkında oldukça iyi bir anlayışa sahipsinizdir: farklı temel parçalar, süreçler, iş akışları, çıktılar ve parametreler. Bunları komut satırından temel geçersiz kılmalarla nasıl yapılandıracağınızı, dinamik bir betikle dinamik bir girdi bloğunu nasıl oluşturacağınızı ve tüm iş yükü yürütmelerinizi nasıl yöneteceğinizi biliyorsunuz: neleri çalıştırdığınızı görme, devam ettirme, temizleme. Çok şey var. Uzun bir yol kat ettiniz. Dolayısıyla bir mola vermek, kısa bir yürüyüş yapmak ve bir fincan çay içmek istiyorsanız, şimdi muhtemelen iyi bir zaman. Hak ettiniz.

Buradan itibaren temelde bu temel üzerine inşa ediyoruz. Bunu nasıl daha karmaşık, daha güçlü hâle getirebiliriz? Nasıl daha esnek hâle getirebiliriz? Analizimizi ölçekte yapmak istediğimiz şeyleri nasıl yapabiliriz?

## Sınav

Şimdi web sayfasında birinci bölüme, hello world'e kaydırırsanız küçük bir sınav göreceksiniz; bu, Nextflow eğitiminin bu sürümü için yaptığımız yeni bir şey. Geçip bu bölümde yaptığımız tüm materyali anladığınızı kontrol etmek için kendinizi sınayabilirsiniz.

Bu bize gönderilmiyor veya herhangi bir yere kaydedilmiyor; yalnızca tarayıcınızda saklanıyor. Dolayısıyla cevaplarınızın ne olduğunu bilmiyoruz; ancak hiçbir şeyi kaçırmadığınızı veya yanlış anlamadığınızı kontrol etmek için küçük bir öz değerlendirme. İstediğiniz kadar deneyebilirsiniz.

Benim gibi VS Code örneğinizde terminalde kalmak istiyorsanız, _quiz_ komutunu yazıp hangi bölümde olduğunuzu söyleyebilirsiniz. _"Hello World"_ yapıyoruz ve ardından web tarayıcısındakiyle tamamen aynı sınav sorularını, ancak yalnızca terminalinizde yapabilirsiniz.

Harika. Tamam. Umarım bundan keyif alırsınız. Biraz eğlenin ve Nextflow kanalları hakkında konuşmak için bir dakika içinde bir sonraki bölümde görüşürüz.
