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

O zaman gerçekten basit bir yerden başlayalım. Terminale bir şey yazdırmak için "echo" ile başlayalım. "Hello World". Enter'a basıyorum ve bu terminale gidiyor. Hello World. Umarım bu kursu izleyen herkes için sürpriz olmamıştır.

Tamam, bununla bir şeyler yapalım. Sadece terminale yazdırmak yerine bir dosyaya yazalım. Klavyemde yukarı ok tuşuna basıyorum; bu Bash geçmişini döngüsel olarak gösteriyor ve son komutumu veriyor. Sonuna küçük bir büyüktür işareti ekleyeceğim; bu, komutun çıktısını bir dosyaya yönlendiriyor. Dosyayı output.txt olarak adlandıracağım.

Komutu çalıştırmak için tekrar Enter'a basıyorum; bu sefer terminalde hiçbir şey görünmüyor, ancak sol tarafta output.txt adında yeni bir dosyanın belirdiğini görebiliyoruz.

Bunu terminalde cat gibi bir komutla görüntüleyebiliriz. cat output.txt yazıyorum ve gerçekten "Hello World" yazıyor. Ayrıca üzerine çift tıklayarak VS Code'daki kod düzenleyicide açabiliriz.

## 1.1. Kodu İnceleme

Pekâlâ. Basit olduğunu söylemiştim. Sırada ne var? Bu süreci alıp tekrar deneyelim; ancak bu sefer Nextflow içinde yapalım.

Dediğim gibi, bu kursun tüm farklı bölümleri bir betikle başlıyor ve bu bölümün betiği Hello World olarak adlandırılıyor. O yüzden Hello World'ü bulacağım. Tek tıklayınca önizleme gösteriyor; düzenleyicide açmak için çift tıklıyorum. Terminali de hızlıca kapatıyorum.

Bu çok basit bir betik; olabildiğince sade. Yalnızca 22 satır uzunluğunda ve temelde aynı şeyi yapıyor. Hatta bir kısmı tanıdık gelecektir. Az önce yazdığımız şey bu. Bash komutunun bir dosyaya yönlendirildiğini görebiliyoruz.

Tamam. Başka ne var? Ayrıca bu dosyada Nextflow'un temel kavramlarından bazılarını görmeye başlayabiliriz. Burada kırmızıyla bir process ve bir workflow var. Bunlar Nextflow'daki özel anahtar kelimeler ve özel terminoloji.

## 1.1.1. Process Tanımı

Bir iş akışı içindeki farklı süreçler, iş akışınızın farklı mantıksal birimlerini kapsar. Her süreç tek bir şey yapar.

Çalıştırdığımızda bir görev ya da birden fazla görev oluşturur; bunlar bir pipeline'ın gerçek yürütme adımlarıdır. Tüm süreçler daha sonra altta gördüğümüz bir workflow bloğu içinde düzenlenir; bu örnekte yalnızca o tek süreci çalıştırır.

Süreç adı bu anahtar kelimeyi takip eder ve bu ad neredeyse her şey olabilir. Sürecin içeriği ise bu süslü parantezlerin içindedir.

Bir süreç için gerçekten tek bir gereklilik vardır: bir tür script veya exec bloğu içermesi gerekir. Bu, üçlü tırnak işaretleri içindedir ve pipeline'ı çalıştırdığımızda çalışma dizinine yazılan Bash betiğidir; bilgisayarınızda veya sunucunuzda gerçekten çalışan şey budur.

Bu genellikle Bash'tir, ancak en üste farklı bir shebang ekleyebilirsiniz; Python betiği veya R betiği de olabilir. Fark etmez. Bu betikteki her şey yürütülecektir.

Bu sürece eklediğimiz bir diğer şey ise çıktı bildirimidir. Bu, Nextflow'a bu sürecin output.txt adlı bir çıktı dosyası beklediğini söyler. Bunun bir path olduğunu belirtir; dolayısıyla bir dosya gibi ele alınmalıdır. Örneğin bu val olsaydı, bir değişken veya değer gibi ele alınacaktı.

Bu dosyayı oluşturmadığını unutmayın. Gerçekte onu oluşturan aşağıdaki betiktir. Bu yalnızca Nextflow'a bu dosya adıyla bir çıktı dosyası beklemesini söyler.

## 1.1.2. Workflow Tanımı

Tamam. Altta bir workflow var ve yine bir bildirim var. Bu Main olarak adlandırılmış. Bu, bir betik bloğunun iş akışı eşdeğeridir diyebiliriz. İş akışının bir şey yaptığı kısımdır. Bu örnekte sayHello adlı süreci çağırıyoruz.

Elbette normalde pipeline'ınız bundan çok daha karmaşık görünecektir. Muhtemelen birden fazla süreciniz olacak ve aralarındaki veri akışını düzenlemek için kanallar kullanacaksınız. Bunu bu kursun sonraki bölümlerinde ele alacağız; ancak şimdilik bu yeterli. Bu geçerli bir pipeline ve çalışması gerekiyor.

VS Code'da buradaki DAG önizlemesine bile tıklayabilirim. DAG veya DAG, pipeline'daki veri akışı yapısının bir temsilidir ve yan tarafta bir mermaid diyagramı olarak render edildiğini görebiliriz. Bu örnekte çok basit: bir kutu var, bu iş akışı; bir de sayHello adlı süreç. Ancak ilerledikçe bu daha ilginç görünebilir.

## 1.2. İş Akışını Çalıştırma

Tamam, bu iş akışını çalıştırmayı deneyelim ve ne olduğunu görelim.

Altta terminali tekrar açacağım, çıktıyı temizleyeceğim ve Nextflow Run yazacağım. Ardından betik adını yazacağım: hello-world.nf. Enter'a basıyorum.

Tamam, üstte bazı standart bilgiler var; bunlar Nextflow'un çalıştığını, hangi sürümün çalıştığını, betik adını ve diğer bilgileri gösteriyor.

Gerçekten aradığımız önemli şey _burada_: yürütülen farklı görevlerin özeti.

Sizinkinde küçük yeşil bir tik işareti varsa, tebrikler. İlk pipeline'ınızı çalıştırdınız. Harika.

Burada çalışan sürecin adı yazıyor: Say Hello. Bir kez çalıştığını ve başarılı olduğunu söylüyor. Bu ilerledikçe güncelleniyor; dolayısıyla daha büyük bir pipeline çalıştırırken ilerlemeyi burada göreceksiniz. Ancak bu çok küçük olduğu için neredeyse anında çalışıyor.

## 1.2.2. Çıktıyı ve Günlükleri Work Dizininde Bulma

Bir Nextflow pipeline'ı çalıştırdığınızda, bu süreçlerin her biri birbirine bağlanır ve her süreç, daha önce söylediğim gibi, bir veya birden fazla görev oluşturabilir. Bu örnekte bu süreçten tek bir görev vardı. Yalnızca bir kez çalıştı ve bu görev _hash_'i altında tamamlandı.

Nextflow, çalışma dizininizdeki dosyalarla doğrudan ilgilenmez; work adlı özel bir klasör oluşturur. "ls" yaparsam burada belirdiğini göreceğiz: _work_. İçinde çalışan her görev için alt dizinler var. Bu hash ile eşleşiyor. "ls work/c4" yaparsam ve ardından kısaltılmış ama 203 ile başlayan kısmı görürsem, bu pipeline'ı çalıştırdığımızda bu süreç tarafından oluşturulan çalışma dizinidir. Bunu yan tarafta da görebilirsiniz.

Bu dosyaları listelediğimde output.txt dosyasının oluşturulduğunu görebilirsiniz. Burada da görebilirsiniz. Ayrıca normal "ls" ile görünmeyen bir sürü gizli dosya var.

output.txt'ye tıklarsam, çıktımız orada. Harika. Pipeline çalıştı.

Temelde tek satırlık bir Bash betiği çalıştırmak için oldukça fazla şablon kod gibi görünebilir; ancak süreçlerimiz daha karmaşık hale geldikçe daha mantıklı gelecek. Nextflow'daki bu work dizini ve oluşturulan bu dosyalar, Nextflow'u bu kadar güçlü kılan şeyin gerçek omurgasıdır.

Her görev, pipeline'ın her öğesi diğer tüm görevlerden izole edilmiştir. Tekrarlanabilir. Birbirleriyle çakışmıyorlar ve her şey paralel çalışabilir. Bu izolasyon sayesinde alıştığınızda gerçekten güzel bir yöntem; tek bir göreve girip tam olarak ne olduğunu görebilir ve hata ayıklayabilirsiniz.

Work dizinindeki bu diğer dosyalara hızlıca bakalım. Yukarıdan aşağıya, _.command.begin_ adlı bir dosya var. Bu boş. Nextflow'un "tamam, görevi başlatıyorum" dediği, sentinel dosya olarak adlandırılan bir dosya. Burada ilginç bir şey yok.

Ardından _.command.error_, _.command.log_ ve _.command.out_ var. Bunların hepsi çalışan Bash komutundan veya bu betikten gelen çıktılar. Bu standart hata. Bu standart çıktı. Bu ikisi çıktıkları sırayla birleştirilmiş hali. Böylece mantıksal sırayı görürsünüz.

Tamam, bunlar da bu örnekte boştu; pek ilginç değil. Ancak _.command.run_'a gelince işler daha ilginç hale geliyor.

Bu genellikle çok uzun bir betiktir. Nextflow'un gerçekte yürüttüğü şey budur. Buraya girerseniz Nextflow'un tüm iç mantığını görmeye başlarsınız; ne yaptığını ve sürecinizi nasıl yürüttüğünü anlarsınız. Bu, yerel olarak mı çalıştığınıza yoksa SLURM'a bir iş olarak mı gönderdiğinize bağlı olarak değişir; ikinci durumda üstte SLURM başlıkları olacaktır. Tüm bu farklı kurulumlar.

Genel olarak bu dosyaya hiçbir zaman bakmanız gerekmez. Nextflow tarafından otomatik olarak oluşturulur ve pipeline'ınıza özgü gerçekten özel bir şey içermez. Ama çalışan şeyin özü budur.

Bir sonraki dosya çok daha ilginç. _.command.sh_, sürecinizden gelen oluşturulmuş betiktir. Burada Nextflow'un Bash başlığını eklediğini ve ardından betik bloğumuzda bulunan komutumuzu yürüttüğünü görebilirsiniz.

_.command.run_ dosyasının yaptığı tek şey bu _.command.sh_ dosyasını çalıştırmaktır.

Bu gerçekten kullanışlı bir dosya; bir şeyi hata ayıklamaya çalışırken ve Nextflow pipeline'ınızın mantığının beklediğiniz şeyi yapıp yapmadığını kontrol ederken en çok baktığınız dosya genellikle budur.

Son olarak _.exitcode_ adlı bir dosya var; bu, görevden çıkış kodunu yakalar. Bu örnekte başarılıydı, dolayısıyla çıkış kodu sıfırdı.

Bir şeyler ters giderse, bellek yetersizliği veya başka bir sorun nedeniyle başarısız olursa, ne yanlış gittiğini anlamak için bu çok kullanışlıdır.

## 1.3. İş Akışını Tekrar Çalıştırma

Work dizinleri hakkında anlaşılması gereken bir şey daha var: bu pipeline'ı tekrar tekrar çalıştırmaya devam edersem, yani _"nextflow run hello-world.nf"_ yaparsam, tam olarak aynı şeyi yapacak; ancak bu sefer yeni bir görev kimliği olacak. Bu hash'in farklı olduğunu görebilirsiniz ve şimdi work'e bakarsam iki hash dizini var. Bunlar yine birbirinden ayrı.

Dolayısıyla bir Nextflow iş akışını her çalıştırdığınızda, cache kullanan resume'u kullanmadığınız sürece (buna daha sonra değineceğiz), bu süreçleri birbirinden ayrı yeni work dizinlerinde yeniden çalıştıracak. Dosya adı çakışması yaşamazsınız, bu tür sorunlarla karşılaşmazsınız. Her şey izole ve temiz.

Bu dizine girersek tüm aynı dosyaları ve sıfırdan yeniden oluşturulmuş aynı _output.txt_'yi görebilirsiniz.

## 2. Çıktıları Yayımlama

Tamam, bu Nextflow'un kendi içinde harika; pipeline'ınızı çalıştırırken her şeyin birbirinden ayrı, temiz ve yönetilebilir olmasını sağlıyor.

Ancak sonuçlarınızı incelemeye çalışan biri için pek kullanışlı değil. Sonuç dosyalarınızı bulmak için binlerce farklı work dizinini karıştırmak istemezsiniz. Zaten bunu yapmanız da gerekmiyor. Work dizinleri, dosyalarınızın oluşturulduğu nihai konum olarak tasarlanmamıştır.

Bunu dosyalarımızı yayımlayarak yapıyoruz.

## 2.1.1. sayHello Sürecinin Çıktısını Bildirme

Betiğimize geri dönersek, iş akışı bloğumuzda çalışacağız. Hangi dosyaları beklediğimizi, hangi dosyalarla ilgilendiğimizi söyleyeceğiz ve ardından altında output bloğu adlı yeni bir blok oluşturacağız.

Bu, Nextflow'un 26.04 sürümünde varsayılan olarak gelen sözdizimi ayrıştırıcısıyla birlikte gelen yeni sözdizimidir. Daha önce biraz Nextflow kullandıysanız, bu yeni olan şeylerden biridir.

main bloğumuz var ve şimdi publish diyeceğim ve Nextflow'a yayımlamadan ne bekleyeceğini söyleyeceğim. Buna _first_output_ diyeceğiz ve _sayHello.out_ olarak adlandıracağız.

Orada yanlışlıkla bir yazım hatası yaptım, ancak bu aynı zamanda Nextflow VS Code uzantısının bazı özelliklerini de vurgulamak için iyi bir fırsat. Hemen altında küçük kırmızı dalgalı bir çizgi belirdi; bir şeylerin yanlış olduğunu söylüyor. Üzerine gelince bu değişkenin tanımlı olmadığını söylüyor.

Bu örnekte oldukça açık; yazım hatası yaptım. sayHello yazmak istedim ve dalgalı çizgi kayboluyor.

Şimdi mor renkte. Nextflow sözdizimi ayrıştırıcısı bunun bir süreç olduğunu biliyor ve üzerine gelince bu sürecin nasıl göründüğünün kısaltılmış bir temsilini veriyor. Böylece bir bakışta hiçbir girdi almadığını ve bu çıktıyı verdiğini hızlıca görebiliyorum. VS Code'da bu uzantıyla çalışmak, kod yazarken size çok sayıda bağlamsal bilgi sağlıyor.

Bu sürecin çıktısına _.out_ sözdizimi ile başvurabileceğimizi unutmayın. Şu an buna istediğimiz adı verebiliriz; bu yalnızca rastgele bir değişken adıdır.

## 2.1.2. Betiğe Bir output: Bloğu Ekleme

Önemli hale geldiği yer, buradaki yeni bloğumuzu yaptığımızda; bu artık workflow bloğunun altında, workflow'un içinde değiliz. Yine süslü parantezler. Burada Nextflow'a iş akışı tarafından oluşturulan tüm dosyaları nereye koyacağını söylüyoruz.

Şimdi burada oluşturduğum bu değişken adını alacağım, oraya koyacağım ve bunun için süslü parantezler ekleyeceğim. Nextflow'a bir path kullanmasını söyleyeceğim. Hata. Tırnak içinde path. Ve nokta kullanacağım. Bu yalnızca Nextflow'a dosyayı sonuçlar dizininin köküne koymasını söylüyor. Herhangi bir alt dizin veya benzeri bir şey yok.

İş akışımızı tekrar çalıştırmayı deneyelim. _"nextflow run hello-world.nf"_ yaparsam, umarım temelde tamamen aynı görünecektir. Nextflow açısından gerçekten hiçbir şey değişmedi. Aynı şeyleri çalıştırıyor. Sadece bunları yine work dizinlerinde yapıyor.

Ancak şimdi _"ls results/"_ yaparsam, burada results adında yeni bir dizin oluşturulduğunu göreceksiniz; bu, iş akışı yayımlaması için varsayılan temel dizindir. İçinde output.txt adlı bir dosya var.

_"ls -l results"_ yaparsam, bunun aslında work dizinine yumuşak bağlantılı olduğunu göreceksiniz. Yani bu gerçek bir dosya değil; work dizinine bağlantılı ve tüm dosyaları orada toplamış.

## 2.2. Özel Bir Konum Belirleme

"Results" bu path için varsayılan addır. İş akışını tekrar çalıştırırsam ve bu sefer tek tire ile _"-output-dir **my**results"_ yaparsam (çünkü bu bir çekirdek Nextflow seçeneği), kısaca _"-o"_ da yapılabilir. Bu, dosyaların depolandığı yer için farklı bir temel dizin belirleyecek ve bir kez daha, _myresults/_'da artık bir _output.txt_ var.

Bu harika, ancak muhtemelen tüm dosyaların kökte olmasını istemiyoruz. Biraz organizasyon istiyoruz; bu yüzden burada istediğimiz adda bir alt dizin de oluşturabiliriz. Diyelim ki _"path 'hello_world'"_ ve bunu tekrar çalıştırıyorum. _"nextflow run hello-world.nf"_. Sonuçlar dizinine bir alt dizine gitmeli ve gerçekten de şimdi üstteki results altında _hello_world/_ ve _output.txt_ var.

Dikkat edilmesi gereken önemli bir şey: eski _output.txt_ dosyası hâlâ orada. Bunu yaptığınızda sonuçlar dizini silinmiyor. Sadece yeni dosyalar oraya kopyalanıyor. Aynı dosya adına sahip dosyaların üzerine yazacaklar, ancak eskileri temizlemeyecekler. Bu nedenle pipeline'ları yeniden çalıştırırken biraz dikkatli olmanız gerekiyor. Zaten orada olan dosyaların üzerine yazmak istemiyorsanız boş bir dizin kullandığınızdan emin olun.

## 2.3. Yayımlama Modunu Kopyalama Olarak Ayarlama

Tamam, bu dosyaların yumuşak bağlantılar olduğundan bahsetmiştim. _"ls -l results/hello_world/"_ yaparsam, work dizinine yumuşak bağlantı verdiğini görebilirsiniz. Bu, HPC gibi bir ortamda çalışıyorsanız ve bunlar gerçekten büyük dosyalarsa ve bunları kopyalamak istemiyorsanız genellikle iyi bir şeydir; çünkü dosyaların dosya sisteminde yalnızca bir kez depolandığı anlamına gelir.

Ancak work dizinini silerseniz: _"rm -r work"_ yapıp oluşturulan tüm ara dosyaları temizlersem. Şimdi bu dosyayı okumaya çalışırsam _"results/hello_world/"_. Artık var olmayan bir dosyaya yumuşak bağlantı gösterecek ve veriler sonsuza dek kaybolacak ve geri alınamayacak; bu pek iyi olmayabilir.

Bu nedenle genel olarak, yapabiliyorsanız yumuşak bağlantı yerine dosyaları kopyalamanın iyi bir uygulama olduğunu söylüyorum; çünkü daha güvenli. Sadece bu work dizinlerini silmediğiniz sürece iki kat daha fazla disk alanı kullanacağını unutmayın.

Bunu output bloğuyla yapmak için, buradaki ilk çıktıya gideceğim. Daha önce path'i ayarladım ve şimdi mode'u ayarlayacağım. Yazarken VS Code uzantısının öneriler sunduğunu görebilirsiniz; bunun bir output yönergesi olduğunu biliyor. copy diyeceğim. Kaydediyorum.

İş akışını yeniden çalıştıralım. Dosyaları yeniden oluşturacak, yeni work dizini.

Şimdi _"ls -l results/hello_world/"_ yaparsam, bunun gerçek bir dosya olduğunu ve artık yumuşak bağlantı olmadığını görebilirsiniz; Nextflow bunu kopyaladı. Bunu bilmek güzel. Dolayısıyla path ve mode, kendinizi oldukça sık yazarken bulacağınız şeyler.

Elbette bu çok basit. Bunu ilerledikçe daha karmaşık ve güçlü hale getireceğiz; bu şeyleri nasıl dinamik ve çok ayrıntılı olmayan bir şekilde yapacağınızı göreceksiniz.

## 2.4. Süreç Düzeyindeki publishDir Yönergeleri Hakkında Not

Şimdi, buna başlarken söylediğim gibi, bu oldukça yeni bir sözdizimi biçimi. Bunu kaydederken yalnızca Nextflow'un en son sürümlerinde mevcut ve Workflow Outputs olarak adlandırılıyor.

Bunu kullanırsanız harika. Nextflow içinde Nextflow Lineage gibi pek çok başka harika özelliğin kilidini açıyor; bu dosyalar oluşturulurken kökenlerini takip etmeye yardımcı oluyor. Yakında 26.04'te varsayılan olacak. Ve gelecekte belirli bir tarihte, iş akışlarınızı yazmanın tek yolu bu olacak.

Ancak şu anda bu geçiş aşamasındayken, kullandığınız gerçek dünya pipeline'larında publishDir adlı bir şey görebilirsiniz; bu eski yöntemdir ve bu, iş akışı ve çıktı düzeyinde değil, süreç düzeyinde tanımlanır.

Bu bildirim temelde aynı şeyi söylüyor. Sonuç dosyalarını results adlı bir dizine yayımla ve kopyalama modunu kullan diyor. Sözdiziminin çok benzer olduğunu görebilirsiniz. Ancak artık yeni pipeline'lar yazarken, yapay zeka sonuçlarında, belgelerde veya diğer pipeline'larda göreseniz bile bu publishDir yönergesini kullanmamaya çalışın; çünkü bu eski yöntemdir.

2026'da hepimiz workflow outputs kullanıyor olmalıyız.

Bunu yapıyorsanız ve daha önce Nextflow kullandıysanız, Nextflow belgelerine gidebilirsiniz: nextflow.io/docs/. Eğiticiler bölümüne kaydırırsam, _Migrating to Workflow Outputs_ adlı bir eğitici var.

Gerçekten iyi. Tüm sözdizimini, eski sözdizimle nasıl eşdeğer olduğunu, neden değiştirdiğimizi ele alıyor; bir zaman çizelgesi ve her şey var. Çok sayıda örnekle tüm farklı senaryoları ele alıyor. Böylece mevcut Nextflow kodunu yeni sözdizimine kolayca dönüştürebilirsiniz.

## 3.1. sayHello Sürecini Değişken Girdi Bekleyecek Şekilde Değiştirme

Tamam, basit betiğimiz var; bir süreç çalıştırıyor, bir dosya oluşturuyor, Nextflow'a bunun bir çıktı olduğunu söylüyor ve ardından Nextflow'a bu dosyayı nereye kaydedeceğini söylüyoruz. Bu iyi bir başlangıç.

Ancak her şey sabit kodlanmış olmasaydı daha ilginç olurdu. Şimdi, bu sürecin bir değişken girdi alabileceğini Nextflow'a nasıl söyleyeceğimizi düşünelim; bu, bir iş akışı başlattığımızda çalışma zamanında kontrol edebileceğimiz bir şeydir.

Bunu gerçekleştirmek için birkaç farklı şey yapmamız gerekiyor.

Öncelikle, bu sürece bir girdi değişkeni kabul edebileceğini söylememiz gerekiyor; yeni bir bildirim bloğu olarak buraya _input_ yazıyoruz. Buna _"val greeting"_ diyeceğiz.

val kısmı, aşağıdaki path'in eşdeğeridir. Nextflow'a bunun bir değişken olduğunu söylüyor; bu örnekte bir string gibi. Üzerine tekrar gelirseniz, uzantı bunun ne anlama geldiğini söylüyor.

Ardından Nextflow'a bununla ne yapacağını söyleyeceğiz. Bir değişken olduğunu söylemek yeterli değil. Betik içinde bu değişkeni nasıl kullanacağınızı belirtmeniz gerekiyor. Bu yüzden buradaki sabit kodlanmış string'i kaldıracağım ve bir değişken koyacağım.

Bunu süslü parantezler olmadan hızlıca yapacağım; bunun izin verilen ve eski stil bir yöntem olduğunu göstermek için. Ancak yeni sözdizimi ile bunu gerçekten süslü parantezlerin içine koymayı öneriyoruz; bu, Nextflow tarafından burada enterpolasyon yapıldığını çok net hale getiriyor.

Harika. Yani _"input greeting"_, _$\{greeting\}_ içine gidiyor. Son olarak, iş akışı düzeyinde Nextflow'a bu sürecin artık bir girdi aldığını söylememiz gerekiyor. Bunu yapmak için temelde ona bir değişken vereceğiz.

## 3.2. Kullanıcı Girdisini Yakalamak İçin Bir Komut Satırı Parametresi Ayarlama

Bunu yeniden Hello World gibi sabit kodlayabiliriz ve bu gayet iyi çalışır; ancak açıkçası bize hiçbir avantaj sağlamaz. Bunu çalışma zamanında yapılandırabilmek istedik; Nextflow'u başlattığınızda CLI'da yapabilmek istiyoruz.

Bunu yapmanın yolu, _params_ adlı özel bir Nextflow kavramıdır. Buna _params.input_ diyeceğiz.

Bu, bu girdi değişkenini CLI'da açığa çıkarır ve Nextflow'u başlatırken çift tire kullandığımız yer burasıdır.

Buna istediğim adı verebilirim; _hello_, _greeting_ diyebilirim. Fark etmez. Oraya ne yazarsam, bir pipeline başlatırken CLI seçeneği olarak açığa çıkacak. Bu Nextflow'un gerçek bir sihir numarası; çünkü bu parametrelerle iş akışı betiğinizi çok hızlı bir şekilde oluşturabileceğiniz anlamına geliyor. Temelde pipeline'ınız için özel bir CLI oluşturuyorsunuz; başlatırken farklı seçenekleri anında özelleştirmeyi gerçekten kolaylaştırıyor.

O zaman deneyelim. Terminale geri dönelim. _"nextflow run"_ komutumuz burada. Şimdi daha önce gördüğümüz _"params.input"_ ile eşleşen _"--input"_ yapacağım. Sanırım belgelerde Fransızca. Geraldine Fransızca konuşmayı seviyor. Ben İsveç'te yaşadığım için İsveççe yapacağım. "_Hej Världen_" diyeceğim ve Enter'a basacağım.

Tek tırnak veya çift tırnak kullanabilirsiniz; bu yalnızca Bash'in bunu nasıl yorumladığını etkiler.

Nextflow pipeline'ını tamamen aynı şekilde çalıştırıyor. Çalışma dizinini ve her şeyin aynı olduğunu görebilirsiniz. Ancak şimdi _"results/hello_world/output"_'a gidersem. Bunun yerine güzel İsveçcemizi görebiliyoruz.

Yani CLI'dan bir parametreye dinamik olarak bir girdi ilettik. Bunu sürece girdi olarak ilettik ve süreç bunu yorumlayarak bir betik bloğuna koydu; bu da o betik sonucunun çıktısını dinamik olarak değiştirdi. Oldukça harika.

Burada çok az sözdizimi ile oldukça karmaşık bir mantık. Umarım bunun artık nasıl ölçeklendiğini görebiliyorsunuzdur. Pipeline'larımızın mantığını ve özelleştirilebilirliğini Nextflow betiğine bu şekilde gerçekten yerleştiriyoruz.

## 3.4. Komut Satırı Parametreleri İçin Varsayılan Değerleri Kullanma

Tamam, bu harika. Ancak şu an sorun şu: bu pipeline'ı her çalıştırdığımda, çalışması için --input yapmam gerekiyor.

Bu parametre olmadan çalıştırmayı denersem, Nextflow bu parametreye ihtiyaç duyduğunu ancak ayarlanmadığını söyleyen bir hata verecek. Ne yapacağını bilmiyordu.

Bu arada bu harika yeni bir şey. Geçmişte Nextflow boş bir string ile çalışırdı ve anlaşılması güç tüm türlerde garip hatalar alırdınız. Ancak yeni Nextflow sözdizimi ayrıştırıcısında biraz daha dikkatli davranıyor ve hemen söylüyor.

Bu nedenle her zaman her seçeneği belirtmek istemeyiz. Makul varsayılanlar belirlemek iyi bir uygulamadır. Peki bunu betiğimizde nasıl yapıyoruz?

Bunu yazarken _params.input_'i doğrudan kullandığımız yere koyduğumuzu fark edeceksiniz. Dolayısıyla açık çözüm bir varsayılan tanımlamaktır ve bunu iş akışındaki özel bir params bloğunda betiğin en üstünde yapıyoruz. Bu iş akışı betiğinde.

Yine burada bazı yeni sözdizimi var; dikkat edin. Bu gerçekten harika şeyler. Burada beklenen parametrenin adı var.

Bu iki nokta üst üste karakterinden sonra değişkenin türünü tanımlıyoruz. Bunu yapmak zorunda değilsiniz; boş bırakabilirsiniz, ancak gerçekten güzel. Nextflow'a bir string beklediğimizi ve buna göre davranmasını söylüyor.

Bunun yerine bir sayı istiyorsak, örneğin float yazabiliriz; bu kayan noktalı bir sayı istediğimizi söyler. Bununla çalıştırmayı denersek, float olmayan bir string verirsek hata verecektir. String olarak geçirirsek, string olduğunu bilir. Baştaki sıfırları olsa ve tamamen sayısal olsa bile, yine de gerçek bir string olarak geçirecektir.

Bu tür güvenliği Nextflow'un çok yeni bir özelliği; ancak kodunuzu yazmayı ve çalıştırmayı daha güvenli hale getirmek için gerçekten güçlü.

Ardından eşittir sembolü ve ardından buradaki varsayılan değer var. Nextflow başlangıçta Barselona'da yazıldı; bu nedenle varsayılan olarak biraz İspanyolca olması uygun görünüyor: _"Holà mundo!"_.

Tamam, o betiği kaydediyorum, geri dönüyorum, betiği _--input_ olmadan tekrar çalıştırıyorum. Bu sefer çalışmalı ve _results_'ta yeni dosyamızı oluşturmalı. Bu dosyada artık _"Holà mundo!"_ yazıyor.

Ancak bu yalnızca bir varsayılan; dolayısıyla daha önce yaptığımız şeyi hâlâ yapamayacağımız anlamına gelmiyor. Geri dönüp eski betiğimi bulursam, _"Hej Världen"_, komut satırında _--input_ yaptığım için bu varsayılanın üzerine yazacak ve output.txt dosyasında yeniden kullanacak.

Yani betikteki bu yalnızca ayarladığım varsayılan değerdir.

İş akışımızı daha karmaşık hale getirip daha fazla parametre ekledikçe, betiğin en üstündeki bu params bloğu hepsini tek bir yerde toplamaya başlayacak.

Betiğinizde oldukça güzel bir simetri elde ediyorsunuz; burada tüm iş akışı girdileriniz ve altta iş akışı çıktılarınız var. İş akışınızın dış dünyayla arayüzünün ne olduğu çok net. Böylece yeni sözdizimi ile yeni bir pipeline'ı çok hızlı alıp nasıl kullanacağınızı anlayabilirsiniz.

Son bir harika şey daha. Bununla bir varsayılan değer ayarlamak zorunda değiliz. params input yapıp varsayılan değer ayarlamazsak, Nextflow'a bu parametrenin zorunlu olduğunu söylüyor; yine pipeline, bu olmadan çalışmayı reddedecek, ancak null hakkında bir şey yerine daha kullanışlı bir hata mesajı verecek.

Yani girdinin zorunlu olduğunu ancak komut satırında belirtilmediğini söylüyor. Çok güzel.

Tamam, umarım artık Nextflow pipeline'ınızı değişken girdiler ve parametrelerle nasıl kuracağınız, varsayılanı nasıl ayarlayacağınız, türleri nasıl ayarlayacağınız (Boolean doğru/yanlış bayrağı, integer veya farklı türler olabilir), bunları iş akışınıza nasıl ileteceğiniz, nereye gittiği ve ardından sürecinize nasıl enterpolasyon yapıldığı konusunda net bir fikriniz var. Ayrıca Nextflow'u başlatırken bunları komut satırında nasıl özelleştireceğinizi de biliyorsunuz. Bu, basit Bash komutumuza kıyasla daha ilginç görünmeye başlıyor.

## 4. İş Akışı Yürütmelerini Yönetme

Tamam. Sırada ne var? Bu bölümün son kısmında, tüm farklı iş akışı yürütmelerini nasıl yöneteceğimiz hakkında biraz konuşacağız. Kenar çubuğuma ve Explorer'da work'ün altına bakarsanız, bir sürü farklı pipeline çalıştırdığımı ve bu work dizinlerinin oldukça uzadığını, çok sayıda olduğunu göreceksiniz.

Diğer bir şey de, daha önce söylediğim gibi, bu pipeline'ı her yeniden çalıştırdığımda yeni bir work dizinleri seti oluşturuyor ve tüm süreçleri sıfırdan yeniden çalıştırıyor. Bu iyi bir şey. Amaçlanan davranış bu. Tekrarlanabilir ve her şeyi taze bir şekilde yeniden oluşturuyor. Ancak açıkçası, çok uzun süren süreçler çalıştırıyorsanız, yarı yolda çöktüyse veya pipeline'ın sonunda bir şeyi değiştirdiyseniz, her zaman pipeline'ın başından başlamak zorunda kalmak can sıkıcı.

## 4.1. Bir İş Akışını -resume ile Yeniden Başlatma

Neyse ki Nextflow, daha önce neyin çalıştırıldığını ve neyin mevcut olduğunu bilmekte gerçekten iyidir; eski sonuçları yeniden kullanmak çok basit. Komutun sonuna yeni bir bayrak ekliyoruz: _"-resume"_.

Şimdi, input'ta iki tire var çünkü bu bir parametredir. resume'da yalnızca bir tire var çünkü bu bir çekirdek Nextflow seçeneğidir.

Uzun süredir Nextflow kullansanız bile insanları sürekli şaşırtıyor. Bu nedenle her zaman bir veya iki tire olduğunu hatırlayın. Çekirdek bir Nextflow seçeneği mi olduğuna bağlı.

Tamam, şimdi _-resume_ yapıyorum ve tam olarak aynı iş akışını tekrar çalıştırıyorum. Bu sefer tek önemli farkla neredeyse tamamen aynı görünmeli.

Çıktıda sonuçların önbelleğe alındığını görebilirsiniz. Aslında bu görev hash'i önceki çalıştırmayla tamamen aynı ve o work dizinini bütünüyle yeniden kullandı. Girdiler, çıktılar ve betik değiştirilmemişti. Bu nedenle o dosyayı oradan alıyor ve süreçte aşağı akış adımları varsa, bunları pipeline'daki bir sonraki adıma iletirdi.

Yani pipeline'ı baştan sona çalıştırıyor, ancak yapabildiği her görev için önbelleğe alınmış sonuçları kullanıyor.

_-resume_ yaptığınızda, çalışma dizininizdeki son pipeline çalıştırmasını devam ettirir; ne olursa olsun. Ancak aslında orada yaptığınız herhangi bir önceki çalıştırmadan devam edebilirsiniz. Şimdiye kadar oldukça fazla yaptık.

## 4.2. Geçmiş Yürütmelerin Günlüğünü İnceleme

Hepsine bakmak için _"nextflow run"_ yerine _"nextflow log"_ yapabiliriz; bu bize tüm bu farklı çalıştırmaları gösteren güzel bir çıktı verecek. Ekranımı biraz küçültmem gerekiyor ki görebilelim; tüm bu farklı çalıştırmalar, ne zaman yaptığımız, oturum kimliği, komut ve her şey.

Buraya bakabilir ve bunlardan herhangi birinin çalıştırma adını alıp o belirli olanlardan birini devam ettirebiliriz. Geri gidip _hungry_ekeblad_ adlı olanı devam ettirebilirim. Bunu _resume_'dan sonra koyuyorum.

Bu arada merak ediyorsanız, tüm bu sıfatlar ve bilim insanı adları Nextflow kaynak kodunda. Gidip bularak ve favori bilim insanınızı ekleyerek Nextflow'a ilk pull request'inizi göndermek için gerçekten iyi bir yol.

Her neyse, bunu yaptım ve geri gidip bu iş akışı çalıştırmasından önbelleğe alınmış sonuçlara baktı, hâlâ yeniden kullanabileceğini fark etti ve kullandı. Yani önbelleğe alınmış sonuçları tekrar aldım.

## 4.3. Eski Work Dizinlerini Silme

Harika. Peki bu work dizinlerini temizlemek istersem? Burada çok sayıda var. Çok sayıda dosya var. Belki son birkaç pipeline çalıştırmasından devam etmek istediğimi kesinlikle biliyorum, ancak öncekilerle ilgilenmiyorum.

O zaman buradan birini seçebilir ve başka bir Nextflow komutu kullanabilirim: _"nextflow clean"_. _"nextflow clean"_ yapacağım, _"-before"_ yapacağım ve bu örnekte _reverent_pike_ olan belirli çalıştırma adını belirteceğim. _"-n"_ yapacağım; bu Nextflow'a yalnızca kuru çalıştırma yapmasını söylüyor. Yani gerçekten hiçbir şey yapmadan neyi sileceğini söylüyor. Bu work dizinlerini kaldıracak.

Mantıklı görünüyor. Aynı komutu tekrar yapacağım, ancak _"-n"_ yerine gerçekten temizleme yapmak için _"-f"_ kullanacağım. Bu sefer tüm bu dizinleri gerçekten kaldırdı. Work dizinlerine bakarsam, artık çok daha hafif görünüyor. Harika.

İşte yerel work dizinlerinizin tümünü önbelleği tamamen yok etmeden oldukça güvenli bir şekilde temizlemenin yolu. Yani isterseniz hâlâ devam edebilirsiniz.

Bu bayrakların ne olduğunu unutursanız, her Nextflow komutu için _"nextflow help"_ ve ardından komutun adını yapabilirsiniz. _"nextflow help clean"_ yaparsam, tüm farklı seçenekleri görebilirsiniz: _-after, -before, -but_; bu temizleme davranışını yapılandırmanın tüm farklı yolları. Oldukça harika.

## Özetle

Tamam, Hello Nextflow'un birinci bölümünün sonu bu. Kursa oldukça yoğun bir başlangıç, ancak umarım artık bir Nextflow betiğinin nasıl göründüğü hakkında oldukça iyi bir anlayışa sahipsinizdir: farklı temel parçalar, süreçler, iş akışları, çıktılar ve parametreler. Bunları komut satırından temel geçersiz kılmalarla nasıl yapılandıracağınızı, dinamik bir betikle dinamik bir girdi bloğunu nasıl oluşturacağınızı ve tüm iş yükü yürütmelerinizi nasıl yöneteceğinizi biliyorsunuz: daha önce ne çalıştırdığınızı görme, devam ettirme, temizleme. Çok şey var. Uzun bir yol kat ettiniz. Bu nedenle kısa bir mola verip etrafta biraz dolaşmak ve bir fincan çay içmek istiyorsanız, şimdi tam zamanı. Hak ettiniz.

Buradan itibaren temelde bu temel üzerine inşa ediyoruz. Bunu nasıl daha karmaşık, daha güçlü hale getirebiliriz? Nasıl daha esnek hale getirebiliriz? Analizimizi ölçekte yapmak istediğimiz şeyleri nasıl yapabiliriz?

## Quiz

Şimdi web sayfasında birinci bölüme, hello world'e kaydırırsanız küçük bir quiz göreceksiniz; bu, Nextflow eğitiminin bu sürümü için yaptığımız yeni bir şey. Üzerinden geçip bu bölümde yaptığımız tüm materyali anlayıp anlamadığınızı kontrol edebilirsiniz.

Bu bize gönderilmiyor veya herhangi bir yere kaydedilmiyor; yalnızca tarayıcınızda saklanıyor. Yani cevaplarınızı bilmiyoruz, ancak hiçbir şeyi kaçırıp kaçırmadığınızı veya yanlış anlayıp anlamadığınızı kontrol etmek için küçük bir öz değerlendirme. İstediğiniz kadar deneyebilirsiniz.

Benim gibi VS Code örneğinizdeki terminalde kalmak istiyorsanız, _quiz_ komutunu yazıp hangi bölümde olduğunuzu söyleyebilirsiniz. _"Hello World"_ yapıyoruz ve ardından web tarayıcısındaki quiz sorularının tamamen aynısını, ancak yalnızca terminalinizde yapabilirsiniz.

Harika. Tamam. Umarım keyif alırsınız. Biraz eğlenin ve kanallar hakkında konuşmak için bir dakika içinde bir sonraki bölümde görüşürüz.
