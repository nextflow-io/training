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

Umarım artık Codespaces'te ya da eşdeğer bir ortamda VS Code çalışır durumdadır ve Explorer'da Hello Nextflow klasörünüz tüm bu farklı dosyalarla birlikte çalışma alanınızda görünüyordur.

Terminalde Bash kullanarak çok basit şeylerle başlayacağız; ardından aynı şeyleri Nextflow içinde yapıp yapamayacağımıza bakacağız. Böylece sözdiziminin nasıl göründüğü hakkında bir fikir edineceksiniz.

## 0. Isınma

O zaman gerçekten basit bir yerden başlayalım. Terminale bir şey yazdırmak için "echo" ile başlayalım. "Hello World". Enter'a basıyorum ve bu terminale gidiyor. Hello World. Umarım bu kursu izleyen herkes için sürpriz olmamıştır.

Tamam, bununla bir şeyler yapalım. Sadece terminale yazdırmak yerine bir dosyaya yazalım. Klavyemde yukarı ok tuşuna basıyorum; bu Bash geçmişini döngüsel olarak gösteriyor ve son komutumu veriyor. Sonuna küçük bir büyüktür sembolü ekliyorum; bu komutun çıktısını bir dosyaya yönlendiriyor. Dosyayı output.txt olarak adlandıracağım.

Komutu çalıştırmak için tekrar Enter'a basıyorum. Bu sefer terminalde hiçbir şey görünmüyor, ancak sol tarafta output.txt adında yeni bir dosyanın belirdiğini görebiliyoruz.

Bunu terminalde cat gibi bir komutla görüntüleyebiliriz. cat output.txt yazıyorum ve gerçekten "Hello World" yazıyor. Ayrıca üzerine çift tıklayarak VS Code'daki kod düzenleyicide açabiliriz.

## 1.1. Kodu İnceleme

Pekâlâ. Basit olduğunu söylemiştim. Sırada ne var? Bu süreci alıp tekrar deneyelim; ama bu sefer Nextflow içinde yapalım.

Dediğim gibi, bu kursun tüm farklı bölümleri bir betikle başlıyor ve bu bölümün betiği Hello World olarak adlandırılmış. O yüzden Hello World'ü bulacağım. Tek tıkladığımda önizleme gösteriyor; çift tıklayarak düzenleyicide açıyorum. Terminali de hızlıca kapatıyorum.

Bu çok basit bir betik; olabildiğince sade. Yalnızca 22 satır uzunluğunda ve temelde aynı şeyi yapıyor. Hatta bir kısmı tanıdık gelecektir. Az önce yazdığımız şey bu; bir dosyaya yönlendiren Bash komutunu görebiliyoruz.

Tamam. Başka ne var? Bu dosyada aynı zamanda Nextflow'un temel kavramlarından bazılarını görmeye başlayabiliriz. Burada kırmızıyla bir süreç (process) ve bir iş akışı (workflow) var. Bunlar Nextflow'daki özel anahtar kelimeler ve özel terminoloji.

## 1.1.1. Süreç Tanımı

Bir iş akışı içindeki farklı süreçler, iş akışınızın farklı mantıksal birimlerini kapsar. Her süreç tek bir şey yapar.

Çalıştırdığımızda bir görev ya da birden fazla görev oluşturur; bunlar bir pipeline'ın gerçek yürütme adımlarıdır. Tüm süreçler daha sonra altta gördüğümüz bir iş akışı bloğu içinde düzenlenir; bu örnekte yalnızca o tek süreci çalıştırır.

Süreç adı bu anahtar kelimenin ardından gelir ve bu ad neredeyse her şey olabilir. Sürecin içeriği ise bu süslü parantezlerin içindedir.

Bir süreç için gerçekten tek bir gereklilik vardır: bir tür betik (script) ya da exec bloğu içermesi gerekir. Bu, üçlü tırnak işaretleri içindedir ve pipeline'ı çalıştırdığımızda çalışma dizinine yazılan Bash betiğidir; bilgisayarınızda ya da sunucunuzda gerçekten çalışan şey budur.

Bu genellikle Bash'tir, ancak en üste farklı bir shebang ekleyebilirsiniz; Python betiği ya da R betiği olabilir. Fark etmez. Bu betik içindeki her şey çalıştırılacaktır.

Bu sürece eklediğimiz bir şey daha var: çıktı bildirimi. Bu, Nextflow'a bu sürecin output.txt adında bir çıktı dosyası beklediğini söyler. Bunun bir path olduğunu belirtir; yani bir dosya gibi ele alınmalıdır. Örneğin bu val olsaydı, bir değişken ya da değer gibi ele alınacağını söylerdi.

Bunun bu dosyayı oluşturmadığını unutmayın. Aslında onu oluşturmaz. Bunu aşağıdaki betik yapar. Bu yalnızca Nextflow'a bu dosya adıyla bir çıktı dosyası beklemesini söyler.

## 1.1.2. İş Akışı Tanımı

Tamam. Altta bir iş akışımız var ve yine bir bildirim var. Bu Main olarak adlandırılmış. Bu, bir betik bloğunun iş akışı eşdeğeridir diyebiliriz. İş akışının bir şey yaptığı kısımdır. Bu örnekte sayHello adlı süreci çağırıyoruz.

Tabii ki normalde pipeline'ınız bundan çok daha karmaşık görünecektir. Muhtemelen birden fazla süreciniz olacak ve aralarındaki veri akışını düzenlemek için kanallar kullanacaksınız. Bunu bu kursun sonraki bölümlerinde ele alacağız; ancak şimdilik bu yeterli. Bu geçerli bir pipeline ve çalışması gerekiyor.

VS Code'da buradaki DAG önizlemesine bile tıklayabilirim. DAG ya da DAG, pipeline'daki veri akışı yapısının bir temsilidir ve yan tarafta mermaid diyagramı olarak render edildiğini görebiliriz. Bu örnekte çok basit: bir iş akışını temsil eden bir kutu ve sayHello adlı bir süreç var; ama ilerledikçe daha ilginç görünebilir.

## 1.2. İş Akışını Çalıştırma

Tamam, bu iş akışını çalıştırmayı deneyelim ve ne olduğunu görelim.

Altta terminali tekrar açacağım, çıktıyı temizleyeceğim ve Nextflow Run yazacağım. Ardından betik adını yazacağım: hello-world.nf. Enter'a basıyorum.

Tamam, üstte Nextflow'un çalıştığını, hangi sürümün çalıştığını, betik adını ve diğer bilgileri gösteren standart bilgiler var.

Gerçekten aradığımız önemli şey _şurada_: çalıştırılan farklı görevlerin özeti.

Sizinkinde küçük yeşil bir tik işareti varsa, tebrikler. İlk pipeline'ınızı çalıştırdınız. Harika.

Burada çalışan sürecin adı yazıyor: Say Hello. Bir kez çalıştığını ve başarılı olduğunu söylüyor. Bu ilerledikçe güncellenir; daha büyük bir pipeline çalıştırırken ilerlemeyi burada göreceksiniz. Ama bu çok küçük olduğu için neredeyse anında çalışıyor.

## 1.2.2. Work Dizininde Çıktıyı ve Günlükleri Bulma

Bir Nextflow pipeline'ı çalıştırdığınızda, bu süreçlerin her biri birbirine bağlanır ve her süreç, daha önce söylediğim gibi, bir ya da birden fazla görev oluşturabilir. Bu örnekte bu süreçten tek bir görev vardı. Yalnızca bir kez çalıştı ve bu görev _hash_'i altında gerçekleşti.

Nextflow, çalışma dizininizdeki dosyalarla doğrudan ilgilenmez; work adında özel bir klasör oluşturur. "ls" yaparsam burada belirdiğini göreceğiz: _work_. İçinde çalışan her görev için alt dizinler var. Bu hash ile eşleşiyor. "ls work/c4" yaparsam ve ardından kısaltılmış ama 203 ile başlayan kısmı görürsem, bu pipeline'ı çalıştırdığımızda bu süreç tarafından oluşturulan çalışma dizinidir. Bunu yan tarafta da görebilirsiniz.

Bu dosyaları listelediğimde output.txt dosyasının oluşturulduğunu görebilirsiniz. Bunu burada da görebilirsiniz. Ayrıca normal "ls" ile görünmeyen bir sürü gizli dosya var.

output.txt'ye tıkladığımda, gerçekten çıktımız var. Harika. Pipeline çalıştı.

Temelde tek satırlık bir Bash betiği çalıştırmak için oldukça fazla şablon kod gibi görünebilir; ancak süreçlerimiz daha karmaşık hale geldikçe daha anlamlı olacak. Nextflow'daki bu work dizini ve oluşturulan bu dosyalar, Nextflow'u bu kadar güçlü kılan şeyin gerçek omurgasıdır.

Her görev, pipeline'ın her öğesi diğer tüm görevlerden izole edilmiştir. Tekrarlanabilir. Birbirleriyle çakışmıyorlar ve her şey paralel çalışabilir. Bu izolasyon sayesinde alıştığınızda gerçekten güzel bir yöntem; tek bir göreve girip tam olarak ne olduğunu görebilir ve hata ayıklayabilirsiniz.

Work dizinindeki bu diğer dosyalara hızlıca bakalım. Yukarıdan aşağıya, _.command.begin_ adlı bir dosya var. Bu boş. Nextflow'un "tamam, görevi başlatıyorum" dediği, sentinel dosyası olarak adlandırılan bir dosya. Burada ilginç bir şey yok.

Ardından _.command.error_, _.command.log_ ve _.command.out_ var. Bunların hepsi çalışan Bash komutunun ya da betiğin çıktıları. Bu standart hata. Bu standart çıktı ve bu ikisi çıktıkları sırayla birleştirilmiş. Böylece mantıksal sırayı elde ediyorsunuz.

Tamam, bunlar da bu örnek için boştu; pek ilginç değil. Ama _.command.run_'a gelince işler daha ilginç hale geliyor.

Bu genellikle çok uzun bir betiktir. Nextflow'un gerçekte çalıştırdığı şey budur. Buraya girerseniz Nextflow'un tüm iç mantığını görmeye başlarsınız; ne yaptığını ve sürecinizi nasıl çalıştırdığını anlarsınız. Bu, yerel olarak mı çalıştırdığınıza yoksa SLURM'a bir iş olarak mı gönderdiğinize bağlıdır; ikinci durumda üstte SLURM başlıkları olacaktır. Tüm bu farklı kurulumlar.

Genel olarak bu dosyaya hiç bakmanıza gerek yoktur. Nextflow tarafından otomatik olarak oluşturulur ve pipeline'ınıza özgü gerçekten özel bir şey içermez. Ama çalışan şeyin özü budur.

Bir sonraki dosya çok daha ilginç. _.command.sh_, sürecinizden gelen oluşturulmuş betiktir ve burada Nextflow'un Bash başlığını eklediğini ve ardından betik bloğumuzdaki komutumuzun çalıştırıldığını görebilirsiniz.

_.command.run_ dosyasının yaptığı tek şey bu _.command.sh_ dosyasını çalıştırmaktır.

Bu gerçekten kullanışlı bir dosya; bir şeyi hata ayıklamaya çalışırken ve Nextflow pipeline'ınızın mantığının beklediğiniz şeyi yapıp yapmadığını kontrol ederken en çok baktığınız dosya genellikle budur.

Son olarak _.exitcode_ adlı bir dosya var; bu, görevden çıkış kodunu yakalar. Bu örnekte başarılıydı, dolayısıyla çıkış kodu sıfırdı.

Bir şeyler ters giderse, bellek yetersizliği ya da başka bir sorun nedeniyle başarısız olursa, ne yanlış gittiğini anlamak için bu dosya çok kullanışlıdır.

## 1.3. İş Akışını Tekrar Çalıştırma

Work dizinleri hakkında anlaşılması gereken bir şey daha var: bu pipeline'ı tekrar tekrar çalıştırmaya devam edersem, yani _"nextflow run hello-world.nf"_ yaparsam, tam olarak aynı şeyi yapacak; ancak bu sefer yeni bir görev kimliği olacak. Bu hash'in farklı olduğunu görebilirsiniz ve şimdi work'e bakarsam iki hash dizini var. Bunlar yine birbirinden ayrı.

Yani bir Nextflow iş akışını her çalıştırdığınızda, önbelleği kullanan resume'u kullanmadığınız sürece (buna daha sonra değineceğiz), bu süreçleri birbirinden ayrı yeni work dizinlerinde yeniden çalıştıracaktır. Dosya adı çakışması yaşamazsınız, bu tür sorunlarla karşılaşmazsınız. Her şey izole ve temiz.

Bu dizine girersek aynı dosyaların ve sıfırdan yeniden oluşturulmuş aynı _output.txt_'nin orada olduğunu görebilirsiniz.

## 2. Çıktıları Yayımlamak

Tamam, bu Nextflow'un kendi içinde harika; pipeline'ınızı çalıştırırken her şeyin birbirinden ayrı, temiz ve yönetilebilir olmasını sağlıyor.

Ancak sonuçlarınızı incelemeye çalışan biri için pek kullanışlı değil. Sonuç dosyalarınızı bulmak için binlerce farklı work dizinini karıştırmak istemezsiniz. Zaten bunu yapmanız da gerekmiyor. Work dizinleri, dosyalarınızın oluşturulduğu nihai konum olarak tasarlanmamıştır.

Bunu dosyalarımızı yayımlayarak yapıyoruz.

## 2.1.1. sayHello Sürecinin Çıktısını Bildirme

Betiğimize geri dönersek, iş akışı bloğumuzda çalışacağız. Hangi dosyaları beklediğimizi, hangi dosyalarla ilgilendiğimizi söyleyeceğiz ve ardından output bloğu adında yeni bir blok oluşturacağız.

Bu, Nextflow'un 26.04 sürümünde varsayılan olarak gelen sözdizimi ayrıştırıcısıyla birlikte gelen yeni sözdizimi. Dolayısıyla daha önce biraz Nextflow kullandıysanız, bu yeni olan şeylerden biri.

Main bloğumuz var ve şimdi publish yazacağım; Nextflow'a yayımlama hakkında ne bekleyeceğini söyleyeceğim. Buna _first_output_ diyeceğiz ve _sayHello.out_ olarak adlandıracağız.

Orada yanlışlıkla bir yazım hatası yaptım, ama bu aynı zamanda Nextflow VS Code uzantısının özelliklerini de vurgulamak için iyi bir fırsat. Hemen altında kırmızı dalgalı bir çizgi belirdi; bir şeylerin yanlış olduğunu söylüyor. Üzerine geldiğimde bu değişkenin tanımlı olmadığını söylüyor.

Bu durumda oldukça açık; yazım hatası yaptım. sayHello yazmak istedim ve dalgalı çizgi kayboluyor.

Şimdi mor renkte. Nextflow sözdizimi ayrıştırıcısı bunun bir süreç olduğunu biliyor ve üzerine geldiğimde bu sürecin nasıl göründüğünün kısaltılmış bir temsilini veriyor. Böylece bir bakışta hiçbir girdi almadığını ve bize bu çıktıyı verdiğini hızlıca görebiliyorum. VS Code'da bu uzantıyla çalışmak, kod yazarken size çok sayıda bağlamsal bilgi sağlıyor.

Bu sürecin çıktısına _.out_ sözdizimi ile başvurabileceğimizi unutmayın. Şu an buna istediğimiz her şeyi diyebiliriz; bu yalnızca keyfi bir değişken adı.

## 2.1.2. Betiğe Bir output: Bloğu Ekleme

Önemli hale geldiği yer, buradaki yeni bloğumuzu yaptığımızda; bu artık iş akışı bloğunun altında, artık workflow'un içinde değiliz. Yine süslü parantezler. Ve burada Nextflow'a iş akışı tarafından oluşturulan tüm dosyaları nereye koyacağını söylüyoruz.

Şimdi burada oluşturduğum değişken adını alacağım, oraya koyacağım ve bunun için süslü parantezler ekleyeceğim. Nextflow'a bir path kullanmasını söyleyeceğim. Hata. Tırnak içinde path. Ve nokta kullanacağım. Bu yalnızca Nextflow'a dosyayı sonuçlar dizininin köküne koymasını söylüyor. Herhangi bir alt dizin ya da başka bir şey değil.

İş akışımızı tekrar çalıştırmayı deneyelim. _"nextflow run hello-world.nf"_ yaparsam, umarım temelde tamamen aynı görünecektir. Nextflow açısından gerçekten hiçbir şey değişmedi. Aynı şeyleri çalıştırıyor. Sadece bunları yine work dizinlerinde yapıyor.

Ama şimdi _"ls results/"_ yaparsam, burada results adında yeni bir dizin oluşturulduğunu göreceksiniz; bu, iş akışı yayımlaması için varsayılan temel dizindir. İçinde output.txt adında bir dosya var.

_"ls -l results"_ yaparsam, bunun aslında work dizinine yumuşak bağlantılı (soft link) olduğunu göreceksiniz. Yani bu gerçek bir dosya değil; work dizinine bağlantılı ve tüm dosyaları orada toplamış.

## 2.2. Özel Bir Konum Belirleme

"Results" bu path için varsayılan addır. İş akışını tekrar çalıştırırsam ve bu sefer tek tire ile _"-output-dir **my**results"_ yaparsam (çünkü bu bir çekirdek Nextflow seçeneği), kısaca _"-o"_ da yapılabilir. Bu, dosyaların depolandığı yer için farklı bir temel dizin belirleyecek ve yine _myresults/_ altında bir _output.txt_ olacak.

Bu harika, ama muhtemelen tüm dosyaların kök dizinde olmasını istemiyoruz. Biraz organizasyon istiyoruz; bu yüzden istediğimiz adda bir alt dizin de oluşturabiliriz. Diyelim ki _"path 'hello_world'"_ ve bunu tekrar çalıştırıyorum. _"nextflow run hello-world.nf"_. Sonuçlar dizinine bir alt dizine gitmeli ve gerçekten de şimdi üstteki results altında _hello_world/_ ve _output.txt_ var.

Dikkat edilmesi gereken önemli bir şey: eski _output.txt_ dosyası hâlâ orada. Bunu yaptığınızda sonuçlar dizini silinmiyor. Sadece yeni dosyalar oraya kopyalanıyor. Aynı dosya adına sahip dosyaların üzerine yazacaklar, ancak eskileri temizlemeyecekler. Bu yüzden pipeline'ları yeniden çalıştırırken biraz dikkatli olmanız gerekiyor. Zaten orada olan dosyaların üzerine yazmak istemiyorsanız, boş bir dizin kullandığınızdan emin olun.

## 2.3. Yayımlama Modunu Kopyalama Olarak Ayarlama

Tamam, bu dosyaların yumuşak bağlantılar olduğundan bahsetmiştim. _"ls -l results/hello_world/"_ yaparsam, work dizinine yumuşak bağlantı verdiğini görebilirsiniz. HPC gibi bir ortamda çalışıyorsanız ve bunlar gerçekten büyük dosyalardır ve bunları kopyalamak istemiyorsanız bu genellikle iyi bir şeydir; çünkü dosyaların dosya sisteminde yalnızca bir kez depolandığı anlamına gelir.

Ancak work dizinini silerseniz: _"rm -r work"_ yapıp oluşturulan tüm ara dosyaları temizlersem. Şimdi bu dosyayı okumaya çalışırsam _"results/hello_world/"_. Artık var olmayan bir dosyaya yumuşak bağlantı gösterecek ve veriler sonsuza dek kaybolacak ve geri alınamayacak; bu pek de iyi olmayabilir.

Bu yüzden genel olarak, yapabiliyorsanız yumuşak bağlantı yerine dosyaları kopyalamanın iyi bir uygulama olduğunu söylüyorum; çünkü daha güvenli. Sadece bu work dizinlerini silmediğiniz sürece iki kat daha fazla disk alanı kullanacağının farkında olun.

Bunu output bloğuyla yapmak için, buradaki ilk çıktıya gideceğim. Daha önce path'i ayarladım ve şimdi mode'u ayarlayacağım; yazarken VS Code uzantısının öneriler sunduğunu görebilirsiniz, bunun bir output yönergesi olduğunu biliyor. copy diyeceğim. Kaydet'e basıyorum.

İş akışını yeniden çalıştıralım. Dosyaları yeniden oluşturacak, yeni work dizini.

Şimdi _"ls -l results/hello_world/"_ yaparsam, bunun gerçek bir dosya olduğunu ve artık yumuşak bağlantı olmadığını görebilirsiniz; Nextflow bunu kopyaladı. Bunu bilmek güzel. Yani path ve mode, kendinizi oldukça sık yazarken bulacağınız şeyler.

Tabii ki bu çok basit. Bunu ilerledikçe daha karmaşık ve güçlü hale getireceğiz; bu şeyleri nasıl dinamik ve çok ayrıntılı olmayan hale getireceğinizi göreceksiniz.

## 2.4. Süreç Düzeyindeki publishDir Yönergeleri Hakkında Not

Şimdi, buna başlarken söylediğim gibi, bu oldukça yeni bir sözdizimi biçimi. Bunu kaydederken yalnızca Nextflow'un en son sürümlerinde mevcut ve Workflow Outputs olarak adlandırılıyor.

Bunu kullanırsanız harika. Nextflow içinde Nextflow Lineage gibi pek çok başka harika özelliğin kilidini açıyor; bu dosyalar oluşturulurken kökenlerini takip etmeye yardımcı oluyor ve yakında 26.04'te varsayılan olacak. Gelecekte ilerleyen bir tarihte, iş akışlarınızı yazmanın tek yolu bu olacak.

Ancak şu an bu geçiş aşamasındayken, kullandığınız vahşi doğadaki pipeline'larda publishDir adlı bir şey görebilirsiniz; bu eski yöntemdir ve bu, iş akışı ve output düzeyinde değil, süreç düzeyinde tanımlanır.

Bu bildirim temelde aynı şeyi söylüyor. Sonuç dosyalarını results adlı bir dizine yayımla ve kopyalama modunu kullan diyor. Sözdiziminin çok benzer olduğunu görebilirsiniz. Ancak artık yeni pipeline'lar yazarken, AI sonuçlarında, belgelerde ya da diğer pipeline'larda göreseniz bile bu publishDir yönergesini kullanmamaya çalışın; çünkü bu eski yöntemdir.

2026'da hepimiz workflow outputs kullanıyor olmalıyız.

Bunu yapıyorsanız ve daha önce Nextflow kullandıysanız, buradan Nextflow belgelerine gidebilirsiniz: nextflow.io/docs/. Eğiticiler bölümüne kaydırırsam, _Migrating to Workflow Outputs_ adlı bir eğitici var.

Gerçekten iyi. Tüm sözdizimini, eski sözdizimle nasıl eşdeğer olduğunu, neden değiştirdiğimizi ele alıyor; bir zaman çizelgesi ve her şey var. Çok sayıda örnekle tüm farklı senaryoları ele alıyor. Böylece mevcut Nextflow kodunu yeni sözdizimine kolayca dönüştürebilirsiniz.

## 3.1. sayHello Sürecini Değişken Girdi Bekleyecek Şekilde Değiştirme

Tamam, basit betiğimiz var; bir süreç çalıştırıyor, bir dosya oluşturuyor, Nextflow'a bunun bir çıktı olduğunu söylüyor ve ardından Nextflow'a bu dosyayı nereye kaydedeceğini söylüyoruz. Bu iyi bir başlangıç.

Ama her şey sabit kodlanmış olmasaydı daha ilginç olurdu. Öyleyse şimdi, Nextflow'a bu sürecin bir iş akışı başlattığımızda çalışma zamanında kontrol edebileceğimiz değişken bir girdi alabileceğini nasıl söyleyeceğimizi düşünelim.

Bunu gerçekleştirmek için birkaç farklı şey yapmamız gerekiyor.

Öncelikle, bu sürece bir girdi değişkeni kabul edebileceğini söylememiz gerekiyor; yeni bir bildirim bloğu olarak buraya _input_ yazıyoruz. Buna _"val greeting"_ diyeceğiz.

val kısmı, aşağıdaki path'in eşdeğeridir. Nextflow'a bunun bir değişken olduğunu söylüyor; bu durumda bir string gibi. Üzerine tekrar gelirseniz, uzantı bunun ne anlama geldiğini söylüyor.

Ardından Nextflow'a bununla ne yapacağını söyleyeceğiz. Bir değişken olduğunu söylemek yeterli değil. Betikte bu değişkeni nasıl kullanacağınızı söylemeniz gerekiyor. Bu yüzden buradaki sabit kodlanmış string'i kaldıracağım ve bir değişken koyacağım.

Bunu süslü parantezler olmadan hızlıca yapacağım; buna izin verildiğini ve bunun eski stil yöntem olduğunu göstermek için. Ama şimdi yeni sözdizimi ile, bunu gerçekten süslü parantezlerin içine koymayı öneriyoruz; bu, bunun burada Nextflow tarafından enterpolasyon yapıldığını çok açık hale getiriyor.

Harika. Yani _"input greeting"_, _$\{greeting\}_ içine gidiyor. Son olarak, iş akışı düzeyinde Nextflow'a bu sürecin artık bir girdi aldığını söylememiz gerekiyor. Bunu yapmak için temelde ona bir değişken vereceğiz.

## 3.2. Kullanıcı Girdisini Yakalamak İçin Bir Komut Satırı Parametresi Ayarlama

Bunu yine Hello World gibi sabit kodlayabiliriz ve bu gayet iyi çalışır; ancak açıkça bize herhangi bir avantaj sağlamaz. Bunu çalışma zamanında yapılandırabilmek istedik; yani Nextflow'u başlattığınızda CLI'da yapabilmek istiyoruz.

Bunu yapmanın yolu, _params_ adlı özel bir Nextflow kavramıdır. Buna _params.input_ diyeceğiz.

Bunun yaptığı şey, bu girdi değişkenini CLI'da açığa çıkarmaktır ve Nextflow'u başlatırken çift tire kullandığımız yer burasıdır.

Buna istediğim her şeyi diyebilirim; _hello, greeting_ diyebilirim. Fark etmez. Orada ne yaparsam, bir pipeline başlatırken CLI seçeneği olarak açığa çıkacak. Bu Nextflow'un gerçek bir sihir numarası; çünkü bu parametrelerle iş akışı betiğinizi çok hızlı oluşturabileceğiniz anlamına geliyor ve esasen pipeline'ınız için özel bir CLI oluşturuyorsunuz; başlatırken farklı seçenekleri anında özelleştirmeyi gerçekten kolaylaştırıyor.

Öyleyse deneyelim. Terminale geri dönelim. _"nextflow run"_ komutumuz burada. Şimdi daha önce gördüğümüz _"params.input"_ ile eşleşen _"--input"_ yapacağım. Sanırım belgelerde Fransızca. Geraldine Fransızca konuşmayı seviyor. Ben İsveç'te yaşadığım için İsveççe yapacağım. "_Hej Världen_" diyeceğim ve Enter'a basacağım.

Tek tırnak ya da çift tırnak kullanabilirsiniz; bu yalnızca Bash'in bunu nasıl yorumladığını etkiler.

Nextflow pipeline'ını tamamen aynı şekilde çalıştırıyor. Çalışma dizinini ve her şeyin aynı olduğunu görebilirsiniz. Ama şimdi _"results/hello_world/output"_'a gidersem. Bunun yerine güzel İsveçcemizi görebiliyoruz.

Yani CLI'dan bir parametreye dinamik olarak bir girdi ilettik. Bunu sürece girdi olarak ilettik ve süreç bunu yorumlayıp bir betik bloğuna koydu; bu da o betik sonucunun çıktısını dinamik olarak değiştirdi. Oldukça harika.

Burada çok az sözdizimi ile oldukça karmaşık mantık. Ve umarım artık bunun nasıl ölçeklendiğini görmeye başlıyorsunuzdur. Pipeline'larımızın mantığını ve özelleştirilebilirliğini Nextflow betiğine bu şekilde gerçekten inşa ediyoruz.

## 3.4. Komut Satırı Parametreleri İçin Varsayılan Değerleri Kullanma

Tamam, bu harika. Ancak şu an sorun şu: bu pipeline'ı her çalıştırdığımda, çalışması için --input yapmam gerekiyor.

Bu parametresiz çalıştırmayı denersem, Nextflow bu parametreye ihtiyaç duyduğunu ve ayarlanmadığını söyleyen bir hata verecek. Ve ne yapacağını bilmiyordu.

Bu arada harika yeni bir şey. Geçmişte Nextflow boş bir string ile çalışırdı ve anlaşılması güç tüm türden garip hatalar alırdınız. Ama yeni Nextflow sözdizimi ayrıştırıcısında biraz daha dikkatli ve hemen söylüyor.

Yani her zaman her seçeneği belirtmek istemiyoruz. Makul varsayılanlar belirtmek iyi bir uygulamadır. Peki bunu betiğimizde nasıl yapıyoruz?

Bunu yazdığımızda, _params.input_'i doğrudan kullandığımız yere koyduğumuzu fark edeceksiniz. Yani açık çözüm bir varsayılan tanımlamaktır ve bunu iş akışındaki özel bir params bloğunda betiğin en üstünde yapıyoruz. Bu buradaki iş akışı betiğinde.

Yine burada bazı yeni sözdizimi var; dikkat edin. Bu gerçekten harika şeyler. Burada beklenen parametrenin adı var.

Ve bu iki nokta üst üste karakterinden sonra değişkenin türünü tanımlıyoruz. Bunu yapmak zorunda değilsiniz, boş bırakabilirsiniz; ama gerçekten güzel. Nextflow'a bir string beklediğimizi ve bunu böyle ele almasını söylüyor.

Bunun yerine bir sayı istiyorsak, örneğin float yazabiliriz; bu kayan noktalı bir sayı istediğimizi söyler. Bununla çalıştırmayı denersek, float olmayan bir string verirsek hata verecektir. Ve bunu böyle iletecektir. String yaparsak, bunun bir string olduğunu biliyor. Ve baştaki sıfırları olsa ve tamamen sayısal olsa bile, yine de gerçek bir string olarak iletecektir.

Bu tür güvenliği Nextflow'un çok yeni bir özelliği; ama kodunuzu yazmayı ve çalıştırmayı daha güvenli hale getirmek için gerçekten güçlü.

Ardından eşittir sembolü ve ardından buradaki varsayılan değer var. Nextflow başlangıçta Barselona'da yazıldı; bu yüzden burada biraz İspanyolca olması uygun görünüyor: varsayılan olarak _"Holà mundo!"_.

Tamam, o betiği kaydedeceğim, geri döneceğim, betiği _--input_ olmadan tekrar çalıştıracağım. Bu sefer çalışmalı ve _results_'ta yeni dosyamızı oluşturacak. Bu dosyada şimdi _"Holà mundo!"_ yazıyor.

Ama bu yalnızca bir varsayılan; yani daha önce yaptığımız şeyi hâlâ yapamayacağımız anlamına gelmiyor. Geri dönüp eski betiğimi bulursam, _"Hej Världen"_, komut satırında _--input_ yaptığım için bu varsayılanın üzerine yazacak ve output.txt dosyasında yine onu kullanacak.

Yani betikteki bu yalnızca ayarladığım varsayılan değer.

İş akışımızı daha karmaşık hale getirip daha fazla parametre ekledikçe, betiğin en üstündeki bu params bloğu hepsini tek bir yerde toplamaya başlayacak.

Ve betiğinizde oldukça güzel bir simetri elde ediyorsunuz; burada tüm iş akışı girdileriniz ve altta iş akışı çıktılarınız var. Ve iş akışınızın dış dünyayla arayüzünün ne olduğu çok açık. Böylece yeni sözdizimi ile yeni bir pipeline'ı çok hızlı alıp nasıl kullanacağınızı anlayabilirsiniz.

Son harika bir şey. Bununla bir varsayılan değer ayarlamak zorunda değiliz. params input yaparsak ama varsayılan değer ayarlamazsak, Nextflow'a bu parametrenin zorunlu olduğunu söylüyor ve yine pipeline onsuz çalışmayı reddedecek; ama null hakkında bir şey yerine daha kullanışlı bir hata mesajı verecek.

Yani girdisinin zorunlu olduğunu beklediğimizi söylüyor, ama komut satırında belirtilmedi. Çok güzel.

Tamam, umarım artık Nextflow pipeline'ınızı değişken girdiler ve parametrelerle nasıl kuracağınız, varsayılanı nasıl ayarlayacağınız, türleri nasıl ayarlayacağınız (Boolean doğru/yanlış bayrağı ya da integer ya da farklı türler olabilir), bunları iş akışınıza nasıl ileteceğiniz, nereden geçtiği ve ardından sürecinize nasıl enterpolasyon yapıldığı konusunda net bir fikriniz var. Ayrıca Nextflow'u başlatırken bunları komut satırında nasıl özelleştireceğinizi de biliyorsunuz. Bu, basit Bash komutumuzdan daha ilginç görünmeye başlıyor.

## 4. İş Akışı Çalıştırmalarını Yönetme

Tamam. Sırada ne var? Bu bölümün son kısmında, tüm farklı iş akışı çalıştırmalarını nasıl yöneteceğimiz hakkında biraz konuşacağız. Kenar çubuğuma ve Explorer'da work'ün altına bakarsanız, bir sürü farklı pipeline çalıştırdığımı ve bu work dizinlerinin oldukça uzadığını, çok sayıda olduğunu göreceksiniz.

Öte yandan, daha önce söylediğim gibi, bu pipeline'ı her yeniden çalıştırdığımda yeni bir work dizinleri seti oluşturuyor ve tüm süreçleri sıfırdan yeniden çalıştırıyor; bu iyi bir şey. Bu amaçlanan davranış. Tekrarlanabilir ve her şeyi taze bir şekilde yeniden oluşturuyor. Ama açıkça, çok uzun süren süreçler çalıştırıyorsanız, yarı yolda çöktüyse ya da pipeline'ın sonunda bir şeyi değiştirdiyseniz, her zaman pipeline'ın başından başlamak zorunda kalmak can sıkıcı.

## 4.1. Bir İş Akışını -resume ile Yeniden Başlatma

Neyse ki Nextflow, daha önce neyin çalıştırıldığını ve neyin mevcut olduğunu bilmekte gerçekten iyidir; eski sonuçları yeniden kullanmak çok basit. Komutun sonuna yeni bir bayrak ekliyoruz: _"-resume"_.

Şimdi, input'ta iki tire olduğuna dikkat edin; çünkü bu parametredir. resume'da yalnızca bir tire var; çünkü bu bir çekirdek Nextflow seçeneği.

Uzun süredir Nextflow kullansanız bile insanları sürekli şaşırtıyor. Bu yüzden her zaman bir ya da iki tire olduğunu hatırlayın. Çekirdek Nextflow seçeneği mi olduğuna bağlı.

Tamam, şimdi _-resume_ yapıyorum ve tam olarak aynı iş akışını tekrar çalıştırıyorum. Bu sefer bir temel farkla neredeyse tamamen aynı görünmeli.

Buradaki çıktıda sonuçların önbellekten alındığını görebilirsiniz. Aslında bu görev hash'i önceki çalıştırmayla tamamen aynı ve o work dizinini bütünüyle yeniden kullandı. Girdiler, çıktılar ve betik değiştirilmemişti. Bu yüzden sadece o dosyayı oradan alıyor ve süreçte aşağı akış adımları varsa, bunları pipeline'ın bir sonraki adımına iletirdi.

Yani pipeline'ı baştan sona hâlâ çalıştırıyor; ama yapabildiği her görev için önbelleğe alınmış sonuçları kullanıyor.

_-resume_ yaptığınızda, çalışma dizininizdeki son pipeline çalıştırmasından devam eder; ne olursa olsun. Ama aslında orada yaptığınız herhangi bir önceki çalıştırmadan devam edebilirsiniz. Ve şimdiye kadar oldukça fazla yaptık.

## 4.2. Geçmiş Çalıştırmaların Günlüğünü İnceleme

Hepsine bakmak için _"nextflow run"_ yerine _"nextflow log"_ yapabiliriz; bu bize tüm bu farklı çalıştırmaları gösteren güzel bir çıktı verecek. Ekranımı biraz küçültmem gerekiyor ki görebilelim; tüm bu farklı çalıştırmalar, ne zaman yaptığımız, oturum kimliği, komut ve her şey.

Buraya bakabilir ve bunlardan herhangi birinin çalıştırma adını alıp ardından o belirli olanlardan birini devam ettirebiliriz. Geri gidip _hungry_ekeblad_ adlı olanı devam ettirebilirim. Bunu _resume_'dan sonra koyuyorum.

Bu arada merak ediyorsanız, tüm bu sıfatlar ve bilim insanı adları Nextflow kaynak kodunda. Gidip bularak ve favori bilim insanınızı ekleyerek Nextflow'a ilk pull request'inizi yapmak için gerçekten iyi bir yol.

Her neyse, bunu yaptım ve geri gidip bu iş akışı çalıştırmasından önbelleğe alınmış sonuçlara baktı, hâlâ yeniden kullanabileceğini fark etti ve kullandı. Yani önbelleğe alınmış sonuçları tekrar aldım.

## 4.3. Eski Work Dizinlerini Silme

Harika. Peki bu work dizinlerini temizlemek istersem? Burada çok sayıda var. Çok sayıda dosya var. Belki son birkaç pipeline çalıştırmasından devam etmek istediğimi kesinlikle biliyorum, ama öncekilerle ilgilenmiyorum.

O zaman buradan birini seçebilir ve başka bir Nextflow komutu kullanabilirim: _"nextflow clean"_. _"nextflow clean"_ yapacağım, _"-before"_ yapacağım ve bu durumda _reverent_pike_ olan belirli çalıştırma adını, ve _"-n"_ yapacağım; bu Nextflow'a yalnızca kuru çalıştırma yapmasını söylüyor. Yani gerçekten bir şey yapmadan neyi sileceğini söylüyor. Bu work dizinlerini kaldıracaktı.

Bu mantıklı görünüyor. Yani aynı komutu tekrar yapacağım; ama _"-n"_ yerine gerçekten temizleme yapmak için _"-f"_ yapacağım. Bu sefer tüm bu dizinleri gerçekten kaldırdı. Work dizinlerine girip bakarsam, artık çok daha hafif görünüyor. Harika.

İşte yerel work dizinlerinizin hepsini önbelleği tamamen yok etmeden oldukça güvenli bir şekilde temizlemenin yolu. Yani isterseniz hâlâ devam edebilirsiniz.

Bu bayrakların ne olduğunu unutursanız, her Nextflow komutu için _"nextflow help"_ ve ardından komutun adını yapabilirsiniz. _"nextflow help clean"_ yaparsam, tüm farklı seçenekleri görebilirsiniz: _-after, -before, -but_, bu temizleme davranışını yapılandırmanın tüm farklı yolları. Oldukça harika.

## Özetle

Tamam, Hello Nextflow'un birinci bölümünün sonu bu. Kursa oldukça yoğun bir başlangıç; ama umarım artık bir Nextflow betiğinin nasıl göründüğü hakkında oldukça iyi bir anlayışınız var: farklı temel parçalar, süreçler, iş akışları, çıktılar ve parametreler. Bunları komut satırından temel geçersiz kılmalarla nasıl yapılandıracağınızı, dinamik bir betikle dinamik bir girdi bloğunu nasıl oluşturacağınızı ve tüm iş yükü çalıştırmalarınızı nasıl yöneteceğinizi biliyorsunuz: daha önce ne çalıştırdığınızı görme, devam ettirme, temizleme. Çok şey var. Uzun bir yol kat ettiniz. Bu yüzden kısa bir mola verip etrafta biraz dolaşmak ve bir fincan çay içmek istiyorsanız, şimdi tam zamanı. Hak ettiniz.

Buradan itibaren temelde bu temel üzerine inşa ediyoruz. Bunu nasıl daha karmaşık, daha güçlü hale getirebiliriz? Nasıl daha esnek hale getirebiliriz? Analizimizi ölçekte yapmak istediğimiz şeyleri nasıl yapabiliriz?

## Quiz

Şimdi web sayfasında birinci bölüme, hello world'e aşağı kaydırırsanız küçük bir quiz göreceksiniz; bu Nextflow eğitiminin bu sürümü için yaptığımız yeni bir şey. Üzerinden geçip bu bölümde yaptığımız tüm materyali anladığınızı kontrol etmek için kendinizi sınayabilirsiniz.

Bu bize gönderilmiyor ya da başka bir şey değil; yalnızca tarayıcınızda saklanıyor. Yani cevaplarınızın ne olduğunu bilmiyoruz; ama bir şeyi kaçırıp kaçırmadığınızı ya da yanlış anlayıp anlamadığınızı kontrol etmek için küçük bir öz değerlendirme. Ve istediğiniz kadar deneyebilirsiniz.

Benim gibi belki VS Code örneğinizde terminalde kalmak istiyorsunuzdur; bu durumda _quiz_ komutunu yazıp hangi bölümde olduğunuzu söyleyebilirsiniz. _"Hello World"_ yapıyoruz ve ardından web tarayıcısındaki quiz sorularının tamamen aynısını, ama yalnızca terminalinizde yapabilirsiniz.

Harika. Tamam. Umarım bunun tadını çıkarırsınız. Biraz eğlenin ve bir dakika sonra Nextflow kanalları hakkında konuşmak için bir sonraki bölümde görüşürüz.
