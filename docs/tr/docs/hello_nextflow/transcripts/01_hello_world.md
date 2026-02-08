# Bölüm 1: Hello World - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [ders materyaline](../01_hello_world.md) geri dönünüz.

    Transkriptte gösterilen bölüm numaraları yalnızca yönlendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba ve tekrar hoş geldiniz.

Şu anda "Hello Nextflow" kursunun "Hello World" adlı birinci bölümündesiniz. Bu bölümde, Nextflow'un en temel kavramları hakkında biraz anlayış geliştirmeye başlayacağız.

Umarım artık Codespaces'te veya eşdeğer bir yerde VS Code çalışır durumda kurulmuşsunuzdur ve Workspace'inizde Explorer'da tüm bu farklı dosyalarla birlikte Hello Nextflow klasörünüz vardır.

Terminalde Bash kullanarak çok basit şeyler yaparak başlayacağız, sonra aynı şeyleri Nextflow içinde yapıp yapamayacağımızı göreceğiz, böylece sözdiziminin nasıl göründüğüne dair bir fikir edineceksiniz.

## 0. Isınma

O halde gerçekten basit bir şekilde başlayalım. Terminale bir şey yazdırmak için sadece "echo" ile başlayalım. "Hello World". Enter'a basıyorum ve bu terminale gidiyor. Hello World. Umarım bu kursu izleyen herkes için sürpriz değildir.

Tamam, bununla bir şeyler yapalım. Bunu terminale yazdırmak yerine, bir dosyaya yazalım. Klavyemdeki yukarı ok tuşuna basacağım, bu Bash geçmişinde dolaşır, böylece bana son komutu verir ve bunun sonuna küçük büyüktür sembolü ekleyeceğim, bu da çıktıyı bu komuttan bir dosyaya yönlendirir ve buna output.txt adını vereceğim.

Tekrar Enter'a basıyorum, bu komutu çalıştırmak için, bu sefer terminalde hiçbir şey yok, ama sol tarafta output.txt adlı yeni dosyanın burada ortaya çıktığını görebiliriz.

Bunu terminalde cat gibi bir şeyle görüntüleyebiliriz. Yani cat output.txt ve gerçekten de "Hello World" diyor. Ayrıca çift tıklayarak VS Code'da kod düzenleyicide açabiliriz.

## 1.1. Kodu inceleyin

Pekala. Size söylemiştim, basitti. Sırada ne var? Bu süreci almayı ve tekrar yapmayı deneyelim, ama bu sefer, bunu Nextflow içinde yapalım.

Söylediğim gibi, bu kurstaki tüm farklı bölümler bir betikle başlar ve bu bölümün adı Hello World. O halde Hello World'ü bulacağım. Tek tıkladığımda önizleme yapıyor, burada düzenleyicide açmak için çift tıklayacağım. Ve terminali hızlıca kaldıracağım.

Şimdi bu çok basit bir betik, yani olabildiğince basit. Sadece 22 satır uzunluğunda ve temelde aynı şeyi yapıyor. Aslında. Bunun bir kısmı tanıdık gelmelidir. Az önce yazdığımız bash komutunu görebiliriz. Bash komutumuzu orada bir dosyaya yönlendiren.

Peki, başka ne var? Ayrıca, bu dosyada, Nextflow'un bazı temel kavramlarını görmeye başlayabiliriz. Burada kırmızı renkte bir süreç (process) ve bir iş akışı (workflow) var. Bunlar Nextflow'da özel anahtar kelimeler ve özel terminolojidir.

## 1.1.1. Süreç tanımı

Bir iş akışı içindeki farklı süreçler, iş akışınızın farklı mantıksal birimlerini sarar. Her süreç bir şey yapar.

Çalıştırdığımızda, bir görev veya birden fazla görev oluşturur, bunlar bir pipeline'ın gerçek yapma adımlarıdır. Tüm süreçler daha sonra alt kısımda gördüğümüz workflow bloğu içinde düzenlenir ve bu durumda sadece o tek süreci çalıştırır.

Süreç adı bu anahtar kelimenin ardından gelir ve bu temelde herhangi bir şey olabilir. Ve sonra sürecin içerikleri bu kıvırcık parantezlerin içindedir.

Bir süreç için gerçekten sadece bir gereksinim var, o da bir tür script veya exec bloğu içermesidir. Bu burada üçlü tırnak içindedir ve bu, pipeline'ı çalıştırdığımızda çalışma dizinine yazılan bash betiğidir ve gerçekte bilgisayarınızda veya sunucunuzda çalışan şeydir.

Bu tipik olarak bash'tir, ama burada en üstte farklı bir hash bang da koyabilirsiniz ve bu bir Python betiği veya bir R betiği olabilir. Farketmez. Bu betikteki her şey yürütülecektir.

Bu sürece eklediğimiz bir şey daha var, bu da output bildirimidir. Bu, Nextflow'a bu sürecin output.txt adlı bir çıktı dosyası beklediğini söyler. Path olduğunu söylüyor, yani bir dosya gibi ele alınmalı, mesela eğer bu val olsaydı, bir değişken veya değer gibi olduğunu söylerdi.

Bunun bu dosyayı oluşturmadığını unutmayın. Aslında onu üretmiyor. Bu, burada aşağıdaki betik tarafından yapılır. Sadece Nextflow'a bu dosya adıyla bir çıktı dosyası beklemesini söylüyor.

## 1.1.2. İş akışı tanımı

Tamam. Ve sonra altta bir workflow var ve yine bir bildirimimiz var. Bu Main olarak adlandırılmış. Bu, isterseniz betik bloğunun workflow eşdeğeridir. Bir şey yapan workflow'un parçasıdır. Ve bu durumda, sayHello adlı süreci çağırdığımızı söylüyoruz.

Normalde, tabii ki, pipeline'ınız bundan çok daha karmaşık görünecektir. Muhtemelen birden fazla süreciniz olacak ve aralarındaki veri akışını düzenlemek için kanalları kullanacaksınız. Bunun üzerine bu kursun sonraki bölümlerinde geleceğiz, ama şimdilik bu yeterli. Bu geçerli bir pipeline ve çalışması gerekir.

VS Code'da burada DAG önizleme'ye bile tıklayabilirim. DAG veya DAG, pipeline'daki veri akışı yapısının bir temsilidir ve yanda bir mermaid diyagramı olarak render edildiğini görebiliriz. Bu durumda çok basit. Bir kutu var, bu workflow ve bir süreç var, sayHello adlı, ama ilerledikçe bu daha ilginç görünebilir.

## 1.2. İş akışını çalıştırın

Tamam, bu workflow'u çalıştırmayı deneyelim ve ne olduğunu görelim.

Terminali tekrar altta getireceğim, çıktıyı temizleyeceğim ve Nextflow Run yazacağım. Ve sonra sadece betik adını yazacağım, bu hello-world.nf. Ve Enter'a basacağım.

Tamam, üstte bize Nextflow'un çalıştığını ve hangi sürümün çalıştığını ve betik adının ne olduğunu ve her şeyi anlatan standart bir şeyler var.

Ve gerçekten aradığımız önemli şey _burada_, yürütülen farklı görevlerin bir özeti.

Eğer sizinki de bunun gibi küçük bir yeşil tik işaretiyle görünüyorsa, o zaman aferin. İlk pipeline'ınızı çalıştırdınız. Harika.

Bize burada Ran adlı sürecin adını söylüyor, Say Hello olarak adlandırılmıştı ve bir kez çalıştığını ve başarılı olduğunu söyledi. Bu ilerlerken güncellenir, yani daha büyük bir pipeline çalıştırırken, ilerlemeyi burada temsil edilmiş olarak göreceksiniz. Ama bu çok küçük olduğu için, temelde hemen çalışır.

## 1.2.2. Çıktıyı ve günlükleri work dizininde bulun

Şimdi bir Nextflow pipeline'ı çalıştırdığınızda, bu süreçlerin her biri birbirine dikilir ve daha önce söylediğim gibi her süreç görevler oluşturabilir, bir veya birden fazla. Yani bu durumda, bu süreçten tek bir görevimiz vardı. Sadece bir kez çalıştı ve bu, bu görev _hash_'i altında yapıldı.

Nextflow, çalışma dizininizdeki dosyalarla doğrudan ilgilenmez, work adlı özel bir klasör oluşturur. Ve "ls" yaparsam, burada ortaya çıktığını göreceğiz: _work_, ve bunun içinde çalışan her tek görev için alt dizinler var. Ve bu, bu hash ile eşleşir. Yani "ls work/c4" yaparsam, ve sonra kısaltılmış, ama 203 ile başlıyor ve bu, pipeline'ı çalıştırdığımızda bu süreç tarafından oluşturulan çalışma dizinidir. Ve bunu yanda da görebilirsiniz.

O dosyaları listelediğimde, output.txt dosyasının oluşturulduğunu görebilirsiniz. Burada da görebilirsiniz. Ve normal "ls"'imde görünmeyen bir grup gizli dosya var.

output.txt'ye tıklarsam, gerçekten de çıktımız var. Harika. Yani pipeline çalıştı.

Esasen tek satırlık bir bash betiği çalıştırmak için oldukça fazla standart kod gibi görünebilir, ama süreçlerimiz daha karmaşık hale geldikçe daha mantıklı olacak. Ve Nextflow ile bu work dizini ve oluşturulan bu dosyalar, gerçekten Nextflow'u bu kadar güçlü yapan şeyin bel kemiğidir.

Her görev, bir pipeline'ın her öğesi diğer tüm görevlerden izole edilmiştir. Tekrarlanabilir. Birbirleriyle çakışmazlar ve her şey paralel olarak çalışabilir. Aslında, alıştığınızda gerçekten güzel bir yol, çünkü bu izolasyon sayesinde girip tam olarak tek bir görev için ne olduğunu görebilir ve hata ayıklayabilirsiniz.

Work dizinindeki bu diğer dosyalara hızlıca bakalım. Yukarıdan aşağıya, _.command.begin_ adlı bir dosyamız var. Bu boş. Bu sadece Nextflow'un oluşturduğu, tamam, görevi başlatıyorum diyen bir sentinel dosyasıdır. Orada ilginç bir şey yok.

Sonra _.command.error_, _.command.log_ ve _.command.out_ var. Bunların hepsi çalışan bash komutundan veya bu betikten gelen çıktılardır. Bu standart hata. Bu standart çıktı, ve bu ikisinin çıktıkları gibi birleşimi. Yani mantıksal sırayı alırsınız.

Tamam, bunlar da bu için boştu, o yüzden çok ilginç değil, ama _.command.run_'a geldiğinizde işler daha ilginç hale geliyor.

Bu tipik olarak çok uzun bir betiktir. Ve Nextflow'un aslında yürüttüğü budur. Buraya girerseniz, Nextflow'un tüm iç mantığını görmeye başlayacaksınız ve ne yaptığını ve sürecinizi nasıl yürüttüğünü göreceksiniz. Bu, nerede çalıştırdığınıza bağlı olacaktır, yerel olarak mı çalıştırıyoruz yoksa SLURM'a bir görev olarak mı gönderiyoruz, bu durumda üstte SLURM başlıkları olacak. Tüm bu farklı kurulumlar.

Genel olarak, bu dosyaya gerçekten bakmaya ihtiyacınız yok. Nextflow tarafından otomatik oluşturulur ve içinde pipeline'ınıza özgü gerçekten özellikle benzersiz bir şey yok. Ama bu gerçekten çalışanın özü.

Bir sonraki çok daha ilginç. _.command.sh_, sürecinizden gelen oluşturulmuş betiktir ve burada Nextflow'un Bash başlığını eklediğini görebilirsiniz, ve sonra betik bloğumuzdaki komutumuz yürütüldü.

Ve _.command.run_ dosyasının yaptığı tek şey bu _.command.sh_ dosyasını çalıştırmaktır.

Bu gerçekten kullanışlı bir tanesidir, genellikle bir şeylerin hatasını ayıklamaya çalışırken en çok baktığınız ve Nextflow pipeline'ınızın mantığının beklediğiniz şeyi yaptığını kontrol ettiğiniz şeydir.

Son olarak, _.exitcode_ adlı bir dosyamız var ve bu sadece bir görevden çıkış kodunu yakalar, bu durumda başarılıydı. Yani çıkış kodu sıfırdı.

Eğer bir şeyler ters giderse, belleğiniz biterse veya başka bir şey ve başarısız olursa, o zaman bu ne yanlış gittiğini anlamak için çok kullanışlıdır.

## 1.3. İş akışını tekrar çalıştırın

Work dizinleri hakkında anlaşılması gereken bir şey daha var, eğer bu pipeline'ı tekrar tekrar çalıştırmaya devam edersem, yani _"nextflow run hello-world.nf"_ yaparsam, tamamen aynı şeyi yapacak, ama bu sefer yeni bir görev kimliği olacak. Buradaki bu hash'in farklı olduğunu görebilirsiniz ve şimdi work'e bakarsam, iki hash dizini var. Ve bunlar yine birbirinden ayrıdır.

Yani her Nextflow workflow'u çalıştırdığınızda, önbelleği kullanan resume'u kullanmadığınız sürece, daha sonra değineceğiz, bu süreçleri birbirinden ayrı olan yeni work dizinlerinde yeniden çalıştıracaktır. Herhangi bir dosya adı çakışması yaşamayacaksınız, bunun gibi hiçbir probleminiz olmayacak. Her şey izole ve temiz.

Ve bu dizine girersek, tüm aynı dosyaları ve sıfırdan yeniden oluşturulmuş aynı _output.txt_'yi görebilirsiniz.

## 2. Çıktıları yayınlayın

Tamam, bu Nextflow'un kendisi için harika, pipeline'ınızı çalıştırırken tüm şeyler birbirinden ayrı ve temiz olabilsin ve yönetilebilsin diye.

Ama sonuçlarınızı keşfetmeye çalışan bir kişiyseniz, bu süper kullanışlı değil. Sonuç dosyalarınızı bulmaya çalışarak binlerce farklı work dizininde dolaşmak gerçekten istemezsiniz. Ve gerçekten de öyle olmamalısınız. Work dizinleri, dosyalarınızın oluşturulduğu son durum olması anlamına gelmez.

Bunu dosyalarımızı yayınlayarak yaparız.

## 2.1.1. sayHello sürecinin çıktısını bildirin

O halde betiğimize geri dönersem, workflow bloğumuzda çalışacağız. Hangi dosyaları beklediğini, hangi dosyalarla ilgilendiğimizi söyleyeceğiz ve sonra altında output block adlı yeni bir blok oluşturacağız.

Bu, Nextflow'un sözdizimi ayrıştırıcısıyla birlikte gelen ve Nextflow'un 26.04 sürümünde varsayılan olan yeni sözdizimdir. Yani daha önce Nextflow'u biraz kullandıysanız, bu yeni olan şeylerden biridir.

Yani main bloğunu aldık ve sonra publish diyeceğim ve Nextflow'a yayınlamadan ne bekleyeceğini söyleyeceğim. Buna _first_output_ diyeceğiz ve _sayHello.out_ diyeceğiz.

Yanlışlıkla orada bir yazım hatası yaptım, ama bu aynı zamanda Nextflow VS Code uzantısının bazı özelliklerine de dikkat çekmenin iyi bir fırsatı. Hemen bunun altında küçük dalgalı bir kırmızı çizgi verdiğini görebilirsiniz, bir şeylerin yanlış olduğunu söylüyor. Ve üzerine gelirseniz, bu değişken tanımlanmamış diyor. Ne olduğunu bilmiyorum.

Bu durumda oldukça açık, bir yazım hatası yaptım. sayHello yazmak istedim ve sonra dalgalı çizgi kayboldu.

Şimdi mor. Nextflow sözdizimi ayrıştırıcısı bunun bir süreç olduğunu biliyor ve üzerine geldiğimde, bu sürecin nasıl göründüğünün küçültülmüş bir temsilini veriyor. Yani çok hızlı bir şekilde bir bakışta hiçbir girdi almadığını ve bize bu çıktıyı verdiğini görebilirim. Yani VS Code'da bu uzantıyla çalışmak, kod yazarken size çok fazla bağlamsal bilgi veriyor.

Bu süreçten çıktıya _.out_ sözdizimi ile başvurabileceğimizi unutmayın. Ve şu anda buna istediğimiz her şeyi çağırabiliriz, bu sadece keyfi bir değişken adı.

## 2.1.2. Betiğe bir output: bloğu ekleyin

Önemli olan, bunu burada yeni bloğumuzu yaptığımızda olur ve bu artık workflow bloğunun altında, artık workflow içinde değiliz. Tekrar kıvırcık parantezler. Ve bu sadece Nextflow'a workflow tarafından oluşturulan tüm dosyaları nereye koyacağını söylediğimiz yer.

Şimdi burada oluşturduğum bu değişken adını alacağım ve oraya koyacağım ve bunun için bazı kıvırcık parantezler koyacağım. Ve Nextflow'a bir path kullanmasını söyleyeceğim. Oops. Path, tırnak işaretleri içinde. Ve nokta kullanacağım. Bu sadece Nextflow'a dosyayı results dizininin köküne koymasını söyler. Yani hiçbir alt dizin veya başka bir şey değil.

Workflow'umuzu tekrar çalıştırmayı deneyelim. _"nextflow run hello-world.nf"_ yaparsam, umarım temelde tamamen aynı görünmelidir. Nextflow ile gerçekten hiçbir şey değişmedi. Aynı şeyleri çalıştırıyor. Sadece onları tekrar work dizinlerinde yapıyor.

Ama şimdi _"ls results/"_ yaparsam, burada workflow yayınlama için varsayılan temel dizin olan results adlı yeni bir dizin oluşturulduğunu göreceksiniz. Ve içinde _output.txt_ adlı bir dosya var.

_"ls -l results"_ yaparsam, bunun aslında work dizinine yumuşak bağlantılı olduğunu göreceksiniz. Yani bu gerçek bir dosya değil, work dizinine bağlantılı ve oradaki tüm dosyaları bizim için topladı.

## 2.2. Özel bir konum ayarlayın

"Results" bu path için varsayılan addır. Workflow'u tekrar çalıştırırsam ve bu sefer _dash_ tek tire yaparsam, bu bir çekirdek Nextflow seçeneği olduğu için. _"-output-dir **my** results"_. Kısaca sadece _"-o"_ da yapabilirdi. O zaman dosyaların nerede saklandığı için farklı bir temel dizin ayarlayacak ve bir kez daha, burada _myresults/_'ta, şimdi bir _output.txt_ var.

Bu harika, ama muhtemelen tüm dosyaların sadece kökte olmasını istemiyoruz. Biraz organizasyon istiyoruz, bu yüzden burada istediğimiz her şey olarak adlandırılmış bir alt dizin de oluşturabiliriz. Diyelim ki _"path 'hello_world'"_ ve bunu tekrar çalıştırıyorum. _"nextflow run hello-world.nf"_. Results dizinine bir alt dizine gitmeli ve gerçekten de, şimdi burada üstte results altında _hello_world/_ var ve _output.txt_ var.

Fark edilmesi gereken önemli bir şey, eski _output.txt_ dosyası hala orada. Bunu yaptığınızda results dizini silinmez. Sadece yeni dosyalar oraya kopyalanır. Aynı dosya adına sahiplerse zaten orada olan dosyaların üzerine yazarlar, ama eskileri temizlemezler. Bu yüzden pipeline'ları yeniden çalıştırdığınızda biraz dikkatli olmanız gerekir. Eğer zaten orada olan dosyaların üzerinde olmalarını istemiyorsanız. Boş bir dizin kullandığınızdan emin olun.

## 2.3. Yayınlama modunu kopyala olarak ayarlayın

Tamam, bu dosyaların yumuşak bağlantılar olduğunu söyledim, yani _"ls -l results/hello_world/"_ yaparsam, work dizinine yumuşak bağlantı yaptığını görebilirsiniz. Bu genellikle HPC gibi bir şey üzerinde çalışıyorsanız iyi bir şeydir ve bunlar gerçekten büyük dosyalardır ve onları çoğaltmak istemezsiniz, çünkü bu, dosyaların dosya sisteminde yalnızca bir kez saklandığı anlamına gelir.

Ancak, bu work dizinini silerseniz demektir: _"rm -r work"_ yaparsam ve oluşturulan tüm ara dosyaları temizlersem. Şimdi, bu dosyayı okumaya çalışırsam _"results/hello_world/"_. Artık var olmayan bir dosyaya yumuşak bağlantı olarak işaret edecek ve veriler sonsuza kadar gitti ve geri alınamaz, bu belki pek hoş değil.

Bu yüzden genellikle, yapabilirseniz yumuşak bağlantı yerine dosyaları kopyalamanın iyi bir uygulama olduğunu söylerim, çünkü daha güvenlidir. Sadece bu work dizinlerini silmediğiniz sürece iki kat daha fazla disk alanı kullanacağının farkında olun.

Bunu output bloğuyla yapmak için, burada ilk output'a gideceğim. Daha önce path'i ayarladım ve şimdi mode'u ayarlayacağım ve yazdığım gibi görebilirsiniz, VS Code uzantısı önerilerde bulunuyor, bunun burada bir output yönergesi olduğunu biliyor. Ve copy diyeceğim. Kaydet'e basıyorum.

Workflow'u yeniden çalıştıralım. Dosyaları tekrar oluşturacak, yeni work dizini.

Şimdi, _"ls -l results/hello_world/"_'a gidersem bunun gerçek bir dosya olduğunu görebilirsiniz ve artık yumuşak bir bağlantı değil ve Nextflow bunu kopyaladı. Bilmek iyi. Yani path ve mode oldukça sık yazacağınız şeylerdir.

Şimdi, tabii ki, bu çok basit. Bunu ilerledikçe daha karmaşık ve güçlü hale getireceğiz ve bu şeyleri nasıl dinamik yapacağınızı ve çok ayrıntılı olmayacağınızı göreceksiniz.

## 2.4. Süreç düzeyinde publishDir yönergeleri hakkında not

Şimdi, buna başlarken söylediğim gibi, bu oldukça yeni bir sözdizimi biçimidir. Bunu kaydettiğim sırada yalnızca Nextflow'un en son sürümlerinde mevcuttur ve buna Workflow Outputs denir.

Bunu kullanırsanız, harika. Nextflow içinde, Nextflow Lineage gibi, oluşturuldukları sırada bu dosyaların mirasını takip etmeye yardımcı olmak için başka birçok havalı özelliğin kilidini açar ve yakında 26.04'te varsayılan olacak. Ve gelecekte daha sonraki bir tarihte, bu workflow'larınızı yazmanın tek yolu olacaktır.

Ancak, şu anda bu geçiş aşamasındayken, belki de publishDir adı verilen bir şey kullanan vahşi doğada pipeline'lar görebilirsiniz, bu bunu yapmanın eski yoludur ve bu workflow ve output düzeyinde değil, süreç düzeyinde tanımlanır.

Ve bu bildiri temelde aynı şeyi söyler. Sonuç dosyalarını results adlı bir dizine yayınla diyor ve bir kopyalama modu kullan. Yani sözdizimin çok benzer olduğunu görebilirsiniz. Ama şimdi yeni pipeline'lar yazarken, AI sonuçlarında veya belgelerde veya diğer pipeline'larda görsem bile, bu publishDir yönergesini kullanmamaya çalışın, çünkü bu bunu yapmanın eski yolu.

2026'da hepimiz workflow outputs kullanmalıyız.

Bütün bunlar belgelenmiştir, bunu yapıyorsanız ve daha önce Nextflow kullandıysanız, burada Nextflow dokümanlarına gidebilirsiniz, nextflow.io/docs/. Ve öğreticilere kaydırırsam, _Workflow Outputs'a Geçiş_ adlı bir öğretici var.

Gerçekten çok iyi. Tüm sözdizimini, eski sözdizime nasıl eşdeğer olduğunu, neden değiştirdiğimizi ve bir zaman çizelgesi ve her şeyi anlatıyor. Ve bir sürü örnekle birlikte tüm farklı senaryoları ele alıyor. Böylece mevcut Nextflow kodunu yeni sözdizime kolayca dönüştürebilirsiniz.

## 3.1. sayHello sürecini bir değişken girdisi bekleyecek şekilde değiştirin

Tamam, basit betiğimizi aldık, bir süreç çalıştırıyor, bir dosya oluşturuyor, Nextflow'a bunun bir çıktı olduğunu söylüyor ve sonra Nextflow'a o dosyayı nereye kaydedeceğini söylüyoruz. Bu iyi bir başlangıç.

Ama hepsi sabit kodlanmış olmasaydı daha ilginç olurdu. Yani sonra, bu sürece Nextflow'a, çalışma zamanında bir workflow başlattığımızda kontrol edebileceğimiz bir değişken girdisi alabileceğini nasıl söyleyeceğimizi düşünelim.

Bunun gerçekleşmesi için birkaç farklı şey yapmamız gerekiyor.

İlk olarak, bu sürece bir girdi değişkeni kabul edebileceğini söylememiz gerekiyor ve burada yeni bir bildirim bloğu olarak _input_ yazıyoruz. Ve buna _"val greeting"_ diyeceğiz.

Val kısmı aşağıdaki path'in eşdeğeridir. Nextflow'a bunun bir değişken olduğunu, bu durumda bir string gibi olduğunu söyler. Ve üzerine gelirseniz, bu ne anlama geldiğinin uzantıdan yine size söyler.

Sonra Nextflow'a bununla ne yapacağını söyleyeceğiz. Sadece bir değişken olduğunu söylemek yeterli değildir. Betikteki o değişkeni nasıl kullanacağınızı söylemeniz gerekir. Ve bu yüzden buradaki bu sabit kodlanmış string'i kaldıracağım ve bir değişken koyacağım.

Kıvırcık parantezler olmadan hızlıca yapacağım, sadece bunun izin verildiğini göstermek için ve bu bunu yapmanın eski, stil yolu. Ama şimdi yeni sözdizimi ile, bunu gerçekten bunun gibi kıvırcık parantezlerin içine koymayı öneriyoruz ve bunun Nextflow tarafından burada interpolate edildiğini çok açık hale getirir.

Harika. Yani _"input greeting"_ _$\{greeting\}_ içine girer. Son şey, Nextflow'a workflow düzeyinde bu sürecin artık bir girdi aldığını söylememiz gerekiyor. Ve bunu yapmak için, temelde ona bir değişken vereceğiz.

## 3.2. Kullanıcı girdisini yakalamak için bir komut satırı parametresi ayarlayın

Bunu tekrar Hello World gibi sabit kodlayabiliriz ve bu gayet iyi çalışır, ama açıkçası bize gerçekten hiçbir avantaj sağlamaz. Bunu çalışma zamanında yapılandırabilmeyi istedik, böylece Nextflow'u başlattığınızda CLI'da yapabilmek istiyoruz.

Ve bunu yapma şeklimiz, _params_ adlı özel bir Nextflow kavramıdır. Buna _params.input_ diyeceğiz.

Bunun yaptığı şey, bu input değişkenini CLI'da açığa çıkarmaktır ve orada Nextflow'u başlattığımızda çift tire kullandığımız yerdir.

Buna istediğim her şeyi çağırabilirim, _hello, greeting_ diyebilirim. Farketmez. Orada yaptığım her şey, bir pipeline başlattığımızda CLI seçeneği olarak açığa çıkacaktır. Ve bu Nextflow'un gerçek bir büyüsüdür çünkü bu, workflow betiğinizi bu parametrelerle çok hızlı bir şekilde oluşturabileceğiniz anlamına gelir ve esasen başlattığınızda farklı seçenekleri anında özelleştirmeyi gerçekten kolaylaştırarak pipeline'ınız için özel bir CLI oluşturuyorsunuz.

O halde. Deneyelim. Terminalimize geri dönelim. Burada _"nextflow run"_ komutumuz var. Ve şimdi daha önce gördüğümüz _"params.input"_ ile eşleşen _"--input"_ yapacağım. Sanırım dokümanlarda Fransızca. Geraldine Fransızca konuşmayı sever. Ben İsveççe yapacağım çünkü İsveç'te yaşıyorum. Bu yüzden, "_Hej Världen_" diyeceğim ve Enter'a basacağım.

Tek tırnak veya çift tırnak kullanabilirsiniz, sadece Bash'in onu nasıl yorumladığını etkiler.

Nextflow pipeline'ını tamamen aynı şekilde çalıştırır. Çalışma dizininin ve her şeyin aynı olduğunu görebilirsiniz. Ama şimdi _"results/hello_world/output"_'a gidersem. Burada güzel İsveççe'mizi görebiliriz.

Yani CLI'dan bir parametreye dinamik olarak bir girdi geçirdik. Bunu sürece bir girdi olarak geçirdik ve süreç bunu yorumladı ve bir betik bloğuna koydu, bu da o betik sonucunun çıktısını dinamik olarak değiştirdi. Oldukça havalı.

Burada çok, az sözdizimi ile oldukça karmaşık mantık. Ve umarım bunun artık nasıl ölçeklendiğini görebilirsiniz. Ve bu, pipeline'larımızın mantığını ve özelleştirilebilirliğini Nextflow betiğine gerçekten nasıl inşa ettiğimizdir.

## 3.4. Komut satırı parametreleri için varsayılan değerler kullanın

Tamam, bu harika. Ama şimdi problem, bu pipeline'ı her çalıştırdığımda, çalışması için tire, input yapmam gerekiyor.

Bu parametre olmadan çalıştırmaya çalışırsam, şimdi Nextflow bu parametreye ihtiyacı olduğunu ve ayarlanmadığını söyleyerek bir hata verecek. ve bu yüzden ne yapacağını bilmiyordu.

Bu arada havalı yeni bir şey. Geçmişte, Nextflow sadece boş bir string ile çalışırdı ve anlaşılması zor olan her türlü garip hata alırdınız. Ama yeni Nextflow sözdizimi ayrıştırıcısında, biraz daha dikkatli ve size hemen söylüyor.

Yani her seçeneği belirtmek her zaman istemeyiz. Makul varsayılanlar belirtmek iyi bir uygulamadır. Peki bunu betiğimizde nasıl yaparız?

Bunu yazdığımızda, _params.input_'u doğrudan kullandığımız yere koyduğumuzu fark edeceksiniz. Yani bariz çözüm bir varsayılan tanımlarız ve bunu workflow'daki özel params bloğunda, betiğin en üstünde yaparız. Bu burada workflow betiğindedir.

Yine, burada bazı yeni sözdizimler var, bu yüzden dikkat edin. Bu gerçekten havalı şeyler. Burada beklenecek parametrenin adını aldık.

Ve sonra bu iki nokta karakterinden sonra, değişkenin tipini tanımlıyoruz. Bunu yapmak zorunda değilsiniz, sadece boş bırakabilirsiniz, ama gerçekten güzel. Nextflow'a bir string beklediğimizi söyler ve bunu öyle davran.

Örneğin, bunun yerine bir sayı istersek, float yazabilirdik ve bu da kayan noktalı bir sayı istediğimizi söylerdi. Ve bununla çalıştırmaya çalışırsak, o zaman bir hata verecektir. Eğer ona bir float olmayan bir string verirsek. Ve aynı şekilde geçirecektir. String yaparsak, o zaman bunun bir string olduğunu bilir. Ve baştaki sıfırları olsa ve tamamen sayısal olsa bile, onu yine de gerçek bir string olarak geçecektir.

Yani bu tip güvenliği Nextflow'un çok yeni bir özelliğidir, ama kodunuzu yazarken ve çalıştırırken daha güvenli hale getirmek için gerçekten güçlüdür.

Sonra bundan sonra eşittir sembolünü aldık ve sonra burada varsayılan değer var. Nextflow orijinal olarak Barcelona'da yazıldı, bu yüzden burada biraz, İspanyolca'mız olması uygun görünüyor, varsayılan olarak _"Holà mundo!"_.

Tamam betiği kaydedeceğim, geri döneceğim, betiği _--input_ olmadan tekrar çalıştıracağım. Ve bu sefer çalışmalı ve _results_'ta yeni dosyamızı oluşturacak. Ve bu dosyada şimdi _"Holà mundo!"_ diyor.

Bu sadece bir varsayılan, yani daha önce yaptığımız aynı şeyi hala yapamayacağımız anlamına gelmez. Geri dönüp burada eski betiğimi bulursam, _"Hej Världen"_, komut satırında _--input_ yaptığım için, bu o varsayılanın üzerine yazacak ve onu output.txt dosyasında tekrar kullanacaktır.

Yani betikteki bu sadece ayarladığım varsayılan değerdir.

Workflow'umuzu daha karmaşık hale getirip daha fazla parametre içerecek şekilde oluşturdukça, betiğin üstündeki bu params bloğu hepsini bir yerde toplamaya başlayacak.

Ve betiğinizde oldukça güzel bir simetri elde edersiniz, burada etkin bir şekilde tüm workflow girdileriniz burada ve workflow çıktılarınız en altta. Ve workflow'unuzun dış dünyaya arayüzünün ne olduğu çok açık. Böylece yeni sözdizimi ile yeni bir pipeline'ı çok hızlı bir şekilde alabilir ve nasıl kullanacağınızı anlayabilirsiniz.

Son bir havalı şey. Bununla bir varsayılan değer ayarlamak zorunda değiliz. Eğer params input yaparsak ama bir varsayılan değer ayarlamazsak, o zaman Nextflow'a bu parametrenin gerekli olduğunu söyler ve yine, pipeline onsuz çalışmayacak, ama null hakkında bir şey yerine daha kullanışlı bir hata mesajı verecektir.

Yani bir girdinin gerekli olduğunu beklediğimizi söylüyor, ama komut satırında belirtilmedi. Çok güzel.

Tamam, umarım şimdi Nextflow pipeline'ınızı değişken girdiler ve parametrelerle nasıl kuracağınız, varsayılanı nasıl ayarlayacağınız, tipleri nasıl ayarlayacağınız, Boolean doğru yanlış bayrağı veya bir tamsayı veya burada farklı tipler olabilir konusunda açıktır. Bunları workflow'unuza nasıl geçireceğiniz, nereye gittiği ve sonra sürecinize interpolate edileceği. Ve ayrıca Nextflow'u başlattığınızda komut satırında bunları nasıl özelleştireceğinizi de biliyorsunuz. Bu basit bash komutumuza göre daha ilginç görünmeye başlıyor.

## 4. İş akışı yürütmelerini yönetin

Tamam. Sırada ne var? Bu bölümün son kısmı için, tüm farklı workflow yürütmelerini nasıl yöneteceğimiz hakkında biraz konuşacağız. Kenar çubuğuma ve work altındaki Explorer'a bakarsanız, bir sürü farklı pipeline çalıştırdım ve bu work dizinleri oldukça uzunlaşıyor, bir sürü var.

Ve diğer şey, daha önce söylediğim gibi, bu pipeline'ı her yeniden çalıştırdığımda, yeni bir work dizinleri seti oluşturuyor ve tüm süreçleri sıfırdan yeniden çalıştırıyor, bu iyi bir şey. Bu amaçlanan davranıştır. Tekrarlanabilir ve her şeyi taze yeniden üretiyor. Ama açıkçası, çok uzun süren süreçler çalıştırıyorsanız, pipeline'ınız ortada çöktüyse veya pipeline'ın sonunda bir şeyi değiştirdiyseniz, pipeline'ınızı her zaman baştan başlamak zorunda kalmak sinir bozucu.

## 4.1. Bir workflow'u -resume ile yeniden başlatın

Neyse ki, Nextflow daha önce ne çalıştırıldığını ve ne mevcut olduğunu bilmede gerçekten, iyidir ve o eski sonuçları yeniden kullanmak için çok, basittir. Komutun sonuna sadece yeni bir bayrak ekliyoruz _"-resume"_.

Şimdi, input'ta iki tire olduğuna dikkat edin çünkü bu parametre. Resume'da sadece bir tire var çünkü bu bir çekirdek Nextflow seçeneği.

İnsanları her zaman yanıltır, Nextflow'u uzun süredir kullansanız bile. Bu yüzden her zaman bir veya iki tire hatırlayın. Çekirdek bir Nextflow seçeneği olup olmadığına bağlı.

Tamam, şimdi _-resume_ yapıyorum ve tamamen aynı workflow'u tekrar çalıştırıyorum. Ve bu sefer, bir önemli farkla oldukça aynı görünmelidir.

Buradaki çıktıda, sonuçların önbelleğe alındığını görebilirsiniz. Ve aslında, buradaki bu görev hash'i önceki çalıştırmayla tamamen aynı ve sadece o work dizinini bütünüyle yeniden kullandı. Girdiler ve çıktılar ve betik hepsi değiştirilmemişti. Ve bu yüzden sadece o dosyayı oradan aldı ve pipeline'da aşağı akış adımları olsaydı, onları pipeline'daki bir sonraki adıma geçirirdi.

Yani hala tüm pipeline'ı baştan sona çalıştırıyor, ama yapabildiği her bir görev için önbelleğe alınmış sonuçları kullanıyor.

Şimdi, _-resume_ yaptığınızda, sadece çalışma dizininizdeki son pipeline çalıştırmasını sürdürür, o her neyse. Ama aslında orada yaptığınız herhangi bir önceki çalıştırmadan sürdürebilirsiniz. Ve şimdi oldukça fazla yaptık.

## 4.2. Geçmiş yürütmelerin günlüğünü inceleyin

Hepsine bakmak için, _"nextflow run"_ yerine _"nextflow log"_ yapabiliriz ve bu bize tüm bunların güzel bir çıktısını verecek.. Görebilmemiz için ekranımı biraz daha küçük yapmam gerekiyor, ne zaman yaptığımız tüm bu farklı çalıştırmalar, oturum kimliği, komut ve her şey.

Ve buraya bakabiliriz ve bunlardan herhangi birinin çalıştırma adını alabiliriz ve sonra o spesifik olanlardan birini sürdürebiliriz. Yani geri dönebilirim ve _hungry_ekeblad_ adlı olanı sürdürebilirim. Ve bunu sadece resume'dan sonra koyuyorum.

Bu arada merak ediyorsanız, tüm bu sıfatlar ve bilim insanı adları Nextflow kaynak kodundadır. Gidip onu bularak ve favori bilim insanınızı ekleyerek Nextflow'a ilk pull isteğinizi almanın gerçekten iyi bir yolu.

Ve neyse, bunu yaptım ve geri gitti ve bu workflow çalıştırmasından önbelleğe alınmış sonuçlara baktı, hala onları yeniden kullanabileceğini anladı ve yaptı. Yani tekrar önbelleğe alınmış sonuçları aldım.

## 4.3. Eski work dizinlerini silin

Bu harika. Peki ya bu work dizinlerini temizlemek istersem? Burada bir sürü var. Bir sürü dosya var. Belki son birkaç pipeline çalıştırmasından sürdürmek istediğimi kesin olarak biliyorum, ama bunlardan önceki tüm olanları umursamıyorum.

O zaman burada birini seçebilirim ve başka bir Nextflow komutu kullanabilirim, bu _"nextflow clean"_ ve _"nextflow clean"_ yapabilirim, _"-before"_ yapacağım ve bu durumda _reverent_pike_ olan belirli çalıştırma adı ve _"-n"_ yapacağım, bu Nextflow'a sadece kuru bir çalıştırma yapmasını söyler. Yani sadece neyi sileceğini söyler. Aslında hiçbir şey yapmadan, yani bu work dizinlerini kaldıracaktır.

Bu mantıklı görünüyor. Yani aynı komutu tekrar yapacağım, ama _"-n"_ yerine temizliği gerçekten yapmak için _"-f"_ yapacağım. Ve bu sefer gerçekten tüm bu dizinleri kaldırdı. Ve girip work dizinlerine bakarsam, şimdi çok daha hafif görünüyor. Harika.

Yani tüm yerel work dizinlerinizi önbelleği tamamen yok etmeden oldukça güvenli bir şekilde nasıl temizleyeceğiniz budur. Böylece isterseniz hala sürdürebilirsiniz.

Bu bayrakların ne olduğunu unutursanız, her Nextflow komutu için _"nextflow help"_ yapabilirsiniz ve sonra komutun adını. Yani _"nextflow help clean"_ yaparsam, tüm farklı seçenekleri görebilirsiniz: _-after, -before, -but_, bu temizleme davranışını yapılandırmanın tüm farklı yolları. Oldukça havalı.

## Özet

Tamam, bu Hello Nextflow'un birinci bölümünün sonu. Kursa oldukça yoğun bir başlangıç, ama umarım şimdi bir Nextflow betiğinin nasıl göründüğü hakkında oldukça iyi bir anlayışınız vardır; farklı ana parçalarla, süreçler, iş akışları, çıktılar ve parametreler. Bunları komut satırından temel geçersiz kılmalarla nasıl yapılandıracağınızı, dinamik bir betikle dinamik bir girdi bloğu nasıl yapacağınızı biliyorsunuz ve tüm iş yükü yürütmelerinizi nasıl yöneteceğinizi biliyorsunuz: daha önce ne çalıştırdığınızı görmek, sürdürmek, temizlemek. Çok şey var. Uzun bir yol kat ettiniz. O halde bir mola verip kısa bir yürüyüş ve bir fincan çay içmek isterseniz, muhtemelen şimdi iyi bir zaman. Bunu hak ettiniz.

Buradan itibaren, temelde bu temel üzerine inşa ediyoruz. Bunu nasıl daha karmaşık, daha güçlü hale getirebiliriz? Nasıl daha esnek yapabiliriz? Analizimizi ölçekte yapmak istediğimiz şeyleri yapmak.

## Quiz

Şimdi web sayfasındaki birinci bölüm, hello world'e kadar aşağı kaydırırsanız küçük bir sınav göreceksiniz ve bu, Nextflow eğitiminin bu versiyonu için yaptığımız yeni bir şey. Ve bu bölümde yaptığımız tüm materyali anladığınızı kontrol etmek için kendinizi sınayabilirsiniz.

Bu bize gönderilmez veya başka bir şey değildir, sadece tarayıcınızda saklanır. Yani cevaplarınızın ne olduğunu bilmiyoruz, ama bu sadece bir şeyi kaçırmadığınızdan veya yanlış anlamadığınızdan emin olmak için küçük bir kendi kendini kontrol. Ve istediğiniz kadar deneyebilirsiniz.

Eğer benim gibiyseniz, belki VS Code örneğinizdeki terminalde kalmak istersiniz, bu durumda _quiz_ komutunu yazabilir ve sonra hangi bölümde olduğunuzu söyleyebilirsiniz. Yani _"Hello World"_ yapıyoruz ve sonra web tarayıcısındaki ile tamamen aynı, sınav sorularını yapabilirsiniz, ama sadece terminalinizde.

Havalı. Tamam. Umarım bunun tadını çıkarırsınız. Biraz eğlenin ve, bir dakika sonra kanallar hakkında konuşmak için sonraki bölümde görüşürüz.
