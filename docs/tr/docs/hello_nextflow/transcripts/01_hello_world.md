# Bölüm 1: Hello World - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirme önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../01_hello_world.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba ve tekrar hoş geldiniz.

Şu anda "Hello Nextflow" kursunun "Hello World" adlı Birinci Bölümündesiniz. Bu bölümde, Nextflow'un en temel kavramlarını anlamaya başlayacağız.

Umarım artık Codespaces'te veya eşdeğer bir yerde VS Code çalışır durumdadır ve Explorer'da çalışma alanınızda Hello Nextflow klasörünüz buradaki tüm farklı dosyalarla birlikte bulunmaktadır.

Terminalde Bash kullanarak çok temel şeyler yaparak başlayacağız ve ardından aynı şeyleri Nextflow içinde yapıp yapamayacağımızı göreceğiz, böylece sözdiziminin nasıl göründüğüne dair bir fikir edineceksiniz.

## 0. Isınma

Hadi gerçekten basit bir şeyle başlayalım. Terminale bir şey yazdırmak için sadece "echo" ile başlayalım. "Hello World". Enter'a basıyorum ve bu terminale gidiyor. Hello World. Umarım bu kursu izleyen kimse için sürpriz değildir.

Tamam, bununla bir şeyler yapalım. Bunu sadece terminale yazdırmak yerine, bir dosyaya yazalım. Klavyemdeki yukarı ok tuşuna basacağım, bu Bash geçmişinde gezinir, böylece bana son komutumu verir ve sonuna küçük bir büyüktür sembolü ekleyeceğim, bu da bu komutun çıktısını bir dosyaya yönlendirir ve buna output.txt adını vereceğim.

Tekrar Enter, bu komutu çalıştırmak için, bu sefer terminalde hiçbir şey yok, ama sol tarafta, burada output.txt adlı yeni dosya belirdi.

Bunu terminalde cat gibi bir şeyle görüntüleyebiliriz. Yani cat output.txt ve gerçekten de "Hello World" diyor. Ayrıca çift tıklayabiliriz ve VS Code'daki kod düzenleyicide açılır.

## 1.1. Kodu İnceleyin

Pekala. Size basit olduğunu söylemiştim. Sırada ne var? Hadi bu süreci almayı ve tekrar yapmayı deneyelim, ama bu sefer, bunu Nextflow içinde yapalım.

Dediğim gibi, bu kurstaki farklı bölümlerin hepsi bir betikle başlar ve bu bölümün adı Hello World. Bu yüzden Hello World'ü bulacağım. Tek tıkladığımda önizleme yapıyor, burada düzenleyicide açmak için çift tıklayacağım. Ve terminali hızlıca kaldıracağım.

Şimdi bu çok basit bir betik, yani olabildiğince basit. Sadece 22 satır uzunluğunda ve temelde aynı şeyi yapıyor. Aslında. Bunun bir kısmı tanıdık gelmeli. Az önce yazdığımız bash komutunu görebiliriz. Orada bir dosyaya yönlendirme yapan bash komutumuz var.

Tamam. Başka ne var? Ayrıca, bu dosyada, Nextflow'un bazı temel kavramlarını görmeye başlayabiliriz. Burada kırmızıyla bir süreç ve bir iş akışı var. Bunlar Nextflow'da özel anahtar kelimeler ve özel terminolojidir.

## 1.1.1. Süreç tanımı

Bir iş akışı içindeki farklı süreçler, iş akışınızın farklı mantıksal birimlerini sarar. Her süreç bir şey yapar.

Çalıştırdığımızda, bir görev veya birden fazla görev oluşturur, bunlar bir boru hattının gerçek yapma adımlarıdır. Tüm süreçler daha sonra altta gördüğümüz bir workflow bloğu içinde düzenlenir ve bu durumda sadece o bir süreci çalıştırır.

Süreç adı bu anahtar kelimeyi takip eder ve bu temelde herhangi bir şey olabilir. Ve ardından sürecin içeriği bu süslü parantezler içindedir.

Bir süreç için gerçekten sadece bir gereklilik vardır, o da bir tür betik veya exec bloğu içermesidir. Bu burada üçlü tırnak içindedir ve bu, boru hattını çalıştırdığımızda çalışma dizinine yazılan bash betiğidir ve bilgisayarınızda veya sunucunuzda gerçekten çalışan şey budur.

Bu tipik olarak bash'tir, ancak burada en üstte farklı bir shebang da koyabilirsiniz ve bu bir Python betiği veya R betiği olabilir. Önemli değil. Bu betikteki her şey çalıştırılacaktır.

Bu sürece eklediğimiz bir şey daha var, o da çıktı bildirimi. Bu, Nextflow'a bu sürecin output.txt adlı bir çıktı dosyası beklediğini söyler. Path olduğunu söylüyor, yani bir dosya gibi işlenmeli, örneğin bu val olsaydı, bir değişken veya değer gibi olduğunu söylerdi.

Bunun bu dosyayı oluşturmadığını unutmayın. Aslında onu üretmiyor. Bu aşağıdaki betik tarafından yapılır. Sadece Nextflow'a bu dosya adıyla bir çıktı dosyası beklemesini söylüyor.

## 1.1.2. İş akışı tanımı

Tamam. Ve sonra altta burada bir iş akışımız var ve yine bir bildirimiz var. Bunun adı Main. Bu, isterseniz bir betik bloğunun iş akışı eşdeğeridir. İş akışının bir şey yapan kısmıdır. Ve bu durumda, sayHello adlı süreci çağır diyoruz.

Normalde, elbette, boru hattınız bundan çok daha karmaşık görünecektir. Muhtemelen birden fazla süreciniz olacak ve aralarındaki veri akışını düzenlemek için kanallar kullanacaksınız. Buna bu kursun sonraki bölümlerinde geleceğiz, ama şimdilik bu yeterli. Bu çalışması gereken geçerli bir boru hattıdır.

VS Code'da burada önizleme DAG'ına bile tıklayabilirim. DAG veya DAG, boru hattındaki bir veri akışı yapısının temsilidir ve onu yan tarafta bir mermaid diyagramı olarak işlenmiş görebiliriz. Bu durumda çok basit. İş akışı olan bir kutu ve sayHello adlı bir süreç var, ama ilerledikçe bu daha ilginç görünebilir.

## 1.2. İş akışını çalıştırın

Tamam, bu iş akışını çalıştırmayı deneyelim ve ne olduğunu görelim.

Terminali altta tekrar getireceğim, çıktıyı temizleyeceğim ve Nextflow Run yazacağım. Ve sonra sadece betik adını yazacağım, bu hello-world.nf. Ve enter'a basacağım.

Tamam, en üstte bize Nextflow'un çalıştığını ve hangi sürümün çalıştığını ve betik adının ne olduğunu ve her şeyi söyleyen bazı standart şeyler var.

Ve gerçekten aradığımız önemli şey _burada_, yürütülen farklı görevlerin bir özeti.

Sizinki böyle küçük yeşil bir tik işaretiyle görünüyorsa, tebrikler. İlk boru hattınızı çalıştırdınız. Harika.

Burada bize çalışan sürecin adını söylüyor, Say Hello adındaydı ve bize bir kez çalıştığını ve başarılı olduğunu söyledi. Bu ilerledikçe güncellenir, bu yüzden daha büyük bir boru hattı çalıştırdığınızda, ilerlemeyi burada temsil edildiğini göreceksiniz. Ama bu çok küçük olduğu için, temelde anında çalışır.

## 1.2.2. Çıktıyı ve günlükleri work dizininde bulun

Şimdi bir Nextflow boru hattı çalıştırdığınızda, bu süreçlerin her biri birbirine dikilir ve her süreç, daha önce söylediğim gibi, görevler oluşturabilir, bir veya birden fazla. Yani bu durumda, bu süreçten tek bir görevimiz vardı. Sadece bir kez çalıştı ve bu, bu görev _hash_'i altında yapıldı.

Nextflow, çalışma dizininizdeki dosyalarla doğrudan ilgilenmez, work adında özel bir klasör oluşturur. Ve "ls" yaparsam, burada göründüğünü göreceğiz: _work_, ve bunun içinde çalışan her görev için alt dizinler var. Ve bu, bu hash ile eşleşir. Yani "ls work/c4" yaparsam görebilirsiniz, ve sonra kısaltılmış, ama 203 ile başlıyor ve bu, boru hattını çalıştırdığımızda bu süreç tarafından oluşturulan çalışma dizinidir. Ve onu yan tarafta da görebilirsiniz.

Bu dosyaları listelediğimde, output.txt dosyasının oluşturulduğunu görebilirsiniz. Onu burada da görebilirsiniz. Ve normal "ls" ile gösterilmeyen bir sürü gizli dosya var.

output.txt'ye tıklarsam, gerçekten de çıktımız var. Harika. Yani boru hattı çalıştı.

Esasen tek satırlık bir bash betiği çalıştırmak için oldukça fazla şablon gibi görünebilir, ama süreçlerimiz daha karmaşık hale geldikçe daha mantıklı olacak. Ve Nextflow ile bu work dizini ve oluşturulan bu dosyalar, gerçekten Nextflow'u bu kadar güçlü yapan şeyin omurgasıdır.

Her görev, bir boru hattının her öğesi diğer tüm görevlerden izole edilmiştir. Tekrarlanabilir. Birbirleriyle çakışmazlar ve her şey paralel çalışabilir. Aslında alıştığınızda gerçekten güzel bir yoldur çünkü bu izolasyon sayesinde içeri girip tek bir görev için tam olarak ne olduğunu görebilir ve hata ayıklayabilirsiniz.

Work dizinindeki bu diğer dosyalara hızlıca bir göz atalım. Yukarıdan aşağıya, _.command.begin_ adlı bir dosyamız var. Bu boş. Bu sadece Nextflow tarafından oluşturulan, tamam, görevi başlatıyorum diyen bir sentinel dosyasıdır. Orada ilginç bir şey yok.

Sonra _.command.error_, _.command.log_ ve _.command.out_ var. Bunların hepsi çalışan bash komutundan veya bu betikten çıktılardır. Bu standart hatadır. Bu standart çıktıdır ve bu ikisinin çıktıkları gibi birleştirilmiş halidir. Yani mantıksal sırayı alırsınız.

Tamam, bunlar da bu için boştu, yani çok ilginç değil, ama _.command.run_'a geldiğinizde işler daha ilginç hale geliyor.

Bu tipik olarak çok uzun bir betiktir. Ve Nextflow'un aslında çalıştırdığı şey budur. Buraya girerseniz, Nextflow'un tüm iç mantığını görmeye başlayacaksınız ve ne yaptığını ve sürecinizi nasıl çalıştırdığını göreceksiniz. Bu, nerede çalıştığınıza bağlı olacaktır, yerel olarak mı çalıştırıyoruz yoksa SLURM'a bir iş olarak mı gönderiyoruz, bu durumda en üstte SLURM başlıkları olacak. Tüm bu farklı kurulumlar.

Genel olarak, bu dosyaya gerçekten hiç bakmanıza gerek yok. Nextflow tarafından otomatik olarak oluşturulur ve içinde boru hattınıza özgü gerçekten özel bir şey yok. Ama bu gerçekten çalışanın özüdür.

Bir sonraki çok daha ilginç. _.command.sh_, sürecinizden gelen oluşturulmuş betiktir ve burada Nextflow'un Bash başlığını eklediğini görebilirsiniz ve ardından betik bloğumuzda olan komutumuz çalıştırıldı.

Ve _.command.run_ dosyasının yaptığı tek şey bu _.command.sh_ dosyasını çalıştırmaktır.

Bu gerçekten yararlı olanıdır, bir şeylerin hatasını ayıklamaya çalışırken ve Nextflow boru hattınızın mantığının beklediğiniz şeyi yaptığını kontrol ederken genellikle en çok baktığınız olandır.

Son olarak, _.exitcode_ adlı bir dosyamız var ve bu sadece bir görevden çıkış kodunu yakalar, bu durumda başarılıydı. Yani çıkış kodu sıfırdı.

Bir şeyler ters giderse, bellek tükenirse veya başka bir şey olursa ve başarısız olursa, o zaman bu neyin yanlış gittiğini anlamak için çok yararlıdır.

## 1.3. İş akışını tekrar çalıştırın

Work dizinleri hakkında anlaşılması gereken bir şey daha var, eğer bu boru hattını tekrar tekrar çalıştırmaya devam edersem, yani _"nextflow run hello-world.nf"_ yaparsam, tamamen aynı şeyi yapacak, ama bu sefer yeni bir görev kimliği olacak. Buradaki bu hash'in farklı olduğunu görebilirsiniz ve şimdi work'e bakarsam, iki hash dizini var. Ve bunlar yine birbirinden ayrıdır.

Yani her Nextflow iş akışını çalıştırdığınızda, daha sonra değineceğimiz önbelleği kullanan resume'u kullanmadığınız sürece, bu süreçleri birbirinden ayrı olan yeni work dizinlerinde yeniden çalıştıracaktır. Herhangi bir dosya adı çakışması olmayacak, bunun gibi herhangi bir sorununuz olmayacak. Her şey izole ve temiz.

Ve bu dizine girersek, aynı dosyaların hepsini ve sıfırdan yeniden oluşturulmuş aynı _output.txt_'yi görebilirsiniz.

## 2. Çıktıları yayınlayın

Tamam, bu Nextflow'un kendisi için harika, boru hattınızı çalıştırırken böylece tüm şeyler birbirinden ayrı ve temiz ve yönetilebilir.

Ama sonuçlarınızı keşfetmeye çalışan bir kişiyseniz çok kullanışlı değil. Sonuç dosyalarınızı bulmaya çalışarak binlerce farklı work dizininde dolaşmak istemezsiniz. Ve gerçekten de öyle olması amaçlanmamıştır. Work dizinleri, dosyalarınızın oluşturulduğu son durum olması amaçlanmamıştır.

Bunu dosyalarımızı yayınlayarak yaparız.

## 2.1.1. sayHello sürecinin çıktısını bildirin

Yani betiğimize geri dönersem, burada workflow bloğumuzda çalışacağız. Hangi dosyaları beklediğini, hangi dosyaları önemsediğimizi söyleyeceğiz ve ardından altında output bloğu adında yeni bir blok oluşturacağız.

Bu, Nextflow'un sözdizimi ayrıştırıcısıyla gelen ve 26.04 sürümünde varsayılan olan yeni sözdizimdir. Yani Nextflow'u daha önce biraz kullandıysanız, bu yeni olan şeylerden biridir.

Yani main bloğuna sahibiz ve sonra publish diyeceğim ve Nextflow'a yayınlamadan ne bekleyeceğini söyleyeceğim. Buna _first_output_ diyeceğiz ve _sayHello.out_ diyeceğiz.

Orada yanlışlıkla bir yazım hatası yaptım, ama bu aynı zamanda Nextflow VS Code uzantısının bazı özelliklerini de belirtmek için iyi bir fırsat. Hemen bunun altında küçük dalgalı kırmızı bir çizgi verdiğini görebilirsiniz, bir şeylerin yanlış olduğunu söylüyor. Ve üzerine gelirseniz, bana bu değişkenin tanımlanmadığını söyleyecek. Ne olduğunu bilmiyorum.

Bu durumda oldukça açık, bir yazım hatası yaptım. sayHello yazmak istedim ve sonra dalgalı çizgi kaybolur.

Şimdi mor. Nextflow sözdizimi ayrıştırıcısı bunun bir süreç olduğunu biliyor ve üzerine geldiğimde, bana bu sürecin nasıl göründüğünün küçültülmüş bir temsilini veriyor. Yani çok hızlı bir şekilde bir bakışta herhangi bir girdi almadığını ve bize bu çıktıyı verdiğini görebilirim. Yani VS Code'da bu uzantıyla çalışmak, kod yazarken size çok fazla bağlamsal bilgi verir.

Bu süreçten çıktıya _.out_ sözdizimi ile başvurabileceğimizi unutmayın. Ve şu anda buna istediğimiz gibi adlandırabiliriz, bu sadece keyfi bir değişken adıdır.

## 2.1.2. Betiğe bir output: bloğu ekleyin

Önemli hale geldiği yer, burada yeni bloğumuzu yaptığımızda ve bu artık workflow bloğunun altında, artık workflow içinde değiliz. Tekrar süslü parantezler. Ve bu, Nextflow'a iş akışı tarafından oluşturulan tüm dosyaları nereye koyacağını söylediğimiz yerdir.

Şimdi burada oluşturduğum bu değişken adını alacağım ve onu oraya koyacağım ve bunun için bazı süslü parantezler koyacağım. Ve Nextflow'a bir path kullanmasını söyleyeceğim. Hata. Path, tırnak işaretleri içinde. Ve nokta kullanacağım. Bu sadece Nextflow'a dosyayı results dizininin kökünde koymasını söyler. Yani herhangi bir alt dizin veya başka bir şey değil.

İş akışımızı tekrar çalıştırmayı deneyelim. _"nextflow run hello-world.nf"_ yaparsam, o zaman umarım temelde tamamen aynı görünmeli. Burada Nextflow ile gerçekten hiçbir şey değişmedi. Aynı şeyleri çalıştırıyor. Sadece onları tekrar work dizinlerinde yapıyor.

Ama şimdi _"ls results/"_ yaparsam, burada results adında yeni bir dizin oluşturulduğunu göreceksiniz, bu iş akışı yayınlama için varsayılan temel dizindir. Ve içinde _output.txt_ adlı bir dosya var.

_"ls -l results"_ yaparsam, bunun aslında work dizinine soft link olduğunu göreceksiniz. Yani bu gerçek bir dosya değil, work dizinine bağlı ve bizim için oradaki tüm dosyaları topladı.

## 2.2. Özel bir konum ayarlayın

"Results" bu path için varsayılan addır. İş akışını tekrar çalıştırırsam ve bu sefer _tire_ tek tire yaparsam, bu bir çekirdek Nextflow seçeneği olduğu için. _"-output-dir **my**results"_. Kısaca _"-o"_ da yapabilirdim. O zaman dosyaların depolandığı yer için farklı bir temel dizin ayarlayacak ve bir kez daha, burada _myresults/_'ta, şimdi bir _output.txt_'miz var.

Bu harika, ama muhtemelen tüm dosyaların sadece kökte olmasını istemiyoruz. Biraz organizasyon istiyoruz, bu yüzden burada istediğimiz gibi adlandırabileceğimiz bir alt dizin de oluşturabiliriz. Diyelim ki _"path 'hello_world'"_ ve bunu tekrar çalıştırıyorum. _"nextflow run hello-world.nf"_. Results dizinine bir alt dizine gitmeli ve gerçekten de, şimdi en üstte results altında _hello_world/_ var ve _output.txt_'miz var.

Dikkat edilmesi gereken önemli şey, eski _output.txt_ dosyası hala orada. Bunu yaptığınızda results dizini silinmez. Sadece yeni dosyalar oraya kopyalanır. Aynı dosya adına sahiplerse zaten orada olan dosyaların üzerine yazarlar, ama eskileri temizlemezler. Bu yüzden boru hatlarını yeniden çalıştırdığınızda biraz dikkatli olmanız gerekir. Zaten orada olan dosyaların üzerine gelmelerini istemiyorsanız. Boş bir dizin kullandığınızdan emin olun.

## 2.3. Yayınlama modunu kopyalama olarak ayarlayın

Tamam, bu dosyaların soft link olduğunu söyledim, yani _"ls -l results/hello_world/"_ yaparsam, work dizinine soft link yaptığını görebilirsiniz. Bu genellikle HPC gibi bir şey üzerinde çalışıyorsanız iyi bir şeydir ve bunlar gerçekten büyük dosyalardır ve onları çoğaltmak istemezsiniz, çünkü bu dosyaların dosya sisteminde yalnızca bir kez depolandığı anlamına gelir.

Ancak, bu work dizinini silerseniz: _"rm -r work"_ yaparsam ve oluşturulan tüm ara dosyaları temizlersem. Şimdi, bu dosyayı okumaya çalışırsam _"results/hello_world/"_. Artık var olmayan bir dosyaya soft link olarak işaret edecek ve veriler sonsuza kadar gitti ve geri alınamaz, bu belki de pek iyi değil.

Bu yüzden genellikle, yapabiliyorsanız dosyaları soft link yerine kopyalamanın iyi bir uygulama olduğunu söylerim, çünkü daha güvenlidir. Sadece bu work dizinlerini silmediğiniz sürece iki kat daha fazla disk alanı kullanacağını unutmayın.

Bunu output bloğuyla yapmak için, burada first output'a gideceğim. Daha önce path'i ayarladım ve şimdi mode'u ayarlayacağım ve yazarken görebilirsiniz, VS code uzantısı şeyler öneriyor, bunun burada bir output yönergesi olduğunu biliyor. Ve copy diyeceğim. Kaydet'e basıyorum.

İş akışını yeniden çalıştıralım. Dosyaları tekrar oluşturacak, yeni work dizini.

Şimdi, _"ls -l results/hello_world/"_'e gidersem bunun gerçek bir dosya olduğunu ve artık soft link olmadığını görebilirsiniz ve Nextflow bunu kopyaladı. Bilmek iyi. Yani path ve mode oldukça sık yazacağınız şeylerdir.

Şimdi, elbette, bu çok basit. Bunu ilerledikçe daha karmaşık ve güçlü hale getireceğiz ve bu şeyleri nasıl dinamik yapacağınızı ve çok ayrıntılı olmayacağınızı göreceksiniz.

## 2.4. Süreç düzeyinde publishDir yönergeleri hakkında not

Şimdi, buna başlarken söylediğim gibi, bu oldukça yeni bir sözdizimi biçimidir. Bunu kaydederken sadece Nextflow'un en son sürümlerinde mevcuttur ve buna Workflow Outputs denir.

Bunu kullanırsanız, harika. Nextflow içinde, bu dosyaların oluşturulurken mirasını izlemeye yardımcı olmak için Nextflow Lineage gibi birçok başka harika özelliğin kilidini açar ve yakında 26.04'te varsayılan olacak. Ve gelecekte daha sonraki bir tarihte, bu iş akışlarınızı yazmanın tek yolu olacak.

Ancak, şu anda bu geçiş aşamasındayken, publishDir adı verilen bir şey kullanan doğada boru hatları görebilirsiniz, bu bunu yapmanın eski yoludur ve bu workflow ve output düzeyinde değil, süreç düzeyinde tanımlanır.

Ve bu bildirim temelde aynı şeyi söyler. Results adlı bir dizine sonuç dosyalarını yayınla ve bir kopyalama modu kullan diyor. Yani sözdiziminin çok benzer olduğunu görebilirsiniz. Ama şimdi yeni boru hatları yazarken, yapay zeka sonuçlarında veya dokümantasyonda veya diğer boru hatlarında görsенiz bile bu publishDir yönergesini kullanmamaya çalışın, çünkü bu bunu yapmanın eski yoludur.

2026'da hepimiz workflow outputs kullanıyor olmalıyız.

Bu tamamen belgelenmiştir, bunu yapıyorsanız ve Nextflow'u daha önce kullandıysanız, Nextflow dokümanlarına buradan gidebilirsiniz, nextflow.io/docs/. Ve eğer tutorials'a kadar aşağı kaydırırsam, _Migrating to Workflow Outputs_ adlı bir eğitim var.

Gerçekten iyidir. Tüm sözdizimini, eski sözdizime nasıl eşdeğer olduğunu, neden değiştirdiğimizi anlatır ve bir zaman çizelgesi ve her şey vardır. Ve tonlarca örnekle tüm farklı senaryoları ele alır. Böylece mevcut Nextflow kodunu kolayca yeni sözdizime dönüştürebilirsiniz.

## 3.1. sayHello sürecini değişken bir girdi bekleyecek şekilde değiştirin

Tamam, yani basit betiğimiz var, bir süreç çalıştırıyor, bir dosya oluşturuyor, Nextflow'a bunun bir çıktı olduğunu söylüyor ve sonra Nextflow'a bu dosyayı nereye kaydedeceğini söylüyoruz. Bu iyi bir başlangıç.

Ama hepsi sabit kodlanmış olmasaydı daha ilginç olurdu. Yani sonra, Nextflow'a bu sürecin, bir iş akışını başlattığımızda çalışma zamanında kontrol edebileceğimiz değişken bir girdi alabileceğini nasıl söyleyeceğimizi düşünelim.

Bunun gerçekleşmesi için birkaç farklı şey yapmamız gerekiyor.

İlk olarak, bu sürece bir girdi değişkeni kabul edebileceğini söylememiz gerekiyor ve burada yeni bir bildirim bloğu olarak _input_ yazıyoruz. Ve buna _"val greeting"_ diyeceğiz.

Val kısmı aşağıdaki path'in eşdeğeridir. Nextflow'a bunun bir değişken olduğunu, bu durumda bir string gibi olduğunu söyler. Ve üzerine gelirseniz, uzantıdan bunun ne anlama geldiğini söyler.

Sonra Nextflow'a bununla ne yapacağını söyleyeceğiz. Sadece bir değişken olduğunu söylemek yeterli değil. Betikteki değişkeni nasıl kullanacağınızı söylemeniz gerekir. Ve bu yüzden buradaki bu sabit kodlanmış string'den kurtulacağım ve bir değişken koyacağım.

Hızlıca süslü parantezler olmadan yapacağım, sadece bunun izin verildiğini göstermek için ve bu eski stil yapma yoludur. Ama şimdi yeni sözdizimi ile, bunu gerçekten böyle süslü parantezler içine koymayı öneriyoruz ve bunun Nextflow tarafından burada interpolate edildiğini çok açık hale getirir.

Harika. Yani _"input greeting"_ _$\{greeting\}_'e gider. Son şey, Nextflow'a workflow düzeyinde bu sürecin artık bir girdi aldığını söylememiz gerekiyor. Ve bunu yapmak için, temelde ona bir değişken vereceğiz.

## 3.2. Kullanıcı girdisini yakalamak için bir komut satırı parametresi ayarlayın

Bunu tekrar Hello World gibi sabit kodlayabilirdik ve bu gayet iyi çalışırdı, ama açıkçası bize gerçekten bir avantaj sağlamaz. Bunu çalışma zamanında yapılandırabilmek istedik, bu yüzden Nextflow'u başlattığınızda CLI'da yapabilmek istiyoruz.

Ve bunu yapma şekli, _params_ adı verilen özel bir Nextflow kavramıdır. Buna _params.input_ diyeceğiz.

Bunun yaptığı şey, bu input değişkenini CLI'da açığa çıkarmaktır ve Nextflow'u başlattığımızda çift tire kullandığımız yerdir.

Buna istediğim gibi adlandırabilirim, _hello, greeting_ diyebilirim. Önemli değil. Orada ne yaparsam yapayım, bir boru hattını başlattığınızda CLI seçeneği olarak açığa çıkacaktır. Ve bu Nextflow'un gerçek bir sihirli numarasıdır çünkü bu, iş akışı betiğinizi bu parametrelerle çok hızlı bir şekilde oluşturabileceğiniz anlamına gelir ve esasen boru hattınız için özel bir CLI oluşturuyorsunuz, başlattığınızda farklı seçenekleri anında özelleştirmeyi gerçekten kolaylaştırıyor.

Yani. Hadi deneyelim. Terminalimize geri dönelim. Burada _"nextflow run"_ komutumuz var. Ve şimdi daha önce gördüğümüz _"params.input"_ ile eşleşen _"--input"_ yapacağım. Sanırım dokümanlarda Fransızca. Geraldine Fransızca konuşmayı sever. Ben İsveç'te yaşadığım için İsveççe yapacağım. Bu yüzden "_Hej Världen_" diyeceğim ve enter'a basacağım.

Tek tırnak veya çift tırnak kullanabilirsiniz, bu sadece Bash'in onu nasıl yorumladığını etkiler.

Nextflow boru hattını tamamen aynı şekilde çalıştırır. Çalışma dizininin ve her şeyin aynı olduğunu görebilirsiniz. Ama şimdi _"results/hello_world/output"_'a gidersem. Burada güzel İsveççemizi görebiliriz.

Yani CLI'dan bir parametreye dinamik olarak bir girdi geçirdik. Bunu bir sürece girdi olarak geçirdik ve süreç bunu yorumladı ve bir betik bloğuna koydu, bu da o betik sonucunun çıktısını dinamik olarak değiştirdi. Oldukça havalı.

Burada çok az sözdizimi ile oldukça karmaşık mantık. Ve umarım bunun şimdi nasıl ölçeklendiğini görebilirsiniz. Ve bu gerçekten boru hatlarımızın mantığını ve özelleştirilebilirliğini Nextflow betiğine nasıl inşa ettiğimizdir.

## 3.4. Komut satırı parametreleri için varsayılan değerler kullanın

Tamam, bu harika. Ama şimdi sorun şu ki, bu boru hattını her çalıştırdığımda, çalışması için tire, input yapmam gerekiyor.

Bu parametre olmadan çalıştırmaya çalışırsam, şimdi Nextflow bu parametreye ihtiyacı olduğunu ve ayarlanmadığını söyleyerek bir hata verecek. ve bu yüzden ne yapacağını bilmedi.

Bu arada havalı yeni bir şey. Geçmişte, Nextflow sadece boş bir string ile çalışırdı ve anlaşılması zor olacak her türlü garip hatanız olurdu. Ama yeni Nextflow sözdizimi ayrıştırıcısında, biraz daha dikkatli ve size hemen söylüyor.

Yani her seçeneği her zaman belirtmek istemiyoruz. Mantıklı varsayılanları belirtmek iyi bir uygulamadır. Peki bunu betiğimizde nasıl yaparız?

Bunu yazdığımızda, _params.input_'u doğrudan kullandığımız yere koyduğumuzu fark edeceksiniz. Yani açık çözüm bir varsayılan tanımlarız ve bunu burada iş akışındaki betiğin en üstünde özel bir params bloğunda yaparız. Bu iş akışı betiğinde.

Yine, burada bazı yeni sözdizimler var, bu yüzden dikkat edin. Bu gerçekten havalı şeyler. Burada beklenecek parametrenin adını aldık.

Ve sonra bu iki nokta karakterinden sonra, değişkenin bir türünü tanımlıyoruz. Bunu yapmak zorunda değilsiniz, sadece boş bırakabilirsiniz, ama gerçekten güzel. Nextflow'a bir string beklediğimizi ve onu öyle ele almasını söyler.

Örneğin bunun yerine bir sayı istersek, float yazabiliriz ve bu bir kayan noktalı sayı istediğimizi söyler. Ve bununla çalıştırmaya çalışırsak, o zaman bir hata verecektir. Ona float olmayan bir string verirsek. Ve ayrıca onu öyle geçirecek. String yaparsak, o zaman bunun bir string olduğunu bilir. Ve önde sıfırlar olsa ve tamamen sayısal olsa bile, yine de gerçek bir string olarak geçirecektir.

Yani bu tür güvenliği Nextflow'un çok yeni bir özelliğidir, ama kodunuzu yazmayı ve çalıştırmayı daha güvenli hale getirmek için gerçekten güçlüdür.

Sonra bundan sonra bir eşittir sembolü var ve sonra burada varsayılan değer. Nextflow başlangıçta Barcelona'da yazıldı, bu yüzden burada varsayılan olarak biraz İspanyolca olması uygun görünüyor, _"Holà mundo!"_.

Tamam, o betiği kaydedeceğim, geri döneceğim, betiği _--input_ olmadan tekrar çalıştıracağım. Ve bu sefer çalışmalı ve _results_'ta yeni dosyamızı oluşturacak. Ve bu dosyada şimdi _"Holà mundo!"_ diyor.

Bu sadece bir varsayılan olsa da, hala daha önce yaptığımız aynı şeyi yapamayacağımız anlamına gelmez. Geri dönüp burada eski betiğimi bulursam, _"Hej Världen"_, çünkü komut satırında _--input_ yapıyorum, bu varsayılanın üzerine yazacak ve onu output.txt dosyasında tekrar kullanacak.

Yani betikteki bu sadece ayarladığım varsayılan değerdir.

İş akışımızı daha karmaşık hale getirip daha fazla parametre içerecek şekilde oluşturdukça, betiğin en üstündeki bu params bloğu hepsini tek bir yerde toplamaya başlayacak.

Ve betiğinizde bu oldukça güzel bir simetri elde edersiniz, burada etkili bir şekilde tüm iş akışı girdileriniz var ve altta iş akışı çıktılarınız. Ve iş akışınızın dış dünyaya arayüzünün ne olduğu çok açık. Böylece yeni bir boru hattını yeni sözdizimi ile çok hızlı bir şekilde alabilir ve nasıl kullanılacağını anlayabilirsiniz.

Son bir havalı şey. Bununla bir varsayılan değer ayarlamak zorunda değiliz. Params input yaparsak ama bir varsayılan değer ayarlamazsak, o zaman Nextflow'a bu parametrenin gerekli olduğunu söyler ve yine, boru hattı onsuz çalışmayı başaramaz, ama size null olması hakkında bir şey yerine daha kullanışlı bir hata mesajı verecektir.

Yani bir input beklediğimizi söylüyor, gereklidir, ama komut satırında belirtilmedi. Çok güzel.

Tamam, umarım şimdi Nextflow boru hattınızı değişken girdiler ve parametrelerle nasıl kuracağınız, varsayılanı nasıl ayarlayacağınız, türleri nasıl ayarlayacağınız, bir Boolean doğru yanlış bayrağı veya bir tamsayı veya burada farklı türler olabilir konusunda açıktır. Bunları iş akışınıza nasıl geçireceğiniz, nereden geçtiği ve sonra sürecinize interpolate olur. Ve ayrıca Nextflow'u başlattığınızda komut satırında bunları nasıl özelleştireceğinizi biliyorsunuz. Bu basit bash komutumuzdan daha ilginç görünmeye başlıyor.

## 4. İş akışı yürütmelerini yönetin

Tamam. Sırada ne var? Bu bölümün son kısmı için, tüm farklı iş akışı yürütmelerini nasıl yöneteceğimiz hakkında biraz konuşacağız. Kenar çubuğuma ve work altındaki Explorer'a bakarsanız, bir sürü farklı boru hattı çalıştırdığımı göreceksiniz ve bu work dizinleri oldukça uzuyor, çok sayıda var.

Ve diğer şey, daha önce söylediğim gibi, bu boru hattını her yeniden çalıştırdığımda, yeni bir work dizinleri seti oluşturuyor ve tüm süreçleri sıfırdan yeniden çalıştırıyor, bu iyi bir şey. Bu amaçlanan davranıştır. Tekrarlanabilir ve her şeyi taze bir şekilde yeniden üretiyor. Ama açıkçası, çok uzun süren süreçler çalıştırıyorsanız, yarı yolda çöktüyse veya boru hattının sonunda bir şeyi değiştirdiyseniz, boru hattınızı her zaman baştan başlatmak zorunda kalmak sinir bozucu.

## 4.1. Bir iş akışını -resume ile yeniden başlatın

Neyse ki, Nextflow daha önce neyin çalıştırıldığını ve neyin mevcut olduğunu bilmekte gerçekten iyidir ve bu eski sonuçları yeniden kullanmak çok basittir. Komutun sonuna sadece yeni bir bayrak ekliyoruz _"-resume"_.

Şimdi, input'ta iki tire var çünkü bu parametre. Resume'da sadece bir tire var çünkü bu bir çekirdek Nextflow seçeneği.

İnsanları her zaman tökezletiyor, Nextflow'u uzun süredir kullanıyor olsanız bile. Bu yüzden her zaman bir veya iki tire hatırlayın. Çekirdek bir Nextflow seçeneği olup olmadığına bağlıdır.

Tamam, şimdi _-resume_ yapıyorum ve tamamen aynı iş akışını tekrar çalıştırıyorum. Ve bu sefer bir temel farkla hemen hemen tamamen aynı görünmeli.

Buradaki çıktıda, sonuçların önbelleğe alındığını görebilirsiniz. Ve aslında, buradaki bu görev hash'i önceki çalıştırmayla tamamen aynı ve sadece o work dizinini bütünüyle yeniden kullandı. Girdiler ve çıktılar ve betik hepsi değiştirilmemişti. Ve bu yüzden sadece o dosyayı oradan aldı ve boru hattında aşağı akış adımları olsaydı, onları bir sonraki adıma geçirirdi.

Yani hala tüm boru hattını baştan sona çalıştırıyor, ama yapabildiği her bir görev için önbelleğe alınmış sonuçları kullanıyor.

Şimdi, _-resume_ yaptığınızda, çalışma dizininizdeki son boru hattı çalıştırmasını, her ne olursa olsun devam ettirir. Ama aslında orada yaptığınız herhangi bir önceki çalıştırmadan devam edebilirsiniz. Ve şimdi oldukça fazla yaptık.

## 4.2. Geçmiş yürütmelerin günlüğünü inceleyin

Hepsine bakmak için, _"nextflow run"_ yerine _"nextflow log"_ yapabiliriz ve bu bize tüm bu farklı.. Görebilmemiz için ekranımı biraz küçültmem gerekiyor, tüm bu farklı çalıştırmaları ne zaman yaptığımızı, oturum kimliğini, komutu ve her şeyi gösteren güzel bir çıktı verecek.

Ve buraya bakabiliriz ve bunlardan herhangi birinin çalıştırma adını alabiliriz ve sonra o belirli olanlardan birini devam ettirebiliriz. Yani geri dönebilirim ve _hungry_ekeblad_ adlı olanı devam ettirebilirim. Ve bunu resume'dan sonra koyuyorum.

Bu arada merak ediyorsanız, tüm bu sıfatlar ve bilim insanları isimleri Nextflow kaynak kodundadır. Gidip bularak ve favori bilim insanınızı ekleyerek Nextflow'a ilk pull request'inizi almanın gerçekten iyi bir yolu.

Ve neyse, bunu yaptım ve geri gitti ve bu iş akışı çalıştırmasından önbelleğe alınmış sonuçlara baktı, hala onları yeniden kullanabileceğini fark etti ve yaptı. Yani önbelleğe alınmış sonuçları tekrar aldım.

## 4.3. Eski work dizinlerini silin

Bu harika. Peki ya bu work dizinlerini temizlemek istersem? Burada tonlarca var. Tonlarca dosya var. Belki kesin olarak biliyorum ki son birkaç boru hattı çalıştırmasından devam etmek istiyorum, ama öncekilerle ilgilenmiyorum.

O zaman burada birini seçebilirim ve başka bir Nextflow komutu kullanabilirim, bu _"nextflow clean"_ ve _"nextflow clean"_ yapabilirim, _"-before"_ yapacağım ve belirli çalıştırma adı, bu durumda _reverent_pike_ ve _"-n"_ yapacağım, bu Nextflow'a sadece kuru bir çalıştırma yapmasını söyler. Yani bana sadece neyi sileceğini söyler. Aslında hiçbir şey yapmadan, bu work dizinlerini kaldırırdı.

Bu mantıklı görünüyor. Bu yüzden aynı komutu tekrar yapacağım, ama _"-n"_ yerine temizliği gerçekten yapmak için _"-f"_ yapacağım. Ve bu sefer gerçekten tüm bu dizinleri kaldırdı. Ve girip work dizinlerine bakarsam, şimdi çok daha hafif görünüyor. Harika.

Yani tüm yerel work dizinlerinizi önbelleği tamamen yok etmeden oldukça güvenli bir şekilde nasıl temizleyeceğiniz budur. Böylece isterseniz hala devam edebilirsiniz.

Bu bayrakların ne için olduğunu unutursanız, her Nextflow komutu için _"nextflow help"_ yapabilirsiniz ve ardından komutun adını yazabilirsiniz. Yani _"nextflow help clean"_ yaparsam, tüm farklı seçenekleri görebilirsiniz: _-after, -before, -but_, bu temizleme davranışını yapılandırmanın tüm farklı yolları. Oldukça havalı.

## Özet

Tamam, bu Hello Nextflow'un birinci bölümünün sonu. Kursun oldukça yoğun bir başlangıcı, ama umarım şimdi bir Nextflow betiğinin nasıl göründüğü hakkında oldukça iyi bir anlayışa sahipsiniz; farklı temel parçalarla, süreçler, iş akışları, çıktılar ve parametreler. Bunları komut satırından temel geçersiz kılmalarla nasıl yapılandıracağınızı, dinamik bir betikle dinamik bir girdi bloğu nasıl yapacağınızı biliyorsunuz ve tüm iş yükü yürütmelerinizi nasıl yöneteceğinizi biliyorsunuz: zaten ne çalıştırdığınızı görmek, devam ettirmek, temizlemek. Çok fazla şey var. Uzun bir yol kat ettiniz. Bu yüzden bir mola vermek ve hızlı bir yürüyüş ve bir fincan çay içmek istiyorsanız, şimdi muhtemelen iyi bir zaman. Hak ettiniz.

Buradan itibaren, temelde bu temel üzerine inşa ediyoruz. Bunu nasıl daha karmaşık, daha güçlü hale getirebiliriz? Nasıl daha esnek hale getirebiliriz? Analizimizi ölçekte yapmak istediğimiz şeyleri yapın.

## Quiz

Şimdi web sayfasında birinci bölüm, hello world'e kadar aşağı kaydırırsanız küçük bir quiz göreceksiniz ve bu, Nextflow eğitiminin bu sürümü için yaptığımız yeni bir şey. Ve gidip bu bölümde yaptığımız tüm materyali anladığınızı kontrol etmek için kendinizi sınayabilirsiniz.

Bu bize gönderilmiyor veya başka bir şey, sadece tarayıcınızda saklanıyor. Yani cevaplarınızın ne olduğunu bilmiyoruz, ama bu sadece bir şeyi kaçırmadığınızdan veya yanlış anlamadığınızdan emin olmak için küçük bir kendi kendine kontrol. Ve istediğiniz kadar deneyebilirsiniz.

Eğer benim gibiyseniz, belki VS Code örneğinizde terminalde kalmak istiyorsunuz, bu durumda _quiz_ komutunu yazabilirsiniz ve sonra sadece hangi bölümde olduğunuzu söyleyebilirsiniz. Yani _"Hello World"_ yapıyoruz ve sonra web tarayıcısında olan tamamen aynı quiz sorularını yapabilirsiniz, ama sadece terminalinizde.

Havalı. Tamam. Umarım bundan keyif alırsınız. Biraz eğlenin ve bir sonraki bölümde sadece bir dakika içinde Nextflow kanalları hakkında konuşmak için görüşürüz.
