# Bölüm 1: Hello World - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../01_hello_world.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba ve tekrar hoş geldiniz.

Şu anda "Hello Nextflow" kursunun "Hello World" adlı Birinci Bölümündesiniz. Bu bölümde, Nextflow'un en temel kavramlarını anlamaya başlayacağız.

Umarım artık Codespaces'te veya eşdeğer bir ortamda VS Code çalışır durumdadır ve Explorer'da çalışma alanınızda Hello Nextflow klasörünüz ve buradaki tüm farklı dosyalarınız bulunmaktadır.

Terminalde Bash kullanarak çok temel şeyler yaparak başlayacağız, ardından aynı şeyleri Nextflow içinde yapıp yapamayacağımıza bakacağız, böylece söz diziminin nasıl göründüğüne dair bir fikir edineceksiniz.

## 0. Isınma

Hadi gerçekten basit bir şeyle başlayalım. Sadece "echo" ile başlayalım, terminale bir şeyler yazdırmak için. "Hello World". Enter'a basıyorum ve bu terminale gidiyor. Hello World. Umarım bu kursu izleyen hiç kimse için sürpriz değildir.

Tamam, bununla bir şeyler yapalım. Bunu sadece terminale yazdırmak yerine, bir dosyaya yazalım. Klavyemdeki yukarı ok tuşuna basacağım, bu Bash geçmişinde gezinir, böylece bana son komutumu verir ve sonuna küçük bir büyüktür sembolü ekleyeceğim, bu da bu komutun çıktısını bir dosyaya yönlendirir ve buna output.txt adını vereceğim.

Tekrar Enter, bu komutu çalıştırmak için, bu sefer terminalde hiçbir şey yok, ama sol tarafta yeni dosyanın burada göründüğünü görebiliriz, output.txt adında.

Bunu terminalde cat gibi bir şeyle görüntüleyebiliriz. Yani cat output.txt ve gerçekten de "Hello World" diyor. Ayrıca çift tıklayabiliriz ve VS Code'daki kod düzenleyicide açılır.

## 1.1. Kodu İnceleyin

Pekala. Size basit olduğunu söylemiştim. Sırada ne var? Hadi bu süreci alıp tekrar yapalım, ama bu sefer, bunu Nextflow içinde yapalım.

Dediğim gibi, bu kurstaki farklı bölümlerin hepsi bir betikle başlar ve bu bölümün adı Hello World. Bu yüzden Hello World'ü bulacağım. Tek tıkladığımda önizleme yapıyor, burada düzenleyicide açmak için çift tıklayacağım. Ve terminali hızlıca kaldıracağım.

Şimdi bu çok basit bir betik, yani olabildiğince basit. Sadece 22 satır uzunluğunda ve temelde aynı şeyi yapıyor. Aslında. Bunun bir kısmı tanıdık gelmeli. Az önce yazdığımız bash komutu bu. Orada bir dosyaya yönlendirme yapan bash komutunu görebiliriz.

Tamam. Başka ne var? Ayrıca, bu dosyada Nextflow'un bazı temel kavramlarını görmeye başlayabiliriz. Burada kırmızıyla bir process ve bir workflow var. Bunlar Nextflow'da özel anahtar kelimeler ve özel terminolojidir.

## 1.1.1. Süreç tanımı

Bir iş akışı içindeki farklı süreçler, iş akışınızın farklı mantıksal birimlerini sarar. Her süreç bir şey yapar.

Çalıştırdığımızda, bir görev veya birden fazla görev oluşturur, bunlar bir boru hattının gerçek yapma adımlarıdır. Tüm süreçler daha sonra altta gördüğümüz workflow bloğu içinde düzenlenir ve bu durumda sadece o tek süreci çalıştırır.

Süreç adı buradaki bu anahtar kelimeyi takip eder ve bu temelde herhangi bir şey olabilir. Ve sonra sürecin içeriği bu süslü parantezler içindedir.

Bir süreç için gerçekten sadece bir gereklilik vardır, o da bir tür script veya exec bloğu içermesidir. Bu burada üçlü tırnak içindedir ve bu, boru hattını çalıştırdığımızda çalışma dizinine yazılan bash betiğidir ve bilgisayarınızda veya sunucunuzda gerçekten çalışan şey budur.

Bu genellikle bash'tir, ancak burada en üstte farklı bir shebang da koyabilirsiniz ve bu bir Python betiği veya R betiği olabilir. Önemli değil. Bu betikteki her ne varsa çalıştırılacaktır.

Bu sürece eklediğimiz bir şey daha var, o da output bildirimi. Bu, Nextflow'a bu sürecin output.txt adlı bir çıktı dosyası beklediğini söyler. Path olduğunu söylüyor, yani bir dosya gibi işlenmeli, mesela bu val olsaydı, bir değişken veya değer gibi olduğunu söylerdi.

Bunun bu dosyayı oluşturmadığını unutmayın. Aslında onu üretmiyor. Bu aşağıdaki betik tarafından yapılıyor. Sadece Nextflow'a bu dosya adıyla bir çıktı dosyası beklemesini söylüyor.

## 1.1.2. İş akışı tanımı

Tamam. Ve sonra altta bir workflow var ve yine bir bildirim var. Bunun adı Main. Bu, isterseniz bir script bloğunun workflow eşdeğeridir. İş akışının bir şey yapan kısmıdır. Ve bu durumda, sayHello adlı süreci çağır diyoruz.

Normalde, elbette, boru hattınız bundan çok daha karmaşık görünecektir. Muhtemelen birden fazla süreciniz olacak ve aralarındaki veri akışını düzenlemek için kanallar kullanacaksınız. Buna bu kursun sonraki bölümlerinde geleceğiz, ama şimdilik bu yeterli. Bu çalışması gereken geçerli bir boru hattıdır.

VS Code'da burada DAG önizlemesine bile tıklayabilirim. DAG veya DAG, boru hattındaki bir veri akışı yapısının temsilidir ve onu yan tarafta bir mermaid diyagramı olarak işlenmiş görebiliriz. Bu durumda çok basit. Bir kutu var, bu workflow ve bir süreç var, adı sayHello, ama ilerledikçe bu daha ilginç görünebilir.

## 1.2. İş akışını çalıştırın

Tamam, hadi bu iş akışını çalıştırmayı deneyelim ve ne olduğunu görelim.

Terminali tekrar altta getireceğim, çıktıyı temizleyeceğim ve Nextflow Run yazacağım. Ve sonra sadece betik adını yazacağım, bu hello-world.nf. Ve Enter'a basacağım.

Tamam, üstte bize Nextflow'un çalıştığını ve hangi sürümün çalıştığını ve betik adının ne olduğunu ve her şeyi söyleyen bazı standart şeyler var.

Ve gerçekten aradığımız önemli şey _burada_, bu çalıştırılan farklı görevlerin bir özetidir.

Eğer sizinki böyle küçük yeşil bir tik işaretiyle görünüyorsa, tebrikler. İlk boru hattınızı çalıştırdınız. Harika.

Burada bize çalışan sürecin adını söylüyor, adı Say Hello ve bir kez çalıştığını ve başarılı olduğunu söylüyor. Bu ilerledikçe güncellenir, bu yüzden daha büyük bir boru hattı çalıştırdığınızda, ilerlemeyi burada temsil edildiğini göreceksiniz. Ama bu çok küçük olduğu için, temelde anında çalışır.

## 1.2.2. Çıktıyı ve günlükleri work dizininde bulun

Şimdi bir Nextflow boru hattı çalıştırdığınızda, bu süreçlerin her biri birbirine bağlanır ve her süreç, daha önce söylediğim gibi, görevler oluşturabilir, bir veya birden fazla. Yani bu durumda, bu süreçten tek bir görevimiz vardı. Sadece bir kez çalıştı ve bu, bu görev _hash_'i altında yapıldı.

Nextflow, çalışma dizininizdeki dosyalarla doğrudan ilgilenmez, work adında özel bir klasör oluşturur. Ve "ls" yaparsam, burada göründüğünü göreceğiz: _work_, ve bunun içinde çalışan her görev için alt dizinler var. Ve bu, bu hash ile eşleşir. Yani "ls work/c4" yaparsam görebilirsiniz, ve sonra kısaltılmış, ama 203 ile başlıyor ve bu, boru hattını çalıştırdığımızda bu süreç tarafından oluşturulan çalışma dizinidir. Ve onu yan tarafta da görebilirsiniz.

Bu dosyaları listelediğimde, output.txt dosyasının oluşturulduğunu görebilirsiniz. Bunu burada da görebilirsiniz. Ve normal "ls" ile gösterilmeyen bir sürü gizli dosya var.

output.txt'ye tıklarsam, gerçekten de çıktımız var. Harika. Yani boru hattı çalıştı.

Esasen tek satırlık bir bash betiği çalıştırmak için oldukça fazla şablon gibi görünebilir, ama süreçlerimiz daha karmaşık hale geldikçe daha mantıklı olacak. Ve Nextflow ile bu work dizini ve oluşturulan bu dosyalar, gerçekten Nextflow'u bu kadar güçlü yapan şeyin omurgasıdır.

Her görev, bir boru hattının her öğesi diğer tüm görevlerden izole edilmiştir. Tekrarlanabilir. Birbirleriyle çakışmazlar ve her şey paralel çalışabilir. Aslında alıştığınızda gerçekten güzel bir yoldur çünkü bu izolasyon sayesinde içeri girip tek bir görev için tam olarak ne olduğunu görebilir ve hata ayıklayabilirsiniz.

Work dizinindeki bu diğer dosyalara hızlıca bir göz atalım. Yukarıdan aşağıya, _.command.begin_ adlı bir dosyamız var. Bu boş. Bu sadece Nextflow tarafından oluşturulan, tamam, görevi başlatıyorum diyen bir sentinel dosyasıdır. Orada ilginç bir şey yok.

Sonra _.command.error_, _.command.log_ ve _.command.out_ var. Bunların hepsi çalışan bash komutundan veya bu betikten çıktılardır. Bu standart hata. Bu standart çıktı ve bu ikisinin çıktıkları gibi birleştirilmiş halidir. Yani mantıksal sırayı alırsınız.

Tamam, bunlar da bu için boştu, yani çok ilginç değil, ama _.command.run_'a geldiğinizde işler daha ilginç hale geliyor.

Bu genellikle çok uzun bir betiktir. Ve Nextflow'un aslında çalıştırdığı şey budur. Buraya girerseniz, Nextflow'un tüm iç mantığını görmeye başlarsınız ve ne yaptığını ve sürecinizi nasıl çalıştırdığını görürsünüz. Bu, nerede çalıştırdığınıza bağlı olacaktır, yerel olarak mı çalıştırıyoruz yoksa SLURM'a bir iş olarak mı gönderiyoruz, bu durumda üstte SLURM başlıkları olacaktır. Tüm bu farklı kurulumlar.

Genel olarak, bu dosyaya gerçekten hiç bakmanıza gerek yok. Nextflow tarafından otomatik olarak oluşturulur ve içinde boru hattınıza özgü gerçekten özel bir şey yok. Ama bu gerçekten çalışanın özüdür.

Bir sonraki çok daha ilginç. _.command.sh_, sürecinizden gelen oluşturulmuş betiktir ve burada Nextflow'un Bash başlığını eklediğini görebilirsiniz ve sonra script bloğumuzda olan komutumuz çalıştırıldı.

Ve _.command.run_ dosyasının yaptığı tek şey bu _.command.sh_ dosyasını çalıştırmaktır.

Bu gerçekten yararlı olanıdır, bir şeylerin hatasını ayıklamaya çalışırken ve Nextflow boru hattınızın mantığının beklediğiniz şeyi yaptığını kontrol ederken genellikle en çok baktığınız olandır.

Son olarak, _.exitcode_ adlı bir dosyamız var ve bu sadece bir görevden çıkış kodunu yakalar, bu durumda başarılıydı. Yani çıkış kodu sıfırdı.

Bir şeyler ters giderse, bellek tükenirse veya başka bir şey olur ve başarısız olursa, o zaman bu neyin yanlış gittiğini anlamak için çok yararlıdır.

## 1.3. İş akışını tekrar çalıştırın

Work dizinleri hakkında anlaşılması gereken bir şey daha var, eğer bu boru hattını tekrar tekrar çalıştırmaya devam edersem, yani _"nextflow run hello-world.nf"_ yaparsam, tam olarak aynı şeyi yapacak, ama bu sefer yeni bir görev id'si olacak. Buradaki bu hash'in farklı olduğunu görebilirsiniz ve şimdi work'e bakarsam, iki hash dizini var. Ve bunlar yine birbirinden ayrıdır.

Yani bir Nextflow iş akışını her çalıştırdığınızda, önbelleği kullanan resume'u kullanmadığınız sürece, daha sonra değineceğiz, bu süreçleri birbirinden ayrı olan yeni work dizinlerinde yeniden çalıştıracaktır. Herhangi bir dosya adı çakışması olmayacak, bunun gibi herhangi bir sorununuz olmayacak. Her şey izole ve temiz.

Ve bu dizine girersek, tüm aynı dosyaları ve sıfırdan yeniden oluşturulmuş aynı _output.txt_'yi görebilirsiniz.

## 2. Çıktıları yayınlayın

Tamam, bu Nextflow'un kendisi için harika, boru hattınızı çalıştırırken böylece tüm şeyler birbirinden ayrı ve temiz ve yönetilebilir.

Ama sonuçlarınızı keşfetmeye çalışan bir kişiyseniz çok kullanışlı değil. Sonuç dosyalarınızı bulmaya çalışarak binlerce farklı work dizininde dolaşmak istemezsiniz. Ve gerçekten de öyle olması amaçlanmamıştır. Work dizinleri, dosyalarınızın oluşturulduğu son durum olması amaçlanmamıştır.

Bunu dosyalarımızı yayınlayarak yaparız.

## 2.1.1. sayHello sürecinin çıktısını bildirin

Yani betiğimize geri dönersem, burada workflow bloğumuzda çalışacağız. Hangi dosyaları beklediğini, hangi dosyaları önemsediğimizi söyleyeceğiz ve sonra altında output bloğu adında yeni bir blok oluşturacağız.

Bu, Nextflow'un söz dizimi ayrıştırıcısıyla gelen ve 26.04 sürümünde varsayılan olan yeni söz dizimidir. Yani Nextflow'u daha önce biraz kullandıysanız, bu yeni olan şeylerden biridir.

Yani main bloğumuz var ve sonra publish diyeceğim ve Nextflow'a yayınlamadan ne bekleyeceğini söyleyeceğim. Buna _first_output_ diyeceğiz ve _sayHello.out_ diyeceğiz.

Orada yanlışlıkla bir yazım hatası yaptım, ama bu aynı zamanda Nextflow VS Code uzantısının bazı özelliklerini de göstermek için iyi bir fırsat. Hemen bunun altında küçük dalgalı kırmızı bir çizgi verdiğini görebilirsiniz, bir şeylerin yanlış olduğunu söylüyor. Ve üzerine gelirseniz, bu değişken tanımlı değil diyecektir. Ne olduğunu bilmiyorum.

Bu durumda oldukça açık, bir yazım hatası yaptım. sayHello yazmak istedim ve sonra dalgalı çizgi kaybolur.

Şimdi mor. Nextflow söz dizimi ayrıştırıcısı bunun bir süreç olduğunu biliyor ve üzerine geldiğimde, bana bu sürecin nasıl göründüğünün küçültülmüş bir temsilini veriyor. Yani çok hızlı bir bakışta herhangi bir girdi almadığını ve bize bu çıktıyı verdiğini görebilirim. Yani VS Code'da bu uzantıyla çalışmak, kod yazarken size çok fazla bağlamsal bilgi verir.

Bu süreçten çıktıya _.out_ söz dizimiyle başvurabileceğimizi unutmayın. Ve şu anda buna istediğimiz gibi adlandırabiliriz, bu sadece keyfi bir değişken adıdır.

## 2.1.2. Betiğe bir output: bloğu ekleyin

Önemli hale geldiği yer, burada yeni bloğumuzu yaptığımızda ve bu artık workflow bloğunun altında, artık workflow içinde değiliz. Tekrar süslü parantezler. Ve bu, Nextflow'a iş akışı tarafından oluşturulan tüm dosyaları nereye koyacağını söylediğimiz yerdir.

Şimdi burada oluşturduğum bu değişken adını alacağım ve onu oraya koyacağım ve bunun için bazı süslü parantezler koyacağım. Ve Nextflow'a bir path kullanmasını söyleyeceğim. Hata. Path, tırnak içinde. Ve nokta kullanacağım. Bu sadece Nextflow'a dosyayı results dizininin kökünde koymasını söyler. Yani herhangi bir alt dizin veya başka bir şey değil.

İş akışımızı tekrar çalıştırmayı deneyelim. _"nextflow run hello-world.nf"_ yaparsam, o zaman umarım temelde tamamen aynı görünmeli. Burada Nextflow ile gerçekten hiçbir şey değişmedi. Aynı şeyleri çalıştırıyor. Sadece onları tekrar work dizinlerinde yapıyor.

Ama şimdi _"ls results/"_ yaparsam, burada results adında yeni bir dizin oluşturulduğunu göreceksiniz, bu workflow yayınlama için varsayılan temel dizindir. Ve içinde _output.txt_ adlı bir dosya var.

_"ls -l results"_ yaparsam, bunun aslında work dizinine soft link olduğunu göreceksiniz. Yani bu gerçek bir dosya değil, work dizinine bağlı ve bizim için oradaki tüm dosyaları topladı.

## 2.2. Özel bir konum ayarlayın

"Results", bu path için varsayılan addır. İş akışını tekrar çalıştırırsam ve bu sefer _tire_ tek tire yaparsam, bu bir çekirdek Nextflow seçeneği olduğu için. _"-output-dir **my**results"_. Kısa olarak sadece _"-o"_ da yapabilirdim. O zaman dosyaların depolandığı yer için farklı bir temel dizin ayarlayacak ve bir kez daha, burada _myresults/_ içinde, şimdi bir _output.txt_ var.

Bu harika, ama muhtemelen tüm dosyaların sadece kök dizinde olmasını istemiyoruz. Biraz organizasyon istiyoruz, bu yüzden burada istediğimiz gibi adlandırabileceğimiz bir alt dizin de oluşturabiliriz. Diyelim ki _"path 'hello_world'"_ ve bunu tekrar çalıştırıyorum. _"nextflow run hello-world.nf"_. Results dizinine bir alt dizine gitmeli ve gerçekten de, şimdi en üstte results altında _hello_world/_ var ve _output.txt_ var.

Dikkat edilmesi gereken önemli şey, eski _output.txt_ dosyası hala orada. Bunu yaptığınızda results dizini silinmez. Sadece yeni dosyalar oraya kopyalanır. Aynı dosya adına sahiplerse zaten orada olan dosyaların üzerine yazarlar, ama eskileri temizlemezler. Bu yüzden boru hatlarını yeniden çalıştırdığınızda biraz dikkatli olmanız gerekir. Zaten orada olan dosyaların üzerine gelmelerini istemiyorsanız. Boş bir dizin kullandığınızdan emin olun.

## 2.3. Yayınlama modunu kopyalama olarak ayarlayın

Tamam, bu dosyaların soft link olduğunu söylemiştim, yani _"ls -l results/hello_world/"_ yaparsam, work dizinine soft link yaptığını görebilirsiniz. Bu genellikle HPC gibi bir şey üzerinde çalışıyorsanız iyi bir şeydir ve bunlar gerçekten büyük dosyalardır ve onları çoğaltmak istemezsiniz, çünkü bu dosyaların dosya sisteminde yalnızca bir kez depolandığı anlamına gelir.

Ancak, work dizinini silerseniz: _"rm -r work"_ yaparsam ve oluşturulan tüm ara dosyaları temizlersem. Şimdi, bu dosyayı okumaya çalışırsam _"results/hello_world/"_. Artık var olmayan bir dosyaya soft link olarak işaret edecek ve veriler sonsuza kadar gitti ve geri alınamaz, bu belki de pek iyi değil.

Bu yüzden genellikle, yapabiliyorsanız dosyaları soft link yerine kopyalamanın iyi bir uygulama olduğunu söylerim, çünkü daha güvenlidir. Sadece iki kat daha fazla disk alanı kullanacağının farkında olun, bu work dizinlerini silmediğiniz sürece.

Bunu output bloğuyla yapmak için, burada first output'a gideceğim. Daha önce path'i ayarladım ve şimdi mode'u ayarlayacağım ve yazdıkça VS code uzantısının önerdiği şeyleri görebilirsiniz, burada bir output yönergesi olduğunu biliyor. Ve copy diyeceğim. Kaydet'e basıyorum.

İş akışını yeniden çalıştıralım. Dosyaları tekrar oluşturacak, yeni work dizini.

Şimdi, _"ls -l results/hello_world/"_'e gidersem bunun gerçek bir dosya olduğunu görebilirsiniz ve artık soft link değil ve Nextflow bunu kopyaladı. Bilmek iyi. Yani path ve mode oldukça sık yazacağınız şeylerdir.

Şimdi, elbette, bu çok basit. Bunu daha karmaşık ve güçlü hale getireceğiz ve ilerledikçe bu şeyleri nasıl dinamik yapacağınızı ve çok ayrıntılı olmayacağınızı göreceksiniz.

## 2.4. Süreç düzeyinde publishDir yönergeleri hakkında not

Şimdi, buna başlarken söylediğim gibi, bu oldukça yeni bir söz dizimi biçimidir. Bunu kaydederken yalnızca Nextflow'un en son sürümlerinde mevcuttur ve buna Workflow Outputs denir.

Bunu kullanırsanız, harika. Nextflow içinde başka birçok harika özelliğin kilidini açar, örneğin, bu dosyalar oluşturulurken miraslarını izlemeye yardımcı olmak için Nextflow Lineage ve yakında 26.04'te varsayılan olacak. Ve gelecekte daha sonraki bir tarihte, bu iş akışlarınızı yazmanın tek yolu olacak.

Ancak, şu anda bu geçiş aşamasındayken, publishDir adında bir şey kullanan doğada boru hatları görebilirsiniz, bu bunu yapmanın eski yoludur ve bu workflow ve output düzeyinde değil, süreç düzeyinde tanımlanır.

Ve bu bildirim temelde aynı şeyi söyler. Sonuç dosyalarını results adlı bir dizine yayınla ve bir copy modu kullan diyor. Yani söz diziminin çok benzer olduğunu görebilirsiniz. Ama şimdi yeni boru hatları yazarken, bu publishDir yönergesini kullanmamaya çalışın, AI sonuçlarında veya dokümantasyonda veya diğer boru hatlarında görsенiz bile, çünkü bu bunu yapmanın eski yoludur.

2026'da hepimiz workflow outputs kullanıyor olmalıyız.

Bunların hepsi belgelenmiştir, bunu yapıyorsanız ve daha önce Nextflow kullandıysanız, burada Nextflow dokümanlarına gidebilirsiniz, nextflow.io/docs/. Ve aşağı kaydırırsam tutorials'a, _Migrating to Workflow Outputs_ adlı bir eğitim var.

Gerçekten iyidir. Tüm söz dizimini, eski söz dizimine nasıl eşdeğer olduğunu, neden değiştirdiğimizi ve bir zaman çizelgesi ve her şeyi anlatır. Ve tonlarca örnekle tüm farklı senaryoları anlatır. Böylece mevcut Nextflow kodunu kolayca yeni söz dizimine dönüştürebilirsiniz.

## 3.1. sayHello sürecini değişken bir girdi bekleyecek şekilde değiştirin

Tamam, yani bir süreç çalıştıran, bir dosya oluşturan, Nextflow'a bunun bir çıktı olduğunu söyleyen ve sonra Nextflow'a bu dosyayı nereye kaydedeceğini söyleyen basit betiğimiz var. Bu iyi bir başlangıç.

Ama hepsi sabit kodlanmış olmasaydı daha ilginç olurdu. Yani sonra, bu sürece bir iş akışını başlattığımızda çalışma zamanında kontrol edebileceğimiz değişken bir girdi alabileceğini Nextflow'a nasıl söyleyeceğimizi düşünelim.

Bunun gerçekleşmesi için birkaç farklı şey yapmamız gerekiyor.

İlk olarak, bu sürece bir girdi değişkeni kabul edebileceğini söylememiz gerekiyor ve burada yeni bir bildirim bloğu olarak _input_ yazıyoruz. Ve buna _"val greeting"_ diyeceğiz.

val kısmı, aşağıdaki path'in eşdeğeridir. Nextflow'a bunun bir değişken olduğunu söyler, bu durumda bir string gibi. Ve üzerine gelirseniz, uzantıdan bunun ne anlama geldiğini söyler.

Sonra Nextflow'a bununla ne yapacağını söyleyeceğiz. Sadece bir değişken olduğunu söylemek yeterli değil. Betikteki değişkeni nasıl kullanacağınızı söylemeniz gerekir. Ve bu yüzden buradaki bu sabit kodlanmış string'i kaldıracağım ve bir değişken koyacağım.

Bunu süslü parantezler olmadan hızlıca yapacağım, sadece bunun izin verildiğini göstermek için ve bu bunu yapmanın eski tarz yoludur. Ama şimdi yeni söz dizimiyle, bunu gerçekten böyle süslü parantezler içine koymayı öneriyoruz ve bunun burada Nextflow tarafından interpolate edildiğini çok açık hale getirir.

Harika. Yani _"input greeting"_ _$\{greeting\}_ içine gider. Son şey, Nextflow'a workflow düzeyinde bu sürecin artık bir girdi aldığını söylememiz gerekiyor. Ve bunu yapmak için, temelde ona bir değişken vereceğiz.

## 3.2. Kullanıcı girdisini yakalamak için bir komut satırı parametresi ayarlayın

Bunu tekrar Hello World gibi sabit kodlayabilirdik ve bu gayet iyi çalışırdı, ama açıkçası bize gerçekten bir avantaj sağlamaz. Bunu çalışma zamanında yapılandırabilmek istedik, bu yüzden Nextflow'u başlattığınızda CLI'da yapabilmek istiyoruz.

Ve bunu yapmanın yolu, _params_ adlı özel bir Nextflow kavramıdır. Buna _params.input_ diyeceğiz.

Bunun yaptığı şey, bu input değişkenini CLI'da açığa çıkarmaktır ve Nextflow'u başlattığımızda çift tire kullandığımız yer burasıdır.

Buna istediğim gibi adlandırabilirim, _hello, greeting_ diyebilirim. Önemli değil. Orada ne yaparsam yapayım, bir boru hattı başlattığınızda CLI seçeneği olarak açığa çıkacaktır. Ve bu Nextflow'un gerçek bir sihirli numarasıdır çünkü bu parametrelerle iş akışı betiğinizi çok hızlı bir şekilde oluşturabileceğiniz anlamına gelir ve esasen boru hattınız için özel bir CLI oluşturuyorsunuz, başlattığınızda farklı seçenekleri anında özelleştirmeyi gerçekten kolaylaştırıyorsunuz.

Yani. Hadi deneyelim. Terminalimize geri dönelim. Burada _"nextflow run"_ komutumuz var. Ve şimdi _"--input"_ yapacağım, bu daha önce gördüğümüz _"params.input"_ ile eşleşiyor. Sanırım dokümanlarda Fransızca. Geraldine Fransızca konuşmayı seviyor. Ben İsveççe yapacağım çünkü İsveç'te yaşıyorum. Yani, "_Hej Världen_" diyeceğim ve enter'a basacağım.

Tek tırnak veya çift tırnak kullanabilirsiniz, bu sadece Bash'in bunu nasıl yorumladığını etkiler.

Nextflow boru hattını tamamen aynı şekilde çalıştırır. Çalışma dizininin ve her şeyin aynı olduğunu görebilirsiniz. Ama şimdi _"results/hello_world/output"_'a gidersem. Burada yerine güzel İsveççemizi görebiliriz.

Yani CLI'dan bir parametreye dinamik olarak bir girdi geçirdik. Bunu bir sürece girdi olarak geçirdik ve süreç bunu yorumladı ve bir script bloğuna koydu, bu da o betik sonucunun çıktısını dinamik olarak değiştirdi. Oldukça havalı.

Burada çok az söz dizimiyle oldukça karmaşık mantık. Ve umarım bunun şimdi nasıl ölçeklendiğini görebilirsiniz. Ve boru hatlarımızın mantığını ve özelleştirilebilirliğini Nextflow betiğine gerçekten bu şekilde inşa ediyoruz.

## 3.4. Komut satırı parametreleri için varsayılan değerler kullanın

Tamam, bu harika. Ama şimdi sorun şu ki, bu boru hattını her çalıştırdığımda, çalışması için tire, input yapmam gerekiyor.

Bu parametre olmadan çalıştırmaya çalışırsam, şimdi Nextflow bu parametreye ihtiyacı olduğunu ve ayarlanmadığını söyleyerek bir hata verecek. ve bu yüzden ne yapacağını bilmiyordu.

Bu arada havalı yeni bir şey. Geçmişte, Nextflow sadece boş bir string ile çalışırdı ve anlaşılması zor olacak her türlü garip hatanız olurdu. Ama yeni Nextflow söz dizimi ayrıştırıcısında, biraz daha dikkatli ve size hemen söylüyor.

Yani her seçeneği her zaman belirtmek istemeyiz. Mantıklı varsayılanlar belirtmek iyi bir uygulamadır. Peki bunu betiğimizde nasıl yaparız?

Bunu yazdığımızda, _params.input_'u doğrudan kullandığımız yere koyduğumuzu fark edeceksiniz. Yani açık çözüm bir varsayılan tanımlarız ve bunu burada workflow'daki betiğin en üstünde özel bir params bloğunda yaparız. Bu workflow betiğinde.

Yine, burada bazı yeni söz dizimi var, bu yüzden dikkat edin. Bu gerçekten havalı şeyler. Burada beklenecek parametrenin adını aldık.

Ve sonra bu iki nokta karakterinden sonra, değişkenin bir türünü tanımlıyoruz. Bunu yapmak zorunda değilsiniz, sadece boş bırakabilirsiniz, ama gerçekten güzel. Nextflow'a bir string beklediğimizi söyler ve buna göre davranır.

Örneğin, bunun yerine bir sayı istersek, float yazabiliriz ve bu bir kayan noktalı sayı istediğimizi söyler. Ve bununla çalıştırmaya çalışırsak, o zaman bir hata verecektir. Ona bir float olmayan bir string verirsek. Ve ayrıca bunu öyle geçirecek. String yaparsak, o zaman bunun bir string olduğunu bilir. Ve önünde sıfırlar olsa ve tamamen sayısal olsa bile, yine de gerçek bir string olarak geçirecektir.

Yani bu tür güvenliği Nextflow'un çok yeni bir özelliğidir, ama kodunuzu yazmayı ve çalıştırmayı daha güvenli hale getirmek için gerçekten güçlüdür.

Sonra bundan sonra bir eşittir sembolümüz var ve sonra burada varsayılan değer. Nextflow başlangıçta Barcelona'da yazıldı, bu yüzden burada biraz İspanyolca olması uygun görünüyor, varsayılan olarak _"Holà mundo!"_.

Tamam, o betiği kaydedeceğim, geri döneceğim, betiği _--input_ olmadan tekrar çalıştıracağım. Ve bu sefer çalışmalı ve _results_'ta yeni dosyamızı oluşturacak. Ve bu dosyada şimdi _"Holà mundo!"_ diyor.

Bu sadece bir varsayılan olsa da, yine de daha önce yaptığımız aynı şeyi yapamayacağımız anlamına gelmez. Geri dönüp burada eski betiğimi bulursam, _"Hej Världen"_, çünkü komut satırında _--input_ yapıyorum, bu varsayılanın üzerine yazacak ve bunu output.txt dosyasında tekrar kullanacak.

Yani betikteki bu sadece ayarladığım varsayılan değerdir.

İş akışımızı daha karmaşık hale getirip daha fazla parametre içerecek şekilde oluşturdukça, betiğin en üstündeki bu params bloğu hepsini tek bir yerde toplamaya başlayacak.

Ve betiğinizde bu oldukça güzel bir simetri elde edersiniz, burada etkili bir şekilde tüm iş akışı girdileriniz burada ve iş akışı çıktılarınız en altta. Ve iş akışınızın dış dünyaya arayüzünün ne olduğu çok açık. Böylece yeni söz dizimiyle yeni bir boru hattını çok hızlı bir şekilde alabilir ve nasıl kullanılacağını anlayabilirsiniz.

Son bir havalı şey. Bununla bir varsayılan değer ayarlamak zorunda değiliz. Params input yaparsak ama bir varsayılan değer ayarlamazsak, o zaman Nextflow'a bu parametrenin gerekli olduğunu söyler ve yine, boru hattı onsuz çalışmayı başaramaz, ama size null olması hakkında bir şey yerine daha kullanışlı bir hata mesajı verecektir.

Yani bir input beklediğimizi söylüyor, gerekli ama komut satırında belirtilmedi. Çok güzel.

Tamam, yani umarım şimdi Nextflow boru hattınızı değişken girdiler ve parametrelerle nasıl kuracağınız, varsayılanı nasıl ayarlayacağınız, türleri nasıl ayarlayacağınız, bir Boolean doğru yanlış bayrağı veya bir integer veya burada farklı türler olabilir konusunda açıktır. Bunları iş akışınıza nasıl geçireceğiniz, nereden geçtiği ve sonra sürecinize interpolate olduğu. Ve ayrıca Nextflow'u başlattığınızda komut satırında bunları nasıl özelleştireceğinizi biliyorsunuz. Bu basit bash komutumuzdan daha ilginç görünmeye başlıyor.

## 4. İş akışı yürütmelerini yönetin

Tamam. Sırada ne var? Bu bölümün son kısmı için, tüm farklı iş akışı yürütmelerini nasıl yöneteceğimiz hakkında biraz konuşacağız. Kenar çubuğuma ve Explorer'da work'ün altına bakarsanız, bir sürü farklı boru hattı çalıştırdığımı göreceksiniz ve bu work dizinleri oldukça uzun, çok sayıda var.

Ve diğer şey, daha önce söylediğim gibi, bu boru hattını her yeniden çalıştırdığımda, yeni bir work dizinleri seti oluşturuyor ve tüm süreçleri sıfırdan yeniden çalıştırıyor, bu iyi bir şey. Bu amaçlanan davranıştır. Tekrarlanabilir ve her şeyi taze bir şekilde yeniden üretiyor. Ama açıkçası, çok uzun süren süreçler çalıştırıyorsanız, yarı yolda çöktüyse veya boru hattının sonunda bir şeyi değiştirdiyseniz, boru hattınızı her zaman baştan başlatmak zorunda kalmak sinir bozucu.

## 4.1. Bir iş akışını -resume ile yeniden başlatın

Neyse ki, Nextflow daha önce neyin çalıştırıldığını ve neyin mevcut olduğunu bilmekte gerçekten iyidir ve bu eski sonuçları yeniden kullanmak çok basittir. Komutun sonuna sadece yeni bir bayrak ekliyoruz _"-resume"_.

Şimdi, input'ta iki tire var çünkü bu parametre. Resume'da sadece bir tire var çünkü bu bir çekirdek Nextflow seçeneği.

Bu insanları her zaman yanıltır, Nextflow'u uzun süredir kullanıyor olsanız bile. Bu yüzden her zaman bir veya iki tire hatırlayın. Çekirdek bir Nextflow seçeneği olup olmadığına bağlıdır.

Tamam, şimdi _-resume_ yapıyorum ve tamamen aynı iş akışını tekrar çalıştırıyorum. Ve bu sefer bir önemli farkla oldukça aynı görünmeli.

Buradaki çıktıda, sonuçların önbelleğe alındığını görebilirsiniz. Ve aslında, buradaki bu görev hash'i önceki çalıştırmayla tamamen aynı ve sadece o work dizinini bütünüyle yeniden kullandı. Girdiler ve çıktılar ve betik hepsi değiştirilmemişti. Ve bu yüzden sadece o dosyayı oradan aldı ve boru hattında aşağı akış adımları olsaydı, onları boru hattındaki bir sonraki adıma geçirirdi.

Yani hala tüm boru hattını baştan sona çalıştırıyor, ama yapabildiği her bir görev için önbelleğe alınmış sonuçları kullanıyor.

Şimdi, _-resume_ yaptığınızda, çalışma dizininizdeki son boru hattı çalıştırmasını devam ettirir, her ne olursa olsun. Ama aslında orada yaptığınız herhangi bir önceki çalıştırmadan devam edebilirsiniz. Ve şimdi oldukça çok yaptık.

## 4.2. Geçmiş yürütmelerin günlüğünü inceleyin

Hepsine bakmak için, _"nextflow run"_ yerine _"nextflow log"_ yapabiliriz ve bu bize tüm bu farklı.. Görebilmemiz için ekranımı biraz küçültmem gerekiyor, tüm bu farklı çalıştırmaları ne zaman yaptığımızı, oturum id'sini, komutu ve her şeyi gösteren güzel bir çıktı verecek.

Ve buraya bakabiliriz ve bunlardan herhangi birinin çalıştırma adını alabiliriz ve sonra o belirli olanlardan birini devam ettirebiliriz. Yani geri dönebilirim ve _hungry_ekeblad_ adlı olanı devam ettirebilirim. Ve bunu sadece _resume_'dan sonra koyuyorum.

Bu arada merak ediyorsanız, tüm bu sıfatlar ve bilim insanları isimleri Nextflow kaynak kodundadır. Gidip bularak ve favori bilim insanınızı ekleyerek Nextflow'a ilk pull request'inizi almanın gerçekten iyi bir yolu.

Ve neyse, bunu yaptım ve geri gitti ve bu iş akışı çalıştırmasından önbelleğe alınmış sonuçlara baktı, hala onları yeniden kullanabileceğini fark etti ve yaptı. Yani önbelleğe alınmış sonuçları tekrar aldım.

## 4.3. Eski work dizinlerini silin

Bu harika. Peki ya bu work dizinlerini temizlemek istersem? Burada tonlarca var. Tonlarca dosya var. Belki kesin olarak biliyorum ki son birkaç boru hattı çalıştırmasından devam etmek istiyorum, ama öncekilerle ilgilenmiyorum.

O zaman burada birini seçebilirim ve başka bir Nextflow komutu kullanabilirim, bu _"nextflow clean"_ ve _"nextflow clean"_ yapabilirim, _"-before"_ yapacağım ve belirli çalıştırma adı, bu durumda _reverent_pike_ ve _"-n"_ yapacağım, bu Nextflow'a sadece kuru bir çalıştırma yapmasını söyler. Yani bana sadece neyi sileceğini söyler. Aslında hiçbir şey yapmadan, bu work dizinlerini kaldırırdı.

Bu mantıklı görünüyor. Yani aynı komutu tekrar yapacağım, ama _"-n"_ yerine temizliği gerçekten yapmak için _"-f"_ yapacağım. Ve bu sefer gerçekten tüm bu dizinleri kaldırdı. Ve girip work dizinlerine bakarsam, şimdi çok daha hafif görünüyor. Harika.

Yani tüm yerel work dizinlerinizi önbelleği tamamen yok etmeden oldukça güvenli bir şekilde temizlemenin yolu budur. Böylece isterseniz hala devam edebilirsiniz.

Bu bayrakların ne için olduğunu unutursanız, her Nextflow komutu için _"nextflow help"_ yapabilirsiniz ve sonra komutun adını. Yani _"nextflow help clean"_ yaparsam, tüm farklı seçenekleri görebilirsiniz: _-after, -before, -but_, bu temizleme davranışını yapılandırmanın tüm farklı yolları. Oldukça havalı.

## Özet

Tamam, bu Hello Nextflow'un birinci bölümünün sonu. Kursa oldukça yoğun bir başlangıç, ama umarım şimdi bir Nextflow betiğinin nasıl göründüğü konusunda oldukça iyi bir anlayışa sahipsiniz; farklı anahtar parçalarla, süreçler, iş akışları, çıktılar ve parametreler. Bunları komut satırından temel geçersiz kılmalarla nasıl yapılandıracağınızı, dinamik bir script ile dinamik bir input bloğu nasıl yapacağınızı biliyorsunuz ve tüm iş yükü yürütmelerinizi nasıl yöneteceğinizi biliyorsunuz: zaten ne çalıştırdığınızı görmek, devam ettirmek, temizlemek. Çok fazla şey var. Uzun bir yol kat ettiniz. Yani bir mola vermek ve kısa bir yürüyüş ve bir fincan çay içmek istiyorsanız, şimdi muhtemelen iyi bir zaman. Bunu hak ettiniz.

Buradan itibaren, temelde bu temel üzerine inşa ediyoruz. Bunu nasıl daha karmaşık, daha güçlü hale getirebiliriz? Nasıl daha esnek hale getirebiliriz? Analizimizi ölçekte yapmak istediğimiz şeyleri yapın.

## Quiz

Şimdi web sayfasında birinci bölüm, hello world'e kaydırırsanız küçük bir quiz göreceksiniz ve bu, Nextflow eğitiminin bu versiyonu için yaptığımız yeni bir şey. Ve gidip bu bölümde yaptığımız tüm materyali anladığınızı kontrol etmek için kendinizi sınayabilirsiniz.

Bu bize gönderilmiyor veya başka bir şey, sadece tarayıcınızda saklanıyor. Yani cevaplarınızın ne olduğunu bilmiyoruz, ama bu sadece bir şeyi kaçırmadığınızdan veya yanlış anlamadığınızdan emin olmak için küçük bir kendi kendine kontrol. Ve istediğiniz kadar deneyebilirsiniz.

Eğer benim gibiyseniz, belki VS Code örneğinizde terminalde kalmak istiyorsunuz, bu durumda _quiz_ komutunu yazabilirsiniz ve sonra sadece hangi bölümde olduğunuzu söylersiniz. Yani _"Hello World"_ yapıyoruz ve sonra web tarayıcısında olan tamamen aynı quiz sorularını yapabilirsiniz, ama sadece terminalinizde.

Harika. Tamam. Umarım bundan keyif alırsınız. Biraz eğlenin ve bir dakika içinde Nextflow kanalları hakkında konuşmak için bir sonraki bölümde görüşürüz.
