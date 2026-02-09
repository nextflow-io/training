# Bölüm 3: Hello Workflow - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../03_hello_workflow.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş geldiniz ve özet

Merhaba ve Hello Nextflow'un üçüncü bölümüne tekrar hoş geldiniz. Bu bölümün adı Hello Workflow ve bu kurs bölümünde pipeline veya iş akışı ismini gerçekten haklı çıkarmaya başlıyoruz.

Şimdiye kadar tek bir süreçten oluşan basit pipeline betiğimizi alacağız ve ek süreçler eklemeye başlayacağız ve Nextflow'un bu orkestrasyon ve pipeline boyunca veri akışını nasıl yönettiğini göreceğiz.

Kod alanlarımıza geri dönelim. Temiz tutmaya çalışmak için tüm .nextflow\* dizinlerimi, work dizinlerini ve her şeyi sildiğimi göreceksiniz. Kursun önceki bölümlerinden kalan bu dosyalar hala etrafta dolaşıyorsa endişelenmeyin.

hello-workflow.nf adlı bir dosyadan çalışacağız. Daha önce olduğu gibi, bu temelde şu ana kadar oluşturduğumuz betiği temsil ediyor ve bize temiz bir başlangıç noktası veriyor. Ve yine çıktıda yolun artık hello_workflow olduğunu görebiliriz. Yani yayınlanan dosyalar results klasörünüzdeki farklı bir alt dizine gitmelidir.

Şu ana kadar nerede olduğumuzu özetlemek gerekirse, burada tek bir sürecimiz var, bir girdi selamlaması, bir çıktı selamlama dosyası. Ve sonra sadece bir dosyaya echo komutu yapan basit Bash betiği.

Tek bir iş akışı girdimiz var, buradaki params bloğu, bir yol beklediğini ve varsayılanın data/greetings.csv olduğunu söylüyor, bu da buradaki dosya.

Sonra iş akışının kendisinde, bir main bloğumuz var. Bir kanal oluşturuyoruz. CSV'yi satırlara ayrıştırıyoruz ve sonra her dizinin ilk elemanını alıyoruz ve bu kanalı o sürece geçiriyoruz, bu da üç görev oluşturuyor ve iş akışından, o sürecin çıktılarını yayınlıyoruz.

Ve son olarak, output bloğunda, Nextflow'a bu dosyaları bu kanaldan hello_workflow adlı dizine yayınlamasını söylüyoruz. Ve bu dosyaları soft link yapmak yerine kopyalamasını istiyoruz.

## 1. İş akışına ikinci bir adım ekleyin

Tamam, bu bölümde iş akışımıza ikinci bir süreç ekleyeceğiz. sayHello sürecinin çıktılarını alacağız ve bunları ikinci bir adımda işleyeceğiz, bu da bu dosyaların içindeki tüm harfleri convertToUppercase ile büyük harfe dönüştürecek.

Bu sadece aptalca bir örnek, yine sadece basit bir string işleme, ama size iş akışı içindeki mantığı nasıl alabileceğimizi gösteriyor.

Bunun için "tr" adında bir bash komutu kullanacağız, bu translate'in kısaltması. Sonsuza kadar var olan bir Unix komutu. Eğer aşina değilseniz, sizi suçlamam. Sanırım eğitimden önce hiç kullanmadım, ama terminalde çok hızlı bir şekilde deneyebilirsiniz. Eğer "echo 'hello world'" yaparsam ve sonra 'tr'ye pipe edersem ve sonra tırnak içinde karakter aralığı söylerseniz, yani A'dan Z'ye, küçük harf, ve sonra A'dan Z'ye büyük harf yapmak istiyorsunuz. Ve sadece şu harfleri şu harflere çevir diyor.

Ve enter'a bastığımda, her şeyi büyük harfe çevirdiğini görebilirsiniz. İnsanlara bağırmayı seviyorsanız çok güzel.

Yani bu, ikinci sürecimizde kullanacağımız çok basit bir bash komutu stili.

## 1.2. Büyük harfe çevirme adımını Nextflow süreci olarak yazın

Eğer betiğime geri dönersem, biraz hile yapacağım ve kodu eğitim dokümanlarından kopyalayacağım. Ama tam olarak neler olduğunu görebilirsiniz.

Burada yeni bir sürecimiz var. Buna convertToUpper dedik, ama istediğimiz her şeyi çağırabiliriz.

Daha önce yaptığımız gibi tek bir girdi yolumuz var. Bu bir value kanalı değil, bir path kanalı. Ve sonra tek bir çıktı.

Script bloğunda girdi dosyasında "cat" yapıyoruz. Ve bunu istediğimiz takdirde süslü parantezlere koyabiliriz. ve bu değişkeni alır. Ve aynı bash komutunu pipe'da çalıştırıyoruz ve sonuçları bu dosya adıyla bir dosyaya yazıyoruz ve bu çıktı yolu tarafından alınıyor.

Şimdi bu yeni süreçle bir şeyler yapmamız gerekiyor. Bu yüzden iş akışının farklı mantığını oluşturduğumuz yere gideceğiz ve o ilk süreçten sonra, ikinci sürecimizi çalıştıracağız. Yani convertToUpper burada sürecin adı.

Bir girdi alıyor, bu yüzden onu tek başına çağıramayız. İlk sürecin çıktısını işlemek istiyoruz. Yani tıpkı bununla yaptığımız gibi, bu sonuçları yayınladığımız sayHello out. Aynı sonuçları burada girdi olarak kullanmak istiyoruz, bu yüzden bunları kopyalayıp oraya koyabiliriz.

sayHello sürecini ".out" istiyoruz ve Nextflow bunun burada basit bir tek çıktı kaydı anlamına geldiğini biliyor, bu da bu dosya. Yani bu daha sonra ikinci bir sürece girdi olarak geçirilecek.

## 1.5. İş akışı çıktı yayınlamasını ayarlayın

Tamam. Ve son olarak, bu ikinci sürecin sonuçlarını gerçekten kaydetmek için, bunları iş akışından da yayınlamamız ve sonra output bloğunda tanımlamamız gerekiyor, daha önce olduğu gibi aynı sözdizimi. Yani bunu kopyalayabilir ve ikinci çıktılar diyebiliriz, ya da ne derseniz deyin.

İlgilendiğimiz süreç adını alın, convertToUpper out, ve sonra burada output bloğunda. Bunu ekleyin ve burada aynı özellikleri yapabiliriz. Yani bu dosyaların da Hello Workflow alt dizininde olmasını istiyoruz ve bunları da kopyalamak istiyoruz.

Harika. Hadi çalıştırmayı deneyelim. Eğer terminali açarsam ve "nextflow run hello-workflow.nf" yaparsam, ne yaptığını göreceğiz. Önceki bölümlerden farklı görünüp görünmediğine bakalım.

Yani Nextflow'u başlatıyor. Dokümanlarda bunu "-resume" ile yapmayı söylüyor, ama ben tüm work dizinlerimi sildim, bu yüzden burada bir fark yaratmazdı. Ama yapsaydınız, o zaman bu da işe yarardı.

Ve neredeyse tamamen aynı görünüyor. Ama şimdi burada ikinci bir çıktı satırı görebilirsiniz, az önce eklediğimiz ikinci sürecin adını görebilirsiniz. Ve gerçekten de, üç kez başarıyla çalıştığını görebilirsiniz.

Harika. Eğer önceki work dizinlerim etrafta olsaydı ve bunu "-resume" ile yapsaydım, bunlar önbelleğe alınmış olurdu, sadece pipeline'daki ilk adım. Çünkü bu çıktılar tamamen aynıydı, bu yüzden Nextflow bunları tekrar kullanmayı bilirdi.

Ve böylece gerekirse iş akışınızı adım adım yinelemeli olarak oluşturmak için -resume'u nasıl kullanabileceğinizi görebilirsiniz.

Tamam, işe yarayıp yaramadığını görmek için buradaki results dizinine bakalım. Burada daha fazla dosya olduğunu görebiliriz. Daha önce olduğu gibi ilk süreçten orijinal dosyalarımızı aldık. Ve gerçekten de, upper dosyalarımız var ve harfler tamamen büyük harf, yani işe yaramış. Görmek gerçekten çok güzel.

Bu work dizinlerinin içini kontrol etmek de ilginç. Daha önce olduğu gibi, buradaki hash work dizinlerine karşılık geliyor. Yani "ls work"e bakarsam, ve sonra onu genişletirsem, burada farklı dosyaları göreceğiz.

İlk süreçten çıktı dosyasını görüyoruz, bu buraya girdi olarak çekilmiş. Ve oluşturulan yeni çıktı dosyasını görebiliriz.

Şimdi bunu "-la" ile yaparsam, listelemek ve tüm dosyaları göstermek için, birkaç şey daha göreceğiz. İlk olarak, bu dosyanın aslında ilk sürece bir soft link olduğunu göreceksiniz. Dosya alanından tasarruf etmek için bu temelde her zaman mümkünse bir soft link'tir. Dosyaları burada yayınlamıyoruz ve sadece ilk görevden o dosyaya referans veriyor, ikinci bir göreve, böylece her şey bir çalışma dizini içinde kapsüllenmiş ve güvenli ve diğer her şeyden izole edilmiş.

Ve orada olması gerekiyor çünkü .command.sh dosyasına bakarsak, yani "cat work/b8/56\*" yaparsam, buradaki dosya parçalarının göreceli olduğunu görebilirsiniz, yani aynı çalışma dizinine soft link yapılmış o girdi dosyasını cat yapıyor.

Yani her work dizini böyle görünecek. Nextflow'da baktığınızda, orada aynı çalışma dizinine sahnelenmiş tüm girdi dosyalarına sahip olacaksınız. Ve sonra oluşturulan herhangi bir çıktı dosyanız da olacak. Yani bu harika. Beklediğimiz gibi görünüyor.

## 2.1. Toplama komutunu tanımlayın ve terminalde test edin

Tamam, iş akışımıza geri dönelim. Yapmak istediğimiz bir sonraki adım ne?

Şimdi iki sürecimiz var ve bu tek CSV dosyasını alıyorlar, ayrıştırıyorlar ve bölüyorlar. Ve sonra bu süreçlerin her biri için üç görevimiz var ve Nextflow tüm bunların paralelleştirilmesini yönetiyor, bu yüzden mümkün olduğunda yan yana çalışıyor.

İşi paralel çalıştırmak için böyle bölme şekli çok yaygın. Ve bunun tersi de her şeyi geri toplamak. Yani iş akışımızdaki son sürecimizle yapacağımız şey bu, burada üçüncü bir tane olacak, bu üç farklı çıktıyı alıp hepsini tek bir dosyada birleştirecek.

Bunu terminalde oldukça basit bir şekilde yapabiliriz, bunun nasıl görüneceğine dair bir fikir edinmek için.

Results klasörüne gidersem. Yani, "cd results/hello_workflow/", ve burada tüm UPPER dosyalarımız var. Sadece "cat" kullanabilirim, bunu o dosyanın içeriğini yazdırmak için kullanıyoruz ve "cat"e birden fazla dosya verebilirsiniz ve bunları birbiri ardına okuyacaktır.

Yani "UPPER-\*" diyebilirim, bu bana Bash genişletmesiyle aynı üç dosya adı listesini verir. Ve combined.txt diyebilirim. Sanırım dokümanlarda tam dosya adlarını listeler, ama aynı şeyi yapıyor.

Şimdi, "cat combined.txt" kullanırsam, bu üç dosyanın tümünün dosya içeriğine sahip olduğumuzu görebiliriz.

Yani temelde bu sürecin yapacağı şey bu, önceki bir süreçten tüm farklı çıktı dosyalarını tek bir süreç görevinde vermeye çalışacağız ve sonra bunları "cat" ile birleştirecek ve çıktı dosyasını kaydedeceğiz.

## 2.2. Toplama adımını yapmak için yeni bir süreç oluşturun

Tamam, yeni sürecimizi ekleyelim. Bunu eğitim materyallerinden yapıştıracağım ve burada okuyucu için biraz alıştırma bıraktığını bu soru işaretleriyle görebilirsiniz. Ama sürecin genel hatlarının temelde terminalde yaptığımız şey olduğunu görebilirsiniz, bir grup girdi dosyasının "cat"ini yapıyoruz ve bunu burada collected adlı bir çıktı dosyasına yazıyoruz ve sonra çıktı yine o tek yolu bekliyor.

Yani burada bir tür girdiye ihtiyacımız var ve bunlar bir dizi yol olacak. Yani yine, bir girdi yol kanalı tanımlıyoruz ve buna input_files diyelim. Şimdi, bu daha önce bize burada tek bir yol verdi, ama bir yol burada birden fazla dosyaya da sahip olabilir, hala tek bir bildirim olsa bile.

Bunu buraya kopyalayacağım çünkü bu dosyaları "cat" yapmak istiyoruz. Ve burada bir diziyi yazdırmak gibi bazı sorunlarımız olabileceğini düşünebilirsiniz, ama Nextflow genellikle bu konuda oldukça mantıklı. Ve eğer bunun gibi içinde birden fazla dosya olan bir kanal verilirse, hepsini boşluk ayırıcılarla bir araya getirecektir. Yani bu bize doğru sözdizimini verecek.

Bu harika. Şimdi yeni sürecimizi bağlayalım. Workflow'a gidiyorum. Çıktıları birleştir diyeceğim, yeni süreç adı ve daha önce olduğu gibi. Bu önceki süreci alacağım, convertToUpper ve ".out" yapacağım.

Harika. Hadi deneyelim ve terminalde işe yarayıp yaramadığına bakalım. Sadece birkaç dizin yukarı gidersem ve sonra Nextflow komutunu yeniden çalıştırırsam, ne olacağını göreceğiz.

Yani iş akışı başlatıldı ve şimdi üç farklı süreç adımız olduğunu görebilirsiniz, bu harika. İlk ikisi daha önce olduğu gibi görünüyor ve üçüncü yeni olan çalışıyor, bu iyi.

Ancak, burada biraz garip bir şey var. Bu çıktı dosyalarını tek bir dosyada birleştirmek istedik, ama bu sürecin bir kez değil üç kez çalıştığını görebiliriz.

Gerçekten de, bu work dizinlerinden birine girersek. Ve "cat work/" "collected" yaparsak, göreceğiz. Burada sadece tek bir kelime var, üç değil.

Ve olan şey, Nextflow'un önceki adımlarda yaptığı gibi bu paralelleştirmeyi sürdürmesi oldu. Ve bu süreç bize üç elemanlı bir kanal verdi ve bu üç kanal elemanı bizim downstream sürecimize geçirildi, bu da üç süreç görevi oluşturdu.

Temelde üç ayrı kez toplamaya çalıştı ve her seferinde sadece tek bir dosyası vardı, bu yüzden sadece tek dosyayı bir çıktıya cat yaptı ve aslında bunu .command.sh dosyasında da görebiliriz.

Eğer .command.sh yaparsam, burada sadece tek bir dosya adı olduğunu görebiliriz ve bu çalışma dizinine sadece tek bir dosya sahnelendi.

## 2.3. Toplama adımını iş akışına ekleyin

Yani bir şekilde Nextflow'a önceki bir süreçten tüm bu çıktıları bir araya getirmesini ve bunları bu downstream sürece üç yerine tek bir kanal elemanı olarak vermesini söylememiz gerekiyor.

Bunu _collect_ adlı bir kanal operatörü ile yapıyoruz.

Bu, Nextflow pipeline'larında her zaman göreceğiniz süper kullanışlı bir operatör. Bu burada bir kanal, bu çıktı kanalı, tıpkı en üstte oluşturduğumuz gibi. Ve böylece daha önce yaptığımız gibi ona kanal operatörleri ekleyebiliriz. Sadece nokta yapabiliriz ve sonra bu durumda, collect, parantezler.

Ve ihtiyacımız olan tek şey bu. Bu daha sonra bu kanalı bu sürece geçirilmeden önce manipüle edecek.

Ona ne olduğunu görmek istiyorsanız, bunu burada da görüntüleyebiliriz. Yani burada, bu bu süreci çalıştırmakla hiç ilgili değil, bu yüzden o süreci çalıştırdıktan sonra herhangi bir noktaya koyabilirim. Ama aynı çıktı kanalını alıyoruz ve .view ile bakıyoruz ve sonra .collect.view ile tekrar bakıyoruz.

Ve bunu çalıştırdığımızda, bize o kanalın iki farklı yapısını gösterecek, collect'ten önce ve sonra. Hadi şimdi deneyelim. Tamam, bazı çıktılar oldukça uzun olduğu için biraz uzaklaştırdım, ama pipeline'ı çalıştırırsam, işe yarayıp yaramadığını göreceğiz.

Üçüncü sürecin sadece bir kez çalışmasını umuyorum, çünkü çıktıları topluyor ve gerçekten de, collectGreetings'i bir'den bir olarak görebilirsiniz. Yani bu sadece bir görev çalıştırdı.

Ve sonra view ifadelerine bakarsak, önce'nin üç elemanı için üç view ifademiz var, her birinde bir dosya yolu var.

Ve sonra o collect ifadesinden sonra, bu sadece bir kez tetiklendi çünkü o kanalda tek bir eleman var. Ve şimdi üç farklı dosya yolunun bu listesine sahibiz.

Tam olarak umduğumuz şey bu. Ve umarım görebilirsiniz, bu temelde CSV dizilerinden ayrı kanal elemanlarına gitmek için yaptığımız o "map" operatörünün tersi. Şimdi ayrı kanal elemanlarını alıyoruz ve tek bir diziye geri koyuyoruz.

Harika, bu view ifadelerini temizleyebiliriz. Bunlara artık ihtiyacımız yok. Bir sonraki adıma geçebiliriz.

Daha ileri gitmeden önce ve unutmadan önce, burada yeni bir publish ifadesi ekleyeceğim. Üçüncü çıktı. Bunu iş akışınızda daha anlamsal ve açıklayıcı bir şey olarak adlandırabilirsiniz. Ve sonra bunu output bloğuna tekrar ekleyeceğim ve path 'hello_workflow' mode 'copy' diyeceğim. Böylece bu süreç tarafından oluşturulan çıktı dosyası buradaki results klasörümüze kaydedilecek.

Bunun işe yarayıp yaramadığını hızlıca kontrol etmek için. Şimdi biraz daha temiz olmalı çünkü o view ifadelerimiz yok. Ve burada yeni çıktı dosyamızı alıp almadığımızı göreceğiz. Bir'den bir görev çalıştı, collected adında yeni bir dosya aldık ve şimdi bu üç kelimenin hepsine sahibiz. Harika. Sırada ne var?

## 3. Bir sürece ek parametreler geçirin

Tamam. Sonra tek bir sürece birden fazla girdinin nasıl işleneceğine bakacağız. Şimdiye kadar tüm süreçlerimizin sadece bir şeyi girdi olarak aldığını görebilirsiniz. Hepsinin girdileri altında tek bir satırı var.

Bunu Nextflow'un farklı bir batch tanımlayıcısı belirtmesine izin vererek göstereceğiz, böylece belki bu iş akışını birden çok kez çalıştırırsınız ve her seferinde farklı bir batch ID verebilirsiniz.

Basitçe collectGreetings için girdide ikinci bir satır ekleyeceğim. Ve buna "val" diyeceğim, çünkü bu bir string. Şimdi bu bir değer, bir yol değil ve buna "batch_name" diyeceğim.

Sonra bu değişkeni kullanmak için aşağıdaki betiği düzenleyeceğim ve bunu eğitim materyaliyle aynı yere koymaya çalışacağım. Yani bunu bu dosya yolunun ortasına koyuyorum COLLECTED-$\{batch_name\}-output.

Henüz bitmedi. Nextflow'a çıktı dosya adlarının ne olacağını söylememiz gerektiğini unutmayın. Yani burada da aynı şeyi yapmamız gerekiyor: COLLECTED-$\{batch_name\}-output.txt".

Harika. Nextflow şimdi ikinci bir değişken girdisi alıyor ve bunu betiğe ve çıktıya interpolasyon yapıyor.

Son bir şey, şimdi bunun nerede çağrıldığını bulmamız ve sürece ikinci girdiye geçmemiz gerekiyor. Bu, başka herhangi bir dildeki bir fonksiyona başka herhangi bir girdi gibi.

Eğitimde daha önce yaptığımız gibi, burada özel "params"ı kullanacağım ve buna "params.batch" diyeceğiz, böylece bir -- batch CLI seçeneğimiz olabilir. Ve şimdi sürecimizin burada virgülle ayrılmış iki ayrı girdisi olduğunu görebilirsiniz, bunlar geçiriliyor.

Sırayı doğru yapmak gerçekten önemli, yani burada kanal için argümanların sırası ve sonra param eşleşmelidir. Kanal ve oradaki batch adı. Bu sadece konumsal eşleştirme.

Tamam. Bu pipeline'ı şimdi doğrudan --batch ile çalıştırabilirim, ama önce doğru olanı yapalım ve bunu burada Params'ta tanımlayalım. Yani bunu batch'e ekleyeceğim ve sonra bunun bir string olduğunu söyleyeceğiz ve bir varsayılan verelim. Yani sadece batch diyelim. Tamam? Şimdi iş akışını çalıştırmayı deneyelim.

--batch Trio. Sanırım eğitim materyalinde öyle diyor, ama orada istediğimiz herhangi bir string kullanabiliriz. Ve umarım bu sonuç çıktı dosyasının buraya geldiğini göreceğiz.

Ve gerçekten de, COLLECTED-trio-output - bu düzgün çalıştı. Dosyamızı yeniden adlandırdı. Ve şimdi bunun yararlı olduğunu hayal edebilirsiniz çünkü bunu farklı bir batch adıyla tekrar çalıştırırsam, replicate_two gibi, o zaman bize burada farklı bir batch adı verecek.

Ve bu durumda çıktı dosyalarının üzerine yazmayacak. Yani bu güzel.

## 4. Toplayıcı adıma bir çıktı ekleyin

Tamam, şimdi sürecimize birden fazla girdimiz var. Ama birden fazla çıktı oluşturmak istersek ne olur? Buradaki örneğimiz, bu süreç için bir rapor oluşturacağız, sadece kaç dosyanın toplandığını söyleyerek.

Ve bunu burada bir echo komutuyla yapacağız. Yani echo diyebiliriz. Vardı, bunu eğitim materyalinden kopyalayacağım, böylece beni yazarken izlemek zorunda kalmazsınız.

Bu batch'te $\{count_greetings\} selamlama vardı ve bunu şimdi $\{batch_name\} adlı yeni bir dosyaya kaydedin, yani aynı değişken, bunu istediğimiz kadar yeniden kullanabiliriz, report.txt.

## 4.1.1. Toplanan selamlama sayısını sayın

Bunu bir şekilde hesaplamamız gerekiyor. Bu mantığı istersek Bash betiğinde yapabiliriz, Bash mantığını kullanarak. Ancak, script bloğu içinde ve alıntılanan bölümün üstünde olduğu sürece, doğrudan Nextflow kodu içinde de betik yazabiliriz.

Buradaki hiçbir şey son render edilen betiğe dahil edilmeyecek ve sadece bir görevi render ettiğinde Nextflow tarafından yürütülecek.

Yani burada sadece biraz mantık yapıyoruz. count_greetings adında yeni bir değişken oluşturuyoruz. Burada input files kanalını alıyoruz ve üzerinde .size() çağırıyoruz.

Tamam, o fonksiyon bana bu değişkene bir sayı verecek ve şimdi uyarımız gitti çünkü bu değişken tanımlanıyor.

Tamam, yani work dizininde o ikinci dosyayı oluşturuyoruz, ama Nextflow'a bunu bu sürecin yayınlanmış bir çıktısı olarak beklemesini söylememiz gerekiyor. Bunu ilk dosya için yaptığımız gibi tamamen aynı sözdizimi ile yapıyoruz.

Path diyoruz çünkü bu, yine, istersek burada bir değişken yayınlayabiliriz "val" ile, ama "path" diyeceğiz. Ve sonra beklenen dosya adı. Burada vurgulanmadığına dikkat edin. Bunun nedeni tek tırnak kullanmam. Çift tırnak kullanmam gerekiyor.

## 4.1.2. Rapor dosyasını yayınlayın ve çıktıları adlandırın

Tamam, bu harika. Ve şimdi bu çıktılara burada tıpkı burada yaptığım gibi erişebiliriz. Ama şimdi bu farklı nesnelerin bir dizisi, bu yüzden ilkini almak için collectGreetings.out[0] yapabilirim, ya da yeni raporumuz olan ikincisini almak için bir.

Ama bunu yapmayı pek sevmiyorum çünkü dizin sayımını karıştırmak oldukça kolay. Ve orada oturur çok sayar ve yeni bir çıktı eklersiniz ve aniden her şey bozulur. Yani

her şeye isimle referans vermek çok daha güzel. Ve bunu burada "emit" adlı özel bir anahtarla yapabiliriz.

Yani buna istediğimizi çağırabiliriz. Buna emit outfile diyelim ve emit reports. Bunları tanımlarsanız ve bir veya birçoğunda yapabilirsiniz, size kalmış. Şimdi buraya gidebilirim ve bunun yerine dot out dot reports gidebilirim ve sadece isimle çağırabilirim, bu kodunuzu okuduğunuzda anlamak çok daha kolay ve koddaki değişikliklere karşı daha güvenli.

Burada .out.report ekledim, ama aslında iki farklı çıktının yayınlanması gerekiyor. Bu yüzden bunu collected ve report gibi daha ilginç bir şey olarak yeniden adlandıracağım ve ona ne dedim? Out file dedim, üzgünüm. Yani burada emit adı outfile ve report. çünkü iki farklı çıktı kanalı yayınlıyoruz ve bu yüzden publish bloğunda her ikisine de referans vermemiz gerekiyor.

Sonra bunları output bloğunda da tanımlamamız gerekiyor. Yani bunu collected olarak yeniden adlandırdım ve yine, reports için, burada biraz ayrıntılı, ama yeni bir iş akışını okumaya geldiğinizde, burada yan yana listelenen tüm farklı çıktıları, tüm farklı kanalları görmek gerçekten yararlı ve bunu daha az ayrıntılı yapmanın yolları var, bunlara daha sonra değineceğiz.

Tamam, hadi deneyelim ve iş akışımızı çalıştıralım ve ne olacağını görelim.

Umarım şimdi temelde daha önce olduğu gibi çalışmalı. Ve burada replicate_two adında yeni bir çıktı dosyası alacağız, report. Ve işte burada. Açıldı ve bu batch'te üç selamlama var diyor, ki bu beklediğimiz şey, yani mükemmel.

Eğer work dizinine girersem, sadece Nextflow kodunda yürütüldüğünü size kanıtlamak için, bash betiği yerine, cat work/ command.sh'ye gidebilirim ve burada sadece bu stringi doğrudan echo ettiğini göreceksiniz. Bu batch'te üç selamlama vardı ve bu yüzden o değişken Nextflow tarafından interpolasyon yapıldı. .command.sh dosyasını yazmadan önce script bloğunda hesaplandı. Yani sonuçta ortaya çıkan değişken hesaplaması, bu durumda hesaplama ortamınızda yürütülmeden önce temelde buna sabit kodlanmış.

Ve böylece script arasındaki ayrımı görebilirsiniz. Buradaki blok ve üstündeki herhangi bir şey. Umarım bu mantıklıdır.

## Özet ve quiz

Tamam, Hello Nextflow'un bu bölümünün sonu. Yani daha önce olduğu gibi, gidip quiz'e bakın. Web sayfasında veya CLI'da yapın, bazı sorulardan geçin ve ele aldığımız materyallerden bazılarını anlayıp anlamadığınızı kontrol edin. Anlamadığınız bir şeyi vurgulayan bir şey olup olmadığına bakın. Çok fazla soru yok. Yapması güzel ve kolay. Ya da bunu burada web sayfasında da yapabilirsiniz.

Ve küçük bir mola verin, biraz dolaşın ve geri gelin ve Hello Nextflow'un dördüncü bölümünde bize katılın, modüllerden bahsedeceğiz. Çok teşekkür ederim.
