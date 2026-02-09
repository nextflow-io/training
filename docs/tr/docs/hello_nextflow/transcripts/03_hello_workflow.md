# Bölüm 3: Hello Workflow - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../03_hello_workflow.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş geldiniz ve özet

Merhaba ve Hello Nextflow'un üçüncü bölümüne tekrar hoş geldiniz. Bu bölümün adı Hello Workflow ve bu kurs bölümünde pipeline veya iş akışı adını gerçekten hak etmeye başlıyoruz.

Şu ana kadarki tek süreçli basit pipeline betiğimizi alacağız ve ek süreçler eklemeye başlayıp Nextflow'un bu orkestrasyon ve pipeline boyunca veri akışını nasıl yönettiğini göreceğiz.

Code spaces'e geri dönelim. Temiz tutmak için tüm .nextflow\* dizinlerini, work dizinlerini ve her şeyi sildim. Kursun önceki bölümlerinden bu dosyalar hala etrafta dolanıyorsa endişelenmeyin.

hello-workflow.nf adlı bir dosyadan çalışacağız. Daha önceki gibi, bu temelde bu noktaya kadar oluşturduğumuz betiği temsil ediyor ve bize temiz bir başlangıç noktası veriyor. Ayrıca, aşağıdaki çıktıda yolun şimdi hello_workflow olduğunu görebiliriz. Yani yayınlanan dosyalar results klasörünüzdeki farklı bir alt dizine gitmeli.

Şu ana kadar nerede olduğumuzu özetlemek gerekirse, burada tek bir girdi selamlama, tek bir çıktı selamlama dosyası olan tek bir sürecimiz var. Ardından bir dosyaya sadece bir echo komutu yapan basit Bash betiği.

Tek bir iş akışı girdimiz var, burada params bloğu, bir path beklediğini ve varsayılanın data/greetings.csv olduğunu söylüyoruz, bu da buradaki dosya.

Sonra iş akışının kendisinde bir main bloğumuz var. Bir kanal oluşturuyoruz. CSV'yi satırlara ayrıştırıyoruz ve ardından her dizinin ilk elemanını alıyoruz ve bu kanalı o sürece aktarıyoruz, bu da üç görev oluşturuyor ve bu sürecin çıktılarını iş akışından yayınlıyoruz.

Ve son olarak, output bloğunda, Nextflow'a bu dosyaları bu kanaldan hello_workflow adlı dizine yayınlamasını söylüyoruz. Ve bu dosyaları soft link yerine kopyalamasını istiyoruz.

## 1. İş akışına ikinci bir adım ekleyin

Tamam, bu bölümde iş akışımıza ikinci bir süreç ekleyeceğiz. sayHello sürecinin çıktılarını alacağız ve bunları bu dosyaların içindeki tüm harfleri büyük harfe dönüştürecek convertToUppercase olan ikinci bir adımda işleyeceğiz.

Bu sadece aptalca bir örnek, yine sadece basit bir string işleme, ama size iş akışı içindeki mantığı nasıl alabileceğimizi gösteriyor.

Bunun için "tr" (translate'in kısaltması) adında bir bash komutu kullanacağız. Sonsuza kadar var olan bir Unix komutu. Eğer buna aşina değilseniz, sizi suçlamam. Sanırım bunu eğitimden önce hiç kullanmadım, ama terminalde çok hızlı bir şekilde deneyebilirsiniz. Eğer "echo 'hello world'" yaparsam ve ardından 'tr'ye pipe edersem ve ardından tırnak içinde karakter aralığı söylerseniz, yani küçük harf A'dan Z'ye, ve ardından büyük harf A'dan Z'ye yapmak istersiniz. Ve sadece bu harfleri bu harflere çevir der.

Ve enter'a bastığımda, her şeyi büyük harfe çevirdiğini görebilirsiniz. İnsanlara bağırmayı seviyorsanız çok güzel.

Bu ikinci sürecimizde kullanacağımız çok basit bir bash komutu stili.

## 1.2. Büyük harfe çevirme adımını Nextflow süreci olarak yazın

O yüzden betiğime geri dönersem, biraz hile yapacağım ve kodu eğitim dokümanlarından kopyalayacağım. Ama tam olarak neler olduğunu görebilirsiniz.

Burada yeni bir sürecimiz var. Buna convertToUpper dedik, ama istediğimiz her şeyi çağırabiliriz.

Daha önce yaptığımız gibi tek bir girdi path'imiz var. Bir value kanalı değil, bir path kanalı. Ve sonra tek bir çıktı.

Script bloğunda girdi dosyasında "cat" yapıyoruz. Ve istarsak bunu kıvırcık parantezlere koyabiliriz. ve bu o değişkeni alır. Ve aynı bash komutunu pipe'da çalıştırıyoruz ve sonuçları bu dosya adıyla bir dosyaya yazıyoruz ve bu çıktı path tarafından alınır.

Şimdi bu yeni süreçle bir şeyler yapmamız gerekiyor. O yüzden iş akışının farklı mantığını oluşturduğumuz yere gideceğim ve bu ilk süreçten sonra ikinci sürecimizi çalıştıracağız. Yani convertToUpper burada sürecin adı.

Bir girdi alıyor, bu yüzden onu tek başına çağıramayız. İlk sürecin çıktısını işlemek istiyoruz. Yani tıpkı bununla yaptığımız gibi, bu sonuçları yayınladığımız sayHello out. Bu aynı sonuçları burada girdi olarak kullanmak istiyoruz, bu yüzden bunları kopyalayıp oraya koyabiliriz.

sayHello sürecini ".out" istiyoruz ve Nextflow bunun burada bu dosya olan basit tek bir çıktı kaydı anlamına geldiğini biliyor. Yani bu daha sonra ikinci bir sürece girdi olarak aktarılacak.

## 1.5. İş akışı çıktı yayınlamasını ayarlayın

Tamam. Ve son olarak, bu ikinci sürecin sonuçlarını gerçekten kaydetmek için, onları da iş akışından yayınlamamız ve ardından output bloğunda tanımlamamız gerekiyor, daha öncekiyle aynı sözdizimi. Bu yüzden bunu kopyalayabilir ve second outputs diyebiliriz, veya ne derseniz deyin.

İlgilendiğimiz süreç adını alın, convertToUpper out, ve sonra burada output bloğunda. Bunu ekleyin ve burada aynı nitelikleri yapabiliriz. Yani bu dosyaların da Hello Workflow alt dizininde olmasını istiyoruz ve onları da kopyalamak istiyoruz.

Harika. Çalıştırmayı deneyelim. Terminali açarsam ve "nextflow run hello-workflow.nf" yaparsam, ne yapacağını göreceğiz. Önceki bölümlerden farklı görünüp görünmediğini görelim.

Yani Nextflow'u başlatıyor. Dokümanlarda bunu "-resume" ile yapmayı söylüyor, ama ben tüm work dizinlerimi sildim, bu yüzden burada bir fark yaratmazdı. Ama yapsaydınız, o zaman bu da çalışırdı.

Ve neredeyse tamamen aynı görünüyor. Ama şimdi burada ikinci bir çıktı satırı görebilirsiniz, burada az önce eklediğimiz ikinci sürecin adını görebilirsiniz. Ve gerçekten de, üç kez başarıyla çalıştığını görebilirsiniz.

Harika. Önceki work dizinlerim etrafta olsaydı ve bunu "-resume" ile yapmış olsaydım, bunlar sadece pipeline'daki ilk adımı önbelleğe alınmış olurdu. Çünkü bu çıktılar tamamen aynıydı, bu yüzden Nextflow bunları tekrar kullanmayı bilirdi.

Ve böylece gerekirse iş akışınızı adım adım yinelemeli olarak oluşturmak için -resume'u nasıl kullanabileceğinizi görebilirsiniz.

Tamam, işe yarayıp yaramadığını görmek için buradaki results dizinine bakalım. Burada birkaç dosya daha var. İlk süreçten daha önce yaptığımız gibi orijinal dosyalarımız var. Ve gerçekten de, upper dosyalarımız var ve harfler tamamen büyük harf, yani işe yaramış. Görmek gerçekten güzel.

Bu work dizinlerinin içine bakmak da ilginç. Daha önce olduğu gibi, buradaki hash work dizinlerine karşılık geliyor. Bu yüzden "ls work"e bakarsam ve sonra onu genişletirsem, burada farklı dosyaları göreceğiz.

İlk süreçten, buraya girdi olarak çekilen çıktı dosyasını görüyoruz. Ve oluşturulan yeni çıktı dosyasını görebiliriz.

Şimdi bunu "-la" ile tüm dosyaları listele ve göster yaparsam, birkaç şey daha göreceğiz. İlk olarak, bu dosyanın aslında ilk sürece bir soft link olduğunu göreceksiniz. Bu, dosya alanından tasarruf etmek için mümkünse her zaman temelde bir soft link'tir. Dosyaları burada yayınlamıyoruz ve bu sadece o dosyaya ilk görevden ikinci göreve referans veriyor, böylece her şey bir çalışma dizininde kapsüllenmiş ve güvenli ve diğer her şeyden izole edilmiş.

Ve bu orada olması gerekiyor çünkü .command.sh dosyasına bakarsak, yani "cat work/b8/56\*" yaparsam, buradaki dosya parçalarının göreceli olduğunu görebilirsiniz, bu yüzden aynı çalışma dizinine soft link edilmiş olan o girdi dosyasını catliyor.

Yani her work dizini böyle görünecek. Nextflow'da baktığınızda, o work dizinine staged edilmiş tüm girdi dosyalarına sahip olacaksınız. Ve sonra oluşturulan herhangi bir çıktı dosyanız da olacak. Yani bu harika. Beklediğimiz gibi görünüyor.

## 2.1. Toplama komutunu tanımlayın ve terminalde test edin

Tamam, iş akışımıza geri dönelim. Yapmak istediğimiz bir sonraki adım ne?

Şimdi iki sürecimiz var ve bu bir CSV dosyasını alıyor, ayrıştırıyor ve bölüyor. Ve sonra bu süreçlerin her biri için üç görevimiz var ve Nextflow tüm bunların paralelleştirilmesini yönetiyor, böylece mümkün olduğunda yan yana çalışıyor.

İşleri paralel çalıştırmak için böyle bölme yolu çok yaygın. Ve bunun tersi her şeyi geri toplamaktır. İş akışındaki son sürecimizle yapacağımız şey budur, burada üçüncü bir tane olacak, bu üç farklı çıktıyı alıyor ve hepsini tek bir dosyada birleştiriyor.

Bunu terminalde oldukça basit bir şekilde yapabiliriz, bunun nasıl görüneceğine dair bir fikir edinmek için.

Results klasörüne gidersem. Yani, "cd results/hello_workflow/", ve burada tüm UPPER dosyalarımız var. Sadece "cat" kullanabilirim, bunu o dosyanın içeriğini yazdırmak için kullanıyoruz ve "cat"e birden fazla dosya verebilirsiniz ve birini diğerinin ardından okuyacak.

Yani "UPPER-\*" diyebilirim, bu bana Bash genişletmesiyle aynı üç dosya adı listesini veriyor. Ve combined.txt diyebilirim. Sanırım dokümanlarda tam dosya adlarını listeler, ama aynı şeyi yapıyor.

Şimdi, "cat combined.txt" kullanırsam, o üç dosyanın tümünün dosya içeriklerine sahip olduğumuzu görebiliriz.

Yani temelde bu sürecin yapacağı şey budur, ona önceki bir süreçten gelen tüm farklı çıktı dosyalarını tek bir süreç görevinde vermeye çalışacağız ve sonra onları "cat" ile birleştireceğiz ve çıktı dosyasını kaydedeceğiz.

## 2.2. Toplama adımını yapmak için yeni bir süreç oluşturun

Tamam, yeni sürecimizi ekleyelim. Bunu eğitim materyallerinden yapıştıracağım ve burada bize bu soru işaretleriyle okuyucu için biraz alıştırma bıraktığını görebilirsiniz. Ama sürecin genel hatlarının temelde terminalde yaptığımız şey olduğunu görebilirsiniz, bir grup girdi dosyasının "cat"ini yapıyoruz ve burada collected adlı bir çıktı dosyasına yazıyoruz ve sonra çıktı yine o tek path'i bekliyor.

Bu yüzden burada bir çeşit girdiye ihtiyacımız var ve bunlar bir dizi path olacak. Yine, bir girdi path kanalı tanımlıyoruz ve buna input_files diyelim. Şimdi, bu daha önce bize burada tek bir path verdi, ama bir path aynı zamanda burada birden fazla dosyaya da sahip olabilir, tek bir bildirim olmasına rağmen.

Bunu buraya kopyalayacağım çünkü bu dosyaları catlamak istiyoruz. Ve burada dizi yazdırmak gibi bazı sorunlarımız olabileceğini düşünebilirsiniz, ancak Nextflow genellikle bu konuda oldukça mantıklı. Ve eğer bunun gibi içinde birden fazla dosya olan bir kanal verilirse, hepsini boşluk ayırıcılarla bir araya getirecektir. Yani bu bize doğru sözdizimini verecek.

Bu harika. Şimdi yeni sürecimizi bağlayalım. Workflow'a gidiyorum. Çıktıları birleştir diyeceğim, yeni süreç adı, ve daha önceki gibi. Bu önceki süreci alacağım, convertToUpper ve ".out" yapacağım.

Harika. Deneyelim ve terminalde çalışıp çalışmadığını görelim. Sadece birkaç dizin yukarı çıkıp sonra Nextflow komutunu yeniden çalıştırırsam, ne olacağını göreceğiz.

Yani iş akışı başlatıldı ve şimdi üç farklı süreç adımız olduğunu görebilirsiniz, bu harika. İlk ikisi daha önce olduğu gibi görünüyor ve üçüncü yeni olan çalışıyor, bu iyi.

Ancak, burada biraz garip bir şey var. Bu çıktı dosyalarını tek bir dosyada birleştirmek istedik, ama bu sürecin bir kez değil üç kez çalıştığını görebiliriz.

Gerçekten de, bu work dizinlerinden birine girersek. Ve "cat work/" "collected" yaparsak, göreceğiz. Burada sadece tek bir kelime var, üç değil.

Ve olan şey, Nextflow'un önceki adımlarda olduğu gibi bu paralelleştirmeye devam etmesi oldu. Ve bu süreç bize üç elemanlı bir kanal verdi ve bu üç kanal elemanı bizim downstream sürecimize aktarıldı, bu da üç süreç görevi oluşturdu.

Temelde üç ayrı kez toplamaya çalıştı ve her seferinde sadece tek bir dosyası vardı, bu yüzden sadece cat tek dosya çıktıya yaptı ve aslında bunu .command.sh dosyasında da görebiliriz.

Eğer .command.sh yaparsam, burada sadece tek bir dosya adı olduğunu görebiliriz ve o çalışma dizinine yalnızca tek bir dosya staged edildi.

## 2.3. Toplama adımını iş akışına ekleyin

Bu yüzden bir şekilde Nextflow'a önceki bir süreçten tüm bu çıktıları bir araya getirmesini ve bunları bu downstream sürece üç yerine tek bir kanal elemanı olarak vermesini söylememiz gerekiyor.

Bunu _collect_ adlı bir kanal operatörü ile yapıyoruz.

Bu, Nextflow pipeline'larında her zaman göreceğiniz süper kullanışlı bir operatör. Bu burada bir kanal, bu çıktı kanalı, tıpkı yukarıda oluşturduğumuz gibi. Ve bu yüzden daha önce yaptığımız gibi ona kanal operatörleri ekleyebiliriz. Sadece nokta yapabiliriz ve sonra bu durumda, collect, parantezler.

Ve ihtiyacımız olan tek şey bu. Bu daha sonra bu kanalı bu sürece aktarılmadan önce manipüle edecek.

Ona ne olduğunu görmek istiyorsanız, onu burada da görüntüleyebiliriz. Yani burada, bu süreç çalıştırmakla hiç ilgili değil, bu yüzden bunu o süreci çalıştırdıktan sonra herhangi bir noktaya koyabilirim. Ama aynı çıktı kanalını alıyoruz ve .view ile bakıyoruz ve sonra .collect.view ile tekrar bakıyoruz.

Ve bunu çalıştırdığımızda, bize o kanalın iki farklı yapısını gösterecek, collect'ten önce ve sonra. Şimdi bunu deneyelim. Tamam, bazı çıktılar oldukça uzun olduğu için biraz zoom out yaptım, ama pipeline'ı çalıştırırsam, işe yarayıp yaramayacağını göreceğiz.

Üçüncü sürecin sadece bir kez çalışmasını umuyorum, çünkü çıktıları topluyor ve gerçekten de, collectGreetings'i bir tanesi bir olarak görebilirsiniz. Yani bu sadece bir görev çalıştırdı.

Ve sonra view ifadelerine bakarsak, öncesinin üç elemanı için üç view ifademiz var, her birinde bir dosya yolu var.

Ve sonra o collect ifadesinden sonra, bu sadece bir kez tetiklendi çünkü o kanalda tek bir eleman var. Ve şimdi üç farklı dosya yolunun bu listesine sahibiz.

Tam olarak umduğumuz şey bu. Ve umarım görebilirsiniz, bu temelde CSV dizilerinden ayrı kanal elemanlarına gitmek için yaptığımız o "map" operatörünün tersi. Şimdi ayrı kanal elemanlarını alıyor ve tek bir diziye geri koyuyoruz.

Harika, bu view ifadelerini temizleyebiliriz. Bunlara artık ihtiyacımız yok. Bir sonraki adıma geçebiliriz.

Daha fazla ilerlemeden önce ve unutmadan önce, burada yeni bir publish ifadesi ekleyeceğim. Third output. İş akışınızda buna daha anlamlı ve açıklayıcı bir şey adlandırabilirsiniz. Ve sonra bunu output bloğuna tekrar ekleyeceğim ve path 'hello_workflow' mode 'copy' diyeceğim. Böylece bu süreç tarafından oluşturulan çıktı dosyası buradaki results klasörümüze kaydedilecek.

Çalışıp çalışmadığını hızlıca kontrol etmek için. Şimdi biraz daha temiz olmalı çünkü o view ifadelerimiz yok. Ve burada yeni çıktı dosyamızı alıp almadığımızı göreceğiz. Birinin, bir görev çalıştı, collected adında yeni bir dosya aldık ve şimdi o üç kelimenin hepsine sahibiz. Harika. Sırada ne var?

## 3. Bir sürece ek parametreler aktarın

Tamam. Sonra tek bir sürece birden fazla girdinin işlenmesine bakacağız. Şimdiye kadar görebilirsiniz ki tüm süreçlerimiz girdi olarak sadece bir şey alıyor. Hepsinin girdisi altında tek bir satırı var.

Bunu Nextflow'un farklı bir batch tanımlayıcısı belirtmesine izin vererek göstereceğiz, böylece belki bu iş akışını birden çok kez çalıştırırsınız ve her seferinde farklı bir batch ID verebilirsiniz.

Sadece collectGreetings için girdiye ikinci bir satır ekleyeceğim. Ve bunu "val" çağıracağım, çünkü bu bir string. Şimdi bu bir değer, bir path değil ve buna "batch_name" diyeceğim.

Sonra bu değişkeni kullanmak için aşağıdaki script'i düzenleyeceğim ve bunu eğitim materyaliyle aynı yere koymaya çalışacağım. Yani bunu bu dosya yolunun ortasına koyuyorum COLLECTED-$\{batch_name\}-output.

Henüz tamamlamadık. Unutmayın ki çıktı dosya adlarının ne olacağını Nextflow'a söylememiz gerekiyor. Bu yüzden aynı şeyi burada da yapmamız gerekiyor: COLLECTED-$\{batch_name\}-output.txt".

Harika. Nextflow şimdi ikinci bir değişken girdisi alıyor ve bunu script'e ve çıktıya interpolasyon yapıyor.

Son bir şey, şimdi bunun nerede çağrıldığını bulmamız gerekiyor ve sürece ikinci girdiyi aktarmamız gerekiyor. Bu, başka herhangi bir dildeki bir fonksiyona girdi yapmak gibi.

Eğitimde daha önce yaptığımız gibi, burada özel "params"ı kullanacağım ve buna "params.batch" diyeceğiz, böylece -- batch CLI seçeneğimiz olabilir. Ve şimdi sürecimizin burada aktarılan iki ayrı girdisi olduğunu görebilirsiniz, virgülle ayrılmış.

Sırayı doğru yapmak gerçekten önemli, bu yüzden burada kanal ve ardından param için argüman sırası eşleşmeli. Kanal ve oradaki batch adı. Bu sadece konumsal eşleştirme.

Tamam. Bu pipeline'ı şimdi doğrudan --batch ile çalıştırabilirim, ama önce doğru şeyi yapalım ve bunu burada Params'ta tanımlayalım. Yani bunu batch'e ekleyeceğim ve sonra bunun bir string olduğunu söyleyeceğiz ve ona bir varsayılan verelim. Yani sadece ona batch diyelim. Tamam? Şimdi iş akışını çalıştırmayı deneyelim.

--batch Trio. Sanırım eğitim materyalinde bunu söylüyor, ama orada istediğimiz herhangi bir string kullanabiliriz. Ve umarım o sonuç çıktı dosyasının buraya geldiğini göreceğiz.

Ve gerçekten de, COLLECTED-trio-output - bu düzgün çalıştı. Dosyamızı yeniden adlandırdı. Ve şimdi bunun yararlı olduğunu hayal edebilirsiniz çünkü bunu replicate_two gibi farklı bir batch adıyla tekrar çalıştırırsam, o zaman bize burada farklı bir batch adı verecek.

Ve bu durumda çıktı dosyalarının üzerine yazmayacak. Yani bu güzel.

## 4. Toplayıcı adıma bir çıktı ekleyin

Tamam, şimdi sürecimize birden fazla girdi aldık. Ama birden fazla çıktı oluşturmak istersek ne olur? Buradaki örneğimiz o zaman bu süreç için kaç dosyanın toplandığını söyleyen bir rapor oluşturacağız.

Ve bunu burada bir echo komutuyla yapacağız. Yani echo diyebiliriz. There were, bunu eğitim materyalinden kopyalayacağım, böylece bunu yazarken sizi izlemek zorunda kalmazsınız.

There were $\{count_greetings\} greetings in this batch, ve bunu şimdi $\{batch_name\} adlı yeni bir dosyaya kaydet, yani aynı değişken, bunu istediğimiz kadar yeniden kullanabiliriz, report.txt.

## 4.1.1. Toplanan selamlamaların sayısını sayın

Bunu bir şekilde hesaplamamız gerekiyor. İstersek o mantığı Bash script'inde yapabiliriz, Bash mantığını kullanarak. Ancak, süreçte script bloğunun içinde ve alıntılanan bölümün üstünde olduğu sürece, doğrudan Nextflow kodu içinde betik yazabiliriz.

Buradaki hiçbir şey son render edilmiş script'e dahil edilmeyecek ve sadece bir görevi render ettiğinde Nextflow tarafından yürütülecek.

Yani burada sadece biraz mantık yapıyoruz. count_greetings adında yeni bir değişken oluşturuyoruz. Burada input files kanalını alıyoruz ve üzerinde .size() çağırıyoruz.

Tamam, o fonksiyon bana bu değişkende bir sayı verecek ve şimdi uyarımız kayboldu çünkü bu değişken tanımlanıyor.

Tamam, work dizininde o ikinci dosyayı oluşturuyoruz, ama Nextflow'a bunu bu sürecin yayınlanmış çıktısı olarak beklemesini söylememiz gerekiyor. Bunu ilk dosya için yaptığımızla tamamen aynı sözdizimi ile yapıyoruz.

"path" diyoruz çünkü, yine, istersek burada bir değişken yayınlayabiliriz "val" ile, ama "path" diyeceğiz. Ve sonra beklenen dosya adı. Burada vurgulanmadığına dikkat edin. Bunun nedeni tek tırnak kullanmam. Çift tırnak kullanmam gerekiyor.

## 4.1.2. Rapor dosyasını yayınlayın ve çıktıları adlandırın

Tamam, bu harika. Ve şimdi bu çıktılara burada tıpkı burada yaptığım gibi erişebiliriz. Ama şimdi bu farklı nesnelerin bir dizisi, bu yüzden ilkini almak için collectGreetings.out[0] veya bizim yeni raporumuz olan ikincisini almak için bir yapabilirim.

Ama bunu yapmayı pek sevmiyorum çünkü indeks saymayı karıştırmak oldukça kolay. Ve orada oturmuş çok fazla satır sayıyorsunuz ve yeni bir çıktı ekliyorsunuz ve aniden her şey bozuluyor. Yani

her şeyi isme göre referans vermek çok daha güzel. Ve bunu burada "emit" adlı özel bir anahtar ile yapabiliriz.

Buna istediğimiz her şeyi çağırabiliriz. Emit outfile, ve emit reports diyelim. Bunları tanımlarsanız ve bunu bir ya da birçok üzerinde yapabilirsiniz, size kalmış. Şimdi buraya gidebilirim ve bunun yerine dot out dot reports yapabilir ve sadece isme göre çağırabilirim, bu kodunuzu okuduğunuzda anlamak çok daha kolay ve koddaki değişikliklere karşı daha güvenli.

Burada .out.report ekledim, ama aslında yayınlanan iki farklı çıktıya sahip olmam gerekiyor. Bu yüzden collected ve report gibi daha ilginç bir şey olarak yeniden adlandıracağım ve buna ne demiştim? out file dedim, pardon. Yani burada emit adı outfile ve report. çünkü iki farklı çıktı kanalı yayınlıyoruz ve bu yüzden publish bloğunda her ikisine de referans vermemiz gerekiyor.

Sonra bunları output bloğunda da tanımlamamız gerekiyor. Yani bunu collected olarak yeniden adlandırdım, ve yine, reports için, burada biraz ayrıntılı ama yeni bir iş akışını okumaya geldiğinizde gerçekten yararlı, burada yan yana listelenen tüm farklı çıktıları, tüm farklı kanalları görmek, ve bunu daha az ayrıntılı yapmanın yolları var, buna daha sonra değineceğiz.

Tamam, deneyelim ve iş akışımızı çalıştırıp ne olacağını görelim.

Umarım şimdi temelde daha önce olduğu gibi çalışmalı. Ve burada replicate_two adında yeni bir çıktı dosyası alacağız, report. Ve işte burada. Açıldı ve toplu işteki üç selamlama var diyor, bu beklediğimiz şey, yani mükemmel.

Nextflow kodunda yürütüldüğünü kanıtlamak için buraya work dizinine girersem, bash script'inde değil, cat work/ command.sh'a gidebilirim ve burada sadece bu string'i doğrudan yankıladığını göreceksiniz. There were three greetings in this batch, ve bu değişken Nextflow tarafından interpolasyon yapıldı. .command.sh dosyasını yazmadan önce script bloğunda hesaplandı. Yani sonuç değişken hesaplaması temelde bu durumda compute ortamınızda yürütülmeden önce buna hard coded edilmiş durumda.

Ve böylece script arasındaki o ayrımı görebilirsiniz. Bloğu burada ve üzerindeki herhangi bir şey. \_Umarım bu mantıklı gelir.

## Özet ve test

Tamam, bu Hello Nextflow'un bu bölümünün sonu. Yani daha önce olduğu gibi, gidin ve testi kontrol edin. Bunu web sayfasında veya CLI'da yapın, bazı sorulardan geçin ve kapsadığımız materyalden bazılarını anladığınızı kontrol edin. Anlamadığınız herhangi bir şeyi vurgulayan bir şey olup olmadığına bakın. Çok fazla soru değil. Yapmak güzel ve kolay. Veya bunu burada web sayfasında da yapabilirsiniz.

Ve küçük bir mola verin, biraz dolaşın ve geri gelin ve Hello Nextflow'un dördüncü bölümünde bize katılın, burada modüllerden bahsedeceğiz. Çok teşekkür ederim.
