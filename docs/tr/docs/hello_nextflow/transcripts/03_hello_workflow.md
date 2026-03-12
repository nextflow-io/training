# Bölüm 3: Hello Workflow - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../03_hello_workflow.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş geldiniz ve özet

Merhaba ve Hello Nextflow'un üçüncü bölümüne tekrar hoş geldiniz. Bu bölümün adı Hello Workflow ve kursun bu bölümünde pipeline veya iş akışı adını gerçekten hak etmeye başlıyoruz.

Şu ana kadarki tek süreçli basit pipeline betiğimizi alacağız ve ek süreçler eklemeye başlayıp Nextflow'un bu orkestrasyon ile pipeline boyunca veri akışını nasıl yönettiğini göreceğiz.

Code Spaces'e geri dönelim. Temiz tutmak için tüm `.nextflow*` dizinlerini, work dizinlerini ve her şeyi sildim. Kursun önceki bölümlerinden bu dosyalar hâlâ etrafta dolanıyorsa endişelenmeyin.

`hello-workflow.nf` adlı bir dosyadan çalışacağız. Daha önceki gibi, bu temelde bu noktaya kadar oluşturduğumuz betiği temsil ediyor ve bize temiz bir başlangıç noktası veriyor. Ayrıca aşağıdaki çıktıda yolun artık `hello_workflow` olduğunu görebilirsiniz. Yani yayımlanan dosyalar, `results` klasörünüzdeki farklı bir alt dizine gitmeli.

Şu ana kadar nerede olduğumuzu özetlemek gerekirse, burada tek bir girdi selamlama ve tek bir çıktı selamlama dosyası olan tek bir sürecimiz var. Ardından bir dosyaya yalnızca bir echo komutu yapan basit bir Bash betiği.

Tek bir iş akışı girdimiz var; buradaki `params` bloğu, bir path beklediğini ve varsayılanın `data/greetings.csv` olduğunu söylüyor. Bu da yukarıdaki dosya.

Sonra iş akışının kendisinde bir `main` bloğumuz var. Bir kanal oluşturuyoruz. CSV'yi satırlara ayrıştırıyoruz, ardından her dizinin ilk elemanını alıyoruz ve bu kanalı o sürece aktarıyoruz. Bu da üç görev oluşturuyor ve bu sürecin çıktılarını iş akışından yayımlıyoruz.

Son olarak, output bloğunda Nextflow'a bu dosyaları bu kanaldan `hello_workflow` adlı dizine yayımlamasını söylüyoruz. Ayrıca bu dosyaları soft link yerine kopyalamasını istiyoruz.

## 1. İş akışına ikinci bir adım ekleyin

Tamam, bu bölümde iş akışımıza ikinci bir süreç ekleyeceğiz. `sayHello` sürecinin çıktılarını alacağız ve bunları `convertToUppercase` adlı ikinci bir adımda işleyeceğiz; bu adım söz konusu dosyaların içindeki tüm harfleri büyük harfe dönüştürecek.

Bu sadece aptalca bir örnek; yine basit bir string işleme. Ama iş akışı içindeki mantığı nasıl alabileceğimizi gösteriyor.

Bunun için "tr" (translate'in kısaltması) adında bir Bash komutu kullanacağız. Sonsuza kadar var olan bir Unix komutu. Buna aşina değilseniz sizi suçlamam. Sanırım eğitimden önce hiç kullanmamıştım, ama terminalde çok hızlı bir şekilde deneyebilirsiniz. `echo 'hello world'` yaparsanız ve ardından `tr`'ye pipe ederseniz, tırnak içinde karakter aralığı belirtirsiniz: küçük harf A'dan Z'ye, ardından büyük harf A'dan Z'ye. Bu da "bu harfleri şu harflere çevir" anlamına gelir.

Enter'a bastığımda her şeyin büyük harfe çevrildiğini görebilirsiniz. İnsanlara bağırmayı seviyorsanız çok güzel.

İkinci sürecimizde kullanacağımız Bash komutu işte bu kadar basit.

## 1.2. Büyük harfe çevirme adımını Nextflow süreci olarak yazın

Betiğime geri dönersem, biraz hile yapacağım ve kodu eğitim dokümanlarından kopyalayacağım. Ama tam olarak neler olduğunu görebilirsiniz.

Burada yeni bir sürecimiz var. Buna `convertToUpper` dedik, ama istediğimiz her şeyi çağırabiliriz.

Daha önce yaptığımız gibi tek bir girdi path'imiz var. Bir value kanalı değil, bir path kanalı. Ve sonra tek bir çıktı.

Script bloğunda girdi dosyasında `cat` yapıyoruz. İstersek bunu kıvırcık parantezlere koyabiliriz; bu da o değişkeni alır. Aynı Bash komutunu pipe'da çalıştırıyoruz ve sonuçları bu dosya adıyla bir dosyaya yazıyoruz; bu da output path tarafından alınır.

Şimdi bu yeni süreçle bir şeyler yapmamız gerekiyor. İş akışının farklı mantığını oluşturduğumuz yere gidip bu ilk süreçten sonra ikinci sürecimizi çalıştıracağız. Yani `convertToUpper`, burada sürecin adı.

Bir girdi aldığı için onu tek başına çağıramayız. İlk sürecin çıktısını işlemek istiyoruz. Tıpkı bu sonuçları yayımladığımız `sayHello.out` ile yaptığımız gibi, aynı sonuçları burada girdi olarak kullanmak istiyoruz; bu yüzden bunları kopyalayıp oraya koyabiliriz.

`sayHello` sürecinin `.out`'unu istiyoruz ve Nextflow bunun buradaki dosya olan basit tek bir çıktı kaydı anlamına geldiğini biliyor. Yani bu daha sonra ikinci bir sürece girdi olarak aktarılacak.

## 1.5. İş akışı çıktı yayımlamasını ayarlayın

Tamam. Son olarak, bu ikinci sürecin sonuçlarını gerçekten kaydetmek için onları da iş akışından yayımlamamız ve ardından output bloğunda tanımlamamız gerekiyor; daha öncekiyle aynı sözdizimi. Bu yüzden bunu kopyalayıp `second outputs` diyebiliriz ya da ne derseniz deyin.

İlgilendiğimiz sürecin adını alın: `convertToUpper.out`. Ardından aşağıdaki output bloğunda bunu ekleyin ve burada aynı nitelikleri yapabiliriz. Yani bu dosyaların da `Hello Workflow` alt dizininde olmasını istiyoruz ve onları da kopyalamak istiyoruz.

Harika. Çalıştırmayı deneyelim. Terminali açıp `nextflow run hello-workflow.nf` yaparsam ne yapacağını göreceğiz. Önceki bölümlerden farklı görünüp görünmediğine bakalım.

Nextflow başlatılıyor. Dokümanlarda bunu `-resume` ile yapmayı söylüyor, ama tüm work dizinlerimi sildim; bu yüzden burada bir fark yaratmazdı. Ama silmemiş olsaydınız, o zaman da çalışırdı.

Neredeyse tamamen aynı görünüyor. Ama şimdi burada ikinci bir çıktı satırı görebilirsiniz; az önce eklediğimiz ikinci sürecin adını görebilirsiniz. Gerçekten de üç kez başarıyla çalıştığını görebilirsiniz.

Harika. Önceki work dizinlerim etrafta olsaydı ve bunu `-resume` ile yapmış olsaydım, bunlar pipeline'daki yalnızca ilk adım için önbelleğe alınmış olurdu. Çünkü bu çıktılar tamamen aynıydı; bu yüzden Nextflow bunları tekrar kullanmayı bilirdi.

Böylece gerekirse iş akışınızı adım adım yinelemeli olarak oluşturmak için `-resume`'u nasıl kullanabileceğinizi görebilirsiniz.

Tamam, işe yarayıp yaramadığını görmek için buradaki `results` dizinine bakalım. Burada birkaç dosya daha var. İlk süreçten daha önce yaptığımız gibi orijinal dosyalarımız var. Gerçekten de `upper` dosyalarımız var ve harfler tamamen büyük harf; yani işe yaramış. Görmek gerçekten güzel.

Bu work dizinlerinin içine bakmak da ilginç. Daha önce olduğu gibi buradaki hash, work dizinlerine karşılık geliyor. `ls work` yapıp genişletirsem burada farklı dosyaları göreceğiz.

İlk süreçten, buraya girdi olarak çekilen çıktı dosyasını görüyoruz. Ayrıca oluşturulan yeni çıktı dosyasını da görebiliriz.

Şimdi bunu `-la` ile tüm dosyaları listeleyecek şekilde yaparsam birkaç şey daha göreceğiz. İlk olarak, bu dosyanın aslında ilk sürece bir soft link olduğunu göreceksiniz. Dosya alanından tasarruf etmek için mümkünse bu temelde her zaman bir soft link'tir. Dosyaları burada yayımlamıyoruz; bu sadece o dosyaya birinci görevden ikinci göreve referans veriyor. Böylece her şey o çalışma dizininde kapsüllenmiş, güvenli ve diğer her şeyden izole edilmiş oluyor.

Bunun orada olması gerekiyor çünkü `.command.sh` dosyasına bakarsak, yani `cat work/b8/56*` yaparsak, buradaki dosya parçalarının göreli olduğunu görebilirsiniz. Yani aynı çalışma dizinine soft link edilmiş olan o girdi dosyasını `cat`liyor.

Her work dizini böyle görünecek. Nextflow'da baktığınızda, o work dizinine staged edilmiş tüm girdi dosyalarına sahip olacaksınız. Ardından oluşturulan herhangi bir çıktı dosyanız da olacak. Yani bu harika. Beklediğimiz gibi görünüyor.

## 2.1. Toplama komutunu tanımlayın ve terminalde test edin

Tamam, iş akışımıza geri dönelim. Yapmak istediğimiz bir sonraki adım ne?

Şimdi iki sürecimiz var ve bunlar bu bir CSV dosyasını alıyor, ayrıştırıyor ve bölüyor. Ardından bu süreçlerin her biri için üç görevimiz var ve Nextflow tüm bunların paralelleştirilmesini yönetiyor; böylece mümkün olduğunda yan yana çalışıyor.

İşleri paralel çalıştırmak için bu şekilde bölme yolu çok yaygın. Bunun tersi ise her şeyi geri toplamaktır. İş akışındaki son sürecimizle yapacağımız şey budur; burada üçüncü bir süreç olacak ve bu üç farklı çıktıyı alıp hepsini tek bir dosyada birleştirecek.

Bunu terminalde oldukça basit bir şekilde yapabiliriz; nasıl görüneceğine dair bir fikir edinmek için.

`results` klasörüne gidersem: `cd results/hello_workflow/` ve burada tüm `UPPER` dosyalarımız var. Sadece `cat` kullanabilirim; bunu o dosyanın içeriğini yazdırmak için kullanıyoruz. `cat`'e birden fazla dosya verebilirsiniz ve birini diğerinin ardından okuyacak.

Yani `UPPER-*` diyebilirim; bu bana Bash genişletmesiyle aynı üç dosya adı listesini veriyor. Ve `combined.txt` diyebilirim. Sanırım dokümanlarda tam dosya adları listeleniyor, ama aynı şeyi yapıyor.

Şimdi `cat combined.txt` kullanırsam, o üç dosyanın tümünün dosya içeriklerine sahip olduğumuzu görebiliriz.

Yani temelde bu sürecin yapacağı şey budur: önceki bir süreçten gelen tüm farklı çıktı dosyalarını tek bir süreç görevinde vermeye çalışacağız, ardından onları `cat` ile birleştirip çıktı dosyasını kaydedeceğiz.

## 2.2. Toplama adımını yapmak için yeni bir süreç oluşturun

Tamam, yeni sürecimizi ekleyelim. Bunu eğitim materyallerinden yapıştıracağım ve burada okuyucu için bu soru işaretleriyle biraz alıştırma bırakıldığını görebilirsiniz. Ama sürecin genel hatlarının temelde terminalde yaptığımız şey olduğunu görebilirsiniz: bir grup girdi dosyasının `cat`'ini yapıyoruz ve burada `collected` adlı bir çıktı dosyasına yazıyoruz; ardından çıktı yine o tek path'i bekliyor.

Bu yüzden burada bir çeşit girdiye ihtiyacımız var ve bunlar bir dizi path olacak. Yine bir girdi path kanalı tanımlıyoruz ve buna `input_files` diyelim. Şimdi bu daha önce bize burada tek bir path verdi, ama tek bir bildirim olmasına rağmen bir path aynı zamanda birden fazla dosyaya da sahip olabilir.

Bunu aşağıya kopyalayacağım çünkü bu dosyaları `cat`lemek istiyoruz. Bir dizi yazdırmak gibi bazı sorunlarımız olabileceğini düşünebilirsiniz, ancak Nextflow genellikle bu konuda oldukça mantıklı davranır. İçinde birden fazla dosya olan bir kanal verilirse, hepsini boşluk ayırıcılarla bir araya getirecektir. Yani bu bize doğru sözdizimini verecek.

Bu harika. Şimdi yeni sürecimizi bağlayalım. Workflow'a gidiyorum. Çıktıları birleştir diyeceğim; yeni sürecin adı, daha önceki gibi. Bu önceki süreci alacağım: `convertToUpper` ve `.out` yapacağım.

Harika. Deneyelim ve terminalde çalışıp çalışmadığını görelim. Birkaç dizin yukarı çıkıp Nextflow komutunu yeniden çalıştırdığımda ne olacağını göreceğiz.

İş akışı başlatıldı ve şimdi üç farklı süreç adımız olduğunu görebilirsiniz; bu harika. İlk ikisi daha önce olduğu gibi görünüyor ve üçüncü yeni olan çalışıyor; bu iyi.

Ancak burada biraz garip bir şey var. Bu çıktı dosyalarını tek bir dosyada birleştirmek istedik, ama bu sürecin bir kez değil üç kez çalıştığını görebiliriz.

Gerçekten de bu work dizinlerinden birine girip `cat work/` `collected` yaparsak göreceğiz. Burada sadece tek bir kelime var, üç değil.

Olan şu: Nextflow, önceki adımlarda olduğu gibi bu paralelleştirmeye devam etti. Bu süreç bize üç elemanlı bir kanal verdi ve bu üç kanal elemanı downstream sürecimize aktarıldı; bu da üç süreç görevi oluşturdu.

Temelde üç ayrı kez toplamaya çalıştı ve her seferinde sadece tek bir dosyası vardı; bu yüzden sadece tek dosyayı çıktıya `cat`ledi. Bunu `.command.sh` dosyasında da görebiliriz.

`.command.sh` yaparsam, burada sadece tek bir dosya adı olduğunu görebiliriz ve o çalışma dizinine yalnızca tek bir dosya staged edildi.

## 2.3. Toplama adımını iş akışına ekleyin

Bu yüzden bir şekilde Nextflow'a önceki bir süreçten tüm bu çıktıları bir araya getirmesini ve bunları bu downstream sürece üç yerine tek bir kanal elemanı olarak vermesini söylememiz gerekiyor.

Bunu _collect_ adlı bir kanal operatörü ile yapıyoruz.

Bu, Nextflow pipeline'larında her zaman göreceğiniz son derece kullanışlı bir operatör. Bu burada bir kanal; yukarıda oluşturduğumuz çıktı kanalının aynısı. Bu yüzden daha önce yaptığımız gibi ona kanal operatörleri ekleyebiliriz. Sadece nokta yapabiliriz ve ardından bu durumda `collect` ve parantezler.

İhtiyacımız olan tek şey bu. Bu daha sonra bu kanalı sürece aktarılmadan önce manipüle edecek.

Ona ne olduğunu görmek istiyorsanız, burada da görüntüleyebiliriz. Burada bu süreç çalıştırmakla hiç ilgili değil; bu yüzden bunu o süreci çalıştırdıktan sonra herhangi bir noktaya koyabilirim. Ama aynı çıktı kanalını alıyoruz ve `.view` ile bakıyoruz, ardından `.collect.view` ile tekrar bakıyoruz.

Bunu çalıştırdığımızda, bize o kanalın iki farklı yapısını gösterecek: `collect`'ten önce ve sonra. Şimdi bunu deneyelim. Tamam, bazı çıktılar oldukça uzun olduğu için biraz zoom out yaptım, ama pipeline'ı çalıştırdığımda işe yarayıp yaramayacağını göreceğiz.

Üçüncü sürecin sadece bir kez çalışmasını umuyorum çünkü çıktıları topluyor. Gerçekten de `collectGreetings`'i bir tanesi bir olarak görebilirsiniz. Yani sadece bir görev çalıştırdı.

Ardından view ifadelerine bakarsak, öncesi için üç view ifademiz var; her birinde bir dosya yolu olan üç eleman.

O `collect` ifadesinden sonra ise bu yalnızca bir kez tetiklendi çünkü o kanalda tek bir eleman var. Artık üç farklı dosya yolunun bu listesine sahibiz.

Tam olarak umduğumuz şey bu. Umarım görebilirsiniz; bu temelde CSV dizilerinden ayrı kanal elemanlarına gitmek için yaptığımız `map` operatörünün tersi. Şimdi ayrı kanal elemanlarını alıp tek bir diziye geri koyuyoruz.

Harika, bu view ifadelerini temizleyebiliriz. Bunlara artık ihtiyacımız yok. Bir sonraki adıma geçebiliriz.

Daha fazla ilerlemeden ve unutmadan önce, buraya yeni bir publish ifadesi ekleyeceğim. `Third output`. İş akışınızda buna daha anlamlı ve açıklayıcı bir şey adlandırabilirsiniz. Ardından bunu output bloğuna tekrar ekleyeceğim ve `path 'hello_workflow' mode 'copy'` diyeceğim. Böylece bu süreç tarafından oluşturulan çıktı dosyası buradaki `results` klasörümüze kaydedilecek.

Çalışıp çalışmadığını hızlıca kontrol etmek için. Şimdi biraz daha temiz olmalı çünkü o view ifadelerimiz yok. Yeni çıktı dosyamızı alıp almadığımızı göreceğiz. Bir görev çalıştı, `collected` adında yeni bir dosya aldık ve artık o üç kelimenin hepsine sahibiz. Harika. Sırada ne var?

## 3. Bir sürece ek parametreler aktarın

Tamam. Sonra tek bir sürece birden fazla girdi işlenmesine bakacağız. Şimdiye kadar tüm süreçlerimizin girdi olarak sadece bir şey aldığını görebilirsiniz. Hepsinin `input` altında tek bir satırı var.

Bunu Nextflow'un farklı bir batch tanımlayıcısı belirtmesine izin vererek göstereceğiz; böylece bu iş akışını birden çok kez çalıştırıp her seferinde farklı bir batch ID verebilirsiniz.

`collectGreetings` için girdiye sadece ikinci bir satır ekleyeceğim. Buna `val` diyeceğim çünkü bu bir string. Şimdi bu bir değer, bir path değil ve buna `batch_name` diyeceğim.

Ardından bu değişkeni kullanmak için aşağıdaki script'i düzenleyeceğim ve bunu eğitim materyaliyle aynı yere koymaya çalışacağım. Yani bunu bu dosya yolunun ortasına koyuyorum: `COLLECTED-${batch_name}-output`.

Henüz tamamlamadık. Unutmayın ki çıktı dosya adlarının ne olacağını Nextflow'a söylememiz gerekiyor. Bu yüzden aynı şeyi burada da yapmamız gerekiyor: `COLLECTED-${batch_name}-output.txt`.

Harika. Nextflow artık ikinci bir değişken girdi alıyor ve bunu script'e ve çıktıya interpolasyon yapıyor.

Son bir şey: şimdi bunun nerede çağrıldığını bulmamız ve sürece ikinci girdiyi aktarmamız gerekiyor. Bu, başka herhangi bir dildeki bir fonksiyona girdi yapmak gibi.

Eğitimde daha önce yaptığımız gibi, burada özel `params`'ı kullanacağım ve buna `params.batch` diyeceğiz; böylece `--batch` CLI seçeneğimiz olabilir. Artık sürecimizin burada virgülle ayrılmış iki ayrı girdisi olduğunu görebilirsiniz.

Sırayı doğru yapmak gerçekten önemli; bu yüzden burada kanal ve ardından param için argüman sırası eşleşmeli. Kanal ve oradaki batch adı. Bu sadece konumsal eşleştirme.

Tamam. Bu pipeline'ı şimdi doğrudan `--batch` ile çalıştırabilirim, ama önce doğru şeyi yapalım ve bunu burada `params` bloğunda tanımlayalım. Yani bunu `batch`'e ekleyeceğim, ardından bunun bir string olduğunu söyleyeceğiz ve ona bir varsayılan verelim. Sadece `batch` diyelim. Tamam? Şimdi iş akışını çalıştırmayı deneyelim.

`--batch Trio`. Sanırım eğitim materyalinde bunu söylüyor, ama orada istediğimiz herhangi bir string kullanabiliriz. Umarım o sonuç çıktı dosyasının buraya geldiğini göreceğiz.

Gerçekten de, `COLLECTED-trio-output` - bu düzgün çalıştı. Dosyamızı yeniden adlandırdı. Şimdi bunun yararlı olduğunu hayal edebilirsiniz çünkü bunu `replicate_two` gibi farklı bir batch adıyla tekrar çalıştırırsam, bize burada farklı bir batch adı verecek.

Bu durumda çıktı dosyalarının üzerine yazmayacak. Yani bu güzel.

## 4. Toplayıcı adıma bir çıktı ekleyin

Tamam, şimdi sürecimize birden fazla girdi aldık. Ama birden fazla çıktı oluşturmak istersek ne olur? Buradaki örneğimizde bu süreç için kaç dosyanın toplandığını söyleyen bir rapor oluşturacağız.

Bunu burada bir echo komutuyla yapacağız. Yani `echo` diyebiliriz. `There were`, bunu eğitim materyalinden kopyalayacağım; böylece yazarken sizi izlemek zorunda kalmazsınız.

`There were ${count_greetings} greetings in this batch`, ve bunu şimdi `${batch_name}` adlı yeni bir dosyaya kaydet; yani aynı değişken, bunu istediğimiz kadar yeniden kullanabiliriz: `report.txt`.

## 4.1.1. Toplanan selamlamaların sayısını sayın

Bunu bir şekilde hesaplamamız gerekiyor. İstersek o mantığı Bash script'inde yapabiliriz, Bash mantığını kullanarak. Ancak, süreçte script bloğunun içinde ve alıntılanan bölümün üstünde olduğu sürece, doğrudan Nextflow kodu içinde betik yazabiliriz.

Buradaki hiçbir şey son render edilmiş script'e dahil edilmeyecek ve sadece bir görevi render ettiğinde Nextflow tarafından yürütülecek.

Yani burada sadece biraz mantık yapıyoruz. `count_greetings` adında yeni bir değişken oluşturuyoruz. Burada `input_files` kanalını alıyoruz ve üzerinde `.size()` çağırıyoruz.

Tamam, o fonksiyon bana bu değişkende bir sayı verecek ve artık uyarımız kayboldu çünkü bu değişken tanımlanıyor.

Tamam, work dizininde o ikinci dosyayı oluşturuyoruz, ama Nextflow'a bunu bu sürecin yayımlanan çıktısı olarak beklemesini söylememiz gerekiyor. Bunu ilk dosya için yaptığımızla tamamen aynı sözdizimi ile yapıyoruz.

`path` diyoruz çünkü yine, istersek burada bir değişken yayımlayabiliriz `val` ile, ama `path` diyeceğiz. Ardından beklenen dosya adı. Burada vurgulanmadığına dikkat edin. Bunun nedeni tek tırnak kullanmam. Çift tırnak kullanmam gerekiyor.

## 4.1.2. Rapor dosyasını yayımlayın ve çıktıları adlandırın

Tamam, bu harika. Artık bu çıktılara burada tıpkı burada yaptığım gibi erişebiliriz. Ama şimdi bu farklı nesnelerin bir dizisi; bu yüzden ilkini almak için `collectGreetings.out[0]` veya yeni raporumuz olan ikincisini almak için `[1]` yapabilirim.

Ama bunu yapmayı pek sevmiyorum çünkü indeks saymayı karıştırmak oldukça kolay. Oturmuş çok fazla satır sayıyorsunuz, yeni bir çıktı ekliyorsunuz ve aniden her şey bozuluyor. Bu yüzden

her şeyi isme göre referans vermek çok daha güzel. Bunu burada `emit` adlı özel bir anahtar ile yapabiliriz.

Buna istediğimiz her şeyi çağırabiliriz. `emit outfile` ve `emit reports` diyelim. Bunları tanımlarsanız - bunu bir ya da birçok üzerinde yapabilirsiniz, size kalmış - artık buraya gidip bunun yerine `.out.reports` yapabilir ve sadece isme göre çağırabilirim. Bu, kodunuzu okuduğunuzda anlamak çok daha kolay ve koddaki değişikliklere karşı daha güvenli.

Buraya `.out.report` ekledim, ama aslında yayımlanan iki farklı çıktıya sahip olmam gerekiyor. Bu yüzden `collected` ve `report` gibi daha ilginç bir şey olarak yeniden adlandıracağım. Peki buna ne demiştim? `outfile` dedim, pardon. Yani burada emit adları `outfile` ve `report`. Çünkü iki farklı çıktı kanalı yayımlıyoruz ve bu yüzden publish bloğunda her ikisine de referans vermemiz gerekiyor.

Ardından bunları output bloğunda da tanımlamamız gerekiyor. Yani bunu `collected` olarak yeniden adlandırdım ve yine `reports` için. Burada biraz ayrıntılı, ama yeni bir iş akışını okumaya geldiğinizde gerçekten yararlı; burada yan yana listelenen tüm farklı çıktıları, tüm farklı kanalları görmek. Bunu daha az ayrıntılı yapmanın yolları var; buna daha sonra değineceğiz.

Tamam, deneyelim ve iş akışımızı çalıştırıp ne olacağını görelim.

Umarım şimdi temelde daha önce olduğu gibi çalışmalı. Burada `replicate_two` adında yeni bir çıktı dosyası alacağız: `report`. Ve işte burada. Açıldı ve toplu işteki üç selamlama var diyor; beklediğimiz şey bu, yani mükemmel.

Nextflow kodunda yürütüldüğünü kanıtlamak için work dizinine girersem - Bash script'inde değil - `cat work/command.sh`'a gidebilirim ve burada sadece bu string'i doğrudan yankıladığını göreceksiniz. `There were three greetings in this batch` ve bu değişken Nextflow tarafından interpolasyon yapıldı. `.command.sh` dosyasını yazmadan önce script bloğunda hesaplandı. Yani sonuç değişken hesaplaması, bu durumda compute ortamınızda yürütülmeden önce temelde buna hard coded edilmiş durumda.

Böylece script bloğu ile onun üzerindeki herhangi bir şey arasındaki o ayrımı görebilirsiniz. Umarım bu mantıklı gelir.

## Özet ve test

Tamam, Hello Nextflow'un bu bölümünün sonu. Yani daha önce olduğu gibi, gidin ve testi kontrol edin. Bunu web sayfasında veya CLI'da yapın, bazı sorulardan geçin ve kapsadığımız materyalden bazılarını anladığınızı kontrol edin. Anlamadığınız herhangi bir şeyi vurgulayan bir şey olup olmadığına bakın. Çok fazla soru değil. Yapmak güzel ve kolay. Ya da bunu burada web sayfasında da yapabilirsiniz.

Küçük bir mola verin, biraz dolaşın ve geri gelin; Hello Nextflow'un dördüncü bölümünde bize katılın. Orada modüllerden bahsedeceğiz. Çok teşekkür ederim.
