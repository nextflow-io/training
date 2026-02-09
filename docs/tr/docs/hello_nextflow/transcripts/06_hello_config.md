# Bölüm 6: Hello Config - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [eğitim materyaline](../06_hello_config.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba ve Hello Nextflow'un Altıncı Bölümüne tekrar hoş geldiniz. Bu bölüm tamamen yapılandırma dosyaları hakkında ve bu kursun son bölümü.

Nextflow özellikle iki konuda çok iyidir: tekrarlanabilirlik ve taşınabilirlik. Yapılandırma dosyaları, bunlardan ikincisinin gerçekten parladığı yerdir. Bir Nextflow pipeline'ını farklı şekillerde çalışacak ve farklı sistemlerde çalışacak şekilde yapılandırma yeteneği, altta yatan pipeline kodunu düzenlemek zorunda kalmadan.

Bu süper güç, Nextflow pipeline'larının farklı yerlerdeki diğer insanlar tarafından veya kendinizin erişiminiz olabilecek farklı altyapılarda yeniden kullanılmasına olanak tanır.

Bu, pipeline kodunu dizüstü bilgisayarınızda geliştirebileceğiniz, buluta gönderebileceğiniz, HPC'nizde çalıştırabileceğiniz anlamına gelir ve aynı pipeline kodu ve her yerde çalışır.

Bu bölümde birkaç konuyu ele alacağız. Nextflow'un yapılandırma dosyalarını nasıl işlediği, bunları nereden yüklediği, nasıl yazdığınız ve nasıl yapılandırdığınız ve pipeline'ın kendisi ile bir yapılandırma dosyasına ne gitmesi gerektiği arasındaki ayrım ile başlayacağız.

Ardından, çıktı dosyalarının nerede saklandığını değiştirmek ve ayrıca pipeline'ı farklı altyapılarda çalışacak şekilde nasıl yapılandıracağınız gibi bazı yaygın kullanım durumlarına geçeceğiz, hem farklı yazılım paketleme türlerini kullanarak hem de işleri farklı altyapılara göndererek.

## Yapılandırma dosyası hiyerarşileri

Tamam, başlayalım. Yapılandırma dosyalarını yükleme söz konusu olduğunda, Nextflow birçok farklı yerden çekebilir, bu iyi bir şey ve aynı zamanda biraz riskli bir şey olabilir çünkü bazen bir yapılandırma dosyasını nereden aldığını ve hangi sırayla yüklediğini bilmek biraz zor olabilir.

Bu yüzden bu bağlantıya tıklamanızı gerçekten tavsiye ederim, bu bizi Nextflow dokümanlarına götürür. Ve bu yapılandırma sayfasında, yapılandırmanın yüklendiği temel yerleri ve önemlisi, bu şeylerin yüklendiği sırayı listeler.

Görebileceğiniz gibi, Nextflow ana dizininize bir yapılandırma dosyası koyabilirsiniz, bu genellikle ana dizininizde ".nextflow"dur. Ve bu dosya sisteminizde her Nextflow çalıştırması tarafından her zaman yüklenecektir.

Bakılacak bir sonraki yer, pipeline'ınızın kök dizininde veya deposunda "nextflow.config" adlı bir dosyadır.

Bundan sonra, "nextflow.config" adlı başka bir dosya, ancak bu sefer Nextflow'u başlattığınız dizinde: başlatma dizini.

Son olarak, komut satırında "-c" argümanı ile yapılandırma dosyası yolları sağlayabilirsiniz ve bunu birden çok kez yapabilirsiniz. Ve bunlar belirttiğiniz sırayla uygulanır.

İsterseniz tüm bu konumlarda yapılandırma dosyaları sağlayabilirsiniz ve bunlar yinelemeli olarak yüklenecek, her biri yalnızca çakıştıkları yapılandırma kapsamlarında bir öncekinin üzerine yazacaktır.

Bu gerçekten güçlü bir sistemdir çünkü mantıklı varsayılanlar ayarlayabileceğiniz ve ardından o yapılandırmaya doğru daraldıkça giderek daha spesifik hale gelebileceğiniz anlamına gelir.

## 0. Isınma: hello-config.nf'yi çalıştırın

Tamam, bunu kapatalım ve Codespaces'imize atlayalım ve başlayalım. Daha önce olduğu gibi burayı temizledim, önceki sonuçlar dizinlerimi, Nextflow'umu ve çalışma dizinlerimi vb. kaldırdım. Hala bu dosyalar etrafta dolaşıyorsa endişelenmeyin. Sadece çok yakınlaştırdığım için aksi takdirde işler çok hızlı karışıyor.

hello-config.nf ile çalışacağız, dizinimizdeki son dosya ve bu, önceki bölümde kaldığımız yerden devam etmeli.

Yani modül dosyalarından dahil edilen dört farklı sürecimiz var. Pipeline parametrelerimiz, farklı süreçleri çağırdığımız ve kanalları birbirine bağladığımız iş akışı bloğumuz, çıktı kanallarını yayınlama ve ardından bu dosyaların nerede saklanması gerektiğini ve nasıl kopyalanması gerektiğini tanımladığımız alttaki çıktı bloğu var.

Ayrıca son bölümden zaten bir "nextflow.config" dosyamız var, Docker'ı etkinleştirdiğimiz yer ve bugün bu dosyayı oluşturacağız.

Daha önce olduğu gibi, bu ana betikteki çıktı yolunu hello config olarak değiştirdik, böylece daha önce oluşturduğunuz önceki sonuçlarla çakışmaz.

Tamam, her şeyin hala beklediğimiz gibi çalıştığını hızlıca kontrol edelim. Bir terminal açalım ve nextflow run hello-config.nf yapalım. Nextflow yükleniyor. Dört farklı sürecimizi çalıştırmalı. cowpy kullanarak güzel ascii sanatı oluşturmalı ve ardından sonuçlarımızı o dizindeki sonuçlar dosyalarımıza kaydetmeli.

Burada bu dosyaların beklediğimiz gibi göründüğünden emin olmak için hızlıca bakabilirim ve evet, işte dev Hindi'miz. Harika.

## 1.1. Varsayılan değerleri nextflow.config'e taşıyın

Şimdi yapacağımız ilk şey, betiğimizden bazı şeyleri yapılandırma dosyamıza taşımaya başlamak.

Ve bu aşamada çoğunlukla parametrelerle ilgileniyoruz. Varsayılan değerleri yapılandırma dosyasına almak istiyoruz, böylece varsayılanların ne olduğu daha net olur ve insanların bunların üzerine yazması daha kolay olur.

Bu params bloğunu buradan alıp yapılandırma dosyasına koyacağım. Ve burada biraz dikkatli olmamız gerekiyor, çünkü şu anda sözdizimi yapılandırma ve betikler arasında biraz farklı. Yapılandırma dosyası tür bildirimlerini alamaz çünkü bu parametreleri gerçekten tanımlamıyoruz, sadece onlara referans veriyoruz. Bu yüzden bunlardan kurtulacağım.

Ama bunun dışında çok benzer. Bir params bloğumuz var ve ardından farklı girdi parametrelerimiz, batch parametresi, character parametresi var.

Şimdi betiğime geri dönebilirim ve artık bu varsayılanları tanımlamama gerek yok çünkü bu değerler artık Nextflow yapılandırma dosyamda.

Ancak, parametre adlarını ve türlerini bırakıyorum, böylece Nextflow bu bilgiyi biliyor ve hala tüm tür güvenliğini ve her şeyi yapabiliyor.

Tamam. Bu dosyaları kaydedelim ve her şeyin daha önce olduğu gibi çalıştığını hızlıca kontrol edelim. Burada herhangi bir değişiklik olmamalı. Değerleri aynı tuttuk. Sadece tanımlandıkları yeri taşıdık.

Harika.

## 1.2. Çalıştırmaya özgü bir yapılandırma dosyası kullanın

Şimdi, şimdiye kadar Nextflow'u pipeline betiğimizin bulunduğu dizinden başlatıyorduk. Yani başlatma dizinimiz ve pipeline dizinimiz aynı şey gibi.

Farklı başlatma dizinleriyle farklı yapılandırma dosyalarına nasıl sahip olabileceğimizi göstermek için, şimdi yeni bir alt dizin oluşturacağız.

Yani mkdir diyeceğim ve buna tux-run diyeceğiz.

Ve sonra cd yapacağım, tux-run'a dizin değiştireceğim. Ve çalışma dizinimizin artık pipeline betiklerinin bulunduğu dizinde olmadığını unutmayın.

Tamam, yeni bir "nextflow.config" dosyası oluşturalım. Yani touch nextflow config ve VS Code'da açalım. Ayrıca buradaki kenar çubuğunda artık bu alt dizinde olduğumuzu görebilirsiniz.

Şimdi üst seviye nextflow.config'de sahip olduğumuz aynı params bloğunu alabiliriz, bunu kopyalayalım ve şimdi bu değerleri değiştirebiliriz.

İlk olarak, veri artık farklı bir göreceli yol çünkü bir alt dizindeyiz, bu yüzden bunu güncellememiz gerekiyor. Ve sonra batch'i experiment olarak değiştireceğiz ve karakteri Turkey'den tux'a değiştireceğiz.

Şimdi orada kaydet'e tıklayın ve deneyelim. Tıpkı veri gibi, şimdi betiğe ulaşmak için ../ demem gerekiyor. Yani Hello config. Ve enter'a basıyorum.

Pipeline kodu hiç değişmedi, ama şimdi iki yapılandırma seti yüklenecek ve başlatma dizini yapılandırma dosyası, pipeline nextflow.config'de ayarlanan varsayılanların üzerine yazmalı ve farklı sonuç setleri almalıyız.

Elbette, buradaki dizinimizde, tux-run içinde, bir nokta Nextflow dizini ve bir çalışma dizini olduğunu görebilirsiniz ve bunun nedeni bunların her zaman başlatma dizininizde oluşturulmasıdır. Yani bunlar önceki çalıştırmalardan sahip olduğumuz çalışma ve sonuçlardan farklı.

Şimdi, sonuçlara bakarsam, toplananımızı görebiliriz ve işte küçük tux karakterimiz. Yani bu parametrelerin düzgün bir şekilde yorumlandığını görebilirsiniz.

## 1.3. Bir parametre dosyası kullanın

Tamam. Daha önce yüklenebilecek farklı yapılandırma dosyalarından bahsederken, yapılandırmayı alabileceğimiz başka bir yeri kaçırdım.

Bunu komut satırından tire tire parametre adlarıyla gördüğümüz gibi alabiliriz, ancak aynı zamanda sadece parametrelerden oluşan bir YAML veya JSON dosyası da sağlayabiliriz.

Yapılandırma dosyası her türlü farklı kapsama sahip olabilir, ancak bu dosyalar sadece parametrelerdir ve birçok parametreyi aynı anda sağlamanın kullanıcı dostu bir yoludur ve belki de biraz daha tekrarlanabilir bir yoldur çünkü bunları dosyaya yazarsınız, bu yüzden daha sonraki bir aşamada almak kolaydır.

Öyleyse terminalimize geri dönelim ve unutmadan önce, bir dizin yukarı çıktığımızdan emin olalım, böylece artık alt dizinde değilim ve burada sahip olduğumuz test-params.yaml adlı YAML dosyasına bakacağım.

Yani sadece code test-params.yaml yaparsam, bunun sadece normal bir YAML dosyası olduğunu görebilirsiniz. Özel bir şey yok. Anahtarlar parametre adlarımız, YAML biçimlendirmesiyle burada bir iki nokta üst üste ve ardından bir değer.

Bunun Nextflow kodu olmadığını unutmayın, bu yüzden buraya değişkenler gibi şeyler koyamayız. Bunlar sadece statik değerlerdir.

Ayrıca JSON aslında YAML olarak ayrıştırıldığı için, test-params.json dosyamız da olabilir, bu çok benzer görünüyor. Sadece farklı veri formatı.

Yani burada iki farklı test dosyamız var ve biraz farklı değişkenlerimiz var.

Tamam, bunları Nextflow'a nasıl veririz? Çok basit. Nextflow run hello config yapıyoruz, daha önce olduğu gibi. Ve yapılandırma dosyası için "-c" yerine veya bu varsayılan dosya adlarını yüklemek yerine, -params-file yapıyoruz. Tek tire çünkü bu bir çekirdek Nextflow seçeneği.

Ve sonra o dosyanın yolunu geçirin. Yani "-params-file test-params.yaml" yapacağım ve bunların düzgün bir şekilde yüklenip yüklenmediğini göreceğiz.

Tamam. Çalıştı. Bu YAML dosyasında ne olduğunu kendimize hatırlatalım. Yani batch YAML olarak ayarlandı, bu yüzden böyle çağrılmalı ve bir stegosaurus olmalı. Öyleyse gidip sonuçlara bakalım. Ve COLLECTED-yaml'ımız var. Bakalım bir Stegosaurus'umuz var mı. Harika, şapka takan bir Stegosaurus. İşte bunu seviyoruz.

Yani bu gerçekten iyi çalıştı ve JSON dosyasıyla tamamen aynı. Sadece dosya uzantısını burada değiştiriyoruz ve Nextflow bunu nasıl okuyacağını biliyor.

Ve bu durumda, JSON adında bir batch'imiz olmalı ve bir kaplumbağamız olmalı. Bir bakalım. Harika. En sevdiğim CLI araçlarından biri.

## 2.1. Çıktı dizinini -output-dir ile özelleştirin

Tamam, bu çoğunlukla pipeline'a girdiler ve parametreleri değiştirmek hakkındaydı. Peki çıktılar ne olacak?

Şimdi, parametreleri kullanarak alt dizinleri değiştiriyor olsak da, tüm dosyalarımızın hala sonuçlara gittiğini fark etmiş olabilirsiniz.

Tüm dosyaların yayınlandığı temel dizini -output-dir adlı bir komut satırı bayrağıyla değiştirebiliriz. Yani nextflow run hello config yaparsam ve sonra -output-dir yaparsam ve buna "custom-outdir-cli" diyeceğiz. Yazamıyorum. Böylece bu dosyaların nereden geldiğini hatırlayalım.

Bu bir çekirdek Nextflow seçeneği ve çok yeni bir seçenek. Bu yakın zamanda eklendi ve bu, yeni dil ayrıştırıcısı ve her şeyle yapabileceğimiz şeylerden biri.

Yazmak biraz uzun. İsterseniz sadece "-o" da diyebilirsiniz. Yani sadece geri gidersem. Bunu sadece "-o" olarak kısaltabilirim, bu biraz daha basit.

Tamam. Bunu çalıştırıyoruz. Pipeline'ımızda veya bu noktada yapılandırmamızda hiçbir şeyi değiştirmedik ve umarım tüm sonuçlarımızı farklı bir üst düzey dizine kaydetmeli. Ve bunu temel olarak istediğiniz herhangi bir yola ayarlayabileceğinizi hayal edebilirsiniz.

En üstte geldi. custom-outdir-cli'miz var ve tüm dosyalar orada aynı alt dizinleri ve dosya adlarıyla tamamen aynı şekilde düzenlenmiş. Yani bu, pipeline'ın sonuçlarını nereye yayınladığını değiştirmenin gerçekten kolay bir yolu, bu sonuçların nasıl düzenlendiğini çok fazla düşünmeden.

## 2.1.2. Çıktı bloğundan sabit kodlanmış yolları kaldırın

Bu dizine bakarsam, hala Hello Config adında bir alt dizinimiz olduğunu görebiliriz, bu şimdi biraz gereksiz görünüyor.

Öyleyse betiğimizi tekrar yükleyelim ve şimdi bu alt dizini alttaki çıktı bloğundan kaldırabiliriz. Çünkü artık gerçekten ihtiyacımız yok. Yani bunu şimdi yapabiliriz, buradan silebiliriz. Ve eğer sadece buysa, bunu tamamen silebilir veya boş bir dize olarak bırakabilirsiniz. Şimdilik boş bir dize olarak bırakacağım, çünkü geri gelip gelecekte yerine farklı şeyler koyacağız. Ama alt dizinler umurunda değilse, path bildirimini tamamen kaldırmak en temiz olanıdır.

Tamam, kaydet'e basalım. Hızlıca kendimizi kontrol edelim, çünkü görünüşe göre hata yapıyoruz. Orada mevcut dosyalarla karışmayalım diye "custom-outdir-cli" dizinini kaldıracağım. Çünkü unutmayın, bir şeyler yayınladığınızda, orada zaten olan dosyaları kaldırmaz. Sadece yenilerini ekler. O komutu tekrar çalıştıralım, custom-outdir-cli.

Ve şimdi "ls custom-outdir-cli" yaparsanız, artık orada Hello Config adında bir dizin yok.

## 2.2.1. Yapılandırma dosyasında outputDir'i ayarlayın

Tamam, buradaki komut satırı bayrağı, "-o" veya "-output-dir" iyi. Ama bunun için yapılandırmada varsayılanları ayarlamaya ne dersiniz? Bunu nasıl yaparız?

"nextflow.config" dosyasını açıyorum, diğer her şeyi kapatıyorum ve bundan kurtuluyorum. Buraya yeni bir yapılandırma seçeneği ekleyebiliriz, bunu eğitim materyali web sitesinden kopyaladım ve outputDir olarak adlandırılıyor.

Herhangi bir kapsam altında değil. Parametreler veya başka bir şey altında değil. Üst düzey ve bunu bir dizeye ayarlayabiliriz. Şimdi yapılacak basit bir şey, bunu sabit kodlanmış bir dize olarak sonuçlardan başka bir şeye değiştirmektir. Ancak bu bir Nextflow yapılandırma dosyasında olduğu için, burada biraz akıllı olabiliriz ve değişkenleri de dahil edebiliriz.

Ve burada bu dizenin bir parçası olan bir params değişkeni, params.batch dahil ettiğimizi görebilirsiniz. Bu, başka yerlerden gelen değişkenleri yeniden kullanabileceğimiz anlamına gelir. Ve bu durumda, Nextflow Pipeline'ı çalıştırdığımızda --batch yaparsak, batch adının ne olduğuna göre özel yolumuzda bir alt dizin alacağız.

Tamam, bunu deneyelim ve sonuçların nasıl göründüğüne hızlıca bakalım. Yani nextflow run hello config ve --batch my_run yaparsam. Yapılandırmanın nasıl göründüğünü kendimize hatırlatalım. Yani custom-outdir-config.

Tree custom-outdir-config. Ve batch'in my_run olarak adlandırıldığını görebilirsiniz. Ve sonra my_run adında o alt dizinimiz var. Yani o dinamik dosya yolu çalıştı.

Ve sadece bu değil, artık varsayılan bir sonuçlar dizinine gitmedi ve temel dizini değiştirmek için komut satırında hiçbir şey belirtmek zorunda kalmadım. Yani varsayılan outputDir için varsayılan değeri başarıyla sıfırladık.

## 2.2.2. Batch ve süreç adlarıyla alt dizinler

Tamam, bunu biraz daha ileri götürelim. Bu, yapılandırma dosyası içinde dinamik bir değişken. Peki betiğin kendisi ne olacak? Şimdi, şimdiye kadar burada bu yollarımız vardı ve bunlar da dinamik olabilir. Yani sadece bir şeyi sabit kodlamak yerine, bazı kıvrımlı parantezler koyabilir ve dinamik bir şey koyabiliriz.

Örneğin, sayHello adlı süreçlerimiz var. sayHello.name yapabiliriz, bu sürecin bir özelliği, bu durumda sadece "sayHello" olduğu için biraz sıkıcı. Ama değişken.

Yani bu size bir fikir veriyor. Bunu buraya koyabiliriz ve convertToUpper.name, collectGreetings.name, tekrar collectGreetings.name ve cowpy diyebiliriz.

Şimdi çalıştırdığımızda, temel dizin hala custom-outdir-config olacak. Ve params.batch adında bir alt dizinde olacak, ancak bunun altındaki alt dizinler süreç adına göre düzenlenmelidir.

Bunu deneyelim ve çalışıp çalışmadığını görelim. Yani karışmayalım diye önceki dizini kaldıracağım ve tamamen aynı Nextflow Run komutunu kullanacağım.

Aynı şekilde çalışmalı. Biraz daha hızlı yapmak ve önceden hesaplanmış sonuçları kullanmak için bunların hepsinde tire resume kullanabilirdim. Şimdi, tree custom-outdir-config yaparsam, sonuçlarda değil, batch adıyla temel dizinimizde olduğunu görebilirsiniz. Ve tüm sonuçların artık süreçten sonra adlandırılmış alt dizinler içinde düzenlendiğini görebilirsiniz. Yani iki farklı yerde dinamik çıktı yolları tanımlıyoruz.

Tamam. Son şey, daha önce sahip olduğumuz ara klasörleri geri ekleyelim çünkü onlar biraz güzeldi. Intermediates.

Ve ayrıca bu params.batch hakkında biraz düşünebiliriz, belki bir pipeline geliştiricisi olarak alt dizinde bunu gerçekten sevdim, ancak pipeline'ın son kullanıcıları CLI'de "-o" veya -output-dir ayarlıyorsa, bu tüm ifadeyi tamamen üzerine yazıyor ve o alt dizini kaybediyoruz.

Yani yapabileceğimiz şey, o dinamik yolu ezilecek olan outputDir yapılandırmasından çıkarıp, ezilmeyen çıktı yoluna koymaktır.

Yani params.batch slash intermediates slash sayHello.name yapabiliriz ve tüm bunları çift tırnaklı bir dizede yapabiliriz, böylece Nextflow tarafından yorumlanır.

Şimdi kopyalayabilirim, hata. Bunları diğer süreçlere kopyalayın. Hepsini tırnak içine koymayı unutmayın. Ve intermediates'i bu belirli çıktılardan kaldırın.

Tamam? Şimdi biraz daha karmaşık görünüyor, ancak kodumuzdaki çıktı dizini yapısını gerçekten güzel organize etmeye başladığımızı görebilirsiniz.

Ve gerçekten güzel olan şey, kodda geçmeyen bu ekstra karmaşıklığın CLI'ye geçmemesidir. Yani komutumuzu -output-dir ve hangi batch değişkenleriyle çalıştırabiliriz, sadece pipeline'ı nasıl çalıştıracağımızı düşünerek ve koddaki şeyler hakkında çok fazla düşünmeden. Ve çıktı dosyalarımız gerçekten güzel bir şekilde, çok iyi organize edilmiş bir şekilde oluşturulacak, bu da temel olarak pipeline'ı kullanan insanlar için güzel.

Harika. Bunu yazarken bir hata yaptığımı fark ettim. Bakalım biri beni burada yakaladı mı.. collectGreetings.name'imiz var, yani bir şeyler ters gitti. Ve evet, elbette, yanlışlıkla bunları kıvrımlı parantezlere koymayı unuttum.

Yani unutmayın, kodunuzu yazarken dikkatli olun ve Nextflow'a neyin bir değişken olduğunu ve neyin sadece bir dize olduğunu söylediğinizden emin olun. Çünkü size tam olarak söylediğinizi yapacaktır. Ve daha fazlası değil. Tüm iyi bilgisayarlar gibi. Tamam, bu düzeltmeli.

## 2.3. Yayınlama modunu iş akışı düzeyinde ayarlayın

Bu betiğin hala sevmediğim bir kısmı var, o da mode copy'yi tekrar tekrar yazmamız gerçeği ve eğer sevmediğimiz bir şey varsa, o da kendimizi tekrar etmektir.

Yani bunu alıp yapılandırmaya taşıyarak biraz temizleyebiliriz. Ve aslında, bunu tek seferde tüm pipeline için ayarlayabiliriz. Yani bunu birden çok kez söylemek zorunda değiliz.

Yapılandırma dosyamıza gidiyoruz ve burada workflow adında yeni bir kapsamımız var. Ve ya kıvrımlı parantezler yapabiliriz ya da nokta notasyonu yapabiliriz. Hiçbir fark yaratmaz. Eğitim materyali web sitesi nokta notasyonu kullanıyor. Output diyebilirim ve karıştırıp eşleştirebiliriz, yani mode equals copy. Harika.

Ve şimdi buraya geri dönebilir ve bunları silebiliriz. Şimdi bunları yerinde bırakabilirdik. Yapılandırma temel olarak burada yazılanın üzerine yazıyor, ancak pipeline düzeyinde yapılandırmada olduğu için ve bu iki dosya birlikte gönderildiği için, bunu gerçekten iki kez yapmanın bir nedeni yok.

Tamam. Kendimizi kontrol edelim, çünkü görünüşe göre hata yapıyoruz. Bunu tekrar çalıştıralım ve dosyaları yayınlamak için copy modunu doğru kullandığımızı kontrol edelim. Yani betiği tekrar çalıştıracağız ve bu sefer sonuçları config-output-mode adlı bir dizine koyduk, dosyaların orada nasıl göründüğüne bakalım.

Ve sonra batch'e bakmak için "ls -l" yaparsam ve örneğin cowpy'ye bakabiliriz. Ve evet, bunun yumuşak bir bağlantı olmayan uygun bir dosya olduğunu görmeliyiz, yani o yapılandırma özelliği düzgün bir şekilde uygulanmış.

## 3. Bir yazılım paketleme teknolojisi seçin

Tamam. Şimdiye kadar girdilere ve çıktılara, iş akışının çalıştığı dosyalara odaklanıyorduk. Ama altyapı ne olacak? Başlangıçta Nextflow'un aynı pipeline'ı farklı bilgi işlem kurulumlarında çalıştırmanıza izin verdiğini söyledim. Peki bu nasıl görünüyor?

Bunu göstermek için, cowpy'yi çalıştırmak için Docker kullanmaktan geçeceğiz ve bunun yerine aynı şeyi yapmak için Conda kullanacağız.

Bunu çok basit bir şekilde yapabilirim. Code'a gidersem, "nextflow.config". Hatırlarsanız en üstte, daha önce docker.enabled'ı tanımladık ve son bölümde cowpy ile konteyneri kullanabilmemiz için.

Nextflow'a Docker kullanmamasını söyleyeceğim. Bunu false olarak ayarlayın. Ve Conda enabled equals true diyeceğim. Yani Nextflow'a, lütfen Conda kullan diyorum.

Şimdi sadece Conda'yı etkinleştirmek tek başına yeterli değil. Tıpkı Docker'da yaptığımız gibi, Nextflow'a ihtiyaç duyduğu yazılımı nereden alabileceğini söylemeliyiz.

Yani buradaki modüllere atlarsak. Ve cowpy betiğini açarsak. En üstte bir konteyner bildirimimiz olduğunu görebiliriz. Ve konteyner Docker tarafından kullanılır, aynı zamanda Singularity, Apptainer ve diğer birçok yazılım aracı tarafından da kullanılır.

Ancak Conda için kullanılamaz, bu yüzden "conda" adında ayrı bir bildirimimiz var ve sadece "cowpy" yazabiliriz. Ve bu, conda paket çözümlemesinin yerel conda ortamınıza göre bunu çözmenin en iyi yolunu bulmasını sağlar.

Veya eğitim materyali web sitesinin yapmasını söylediği şeyi yapmak iyi bir uygulamadır, bu da çift iki nokta üst üste notasyonuyla belirli bir conda kanalı tanımlamak ve kesinlikle yazılımın belirli bir sürümünü tanımlamaktır, böylece pipeline'ı çalıştıran her kişi aynı sürümü alır.

Konteynırların bu açıdan biraz daha üstün olduğunu unutmayın, çünkü Conda ile bir şey yüklediğinizde, o paket için tüm bağımlılıkları çözmeye devam eder ve bunlar zaman içinde değişebilir. Bağımlılık kayması denir.

Yani konteynerler, tüm yazılım bağımlılıkları yığınını tamamen kilitlediği için, A, çalışacağından ve B, tekrarlanabilir olacağından biraz daha emin olabilirsiniz.

Yani Docker veya Singularity veya Apptainer kullanabiliyorsanız, kesinlikle bunu tavsiye ederim.

Şimdi bunun güzel yanı, pipeline geliştiricisi tarafından yazılan modül dosyasının artık hem Container hem de Conda'ya sahip olması ve bu yüzden bu pipeline'ı çalıştıran kişiye, hangi yazılım paketleme çözümünü kullandığınız umurumuzda değil diyoruz. Hem Docker hem de Conda ile çalışacak ve her iki durumda da yazılımı nereden alacağınız burası.

Terminali açabiliriz ve bunu deneyelim. Yani nextflow run hello config --batch conda. Ve bu conda ile ilk kez çalıştığında, o belirli sürece geldiğinde biraz yavaş olacak, çünkü "conda install" çalıştırması gerekiyor.

Ve sadece bu bir süreç için özel bir conda ortamı oluşturuyor. Yani terminalimde sahip olduğum global conda ortamımı kullanmıyor. Sadece o bir süreç için bir tane oluşturuyor. Bu iyi çünkü iş akışınızdaki farklı süreçler arasında bağımlılık çakışmaları gibi şeylerden kaçınır. Süreçlerinizin Python'un farklı sürümlerine ihtiyaç duyan araçları varsa veya bunun gibi şeyler, bu sorun değil çünkü farklı conda ortamları kullanıyorlar.

Nextflow bu conda ortamlarını yerel olarak önbelleğe alır, size o yolun nerede olduğunu söyler, burada çalışma dizininde. Ve bu yüzden bu betiği Conda ile bir sonraki çalıştırdığımda, çok daha hızlı olacak çünkü mevcut conda ortamını bulacak ve sadece yeniden kullanacak. Ama ilk kez yaptığımızda, gidip getirmesi, çözmesi, tüm bağımlılıkları indirmesi ve her şeyi kurması gerekiyor.

Tamam, harika, çalıştı. Pipeline'ın şu anda kullanmak üzere yapılandırıldığını kendimize hatırlatabiliriz. Yapılandırma dosyasına bakarsak, şu anda benim için "custom-outdir-config" idi. O temel dizine gidersem. Ve --batch conda yaptım. İşte conda alt dizinimiz. Yani çalıştı ve işte cowpy çıktımız.

Yani cowpy'yi getirdi, conda kullanarak yerel sistemime yükledi ve süreci çalıştırdı. Ve harika olan şey, o son kullanıcı olarak, oradaki yazılım yönetimi hakkında hiç düşünmek zorunda kalmadım. Nextflow benim için halletti. Bu sistemde conda kullanmam gerektiğini söyledim. Pipeline geliştiricisi hangi paketlere ihtiyacım olduğunu söyledi. Ve Nextflow gerisini yaptı. Çok güçlü.

Aslında farklı teknolojilerin bir karışımını kullanabileceğinizi unutmayın. Yani belirli süreçler için Docker'ı etkinleştirebilirim ve diğer süreçler için conda veya bazı süreçlerin yüklü olan yerel yazılımımı kullanması gerektiğini söyleyebilirim. Bu oldukça alışılmadık, ancak mümkün ve bazı durumlarda, örneğin Docker'da paketlemesi zor olabilecek belirli yazılımlar kullanıyorsanız, bir kaçış yolunuz var.

## 4. Bir yürütme platformu seçin

Yani bu yazılım paketleme. Diğer sistemlere taşınabilirliğin diğer kısmı, işlerin gerçekte nerede çalıştığıdır. Şu anda, temel olarak dizüstü bilgisayarımda veya bu Codespaces'te çalışıyorum, bu tek bir bilgisayar. Süslü bir şey yok. Nextflow işleri elinden geldiğince paralelleştirme konusunda biraz akıllı davranıyor, ama hepsi bir sistemde.

Şimdi, bir HPC'de çalışıyorsanız, muhtemelen SLURM veya PBS veya bir şey gibi bir tür iş zamanlayıcınız vardır ve işleri o zamanlayıcıya gönderirsiniz ve tüm işleri farklı hesaplama düğümlerine dağıtır.

Çalıştırmanın başka bir yolu bulutta. Yani belki AWS Batch veya Azure Cloud veya Google kullanıyorsunuz. Ve bunların hepsi benzer bir sistemde çalışır, bir zamanlayıcınız var ve işleri gönderirsiniz ve hesaplanmak üzere farklı yerlere gönderilirler.

Şimdi biyoinformatiğe başladığımda uzak geçmişte, herkesin analiz çalıştırmak için yazılımı hesaplama altyapılarına çok bağlıydı, bu da tekrarlamayı neredeyse imkansız hale getirdi.

Ancak Nextflow'daki bu yapılandırma ayrımıyla ve Nextflow'un çok farklı hesaplama altyapısı arka uçlarıyla etkileşim kurma yeteneğiyle, pipeline'ımızı pipeline kodunu hiç değiştirmeden almak ve sadece bunu değiştirmek çok basit.

## 4.1. Farklı bir arka ucu hedefleme

Yani "nextflow.config" dosyamıza gidersek ve şimdi biraz süreç düzeyinde yapılandırma koyabiliriz. Yani en üste süreç kapsamı koyarsam ve yürütücüyü ayarlayabilirim ve burada varsayılan olan local olarak ayarlanmış.

Bunun süreç düzeyinde olduğu için, şeyleri farklı süreçlere hedefleyebileceğimizi unutmayın. Ve böylece aslında yürütücüleri sürece özgü olacak şekilde ayarlayabilir ve bazı işlerin yerel olarak, Nextflow işinin yürütüldüğü her yerde çalışabileceği hibrit bir yürütmeye sahip olabilirsiniz. Bazıları farklı HPC'ye gönderilir ve bazıları buluta gönderilebilir. İstediğiniz kadar akıllı olabilirsiniz.

Şimdi, bunu böyle bir eğitim ortamında göstermek çok zor çünkü göndermek için bir HPC'm yok. Ama yapabileceğim şey slurm yazarsam, biraz hile yapabiliriz ve bunun hissini alabilirsiniz.

Ve bu gerçekten SLURM'da çalışmaya alışkın ve SLURM başlıklarının nasıl göründüğünü bilen insanlar için en ilginç olanı. Ama nextflow run, hello config yaparsam. Başarısız olacak çünkü mevcut olmayan bir kümeye iş göndermeye çalışacak. Yani sbatch'in mevcut olmadığı hakkında bir tür hata alacağız.

Evet, yazılmış. Bu araç. Bu, bir slurm kümesine iş göndermek için kullandığınız CLI aracı. Ama yapabileceğimiz şey, buradaki çalışma dizinimize gidip komuta tıklayıp, o dizini açıp .command.run'a bakabiliriz. Ve .command.run dosyasının en üstünde, teorik bir SLURM kümesine bu iş gönderimini nasıl ele alacağını söyleyen sbatch başlıklarımız var.

Ve böylece Nextflow'un akıllı olduğunu, doğru şeyleri yaptığını görebilirsiniz. Sadece göndermek için bir kümemiz yoktu.

## 5. Hesaplama kaynak tahsislerini kontrol edin

Farklı hesaplama altyapıları arasında başka ne farklı? Başka bir şey, ne kadar kullanılabilir kaynağınız olduğu ve aslında birçok hesaplama ortamında, bir işin kaç CPU'ya ve ne kadar belleğe ihtiyaç duyduğunu belirtmeniz bir gerekliliktir.

Yine, Nextflow bunu bizim için soyutlar, böylece artık tek bir hesaplama ortamı türüne özgü değildir ve burada süreç düzeyinde kapsam yazabiliriz. CPUs equals one, memory equals two gigabytes. Pipeline'ımız çok talepkar değil, bu yüzden bu iyi olmalı.

Şimdi, burada bu sayıları tahmin ettim, ama kullanmak için mantıklı bir kaynak miktarının ne olduğunu nasıl bilirsiniz? Birçok örneğin büyük bir pipeline'ının tüm bu farklı süreçlerini incelemek ve kaynak kullanımının ne olduğunu anlamak oldukça zor bir iş.

Yani bunun için iyi bir yaklaşım, bu değerleri başlangıçta yüksek sayılara ayarlamaktır, böylece pipeline'ınız herhangi bir hata olmadan çalışır ve ardından Nextflow'dan sizin için bir kullanım raporu oluşturmasını isteyin.

Bunu yapmak çok kolay, bu yüzden bir terminale geri döneceğim. Oh, bunu local'e geri ayarlamayı hatırlamam gerekiyor ki pipeline'ım gerçekten çalışsın. Ve nextflow run diyeceğim ve -with-report komut satırı bayrağını kullanacağım.

Ve bunu boş bırakabilirim ve varsayılan bir dosya adı verecek, ama ona belirli bir dosya adı vereceğim, böylece belirli bir yere kaydedilecek.

Enter'a basın ve pipeline tamamen normal şekilde çalışır, ancak bittiğinde, benim için güzel bir HTML raporu oluşturacak.

Yani buradaki kenar çubuğunda, bu HTML dosyam var. Bunu yerel olarak çalıştırıyor olsaydım, sadece açardım. Codespaces'te olduğum için, bunun üzerine sağ tıklayacağım ve indir'e tıklayacağım, bu da yerel bilgisayarıma indirecek. Ve web tarayıcısında kolayca açabilirim.

Nextflow herhangi bir pipeline için böyle bir rapor oluşturabilir ve gerçekten güzel bilgiler var. Yani bunları her zaman kaydetmek iyi bir uygulamadır. Bize ne zaman çalıştırdığımızı, nerede çalıştırdığımızı, başarılı olup olmadığını, hangi parametrelerin kullanıldığını, CLI komutunun ne olduğunu, bunun gibi şeyleri söyler.

Ve ayrıca kaynak kullanımı hakkında bu grafikler var. Yani bize her süreç için hangi CPU çağrılarının yüzdesinin kullanıldığını bir kutu grafiği olarak söyler, çünkü her süreç için birçok görev var, böylece dağılımı görebiliriz.

Burada süreçlerimizi görebilirsiniz, cowpy ve collectGreetings sadece tek bir göreve sahipti, bu yüzden sadece tek bir çizgi. Ve hem CPU hem de bellek ve iş süresi var ve çok hızlıydılar.

Bu arada, Seqera Platform kullanıyorsanız, hiçbir şey yapmadan Platform arayüzüne yerleştirilmiş aynı grafikleri alırsınız. Yani bu bilgiyi her zaman parmaklarınızın ucunda alırsınız.

Tamam, bu raporu kullanabiliriz ve gerçek bir çalıştırmada ve pipeline'ımız tarafından kaç CPU ve ne kadar bellek kullanıldığı hakkında bir fikir edinebiliriz ve geri gelip bu değerleri yapılandırma dosyamıza koyabiliriz, böylece bir dahaki sefere belki o kadar çok talep etmeyiz. Ve biraz daha yalın olabiliriz.

Şimdi pipeline yapılandırma dosyalarını yapılandırma konusunda gerçekten akıllı olabilirsiniz. Ve yine, Seqera Platform kullanıyorsanız, bir ampul gibi görünen küçük bir düğmeye dikkat edin. Çünkü ona tıklarsanız, verilerinize, çalıştırmanıza ve pipeline'ınıza özel olarak uyarlanmış, son derece optimize edilmiş bir yapılandırma dosyası oluşturacaktır. Mümkün olan en verimli şekilde çalıştırmak için.

Ama şimdilik, aslında Nextflow'un verdiği varsayılan CPU sayısının iyi olduğunu ve sadece bir gigabayt belleğe ihtiyacımız olduğunu söyleyeceğim.

## 5.3. Belirli bir süreç için kaynak tahsisleri ayarlayın

Şimdi, gerçek hayatta, pipeline'ınızdaki tüm süreçlerin aynı gereksinimlere ihtiyaç duyması oldukça alışılmadık bir durumdur. Kaynaklar açısından çok az şeye ihtiyaç duyan ve oldukça hızlı çalışan bir raporlama aracı olan MultiQC gibi bir şeyiniz olabilir.

Ve sonra belki bir referans genomu indeksleyen veya bir hizalama yapan veya başka bir iş yapan bir şeyiniz var. Ne olduğu önemli değil, çok fazla kaynak alır. Ve bu yüzden bir zamanlayıcıya bu farklı iş gönderimleri için, farklı miktarlarda kaynak vermek istersiniz.

Bu süreç kapsamı altında, belirli süreçleri farklı şekillerde hedefleyen bir yapılandırma tanımlayabiliriz.

Burada withName kullanıyoruz, etiketleri de kullanabiliriz ve bunlar bir veya birden fazla süreci hedeflemek için bir desen kullanabilir. Burada sadece cowpy adında bir adı olan herhangi bir sürecin iki gigabayt bellek ve iki CPU'ya ayarlandığını söylüyoruz ve bu, üst düzey süreç olanından daha spesifik bir seçici olduğu için, bu durumlarda üzerine yazılır, böylece burada pipeline'ınızdaki tüm farklı süreçlerinizi gerçekten verimli hale getirmek için uyarlayan güzel bir yapılandırma dosyası oluşturabilirsiniz.

## 5.5. Kaynak limitleri ekleyin

Şimdi bir pipeline geliştiricisi olarak, muhtemelen araçları oldukça iyi biliyorum ve her şeyin mümkün olduğunca hızlı ve iyi çalışmasını istiyorum. Bu yüzden bunlardan bazıları için oldukça yüksek sayılar koyabilirim çünkü cowpy'ye 20 CPU verirsem çok daha hızlı çalışacağını biliyorum.

Bu, dizüstü bilgisayarınızda veya GitHub Actions Sürekli Entegrasyon testinde veya belki 20 CPU'nun mevcut olmadığı başka bir sistemde çalıştırmaya gittiğinizde sorun değil.

Şimdi pipeline'ı çalıştırmaya çalıştığınızda, çökecek çünkü Nextflow, bu işi hiçbir yere gönderemem diyecek. Kullanılabilir kaynaklarım yok.

Şimdi bu sert çöküşten kaçınmak için, sistemimize özgü olan, kaynak limitleri adı verilen biraz daha fazla yapılandırma ekleyebiliriz. Ve bu şöyle görünüyor. Yine süreç kapsamı altında.

Ve kaynak limitleri, temel olarak sahip olduğunuz şeyin tavanını belirtebilirsiniz. Burada bir harita ve bu harita içinde belleği, CPU'ları ve zamanı ayarlayabilirsiniz.

Şimdi olan şey, Nextflow bir süreçten bir görev gönderdiğinde, neyin istendiğine bakar ve temel olarak bunun ve bunun arasında bir minimum yapar. Yani 20 CPU talep ettiysek, ancak yalnızca dört tane mevcutsa, dört talep edecektir. Pipeline çökmez ve pipeline geliştiricisi tarafından tasarlandığı şeye mümkün olduğunca yakın kullanır.

## 6. Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profiller kullanın

Tamam. Buradaki kaynak limitlerinin sisteme özgü olabileceğini söyledim ve belki pipeline'ımda bir Nextflow yapılandırma dosyam var ve insanların bunu bir dizi farklı yerde kullanacağını biliyorum. Şimdi, herkesi her seferinde kendi Nextflow yapılandırma dosyalarını oluşturmaya zorlamak yerine, yapabileceğim şey, farklı yapılandırma ön ayarlarını yapılandırma profillerine gruplamaktır.

Burada biraz aşağı kaydıracağım ve sonra parametrelerin hemen ötesinde, çünkü buradaki yapılandırma dosyasının sırası önemli, yapılandırma dosyası sırayla yüklenir, bu yüzden bu profilleri her şeyden sonra koyacağım, böylece önceden tanımlanmış parametrelerin üzerine yazar. Ve bu profilleri eğitim materyalinden yapıştıracağım.

Yani profiller adında yeni bir üst düzey kapsam var. Burada keyfi adlarımız olabilir. Yani my_laptop ve univ_hpc'miz var. Ve burada daha önce yaptığımız aynı yapılandırma parametrelerini ayarladığımızı görebiliriz. Şimdi sadece bir profil içinde. Yani dizüstü bilgisayarımda çalıştırmak için yerel bir yürütücümüz var ve HPC'de bir SLURM kümesine gönderiyorum.

Yerel olarak Docker kullanıyorum, HPC'de conda ve HPC sisteminin çok daha yüksek kaynak limitleri var.

Şimdi pipeline'ı -profile CLI seçeneğiyle çalıştırabilirim, hangi profili kullanmak istediğimi söyleyin. Yani my_laptop kullanacağım ve Nextflow o profil kapsamı içindeki tüm yapılandırmayı uygulayacak. Yani şimdi deneyebilirim. Daha öncekiyle aynı komut. Nextflow run hello config ve dash profile yapıyorum, tek tire çünkü çekirdek Nextflow seçeneği, dash profile my_laptop.

Şimdi tüm o yapılandırma seçeneğini toplu olarak uygulayacak. Oh, ve görebilirsiniz, daha önce bunun olabileceğini söyledim, süreç gereksinimi, dört CPU istedi ve bu Codespaces örneğinde sadece ikisi var.

Yani bu, süreç kaynak limitlerini denemek için iyi bir fırsat ve dizüstü bilgisayarımda veya bu Codespaces'te sadece iki CPU'ya sahip olduğumu söyleyin. Şimdi tekrar çalıştırırsak, o gereksinimi ikiye sınırlamalı ve umarım pipeline çalışacak. Harika.

## 6.2. Test parametrelerinin bir profilini oluşturun

Bu profillerin yalnızca altyapıları hakkında yapılandırmaya sahip olması gerekmediğini unutmayın. Burada parametreler de dahil olmak üzere herhangi bir yapılandırmanın gruplamalarına sahip olabilirsiniz.

Yani insanların pipeline'larında çok sık göreceğiniz başka bir şey, normalde kullanıcı başına göndereceğiniz parametreleri içeren bir test profilidir. Ama burada, test durumlarını çalıştırmak istediğimde temel olarak farklı mantıklı varsayılanlarımız var.

Ve bu harika çünkü gerekli parametreler olabilecek tüm bu şeyleri belirtmek zorunda değilim. Aksi takdirde sadece dash profile test diyebilirim ve kutudan çıkar çıkmaz çalışacak.

Şimdi dikkat edilmesi gereken bir şey, profillerin birden fazla birleştirilebileceğidir. Yani burada profile my_laptop yapabilirim ve sonra test'i de ekleyebilirim. Profile'ı iki kez yapmıyorum. Sadece burada boşluksuz virgülle ayrılmış bir liste yapıyorum. Ve bu profilleri sırayla uygulayacak. Yani my_laptop profilinden yapılandırmayı alacak ve sonra test yapılandırmasını üstüne uygulayacak.

Gerçekten kullanışlı ve pipeline'ınızı çalıştırmayı kolaylaştırmak için burada birçok mantıklı varsayılan grup ayarlayabileceğinizi görebilirsiniz.

## 6.3. Çözümlenmiş yapılandırmayı görmek için nextflow config kullanın

Umarım, Nextflow yapılandırma çözümlemesinin güçlü olduğuna sizi ikna etmişimdir, ancak yapılandırma sağlamanın yaklaşık 20 farklı yolunu söyledikten sonra bu noktada biraz şaşkın oluyorsanız sizi suçlamam ve tüm bu farklı katmanları bir soğan kabuğu gibi verin.

Yani Nextflow için nihai çözümlenmiş yapılandırmanın ne olduğundan emin değilseniz, "nextflow config" adında bir komut olduğunu bilin ve bunu çalıştırabiliriz ve bize mevcut konumumuzdaki çözümlenmiş yapılandırmanın ne olduğunu söyleyecektir.

Yani burada çalıştırdığımda, mevcut çalışma dizinindeki "nextflow.config" dosyasını bulur ve tüm farklı yapılandırmaları işler ve bana çözümlenmiş çıktıyı verir.

Nextflow yapılandırma dosyasının profile CLI seçeneğini de alabileceğini unutmayın. Yani ona my_laptop ve test profillerinde çözümlemesini söylersem ve ayrıca my_laptop yapılandırma seçeneğinden kaynak limitlerini de uyguladığını görebilirsiniz ve ayrıca testte olan parametreleri ayarladı.

Yani bu, yapılandırma çözümlemesinin nasıl çalıştığını keşfetmenin güzel bir yolu, eğer hiç emin değilseniz.

## Özet

Tamam, bu kadar. Bu kısaca Nextflow yapılandırması. Yapılandırma ile çok şey yapabilirsiniz. Gerçekten güçlü. Ama bunlar kendinizi yapacağınız en yaygın kullanım durumlarının çoğu ve bu kavramlar tüm farklı seçeneklere uygulanır.

Kendinize bir alkış verin çünkü bu Hello Nextflow eğitim kursunun sonu. Artık umarım hem sıfırdan kendi Nextflow pipeline'ınızı yazmakta, yapılandırmakta ve çalıştırmakta hem de tüm ayrıntıları ve dikkat edilmesi gereken şeyleri biliyorsunuz.

Yapılandırma eğitim sayfasında deneyebileceğiniz bir test daha var. Yani aşağı inin ve bunu deneyin ve yapılandırma hakkındaki tüm bu bölümleri anladığınızdan emin olun.

Ve, bu eğitim kursundan sonra yapılması iyi olabilecek bazı sonraki adımlar hakkında hızlı bir özet için son videoya katılın.

Bizimle kaldığınız için teşekkürler. Aferin ve bir sonraki videoda görüşürüz.
