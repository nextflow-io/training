# Bölüm 6: Hello Config - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa sadece transkripti göstermektedir. Adım adım talimatların tamamı için [eğitim materyaline](../06_hello_config.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları sadece bilgilendirme amaçlıdır ve materyaldeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba ve Hello Nextflow'un Altıncı Bölümüne tekrar hoş geldiniz. Bu bölüm tamamen yapılandırma dosyaları ile ilgili ve bu kursun son bölümü.

Nextflow özellikle iki konuda çok iyidir: tekrarlanabilirlik ve taşınabilirlik. Yapılandırma dosyaları, bunlardan ikincisinin gerçekten parladığı yerdir. Bir Nextflow pipeline'ını, altta yatan pipeline kodunu düzenlemek zorunda kalmadan farklı şekillerde çalışacak ve farklı sistemlerde çalışacak şekilde yapılandırma yeteneği.

Bu süper güç, Nextflow pipeline'larının farklı yerlerde başka insanlar tarafından veya kendinizin erişiminiz olabilecek farklı altyapılarda yeniden kullanılmasına olanak tanır.

Bu, pipeline kodunu dizüstü bilgisayarınızda geliştirebileceğiniz, buluta gönderebileceğiniz, HPC'nizde çalıştırabileceğiniz ve aynı pipeline kodunun her yerde çalıştığı anlamına gelir.

Bu bölümde birkaç konuyu inceleyeceğiz. Nextflow'un yapılandırma dosyalarını nasıl ele aldığı, nereden yüklediği, nasıl yazdığınız ve nasıl yapılandırdığınız ve pipeline'ın kendisi ile yapılandırma dosyasına ne gitmesi gerektiği arasındaki ayrım ile başlayacağız.

Daha sonra çıktı dosyalarının nerede saklandığını değiştirmek ve pipeline'ı farklı altyapılarda çalıştırmak gibi bazı yaygın kullanım durumlarına geçeceğiz; hem farklı yazılım paketleme türlerini kullanarak hem de işleri farklı altyapılara göndererek.

## Yapılandırma dosyası hiyerarşileri

Tamam, başlayalım. Yapılandırma dosyalarını yükleme konusuna gelince, Nextflow birçok farklı yerden çekebilir; bu iyi bir şey ve aynı zamanda biraz riskli bir şey olabilir çünkü bazen yapılandırma dosyasını nereden aldığını ve hangi sırayla yüklediğini bilmek biraz zor olabilir.

Bu yüzden buraya tıklamanızı gerçekten öneriyorum, bu bizi Nextflow dokümanlarına götürür. Ve bu yapılandırma sayfasında, yapılandırmanın yüklendiği önemli yerleri ve önemlisi bunların yüklendiği sırayı listeler.

Görebileceğiniz gibi, Nextflow ana dizininize bir yapılandırma dosyası koyabilirsiniz; bu genellikle ana dizininizde ".nextflow" olur. Ve bu dosya sisteminizde her Nextflow çalıştırması tarafından her zaman yüklenecektir.

Bakılacak bir sonraki yer, pipeline repository'nizin veya dizininizin kök dizininde "nextflow.config" adlı bir dosyadır.

Bundan sonra, yine "nextflow.config" adlı başka bir dosya, ancak bu sefer Nextflow'u başlattığınız dizinde: başlatma dizini.

Son olarak, komut satırında "-c" argümanı ile yapılandırma dosyası yolları sağlayabilirsiniz ve bunu birden çok kez yapabilirsiniz. Ve belirttiğiniz sırayla uygulanırlar.

İsterseniz bu konumların tümüne yapılandırma dosyaları sağlayabilirsiniz ve bunlar tekrarlı olarak yüklenecek, her biri yalnızca çakıştıkları yapılandırma kapsamlarında bir öncekinin üzerine yazacak.

Bu gerçekten güçlü bir sistemdir çünkü mantıklı varsayılanlar belirleyebilir ve sonra yapılandırmada giderek daha spesifik hale gelebilirsiniz.

## 0. Isınma: hello-config.nf'yi çalıştırın

Tamam, bunu kapatalım ve Codespaces'e atlayalım ve başlayalım. Daha önce olduğu gibi burayı temizledim, önceki sonuç dizinlerimi, Nextflow'umu ve çalışma dizinlerimi vb. kaldırdım. Bu dosyalar hâlâ etrafta dolaşıyor olsa bile endişelenmeyin. Sadece ben çok fazla zoom yaptığım için, aksi halde işler çok çabuk karışıyor.

hello-config.nf ile çalışacağız, dizinimizdeki son dosya ve bu, önceki bölümde bıraktığımız yerden devam etmeli.

Yani modül dosyalarından dahil edilen dört farklı sürecimiz var. Pipeline parametrelerimiz, farklı süreçleri çağırdığımız ve kanalları birbirine bağladığımız workflow bloğumuz, çıktı kanallarını yayınlıyoruz ve ardından bu dosyaların nerede saklanması ve nasıl kopyalanması gerektiğini tanımladığımız en alttaki çıktı bloğu var.

Son bölümden bir "nextflow.config" dosyamız da zaten var, Docker'ı etkinleştirdiğimiz yer ve bugün bu dosyayı geliştireceğiz.

Daha önce olduğu gibi, daha önce oluşturduğunuz sonuçlarla çakışmaması için bu ana betikte çıktı yolunu hello config olarak değiştirdik.

Tamam, her şeyin hâlâ beklediğimiz gibi çalıştığını hızlıca kontrol edelim. Bir terminal açıyorum ve nextflow run hello-config.nf yapıyoruz. Nextflow yükleniyor. Dört farklı sürecimizi çalıştırmalı. Cowpy kullanarak güzel asci sanatı oluşturmalı ve sonra sonuçlarımızı o dizindeki sonuç dosyalarına kaydetmeli.

Bu dosyaların beklediğimiz gibi göründüğünden emin olmak için burada hızlıca bir göz atabilirim ve işte dev Hindi'miz. Harika.

## 1.1. Varsayılan değerleri nextflow.config'e taşıyın

Şimdi yapacağımız ilk şey, betiğimizden bazı şeyleri yapılandırma dosyamıza taşımaya başlamak.

Ve bu aşamada çoğunlukla parametrelerle ilgileniyoruz. Varsayılan değerleri yapılandırma dosyasına almak istiyoruz, böylece varsayılanların ne olduğu daha net olur ve insanların bunların üzerine yazması daha kolay olur.

Bu params bloğunu betikten alıp yapılandırma dosyasına koyacağım. Ve burada biraz dikkatli olmamız gerekiyor, çünkü şu anda yapılandırma ve betikler arasındaki sözdizimi biraz farklı. Yapılandırma dosyası tür bildirimlerini kabul edemez çünkü bu parametreleri gerçekten tanımlamıyoruz, sadece referans veriyoruz. Bu yüzden bunlardan kurtulacağım.

Ama bunun dışında çok benzer. Bir params bloğumuz var ve sonra farklı girdi parametrelerimiz, batch parametresi, character parametresi var.

Şimdi betiğime geri dönebilirim ve bu varsayılanları artık tanımlamama gerek yok çünkü bu değerler artık Nextflow config dosyamda.

Ancak, parametre adlarını ve türlerini bırakıyorum, böylece Nextflow bu bilgiyi biliyor ve yine de tüm tür güvenliğini ve her şeyi yapabiliyor.

Tamam. Bu dosyaları kaydediyoruz ve her şeyin daha önce olduğu gibi çalışıp çalışmadığını hızlıca kontrol ediyoruz. Burada herhangi bir değişiklik olmamalı. Değerleri aynı tuttuk. Sadece nerede tanımlandıklarını taşıdık.

Harika.

## 1.2. Çalıştırmaya özel bir yapılandırma dosyası kullanın

Şimdi, şimdiye kadar Nextflow'u pipeline betiğimizin olduğu aynı dizinden başlatıyorduk. Yani başlatma dizinimiz ve pipeline dizinimiz aynı şey gibi.

Farklı başlatma dizinleriyle nasıl farklı yapılandırma dosyalarına sahip olabileceğimizi göstermek için, şimdi yeni bir alt dizin oluşturacağız.

mkdir diyeceğim ve buna tux-run diyeceğiz.

Sonra cd ile dizini tux-run'a değiştireceğim. Ve artık çalışma dizinimizin artık pipeline betiklerinin bulunduğu dizinle aynı olmadığını unutmayın.

Tamam, yeni bir "nextflow.config" dosyası oluşturalım. touch nextflow config yapıyorum ve VS Code'da açalım. Kenar çubuğunda da artık bu alt dizinde olduğumuzu görebilirsiniz.

Şimdi üst düzey nextflow.config'de sahip olduğumuz aynı params bloğunu alabiliriz, bunu buraya kopyalayabiliriz ve şimdi bu değerleri değiştirebiliriz.

İlk olarak, data artık farklı bir göreceli yol çünkü bir alt dizindeyiz, bu yüzden bunu güncellememiz gerekiyor. Ve sonra batch'i experiment olarak değiştireceğiz ve character'ı Turkey'den tux'a değiştireceğiz.

Şimdi orada kaydet'e tıklayın ve deneyelim. Data'da olduğu gibi, şimdi betiğe ulaşmak için ../ demem gerekiyor. Yani Hello config. Ve enter'a basıyorum.

Pipeline kodu hiç değişmedi, ama şimdi iki yapılandırma seti yüklenecek ve başlatma dizini yapılandırma dosyası, pipeline nextflow.config'de ayarlanan varsayılanların üzerine yazmalı ve farklı sonuç setleri almalıyız.

Nitekim, buradaki dizinimizde, tux-run içinde, bir dot Nextflow dizini ve bir work dizini olduğunu görebilirsiniz ve bunun nedeni bunların her zaman başlatma dizininizde oluşturulması. Yani bunlar daha önceki çalıştırmalardan aldığımız work ve results'lardan farklı.

Şimdi, results'a bakarsam, topladığımızı görebiliriz ve işte küçük tux karakterimiz. Yani bu parametrelerin düzgün bir şekilde yorumlandığını görebilirsiniz.

## 1.3. Bir parametre dosyası kullanın

Tamam. Daha önce yüklenebilecek farklı yapılandırma dosyalarından bahsederken, yapılandırmayı alabileceğimiz başka bir yeri kaçırdım.

Gördüğümüz gibi komut satırından, tire tire parametre adlarıyla alabiliriz, ama aynı zamanda sadece parametrelerin olduğu bir YAML veya JSON dosyası da sağlayabiliriz.

Yapılandırma dosyası birçok farklı kapsam türüne sahip olabilir, ancak bu dosyalar sadece parametrelerdir ve birçok parametreyi aynı anda sağlamanın kullanıcı dostu bir yoludur ve belki de biraz daha tekrarlanabilir bir yoldur çünkü bunları dosyaya yazarsınız, bu yüzden daha sonraki bir aşamada almak kolaydır.

Terminalimize geri dönelim ve unutmadan önce, bir dizin yukarı çıktığımızdan emin olalım, böylece artık alt dizinde değilim ve burada test-params.yaml adlı sahip olduğumuz YAML dosyasına bakacağım.

code test-params.yaml yaparsam, bunun sadece normal bir YAML dosyası olduğunu görebilirsiniz. Özel bir şey yok. Anahtarlar parametre adlarımız, YAML biçimlendirmesiyle yani burada bir iki nokta üst üste, ve sonra bir değer.

Bunun Nextflow kodu olmadığını unutmayın, bu yüzden buraya değişkenler gibi şeyler koyamayız. Bunlar sadece statik değerlerdir.

Ayrıca JSON aslında YAML olarak ayrıştırıldığı için, test-params.json dosyamız da olabilir, bu da çok benzer görünüyor. Sadece farklı veri formatı.

Yani burada iki farklı test dosyamız var ve biraz farklı değişkenlerimiz var.

Tamam, bunları Nextflow'a nasıl veririz? Çok basit. Nextflow run hello config yapıyoruz, daha önceki gibi. Ve yapılandırma dosyası için "-c" yerine veya bu varsayılan dosya adlarını yüklemek yerine, -params-file yapıyoruz. Tek tire çünkü bu bir temel Nextflow seçeneği.

Ve sonra o dosyanın yolunu geçiriyoruz. Yani "-params-file test-params.yaml" yapacağım ve bunların düzgün bir şekilde yüklenip yüklenmediğini göreceğiz.

Tamam. Çalıştı. Bu YAML dosyasında ne olduğunu kendimize hatırlatalım. Batch YAML olarak ayarlandı, yani böyle çağrılmalı ve bir stegosaurus olmalı. Yani results'a bakalım. Ve COLLECTED-yaml'ımız var. Bir Stegosaurus'umuz var mı görelim. Harika, şapkalı bir Stegosaurus. İstediğimiz bu.

Yani bu gerçekten iyi çalıştı ve JSON dosyasıyla tamamen aynı. Sadece burada dosya uzantısını değiştiriyoruz ve Nextflow bunu nasıl okuyacağını biliyor.

Ve bu durumda, JSON adlı bir batch'imiz olmalı ve bir kaplumbağamız olmalı. Bir göz atalım. Harika. En sevdiğim CLI araçlarından biri.

## 2.1. Çıktı dizinini -output-dir ile özelleştirin

Tamam, bu çoğunlukla pipeline'a girişler ve parametreleri değiştirme hakkında düşünmekti. Çıktılar ne olacak?

Şimdi, parametreleri kullanarak alt dizinleri değiştiriyor olsak da, tüm dosyalarımızın hâlâ results'a gittiğini fark etmiş olabilirsiniz.

Tüm dosyaların yayınlandığı bu temel dizini -output-dir adlı bir komut satırı bayrağı ile değiştirebiliriz. Nextflow run hello config yaparsam ve sonra -output-dir yaparsam, buna "custom-outdir-cli" diyeceğiz. Yazmak zahmetli. Sadece bu dosyaların nereden geldiğini hatırlayalım.

Bu bir temel Nextflow seçeneği ve çok yeni. Bu yakın zamanda eklendi ve bu, yeni dil ayrıştırıcısı ve her şeyle yapabileceğimiz şeylerden biri.

Yazmak biraz uzun. İsterseniz sadece "-o" de diyebilirsiniz. Yani sadece geri gidersem. Bunu "-o" olarak kısaltabilirim, bu biraz daha basit.

Tamam. Bunu çalıştırıyoruz. Pipeline'ımızda veya bu noktada yapılandırmamızda hiçbir şeyi değiştirmedik ve umarım tüm sonuçlarımızı farklı bir üst düzey dizine kaydetmeli. Ve hayal edebilirsiniz ki bunu istediğiniz temelde herhangi bir yola ayarlayabilirsiniz.

En üstte geldi. Bir custom-outdir-cli'miz var ve tüm dosyalar tam olarak aynı şekilde organize edilmiş, aynı alt dizinleri ve dosya adlarıyla. Yani bu, pipeline'ın sonuçlarını nasıl organize ettikleri hakkında çok fazla düşünmeden nereye yayınladığını değiştirmenin gerçekten kolay bir yolu.

## 2.1.2. Çıktı bloğundan sabit kodlanmış yolları kaldırın

Bu dizine bakarsam, hâlâ Hello Config adlı bir alt dizinimiz olduğunu görebiliyoruz, bu şimdi biraz gereksiz görünüyor.

Betiğimizi tekrar yükleyelim ve şimdi bu alt dizini en alttaki çıktı bloğundan kaldırabiliriz. Çünkü artık buna gerçekten ihtiyacımız yok. Yani bunu şimdi yapabiliriz, buradan silebiliriz. Ve eğer sadece bu ise, bunu tamamen silebilir veya boş bir dize olarak bırakabilirsiniz. Ben şimdilik boş bir dize olarak bırakacağım, çünkü geri gelip gelecekte bunun yerine bazı farklı şeyler koyacağız. Ama alt dizinlerle ilgilenmiyorsanız, yol bildirimini tamamen kaldırmak en temiz yol.

Tamam, kaydet'e basalım. Kendimizi hızlıca kontrol edelim, çünkü görünüşe göre hata yapıyoruz. Tekrar deneyelim. Aslında "custom-outdir-cli" dizinini kaldıracağım, böylece oradaki mevcut dosyalarla karışmayız. Unutmayın, bir şeyler yayınladığınızda, zaten orada olan dosyaları kaldırmaz. Sadece yenilerini ekler. O komutu tekrar çalıştıralım, custom-outdir-cli.

Ve şimdi "ls custom-outdir-cli" yaparsanız, artık orada Hello Config adlı bir dizin yok.

## 2.2.1. Yapılandırma dosyasında outputDir'i ayarlayın

Tamam, buradaki komut satırı bayrağı, "-o" veya "-output-dir" iyi. Ama bunun için yapılandırmada varsayılanları ayarlamaya ne dersiniz? Bunu nasıl yaparız?

"nextflow.config" dosyasını açıyorum, diğer her şeyi kapatıyorum ve bundan kurtuluyorum. Burada yeni bir yapılandırma seçeneği ekleyebiliriz, bunu eğitim materyali web sitesinden kopyaladım ve adı outputDir.

Herhangi bir kapsamın altında değil. params veya başka bir şeyin altında değil. Üst düzey ve bunu bir dizeye ayarlayabiliriz. Şimdi yapılacak basit bir şey, sadece sabit kodlanmış bir dize olarak results dışında bir şeye değiştirmek. Ama bu bir Nextflow yapılandırma dosyasında olduğu için, burada biraz akıllı olabiliriz ve değişkenler de dahil edebiliriz.

Ve burada params.batch olan bir params değişkeni dahil ettiğimizi görebilirsiniz, bu bu dizenin bir parçası. Bu, başka yerlerden gelen değişkenleri yeniden kullanabileceğimiz anlamına gelir. Ve bu durumda, Nextflow Pipeline'ı çalıştırdığımızda --batch yaparsak, batch adının ne olduğuna göre özel yolumuzda bir alt dizin alacağız.

Tamam, bunu deneyelim ve sonuçların nasıl göründüğüne hızlıca bir göz atalım. Nextflow run hello config yaparsam ve --batch my_run. Yapılandırmanın nasıl göründüğünü kendimize hatırlatalım. Yani custom-outdir-config.

tree custom-outdir-config. Ve batch'in my_run olarak çağrıldığını görebilirsiniz. Ve sonra my_run adlı o alt dizinimiz var. Yani o dinamik dosya yolu çalıştı.

Ve sadece bu da değil, artık varsayılan results dizinine gitmedi ve temel dizini değiştirmek için komut satırında hiçbir şey belirtmek zorunda kalmadım. Yani varsayılan outputDir için varsayılan değeri başarıyla sıfırladık.

## 2.2.2. Batch ve süreç adlarıyla alt dizinler

Tamam, bunu biraz daha ilerletelim. Bu, yapılandırma dosyası içinde dinamik bir değişken. Peki ya betiğin kendisi? Şimdiye kadar burada bu yollarımız vardı ve bunlar da dinamik olabilir. Yani sadece bir şeyi sabit kodlamak yerine, bazı kıvırcık parantezler koyabilir ve dinamik bir şey koyabiliriz.

Örneğin, sayHello adlı süreçlerimiz var. sayHello.name yapabiliriz, bu sürecin bir özelliği, bu durumda sadece "sayHello" olması biraz sıkıcı. Ama değişken.

Yani bu size bir fikir veriyor. Bunu buraya koyabiliriz ve convertToUpper.name, collectGreetings.name, yine collectGreetings.name ve cowpy diyebiliriz.

Şimdi çalıştırdığımızda, temel dizin hâlâ custom-outdir-config olacak. Ve params.batch adlı bir alt dizinde olacak, ama bunun altındaki alt dizinler süreç adına göre organize edilmeli.

Bunu deneyelim ve çalışıp çalışmadığına bakalım. Önceki dizini kaldıracağım ki karışmayalım ve tam olarak aynı Nextflow Run komutunu kullanacağım.

Aynı şekilde çalışmalı. Biraz daha hızlı hale getirmek ve daha önce hesaplanan sonuçları kullanmak için bunların hepsinde tire resume kullanabilirdim. Şimdi, tree custom-outdir-config yaparsam, results'ta değil, özel dizinimizde, batch adıyla. Ve tüm sonuçların şimdi süreç adından sonra adlandırılmış alt dizinler içinde organize edildiğini görebilirsiniz. Yani burada dinamik çıktı yollarını tanımladığımız iki farklı yerimiz var.

Tamam. Son şey, o ara klasörleri geri ekleyelim, çünkü oldukça güzeldi. Intermediates.

Ve ayrıca bu params.batch hakkında biraz düşünebiliriz, belki bir pipeline geliştiricisi olarak alt dizinde olmasını gerçekten sevdim, ama pipeline'ın son kullanıcıları CLI'de "-o" veya -output-dir ayarlıyorsa, bu tamamen bu ifadenin tamamının üzerine yazıyor ve o alt dizini kaybediyoruz.

Yani yapabileceğimiz şey, ezilecek outputDir yapılandırmasından o dinamik yolu çıkarmak ve ezilmeyen çıktı yoluna koymak.

Yani params.batch / intermediates / sayHello.name yapabiliriz ve tüm bunları çift tırnaklı bir dizede yapabiliriz, böylece Nextflow tarafından yorumlanır.

Şimdi kopyalayabilir, hoop. Bunları diğer süreçlere kopyalayın. Hepsini tırnak içine koymayı unutmayın. Ve intermediates'i bu belirli çıktılardan kaldırın.

Tamam? Şimdi biraz daha karmaşık görünüyor, ama kodumuzda gerçekten güzel organize edilmiş çıktı dizini yapısı oluşturmaya başladığımızı görebilirsiniz.

Ve gerçekten güzel olan şey, kodda CLI'ye geçmeyen bu ekstra karmaşıklık. Yani -output-dir ve hangi batch değişkenleriyle komutumuz çalıştırabiliriz, sadece pipeline'ı nasıl çalıştıracağımızı düşünerek ve kodda ne olduğunu çok fazla düşünmeden. Ve çıktı dosyalarımız gerçekten güzel bir şekilde, çok iyi organize bir şekilde oluşturulacak, bu pipeline'ı kullanan insanlar için güzel.

Harika. Bunu yazarken bir hata yaptığımı fark ediyorum. Burada beni yakalayan oldu mu görelim. collectGreetings.name'imiz var, yani bir şeyler ters gitti. Ve evet, kesinlikle, yanlışlıkla bunları kıvırcık parantez içine koymayı unuttum.

Yani unutmayın, kodunuzu yazarken dikkatli olun ve Nextflow'a neyin bir değişken olduğunu ve neyin sadece bir dize olduğunu söylediğinizden emin olun. Çünkü size söylediğiniz şeyi aynen yapacak. Ve daha fazlası değil. Tüm iyi bilgisayarlar gibi. Tamam, bu düzeltmeli.

## 2.3. Yayınlama modunu workflow düzeyinde ayarlayın

Bu betiğin hâlâ sevmediğim bir kısmı var, o da mode copy'yi tekrar tekrar yazmamız. Ve tekrarlamaktan hoşlanmadığımız bir şey varsa, o da kendimizi tekrarlamak.

Bunu biraz temizleyebiliriz, bunu alıp yapılandırmaya taşıyarak. Ve aslında, tüm pipeline için bir kerede ayarlayabiliriz. Yani birden çok kez söylemek zorunda kalmayız.

Yapılandırma dosyamıza gidiyoruz ve burada workflow adında yeni bir kapsamımız var. Ve ya kıvırcık parantezler yapabiliriz ya da nokta notasyonu kullanabiliriz. Fark etmez. Eğitim materyali web sitesi nokta notasyonunu kullanıyor. output diyebilirim ve karıştırıp eşleştirebiliriz, yani mode equals copy. Harika.

Ve şimdi buraya geri dönebiliriz ve bunları silebiliriz. Şimdi bunları yerinde bırakabiliriz. Yapılandırma temelde burada yazılanların üzerine yazıyor, ama pipeline düzeyinde yapılandırmaya sahip olduğumuz için ve bu iki dosya birlikte gönderildiği için, gerçekten iki kez yapmanın bir anlamı yok.

Tamam. Kendimizi kontrol edelim, çünkü görünüşe göre hata yapıyoruz. Tekrar çalıştıralım ve dosyaları yayınlamak için copy modunu doğru kullanıp kullanmadığımızı kontrol edelim. Betiği tekrar çalıştıracağız ve bu sefer sonuçları config-output-mode adlı bir dizine koyduk, oradaki dosyaların nasıl göründüğüne bakalım.

Ve sonra batch'e bakmak için "ls -l" yapabilirim ve örneğin cowpy'ye bakabiliriz. Ve görmeliyiz ki, evet, bu soft link olmayan uygun bir dosya, yani bu yapılandırma özniteliği düzgün bir şekilde uygulanmış.

## 3. Bir yazılım paketleme teknolojisi seçin

Tamam. Şimdiye kadar girişlere ve çıkışlara, workflow'un çalıştığı dosyalara odaklanıyorduk. Ama altyapı ne olacak? Başta söyledim ki Nextflow aynı pipeline'ın farklı bilgi işlem kurulumlarında çalıştırılmasına izin verir. Peki bu nasıl görünüyor?

Bunu göstermek için, cowpy'yi çalıştırmak için Docker kullanmaktan Conda kullanmaya geçeceğiz, aynı şeyi yapmak için.

Bunu çok basit bir şekilde yapabilirim. code'a gidersem, "nextflow.config". Hatırlarsanız en üstte, daha önce docker.enabled'ı tanımladık, geçen bölümde böylece cowpy ile konteyner kullanabilelim.

Nextflow'a Docker kullanmamasını söyleyeceğim. Bunu false olarak ayarlayın. Ve Conda enabled equals true diyeceğim. Yani Nextflow'a, lütfen Conda kullan diyeceğim.

Şimdi sadece Conda'yı etkinleştirmek kendi başına yeterli değil. Tıpkı Docker ile yaptığımız gibi, Nextflow'a ihtiyaç duyduğu yazılımı nereden alabileceğini söylemeliyiz.

Yani buradaki modüllere atlarsak. Ve cowpy betiğini açarsak. En üstte bir container bildirimimiz olduğunu görebiliriz. Ve konteyner Docker tarafından, aynı zamanda Singularity, Apptainer ve diğer birçok yazılım aracı tarafından kullanılır.

Ama Conda için kullanılamaz, bu yüzden "conda" adında ayrı bir bildirimimiz var ve sadece "cowpy" yazabiliriz. Ve bu, çözmeyi yerel conda ortamınıza göre en iyi şekilde çözmek için conda paket çözünürlüğüne bırakacak.

Veya eğitim materyali web sitesinin yapmasını söylediği şeyi yapmak iyi bir uygulamadır, yani çift iki nokta üst üste notasyonuyla belirli bir conda kanalı tanımlamak ve kesinlikle yazılımın belirli bir sürümünü tanımlamak, böylece pipeline'ı çalıştıran her kişi aynı sürümü alacak.

Konteynırların bu açıdan biraz daha üstün olduğunu unutmayın, çünkü Conda ile bir şey yüklediğinizde, yine de o paket için tüm bağımlılıkları çözecek ve bunlar zamanla değişebilir. Buna bağımlılık kayması denir.

Yani konteynırlar, tüm yazılım bağımlılıkları yığınını en alta kadar kilitler, bu yüzden A, çalışacağından ve B, tekrarlanabilir olacağından biraz daha emin olabilirsiniz.

Yani Docker veya Singularity veya Apptainer kullanabiliyorsanız, kesinlikle bunu öneririm.

Şimdi bunun güzel yanı, pipeline geliştiricisi tarafından yazılan modül dosyasının artık hem Container hem de Conda'ya sahip olması ve bu yüzden bu pipeline'ı çalıştıran kişiye, hangi yazılım paketleme çözümünü kullandığınızı umursamıyoruz. Her iki durumda da çalışacak ve Docker ve Conda ile her iki durumda da yazılımı nereden alacağı burası diyoruz.

Terminali açabiliriz ve bunu deneyelim. Nextflow run hello config --batch conda. Ve bu ilk defa conda ile çalıştığında, o belirli sürece geldiğinde biraz yavaş olacak, çünkü "conda install" çalıştırması gerekiyor.

Ve sadece bu bir süreç için özel bir conda ortamı oluşturuyor. Yani terminalimde sahip olduğum global conda ortamımı kullanmıyor. Sadece o bir süreç için bir tane oluşturuyor. Bu iyi çünkü workflow'unuzdaki farklı süreçler arasında bağımlılık çakışmaları gibi şeylerden kaçınıyor. Süreçlerinizin farklı Python sürümlerine ihtiyaç duyan araçları varsa veya bunun gibi şeyler varsa, bu sorun değil çünkü farklı conda ortamları kullanıyorlar.

Nextflow bu conda ortamlarını yerel olarak önbelleğe alır, size bu yolun nerede olduğunu söyler, buradaki work dizininde. Ve bu yüzden bu betiği Conda ile bir sonraki çalıştırdığımda, çok daha hızlı olacak çünkü o mevcut conda ortamını bulacak ve sadece yeniden kullanacak. Ama ilk defa yaptığımızda, gidip alması, çözmesi, tüm bağımlılıkları indirmesi ve her şeyi ayarlaması gerekiyor.

Tamam, harika, çalıştı. Pipeline'ın şu anda kullanmak üzere yapılandırıldığı şeyi kendimize hatırlatabiliriz. Yapılandırma dosyasına bakarsak, şu anda benim için "custom-outdir-config" idi. O temel dizine gidersem. Ve --batch conda yaptım. conda alt dizinimiz var. Yani çalıştı ve cowpy çıktımız var.

Yani cowpy'yi aldı, conda kullanarak yerel sistemime kurdu ve süreci çalıştırdı. Ve harika olan şey, son kullanıcı olarak, orada yazılım yönetimi hakkında hiç düşünmek zorunda kalmadım. Nextflow benim için halletti. Dedim ki, bu sistemde conda kullanmam gerekiyor. Pipeline geliştiricisi hangi paketlere ihtiyacım olduğunu söyledi. Ve Nextflow gerisini yaptı. Çok güçlü.

Aslında farklı teknolojilerin bir karışımını kullanabileceğinizi unutmayın. Yani belirli süreçler için Docker'ı etkinleştirebilirim ve diğer süreçler için conda'yı veya bazı süreçlerin yerel olarak kurulu olan yazılımı kullanması gerektiğini söyleyebilirim. Bu oldukça alışılmadık, ama mümkün ve bazı durumlarda, örneğin Docker'da paketlemesi zor olabilecek belirli yazılımlar kullanıyorsanız, bir kaçış yolunuz var.

## 4. Bir yürütme platformu seçin

Yani bu yazılım paketleme. Diğer sistemlere taşınabilirliğin diğer kısmı, işlerin gerçekte nerede çalıştığıdır. Şu anda, temelde dizüstü bilgisayarımda veya bu Codespaces'te çalışıyorum, bu tek bir bilgisayar. Gösterişli bir şey yok. Nextflow işleri elinden geldiğince paralelleştirme konusunda biraz akıllı davranıyor, ama hepsi bir sistemde.

Şimdi, bir HPC'de çalışıyorsanız, muhtemelen SLURM veya PBS veya başka bir şey gibi bir tür iş zamanlayıcınız var ve bu zamanlayıcıya iş gönderirsiniz ve o tüm işleri farklı hesaplama düğümlerine dağıtır.

Çalışmanın başka bir yolu bulutta. Belki AWS Batch veya Azure Cloud veya Google kullanıyorsunuz. Ve bunların hepsi benzer bir sistemde çalışır, bir zamanlayıcınız var ve iş gönderirsiniz ve farklı yerlere hesaplanmak üzere gönderilir.

Şimdi biyoinformatiğe başladığımda uzak geçmişte, herkesin analiz çalıştırmak için yazılımı hesaplama altyapılarına çok bağlıydı, bu da çoğaltmayı neredeyse imkansız hale getiriyordu.

Ama Nextflow'daki bu yapılandırma ayrımı ile ve Nextflow'un birçok farklı hesaplama altyapısı arka uçlarıyla etkileşim kurma yeteneği ile, pipeline'ımızı pipeline kodunu hiç değiştirmeden almak ve bunu değiştirmek çok basit.

## 4.1. Farklı bir arka ucu hedefleme

"nextflow.config" dosyamıza gidersek ve şimdi süreç düzeyinde yapılandırma koyabiliriz. En üste process kapsamını koyarsam ve executor'ı ayarlayabilirim ve burada varsayılan olan local olarak ayarlanmış.

Bunun süreç düzeyinde olduğunu unutmayın, şeyleri farklı süreçlere hedefleyebiliriz. Ve bu yüzden executor'ları gerçekten süreç spesifik olacak şekilde ayarlayabilir ve hibrit yürütmeye sahip olabilirsiniz, bazı işler Nextflow işinin yürütüldüğü yerde yerel olarak çalışabilir. Bazıları farklı HPC'ye gönderilir ve bazıları buluta gönderilebilir. İstediğiniz kadar akıllı olabilirsiniz.

Şimdi, bunu böyle bir eğitim ortamında göstermek çok zor çünkü göndereceğim bir HPC'im yok. Ama yapabileceğim şey slurm yazarsam, biraz hile yapabiliriz ve bunun bir hissini alabilirsiniz.

Ve bu çoğunlukla SLURM'da çalışmaya alışık olan ve SLURM başlıklarının nasıl göründüğünü bilen insanlar için ilginç. Ama Nextflow run, hello config yaparsam. Başarısız olacak çünkü mevcut olmayan bir kümeye iş göndermeye çalışacak. Yani sbatch'in mevcut olmadığı hakkında bir tür hata alacağız.

Evet, yazıldı. Bu araç. İşleri bir slurm kümesine göndermek için kullandığınız CLI aracı bu. Ama yapabileceğimiz şey, buraya tıklayarak work dizinimize gidip o dizini açmak ve .command.run'a bakmak. Ve .command.run dosyasının en üstünde, teorik bir SLURM kümesine bu iş gönderimini nasıl ele alacağını söyleyen sbatch başlıklarımız var.

Ve Nextflow'un akıllı olduğunu görebilirsiniz, tüm doğru şeyleri yapıyor. Sadece göndereceğimiz bir kümemiz yoktu.

## 5. Hesaplama kaynak tahsislerini kontrol edin

Farklı hesaplama altyapıları arasında başka ne farklı? Başka bir şey, ne kadar kullanılabilir kaynağınız olduğu ve aslında birçok hesaplama ortamında, bir işin kaç CPU'ya ve ne kadar belleğe ihtiyaç duyduğunu belirtmeniz bir gerekliliktir.

Yine, Nextflow bunu bizim için soyutlar, böylece artık tek bir hesaplama ortamı türüne özgü değildir ve burada süreç düzeyi kapsamında yazabiliriz. CPUs equals one, memory equals two gigabytes. Pipeline'ımız çok talepkar değil, bu yüzden sorun olmamalı.

Şimdi, bu sayıları sadece tahmin ettim, ama kullanmak için mantıklı bir kaynak miktarının ne olduğunu nasıl bilirsiniz? Birçok örnekli büyük bir pipeline'ın tüm bu farklı süreçlerini inceleyip kaynak kullanımının ne olduğunu anlamak oldukça zor bir iş.

Bu yüzden bunun için iyi bir yaklaşım, bu değerleri başlamak için yüksek sayılara ayarlamak, böylece pipeline'ınız herhangi bir hata olmadan çalışır ve sonra Nextflow'dan sizin için bir kullanım raporu oluşturmasını istemek.

Bunu yapmak çok kolay, bu yüzden bir terminale geri döneceğim. Oh, bunu tekrar local olarak ayarlamayı hatırlamam gerekiyor ki pipeline'ım gerçekten çalışsın. Ve Nextflow run diyeceğim ve -with-report adlı bir komut satırı bayrağı kullanacağım.

Ve bunu boş bırakabilirim ve varsayılan bir dosya adı verecek, ama ona belirli bir dosya adı vereceğim, böylece belirli bir yere kaydedilecek.

Enter'a basın ve pipeline tam olarak normal şekilde çalışır, ama bittiğinde benim için güzel bir HTML raporu oluşturacak.

Kenar çubuğunda burada bu HTML dosyam var. Bunu yerel olarak çalıştırıyor olsaydım, sadece açardım. Codespaces'te olduğum için, üzerine sağ tıklayıp download'a tıklayacağım, bu da onu yerel bilgisayarıma indirecek. Ve web tarayıcısında kolayca açabilirim.

Nextflow herhangi bir pipeline için böyle bir rapor oluşturabilir ve gerçekten güzel bilgiler var. Yani bunları her zaman kaydetmek iyi bir uygulamadır. Bize ne zaman çalıştırdığımızı, nerede çalıştırdığımızı, başarılı olup olmadığını, hangi parametrelerin kullanıldığını, CLI komutunun ne olduğunu, bunun gibi şeyleri söyler.

Ve ayrıca kaynak kullanımı hakkında bu grafikler var. Yani her süreç için hangi yüzde CPU çağrılarının kullanıldığını burada bir kutu grafiği olarak söyler, çünkü her süreç için birçok görev var, bu yüzden dağılımı görebiliriz.

Burada süreçlerimizi görebilirsiniz, cowpy ve collectGreetings'in sadece tek bir görevi vardı, bu yüzden sadece tek bir çizgi. Ve hem CPU hem de bellek ve iş süresi var ve çok hızlıydılar.

Bu arada, Seqera Platform kullanıyorsanız, herhangi bir şey yapmadan Platform arayüzüne yerleşik aynı grafikleri alırsınız. Yani bu bilgiyi her zaman parmaklarınızın ucunda alırsınız.

Tamam, bu raporu kullanabiliriz ve gerçek bir çalıştırmada, pipeline'ımız tarafından kaç CPU ve ne kadar bellek kullanıldığı hakkında bir his edinebiliriz ve geri gelip bu değerleri yapılandırma dosyamıza koyabiliriz, böylece bir dahaki sefere belki o kadar çok talep etmeyiz. Ve biraz daha yalın olabiliriz.

Şimdi pipeline yapılandırma dosyalarını yapılandırma konusunda gerçekten akıllı olabilirsiniz. Ve yine, Seqera Platform kullanıyorsanız, ampul gibi görünen küçük bir düğmeye dikkat edin. Çünkü buna tıklarsanız, verilerinize, çalıştırmanıza ve pipeline'ınıza özel olarak optimize edilmiş bir yapılandırma dosyası oluşturacak. En verimli şekilde çalıştırmak için.

Ama şimdilik, aslında Nextflow'un verdiği varsayılan CPU sayısının iyi olduğunu ve sadece bir gigabayt belleğe ihtiyacımız olduğunu söyleyeceğim.

## 5.3. Belirli bir süreç için kaynak tahsislerini ayarlayın

Şimdi, gerçek hayatta, pipeline'ınızdaki tüm süreçlerin aynı gereksinimlere ihtiyaç duyması oldukça alışılmadık. Bir raporlama aracı olarak MultiQC gibi bir şeyiniz olabilir, bu kaynaklar açısından çok az şeye ihtiyaç duyar ve oldukça hızlı çalışır.

Ve sonra belki bir referans genomu indeksleyen veya bir hizalama yapan veya başka bir iş yapan bir şeyiniz olabilir. Ne olduğu önemli değil, çok fazla kaynak alır. Ve bu yüzden bir zamanlayıcıya bu farklı iş gönderimleri için, farklı miktarlarda kaynak vermek istersiniz.

Bu process kapsamının altında, belirli süreçleri farklı şekillerde hedefleyen bir yapılandırma tanımlayabiliriz.

Burada withName kullanıyoruz, etiketler de kullanabiliriz ve bunlar bir veya birden çok süreci hedeflemek için bir desen kullanabilir. Burada sadece cowpy adında bir adı olan herhangi bir süreç için two gigabytes memory ve two CPUs ayarlayın diyoruz ve bu üst düzey process bir yerine daha spesifik bir seçici olduğu için, bu durumlarda bunun üzerine yazılır, böylece burada güzel bir yapılandırma dosyası oluşturabilirsiniz, bu da pipeline'ınızdaki tüm farklı süreçlerinizi gerçekten verimli hale getirmek için özelleştirir.

## 5.5. Kaynak sınırları ekleyin

Şimdi bir pipeline geliştiricisi olarak, muhtemelen araçları oldukça iyi biliyorum ve her şeyin mümkün olduğunca hızlı ve iyi çalışmasını istiyorum. Bu yüzden bazıları için oldukça yüksek sayılar koyabilirim çünkü cowpy'ye 20 CPU verdiğimde çok daha hızlı çalışacağını biliyorum.

Bu, dizüstü bilgisayarınızda veya GitHub Actions Continuous Integration testinde veya belki de 20 CPU'nun mevcut olmadığı başka bir sistemde çalıştırmaya gittiğinizde sorun olmaz.

Şimdi pipeline'ı çalıştırmaya çalıştığınızda, çökecektir çünkü Nextflow, bu işi hiçbir yere gönderemem diyecektir. Mevcut kaynaklarım yok.

Şimdi bu sert çöküşten kaçınmak için, sistemimize özgü, kaynak sınırları adı verilen biraz daha yapılandırma ekleyebiliriz. Ve bu şuna benziyor. Yine process kapsamının altında.

Ve kaynak sınırları, temelde mevcut olan tavanı belirtebilirsiniz. Burada bir map var ve bu map içinde bellek, CPU'lar ve zamanı ayarlayabilirsiniz.

Şimdi olan şey, Nextflow bir süreçten bir görev gönderdiğinde, neyin istendiğine bakar ve temelde bunun ve bunun arasında bir minimum yapar. Yani 20 CPU istedik, ama sadece dört tane mevcut ise, dört tane isteyecek. Pipeline çökmez ve pipeline geliştiricisi tarafından tasarlandığı şeye mümkün olduğunca yakın kullanır.

## 6. Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profiller kullanın

Tamam. Buradaki kaynak sınırlarının sisteme özgü olabileceğini söyledim ve belki pipeline'ımda bir Nextflow config dosyam var ve insanların bunu bir dizi farklı yerde kullanacağını biliyorum. Şimdi, herkesi her seferinde kendi Nextflow config dosyasını oluşturmaya zorlamak yerine, yapabileceğim şey farklı yapılandırma ön ayarlarını config profiles içinde bir araya gruplandırmak.

Burada biraz aşağı kaydıracağım ve params'ın hemen sonrasına, çünkü buradaki yapılandırma dosyasının sırası önemli, yapılandırma dosyası sırayla yüklenir, bu yüzden bu profilleri diğer her şeyden sonra koyacağım, böylece daha önce tanımlanan parametrelerin üzerine yazacak. Ve bu profilleri eğitim materyalinden yapıştıracağım.

Yani profiles adında yeni bir üst düzey kapsam var. Burada keyfi adlara sahip olabiliriz. Yani my_laptop ve univ_hpc'miz var. Ve burada daha önce yaptığımız aynı yapılandırma parametrelerini ayarladığımızı görebiliriz. Şimdi sadece bir profile içinde. Yani my_laptop'ta çalıştırmak için local bir executor'ımız var ve HPC'de bir SLURM kümesine gönderiyorum.

Yerel olarak Docker kullanıyorum, HPC'de conda kullanıyorum ve HPC sistemi çok daha yüksek kaynak sınırlarına sahip.

Şimdi -profile CLI seçeneği ile pipeline'ı çalıştırabilirim, hangi profili kullanmak istediğimi söylerim. Yani my_laptop'ı kullanacağım ve Nextflow o profile kapsamı içindeki tüm yapılandırmayı uygulayacak. Bunu şimdi deneyebilirim. Önceki komutla aynı. Nextflow run hello config ve tire profile, tek tire çünkü temel Nextflow seçeneği, tire profile my_laptop.

Şimdi bu yapılandırma seçeneğini toplu olarak uygulayacak. Oh ve gördünüz, daha önce bunun olabileceğini söyledim, süreç gereksinimi dört CPU istedi ve bu Codespaces örneğinde sadece ikisine sahibim.

Yani bu, süreç kaynak sınırlarını denemek için iyi bir fırsat ve my_laptop'ımda veya bu Codespaces'te sadece iki CPU'ya sahip olduğumu söylemek. Şimdi tekrar çalıştırırsak, bu gereksinimi ikiye indirmeli ve umarım pipeline çalışacak. Harika.

## 6.2. Test parametrelerinin bir profilini oluşturun

Bu profillerin sadece altyapıları hakkında yapılandırmaya sahip olmaları gerektiğini unutmayın. Burada parametreler de dahil olmak üzere herhangi bir yapılandırma grubuna sahip olabilirsiniz.

Yani insanların pipeline'larında sıklıkla göreceğiniz başka bir şey, normalde kullanıcı bazında göndereceğiniz parametreleri içeren bir test profilidir. Ama burada, test durumlarını çalıştırmak istediğimde temelde farklı mantıklı varsayılanlarımız var.

Ve bu harika çünkü gerekli olan tüm bu şeyleri belirtmek zorunda değilim. Aksi takdirde sadece tire profile test diyebilirim ve kutunun dışında çalışacak.

Şimdi dikkat edilmesi gereken bir şey, profillerin birden fazla da birleştirilebileceğidir. Yani burada profile my_laptop yapabilirim ve sonra test'i de ekleyebilirim. profile'ı iki kez yapmıyorum. Sadece boşluksuz virgülle ayrılmış bir liste yapıyorum. Ve bu profilleri sırayla uygulayacak. Yani my_laptop profile'dan yapılandırmayı alacak ve sonra test yapılandırmasını üzerine uygulayacak.

Gerçekten uygun ve pipeline'ınızı çalıştırmayı kolaylaştırmak için burada birçok mantıklı varsayılan grup oluşturabileceğinizi görebilirsiniz.

## 6.3. Çözümlenmiş yapılandırmayı görmek için nextflow config kullanın

Umarım, Nextflow yapılandırma çözünürlüğünün güçlü olduğuna sizi ikna etmişimdir, ama yapılandırma sağlamanın 20 farklı yolu ve soğan kabuğu gibi tüm bu farklı katmanları söyledikten sonra biraz şaşkın oluyorsanız size hak veririm.

Bu yüzden Nextflow için nihai çözümlenmiş yapılandırmanın ne olduğundan emin değilseniz, "nextflow config" adında bir komut olduğunu bilin ve bunu çalıştırabiliriz ve mevcut konumumuzda çözümlenmiş yapılandırmanın ne olduğunu bize söyleyecek.

Bunu burada çalıştırdığımda, mevcut çalışma dizininde "nextflow.config" dosyasını bulur ve tüm farklı yapılandırmaları işler ve bana çözümlenmiş çıktıyı verir.

Nextflow config dosyasının profile CLI seçeneğini de alabileceğini unutmayın. Bu yüzden my_laptop ve test profillerinde çözümlemesini söylersem, ve burada my_laptop yapılandırma seçeneğinden kaynak sınırlarını da uyguladığını ve testte olan parametreleri de ayarladığını görebilirsiniz.

Yani bu, yapılandırma çözünürlüğünün nasıl çalıştığını keşfetmek için güzel bir yol, eğer hiç emin değilseniz.

## Tamamlama

Tamam, bu kadar. Bu kısaca Nextflow config. Yapılandırma ile çok şey yapabilirsiniz. Gerçekten güçlü. Ama bunlar bulacağınız en yaygın kullanım durumlarının çoğu ve bu kavramlar tüm farklı seçeneklere uygulanır.

Kendinizi sırtlayın çünkü bu Hello Nextflow eğitim kursunun sonu. Umarım şimdi sıfırdan kendi Nextflow pipeline'ınızı yazmaya, yapılandırmaya ve çalıştırmaya ve tüm incelikleri ve dikkat edilmesi gereken şeyleri bilmeye güveniyorsunuz.

Yapılandırma eğitim sayfasında deneyebileceğiniz bir quiz daha var. Bu yüzden aşağı inin ve deneyin ve yapılandırma hakkındaki tüm bu bölümleri anladığınızdan emin olun.

Ve bu eğitim kursundan sonra yapılması iyi olabilecek bazı sonraki adımlar hakkında hızlı bir özet için son videoda bize katılın.

Bizimle kalın için teşekkürler. Aferin ve bir sonraki videoda görüşürüz.
