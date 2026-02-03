# Bölüm 5: Merhaba Konteynerler - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa sadece transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../05_hello_containers.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba, Merhaba Nextflow eğitim kursunun Beşinci Bölümüne hoş geldiniz.

Bu bölümün adı Merhaba Konteynerler. Nextflow'un Docker ve Singularity gibi araçlarla nasıl entegre olduğundan ve yazılım konteynerlerini kullanarak pipeline'ınızın kullanıcılarına yazılım sağlamaktan bahsedeceğiz.

Bu, insanlar pipeline'ınızı çalıştırdığında, tüm farklı araçları kendilerinin yüklemek zorunda kalmadıkları anlamına gelir. Nextflow bunu onlar için yapacaktır.

Konteynerler son derece güçlü bir teknolojidir ve tekrarlanabilirlik ve kullanım kolaylığı açısından çok önemlidir. Konteynerlerin kendilerine kısa bir giriş yaparak başlayacağız, bazı docker komutlarını manuel olarak çalıştıracağız ve sonra aynı konteynerleri Nextflow pipeline'ımıza koyacağız.

Tamam. Hadi başlayalım.

Daha önce olduğu gibi, eğitim materyalini yükleyerek başlayalım. training.nextflow.io adresine gidin. Merhaba Nextflow, Bölüm Beş, Merhaba Konteynerler.

Codespaces ortamıma atlayacağım ve burada solda hello containers dot nf'i görüyoruz.

Daha önce olduğu gibi, bu dördüncü bölümü bitirdiğimiz aynı betik, bu yüzden tanıdık gelmeli.

Girdi dosyasını ve batch adını belirtmek için komut satırı parametrelerimiz var. Üç modülümüzü dahil ediyoruz ve üç işlemi çalıştırdığımız workflow'umuz var.

## 0. Isınma: hello-containers.nf'yi Çalıştırın

Bu workflow'u tekrar çalıştırmaktan ve beklediğiniz çıktıları ürettiğini iki kez kontrol etmekten çekinmeyin. Şimdilik, aslında kapatacağım ve terminale dalacağım.

## 1. Bir konteyneri 'manuel olarak' kullanın

Bu bölüme başlamak için, konteyner teknolojisi üzerine bir özet yapacağız. Docker veya singularity ya da diğer konteyner teknolojilerine çok alışkınsanız, bunu bir tazeleme olarak düşünün veya tamamen atlamaktan çekinmeyin.

Nextflow birçok farklı konteyner teknolojisini destekler. Buna Docker, Singularity, Podman, Shifter, Charliecloud ve daha fazlası dahildir.

Bu eğitimde Docker'a odaklanacağız. Code spaces ortamında önceden yüklenmiş olarak gelir ve özellikle kendi bilgisayarınızda veya dizüstü bilgisayarınızda geliştirme yapıyorsanız en popüler konteyner teknolojilerinden biridir.

Paylaşımlı bir HPC'de akademik bir ortamda çalışıyorsanız, Singularity'nin mevcut olduğunu ve Docker'ın olmadığını görebilirsiniz. Sorun değil. Tüm kavramlar tamamen aynıdır. Birkaç manuel komut farklıdır, ancak Docker'ı anlıyorsanız, singularity'yi de anlayacaksınız.

Aslında, Singularity da Code Spaces ortamında yüklüdür. İsterseniz, aynı görevleri Docker yerine Singularity kullanarak yapmayı deneyebilirsiniz.

Peki, konteyner teknolojisi nedir? Docker'ın arkasındaki fikir, uzak bir kaynaktan bir imaj çekebilmesi, yerel makinenize çekebilmesi ve ardından bu imaja dayalı bir konteyner oluşturabilmesidir.

Bu çalışan konteyner, bilgisayarınızda çalışan bir sanal makineye biraz benzer. Ortamınızdan izole edilmiştir ve önceden paketlenmiş bir işletim sistemi ve bir dizi mevcut yazılımla birlikte gelir.

## 1.1. Konteyner imajını çekin

Önceden var olan bir imajı çekmek için ihtiyacımız olan sözdizimi "docker pull". Terminalime yazacağım, ancak şimdi oynayacak bir imaja ihtiyacımız var.

İmajları kendiniz oluşturabilirsiniz. Onları Docker Hub veya quay.io gibi genel kayıt defterlerinde bulabilirsiniz. Ancak imajları hızlı bir şekilde almanın gerçekten iyi bir yolu Seqera Containers kullanmaktır.

Bu, 2024'te oluşturduğumuz, giriş yapmadan veya herhangi bir şey olmadan kullanabileceğiniz, kullanımı ücretsiz bir topluluk hizmetidir.

seqera.io/containers adresine giderseniz veya burada üstteki containers'a tıklarsanız, size bir arama arayüzü sunulur ve Conda'da veya Python Package Index'te bulunan herhangi bir aracın adını yazabilirsiniz.

Varsayılan olarak, Bioconda ve Conda Forge kanallarını arar, ancak istiyorsanız herhangi bir Conda kanalını önek olarak ekleyebilirsiniz.

Biraz eğlenmek için cowpy kullanalım. cowpy yazacağım. Bana Python Package Index ve Conda Forge'dan sonuçlar veriyor. Konteynerime eklemek için buna tıklayacağım. İsteseydim buraya birden fazla paket ekleyebilirdim. Docker'ı seçin, linux/amd64'ü seçin ve Get Container'a tıklayın.

Bu, henüz oluşturulmadıysa imajı benim için talep üzerine oluşturur ve kopyalayabileceğim bir URL verir.

İlgileniyorsanız, view Build Details'e tıklayabilirsiniz ve bu sizi kullanılan conda ortam dosyasını ve güvenlik tarama sonuçlarıyla birlikte yapı için tam yapı günlüğünü gösteren bir sayfaya götürür.

Code spaces'e geri dönersem, artık bu konteyner adını yapıştırabilir ve enter tuşuna basabilirim.

Docker şimdi bu konteyner imajı içindeki tüm farklı katmanları indirir ve şimdi bu imajın kullanıma hazır olduğunu söyler.

## Bir Singularity imajı çekme

Singularity kullanıyorsanız, süreç temelde aynıdır. İmaj paketlerimizi seçiyoruz, cowpy'yi seçiyoruz. Şimdi Docker yerine Singularity'yi seçiyoruz ve Get Container'a tıklıyoruz. Bu bize oras:// kullanan bir imaj URL'si verir. Veya isterseniz, bu kutuyu işaretleyerek https:// kullanabilirsiniz. Bu URL'yi kopyalayın. Şimdi Code Spaces'e gidin. Aslında bu alanda Singularity ile aynı olan Apptainer yüklü, ancak birbirlerinin takma adılar. Bu yüzden apptainer pull yapacağım ve sonra ona cowpy sif diyeceğim, ama istediğiniz gibi adlandırabilirsiniz. URL'yi yapıştırın. Ve bu imajı benim için indirecek.

ls -lh yapabilir ve cowpy.sif'i görebilirim

Singularity Docker'dan farklıdır, singularity tüm imajları düz dosyalarda saklar, oysa Docker, tüm katmanları ana makinenizde ayrı ayrı tutan bir kayıt defterine sahiptir ve tüm bunları takip etmek için çalışan bir daemon'a sahiptir.

## 1.2. Cowpy'yi tek seferlik bir komut olarak çalıştırmak için konteyneri kullanın

Tamam, Docker'a geri dönelim. Artık docker run yaparak oluşturduğumuz bu imajı çalıştırmayı deneyebiliriz.

Dash dash rm yapacağım, bu sadece imajın tek seferlik bir çalıştırmasını yapar. Ve imaj URL'sini yapıştıracağım. Ve son olarak, bunu çalıştırmak istediğiniz bir komutla bitirirsiniz.

Oluşturduğumuz imajda cowpy yüklüydü, o yüzden cowpy'yi deneyelim.

İşte. Komutumuzu çalıştırdı. Yerel olarak yüklü cowpy'm yok. Çalıştırmayı denersem, mevcut olmadığını görebilirsiniz. Ancak bu komutta, Docker kullanarak çalıştırdım ve bu çıktıyı doğru bir şekilde üretti.

## 1.3. Cowpy'yi etkileşimli olarak çalıştırmak için konteyneri kullanın

İstersek bundan daha ileri gidebilir ve etkileşimli olarak bir konteyner başlatıp içeride etrafta bakınabiliriz. Yine, "docker run dash dash rm" yapıyorum. Şimdi dash it yapacağım, bu Docker'a etkileşimli bir terminal istediğimizi söyler. İmaj URL'sini tekrar yapıyorum ve bu sefer cowpy yapmak yerine, bin bash yapacağım çünkü çalıştırmak istediğimiz komut bash.

Bu bizi bu çalışan konteynere götürür ve prompt'un şimdi değiştiğini görebilirsiniz.

LS slash yaparsam buradaki dizinlerin farklı olduğunu görebilirsiniz.

Sağ tarafta GitHub Code Spaces'te çalışan ikinci bir terminal açarsam ve LS slash yaparsam, workspaces ve temp gibi dizinlerimiz olduğunu görürsünüz, oysa burada Docker'da farklı.

Bu ortam Docker içinde tamamen ayrı ve ana ortamımdan izole. Bu iyi bir şey, çünkü bu komutun yürütülmesini Docker imajına izole eder ve farklı ana sistemlerde farklı insanlar arasında tekrarlanabilir tutar.

Ana sisteminizdeki verileri Docker imajı içinde kullanmak istiyorsanız, bunu açıkça konteynere monte etmeniz gerekir.

Bunu birazdan yapacağız.

## 1.3.2. İstenen araç komut(lar)ını çalıştırın

Ama önce, cowpy'yi çalıştırıp çalıştıramayacağımızı görelim. Yine, komut artık doğrudan komut satırında mevcut ve daha karmaşık şeyler yapmaya ve argümanlar geçirmeye başlayabiliriz. Hello containers ve inek yerine, tux penguenini yapalım. Başka nelerimiz var görelim.

Cheese yapalım. Harika. Dragon ve Cow'a ne dersiniz? Oldukça iyi.

## 1.3.3. Konteynerden çıkın

Tamam. Bu konteynerde hiç verim olmadığı için daha fazla bir şey yapamam. O halde bu çalışan imajdan çıkalım ve konteynere bazı verileri monte edip edemeyeceğimizi görelim. Bunu control D yaparak veya exit yazarak yapabilirim. Tamam, şimdi normal GitHub code space'ime geri döndüm.

## 1.3.4. Verileri konteynere monte edin

Docker konteynerine bazı verileri monte etmek için dash V kullanmam gerekiyor. Bu yüzden önceki docker komutumu alacağım, başa gidip dash v yapacağım. Mevcut yerel çalışma dizini için "." yapacağım ve sonra bunun ana dizinde nereye monte edilmesi gerektiğini söylemek için bir iki nokta yapıp slash data yapacağım. Yani bu özel dizini konteynerdeki slash data konumuna monte ediyor.

Şimdi LS slash yaparsam data adında yeni bir dizin görebiliriz ve LS data yaparsam, burada yan çubukta sahip olduğumuz tüm dosyaları görebilirsiniz. Harika.

## 1.3.5. Monte edilen verileri kullanın

Şimdi Docker imajı içinde ana sistemde bulunan bazı dosyaları kullanmaya başlayabiliriz. Yani cat data greetings csv diyebilirim. Hatırlarsanız, bu daha önceki farklı selamlamalarımızın bulunduğu CSV dosyamız ve bunu cowpy'ye pipe edebilirim. Harika. Şimdi bir yerlere varıyoruz.

Tamam. Docker'ı etkileşimli olarak çalıştırmak için bu kadar yeter. Umarım artık Docker'ın kabaca ne olduğu ve hem tek seferlik bir şekilde komut çalıştırmak hem de bir imajı etkileşimli olarak kullanmak için nasıl kullanılacağı konusunda bir fikriniz olmuştur. Singularity kullanıyorsanız, komutların hepsi çok benzerdir, sadece apptainer exec veya apptainer run veya singularity exec veya singularity run gibi şeyler yaparsınız.

## 2. Nextflow'da konteynerler kullanın

Sonra Nextflow workflow'umuza geri döneceğiz ve bu teknolojiyi Nextflow pipeline'ı içinde nasıl kullanacağımızı göreceğiz.

Hadi terminali kapatalım ve Hello Containers'ı tekrar açalım.

## 2.1. Bir cowpy modülü yazın

Cowpy örneğimize bağlı kalmak için, workflow'umuzda cowpy kullanan yeni bir işlem oluşturalım. Modüllere gidelim, yeni bir dosya oluşturalım ve ona cowpy nf diyelim. Şimdi biraz hile yapacağım ve eğitim materyalinden kodu kopyalayıp kaydet'e basacağım. Ve bir bakalım.

Yani bu basit bir işlem. Umarım artık bir işlemin yapı taşlarının nasıl göründüğünü anlıyorsunuzdur. publishDir'imiz tekrar var, results'a gidiyor. İki girdimiz var, bir girdi dosyası ve character adında bir string. Bir çıktımız var cowpy input file ve tam olarak Docker imajımızda bir saniye önce manuel olarak çalıştırdığımıza benzeyen bir betiğimiz var: bir dosyayı yazdırmak için cat, bunu cowpy'ye pipe ederek, hangi tür cowpy karakterini kullanmak istediğimizi söyleyerek ve bunu burada çıktı olarak geçirdiğimiz çıktı dosyasına çıktılayarak.

## 2.2. Workflow'a cowpy ekleyin

Tamam, workflow'umuza geri dönelim, bu yeni işlemi içe aktaralım. Yani modules cowpy nf'den cowpy. Hangi karakteri istediğimizi belirtebilmemiz için yeni bir parametre oluşturalım. Varsayılan olarak Turkey diyelim. Ve sonra workflow'un sonunda bu yeni işlemi çağıralım,

cowpy. Ve burada Collect Greetings'ten çıktıyı kullanalım. Yani collect greetings out, burada out file. Ve sonra ihtiyacımız olan ikinci bir argüman var, bu da az önce yaptığımız yeni params. params dot character.

## 2.2.4. Çalıştığını doğrulamak için workflow'u çalıştırın

Tamam, yeni işlemimizin çalışıp çalışmadığını görelim. Nextflow run hello containers. Bu, ilk üç işlemi çalıştırmalı ve sonra sonunda cowpy'yi çalıştırmaya çalışmalı.

Bir hata aldık. Burada söylediği şey, cowpy'de bir hata oldu ve 127 çıkış durumu vardı ve elbette, komut sh cowpy komutu bulunamadı.

Nextflow'a cowpy için mevcut bir Docker imajımız olduğunu söylemedik, bu yüzden onu ana sistemimizde çalıştırmaya çalıştı ve ana sistemimizde yüklü cowpy'miz yok, bu yüzden bir hata tetikledi.

## 2.3. Çalıştırmak için bir konteyner kullanın

Yapmamız gereken şey, Nextflow'a mevcut bir konteynerimiz olduğunu söylememiz. Cowpy işlemimize gidelim ve işlemin üstüne container adında yeni bir yönerge ekleyeceğiz.

Sonra imajımızı buluyoruz, URL'yi kopyalıyoruz ve bunu bir string'e koyuyoruz.

Bu tek başına yeterli değil çünkü bir X Flow pipeline'ı yazılımı belirtmenin birkaç yolu olabilir. Örneğin conda conda-forge cowpy da yapabilirim. Ve Nextflow'un hangi teknolojiyi kullanmak istediğinizi bilmesi gerekir.

## 2.3.2. nextflow.config dosyası aracılığıyla Docker kullanımını etkinleştirin

Bu yüzden Docker etkinken çalıştırmak için, kendimizi biraz ileride atacağız ve bir sonraki bölümde daha ayrıntılı olarak ele alacağımız Nextflow config dosyasını kullanacağız. Bu dizinde Nextflow Config adında bir dosyamız olduğunu görebilirsiniz ve burada zaten docker.enabled False var.

Docker'ı etkinleştirmek için bunu True olarak değiştireceğiz ve sonra workflow'u tekrar çalıştırmayı deneyebiliriz.

## 2.3.3. Docker etkinken workflow'u çalıştırın

Nextflow run hello containers nf ve bu sefer cowpy başarıyla çalıştı. Results'a bakalım. cowpy collected test ve işte Turkey'imiz. Harika.

Yani arka planda, Nextflow o işlem için mevcut bir konteyneri olduğunu biliyordu.

İmajı çekti ve komutları bizim için çalıştırdı.

## 2.3.4. Nextflow'un konteynerli görevi nasıl başlattığını inceleyin

Merak ediyorsanız, aslında tam olarak ne yaptığını work dizinine bakarak görebiliriz. Code work yaparsam, sonra hash ve sonra command run, hatırlarsanız o görev için yürütülen gerçek dosya, içeri girebilir ve NXF launch adlı bir fonksiyon arayabiliriz. Ve burada Nextflow'un kullandığı tam docker komutunu görebilirsiniz, bu daha önce terminalde manuel olarak yaptığımıza çok benziyor. Docker run. Bu ana dizini konteynere bağlama ve sonra konteyner URL'sini belirtme.

Yani burada sihir yok. Sadece Nextflow sizin için ağır işi otomatik olarak yapıyor, pipeline'ınızda konteynerleri kolayca belirtebileceğiniz ve sonra workflow'unuzu çalıştıran diğer herkesin kolayca kullanabileceği bir şekilde. Ve bu insanların artık analiz pipeline'ınızı çalıştırmak için yazılım yönetimi hakkında düşünmelerine gerek kalmıyor.

Çok, çok basit, çok kullanışlı ve aynı zamanda gerçekten tekrarlanabilir. Her açıdan iyi.

Tamam, harika iş. Bu Beşinci Bölümün sonu. Nextflow yapılandırmasından daha ayrıntılı olarak bahsedeceğimiz bu Merhaba Nextflow eğitiminin son bölümü olan altıncı bölüm için bir sonraki videoda bize katılın.

Bir sonraki videoda görüşürüz.

[Sonraki video transkripti :octicons-arrow-right-24:](06_hello_config.md)
