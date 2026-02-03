# Bölüm 6: Merhaba Config - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../06_hello_config.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba, Merhaba Nextflow eğitim kursunun altıncı bölümüne hoş geldiniz.

Bu bölümün adı Merhaba Config ve eğitim kursumuzun son bölümü.

Bu bölümde Nextflow konfigürasyonundan bahsedeceğiz. Nextflow konfigürasyonu gerçekten güçlüdür. Aynı pipeline'ı farklı yazılım sağlama yöntemleriyle, farklı hesaplama altyapılarında ve pipeline'ın kendisinde farklı seçeneklerle çalıştırmamızı sağlar.

Bu, başkaları tarafından oluşturulmuş Nextflow pipeline'larını, tamamen farklı bir altyapı için oluşturulmuş olsalar bile, kendi sisteminizde çalıştırabileceğiniz anlamına gelir. Nextflow'u yapılandırma yeteneği, iş akışlarını gerçekten taşınabilir ve paylaşılabilir hale getirir.

Bu bölümde, önceki bölümlerde oluşturduğumuz iş akışını kullanacağız, ancak iş akışı kodunu hiç düzenlemeyeceğiz. Sadece Nextflow yapılandırma dosyamıza bakacağız ve yapılandırmayı değiştirmenin Nextflow'un çalışma şeklini nasıl değiştirdiğini göreceğiz.

Tamam, hadi başlayalım.

Daha önce olduğu gibi, training.nextflow.io adresine giderek başlayalım. Sol tarafta Merhaba Nextflow ve altıncı bölüm. Merhaba config. Şimdi GitHub Codespaces ortamıma gidip kullanacağımız betiği kontrol edeceğim.

## 0. Isınma: Docker'ın etkin olduğunu kontrol edin ve Hello Config iş akışını çalıştırın

Bu betik Hello Config olarak adlandırılıyor ve daha önce olduğumuz yerden başlıyor. Yani üç parametremizle tamamen aynı görünüyor. CSV dosyası için greetings, çıktı toplu iş adı için batch ve cowpy adı için character. Farklı süreçlerin dört içe aktarımı var ve sonra bunları birbirine zincirleyen bir iş akışımız var.

Aslında bu dosyayı şimdi kapatacağım çünkü bu bölümde Nextflow dosyasına hiç dokunmayacağız. Tamamen yapılandırma dosyası içinde çalışacağız. Önceki beşinci bölümde kısaca baktığımız nextflow.config dosyasına bakarsam, burada tek bir ifade görebiliriz: docker.enabled = true, bu da Nextflow'a bu iş akışını çalıştırırken Docker kullanmasını söylüyor.

Burada, pipeline kökünde, Nextflow'u çalıştırdığımda otomatik olarak yüklenen nextflow.config kullanıyorum. Ancak unutmayın, Nextflow yapılandırma dosyalarını birden fazla yerden yükleyebilir.

Nextflow dokümanlarını kontrol edersem, Configuration'a gidersem, bu yerlerin bir listesini ve yüklendikleri önceliği görebilirsiniz.

Tamam. İş akışımızın beklediğimiz gibi çalıştığını kontrol edelim. Bir terminal açın. nextflow run hello-config yazın ve enter'a basın. Bu dört sürecin çalışması ve bir cowpy komutuyla bitmesi gerekiyor. Evet, bu düzgün çalıştı. Docker etkin olduğundan, Docker'ı çekti ve beşinci bölümün sonunda olduğu gibi benim için cowpy'yi çalıştırdı.

## 1. Hangi yazılım paketleme teknolojisinin kullanılacağını belirleyin

Tamam. Diyelim ki bir HPC üzerinde çalışıyorum ve Docker yüklü değil. Bu senaryoda yapılacak en iyi şey Singularity veya Apptainer kullanmak olurdu. Bunu yapacak olsaydım, cowpy modülüne gidip bu konteyneri Seqera Containers'dan da alabileceğiniz bir oras:// ile singularity imajı kullanacak şekilde değiştirirdim, önceki bölümde gösterdiğim gibi.

Sonra nextflow.config'e gidip docker.enabled'ı false olarak ayarlayıp singularity.enabled = true yapardım. Ya da Apptainer kullanıyorsam, apptainer.enabled = true yapardım ve bu işe yarardı.

Nextflow, konteynerler dışında başka teknolojileri de destekliyor, aşina olabileceğiniz bir şey conda. Burada conda.enabled = true yapıp Docker'ı false olarak ayarlayabiliriz. conda aynı container yönergesini kullanmıyor. Bunun yerine, conda adında yeni bir tane ekleyebiliriz. Sonra kullanmak istediğimiz conda paketini belirtiriz. Pipeline'ı mümkün olduğunca tekrarlanabilir kılmak için mümkün olduğunca spesifik olmak iyi bir pratiktir. Bu yüzden conda kanalını, conda-forge'u ve ardından cowpy'yi ve tam sürümü olan 1.1.5'i belirteceğim.

İsterse sadece cowpy da yazabilirdim, ancak bu pipeline'ın farklı çalıştırmalarında cowpy'nin farklı bir sürümüne çözümlenebilir.

Bunun güzel yanı docker yönergesine hiç dokunmamış olmam. Bu Docker imajı hala orada. Şimdi sadece iki alternatif sağlıyorum ve bunlar yalnızca bir yapılandırma dosyası kullanılarak açılıp kapatılabilir.

## 1.3. Conda kullanabileceğini doğrulamak için iş akışını çalıştırın

Conda artık etkin, o halde deneyelim.

Harika. Çalışıyor ve Nextflow'dan burada bir mesaj görebilirsiniz, Nextflow'un benim için bir conda ortamı oluşturduğunu söylüyor ve bu önbellek konumunu kullanıyor.

Arka planda, Nextflow benim için "conda create" komutları çalıştırarak sadece istediğim paketleri içeren yeni bir izole conda ortamı oluşturuyor ve ardından süreci çalıştırabilmesi için bu conda paketlerini yüklüyor ve indiriyor.

Orada biraz zaman aldığını görebilirsiniz çünkü ilk kez ortamı oluşturuyor ve yazılımı yüklüyordu. Ancak, bu ortamı önbelleğe aldı, bu yüzden aynı Nextflow komutunu tekrar çalıştırırsam, aynı conda ortamını yeniden kullanacağı için çok daha hızlı olmalı.

Bunun harika yanlarından biri, bu yönergelerin sadece tüm iş akışı için değil, süreç seviyesinde de belirtilebilmesidir. Yani istersen, farklı süreçler için hangi teknolojinin kullanıldığını karıştırıp eşleştirebilirsin.

## 2. Süreç yönergeleriyle hesaplama kaynaklarını tahsis edin

Nextflow yapılandırma dosyası sadece yazılım paketlemeden çok daha fazlasını yapabilir. Ayrıca Nextflow'a pipeline'daki adımları gerçekte nasıl çalıştıracağını da söyleyebiliriz. Bir örnek, bir ana sistem sistemine her çalıştırılan göreve hangi kaynakların kullanılabilir hale getirilmesi gerektiğini söylemektir.

Varsayılan olarak, Nextflow çok fazla şey vermez. Her sürece tek bir CPU ve sadece iki gigabayt bellek verir.

Bu muhtemelen değiştirmek isteyeceğimiz bir şeydir, böylece çalışması uzun süren süreçler daha fazla kaynağa sahip olabilir ve daha hızlı çalışabilir, ancak bir sürece ne tahsis edeceğinizi bilmek zor olabilir. Nextflow bunun için size yardımcı olacak bazı güzel hilelere sahip.

## 2.1. Bir kaynak kullanım raporu oluşturmak için iş akışını çalıştırın

İş akışını tekrar çalıştıralım. Bu sefer ek bir argüman ekleyeceğim, -with-report. Bu bir temel Nextflow seçeneği, bu yüzden tek bir tire. Ve sonra istediğim dosya adı. Bu durumda, buna report-config-one.html diyeceğim.

İş akışını tekrar çalıştıracağım. Daha önce olduğu gibi çalışacak, ancak bana ek bir yardımcı rapor verecek, görebilirsiniz şimdi kenar çubuğunda ortaya çıktı.

Bu dosyaya sağ tıklayacağım, download'a tıklayacağım, bu dosyayı GitHub Codespaces'den yerel sistemime indirir, böylece web tarayıcısında kolayca görüntüleyebilirim.

Bu rapor herhangi bir Nextflow çalıştırması için oluşturulabilir ve içinde çok fazla bilgi var. En üstte, hangi komutun kullanıldığı, iş akışının ne zaman çalıştığı, ne kadar sürdüğü hakkında bazı metaverilerle başlar, ancak aşağı kaydırdıkça, pipeline'daki her adım tarafından kullanılan kaynaklar hakkında daha ayrıntılı bilgi alırız.

Her süreç farklı görevler için birden çok kez çalıştığından. Kullandığımız kaynakların her süreç için varyasyonunu gösteren bir kutu grafiğimiz var.

Biraz daha aşağı kaydırırsam, kullanılan bellek ve görev süresi hakkında benzer bilgiler görürüm. Ayrıca disk okuma yazma.

Uzun süren görevlere sahip büyük bir pipeline için, bu, talep ettiğiniz kaynakların konfigürasyonunu ince ayarlamak için, fazla talep etmemeniz ama aynı zamanda hızlı çalışması için yeterli sağlayabilmeniz için çok bilgilendirici olabilir.

Raporu aşağı kaydırmaya devam edersem, ayrıca bir görev tablosu görürüz, bu bize iş akışında çalıştırılan her bir görev hakkında ayrıntılı bilgi gösterir. Bu, çalıştırılan çözümlenmiş betik gibi bilgileri içerir.

Tamam, yapılandırma dosyamıza geri dönelim. İş akışımız için gerçekten fazla bir şeye ihtiyacımız olmadığını gördük, o halde Nextflow'a iş akışındaki her süreç için sadece bir gigabayt belleğe ihtiyacımız olduğunu söyleyelim.

Şimdi bunu bu şekilde process seviyesinde tanımladığımızda, bu pipeline'daki her bir sürece uygulanır.

## 2.3. Tek bir süreç için kaynak tahsisleri ayarlayın

Tartışma uğruna, cowpy'nin gerçekten çok fazla iş yaptığını ve diğer görevlerden daha fazla kaynağa ihtiyaç duyduğunu varsayalım. Burada sadece o sürece uygulanan ekstra bir yapılandırma bloğu tanımlayabiliriz, withName: cowpy kullanarak.

Buna yapılandırma seçici denir ve burada farklı süreçlerle eşleşmek için farklı desenler tanımlayabiliriz. Örneğin, cow\* yapabilirim. Sonra bunu bazı süslü parantezlerle takip ediyorum ve ona bir yerine iki gigabayt bellek verelim ve diyelim ki iki CPU.

Şimdi Nextflow iş akışındaki her sürece bir gigabayt verecek, bu istek hariç, daha spesifik. Yani onu geçersiz kılar. Ve sadece cowpy olarak adlandırılan süreçler için iki gig bellek ve iki CPU alacak.

Nextflow'un kaynak kullanımı konusunda akıllı olduğunu unutmayın. Bu nedenle, bu sayıları daha yüksek değerlere getirmeye başlarsanız, Nextflow'un tüm görevleri paralel olarak çalıştırmak yerine, görev gönderimlerini birbiri ardına kuyruğa almaya başladığını göreceksiniz, böylece mevcut kaynakları fazla talep etmez.

## 2.4. Değiştirilmiş yapılandırmayla iş akışını çalıştırın

İş akışını tekrar çalıştırmayı deneyelim ve bu sefer yeni bir rapor kaydedelim.

Tamam, bu dosyayı indirebilir ve bir göz atabiliriz.

Evet, şaşırtıcı olmayan bir şekilde, temelde tamamen aynı görünüyor çünkü bu gerçek bir şey yapmayan sahte bir iş akışı. Ancak sınırları tanımlamanın ve bu tür raporlamayla gerçek hayattan iş akışları yapmanın bu yinelemeli yaklaşımının, uygun yapılandırma ayarlamak ve sahip olduğunuz hesaplama kaynaklarından gerçekten en iyi şekilde yararlanmak için kanıta dayalı bir yaklaşım yapmanızı nasıl sağladığını hayal edebilirsiniz.

Bu konuda gerçekten akıllı olmaya başlayabilirsiniz. Nextflow'un hataları yeniden deneme yerleşik bir yeteneği vardır ve yapılandırma dosyanızda bir closure kullanarak ve kullanılabilir hale getirilen kaynakları dinamik olarak ayarlayarak bundan yararlanabilirsiniz. Yani burada Nextflow'a bu iki gigabaytı yeniden deneme denemesiyle çarmasını söyledim. Yani ikinci yeniden deneme dört gig alacak, üçüncü yeniden deneme altı gig alacak ve böyle devam edecek. Bu biraz bu eğitim kursunun kapsamının ötesinde, ancak ilgileniyorsanız, dinamik yeniden deneme mantığı hakkında güzel bir bölümü olan Nextflow dokümanlarına göz atın.

## 2.5. Kaynak sınırları ekleyin

Şimdi, bu konuda fark edebileceğiniz bir şey, bu tür şeylerin yanlışlıkla sisteminizdeki mevcut kaynaklardan fazlasını istemek oldukça kolay hale getirebileceğidir. Mevcut olandan daha fazla kaynak talep ederseniz Nextflow yapılandırmanız hakkında bir hata verecek ve çalıştırmayı durduracaktır. Bunu önlemek için kaynak sınırları adı verilen bir şey kullanabilirsiniz.

İş akışımızda process kapsamı altında, bu şekilde bir dizi alan kaynak sınırları tanımlayabiliriz ve bu sistemde mevcut olan maksimum bellek, CPU ve zamanı belirtebiliriz.

Burada yüksek değerler ayarlamak, istenen kaynakların miktarını artırmaz. İsteklerimizde hala bir gigabayt kullanıyor olacağız, ancak bu isteklerden herhangi biri 750'ye ulaşırsa, o tavana çarpacakları ve bundan fazlası istenmeyeceği anlamına gelir, bu da Nextflow'un çalışmaya devam edeceği ve kullanılamayan kaynaklar nedeniyle çökmeyeceği anlamına gelir.

Yani bu kullanılacak güzel bir güvenlik önlemidir, özellikle kaynak tahsisinizle dinamik mantık kullanıyorsanız.

Bunun gerçekten yararlı olduğu diğer durum, herkese açık ve sizin tarafınızdan kontrol edilmeyen pipeline'lar kullanıyorsanızdır. Yapılandırma varsayılanlarıyla gelebilirler ve Nextflow otomatik olarak herhangi bir kaynak talebini sisteminizde çalışacak şekilde eşiklemenin doğru yaklaşımını alacaktır.

Tamam, harika. Yazılımdan bahsettik. Kaynak tahsisinden bahsettik ve tüm süreçler ve belirli süreçler için hem yapılandırmanın farklı kapsamlarını açıkladık.

## 3. İş akışı parametrelerini saklamak için bir parametre dosyası kullanın

Tamam, şimdi dikkatimizi parametrelere çevireceğiz. Yapılandırma dosyasında daha önce Nextflow betiğinde yaptığımız gibi parametreler tanımlayabiliriz. Yani params.greeting = 'hello' veya params kapsamını kullanıp foo = 'bar' ayarlayabiliriz.

Ve bu, iş akışınız için varsayılanları ayarlamak için harika. Ancak, pipeline'ları çalıştırdığınızda, parametreleri bir JSON veya bir YAML dosyasında belirtmek güzel olabilir.

Böyle bir dosya kullanmak, tire tire ile komut satırı seçenekleri belirtmekten çok daha iyidir. Bir iş akışını çalıştırdığınızda, birçok parametre belirtmeniz gerekebilir ve hepsini tek bir CLI'de yazmak sıkıcı ve hataya açık olabilir. Ayrıca, kullandığınız tüm parametreleri hatırlamanız olası değildir, bu yüzden bunu bir dosyaya kodlarsanız, gelecekte aynı parametreleri kullanarak iş akışını tekrar başlatmak daha kolaydır.

Burada test-params adında örnek bir dosyamız var ve bunun iş akışımızda sahip olduğumuz üç parametreyi üç farklı değerle belirttiğini görebilirsiniz. Şahsen, YAML'ı JSON'dan yazmayı daha kolay buluyorum. O yüzden sadece işe yaradığını göstermek için, test.yaml adında yeni bir dosya oluşturacağım ve bunları kopyalayacağım, tırnak işaretlerinden kurtulacağım. Ve kaydet'e basacağım.

Bu JSON ve YAML dosyalarının yazılması daha kolay olabilir çünkü daha tanıdık bir sözdizimi vardır. Ancak bunların yalnızca parametreler için olduğunu ve yalnızca bunun gibi anahtar değer sözdizimi aldıklarını unutmayın.

## 3.1. Bir parametre dosyası kullanarak iş akışını çalıştırın

Hadi deneyelim. Daha önce olduğu gibi aynı komutu yapın. Raporu kaldırın ve -params-file test-params.yaml yapacağım.

Hayır, bu bir temel Nextflow seçeneği, yani tek bir tire.

Tamam. İş akışını çalıştırdı ve hepsini komut satırında belirtmem yerine o YAML dosyasındaki parametreleri kullandı. Sadece bu basit örnek için aşırıya kaçmış gibi görünebilir, ancak 10 veya 20 farklı parametreniz varsa, manuel olarak yazmak acı verici olabilir ve bu, bir kod düzenleyicide düzenlemek ve tekrarlanabilirlik adına saklamak için çok daha kolaydır.

## 3. İşi yapmak için hangi yürütücü(ler)in kullanılması gerektiğini belirleyin

Tamam. Docker ve conda ile yazılım paketlemeden bahsettik. CPU ve bellekle süreç kaynak gereksinimlerinden bahsettik. Ve iş akışlarını çalıştırırken parametrelerin nasıl belirtileceği hakkında biraz konuştuk.

Yapılandırmanın son kısımları gerçekten yürütme, temel hesaplama altyapısının kendisidir ve bu Nextflow'un taçındaki gerçek mücevherdir: bu aynı iş akışını birden fazla farklı hesaplama altyapısında çalıştırabiliriz.

Aslında bir saniyeliğine yazılı eğitim materyaline geçeceğim. Bu eğitim bölümünün altında, farklı yürütücülerin, bu durumda HPC zamanlayıcıların, bir iş göndermek için gereken kaynak gereksinimlerini nasıl tanımladıklarının birkaç farklı örneğini görebiliriz.

Yani Slurm için, --mem ve CPU numarasını tanımlayan bu SBATCH başlıkları var. PBS kullanıyorsanız, farklı başlıklarınız var ve Grid Engine kullanıyorsanız, yine farklı başlıklar var.

AWS Batch, Google Cloud, Azure veya daha fazlası olsun, bulutta çalıştırmak istiyorsanız bunun daha da farklı olduğunu hayal edebilirsiniz.

Bu temel hesaplama altyapılarının her birine yürütücü denir ve Nextflow bu farklı yürütücülerin hepsine doğru sözdizimi ile iş göndermek için nasıl konuşacağını bilir.

İyi haber şu ki bunun hakkında bilmenize gerek yok. Tek yapmanız gereken Nextflow'a hangi yürütücüyü kullanacağını söylemek.

## 3.1. Farklı bir arka ucu hedefleme

Yapılandırma dosyamıza geri dönüyoruz ve process altında executor yapıyoruz ve local yazacağım.

Local aslında varsayılandır, başka bir yürütücü belirtmezseniz, local kullanılacaktır ve bu sadece ana sisteminiz anlamına gelir, Nextflow'u nerede başlattıysanız.

Bunun yerine Slurm belirtebilirim. Ve bu Slurm işleri gönderir ya da AWS Batch diyebilirim ve bu AWS batch'e işler gönderir.

Bazı durumlarda ek yapılandırmaya ihtiyacınız var, örneğin bulutta çalışmak belirli kimlik bilgileri gerektirecektir, ancak gerçekten bunun özü budur ve iş akışınızı tamamen farklı bir hesaplama ortamında çalıştırmak için bir veya iki satır yapılandırma kadar basit olabilir.

Codespaces içinde basit bir sistemde çalışıyor olsak bile, yine de bununla biraz oynayabilir ve Slurm üzerinde çalışıyormuşuz gibi davranabilirim. Eğer sonra iş akışını tekrar başlatırsam, nextflow run hello-config. Slurm'a iş gönderemeyeceği için başarısız olacak. Ancak yine de çalışma dizinlerine gidip Nextflow'un ne yaptığını görebiliriz. Yani bu çalışma dizinine gidersek ve .command.run'a bakarsak. Bu dosyanın en üstünde, Slurm işi için gereken kaynakları belirtmeye çalışan bu sbatch başlık satırlarına sahibiz.

## 4. Önceden ayarlanmış yapılandırmaları seçmek için profiller kullanın

Tamam, neredeyse geldik. Bu bölümün son kısmı yapılandırma profilleri hakkında konuşmak. Pipeline'ınızı birkaç farklı sistemde çalıştırıyorsanız, her seferinde belirtmeniz gereken tüm bu farklı Nextflow yapılandırma dosyalarına sahip olmak can sıkıcı olabilir.

Bunun yerine, Nextflow yapılandırma dosyanızda yapılandırma gruplamalarını kodlayabilir ve bu grupları bir profil bayrağı kullanarak açıp kapatabilirsiniz. Bunun nasıl göründüğüne bakalım.

## 4.1. Yerel geliştirme ve HPC'de yürütme arasında geçiş yapmak için profiller oluşturun

Örneğimizde iki profil oluşturacağız, biri dizüstü bilgisayarım için diğeri daha ağır bir HPC sistemi için. Biraz hile yapacağım ve kodu eğitim materyalinden kopyalayıp buraya yapıştıracağım.

profiles adında yeni bir kapsamımız var ve sonra her profil için, herhangi bir şey olabilecek bir adımız var. Ve bunun içinde, daha önce yazdığımız üst düzey yapılandırmayla tamamen aynı görünen bir yapılandırmamız var. Yani yine, process kapsamımız var. Docker kapsamı.

my_laptop adlı profilde. Local yürütücüsü kullanarak çalıştırmasını, yani ana sistemimde çalışmasını ve Docker kullanmasını söylüyorum.

university_hpc profilinde burada, iş göndermek için Slurm kullanmasını, Docker yerine conda kullanmasını ve kullanıyorum HPC'deki düğümlerin sistem boyutuyla eşleşebilecek farklı kaynak sınırları belirtiyorum.

Varsayılan olarak, Nextflow'u çalıştırdığımda bu yapılandırmanın hiçbiri kullanılmayacak, bu profillerden birini kullanmak istediğimi belirtmem gerekiyor.

## 4.2. Bir profille iş akışını çalıştırın

Hadi nextflow run hello-config yapalım. Ve -profile yapacağım, tek tire çünkü bu bir temel Nextflow seçeneği. Ve sonra verdiğim ad, my_laptop. Nextflow şimdi bu yapılandırma profili içinde belirtilen yapılandırma bloğunu kullanmalı ve Nextflow'u çalıştırdığında onu uygulamalı. Diğer yapılandırma bloğunu kullanmak isterse, sadece o profil adını değiştirmem gerekiyor. Hatırlaması çok daha kolay. Kullanması çok daha kolay.

## 4.3. Bir test profili oluşturun

Not, profiller her türlü yapılandırmaya sahip olabilir, bu yüzden yürütme ortamınızla ilgili olması gerekmez. Örneğin, burada bir dizi parametreye sahip yeni bir profil oluşturalım. Bunu tux olarak değiştirebilir ve my_profile olarak değiştirebiliriz ve şimdi profile test yaptığımızda, iş akışının en üst düzeyinde belirtilen parametrelerin üzerine yazılacak bu parametreleri belirtecektir.

Nextflow'u çalıştırdığınızda, birden fazla profili zincirleyebilirsiniz ve sırayla uygulanacaklardır.

## 4.4. Test profiliyle iş akışını yerel olarak çalıştırın

Yani önceki komutu alıp virgül test yapabilirim. Bu önce my_laptop yapılandırmasını uygulayacak ve sonra test yapılandırmasını uygulayacaktır. Herhangi bir örtüşme varsa, o zaman sağdaki profil önceki profillerdeki herhangi bir yapılandırmanın üzerine yazacaktır. Enter'a basarsam, bakalım ne olacak.

Tamam, burada yeni bir sonuç dosyamız var. Seçeneklerden biri olarak belirttiğim My Profile'ı görebilirsiniz. Ve ayrıca cowpy, my_profile'ı görebiliriz ve elbette orada tux var. Yani bu işe yaradı.

## Özet

Tamam! İnanılmaz. Bu kadar. Kursun sonuna geldiniz. Biraz kutlama konfetisi alıyorsunuz. Bu bölümü bitirdiğiniz için tebrikler.

[Sonraki video transkripti :octicons-arrow-right-24:](07_next_steps.md)
