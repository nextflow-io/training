# Bölüm 4: Hello Modules - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım tüm talimatlar için [ders materyaline](../04_hello_modules.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgi amaçlıdır ve materyallerdeki tüm bölüm numaralarını kapsamayabilir.

## Hoş Geldiniz

Merhaba ve Hello Nextflow'un dördüncü bölümüne tekrar hoş geldiniz. Bu bölüm tamamen modüllerle ilgili ve kursun oldukça kısa bir bölümü. Aslında çok fazla kod yazmayacağız; bu bölüm daha çok boru hattımızdaki kodu nasıl organize ettiğimizle ilgili.

Şimdiye kadar her şeyi tek bir dosyaya koyuyorduk; bu gayet iyi bir yaklaşım ve aslında eskiden Nextflow boru hatlarını bu şekilde oluştururduk.

Ancak boru hattı büyüdükçe betik giderek uzuyor, gezinmesi ve bakımı giderek zorlaşıyor; üstelik kodu paylaşmak da pek mümkün olmuyor.

Nextflow modülleri, süreçleri ana betikten çıkarmamıza ve ardından içe aktarmamıza olanak tanır. Bu sayede kod üzerinde gezinmek kolaylaşır ve modül kodunu farklı boru hatları arasında paylaşabilir hale geliriz.

Belgelerin ana sayfasındaki bu küçük diyagram kavramı güzel bir şekilde gösteriyor. Tek bir büyük betik yerine, farklı modül betiklerinden bu ayrı modül dosyalarını dahil edeceğiz; her şey iş akışına çekilecek ama tamamen aynı şekilde çalışmaya devam edecek.

Hadi GitHub Codespaces'e geçelim ve biraz etrafı inceleyelim. Daha önce olduğu gibi, çalışma alanımı biraz temizledim. Eski Nextflow dizinlerini, work dizinini ve benzerlerini kaldırdım. Ancak bu dosyalar hâlâ duruyorsa sorun değil.

Hello Modules dosyasında çalışmaya başlayacağım; bu dosya temelde önceki bölümün sonunda bıraktığımız hâliyle aynı. Burada üç sürecimiz var. Birkaç params, iş akışı bloğu var; bu blokta üç süreci çalıştırıyor ve kanallarla birbirine bağlıyoruz. Ardından çıktı kanallarını yayımlıyoruz ve bu dosyaların nasıl yayımlanacağını belirten çıktı bloğumuz var.

## 1. Modülleri Depolamak İçin Bir Dizin Oluşturun

Şimdi, dediğim gibi, çok fazla kod yazmayacak ya da düzenlemeyeceğiz. Sadece mevcut kodu taşıyacağız. Nextflow modül dosyaları genellikle tek bir süreç içerir ve kurallara göre bunları normalde `modules` adlı bir dizinde tutarız. Ancak bu dizine istediğiniz adı verebilirsiniz. Ben buradaki depomda bir `modules` dizini oluşturacağım ve ardından her süreç için ayrı bir dosya oluşturacağım. Yani yeni dosya oluşturuyorum: `sayHello.nf`.

## 2. sayHello() için Bir Modül Oluşturun

Şimdi sürecimi alacağım ve bu kodu seçip ana Hello Modules dosyasından kesip buraya yapıştıracağım.

Tabii ki bu tek başına bir şey yapmaz. Ana betiğimiz hâlâ o sürece ihtiyaç duyuyor, dolayısıyla onu bir şekilde geri çekmemiz gerekiyor. Bunu `include` ifadesiyle yapıyoruz.

`include` yazıyorum, ardından süslü parantezler açıyorum ve sürecin adını yazıyorum. Sonra `from` diyorum ve göreli bir dosya yolu veriyorum. Bu betik nerede kayıtlıysa oradan göreli olduğu için `./` ile başlıyor. Yani `modules/sayHello.nf`.

VS Code uzantısının burada oldukça yardımcı olduğuna dikkat edin. Bu dosyayı bulup bulamadığını ve içe aktarmaya çalıştığım süreci bulup bulamadığını bize söylüyor. Kasıtlı olarak bir yazım hatası yaparsam hemen bir hata veriyor ve içe aktarmaya çalıştığım süreci bulamadığını söylüyor. Bu nedenle karşılaştığınız hatalara dikkat edin.

İşte bu kadar. Sürecimiz hâlâ burada. Aşağıda herhangi bir değişiklik yapmaya gerek yok. Süreç aynı ada sahip ve tamamen aynı şekilde çalıştırılıyor. Tek fark, sürecin gerçek kodu artık ayrı bir dosyada.

Nextflow iş akışını yeniden çalıştırabiliriz; tamamen aynı şekilde çalışacak. Kursun bu bölümünün geri kalanı da temelde bu üç süreci kendi dosyalarına taşımaktan ibaret.

Hadi şimdi bunu yapalım. İkinci süreç için hızlıca yeni bir modül dosyası oluşturacağım: `convertToUpper.nf`. O kodu kesip oraya yapıştıracağım. Ardından onu da dahil edeceğim.

Sonra `collectGreetings.nf` için yeni bir dosya oluşturacağım. Onu da kesiyorum.

Çok fazla kesip yapıştırma var.

Şimdi ana iş akışı betiğimiz aniden çok daha kısa, çok daha anlaşılır ve okunması çok daha kolay bir hâle geldi.

Projenin artık farklı dosyalarla nasıl şekillenmeye başladığını görebilirsiniz. İstediğimiz yerlerde ayrıntılara dalabilir, boru hattındaki belirli adımları çok daha kolay bulabilir ve boru hattının ne yaptığına dair hızlıca genel bir bakış elde edebiliriz.

## VS Code ile Modüllerde Gezinme

Tabii ki bunun dezavantajı şu: büyük bir boru hattınız varsa çok sayıda modül dosyanız olacak ve bunlar birden fazla alt dizinde ya da çeşitli yerlerde organize edilmiş olabilir. Yine de burada küçük bir ipucu vermek istiyorum. VS Code uzantısı, kod tabanınızda gezinme ve kod hakkında bilgi verme konusunda oldukça iyidir.

VS Code'un bu süreci anladığını ve üzerine geldiğimde küçük bir özet gösterdiğini görebilirsiniz; böylece kaynak koda gitmek zorunda kalmadan girdilerin ve çıktıların ne olduğunu görebiliyorum; bunlar genellikle bir iş akışında kullanırken en önemli bilgilerdir.

Ayrıca Mac'te Command tuşunu basılı tutup süreç adına tıklarsam dosyayı doğrudan açıyor. Gerçek dosya yollarını düşünmek zorunda kalmadan oraya atlayabiliyorum. Bu özellik her yerde çalışıyor; süreçlerin çağrıldığı yerlerde de kullanabiliyorum. Bu da işleri gerçekten hızlandırıyor.

## 4.4. İş Akışını Çalıştırın

Tamam, boru hattının hâlâ beklediğimiz gibi çalışıp çalışmadığını kontrol edelim. Terminali açalım. `nextflow run hello modules` yazalım ve herhangi bir sorun olmadan çalışıp çalışmadığını görelim.

Umarım tüm bunun amacı boru hattının temelde değişmemiş olması, dolayısıyla önceki çalıştırmadan farklı bir şey görmemeniz gerekiyor. Buradaki çıktı tamamen aynı görünüyor ve aynı dosyaların hepsinin bulunduğu `results` dizinimizi görebiliyoruz; bu harika. Değişiklik olmaması iyi bir şey.

## nf-core/modules Hakkında Bir Not

Konuyu kapatmadan önce, modüller söz konusu olduğunda iş birliğinin gücüne kısaca değinmek istiyorum. Bu dosyalar benim depomda duruyor, dolayısıyla üzerlerinde nasıl iş birliği yapabileceğimiz hemen anlaşılmıyor. Bunu yapmanın birçok farklı yolu var, ancak muhtemelen en büyük ve en iyi bilinen örnek nf-core'dur.

nf-core web sitesine gidip kaynaklar ve modüller bölümüne bakarsanız, nf-core'un devasa bir modül kütüphanesine sahip olduğunu görebilirsiniz; bunu görüntülediğimde neredeyse 1700'e yakın modül var. Favori araçlarımdan herhangi birinin adını yazabilir, başka birinin bunun için zaten bir modül yazıp yazmadığını bulabilir ve tüm girdileri, çıktıları, yazılım container'larını ve diğer bilgileri içeren bu önceden yazılmış modül sürecini görebilirim. Yan tarafta ise bu tek paylaşılan süreci kaç farklı nf-core boru hattının kullandığını görebilirsiniz.

Bu biraz aşırı bir örnek, ancak kodun gerçekten yeniden kullanıldığını görebiliyorsunuz. GitHub kaynağına tıklarsam, yaptığımız şeyle tamamen aynı. Sadece bir dosyadaki büyük bir süreç.

nf-core tarafında, bu dosyaları paylaşmak ve farklı depolara aktarmak için bazı teknikler kullanıyoruz. Bununla ilgili daha fazla bilgi edinmek istiyorsanız, özellikle nf-core kullanımı ve nf-core ile geliştirme hakkındaki kursuma göz atın. Ancak kod yeniden kullanımı kavramının ne kadar güçlü olabileceği konusunda size bir fikir vermek istedim.

## Özet

Modüller için bu kadar. Kursun kısa bir bölümü olduğunu söylemiştim. Sınava bakın, her şeyi anladığınızdan ve her şeyin hâlâ düzgün çalıştığından emin olun. Bir sonraki videoda görüşürüz; o video tamamen yazılım container'larıyla ilgili. Çok teşekkürler.
