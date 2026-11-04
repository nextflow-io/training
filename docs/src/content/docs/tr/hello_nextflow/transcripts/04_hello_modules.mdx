# Bölüm 4: Hello Modules - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım tam talimatlar için [ders materyaline](../04_hello_modules.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgi amaçlıdır ve materyallerdeki tüm bölüm numaralarını kapsamayabilir.

## Hoş Geldiniz

Merhaba ve Hello Nextflow'un dördüncü bölümüne tekrar hoş geldiniz. Bu bölüm tamamen modüllerle ilgili ve kursun oldukça kısa bir bölümü. Aslında çok fazla kod yazmayacağız; bu bölüm daha çok pipeline'ımızdaki kodu nasıl organize ettiğimizle ilgili.

Şimdiye kadar her şeyi tek bir dosyaya koyuyorduk. Bu gayet iyi bir yaklaşım ve aslında eskiden Nextflow pipeline'larını bu şekilde oluştururduk.

Ancak pipeline büyüdükçe betik giderek uzuyor, gezinmesi ve bakımı giderek zorlaşıyor. Üstelik bu yaklaşımda kodu paylaşmak da pek mümkün olmuyor.

Nextflow modülleri, süreçleri ana betikten çıkarıp içe aktarmamıza olanak tanır. Bu sayede kod gezinmesi kolaylaşır ve modül kodunu farklı pipeline'lar arasında paylaşabilir hale geliriz.

Belgelerin ana sayfasındaki bu küçük diyagram kavramı güzel bir şekilde gösteriyor. Tek bir büyük betik yerine, farklı modül betiklerinden bu ayrı modül dosyalarını dahil edeceğiz; her şey iş akışına çekilecek, ancak tam olarak aynı şekilde çalışmaya devam edecek.

Hadi GitHub Codespaces'e geçelim ve biraz inceleyelim. Daha önce olduğu gibi, çalışma alanımı biraz temizledim. Eski Nextflow dizinlerini, work dizinini ve benzerlerini kaldırdım. Ancak bu dosyalar hâlâ elinizde olsa da sorun değil.

Hello Modules dosyasında çalışmaya başlayacağım; bu dosya temelde önceki bölümün sonunda bıraktığımız haldedir. Burada üç sürecimiz var. Birkaç parametre, iş akışı bloğu var; bu blokta üç süreci çalıştırıyor ve kanallarla birbirine bağlıyoruz. Ardından çıktı kanallarını yayımlıyor ve bu dosyaların nasıl yayımlanacağını belirten çıktı bloğuna sahibiz.

## 1. Modülleri Depolamak İçin Bir Dizin Oluşturma

Şimdi, dediğim gibi, çok fazla kod yazmayacak ya da düzenlemeyeceğiz. Sadece mevcut kodu taşıyacağız. Nextflow modül dosyaları genellikle tek bir süreç içerir ve kurala göre bunları normalde `modules` adlı bir dizinde tutarız. Ancak bu dizine istediğiniz adı verebilirsiniz. Ben buradaki depomda bir `modules` dizini oluşturacağım ve ardından her süreç için ayrı bir dosya oluşturacağım. Yeni dosya oluşturuyorum: `sayHello.nf`.

## 2. sayHello() için Bir Modül Oluşturma

Şimdi sürecimi alacağım ve bu kodu seçip ana hello modules dosyasından kesip buraya yapıştıracağım.

Tabii ki bu tek başına bir şey yapmaz. Ana betiğimiz hâlâ o sürece ihtiyaç duyuyor, dolayısıyla onu bir şekilde geri çekmemiz gerekiyor. Bunu `include` ifadesiyle yapıyoruz.

`include` yazıyorum, ardından süslü parantezler açıyorum ve sürecin adını yazıyorum. Sonra `from` diyorum ve göreli bir dosya yolu veriyorum. Bu betik, kaydedildiği konuma göre göreli olduğundan `./` ile başlıyor. Yani `modules/sayHello.nf` oluyor.

VS Code uzantısının burada oldukça yardımcı olduğuna dikkat edin. Bu dosyayı bulup bulamadığını ve içe aktarmaya çalıştığım süreci bulup bulamadığını bize söylüyor. Kasıtlı olarak bir yazım hatası yaparsam hemen bir hata veriyor ve içe aktarmaya çalıştığım süreci bulamadığını söylüyor. Bu nedenle karşılaştığınız hatalara dikkat edin.

İşte bu kadar. Sürecimiz hâlâ burada. Aşağıda herhangi bir değişiklik yapmaya gerek yok. Süreç aynı ada sahip ve tam olarak aynı şekilde çalıştırılıyor. Tek fark, sürecin gerçek kodunun artık ayrı bir dosyada olması.

Nextflow iş akışını yeniden çalıştırabiliriz; tam olarak aynı şekilde çalışacak. Kursun bu bölümünün geri kalanı da temelde bu üç süreci kendi dosyalarına taşımaktan ibaret.

Hadi şimdi bunu yapalım. `convertToUpper.nf` adında ikinci süreç için hızlıca yeni bir modül dosyası oluşturacağım. O kodu kesip oraya yapıştıracağım. Ardından onu da dahil edeceğim.

Sonra `collectGreetings.nf` için yeni bir dosya oluşturacağım. Onu da kesiyorum.

Çok fazla kesme, kopyalama ve yapıştırma var.

Şimdi ana iş akışı betiğimiz aniden çok daha kısa, çok daha anlaşılır ve okunması çok daha kolay bir hale geldi.

Projenin artık farklı dosyalarla nasıl şekillendiğini görebilirsiniz. İstediğimiz yerlerde ayrıntılara dalabilir, pipeline'daki belirli adımları çok daha kolay bulabilir ve pipeline'ın ne yaptığına dair hızlıca genel bir bakış elde edebiliriz.

## VS Code ile Modüllerde Gezinme

Tabii ki bunu yapmanın dezavantajı şu: Büyük bir pipeline'ınız varsa çok sayıda modül dosyanız olacak ve bunlar birden fazla alt dizinde ya da çeşitli yerlerde organize edilmiş olabilir. Yine burada küçük bir ipucu vermek istiyorum. VS Code uzantısı, kod tabanınızda gezinme ve kod hakkında bilgi verme konusunda oldukça iyidir.

VS Code'un bu süreci anladığını ve üzerine geldiğimde küçük bir özet gösterdiğini görebilirsiniz. Böylece kaynak koda gitmek zorunda kalmadan girdilerin ve çıktıların ne olduğunu görebiliyorum; bunlar bir iş akışında kullanırken genellikle en önemli bilgilerdir.

Ayrıca Mac'te Command tuşunu basılı tutup süreç adına tıklarsam dosyayı doğrudan açıyor. Gerçek dosya yollarını düşünmek zorunda kalmadan oraya atlayabiliyorum. Bu özellik her yerde çalışıyor; süreçlerin çağrıldığı yerlerde de kullanabiliyorum. Bu da gezinmeyi gerçekten hızlandırıyor.

## 4.4. İş Akışını Çalıştırma

Tamam, pipeline'ın hâlâ beklediğimiz gibi çalışıp çalışmadığını kontrol edelim. Terminali açalım. `nextflow run hello modules` yazalım ve herhangi bir sorun olmadan çalışıp çalışmadığını görelim.

Umarım tüm bunun amacı pipeline'ın temelde değişmemiş olması, dolayısıyla önceki çalıştırmadan farklı bir şey görmemeniz gerekiyor. Buradaki çıktı tamamen aynı görünüyor ve aynı dosyaların hepsinin bulunduğu sonuçlar dizinimizi görebiliyoruz. Harika. Değişiklik olmaması iyi bir şey.

## nf-core/modules Hakkında Bir Not

Konuyu kapatmadan önce, modüller söz konusu olduğunda iş birliğinin gücüne kısaca değinmek istiyorum. Bu dosyalar benim depomda duruyor, dolayısıyla üzerlerinde nasıl iş birliği yapabileceğimiz hemen anlaşılmıyor. Bunu yapmanın birçok farklı yolu var, ancak muhtemelen en büyük ve en iyi bilinen örnek nf-core'dur.

nf-core web sitesine gidip kaynaklar ve modüller bölümüne bakarsanız, nf-core'un devasa bir modül kütüphanesine sahip olduğunu görebilirsiniz; bunu görüntülediğimde neredeyse 1700'e yakın modül var. Favori araçlarımdan herhangi birinin adını yazabilir, başka birinin bunun için zaten bir modül yazıp yazmadığını bulabilir ve tüm girdileri, çıktıları, yazılım konteynerlerini ve diğer bilgileri içeren bu önceden yazılmış modül sürecini görebilirim. Yan tarafta ise bu tek paylaşılan süreci kullanan kaç farklı nf-core pipeline'ının olduğunu görebilirsiniz.

Bu biraz aşırı bir örnek, ancak kodun gerçekten yeniden kullanıldığını görebiliyorsunuz. GitHub kaynağına tıklarsam, yaptığımız şeyle tamamen aynı. Sadece bir dosyadaki büyük bir süreç.

nf-core tarafında, bu dosyaları paylaşabilmek ve farklı depolara dahil edebilmek için bazı teknikler kullanıyoruz. Bununla ilgili daha fazla bilgi edinmek istiyorsanız, özellikle nf-core kullanımı ve nf-core ile geliştirme hakkındaki kursuma göz atın. Ancak kod yeniden kullanımı kavramının ne kadar güçlü olabileceği konusunda size bir fikir vermek istedim.

## Kapanış

Modüller için bu kadar. Size kısa bir bölüm olduğunu söylemiştim. Sınava göz atın, konuyu anladığınızdan emin olun ve her şeyin hâlâ düzgün çalıştığını kontrol edin. Bir sonraki videoda görüşürüz; o bölüm tamamen yazılım konteynerleriyle ilgili. Çok teşekkürler.
