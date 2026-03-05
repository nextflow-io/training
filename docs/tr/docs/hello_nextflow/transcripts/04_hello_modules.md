# Bölüm 4: Hello Modules - Video Transkripti

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../04_hello_modules.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba ve Hello Nextflow'un dördüncü bölümüne tekrar hoş geldiniz. Bu bölüm tamamen modüllerle ilgili ve kursun oldukça kısa bir bölümü. Aslında pek fazla kod yazmayacağız, daha çok pipeline'ımızdaki kodu nasıl organize ettiğimizle ilgili.

Şimdiye kadar her şeyi tek bir dosyaya koyuyorduk, bu da sorun değil ve aslında eski günlerde Nextflow pipeline'larını böyle oluşturuyorduk.

Ancak pipeline büyüdükçe, betik gittikçe uzuyor ve gezinmesi, bakımı zorlaşıyor ve ayrıca kodun herhangi bir bölümünü gerçekten paylaşamıyoruz.

Nextflow modülleri, süreçleri ana betikten çıkarmamıza ve ardından içe aktarmamıza olanak tanır. Bu, kodun gezinmesinin daha kolay olduğu anlamına gelir ve ayrıca bu modül kodunu farklı pipeline'lar arasında paylaşabileceğimiz anlamına gelir.

Dokümantasyonun ana sayfasındaki bu küçük diyagram konsepti güzel bir şekilde gösteriyor. Tek bir devasa betik yerine, farklı modül betiklerinden bu ayrı modül dosyalarını dahil edeceğiz ve hepsi iş akışına çekilecek, ancak yine de tamamen aynı şekilde çalışacak.

O halde GitHub Codespaces'e atlayalım ve biraz etrafta dolaşalım. Daha önce olduğu gibi, çalışma alanımı burada biraz temizledim. Eski Nextflow dizinlerini ve work dizinini vb. kaldırdım. Ancak bu dosyalar hala etrafta olsa da sorun değil.

hello modules dosyasında çalışmaya başlayacağım, bu temelde önceki bölümün sonunda bıraktığımız yer. Burada üç sürecimiz var. Birkaç parametre var, workflow bloğu var, burada bu üç süreci çalıştırıyoruz ve bunları kanallarla birbirine bağlıyoruz. Ardından çıktı kanallarını yayınlıyoruz ve bu dosyaların nasıl yayınlanacağını söyleyen output bloğumuz var.

## 1. Modülleri saklamak için bir dizin oluşturun

Şimdi, dediğim gibi, gerçekten pek fazla kod yazmayacağız veya düzenlemeyeceğiz. Sadece zaten sahip olduğumuz kodu taşıyacağız. Nextflow modül dosyaları genellikle içlerinde tek bir süreç bulundurur ve geleneksel olarak bunları modules adlı bir dizinde tutarız. Ancak buna istediğiniz ismi verebilirsiniz. Ama burada depomda bir modules dizini tutacağım ve ardından her süreç için bir dosya oluşturacağım. Yani yeni dosya diyeceğim, sayHello.nf.

## 2. sayHello() için bir modül oluşturun

Şimdi sürecimi alacağım ve bu kodu seçeceğim, ana hello modules dosyasından kesip buraya yapıştıracağım.

Açıkçası bu tek başına bir şey yapmaz. Ana betiğimiz hala bu sürece ihtiyaç duyuyor, bu yüzden onu bir şekilde geri çekmemiz gerekiyor. Ve bunu include ifadesiyle yapıyoruz.

Yani include yazıyorum ve bazı süslü parantezler, ve sonra sürecin adını alıyorum. Ve from diyorum, ve sonra ona göreceli bir dosya yolu veriyorum. Yani diyor ki, ./ ile başlıyor çünkü bu betiğin kaydedildiği yerden göreceli. Yani modules sayHello.nf.

VS code uzantısının burada oldukça yardımcı olduğuna dikkat edin. Bize bu dosyayı bulup bulamadığını ve adını verdiğim bir süreç bulup bulamadığını söylüyor. Burada kasıtlı olarak bir yazım hatası yaparsam, hemen bana bir hata veriyor ve içe aktarmaya çalıştığım bu süreci bulamadığını söylüyor. Bu yüzden bulduğunuz hatalara göz kulak olun.

Ve gerçekten bu kadar. Sürecimiz hala burada. Aşağıda hiçbir değişiklik yapılmasına gerek yok. Süreç aynı isme sahip ve tamamen aynı şekilde yürütülüyor. Sadece sürecin gerçek kodu artık ayrı bir dosyada.

Nextflow iş akışını tekrar çalıştırabiliriz, tamamen aynı şekilde çalışacak. Ve bu temelde kursun bu bölümünün geri kalanı sadece bu üç süreci kendi dosyalarına taşımak.

O halde şimdi bunu yapalım. İkinci süreç için hızlıca yeni bir modül dosyası oluşturacağım: convertToUpper.nf. Bu kodu kesip buraya yapıştıracağım. Ve sonra onu dahil edeceğim. hadi, harika.

Ve sonra collectGreetings.nf için yeni bir dosya oluşturacağım. Onu kesip çıkaracağım.

Çok fazla kesme, kesme ve kopyalama ve yapıştırma.

Ve şimdi ana iş akışı betiğimiz aniden çok, çok daha kısa, çok daha ulaşılabilir ve okunması çok daha kolay görünüyor.

Ve projenin artık farklı dosyalarımızla nasıl oluşmaya başladığını görebilirsiniz. İstediğimiz yerlerde detaya dalabiliriz. Pipeline'daki belirli adımları bulmak için çok daha kolay bir şekilde gezinebilir ve pipeline'ın ne yaptığına dair hızlıca bir genel bakış elde edebiliriz.

## VS Code ile modüllerde gezinme

Şimdi, elbette, bunu yapmanın dezavantajı, büyük bir pipeline'ınız varsa, çok sayıda modül dosyanız olacak ve bunlar birden fazla alt dizinde organize edilebilir veya her türlü şey olabilir. Şimdi, yine, burada küçük bir ipucu. VS Code uzantısı kod tabanınızda gezinmede oldukça iyi ve ayrıca size oradaki kod hakkında bilgi veriyor.

VS Code'un bu sürecin ne olduğunu anladığını ve üzerine geldiğimde bana küçük bir genel bakış verdiğini görebilirsiniz, böylece kaynak kodunu bulup gitmek zorunda kalmadan, girdilerin ve çıktıların ne olduğunu görebilirim, bu genellikle bir iş akışında kullanırken en önemli şeydir.

Ve ayrıca command tuşunu basılı tutarsam, Mac kullanıyorum, ve süreç adına tıklarsam, dosyayı doğrudan hemen açar. Çeker. Yani gerçek dosya yollarının ne olduğunu düşünmeden doğrudan oraya atlayabilirim. Ve bu her yerde çalışır, bunu süreçlerin çağrıldığı yerlerde de yapabilirim. Yani bu gerçekten hızlı.

## 4.4. İş akışını çalıştırın

Tamam, pipeline'ın hala beklediğimiz gibi çalıştığını kontrol edelim. Terminali açalım. "nextflow run hello modules" yapalım ve herhangi bir sorun olmadan yürütülüp yürütülmediğini görelim.

Umarım bunun tüm amacı pipeline'ın temelde değişmemiş olmasıdır, bu yüzden gerçekten daha önce çalıştırdığımızdan herhangi bir değişiklik görmemelisiniz. Buradaki çıktı tamamen aynı görünüyor ve tüm aynı dosyalarla results dizinimizi görebilirsiniz, bu harika. Değişiklik olmaması iyi.

## nf-core/modules hakkında bir not

Bitirmeden önce, modüller söz konusu olduğunda işbirliğinin gücüne hızlıca değinmek istiyorum. Bu dosyalar depomda oturuyor, bu yüzden bunlar üzerinde nasıl işbirliği yapabileceğimiz hemen belli değil. Ve bunu yapabileceğiniz birçok farklı yol var, ancak muhtemelen bunun en büyük ve en iyi bilinen örneği nf-core.

nf-core web sitesine gidersem, kaynaklara, ve modüllere gidiyorum. nf-core'un çok büyük bir modül kütüphanesine sahip olduğunu görebilirsiniz, bunu görüntülediğimde neredeyse 1700'ün biraz altında modül var. Ve böylece en sevdiğim araçlardan herhangi birinin adını yazabilirim, başka birinin bunun için zaten bir modül yazıp yazmadığını bulmaya gidebilirim ve burada tüm girdiler, çıktılar, yazılım konteynırları, tüm bu bilgilerle önceden yazılmış bu modül sürecini görebilirim ve yan tarafta burada, kaç farklı nf-core pipeline'ının hepsinin bu tek paylaşılan süreci kullandığını görebilirsiniz.

Bu biraz aşırı bir örnek, ama bunun gerçekten bu kodu yeniden kullandığını görebilirsiniz. Ve bunun için GitHub kaynağına tıklarsam, yaptığımızla tamamen aynı. Sadece bir dosyada büyük bir süreç.

Şimdi nf-core tarafında, bu dosyaları paylaşabilmek ve farklı depolara getirebilmek için bazı numaralar yapıyoruz. Ve bunun hakkında daha fazla bilgi edinmek istiyorsanız, özellikle nf-core ile kullanma ve oluşturma hakkında sahip olduğumuz kursa göz atın. Ancak size bu kod yeniden kullanımı konseptinin ne kadar güçlü olabileceği hakkında bir fikir vermek istedim.

## Özet

Tamam, modüller için bu kadar. Size kursun kısa bir bölümü olduğunu söylemiştim. Testi kontrol edin, anladığınızdan emin olun ve her şeyin hala düzgün çalıştığından emin olun. Ve sizi bir sonraki videoda göreceğim, bu tamamen yazılım konteynırları hakkında. Çok teşekkür ederim.
