# Bölüm 4: Hello Modules - Video Transkripti

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım tam talimatlar için [kurs materyaline](../04_hello_modules.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba ve Hello Nextflow'un dördüncü bölümüne tekrar hoş geldiniz. Bu bölüm tamamen modüllerle ilgili ve kursun oldukça kısa bir bölümü. Aslında fazla kod yazmayacağız, daha çok pipeline'daki kodu nasıl organize edeceğimizle ilgili.

Şimdiye kadar her şeyi tek bir dosyaya koyuyorduk, bu iyidir ve aslında eski günlerde Nextflow pipeline'larını böyle oluşturuyorduk.

Ancak bu pipeline büyüdükçe, betik daha uzun ve daha uzun ve daha uzun hale geliyor ve gezinmesi, bakımını yapmak daha zor oluyor, ayrıca kodun herhangi birini gerçekten paylaşamayacağımız anlamına geliyor.

Nextflow modülleri, süreçleri ana betikten çıkarmamıza ve ardından bunları içe aktarmamıza olanak tanır. Bu, kodun gezinmesinin daha kolay olduğu ve ayrıca bu modül kodunu farklı pipeline'lar arasında paylaşabileceğimiz anlamına gelir.

Dokümantasyonun ana sayfasındaki bu küçük diyagram konsepti güzel bir şekilde gösteriyor. Tek bir devasa betik yerine, farklı modül betiklerinden bu ayrı modül dosyalarını ekleyeceğiz ve hepsi workflow'a çekilecek, ancak hala tamamen aynı şekilde çalışacak.

O halde GitHub Codespaces'e atlayalım ve etrafta biraz bakalım. Daha önceki gibi, buradaki çalışma alanımı biraz temizledim. Eski Nextflow dizinlerini ve work dizinini ve benzeri şeyleri kaldırdım. Ancak bu dosyalar hala etrafta olsa bile önemli değil.

hello modules dosyasında çalışmaya başlayacağım, bu temelde önceki bölümün sonunda bıraktığımız yerdir. Burada üç sürecimiz var. Birkaç params'ımız var, bu üç süreci çalıştırdığımız ve kanallarla birbirine bağladığımız workflow bloğu. Sonra çıktı kanallarını yayınlıyoruz ve bu dosyaların nasıl yayınlanacağını söyleyen output bloğumuz var.

## 1. Modülleri saklamak için bir dizin oluşturun

Şimdi, dediğim gibi, gerçekten çok fazla kod yazmayacağız veya düzenlemeyeceğiz. Sadece zaten sahip olduğumuz kodu etrafta taşıyacağız. Nextflow modül dosyaları tipik olarak tek bir süreç içerir ve geleneksel olarak bunları normalde modules adlı bir dizinde tutarız. Ancak buna istediğiniz şeyi verebilirsiniz. Ama buradaki repomda bir modules dizini tutacağım ve sonra her süreç için bir dosya oluşturacağım. Yani yeni dosya diyeceğim, sayHello.nf.

## 2. sayHello() için bir modül oluşturun

Şimdi sürecimi alacağım ve bu kodu ana hello modules dosyasından seçeceğim, keseceğim ve buraya yapıştıracağım.

Açıkçası bu tek başına bir şey yapmaz. Ana betiğimiz hala bu sürece ihtiyaç duyuyor, bu yüzden onu bir şekilde geri çekmemiz gerekiyor. Ve bunu include ifadesiyle yapıyoruz.

Yani include yazıyorum ve süslü parantezler, ve sonra sürecin adını alıyorum. Ve from diyorum, ve sonra ona göreceli bir dosya yolu veriyorum. Yani bu betiğin kaydedildiği yerden göreceli olduğu için ./ ile başladığını söylüyor. Yani modules sayHello.nf.

VS code uzantısının burada oldukça yardımcı olduğuna dikkat edin. Bize bu dosyayı bulabiliyorsa ve adlandırdığım bir süreç bulabiliyorsa söylüyor. Burada kasıtlı olarak bir yazım hatası koyarsam, hemen bana bir hata veriyor ve içe aktarmaya çalıştığım bu süreci bulamadığını söyleyecek. Bu yüzden bulduğunuz hatalara göz kulak olun.

Ve gerçekten hepsi bu. Sürecimiz hala burada. Aşağıda hiçbir değişiklik gerekmiyor. Süreç aynı ada sahip ve tamamen aynı şekilde yürütülüyor. Sadece sürecin asıl kodu artık ayrı bir dosyada.

Nextflow workflow'unu tekrar çalıştırabiliriz, tamamen aynı şekilde çalışacak. Ve bu temelde kursun bu bölümünün geri kalanı sadece bu üç süreci kendi dosyalarına taşımak.

Şimdi bunu yapalım. İkinci süreç için hızlıca yeni bir modül dosyası oluşturacağım: convertToUpper.nf. O kodu keseceğim, buraya yapıştıracağım. Ve sonra bunu da dahil edeceğim. hadi, harika.

Ve sonra collectGreetings.nf için yeni bir dosya oluşturacağım. Bunu kesin.

Bir sürü kesme, kesme ve kopyalama ve yapıştırma.

Ve şimdi ana workflow betiğimiz aniden çok daha kısa, çok daha ulaşılabilir ve okuması çok daha kolay görünüyor.

Ve projenin artık farklı dosyalarımızla nasıl inşa edilmeye başladığını görebilirsiniz. İstediğimiz yerlerde detaya dalabilir, pipeline'daki belirli adımları bulmak için etrafta çok daha kolay gezinebilir ve pipeline'ın ne yaptığına dair hızlı bir genel bakış alabiliriz.

## VS Code ile modüllerde gezinme

Şimdi, elbette, bunu yapmanın dezavantajı, büyük bir pipeline'ınız varsa, bir sürü modül dosyanız olacak ve bunlar birden fazla alt dizinde organize edilebilir veya her türlü şey olabilir. Şimdi, yine, burada küçük bir ipucu. VS Code uzantısı kod tabanınızda gezinmekte ve ayrıca oradaki kod hakkında size bilgi vermekte oldukça iyidir.

VS Code'un bu sürecin ne olduğunu anladığını ve üzerine geldiğimde bana küçük bir genel bakış verdiğini görebilirsiniz, böylece kaynak kodunu bulup gitmek zorunda kalmadan, bir workflow'da kullanırken tipik olarak en önemli şey olan girdilerin ve çıktıların ne olduğunu görebiliyorum.

Ayrıca command'ı basılı tutarsam, ben Mac'teyim, ve süreç adına tıklarsam, dosyayı doğrudan hemen açar. Onu içeri çeker. Böylece gerçek dosya yollarının ne olduğunu düşünmeden doğrudan oraya atlayabilirim. Ve bu her yerde çalışır, süreçlerin çağrıldığı yerlerde de bunu yapabilirim. Yani bu gerçekten hızlı.

## 4.4. Workflow'u çalıştırın

Tamam, pipeline'ın hala beklediğimiz gibi çalıştığını kontrol edelim. Terminali açalım. "nextflow run hello modules" yapalım ve herhangi bir sorun olmadan çalışıp çalışmadığını görelim.

Umarım bunun amacı pipeline'ın temelde değişmemiş olmasıdır, bu yüzden gerçekten daha önce çalıştırdığımızda gördüğümüz herhangi bir değişikliği görmemelisiniz. Buradaki çıktı tamamen aynı görünüyor ve tüm aynı dosyalarla results dizinimizi görebilirsiniz, bu harika. Değişiklik olmaması iyi.

## nf-core/modules hakkında bir not

Bitirmeden önce, modüller söz konusu olduğunda işbirliğinin gücüne hızlıca değinmek istiyorum. Bu dosyalar benim repomda oturuyor, bu yüzden üzerlerinde nasıl işbirliği yapabileceğimiz hemen belli olmuyor. Ve bunu yapabileceğiniz birçok farklı yol var, ancak muhtemelen bunun en büyük ve en iyi bilinen örneği nf-core.

Eğer nf-core web sitesine gidersem, resources'a, ve modules'e gidiyorum. nf-core'un devasa bir modül kütüphanesine sahip olduğunu görebilirsiniz, bunu görüntülediğimde 1700'ün hemen altında modül var. Ve böylece en sevdiğim araçlardan herhangi birinin adını yazabilirim, başka birinin zaten bunun için bir modül yazıp yazmadığını bulmaya gidebilirim ve burada önceden yazılmış bu modül sürecini tüm girdilerle, çıktılarla, yazılım container'larıyla, tüm bu bilgilerle görebilirim ve burada yanda kaç farklı nf-core pipeline'ının bu tek paylaşılan süreci kullandığını görebilirsiniz.

Bu biraz aşırı bir örnek, ancak bunun gerçekten bu kodu yeniden kullandığını görebilirsiniz. Ve bunun GitHub kaynağına tıklarsam, yaptığımızla tamamen aynı. Sadece bir dosyada büyük bir süreç.

Şimdi nf-core tarafında, bu dosyaları paylaşabilmek ve bunları farklı repolara getirebilmek için bazı numaralar yapıyoruz. Ve bu konuda daha fazla bilgi edinmek istiyorsanız, özellikle nf-core ile kullanma ve oluşturma hakkında sahip olduğumuz kursa göz atın. Ancak size bu kod yeniden kullanımı konseptinin ne kadar güçlü olabileceğine dair bir fikir vermek istedim.

## Özet

Tamam, modüller için bu kadar. Size kursun kısa bir bölümü olduğunu söylemiştim. Testi kontrol edin, anladığınızdan emin olun ve her şeyin hala düzgün çalıştığından emin olun. Ve sizi bir sonraki videoda göreceğim, o video tamamen yazılım container'larıyla ilgili. Çok teşekkür ederim.

I.
