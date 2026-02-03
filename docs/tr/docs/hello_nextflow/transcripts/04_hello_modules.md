# Bölüm 4: Merhaba Modüller - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../04_hello_modules.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyallerdeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba, Merhaba Nextflow eğitim kursunun Dördüncü Bölümüne hoş geldiniz.

Bu bölümün adı Merhaba Modüller ve Nextflow kodunun nasıl modülerleştirileceğinden bahsedeceğiz. Yapacağımız şey, tek bir workflow betiğimizi alıp ayrı dosyalara bölmek.

Bu, workflow'unuz büyüdükçe kodun gezinmesini ve bakımını kolaylaştırır ve ayrıca pipeline'lar arasında modüllerin paylaşılmasını mümkün kılar, böylece aynı aracı kullanan birden fazla pipeline'ınız varsa, o process'i yalnızca bir kez yazmanız yeterli olur.

Bunun klasik bir örneği nf-core modül deposudur; kullanıma hazır modüllerde binlerce farklı aracı vardır ve bunları kurabilir ve workflow'unuzda kullanabilirsiniz.

Nextflow ayrıca alt workflow'larla da çalışabilir; bunlar modüller gibidir, ancak birden fazla process içerirler. Bu eğitimin kapsamı dışındadır, ancak temelde aynı şekilde çalışır.

Tamam. Hadi bir göz atalım.

Her zamanki gibi, training.nextflow.io adresine giderek başlayın.

Kenar çubuğunda "Hello Nextflow"a gidin ve dördüncü bölümü yapıyoruz: "Hello Modules".

Şimdi GitHub Code Spaces ortamıma geçeceğim ve "hello-modules" dosyasına bakacağım.

Daha önce olduğu gibi, önceki bölümün son noktasından başlıyoruz, bu yüzden bu betik tanıdık gelmeli. Üç process'imiz var: say hello, convert to upper ve collect greetings, ve bu üç komutu çalıştıran ve sonunda bir mesaj yayan basit bir workflow. Greeting ve batch adında iki parametremiz var; batch, sonunda toplanan çıktı dosyası için kullanılan adı belirtir.

## 0. Isınma: hello-modules.nf'yi çalıştırın

Bu workflow'un beklediğimiz gibi çalıştığını nextflow run hello, modules yaparak doğrulayabiliriz.

Harika. Her bir process'ten üç görev, bir toplama görevi çalıştırdı ve bu toplu işteki üç karşılama olduğunu söyledi. Results'a gidersek, toplanan test batch çıktısı dahil olmak üzere farklı çıktı dosyalarımızı burada aldık.

## 1. Modülleri depolamak için bir dizin oluşturun

Pekala. Hadi biraz modülerleştirme yapalım.

Modülleri pipeline deponuzda bir alt klasöre koymak genellikle iyi bir fikirdir, sadece işleri düzenli tutmak için. Bunu istediğiniz gibi adlandırabilirsiniz, ancak geleneksel olarak genellikle modules diyoruz.

O halde devam edelim, bir terminale gidelim ve modules'u oluşturalım. VS Code'da kenar çubuğunda belirdiğini görebilirsiniz.

## 2. sayHello() için bir modül oluşturun

Ardından ilk modülüm için yeni bir dosya oluşturacağım. "touch" veya "code" yapabilirsiniz veya bunu kenar çubuğunda yapabilirsiniz, gerçekten önemli değil. Bu yüzden code modules yapacağım ve process'in adını vereceğim. Yani sayHello.nf. NF, Nextflow dosyaları için geleneksel dosya uzantısıdır.

Burada kaydet'e basacağım ve yeni modül dosyamızın ortaya çıktığını göreceğiz.

## 2.2. sayHello process kodunu modül dosyasına taşıyın

Sağda, modül kodunu workflow'dan alacağım. Ayrıca önce burada hash bang'i alacağım ve açıkça bir Nextflow dosyası olması için kopyalayacağım. Sonra bu process'i alacağım ve keseceğim. Yani ana workflow betiğimden kaldıracağım ve bu yeni modüle yapıştıracağım.

Bu modül dosyasının içereceği tüm içerik bu. Sadece tek bir process, workflow yok, mantık yok, sadece tek başına bir process.

Şimdi bu dosyayı kapatabilirim.

## 2.3. workflow bloğundan önce bir import bildirimi ekleyin

Şimdi workflow'umda o ilk process eksik, bu yüzden onu içe aktararak geri getirmemiz gerekiyor. Bunun sözdizimi diğer programlama dillerine çok benzer, bu yüzden tanıdık gelebilir. Include'u süslü parantezlerle yaparız, process'in adı say hello, ve sonra dosya yolu modules, say hello, nf'den. Harika.

Burada birkaç püf noktası var. VS Code uzantısı bu konuda akıllı. Bu dosya yolunu tanır ve üzerine gelip follow link yapabilirsiniz. Ya da Mac'teyim, option tuşuna tıklayabilirim ve bu dosyayı açar. Böylece hızlıca ona geçebiliriz.

Bu process adı şimdi aşağıdaki workflow tarafından kullanılıyor ve burada aynı şeyi yapabiliriz. Bize o process hakkında biraz bilgi gösteriyor ve yine option tuşunu basılı tutabilirim, üzerine tıklayabilirim ve editörde açılacak.

Yani farklı process'leriniz için çok sayıda dosyanız olduğunda VS Code'da kod tabanınızda hızlıca gezinmenin gerçekten hızlı bir yolu.

Tamam. Bu bölüm için temelde bu kadar. Şimdi sadece diğer process'ler için aynı şeyi tekrar yapıyoruz.

## 3. convertToUpper() process'ini modülerleştirin

O halde burada yeni bir dosya oluşturalım. Convert to upper nf olarak adlandıralım. Yine hash bang'i kopyalayın. Ve sonra process'i kesin.

Process adını oraya kopyalayın, yeni process adıyla yeni bir include ifadesi ekleyin.

## 4. collectGreetings() process'ini modülerleştirin

Ve sonra üçüncü process için aynısını yapın. Yeni dosya, connect. Greetings,

hash bang'i yapın. Process'i kesin, process'i yapıştırın ve yeni bir include ifadesi yapın.

Şimdi burada geçersiz include kaynağı diyen bir hata alt çizgisi aldığımı görebilirsiniz. Ve bu aslında biraz fazla hızlı hareket ettiğim için yaptığım gerçek bir hata. Yakından bakarsanız, T'yi kaçırdığımı ve convert to upper'a dönüştürdüğümü görebilirsiniz.

Bu yüzden VS Code çok faydalı bir şekilde bana orada bir hata yaptığımı söyledi. O dosya adını düzeltirsem, hata ortadan kalkar. Bu, VS Code içindeki hata kontrolünün Nextflow kodu yazmak için neden bu kadar yararlı olduğunun iyi bir örneği. Aksi takdirde bunu fark etmezdim ve ancak çok daha sonra workflow'u çalıştırmayı denediğimde öğrenirdim.

Ana pipeline betiğimiz artık çok daha basit görünüyor. İçinde hiçbir process yok, sadece üç include ifademiz ve workflow'umuz var. Workflow'un mantığını değiştirmedik. Process kodunu değiştirmedik, bu yüzden umarım tamamen aynı şekilde çalışmalı.

## 4.4. Daha önce olduğu gibi aynı şeyi yaptığını doğrulamak için workflow'u çalıştırın

Kontrol edelim. Bir terminal açacağım ve daha önce olduğu gibi tamamen aynı komutu çalıştıracağım.

Kesinlikle, process'lerimizi çalıştırdı, say hello, convert to upper collect greetings ve bize tekrar üç karşılama verdi.

Yani kodumuzun yerini değiştirdik, ancak workflow'un nasıl yürütüldüğü hakkında hiçbir şeyi değiştirmedik ve tamamen değişmedi. Tek fark, artık daha temiz kodumuz var, bakımı daha kolay ve başkalarıyla paylaşılması daha kolay.

Ve bu kadar. Kısa bir bölümdü. Basit bir kavram, ancak çok güçlü ve daha karmaşık Nextflow workflow'ları yazma şeklimizin anahtarı. Bu yüzden bunu anlamanız ve kullanma alışkanlığı kazanmanız önemli.

Bir sonraki bölümde, biraz tempomuz değişecek ve Nextflow kodu yazma sözdizimi hakkında bu kadar fazla düşünmeyi bırakacağız ve process'lerin kendilerinde yazılımı nasıl kullandığımız hakkında biraz düşüneceğiz. Bölüm beşte Merhaba Konteynerler için bize katılın.

[Sonraki video transkripti :octicons-arrow-right-24:](05_hello_containers.md)
