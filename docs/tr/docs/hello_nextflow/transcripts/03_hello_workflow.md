# Bölüm 3: Merhaba İş Akışı - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Önemli notlar"

    Bu sayfa yalnızca transkripti göstermektedir. Adım adım talimatların tamamı için [kurs materyaline](../03_hello_workflow.md) geri dönün.

    Transkriptte gösterilen bölüm numaraları yalnızca bilgilendirme amaçlıdır ve materyaldeki tüm bölüm numaralarını içermeyebilir.

## Hoş Geldiniz

Merhaba, "Merhaba Nextflow" eğitim kursunun üçüncü bölümüne hoş geldiniz.

Bu bölümün adı "Merhaba İş Akışı".

İkinci bölümde, tek bir işlemden oluşan basit bir iş akışı oluşturduk, ancak gerçekte, boru hatları birden fazla analiz adımını birbirine zincirleme olarak bağlayabildikleri için faydalıdır.

Bu bölümde, o başlangıç örneğini alıp biraz daha gerçekçi olacak şekilde genişleteceğiz.

Bazı ek adımlar ekleyeceğiz ve bu adımları bağlamak için kanalları nasıl kullandığımıza bakacağız.

Tek bir işlemde birleşebilen birden fazla göreve bakacağız ve birden fazla girdiye ve birden fazla çıktıya sahip olabilen işlemlere göz atacağız.

Pekala, hadi başlayalım.

Baştan başlayalım. Daha önceki gibi. training.nextflow.io'ya gidelim. Merhaba Nextflow, üçüncü bölüm. Merhaba İş Akışı. Ve çalışma alanımızı açalım. Önceki bölümlerden tüm çalışma dosyalarımı temizledim ve Merhaba İş Akışını açacağım.

Şimdi bu, şimdiye kadar üzerinde çalıştığımız aynı dosya, bu yüzden tanıdık gelmelidir. Say hello işlemimiz var. greetings CSV dosyasıyla params.greeting'imiz var ve altta, o CSV dosyasını yükleyen, kanalı oluşturan ve işlemimize ileten iş akışımız var.

## 0. Isınma: hello-workflow.nf'yi Çalıştırma

İsterseniz, bunu deneyebilir ve beklediğimiz gibi çalıştığını iki kez kontrol edebiliriz. nextflow run hello workflow nf için bir terminal açın ve enter'a basın.

Tamam, harika. Üç işlemimiz çalıştı. Üç çıktımızla birlikte results dizinimiz var. Bonjour. Hello. Holà. O yüzden bu dosyaları kapatalım, terminali kapatalım, betiğe geri dönelim.

## 1. İş Akışına İkinci Bir Adım Ekleme

Tamam. Örneğimiz için, temel kalmaya çalışıyoruz ve alan bağımsız olmaya çalışıyoruz. Bu yüzden ikinci işlemimiz sadece bu dizeleri, bu kelimeleri, basit bir şekilde manipüle edecek. Bu dosyaları alıp hepsini büyük harf yapmak için Unix'in translate komutunu kullanacağız. Bunu "tr" komutuyla yapıyoruz.

## 1.1. Büyük Harf Yapma Komutunu Tanımlama ve Terminalde Test Etme

Bunu sadece bash terminalinde deneyebilir ve çalışıp çalışmadığını görebiliriz. Echo, Hello World yapıyorsunuz ve sonra bunu pipe karakteriyle tr'ye iletiyorsunuz, ve ona bir tanıma deseni veriyoruz, a'dan z'ye ve neye dönüştürmesi gerektiği. Büyük harfle A'dan Z'ye.

Bu çok basit çünkü tam anlamıyla A'dan Z karakterlerini yapıyor. Bu yüzden aksanlı veya bunun gibi herhangi bir şey üzerinde çalışmayacaktır. Ama örneğin amaçları için, resmi anlamalısınız.

Enter'a basacağım ve terminale büyük harflerle HELLO WORLD yazdırıyor. Ve daha önce olduğu gibi, istersek bunu bir dosyaya yönlendirebiliriz. Outfile.

Tamam. Bunu temizleyelim.

## 1.1. Büyük Harf Yapma Adımını Nextflow İşlemi Olarak Yazma

Betiğimize geri dönelim ve bu bash komutunu yönetmek için yeni bir işlem yazalım. Önceki işlemi kopyalayacağım, altına yapıştıracağım ve buna convert to upper diyeceğim. Büyük harf için. publishDir results'ı aynı şekilde kullanacağım, ancak burada birkaç değişiklik yapacağım. Bir val almak yerine, path input file alacağım ve çıktı dosyalarımızın üst üste binmemesi için burada bir upper öneki kullanacağım. Ve girdiden değişken adını kullanacağım. Sonra aşağıdaki betiği değiştireceğim ve bunun yerine girdi dosyasında cat kullanacağım ve Bash TR'de yaptığımız gibi, a-z, upper input file .txt. Tamam, kaydet'e tıklayalım.

## 1.2. Workflow Bloğuna Yeni İşlem Çağrısı Ekleme

Şimdi aşağı kaydırırsam, bu işlemi gerçekten çağırmamız gerekiyor. Betiğe sadece işlemi eklemek yeterli değildir. Nextflow'a bu işlemi çalıştırmamız gerektiğini ve nerede yapacağımızı söylemeliyiz.

Buraya convert to upper yazacağım ve

tamam, burada bir argüman beklediğini söyleyen bir hata alıyoruz. Elbette, bu işlemin gerçekten yapacak bir şeyi olması için bu işleme bir şeyler iletmemiz gerekiyor.

## 1.3. İlk İşlemin Çıktısını İkinci İşleme İletme

Yapacağımız şey bu işlemden çıktıyı almak. Yani adı alıyorum, say hello, ve dot out yaptığımda.

Bunun gibi basit bir örnek için, sadece bir çıktısı olan ve bunu yeni bir işleme ilettiğimiz bir işleme sahip olduğumuzda, bir girdisi var, ihtiyacımız olan tek şey bu olmalı. Kaydet'e tıklayacağım, terminali açacağım ve bunu tekrar çalıştırmayı deneyelim.

## 1.4. İş Akışını Tekrar Çalıştırma

Şimdi, bu iş akışını en son çalıştırdığımdan beri work dizinini temizlemedim. Tekrar çalıştıracağım ve bunu kısmi önbelleklemenin nasıl çalıştığını göstermek için bir fırsat olarak kullanacağım. Eğer tek tire resume yaparsam. Umarım, en son çalıştırdığımla tamamen aynı olan o ilk işlemden gelen çıktıları yeniden kullanmalıdır. Ama şimdi daha önce çalışmamış, sıfırdan çalışan yeni bir işlemimiz var. Ve elbette, ilk işlemin önbellek çıktılarını kullandığını ve ikinci çıktının üçte üçünü çalıştırdığını görebilirsiniz. Ayrıca şimdi her iki işlemimizin de burada olduğunu görebilirsiniz, ilk işlemimiz, say hello, üç kez çalıştı ve ikinci işlemimiz convert to upper üç kez çalıştı.

Bunu tekrar çalıştırırsam, hatırlatma olarak, -ansi-log false ile, altı farklı işlem görevinin çalıştığını görmeliyiz, her biri için üçer. Bu tam olarak umduğumuz şeyi yapıyor. İlk işlem üç kez çalışıyor, bu çıktıları ikinci bir işleme iletiyor, o da sonra üç kez çalışıyor.

Hadi work dizininin içine bakalım ve Nextflow'un bu dosya girdilerini nasıl yönettiğini görelim. Eğer buradaki ikinci işlemden bu hash dizinini alırsam, bu dosyalara bakmak için tree komutunu tekrar -a ile kullanabiliriz. Burada girdi dosyamız olan Bonjour-output.txt dosyasını görebilirsiniz ve bu aslında bir symlink. Bu okun bize gösterdiği şey bu ve önceki work dizinindeki dosyaya işaret ediyor.

Bu mantıklı. Nextflow her görevin yürütülmesini kendi kapsüllenmiş dizininde yönetir, böylece tamamen kendi kendine yeterlidir. Ancak, önceki adımlardan gelen dosyaları girdi olarak sağlaması gerekiyor. Work dizininin dışına uzanmak yerine bu dosyaları almak için, Nextflow onları work dizinine aşamalar.

Burada olduğu gibi paylaşımlı bir dosya sistemimiz varsa, bunu ek dosya alanı kullanmaması için bir symlink kullanarak yapar. Farklı konumlarda bucket'lı bulut depolaması kullanırsak, o dosyaları getirir ve gerçekten work dizinine kopyalar.

command sh dosyasına bakalım. Eğer code work, command sh yaparsam, elbette, o dosyaya yerel dizinden eriştiğini görebilirsiniz. Yani her şey çok kendi içinde ve temiz.

Ayrıca results dizinini kontrol edebilir ve bu dosyaların düzgün şekilde çıktılandığından emin olabiliriz. Ve elbette, results'ta, ilk işlemden gelen tüm çıktı dosyalarını ve ikinciden gelen tüm çıktı dosyalarını görebiliriz. Ve umduğumuz gibi hepsi büyük harf.

Nextflow'un gücünün parlamaya başladığı yer burası. Çok az kodla ve Nextflow, bu görevlerin paralel olarak yürütülmesini, ayrı work dizinlerinde temiz kapsülleme ile, girdi ve çıktı dosyalarını aşamalama ve dosya yayımlama işlemlerini otomatik olarak bizim için sadece kutunun dışında yönetti. Yani, bu analiz iş akışlarımızın karmaşıklığını ölçeklendirdikçe, bu işlevselliğin gerçekten çok değerli olduğunu görebilirsiniz.

## 2. Tüm Selamlamaları Toplamak İçin Üçüncü Bir Adım Ekleme

Tamam. Bu adımlar bire-bir idi. İlk işlemden ikinci işlem için bir girdiye giden bir çıktımız vardı. Sonra, bu farklı çıktıları tek bir işlem görevinde nasıl toplayacağımızdan bahsedeceğiz, ki bu da yine çok yaygın bir şey. O yüzden hızlıca terminali açalım ve bunun bir kuru çalıştırmasını yapalım.

## 2.1. Toplama Komutunu Tanımlama ve Terminalde Test Etme

Hile yapacağım ve eğitim materyalinden örnek bash kodunu kopyalayıp enter'a basacağım.

Burada gördüğümüz, bu echo komutunu üç farklı çıktı dosyasına üç kez çalıştırdık, bunu burada görebilirim. Ve sonra bu üç farklı dosyanın her birinin çıktısını yazdırmak için cat komutunu kullandık ve bunu tek bir toplanmış dosyaya yönlendirdik.

Ve eğer "cat COLLECTED-output" yaparsam, bu üç farklı dosyanın içeriğine sahip olduğunu, şimdi tek bir dosyada olduğunu görebilirsiniz.

## 2.2. Toplama Adımını Yapmak İçin Yeni Bir İşlem Oluşturma

Öyleyse aynı şeyi Nextflow boru hattımızda kopyalayıp kopyalayamayacağımıza bakalım.

Yukarı kaydıralım ve üçüncü bir işlem oluşturalım. Bu önceki olanı kopyalayacağım ve bu sefer buna Collect Greetings diyeceğim.

Bash terminalinde, buna collected output txt dedik. Yani burada aynı path output'u söyleyeceğim. Ve aynı şekilde kaydedilmesi için yönlendirmeyi burada yapacağım.

Tamam. O komutun başında ne olduğunu değiştirmemiz gerekiyor ve burada girdi dosyasının ne olduğunu düşünmemiz gerekiyor. Aslında, bu işlem birden fazla girdi dosyası alacak. Path'i tutacağım ve bunu yeni bir değişken olan input files, çoğul, olarak değiştireceğim.

Sonra tekrar, bash betiğimizde yaptığımız gibi cat yapacağım. Ve burada değişkeni kullanacağım.

Şimdi, bunun işe yaramayacağını düşünebilirsiniz. Daha önce bir dizi dizenin veya bir dizi path'in bir işleme iletildiği ve bunun hataya neden olduğu başarısızlıklar gördük. Ama aslında, burada Nextflow bunu bizim için otomatik olarak doğru şekilde yönetecek. Birkaç farklı girdi dosyası alacak ve burada sadece farklı dosya yollarını yazdıracak.

Tabii ki, cat komutunun bunun gibi bir dizi dosya adı alabiliyor olması yardımcı oluyor. Eğer her dosya yolundan önce bir argüman gerektiren veya bunun gibi bir şey gerektiren farklı bir komut kullanıyor olsaydım, bu dosya yollarının yinelemesini işleyebilmek için burada biraz daha kod ve mantığa sahip olmamız gerekirdi. Ama bu durumda, sadece çalışmalı.

## 2.3. İş Akışına Toplama Adımını Ekleme

Tamam, iş akışına gidelim ve yeni işlemimizi ekleyelim. Collect greetings. Ve tekrar, convert to upper out'tan çıktı alalım. Bunu kaydedelim.

Bir deneyin. nextflow run hello workflow.

Tamam, iş akışı çalıştı, ama burada biraz garip bir şey var. İlk adımın üç yürütmesi var, beklediğimiz gibi. İkinci için üç görev, ama sonunda tüm çıktıları birleştirmek için burada sadece tek bir görev beklediğimizde üç görevimiz var.

Results dizinimize girersek. Ayrıca collected output'un üç yerine sadece tek bir değere sahip olduğunu görürüz. Bunun nedeni, o çıktı dosyasının üç farklı değerle üç kez üzerine yazılmış olması.

Bu mantıklı çünkü önceki adımda yaptığımız gibi bir çıktıyı bir girdiye burada iletiyoruz.

## 2.4. Selamlamaları Tek Bir Girdide Toplamak İçin Bir Operatör Kullanma

O yüzden burada üç elemanlı bu kanalı alıp bunları tek bir elemana daraltmak için bir operatöre ihtiyacımız var, böylece o son işlem sadece bir kez çalışır.

Bunu yapmak için, collect operatörünü kullanacağız. Bunu doğrudan iş akışı içinde yapabilirim. .out yapabilirim ve sonunda bir operatöre zincirleyebilirim .collect.

Kaydet'e basın. Ve sonra bu eğitimin amaçları için, daha önce yaptığımız gibi bazı view operatörleri de yapacağım, böylece iş akışını çalıştırdığımızda komut satırında görebiliriz, böylece ne olduğunu anlayabiliriz.

Bu kanalı alacağım, collect'i kaldıracağım ve dot view greetings yapacağım, ve sonra bu satırı kopyalayacağım, collect operatörünü ekleyeceğim. Ve bunu after olarak değiştireceğim.

Bu, bunu çağırdığımız yerden ayrı, ama sorun değil çünkü aynı çıktı kanalında aynı operatör çağrılarını kullanıyoruz.

Tamam, kaydet'e basalım ve terminalde deneyelim. nextflow run yapacağım. Hello, workflow. Betiğimizi yeniden çalıştır.

Tamam. Bu daha iyi görünüyor. Daha önce olduğu gibi ilk iki işlemin üç kez çalıştığını görebiliriz ve şimdi son işlemimiz sadece bir kez çalıştı.

View operatörü tarafından yazdırılanı bakarsak, aşağıda, collect'ten önce dedik, bu buradaki çıktı, ve bu üç kez yazdırılıyor. Ve bunların her biri için tek bir path olduğunu görebilirsiniz. Ve sonra collect'ten sonra, üç path'lik bu diziye sahip olduğumuzu görebilirsiniz. Yani beklediğimiz gibi.

Tamam, results dosyasını kontrol edelim ve bu sefer beklediğimiz gibi olup olmadığına bakalım. Elbette, dosyada şimdi üç satır var - bu üç çıktıyı başarıyla tek bir çıktı dosyasında birleştirdi. Fantastik.

Tamam, temizleyeceğim ve bir sonraki adıma geçelim. Ve işleri temiz tutmak için bu view ifadelerini sileceğim.

## 3. Nihai Çıktı Dosyasını Benzersiz Şekilde Adlandırmak İçin Bir İşleme Birden Fazla Girdi İletme

Tamam. Şimdiye kadar, tüm işlemlerimiz sadece tek bir girdi aldı. Şimdi bunun nasıl çalıştığını görmek için bir işleme birden fazla girdi eklediğimiz bir alıştırma yapacağız. Bunu yapmak için, bu collect greetings örneğini kullanacağız.

İş akışını her çalıştırdığımda, results dizinindeki o dosyanın üzerine yazdı, ki bu istediğimiz şey olmayabilir.

## 3.1. Çıktı Dosyası İçin Kullanıcı Tanımlı Bir Ad Kabul Edecek Şekilde Toplayıcı İşlemi Değiştirme

Bu örnek için, çıktı dosyası adını özelleştirebilmemiz için ek bir parametre ileteceğiz.

Bir işleme ikinci bir girdi eklemek çok basit. Input bloğunda sadece ikinci bir satır eklerim. Bu sefer bir path yerine bir value olacak, çünkü bir dize iletmek istiyoruz ve buna batch underscore name diyeceğim.

Artık bu değişkeni script bloğunda kullanabilirim ve collected tire dollar batch name diyeceğim.

Burada değişken adının etrafında süslü parantezler kullanıyorum. Bu sadece onu dizenin geri kalanından ayrı tutmak için, ve bu durumda muhtemelen gerekli değil, ama bence okunmasını kolaylaştırıyor.

Tamam. Son olarak, output path'i güncellemeyi unutmayın çünkü artık dosya adı değişti, bu yüzden aynı şeyi yapacağım ve beklendiği gibi batch name'i output path'ine koyacağım.

## 3.2. Komut Satırı Batch Parametresi Ekleme

Şimdi bir yerden batch name iletmemiz gerekiyor ve bunu iş akışını çalıştırdığımızda komut satırında yapabilmemiz için ikinci bir parametre oluşturacağım.

Yani params batch name yapacağım ve varsayılan olarak, buna test batch diyelim. Şimdi bu özel parametreler değişkenini aşağıda, işlemi çağırdığımız yerde kullanabilirim.

Ve elbette VS Code bize şimdi bu işleme yeterli argüman olmadığını ve ikinci bir girdi beklediğini söylüyor.

Basitçe virgül yapın ve yeni değişkenimizi iletin ve hata kayboluyor.

Buradaki girdilerin sırası gerçekten önemlidir. İlk işlem girdisi path idi ve ikinci girdi ad. Eğer burada sırayı değiştirirsem, işlemi çağırdığımda da sırayı değiştirmeliyim. Aksi takdirde. Sonraki, yanlış kanalı yanlış girdiye ileteceğiz.

## 3.3. İş Akışını Çalıştırma

Tamam, deneyelim ve çalışıp çalışmadığına bakalım. "nextflow run hello- workflow" yapalım. Tamam, daha önce olduğu gibi çalıştı. Results dizinine bakalım.

Elbette, dosya adımız şimdi "collected test batch output txt" olarak adlandırılıyor. Fantastik.

Ve şimdi tekrar çalıştırarak bunun üzerine yazıp yazamayacağımızı görelim. Bu sefer --batch_name yapacağım, buradaki özel parametre değişken adıyla eşleşecek şekilde. Ve buna demo output diyeceğim.

İş akışını tekrar çalıştırın ve bir şey olup olmadığını göreceğiz.

Tamam, şimdi collected demo output .txt var. Ve bu dosya adı ondan farklı olduğu için, üzerine yazmadı. Her ikisi de şimdi results dizininde mevcut.

## 4. Toplayıcı Adıma Bir Çıktı Ekleme

Tamam, orada bir işleme birden fazla girdi vermeyi gösterdik, ama birden fazla çıktıya ne dersiniz? Bu örnek için, işlenen selamlama sayısını hesaplayacağız ve bunu bu collect greeting adımı için ikincil bir çıktı olarak çıktılayacağız.

## 4.1. Selamlama Sayısını Saymak ve Çıktılamak İçin İşlemi Değiştirme

Burada biraz hile yapacağız. Nextflow işlemleri çok satırlı bir dizeyle bu script bloğuna sahiptir ve bu, dot command dot sh'ye bash çıktısı olarak iletilir. Ama aslında bunun üzerinde herhangi bir özel kod yazabiliriz ve bu, görevin bir parçası olarak yürütülecek ancak bash betiğine dahil edilmeyecek.

Nextflow sözdizimindeki yerleşik işlevlerden biri size olarak adlandırılır. Bu yüzden path input'u alacağım ve count underscore greetings diyeceğim, sadece bir değişken adı tanımlamak için. Input files'ı alacağım ve üzerinde "size" çağıracağım.

Bu işlev, bu girdi kanalının boyutunu sayacak ve bir değişkene atayacak.

Artık bu değişkeni output bloğunun bir parçası olarak döndürebiliriz. Yani val diyoruz, çünkü bir dosya değil value. Ve count greetings.

Şimdi bu kendi başına yeterli ve artık bu işlemden bu farklı çıktılara erişebiliriz. Ancak, onlara konumsal bir şekilde erişmemiz gerekir. Yani sıfır ve bir gibi bir dizin anahtarı kullanarak.

Çıktılara ulaşmayı biraz daha kolaylaştırmak için, onları adlandırabiliriz ve bunu bir emit ifadesi kullanarak yapıyoruz.

Yani virgül emit out file yapıyoruz veya bunu her ne adlandırmak istersem. Ve burada emit count yapıyorum. Bu temelde sadece bir dekoratördür, bu sadece biraz daha temiz kod yazmamıza yardımcı olur, böylece daha sonra workflow bloğunda belirli çıktılara kolayca başvurabiliriz.

## 4.2. İş Akışının Sonunda Çıktıyı Raporlama

Tamam. Eğer workflow bloğuna kaydırırsam, artık collect greetings'in çıktılarını alabilirim, collect greetings, dot out yapın ve VS Code uzantısı tarafından burada önerilen iki adlandırılmış çıktımızı görebiliriz. Çok kullanışlı.

Bu yüzden방금 oluşturduğumuz count değerini almak için dot count yapacağım ve iş akışını çalıştırdığımızda komut satırında yazdırması için view yapacağım. Yani çalıştırdığımızda görebiliriz.

Biraz daha güzel yapmak için closure içine bir şeyler yazalım. num greetings, greetings selamlama vardı.

Ve aslında diğer çıktıyı umursamıyoruz çünkü bunu başka hiçbir işlem için girdi olarak kullanmıyoruz. Ama bunu istersek, aşağıda başka bir işleme girdi olarak kolayca nasıl iletebileceğimizi görebilirsiniz.

## 4.3. İş Akışını Çalıştırma

Kaydet'e tıklayacağız. Terminale bakalım ve deneyelim.

Tamam, fantastik. İşte burada. Üç selamlama var. Bu tamamen doğru.

Tamam, harika. Bu bölümün sonu. Buraya kadar geldiğiniz için hepiniz bitirdiniz. Artık oldukça gerçekçi bir iş akışı oluşturmaya başlıyorsunuz, burada girdileri ve çıktıları ve iş akışımızdaki mantığı yönetebiliyoruz.

Bu iş akışı dosyaları uzadıkça, biraz hantal hale gelmeye başlıyorlar. Bu yüzden bir sonraki bölümde, iş akışı içindeki kodu bulmanın ve sürdürmenin daha kolay olması için Nextflow kodunu ayrı dosyalara nasıl modülerleştirebileceğimize bakacağız.

Bir sonraki videoda, dördüncü bölüm için bize katılın. Merhaba Modüller.

[Sonraki video transkripti :octicons-arrow-right-24:](04_hello_modules.md)
