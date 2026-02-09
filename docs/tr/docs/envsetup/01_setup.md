# GitHub Codespaces

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

GitHub Codespaces, buluttaki sanal makineler tarafından desteklenen, eğitim için önceden yapılandırılmış bir ortam sağlamamıza olanak tanıyan web tabanlı bir platformdur.
Platform, Github (Microsoft'a ait) tarafından işletilmektedir ve Github hesabı olan herkes için ücretsiz olarak (kullanım kotaları ile) erişilebilir durumdadır.

!!! warning "Uyarı"

    Kuruluşlara bağlı hesaplar belirli ek kısıtlamalara tabi olabilir.
    Bu sizin durumunuzsa, bağımsız bir kişisel hesap kullanmanız veya bunun yerine yerel bir kurulum yapmanız gerekebilir.

## GitHub hesabı oluşturma

[GitHub ana sayfasından](https://github.com/) ücretsiz bir GitHub hesabı oluşturabilirsiniz.

## GitHub Codespace'inizi başlatma

GitHub'a giriş yaptıktan sonra, Nextflow eğitim ortamını açmak için tarayıcınızda şu bağlantıyı açın: <https://codespaces.new/nextflow-io/training?quickstart=1&ref=master>

Alternatif olarak, her eğitim kursunda (genellikle Oryantasyon sayfasında) tekrarlanan aşağıda gösterilen düğmeye tıklayabilirsiniz.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Yeni bir GitHub Codespace oluşturabileceğiniz bir sayfa ile karşılaşmalısınız:

![Create a GitHub Codespace](img/codespaces_create.png)

### Yapılandırma

Genel kullanım için herhangi bir şeyi yapılandırmanıza gerek yoktur.
Başladığınız kursta aksi belirtilmedikçe, devam etmek için ana düğmeye tıklamanız yeterlidir.

Ancak, "Change options" düğmesine tıklayarak ortamı özelleştirmek mümkündür.

??? info "Yapılandırma seçenekleri"

    "Change options" düğmesine tıklarsanız, aşağıdakileri özelleştirme seçeneği sunulacaktır:

    #### Branch

    Bu, eğitim materyallerinin farklı bir sürümünü seçmenize olanak tanır.
    `master` dalı genellikle hata düzeltmelerini ve yakın zamanda geliştirilip onaylanmış ancak henüz web sitesinde yayınlanmamış materyalleri içerir.
    Diğer dallar, tam olarak işlevsel olmayabilecek devam eden çalışmaları içerir.

    #### Machine type

    Bu, eğitim boyunca kullanacağınız sanal makineyi özelleştirmenize olanak tanır.

    Daha fazla çekirdeğe sahip bir makine kullanmak, Nextflow'un iş akışı yürütmesini paralelleştirme yeteneğinden daha fazla yararlanmanızı sağlar.
    Ancak, ücretsiz kota tahsisinizi daha hızlı tüketecektir, bu nedenle almayı planladığınız kursun talimatlarında tavsiye edilmedikçe bu ayarı değiştirmenizi önermiyoruz.

    Kotalar hakkında daha fazla ayrıntı için aşağıdaki 'GitHub Codespaces kotaları' bölümüne bakın.

### Başlatma süresi

Yeni bir GitHub Codespaces ortamını ilk kez açmak birkaç dakika sürebilir, çünkü sistemin sanal makinenizi kurması gerekir, bu nedenle bir bekleme süresi olması endişe verici değildir.
Ancak, beş dakikadan fazla sürmemelidir.

## Eğitim arayüzünde gezinme

GitHub Codespaces'iniz yüklendikten sonra, aşağıdakine benzer bir şey görmelisiniz (hesap tercihlerinize bağlı olarak açık modda açılabilir):

![GitHub Codespaces welcome](img/codespaces_welcome.png)

Bu, Nextflow geliştirme için kullanmanızı önerdiğimiz popüler bir kod geliştirme uygulaması olan VSCode IDE'nin arayüzüdür.

- **Ana editör**, Nextflow kodunun ve diğer metin dosyalarının açılacağı yerdir. Kodu burada düzenleyeceksiniz. Codespace'i açtığınızda, bu size `README.md` dosyasının bir önizlemesini gösterecektir.
- Ana editörün altındaki **terminal**, komutları çalıştırmanıza olanak tanır. Kurs talimatlarında verilen tüm komut satırlarını burada çalıştıracaksınız.
- **Kenar çubuğu**, ortamınızı özelleştirmenize ve temel görevleri gerçekleştirmenize (kopyalama, yapıştırma, dosya açma, arama, git vb.) olanak tanır. Varsayılan olarak, deponun içeriğine göz atmanıza olanak tanıyan dosya gezginine açıktır. Gezgindeki bir dosyaya tıklamak, onu ana editör penceresinde açacaktır.

Pencere bölmelerinin göreli oranlarını istediğiniz gibi ayarlayabilirsiniz.

<!-- TODO (future) Link to development best practices side quest? -->

## GitHub Codespaces kullanımı hakkında diğer notlar

### Bir oturuma devam etme

Bir ortam oluşturduktan sonra, kolayca devam edebilir veya yeniden başlatabilir ve kaldığınız yerden devam edebilirsiniz.
Ortamınız 30 dakikalık hareketsizlikten sonra zaman aşımına uğrayacak ve değişikliklerinizi 2 haftaya kadar saklayacaktır.

Bir ortamı <https://github.com/codespaces/> adresinden yeniden açabilirsiniz.
Önceki ortamlar listelenecektir.
Devam etmek için bir oturuma tıklayın.

![List GitHub Codespace sessions](img/codespaces_list.png)

Önceki GitHub Codespaces ortamınızın URL'sini kaydettiyseniz, bunu tarayıcınızda açmanız yeterlidir.
Alternatif olarak, onu ilk etapta oluşturmak için kullandığınız düğmeye tıklayın:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Önceki oturumu görmelisiniz, varsayılan seçenek ona devam etmektir:

![Resume a GitHub Codespace](img/codespaces_resume.png)

### Dosyaları yerel makinenize kaydetme

Gezgin panelinden herhangi bir dosyayı kaydetmek için dosyaya sağ tıklayın ve `Download` seçeneğini seçin.

### GitHub Codespaces kotalarını yönetme

GitHub Codespaces size ayda 15 GB-aylık depolama alanı ve ayda 120 çekirdek-saat verir.
Bu, standart çalışma alanını (2 çekirdek, 8 GB RAM ve 32 GB depolama) kullanarak varsayılan ortam çalışma süresinin yaklaşık 60 saatine eşdeğerdir.

Bunları daha fazla kaynak ile oluşturabilirsiniz (yukarıdaki açıklamaya bakın), ancak bu, ücretsiz kullanımınızı daha hızlı tüketecek ve bu alana daha az saatlik erişiminiz olacaktır.
Örneğin, 2 çekirdekli varsayılan yerine 4 çekirdekli bir makine seçerseniz, kotanız yarı sürede tükenecektir.

İsteğe bağlı olarak, daha fazla kaynağa erişim satın alabilirsiniz.

Daha fazla bilgi için GitHub belgelerine bakın:
[About billing for GitHub Codespaces](https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)
