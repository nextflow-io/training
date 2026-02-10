# GitHub Codespaces

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

GitHub Codespaces, buluttaki sanal makineler tarafından desteklenen, eğitim için önceden yapılandırılmış bir ortam sağlamamıza olanak tanıyan web tabanlı bir platformdur.
Platform, Github (Microsoft'a ait) tarafından işletilmektedir ve Github hesabı olan herkes için ücretsiz olarak (kullanım kotalarıyla) erişilebilir.

!!! warning "Uyarı"

    Kuruluşlara bağlı hesaplar belirli ek kısıtlamalara tabi olabilir.
    Bu sizin durumunuzsa, bağımsız bir kişisel hesap kullanmanız veya yerel kurulum yapmanız gerekebilir.

## GitHub hesabı oluşturma

[GitHub ana sayfasından](https://github.com/) ücretsiz bir GitHub hesabı oluşturabilirsiniz.

## GitHub Codespace'inizi başlatma

GitHub'a giriş yaptıktan sonra, Nextflow eğitim ortamını açmak için tarayıcınızda bu bağlantıyı açın: <https://codespaces.new/nextflow-io/training?quickstart=1&ref=master>

Alternatif olarak, aşağıda gösterilen düğmeye tıklayabilirsiniz; bu düğme her eğitim kursunda (genellikle Oryantasyon sayfasında) tekrarlanmaktadır.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Yeni bir GitHub Codespace oluşturabileceğiniz bir sayfa ile karşılaşmalısınız:

![Create a GitHub Codespace](img/codespaces_create.png)

### Yapılandırma

Genel kullanım için herhangi bir şey yapılandırmanız gerekmez.
Başlattığınız kursta aksi belirtilmedikçe, devam etmek için ana düğmeye tıklamanız yeterlidir.

Ancak "Change options" düğmesine tıklayarak ortamı özelleştirmek mümkündür.

??? info "Yapılandırma seçenekleri"

    "Change options" düğmesine tıklarsanız, aşağıdakileri özelleştirme seçeneği sunulur:

    #### Branch

    Bu, eğitim materyallerinin farklı bir sürümünü seçmenize olanak tanır.
    `master` dalı genellikle hata düzeltmeleri ve yakın zamanda geliştirilmiş ve onaylanmış ancak henüz web sitesinde yayınlanmamış materyalleri içerir.
    Diğer dallar, tam olarak işlevsel olmayabilecek devam eden çalışmaları içerir.

    #### Machine type

    Bu, eğitim üzerinde çalışmak için kullanacağınız sanal makineyi özelleştirmenize olanak tanır.

    Daha fazla çekirdeğe sahip bir makine kullanmak, Nextflow'un iş akışı yürütmeyi paralelleştirme yeteneğinden daha fazla yararlanmanızı sağlar.
    Ancak ücretsiz kota tahsisatınızı daha hızlı tüketir, bu nedenle almayı planladığınız kursun talimatlarında tavsiye edilmedikçe bu ayarı değiştirmenizi önermiyoruz.

    Kotalar hakkında daha fazla ayrıntı için aşağıdaki 'GitHub Codespaces kotaları' bölümüne bakın.

### Başlatma süresi

İlk kez yeni bir GitHub Codespaces ortamı açmak birkaç dakika sürebilir, çünkü sistemin sanal makinenizi kurması gerekir, bu yüzden bir bekleme süresi olursa endişelenmeyin.
Ancak beş dakikadan fazla sürmemelidir.

## Eğitim arayüzünde gezinme

GitHub Codespaces'iniz yüklendikten sonra, aşağıdakine benzer bir şey görmelisiniz (hesap tercihlerinize bağlı olarak açık modda açılabilir):

![GitHub Codespaces welcome](img/codespaces_welcome.png)

Bu, Nextflow geliştirme için kullanmanızı önerdiğimiz popüler bir kod geliştirme uygulaması olan VSCode IDE'nin arayüzüdür.

- **Ana düzenleyici**, Nextflow kodunun ve diğer metin dosyalarının açılacağı yerdir. Burası kodu düzenleyeceğiniz yerdir. Codespace'i açtığınızda, burada `README.md` dosyasının bir önizlemesi gösterilir.
- Ana düzenleyicinin altındaki **terminal**, komutları çalıştırmanıza olanak tanır. Kurs talimatlarında verilen tüm komut satırlarını burada çalıştıracaksınız.
- **Kenar çubuğu**, ortamınızı özelleştirmenize ve temel görevleri (kopyalama, yapıştırma, dosya açma, arama, git vb.) gerçekleştirmenize olanak tanır. Varsayılan olarak, deponun içeriğini görüntülemenize olanak tanıyan dosya gezginine açıktır. Gezginde bir dosyaya tıklamak, dosyayı ana düzenleyici penceresinde açar.

Pencere bölmelerinin göreli oranlarını istediğiniz gibi ayarlayabilirsiniz.

<!-- TODO (future) Link to development best practices side quest? -->

## GitHub Codespaces kullanımı hakkında diğer notlar

### Bir oturuma devam etme

Bir ortam oluşturduğunuzda, kolayca devam edebilir veya yeniden başlatabilir ve kaldığınız yerden devam edebilirsiniz.
Ortamınız 30 dakikalık hareketsizlikten sonra zaman aşımına uğrar ve değişikliklerinizi 2 haftaya kadar kaydeder.

<https://github.com/codespaces/> adresinden bir ortamı yeniden açabilirsiniz.
Önceki ortamlar listelenecektir.
Devam etmek için bir oturuma tıklayın.

![List GitHub Codespace sessions](img/codespaces_list.png)

Önceki GitHub Codespaces ortamınızın URL'sini kaydettiyseniz, tarayıcınızda açabilirsiniz.
Alternatif olarak, ilk başta oluşturmak için kullandığınız düğmeye tıklayın:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Önceki oturumu görmelisiniz, varsayılan seçenek devam etmektir:

![Resume a GitHub Codespace](img/codespaces_resume.png)

### Dosyaları yerel makinenize kaydetme

Gezgin panelinden herhangi bir dosyayı kaydetmek için, dosyaya sağ tıklayın ve `Download` seçeneğini seçin.

### GitHub Codespaces kotalarını yönetme

GitHub Codespaces, ayda 15 GB-ay depolama ve ayda 120 çekirdek-saat verir.
Bu, standart çalışma alanını (2 çekirdek, 8 GB RAM ve 32 GB depolama) kullanarak yaklaşık 60 saatlik varsayılan ortam çalışma süresine eşdeğerdir.

Daha fazla kaynakla oluşturabilirsiniz (yukarıdaki açıklamaya bakın), ancak bu ücretsiz kullanımınızı daha hızlı tüketir ve bu alana daha az erişim saatiniz olur.
Örneğin, varsayılan 2 çekirdekli yerine 4 çekirdekli bir makine seçerseniz, kotanız yarı sürede tükenir.

İsteğe bağlı olarak, daha fazla kaynağa erişim satın alabilirsiniz.

Daha fazla bilgi için GitHub belgelerine bakın:
[GitHub Codespaces için faturalandırma hakkında](https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)
