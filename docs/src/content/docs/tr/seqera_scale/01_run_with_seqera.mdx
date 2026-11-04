# Bölüm 1: Pipeline'ları web arayüzünden başlatma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Scale with Seqera eğitim kursunun bu bölümünde, Seqera Platform'a erişim kurulumunu yapacak ve web arayüzünden üretim ölçeğinde bir pipeline başlatacaksınız.

Çalışma dizininizin [Başlarken](./00_orientation.md) sayfasında belirtildiği gibi `seqera-scale/` olarak ayarlandığından emin olun.

---

## 1. Seqera ile başlarken

Seqera, Nextflow pipeline'larını başlatmak, izlemek ve yönetmek için kapsamlı bir platform sunar.
Bu bölüm, ilk pipeline'ınızı çalıştırmadan önce hesap oluşturma ve platforma alışma süreçlerinde size rehberlik eder.

### 1.1. Ücretsiz hesap oluşturma

[cloud.seqera.io](https://cloud.seqera.io) adresine gidin ve ücretsiz bir hesap oluşturun.
E-posta adresinizi, GitHub veya Google kimlik bilgilerinizi kullanarak kayıt olabilirsiniz.

Ücretsiz hesapla şunlara erişebilirsiniz:

- **Kişisel çalışma alanı**: pipeline ekleyebileceğiniz, hesaplama ortamlarını yapılandırabileceğiniz ve çalıştırmaları yönetebileceğiniz kendi alanınız
- **Community Showcase'e erişim**: önceden yapılandırılmış ayarlar ve örnek çalıştırma verileriyle birlikte sunulan, nf-core ve topluluk pipeline'larından oluşan özenle seçilmiş bir koleksiyon

Hesap kademeleri ve mevcut özellikler hakkında kapsamlı bir genel bakış için [Seqera belgelerine](https://docs.seqera.io) bakın.

### 1.2. Community Showcase'i keşfetme

Kendi pipeline'larınızı başlatmadan önce, birkaç dakikanızı Community Showcase'i keşfetmeye ayırın.
Gerçek pipeline'lar ve verilerle platformun nasıl göründüğüne dair gerçekçi bir ön izleme sunar.

1. [cloud.seqera.io](https://cloud.seqera.io) adresinden giriş yapın.
2. Sol kenar çubuğunda **Showcase** seçeneğine tıklayın.
3. Mevcut pipeline'lara göz atın — Use nf-core kursundan tanıdık birkaç nf-core pipeline'ı göreceksiniz.
4. Bir pipeline'a tıklayarak yapılandırmasını ve başlatma ayarlarını inceleyin.
5. Önceki çalıştırmalara ait görev düzeyindeki ayrıntılar ve raporlar dahil olmak üzere örnek çalıştırma geçmişlerini keşfetmek için **Runs** seçeneğine tıklayın.

Bu salt okunur bir görünümdür; ancak herhangi bir şey çalıştırmadan önce arayüzün nasıl çalıştığını anlamanızı sağlar.

### 1.3. Hesaplama ortamına sahip bir çalışma alanına erişme

Pipeline başlatmak için yapılandırılmış bir hesaplama ortamına sahip bir çalışma alanı gerekir.

Seqera, hesaplama sağlamak için iki yöntemi destekler:

- **Kendi altyapınızı bağlama**: AWS, Azure, Google Cloud ve HPC zamanlayıcıları (SLURM, LSF, PBS ve diğerleri).
  Kurulum kılavuzları için [hesaplama ortamları belgelerine](https://docs.seqera.io) bakın.
- **Seqera Compute**: AWS üzerinde önceden sağlanmış hesaplama ortamları sunan, bulut hesabı kurulumu gerektirmeyen ve ücretli olarak sunulan yönetilen bir hizmet.
  Doğrudan çalışma alanı ayarlarınızdan etkinleştirebilirsiniz.

**Grup eğitimi:**
Bir grup eğitim oturumuna katılıyorsanız, hesaplama ortamı önceden yapılandırılmış bir organizasyon ve çalışma alanına eklenmiş olabilirsiniz.
Eğitmeniniz size organizasyon adını, çalışma alanı adını ve ihtiyaç duyduğunuz diğer ayrıntıları verecektir.

**Bağımsız çalışma:**
Bu eğitimi kendi başınıza tamamlıyorsanız, yukarıdaki seçeneklerden birini kullanarak kişisel çalışma alanınızda bir hesaplama ortamı kurmanız gerekecektir.
Seqera Compute'u denemek için ücretsiz krediler [talep üzerine sunulmaktadır](https://seqera.io/platform/compute/).

!!! note "Not"

    Bu kursun geri kalanında, yapılandırılmış bir hesaplama ortamına sahip bir çalışma alanına erişiminizin olduğu varsayılmaktadır.
    Grup eğitim oturumundaysanız, eğitmeniniz hangi çalışma alanını ve hesaplama ortamını kullanacağınızı onaylayacaktır.

### Özetle

Bir Seqera hesabınız var, Community Showcase'i keşfettiniz ve hesaplama ortamına sahip bir çalışma alanına erişebiliyorsunuz.

### Sırada ne var?

Seqera Cloud web arayüzünden üretim ölçeğinde bir RNA-seq pipeline'ı başlatın.

---

## 2. nf-core/rnaseq'i web arayüzünden başlatma

Use nf-core kursunda ele alındığı üzere, nf-core/rnaseq pipeline'ı toplu RNA dizi verisi analizi için topluluk tarafından denetlenmiş bir pipeline'dır.

Bu bölümde, pipeline'ı çalışma alanınıza ekleyecek, bir çalıştırma başlatacak ve yürütmeyi izleyeceksiniz.

### 2.1. Pipeline'ı çalışma alanınıza ekleme

nf-core/rnaseq, Seqera Pipelines hizmeti aracılığıyla birkaç tıklamayla çalışma alanınıza eklenebilen, özenle seçilmiş bir pipeline koleksiyonunun parçasıdır.

_Kendi pipeline'larınızı nasıl ekleyeceğinizi bu kursun ilerleyen bölümlerinde göstereceğiz._

1. Topluluk koleksiyonuna göz atmak için [**Seqera Pipelines**](https://seqera.io/pipelines) adresine gidin.
2. `rnaseq` için arama yapın ve **nf-core/rnaseq** seçeneğini seçin.
3. **Launch Pipeline** seçeneğine tıklayın veya sayfanın alt kısmındaki **Launch Pipeline** bölümüne gidin.
4. Giriş yaptığınızdan emin olun ve **Organizations**, **Workspace** ile **Compute Environment** açılır menülerinden uygun değerleri seçin.
   **Gruplar için ipucu:** Paylaşılan bir çalışma alanı kullanıyorsanız, pipeline adına benzersiz bir tanımlayıcı (kullanıcı adınız gibi) ekleyin.
5. **Add pipeline to your Seqera account** seçeneğine tıklayın.

**Pipeline added: View Pipeline** mesajını içeren bir kutu görünecektir.
Bağlantıya tıkladığınızda, başlatma panelindeki pipeline girişine yönlendirilirsiniz.

Pipeline artık çalışma alanınızın **Launchpad** panelinde listelenmekte ve başlatılmaya hazır durumdadır.

### 2.2. Pipeline'ı başlatma

**Launchpad** panelinde veya pipeline ayrıntıları sayfasında pipeline'ın **Launch** düğmesine tıklayın.
Bu işlem yapılandırma arayüzünü açar.

Pipeline, `test` profiliyle önceden yapılandırılmıştır; bu nedenle girdi verisi, çıktı dizini ve genom referansı otomatik olarak doldurulmuştur.
Şimdilik diğer parametreleri ve gelişmiş ayarları görmezden gelebilirsiniz.

Çalıştırmayı gerçekten başlatmak için mavi **Launch** düğmesine tıklayın.

### 2.3. Yürütmeyi izleme

Başlatmanın ardından, pipeline'ınızın **Runs** paneline yönlendirilirsiniz.

Çalıştırma görünümü şunları gösterir:

- **Status**: çalıştırmanın mevcut durumu (gönderildi, çalışıyor, başarılı, başarısız)
- **Command line**: Platform'un oluşturup gönderdiği tam `nextflow run` komutu
- **Parameters**: bu çalıştırma için kullanılan tüm parametre değerleri
- **Tasks**: durum, süre ve kaynak kullanımıyla birlikte her süreç çağrısının yer aldığı tablo

Yürütme ayrıntılarını incelemek için herhangi bir görev satırına tıklayın:

- Çalıştırılan `.command.sh` betiği
- stdout ve stderr günlükleri
- CPU, bellek ve G/Ç metrikleri

**Reports** sekmesi, çalıştırma tamamlandığında tüm örneklerdeki kalite kontrol metriklerini bir araya getiren bir MultiQC raporu gösterecektir.

Çalıştırmanın tamamlanması biraz zaman alacağından, şimdilik devam edecek ve çıktılara bakmak için daha sonra geri döneceğiz.

### Özetle

Bir Seqera çalışma alanına pipeline eklemeyi, çalıştırmayı yapılandırıp başlatmayı ve yürütmeyi ölçekte izlemeyi öğrendiniz.

### Sırada ne var?

`tw` CLI'ını kullanarak tüm bunları komut satırından nasıl yapacağınızı öğreneceğiniz [Bölüm 2](./02_launch_from_cli.md)'ye geçin.

---

## Özet

Bu bölümde şunları öğrendiniz:

- Seqera hesabı oluşturma ve Community Showcase'i keşfetme
- Özenle seçilmiş katalogdan pipeline ekleme, üretim ölçeğinde çalıştırma başlatma ve yürütmeyi izleme
