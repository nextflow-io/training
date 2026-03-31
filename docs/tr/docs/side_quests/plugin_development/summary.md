# Özet

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Plugin Geliştirme eğitimini tamamladınız.
Bu sayfa, her bölümde neler geliştirdiğinizi özetlemekte, dağıtım konusunu ele almakta ve bundan sonra nereye gideceğinize dair rehberlik sunmaktadır.

---

## Öğrendikleriniz

### Bölüm 1: Plugin'leri kullanmak

Nextflow plugin'lerinin kullanıcı perspektifinden nasıl çalıştığını keşfettiniz.
nf-schema ve nf-co2footprint'i yüklediniz, bunları `nextflow.config` aracılığıyla yapılandırdınız ve plugin'lerin girdileri nasıl doğrulayabileceğini, fonksiyon ekleyebileceğini ve pipeline yaşam döngüsü olaylarına nasıl bağlanabileceğini gördünüz.

### Bölüm 2: Kurulum

Java 21+ ile bir plugin geliştirme ortamı kurdunuz, `nextflow plugin create` komutunu kullanarak yeni bir plugin projesi oluşturdunuz ve Nextflow'un beklediği proje yapısını öğrendiniz: kaynak dosyalar, derleme yapılandırması ve Makefile iş akışı.

### Bölüm 3: Özel fonksiyonlar

Bir `PluginExtensionPoint` sınıfında `@Function` anotasyonlu metotlar oluşturarak ilk genişletme noktanızı uyguladınız.
`reverseGreeting` ve `decorateGreeting` fonksiyonlarını geliştirdiniz, ardından bunları bir pipeline betiğinden içe aktarıp çağırdınız.

### Bölüm 4: Test etme

Groovy test çerçevesini kullanarak özel fonksiyonlarınız için birim testleri yazdınız.
Testleri `make test` ile nasıl çalıştıracağınızı ve plugin'inizi yüklemeden önce doğru çalıştığını nasıl doğrulayacağınızı öğrendiniz.

### Bölüm 5: Observer'lar

Pipeline yaşam döngüsü olaylarına bağlanmak için `TraceObserver` arayüzünü uyguladınız.
`GreetingObserver` (pipeline başlangıcına ve tamamlanmasına tepki veren) ve `TaskCounterObserver` (tamamlanan görevleri sayan) bileşenlerini geliştirdiniz, ardından bunları bir `TraceObserverFactory` aracılığıyla kaydettiniz.

### Bölüm 6: Yapılandırma

Çalışma zamanında değerleri okumak için `session.config.navigate()` kullanarak plugin'inizi `nextflow.config` üzerinden yapılandırılabilir hale getirdiniz.
"Tanınmayan yapılandırma seçeneği" uyarılarını ortadan kaldırmak ve IDE desteğini etkinleştirmek amacıyla plugin'inizin seçeneklerini resmi olarak bildiren bir `@ConfigScope` sınıfı eklediniz.

---

## Dağıtım

Plugin'iniz yerel ortamda çalışır hale geldikten sonra Nextflow plugin kayıt defteri aracılığıyla başkalarıyla paylaşabilirsiniz.

### Sürüm yönetimi

Sürümleriniz için [anlamsal sürümlemeyi](https://semver.org/) takip edin:

| Sürüm değişikliği         | Ne zaman kullanılır                    | Örnek                                              |
| ------------------------- | -------------------------------------- | -------------------------------------------------- |
| **MAJOR** (1.0.0 → 2.0.0) | Geriye dönük uyumsuz değişiklikler     | Bir fonksiyonu kaldırmak, dönüş türlerini değiştirmek |
| **MINOR** (1.0.0 → 1.1.0) | Yeni özellikler, geriye dönük uyumlu   | Yeni bir fonksiyon eklemek                         |
| **PATCH** (1.0.0 → 1.0.1) | Hata düzeltmeleri, geriye dönük uyumlu | Mevcut bir fonksiyondaki hatayı düzeltmek          |

Her sürümden önce `build.gradle` dosyasındaki sürümü güncelleyin:

```groovy title="build.gradle"
version = '1.0.0'  // Anlamsal sürümlemeyi kullanın: MAJOR.MINOR.PATCH
```

### Kayıt defterine yayımlamak

[Nextflow plugin kayıt defteri](https://registry.nextflow.io/), plugin'leri toplulukla paylaşmanın resmi yoludur.

Yayımlama iş akışı:

1. [Kayıt defterinde](https://registry.nextflow.io/) plugin adınızı talep edin (GitHub hesabınızla giriş yapın)
2. API kimlik bilgilerinizi `~/.gradle/gradle.properties` dosyasında yapılandırın
3. Her şeyin çalıştığını doğrulamak için testleri çalıştırın: `make test`
4. `make release` komutuyla yayımlayın

Adım adım talimatlar için [resmi yayımlama belgelerine](https://www.nextflow.io/docs/latest/guides/gradle-plugin.html#publishing-a-plugin) bakın.

Yayımlandıktan sonra kullanıcılar plugin'inizi herhangi bir yerel kurulum gerektirmeden yükleyebilir:

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting@1.0.0'
}
```

Nextflow, ilk kullanımda plugin'i kayıt defterinden otomatik olarak indirir.

---

## Plugin geliştirme kontrol listesi

- [ ] Java 21+ yüklü
- [ ] `nextflow plugin create <name> <org>` ile proje oluşturuldu
- [ ] `@Function` metotlarıyla genişletme sınıfı uygulandı
- [ ] Birim testleri yazıldı ve `make test` ile çalıştırıldı
- [ ] `make install` ile derlendi ve yüklendi
- [ ] İsteğe bağlı olarak iş akışı olayları için `TraceObserver` uygulamaları eklendi
- [ ] İsteğe bağlı olarak plugin yapılandırması için `ConfigScope` eklendi
- [ ] `nextflow.config` dosyasında `plugins { id 'plugin-id' }` ile etkinleştirildi
- [ ] Fonksiyonlar `include { fn } from 'plugin/plugin-id'` ile içe aktarıldı
- [ ] Sürümlendi ve kayıt defterine yayımlandı

---

## Temel kod kalıpları

**Fonksiyon tanımı:**

```groovy
@Function
String myFunction(String input, String optional = 'default') {
    return input.transform()
}
```

**Plugin yapılandırması:**

```groovy
nextflowPlugin {
    provider = 'my-org'
    className = 'my.org.MyPlugin'
    extensionPoints = ['my.org.MyExtension']
}
```

**İş akışlarında kullanım:**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { item -> myFunction(item) }
        .view()
}
```

---

## Genişletme noktası özeti

| Tür                    | Sınıf/Anotasyon  | Amaç                                                        |
| ---------------------- | ---------------- | ----------------------------------------------------------- |
| Fonksiyon              | `@Function`      | İş akışlarından çağrılabilir                                |
| Trace Observer         | `TraceObserver`  | İş akışı yaşam döngüsü olaylarına bağlanır                  |
| Yapılandırma kapsamı   | `@ScopeName`     | nextflow.config dosyasında plugin yapılandırmasını tanımlar |

---

## Bundan sonra ne yapmalısınız?

Plugin geliştirme yolculuğunuza devam etmek için bazı pratik adımlar aşağıda verilmiştir.

**Gerçek bir şey geliştirin.**
Kendi çalışmanızdan bir kullanım senaryosu seçin: ekibinizin tekrar tekrar kullandığı özel bir fonksiyon, pipeline tamamlandığında Slack bildirimi gönderen bir observer ya da kuruluşunuzun pipeline'larında seçenekleri standartlaştıran bir yapılandırma kapsamı.
Gerçek bir sorundan başlamak, anlayışınızı derinleştirmenin en hızlı yoludur.

**nf-hello'yu referans olarak kullanın.**
[nf-hello](https://github.com/nextflow-io/nf-hello) deposu, resmi minimal plugin örneğidir.
Yeni projeler için iyi bir başlangıç noktasıdır; bir şeyin nasıl yapılandırıldığını kontrol etmeniz gerektiğinde de kullanışlı bir referanstır.

**Resmi belgeleri okuyun.**
Nextflow belgeleri, bu eğitimin ötesindeki konuları kapsamaktadır: kanal fabrikaları, operatör aşırı yükleme ve gelişmiş observer kalıpları bunların arasında sayılabilir.
[Plugin geliştirme](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html) kılavuzu en kapsamlı referans kaynağıdır.

**Mevcut plugin'leri inceleyin.**
[Nextflow plugins deposu](https://github.com/nextflow-io/plugins), nf-schema, nf-wave ve nf-tower gibi resmi plugin'lerin kaynak kodunu içermektedir.
Üretim ortamındaki plugin kodunu okumak, giriş düzeyindeki örneklerin ötesine geçen kalıpları ve kuralları öğrenmenin en iyi yollarından biridir.

---

## Ek kaynaklar

**Resmi belgeler:**

- [Plugin'leri kullanmak](https://www.nextflow.io/docs/latest/plugins/plugins.html): plugin'leri yükleme ve yapılandırmaya yönelik kapsamlı kılavuz
- [Plugin geliştirmek](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html): ayrıntılı plugin geliştirme referansı
- [Yapılandırma kapsamları](https://nextflow.io/docs/latest/developer/config-scopes.html): plugin'ler için yapılandırma kapsamları oluşturma

**Plugin keşfi:**

- [Nextflow Plugin Kayıt Defteri](https://registry.nextflow.io/): mevcut plugin'lere göz atın ve keşfedin
- [Plugin kayıt defteri belgeleri](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html): kayıt defteri belgeleri

**Örnekler ve referanslar:**

- [nf-hello](https://github.com/nextflow-io/nf-hello): basit örnek plugin (başlangıç için idealdir)
- [Nextflow plugins deposu](https://github.com/nextflow-io/plugins): referans için resmi plugin koleksiyonu
