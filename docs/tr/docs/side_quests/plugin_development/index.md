---
title: Plugin Geliştirme
hide:
  - toc
---

# Plugin Geliştirme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow'un plugin sistemi, dili özel fonksiyonlar, izleme kancaları, yürütme arka uçları ve daha fazlasıyla genişletmenize olanak tanır.
Plugin'ler, topluluğun Nextflow'un çekirdeğini değiştirmeden yeni özellikler eklemesini sağlar; bu da onları pipeline'lar arasında yeniden kullanılabilir işlevselliği paylaşmak için ideal kılar.

Bu eğitim boyunca mevcut plugin'leri nasıl kullanacağınızı ve isteğe bağlı olarak kendi plugin'inizi nasıl oluşturacağınızı öğreneceksiniz.

## Hedef kitle ve ön koşullar

1. Bölüm, mevcut plugin'lerin kullanımını kapsar ve tüm Nextflow kullanıcıları için geçerlidir.
   2-6. Bölümler, kendi plugin'inizi oluşturmayı ele alır ve Groovy kodu ile derleme araçlarını içerir.
   Önceden Java veya Groovy deneyimi gerekmez.

**Ön koşullar**

- Bir GitHub hesabı VEYA [burada](../../envsetup/02_local) açıklandığı şekilde yerel kurulum.
- [Hello Nextflow](../../hello_nextflow/index.md) kursunu veya eşdeğerini tamamlamış olmak.
- Java 21 veya üzeri (eğitim ortamına dahildir; yalnızca 2-6. Bölümler için gereklidir).

**Çalışma dizini:** `side-quests/plugin_development`

## Öğrenme hedefleri

Bu eğitimin sonunda şunları yapabileceksiniz:

**Plugin kullanımı (1. Bölüm):**

- Mevcut plugin'leri iş akışlarınıza kurmak ve yapılandırmak
- Plugin fonksiyonlarını içe aktarmak ve kullanmak

**Plugin geliştirme (2-6. Bölümler):**

- Nextflow'un yerleşik proje oluşturucusunu kullanarak yeni bir plugin projesi oluşturmak
- İş akışlarından çağrılabilecek özel fonksiyonlar uygulamak
- Plugin'inizi yerel olarak derlemek, test etmek ve kurmak
- Özel günlükleme veya bildirimler için iş akışı olaylarını (örn. görev tamamlama, pipeline başlangıcı/bitişi) izlemek
- Plugin'leri özelleştirilebilir kılmak için yapılandırma seçenekleri eklemek
- Plugin'inizi dağıtmak

## Ders planı

#### 1. Bölüm: Plugin temelleri

Mevcut plugin'leri bir Nextflow iş akışında kullanın ve davranışlarını yapılandırın.

#### 2. Bölüm: Plugin projesi oluşturma

Yeni bir plugin projesi oluşturun ve yapısını inceleyin.

#### 3. Bölüm: Özel fonksiyonlar

Özel fonksiyonlar uygulayın, plugin'inizi derleyin ve bir iş akışında çalıştırın.

#### 4. Bölüm: Test etme

Spock framework'ü kullanarak birim testleri yazın ve çalıştırın.

#### 5. Bölüm: İş akışı izleme

Görev sayacı oluşturmak için görev tamamlama gibi olaylara yanıt verin.

#### 6. Bölüm: Yapılandırma ve Dağıtım

Plugin'inizi özelleştirilebilir kılmak için `nextflow.config` dosyasından ayarları okuyun, ardından nasıl paylaşacağınızı öğrenin.

Kursa başlamaya hazır mısınız?

[Öğrenmeye başlayın :material-arrow-right:](01_plugin_basics.md){ .md-button .md-button--primary }
