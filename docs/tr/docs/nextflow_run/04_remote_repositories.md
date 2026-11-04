# Bölüm 4: Uzak Pipeline'ları Çalıştırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Şimdiye kadar yerel olarak depolanan iş akışı betiklerini çalıştırdınız.
Pratikte, GitHub gibi uzak depolarda yayımlanan pipeline'ları kendiniz indirmeden çalıştırmak isteyeceksiniz.

Nextflow bunu kolaylaştırır: herhangi bir pipeline'ı doğrudan bir Git deposu URL'sinden çalıştırabilirsiniz.

---

## 1. GitHub'dan bir pipeline çalıştırma

Uzak bir pipeline'ı çalıştırmanın temel sözdizimi `nextflow run <repository>` şeklindedir; burada `<repository>`, `nextflow-io/hello` gibi bir GitHub deposu yolu, tam bir URL ya da GitLab, Bitbucket veya başka bir Git barındırma hizmetine giden bir yol olabilir.

### 1.1. Pipeline'ı başlatma

Resmi Nextflow "hello" demo pipeline'ını çalıştırın.
Bu, bu kursta çalıştırdığınızdan farklı ve çok daha basit bir pipeline'dır: bu eğitim boyunca kullanılan "Hello" pipeline'ından önce gelir ve yalnızca birkaç sabit kodlanmış dil için bir selamlama yazdırır; dolayısıyla alışkın olduğunuz CSV girdisini veya ASCII sanatını beklemeyin.

```bash
nextflow run nextflow-io/hello
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

### 1.2. Pipeline'ın önbelleğe alındığı yeri bulma

Uzak bir pipeline'ı ilk kez çalıştırdığınızda Nextflow onu indirir ve yerel olarak önbelleğe alır.
Sonraki çalıştırmalar, siz açıkça bir güncelleme talep etmediğiniz sürece önbelleğe alınmış sürümü kullanır.

Nextflow, varsayılan olarak indirilen pipeline'ları `$NXF_HOME/assets` dizinine kaydeder.
Belirli bir pipeline'ın nereye kaydedildiğini ve hangi revizyonların mevcut olduğunu öğrenmek için doğrudan Nextflow'a sorabilirsiniz:

```bash
nextflow info nextflow-io/hello
```

??? success "Komut çıktısı"

    ```console
     project name: nextflow-io/hello
     repository  : https://github.com/nextflow-io/hello
     local path  : /workspaces/.nextflow/assets/.repos/nextflow-io/hello
     main script : main.nf
     revisions   :
     > master (default)
       mybranch
       testing
       v1.1 [t]
       v1.2 [t]
       v1.3 [t]
    ```

    Nextflow, yerel olarak zaten kontrol ettiğiniz her revizyonu `>` ile işaretler; geri kalanlar mevcut olmakla birlikte henüz çalışma kopyasına alınmamıştır.

Şimdiye kadar indirdiğiniz tüm pipeline'ları `nextflow list` komutuyla da listeleyebilirsiniz:

```bash
nextflow list
```

??? success "Komut çıktısı"

    ```console
    nextflow-io/hello
    ```

[Use nf-core](../nfcore_use/01_run_demo.md#12-retrieve-the-pipeline-code) kursu, indirilen bir pipeline'ın kaynak koduna nasıl göz atılacağı da dahil olmak üzere bu önbellekleme mekanizmasını daha ayrıntılı ele almaktadır.

### Özetle

Bir pipeline'ı kendiniz indirmeden doğrudan bir GitHub deposundan nasıl çalıştıracağınızı ve ardından yerel olarak nerede bulacağınızı öğrendiniz.

### Sırada ne var?

Tekrarlanabilirlik için uzak bir pipeline'ın belirli bir sürümünü nasıl sabitleyeceğinizi öğrenin.

---

## 2. Tekrarlanabilirlik için sürüm belirleme

Nextflow, varsayılan olarak varsayılan dalın en son revizyonunu çalıştırır.
`-r` bayrağını kullanarak belirli bir sürümü (tag), dalı veya commit'i sabitleyebilirsiniz.

### 2.1. Belirli bir revizyonu sabitleme

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello:v1.3 ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sick_carson] revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Nextflow bu revizyonu ilk kez talep ettiğinizde indirir; bu nedenle `Pulling` ve `downloaded from` satırları görünür. Aynı revizyonu daha sonra tekrar talep ettiğinizde ise doğrudan `Launching` aşamasına geçilir.
Tam bir revizyonu sabitleme, tekrarlanabilirlik açısından son derece önemlidir.
Bu sayede siz ve iş arkadaşlarınız, depoda o tarihten bu yana ne değişmiş olursa olsun, tam olarak aynı pipeline kodunu çalıştırdığınızdan emin olursunuz.

### 2.2. Revizyonlar yalnızca ilgili çalıştırma için geçerlidir

`-r` ile bir revizyonu sabitleme yalnızca belirttiğiniz çalıştırmayı etkiler: sonraki sade bir `nextflow run` komutunun hangi sürümü kullanacağını değiştirmez.
Pipeline'ı `-r` olmadan tekrar çalıştırmayı deneyin:

```bash
nextflow run nextflow-io/hello
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [hungry_maxwell] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Önceki çalıştırma açıkça `v1.3`'ü sabitlemiş olsa da bu çalıştırma doğrudan varsayılan dala (`master`) geri döner.
Nextflow, kullandığınız her revizyon için ayrı bir yerel çalışma kopyası tutar; `nextflow info` komutundaki `>` işaretleri de bunu gösterir. Ancak Nextflow, en son hangisini çalıştırdığınızı asla hatırlamaz.
Tekrarlanabilirlik tamamen size aittir: önceki bir çalıştırmada sabitlediğiniz bir revizyonun hâlâ geçerli olduğunu varsaymak yerine, önemli olduğu her durumda `-r` bayrağını açıkça belirtin.
Bir pipeline'ın varsayılan dalını öğrenmek için `nextflow info <pipeline>` komutunu çalıştırabilirsiniz; varsayılan dal `(default)` olarak işaretlenen daldır.

### Özetle

Uzak bir pipeline'ı tekrarlanabilir çalıştırma için belirli bir sürüme, dala veya commit'e nasıl sabitleyeceğinizi ve bu sabitlemenin yalnızca o tek çalıştırma için geçerli olduğunu, sonraki çalıştırmalar için değil, öğrendiniz.

### Sırada ne var?

Nextflow pipeline'larını çalıştırma ve yönetmenin temellerini öğrendiniz.
Buradan nereye gideceğinizi öğrenmek için [Kurs özeti](next_steps.md) sayfasına bakın.

---

## Özet

Bu bölümde şunları öğrendiniz:

- Bir pipeline'ı kendiniz indirmeden doğrudan bir GitHub deposundan çalıştırmak
- Tekrarlanabilir çalıştırma için uzak bir pipeline'ı belirli bir revizyona sabitlemek ve bu sabitlemenin yalnızca o tek çalıştırma için geçerli olduğunu anlamak
