# Bölüm 1: Temel işlemleri çalıştırın

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow for Bioimaging eğitim kursunun bu ilk bölümünde, temel işlemleri göstermek ve ilgili Nextflow kod bileşenlerine işaret etmek için çok basit, alan-bağımsız bir Hello World örneği kullanacağız.

## 1. Workflow'u çalıştırın

Size `hello-world.nf` adında, `--greeting` adlı bir komut satırı argümanı aracılığıyla girdi alan ve bu selamlamayı içeren bir metin dosyası üreten bir workflow betiği sağlıyoruz.
Henüz koda bakmayacağız; önce çalıştırmanın nasıl göründüğünü görelim.

### 1.1. Workflow'u başlatın ve çalışmayı izleyin

Terminalde aşağıdaki komutu çalıştırın:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

Konsol çıktınız şuna benzer görünmelidir:

```console title="Çıktı" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Tebrikler, ilk Nextflow workflow'unuzu çalıştırdınız!

Buradaki en önemli çıktı son satırdır (satır 6):

```console title="Çıktı" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Bu bize `sayHello` sürecinin bir kez başarıyla çalıştırıldığını söylüyor (`1 of 1 ✔`).

Bu harika, ancak merak ediyor olabilirsiniz: çıktı nerede?

### 1.2. `results` dizininde çıktı dosyasını bulun

Bu workflow, çıktısını `results` adlı bir dizine yayınlamak üzere yapılandırılmıştır.
Mevcut dizininize bakarsanız, workflow'u çalıştırdığınızda Nextflow'un `results` adında yeni bir dizin oluşturduğunu göreceksiniz; bu dizin `output.txt` adlı bir dosya içerir.

```console title="results/" linenums="1"
results
└── output.txt
```

Dosyayı açın; içeriği komut satırında belirttiğiniz selamlama ile eşleşmelidir.

<details>
  <summary>Dosya içeriği</summary>

```console title="results/output.txt" linenums="1"
Hello World!
```

</details>

Harika, workflow'umuz yapması gerekeni yaptı!

Ancak, 'yayınlanmış' sonucun, Nextflow'un workflow'u çalıştırırken ürettiği gerçek çıktının bir kopyası (veya bazı durumlarda bir sembolik bağlantısı) olduğunu unutmayın.

Şimdi, Nextflow'un çalışmayı gerçekte nerede yürüttüğünü görmek için kaputun altına bakacağız.

!!! warning "Uyarı"

    Tüm workflow'lar çıktıları bir results dizinine yayınlayacak şekilde ayarlanmamış olabilir ve/veya dizin adı farklı olabilir.
    Bu bölümde biraz ilerledikten sonra, bu davranışın nerede belirtildiğini nasıl öğreneceğinizi göstereceğiz.

### 1.3. `work/` dizininde orijinal çıktıyı ve günlükleri bulun

Bir workflow çalıştırdığınızda, Nextflow workflow'daki her sürecin her çağrısı için (=pipeline'daki her adım) ayrı bir 'görev dizini' oluşturur.
Her biri için gerekli girdileri hazırlar, ilgili komut(lar)ı çalıştırır ve o dizin içinde çıktıları ve günlük dosyalarını yazar; dizin benzersiz hale getirmek için otomatik olarak bir hash kullanılarak adlandırılır.

Tüm bu görev dizinleri, mevcut dizininizin içinde (komutu çalıştırdığınız yerde) `work` adlı bir dizin altında bulunacaktır.

Bu kafa karıştırıcı gelebilir, o yüzden pratikte bunun nasıl göründüğünü görelim.

Daha önce çalıştırdığımız workflow'un konsol çıktısına geri dönersek, şu satırımız vardı:

```console title="Komut çıktısından alıntı" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Satırın `[a3/7be2fa]` ile nasıl başladığını görüyor musunuz?
Bu, o süreç çağrısı için görev dizini yolunun kısaltılmış bir biçimidir ve `sayHello` süreç çağrısının çıktısını `work/` dizin yolu içinde nerede bulacağınızı söyler.

Tam yolu aşağıdaki komutu yazarak (kendi terminalinizde gördüğünüzle `a3/7be2fa` yerine koyarak) ve tab tuşuna basarak yolu otomatik tamamlayarak veya bir yıldız işareti ekleyerek bulabilirsiniz:

```bash
tree work/a3/7be2fa*
```

Bu, tam dizin yolunu vermelidir: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Orada ne olduğuna bir göz atalım.

!!! Tip "İpucu"

    VSCode dosya gezgininde görev alt dizininin içeriğine göz atarsanız, tüm dosyaları hemen göreceksiniz.
    Ancak, günlük dosyaları terminalde görünmez olarak ayarlanmıştır, bu nedenle bunları görüntülemek için `ls` veya `tree` kullanmak istiyorsanız, görünmez dosyaları görüntülemek için ilgili seçeneği ayarlamanız gerekecektir.

    ```bash
    tree -a work
    ```

Tam alt dizin adları sisteminizde farklı olacaktır.

<details>
  <summary>Dizin içeriği</summary>

```console title="work/"
work
└── a3
    └── 7be2fad5e71e5f49998f795677fd68
        ├── .command.begin
        ├── .command.err
        ├── .command.log
        ├── .command.out
        ├── .command.run
        ├── .command.sh
        ├── .exitcode
        └── output.txt
```

</details>

`output.txt` dosyasını hemen tanımalısınız; aslında bu, `results` dizinine yayınlanan `sayHello` sürecinin orijinal çıktısıdır.
Açarsanız, yine `Hello World!` selamlamasını bulacaksınız.

<details>
  <summary>output.txt dosyasının içeriği</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

Peki ya diğer tüm dosyalar?

Bunlar, Nextflow'un görev yürütmenin bir parçası olarak yazdığı yardımcı ve günlük dosyalarıdır:

- **`.command.begin`**: Görev başlatılır başlatılmaz oluşturulan gözetleyici dosya.
- **`.command.err`**: Süreç çağrısı tarafından yayılan hata mesajları (`stderr`)
- **`.command.log`**: Süreç çağrısı tarafından yayılan tam günlük çıktısı
- **`.command.out`**: Süreç çağrısı tarafından düzenli çıktı (`stdout`)
- **`.command.run`**: Süreç çağrısını çalıştırmak için Nextflow tarafından çalıştırılan tam betik
- **`.command.sh`**: Süreç çağrısı tarafından gerçekten çalıştırılan komut
- **`.exitcode`**: Komuttan kaynaklanan çıkış kodu

`.command.sh` dosyası özellikle yararlıdır çünkü Nextflow'un tüm kayıt tutma ve görev/ortam kurulumu dahil edilmeden çalıştırdığı ana komutu size gösterir.

<details>
  <summary>Dosya içeriği</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip "İpucu"

    Bir şeyler ters gittiğinde ve ne olduğunu sorun gidermeniz gerektiğinde, Nextflow'un workflow talimatlarına, değişken enterpolasyonuna vb. dayalı olarak tam olarak hangi komutu oluşturduğunu kontrol etmek için `command.sh` betiğine bakmak yararlı olabilir.

### 1.4. İsteğe bağlı alıştırma: farklı selamlamalarla yeniden çalıştırın

`--greeting` argümanı için farklı değerlerle workflow'u birkaç kez yeniden çalıştırmayı deneyin, ardından hem `results/` dizininin hem de görev dizinlerinin içeriğine bakın.

İzole görev dizinlerinin çıktılarının ve günlüklerinin nasıl korunduğunu, `results` dizininin içeriğinin ise sonraki çalıştırmaların çıktısıyla nasıl üzerine yazıldığını gözlemleyin.

### Çıkarım

Basit bir Nextflow betiğini nasıl çalıştıracağınızı, çalışmasını nasıl izleyeceğinizi ve çıktılarını nasıl bulacağınızı biliyorsunuz.

### Sırada ne var?

Temel bir Nextflow betiğini nasıl okuyacağınızı ve bileşenlerinin işlevselliğiyle nasıl ilişkili olduğunu nasıl belirleyeceğinizi öğrenin.

---

## 2. Hello World workflow başlangıç betiğini inceleyin

Orada yaptığımız temelde workflow betiğini bir kara kutu gibi ele almaktı.
Artık ne yaptığını gördüğümüze göre, kutuyu açalım ve içine bakalım.

_Buradaki amaç Nextflow kodunun sözdizimini ezberlemek değil, ana bileşenlerin neler olduğu ve nasıl organize edildikleri konusunda temel bir sezgi oluşturmaktır._

### 2.1. Genel kod yapısını inceleyin

`hello-world.nf` betiğini düzenleyici bölmesinde açalım.

<details>
  <summary>Kod</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Bir selamlamayı bir dosyaya yazdırmak için echo kullan
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

workflow {

    // bir selamlama yayınla
    sayHello(params.greeting)
}
```

</details>

Bir Nextflow betiği iki ana temel bileşen türü içerir: bir veya daha fazla **process** ve **workflow**'un kendisi.
Her **process**, pipeline'daki ilgili adımın hangi işlem(ler)i gerçekleştirmesi gerektiğini tanımlarken, **workflow** çeşitli adımları birbirine bağlayan veri akışı mantığını tanımlar.

Önce **process** bloğuna daha yakından bakalım, ardından **workflow** bloğuna bakacağız.

### 2.2. `process` tanımı

Kodun ilk bloğu bir **process** tanımlar.
Process tanımı `process` anahtar kelimesiyle başlar, ardından süreç adı ve son olarak süslü parantezlerle sınırlandırılmış süreç gövdesi gelir.
Süreç gövdesi, çalıştırılacak komutu belirten bir script bloğu içermelidir; bu, komut satırı terminalinde çalıştırabileceğiniz herhangi bir şey olabilir.

Burada `greeting` adlı bir **input** değişkeni alan ve **output**'unu `output.txt` adlı bir dosyaya yazan `sayHello` adında bir **process**'imiz var.

<details>
  <summary>Kod</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Bir selamlamayı bir dosyaya yazdırmak için echo kullan
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}
```

</details>

Bu, sadece bir `input` tanımı, bir `output` tanımı ve çalıştırılacak `script` içeren çok minimal bir süreç tanımıdır.

`input` tanımı, Nextflow'a bir tür değer (bir string, bir sayı, her ne olursa olsun) beklemesini söyleyen `val` niteleyicisini içerir.

`output` tanımı, Nextflow'a bunun bir yol olarak ele alınması gerektiğini söyleyen `path` niteleyicisini içerir (hem dizin yollarını hem de dosyaları içerir).

!!! Tip "İpucu"

    Çıktı tanımı hangi çıktının oluşturulacağını _belirlemez_.
    Sadece beklenen çıktı dosya(lar)ının nerede bulunacağını _bildirir_, böylece Nextflow çalıştırma tamamlandığında bunu arayabilir.

    Bu, komutun başarıyla çalıştırıldığını doğrulamak ve gerekirse çıktıyı aşağı akış süreçlerine aktarmak için gereklidir.
    Çıktı bloğunda bildirilenle eşleşmeyen üretilen çıktı, aşağı akış süreçlerine aktarılmayacaktır.

Gerçek dünyadaki bir pipeline'da, bir süreç genellikle birazdan tanıtacağımız süreç yönergeleri gibi ek bilgiler içerir.

### 2.3. `workflow` tanımı

İkinci kod bloğu **workflow**'un kendisini tanımlar.
Workflow tanımı `workflow` anahtar kelimesiyle başlar, ardından isteğe bağlı bir ad ve süslü parantezlerle sınırlandırılmış workflow gövdesi gelir.

Burada, `--greeting` parametresine verdiğimiz değeri tutan `params.greeting` girdisini alan `sayHello` sürecine bir çağrıdan oluşan bir **workflow**'umuz var.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // bir selamlama yayınla
    sayHello(params.greeting)
}
```

Bu çok minimal bir **workflow** tanımıdır.
Gerçek dünyadaki bir pipeline'da, workflow tipik olarak **channel**'lar tarafından birbirine bağlanan **process**'lere birden fazla çağrı içerir ve değişken girdiler için varsayılan değerler ayarlanmış olabilir.

Bunu kursun 2. Bölümünde nf-core/molkart'ı çalıştırdığımızda uygulamada göreceğiz.

### 2.4. Komut satırı parametrelerinin `params` sistemi

`sayHello()` süreç çağrısına sağladığımız `params.greeting`, Nextflow kodunun düzgün bir parçasıdır ve üzerine fazladan bir dakika harcamaya değer.

Yukarıda belirtildiği gibi, `--greeting` komut satırı parametresinin değerini `sayHello()` süreç çağrısına bu şekilde aktarıyoruz.
Aslında, sadece `params.someParameterName` bildirmek, workflow'a komut satırından `--someParameterName` adında bir parametre vermemizi sağlayacaktır.

!!! Tip "İpucu"

    `params` sistemi kullanılarak bildirilen bu workflow parametreleri her zaman iki tire (`--`) alır.
    Bu, onları sadece bir tire (`-`) alan Nextflow seviyesi parametrelerden ayırır.

### Çıkarım

Artık basit bir Nextflow workflow'unun nasıl yapılandırıldığını ve temel bileşenlerin işlevselliği ile nasıl ilişkili olduğunu biliyorsunuz.

### Sırada ne var?

Workflow çalıştırmalarınızı rahatça yönetmeyi öğrenin.

---

## 3. Workflow çalıştırmalarını yönetin

Workflow'ları nasıl başlatacağınızı ve çıktıları nasıl alacağınızı bilmek harika, ancak yaşamınızı kolaylaştıracak workflow yönetiminin birkaç başka yönü olduğunu hızla fark edeceksiniz.

Burada size aynı workflow'u yeniden başlatmanız gerektiğinde `resume` özelliğinden nasıl yararlanacağınızı, çalıştırma günlüklerini `nextflow log` ile nasıl inceleyeceğinizi ve eski çalışma dizinlerini `nextflow clean` ile nasıl sileceğinizi gösteriyoruz.

### 3.1. `-resume` ile bir workflow'u yeniden başlatın

Bazen, daha önce başlattığınız bir pipeline'ı, zaten başarıyla tamamlanmış herhangi bir işi yeniden yapmadan tekrar çalıştırmak isteyeceksiniz.

Nextflow'un bunu yapmanıza izin veren `-resume` adlı bir seçeneği vardır.
Özellikle, bu modda, tamamen aynı kod, ayarlar ve girdilerle zaten çalıştırılmış olan tüm süreçler atlanacaktır.
Bu, Nextflow'un yalnızca son çalıştırmadan bu yana eklediğiniz veya değiştirdiğiniz süreçleri veya yeni ayarlar veya girdiler sağladığınız süreçleri çalıştıracağı anlamına gelir.

Bunu yapmanın iki temel avantajı vardır:

- Bir pipeline geliştiriyorsanız, değişikliklerinizi test etmek için yalnızca aktif olarak üzerinde çalıştığınız süreç(ler)i çalıştırmanız gerekeceğinden daha hızlı iterasyon yapabilirsiniz.
- Üretimde bir pipeline çalıştırıyorsanız ve bir şeyler ters giderse, çoğu durumda sorunu düzeltebilir ve pipeline'ı yeniden başlatabilirsiniz ve başarısızlık noktasından çalışmaya devam edecektir, bu size çok fazla zaman ve hesaplama tasarrufu sağlayabilir.

Kullanmak için komutunuza `-resume` ekleyin ve çalıştırın:

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

Süreç durum satırına eklenmiş olan `cached:` kısmına bakın (satır 5), bu Nextflow'un bu çalışmayı zaten yaptığını tanıdığı ve önceki başarılı çalıştırmadan sonucu yeniden kullandığı anlamına gelir.

Ayrıca çalışma alt dizini hash'inin önceki çalıştırmayla aynı olduğunu görebilirsiniz.
Nextflow tam anlamıyla size önceki çalıştırmayı gösteriyor ve "Bunu zaten orada yaptım" diyor.

!!! Tip "İpucu"

    Pipeline'ınızı `resume` ile yeniden çalıştırdığınızda, Nextflow daha önce başarıyla çalıştırılan herhangi bir süreç çağrısı tarafından bir `publishDir` dizinine yazılan dosyaların üzerine yazmaz.

### 3.2. Geçmiş çalıştırmaların günlüğünü inceleyin

Bir nextflow workflow'u başlattığınızda, mevcut çalışma dizinindeki `.nextflow` adlı gizli bir dizin altında `history` adlı bir günlük dosyasına bir satır yazılır.

Bu bilgilere erişmenin daha uygun bir yolu `nextflow log` komutunu kullanmaktır.

```bash
nextflow log
```

Bu, günlük dosyasının içeriğini terminale çıktı olarak verir ve mevcut çalışma dizini içinden başlatılmış her Nextflow çalıştırması için zaman damgasını, çalıştırma adını, durumunu ve tam komut satırını gösterir.

### 3.3. Eski çalışma dizinlerini silin

Geliştirme sürecinde, tipik olarak taslak pipeline'larınızı çok sayıda kez çalıştırırsınız, bu da birçok alt dizinde çok sayıda dosyanın birikmesine yol açabilir.
Alt dizinler rastgele adlandırıldığından, hangilerinin daha eski ve hangilerinin daha yeni çalıştırmalar olduğunu adlarından söylemek zordur.

Nextflow, artık umursamadığınız geçmiş çalıştırmaların çalışma alt dizinlerini otomatik olarak silebilen, neyin silineceğini kontrol etmek için çeşitli [seçenekler](https://www.nextflow.io/docs/latest/reference/cli.html#clean) içeren kullanışlı bir `clean` alt komutu içerir.

Zaman damgasına ve/veya komut satırına göre bir çalıştırmayı aramak için Nextflow günlüğünü kullanabilir, ardından önceki çalıştırmalardaki çalışma dizinlerini silmek için `nextflow clean -before <run_name> -f` kullanabilirsiniz.

!!! Warning "Uyarı"

    Geçmiş çalıştırmalardan çalışma alt dizinlerini silmek onları Nextflow'un önbelleğinden kaldırır ve bu dizinlerde depolanan tüm çıktıları siler.
    Bu, Nextflow'un ilgili süreçleri yeniden çalıştırmadan çalıştırmayı sürdürme yeteneğini bozduğu anlamına gelir.

    Önem verdiğiniz veya güvenmeyi planladığınız çıktıları kaydetmek sizin sorumluluğunuzdadır! Bu amaçla `publishDir` yönergesini kullanıyorsanız, `symlink` modunu değil `copy` modunu kullandığınızdan emin olun.

### Çıkarım

Bir pipeline'ı zaten aynı şekilde çalıştırılmış adımları tekrarlamadan nasıl yeniden başlatacağınızı, çalıştırma günlüğünü nasıl inceleyeceğinizi ve eski çalışma dizinlerini temizlemek için `nextflow clean` komutunu nasıl kullanacağınızı biliyorsunuz.

### Sırada ne var?

Artık temel Nextflow işlemlerini anladığınıza göre, nf-core/molkart ile gerçek bir bioimaging pipeline'ı çalıştırmaya hazırsınız.
