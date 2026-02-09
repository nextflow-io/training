# Bölüm 1: Temel işlemleri çalıştırma

Nextflow for Bioimaging eğitim kursunun bu ilk bölümünde, temel işlemleri göstermek ve ilgili Nextflow kod bileşenlerine dikkat çekmek için çok basit, alana özgü olmayan bir Hello World örneği kullanacağız.

## 1. İş akışını çalıştırma

Size `hello-world.nf` adında bir iş akışı betiği sağlıyoruz. Bu betik, `--greeting` adlı bir komut satırı argümanı aracılığıyla bir girdi alır ve bu selamlamayı içeren bir metin dosyası üretir.
Henüz koda bakmayacağız; önce çalıştırmanın nasıl göründüğünü görelim.

### 1.1. İş akışını başlatma ve yürütmeyi izleme

Terminalde aşağıdaki komutu çalıştırın:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

Konsol çıktınız şuna benzer görünmelidir:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Tebrikler, ilk Nextflow iş akışınızı çalıştırdınız!

Buradaki en önemli çıktı son satırdır (satır 6):

```console title="Output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Bu bize `sayHello` sürecinin başarıyla bir kez yürütüldüğünü söyler (`1 of 1 ✔`).

Bu harika, ancak merak ediyor olabilirsiniz: çıktı nerede?

### 1.2. Çıktı dosyasını `results` dizininde bulma

Bu iş akışı, çıktısını `results` adlı bir dizine yayınlayacak şekilde yapılandırılmıştır.
Mevcut dizininize bakarsanız, iş akışını çalıştırdığınızda Nextflow'un `results` adında yeni bir dizin oluşturduğunu ve bu dizinin `output.txt` adlı bir dosya içerdiğini göreceksiniz.

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

Harika, iş akışımız yapması gerekeni yaptı!

Ancak, 'yayınlanan' sonucun, Nextflow'un iş akışını yürüttüğünde ürettiği gerçek çıktının bir kopyası (veya bazı durumlarda bir sembolik bağlantı) olduğunu unutmayın.

Şimdi, Nextflow'un işi gerçekte nerede yürüttüğünü görmek için kaputun altına bakacağız.

!!! warning

    Tüm iş akışları çıktıları bir results dizinine yayınlayacak şekilde ayarlanmayacaktır ve/veya dizin adı farklı olabilir.
    Bu bölümde biraz ileride, bu davranışın nerede belirtildiğini nasıl öğreneceğinizi göstereceğiz.

### 1.3. Orijinal çıktıyı ve günlükleri `work/` dizininde bulma

Bir iş akışı çalıştırdığınızda, Nextflow iş akışındaki her sürecin her çağrısı için (=boru hattındaki her adım) ayrı bir 'görev dizini' oluşturur.
Her biri için gerekli girdileri hazırlar, ilgili komut(lar)ı yürütür ve çıktıları ve günlük dosyalarını o tek dizin içine yazar. Bu dizin, benzersiz olması için otomatik olarak bir hash kullanılarak adlandırılır.

Tüm bu görev dizinleri, mevcut dizininizin (komutu çalıştırdığınız yer) altında `work` adlı bir dizinde bulunacaktır.

Bu kafa karıştırıcı gelebilir, o yüzden pratikte bunun nasıl göründüğünü görelim.

Daha önce çalıştırdığımız iş akışının konsol çıktısına geri dönersek, şu satıra sahiptik:

```console title="Excerpt of command output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Satırın `[a3/7be2fa]` ile başladığını görüyor musunuz?
Bu, o tek süreç çağrısı için görev dizini yolunun kısaltılmış bir biçimidir ve size `sayHello` süreç çağrısının çıktısını `work/` dizin yolu içinde nerede bulacağınızı söyler.

Tam yolu bulmak için aşağıdaki komutu yazabilir (kendi terminalinizde gördüğünüz ile `a3/7be2fa`'yı değiştirerek) ve yolu otomatik tamamlamak için tab tuşuna basabilir veya bir yıldız işareti ekleyebilirsiniz:

```bash
tree work/a3/7be2fa*
```

Bu, tam dizin yolunu vermelidir: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Orada ne olduğuna bir göz atalım.

!!! Tip

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

`output.txt` dosyasını hemen tanımalısınız; bu aslında `results` dizinine yayınlanan `sayHello` sürecinin orijinal çıktısıdır.
Açarsanız, `Hello World!` selamlamasını tekrar bulacaksınız.

<details>
  <summary>output.txt dosya içeriği</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

Peki ya diğer tüm dosyalar?

Bunlar, Nextflow'un görev yürütmesinin bir parçası olarak yazdığı yardımcı ve günlük dosyalarıdır:

- **`.command.begin`**: Görev başlatılır başlatılmaz oluşturulan gözcü dosyası.
- **`.command.err`**: Süreç çağrısı tarafından yayınlanan hata mesajları (`stderr`)
- **`.command.log`**: Süreç çağrısı tarafından yayınlanan tam günlük çıktısı
- **`.command.out`**: Süreç çağrısı tarafından yayınlanan normal çıktı (`stdout`)
- **`.command.run`**: Süreç çağrısını yürütmek için Nextflow tarafından çalıştırılan tam betik
- **`.command.sh`**: Süreç çağrısı tarafından gerçekte çalıştırılan komut
- **`.exitcode`**: Komuttan kaynaklanan çıkış kodu

`.command.sh` dosyası özellikle yararlıdır çünkü Nextflow'un tüm kayıt tutma ve görev/ortam kurulumu dahil edilmeden yürüttüğü ana komutu gösterir.

<details>
  <summary>Dosya içeriği</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip

    Bir şeyler ters gittiğinde ve ne olduğunu araştırmanız gerektiğinde, Nextflow'un iş akışı talimatları, değişken enterpolasyonu vb. temelinde tam olarak hangi komutu oluşturduğunu kontrol etmek için `command.sh` betiğine bakmak yararlı olabilir.

### 1.4. İsteğe bağlı alıştırma: farklı selamlamalarla yeniden çalıştırma

İş akışını `--greeting` argümanı için farklı değerlerle birkaç kez yeniden çalıştırmayı deneyin, ardından hem `results/` dizininin hem de görev dizinlerinin içeriğine bakın.

İzole görev dizinlerinin çıktılarının ve günlüklerinin nasıl korunduğunu, `results` dizininin içeriğinin ise sonraki yürütmelerin çıktısıyla nasıl üzerine yazıldığını gözlemleyin.

### Özet

Basit bir Nextflow betiğini nasıl çalıştıracağınızı, yürütmesini nasıl izleyeceğinizi ve çıktılarını nasıl bulacağınızı biliyorsunuz.

### Sırada ne var?

Temel bir Nextflow betiğini nasıl okuyacağınızı ve bileşenlerinin işlevselliğiyle nasıl ilişkili olduğunu öğrenin.

---

## 2. Hello World iş akışı başlangıç betiğini inceleme

Orada yaptığımız şey temelde iş akışı betiğini bir kara kutu gibi ele almaktı.
Şimdi ne yaptığını gördüğümüze göre, kutuyu açalım ve içine bakalım.

_Buradaki amaç Nextflow kodunun sözdizimini ezberlemek değil, ana bileşenlerin neler olduğu ve nasıl organize edildikleri hakkında bazı temel sezgiler oluşturmaktır._

### 2.1. Genel kod yapısını inceleme

Editör bölmesinde `hello-world.nf` betiğini açalım.

<details>
  <summary>Kod</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Bir dosyaya selamlama yazdırmak için echo kullan
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

    // emit a greeting
    sayHello(params.greeting)
}
```

</details>

Bir Nextflow betiği iki ana temel bileşen türünü içerir: bir veya daha fazla **process** ve **workflow**'un kendisi.
Her **process**, boru hattındaki ilgili adımın hangi işlem(ler)i gerçekleştirmesi gerektiğini tanımlarken, **workflow** çeşitli adımları birbirine bağlayan veri akışı mantığını tanımlar.

Önce **process** bloğuna daha yakından bakalım, sonra **workflow** bloğuna bakacağız.

### 2.2. `process` tanımı

Kodun ilk bloğu bir **process**'i tanımlar.
Süreç tanımı `process` anahtar kelimesiyle başlar, ardından süreç adı ve son olarak süslü parantezlerle sınırlandırılmış süreç gövdesi gelir.
Süreç gövdesi, çalıştırılacak komutu belirten bir script bloğu içermelidir; bu, komut satırı terminalinde çalıştırabileceğiniz herhangi bir şey olabilir.

Burada `greeting` adlı bir **input** değişkeni alan ve **output**'unu `output.txt` adlı bir dosyaya yazan `sayHello` adlı bir **process**'imiz var.

<details>
  <summary>Kod</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Bir dosyaya selamlama yazdırmak için echo kullan
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

Bu, sadece bir `input` tanımı, bir `output` tanımı ve yürütülecek `script`'i içeren çok minimal bir süreç tanımıdır.

`input` tanımı, Nextflow'a bir tür değer (bir dize, bir sayı, her ne olursa olsun) beklemesini söyleyen `val` niteleyicisini içerir.

`output` tanımı, Nextflow'a bunun bir yol olarak ele alınması gerektiğini söyleyen `path` niteleyicisini içerir (hem dizin yollarını hem de dosyaları içerir).

!!! Tip

    Çıktı tanımı hangi çıktının oluşturulacağını _belirlemez_.
    Sadece beklenen çıktı dosya(lar)ının nerede bulunacağını _bildirir_, böylece Nextflow yürütme tamamlandığında onu arayabilir.

    Bu, komutun başarıyla yürütüldüğünü doğrulamak ve gerekirse çıktıyı aşağı akış süreçlerine iletmek için gereklidir.
    Çıktı bloğunda bildirilenle eşleşmeyen üretilen çıktı, aşağı akış süreçlerine iletilmeyecektir.

Gerçek dünya bir boru hattında, bir süreç genellikle birazdan tanıtacağımız süreç yönergeleri gibi ek bilgiler içerir.

### 2.3. `workflow` tanımı

Kodun ikinci bloğu **workflow**'un kendisini tanımlar.
İş akışı tanımı `workflow` anahtar kelimesiyle başlar, ardından isteğe bağlı bir ad ve sonra süslü parantezlerle sınırlandırılmış iş akışı gövdesi gelir.

Burada, `--greeting` parametresine verdiğimiz değeri tutan `params.greeting` girdisini alan `sayHello` sürecine bir çağrıdan oluşan bir **workflow**'umuz var.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // emit a greeting
    sayHello(params.greeting)
}
```

Bu çok minimal bir **workflow** tanımıdır.
Gerçek dünya bir boru hattında, iş akışı tipik olarak **channel**'lar tarafından bağlanan **process**'lere birden fazla çağrı içerir ve değişken girdiler için varsayılan değerler ayarlanmış olabilir.

Bunu kursun 2. Bölümünde nf-core/molkart'ı çalıştırdığımızda uygulamada göreceğiz.

### 2.4. Komut satırı parametrelerinin `params` sistemi

`sayHello()` süreç çağrısına sağladığımız `params.greeting`, düzgün bir Nextflow kodu parçasıdır ve üzerinde ekstra bir dakika harcamaya değer.

Yukarıda belirtildiği gibi, `--greeting` komut satırı parametresinin değerini `sayHello()` süreç çağrısına bu şekilde iletiyoruz.
Aslında, sadece `params.someParameterName` bildirmek, iş akışına komut satırından `--someParameterName` adlı bir parametre vermemizi sağlayacaktır.

!!! Tip

    `params` sistemi kullanılarak bildirilen bu iş akışı parametreleri her zaman iki tire (`--`) alır.
    Bu, onları yalnızca bir tire (`-`) alan Nextflow düzeyindeki parametrelerden ayırır.

### Özet

Artık basit bir Nextflow iş akışının nasıl yapılandırıldığını ve temel bileşenlerin işlevselliğiyle nasıl ilişkili olduğunu biliyorsunuz.

### Sırada ne var?

İş akışı yürütmelerinizi rahatça yönetmeyi öğrenin.

---

## 3. İş akışı yürütmelerini yönetme

İş akışlarını nasıl başlatacağınızı ve çıktıları nasıl alacağınızı bilmek harika, ancak hayatınızı kolaylaştıracak iş akışı yönetiminin birkaç başka yönü olduğunu hızla göreceksiniz.

Burada size aynı iş akışını yeniden başlatmanız gerektiğinde `resume` özelliğinden nasıl yararlanacağınızı, yürütme günlüklerini `nextflow log` ile nasıl inceleyeceğinizi ve eski çalışma dizinlerini `nextflow clean` ile nasıl sileceğinizi gösteriyoruz.

### 3.1. Bir iş akışını `-resume` ile yeniden başlatma

Bazen, daha önce başlattığınız bir boru hattını, zaten başarıyla tamamlanmış herhangi bir işi yeniden yapmadan yeniden çalıştırmak isteyeceksiniz.

Nextflow'un bunu yapmanıza izin veren `-resume` adlı bir seçeneği vardır.
Özellikle, bu modda, tam olarak aynı kod, ayarlar ve girdilerle zaten çalıştırılmış olan tüm süreçler atlanacaktır.
Bu, Nextflow'un yalnızca son çalıştırmadan bu yana eklediğiniz veya değiştirdiğiniz süreçleri veya yeni ayarlar veya girdiler sağladığınız süreçleri çalıştıracağı anlamına gelir.

Bunu yapmanın iki temel avantajı vardır:

- Bir boru hattı geliştiriyorsanız, değişikliklerinizi test etmek için yalnızca aktif olarak üzerinde çalıştığınız süreç(ler)i çalıştırmanız gerektiğinden daha hızlı yineleme yapabilirsiniz.
- Üretimde bir boru hattı çalıştırıyorsanız ve bir şeyler ters giderse, çoğu durumda sorunu düzeltebilir ve boru hattını yeniden başlatabilirsiniz ve hata noktasından çalışmaya devam edecektir, bu da size çok fazla zaman ve hesaplama tasarrufu sağlayabilir.

Kullanmak için, komutunuza `-resume` ekleyin ve çalıştırın:

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

Süreç durum satırına (satır 5) eklenmiş olan `cached:` kısmına bakın; bu, Nextflow'un bu işi zaten yaptığını tanıdığı ve önceki başarılı çalıştırmadan sonucu basitçe yeniden kullandığı anlamına gelir.

Ayrıca çalışma alt dizini hash'inin önceki çalıştırmadakiyle aynı olduğunu görebilirsiniz.
Nextflow kelimenin tam anlamıyla size önceki yürütmeyi işaret ediyor ve "Bunu zaten orada yaptım" diyor.

!!! Tip

    Bir boru hattını `resume` ile yeniden çalıştırdığınızda, Nextflow daha önce başarıyla çalıştırılmış herhangi bir süreç çağrısı tarafından bir `publishDir` dizinine yazılan dosyaların üzerine yazmaz.

### 3.2. Geçmiş yürütmelerin günlüğünü inceleme

Bir nextflow iş akışı başlattığınızda, mevcut çalışma dizinindeki `.nextflow` adlı gizli bir dizin altında `history` adlı bir günlük dosyasına bir satır yazılır.

Bu bilgilere erişmenin daha uygun bir yolu `nextflow log` komutunu kullanmaktır.

```bash
nextflow log
```

Bu, günlük dosyasının içeriğini terminale çıktı olarak verecek ve size mevcut çalışma dizini içinden başlatılmış her Nextflow çalıştırması için zaman damgasını, çalıştırma adını, durumunu ve tam komut satırını gösterecektir.

### 3.3. Eski çalışma dizinlerini silme

Geliştirme süreci boyunca, taslak boru hatlarınızı tipik olarak çok sayıda kez çalıştıracaksınız, bu da birçok alt dizinde çok sayıda dosyanın birikmesine yol açabilir.
Alt dizinler rastgele adlandırıldığından, adlarından hangilerinin daha eski ve hangilerinin daha yeni çalıştırmalar olduğunu söylemek zordur.

Nextflow, artık umursamadığınız geçmiş çalıştırmaların çalışma alt dizinlerini otomatik olarak silebilen, neyin silineceğini kontrol etmek için birkaç [seçenek](https://www.nextflow.io/docs/latest/reference/cli.html#clean) içeren kullanışlı bir `clean` alt komutu içerir.

Zaman damgasına ve/veya komut satırına göre bir çalıştırmayı aramak için Nextflow günlüğünü kullanabilir, ardından daha önceki çalıştırmalardan çalışma dizinlerini silmek için `nextflow clean -before <run_name> -f` kullanabilirsiniz.

!!! Warning

    Geçmiş çalıştırmalardan çalışma alt dizinlerini silmek, bunları Nextflow'un önbelleğinden kaldırır ve bu dizinlerde depolanan tüm çıktıları siler.
    Bu, Nextflow'un ilgili süreçleri yeniden çalıştırmadan yürütmeye devam etme yeteneğini bozar.

    Önem verdiğiniz veya güvenmeyi planladığınız çıktıları kaydetmekten siz sorumlusunuz! Bu amaçla `publishDir` yönergesini kullanıyorsanız, `symlink` modunu değil, `copy` modunu kullandığınızdan emin olun.

### Özet

Zaten aynı şekilde çalıştırılmış adımları tekrarlamadan bir boru hattını nasıl yeniden başlatacağınızı, yürütme günlüğünü nasıl inceleyeceğinizi ve eski çalışma dizinlerini temizlemek için `nextflow clean` komutunu nasıl kullanacağınızı biliyorsunuz.

### Sırada ne var?

Artık temel Nextflow işlemlerini anladığınıza göre, nf-core/molkart ile gerçek bir biyogörüntüleme boru hattı çalıştırmaya hazırsınız.
