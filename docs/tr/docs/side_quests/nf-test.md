# nf-test ile Test Etme

İş akışınızın her bir parçasının yapması gerekeni yaptığını sistematik olarak test edebilmek, tekrarlanabilirlik ve uzun vadeli bakım için kritik öneme sahiptir ve geliştirme sürecinde büyük bir yardımcı olabilir.

Test etmenin neden bu kadar önemli olduğu hakkında biraz konuşalım. Bir iş akışı geliştiriyorsanız, yapacağınız ilk şeylerden biri geçerli olduğunu bildiğiniz ve sonuç üretmesi gereken bazı test verilerini almaktır. Pipeline'a ilk süreci ekler ve çalışması için girdilerinize bağlarsınız. Ardından, her şeyin çalıştığını kontrol etmek için test verileri üzerinde çalıştırırsınız. Bunun işe yaradığını varsayarsak, bir sonraki sürece geçer ve test verilerini tekrar çalıştırırsınız. Pipeline'dan memnun kalana kadar bu süreci tekrarlarsınız.

Ardından, belki de `--skip_process` gibi basit bir doğru ya da yanlış parametresi eklersiniz. Şimdi pipeline'ı beklendiği gibi çalıştığından emin olmak için her parametre ile bir kez olmak üzere iki kez çalıştırmalısınız. Ama bekleyin, `--skip_process`'in gerçekten süreci atladığını nasıl kontrol ederiz? Çıktıları incelemeli veya log dosyalarını kontrol etmeliyiz! Bu zahmetli ve hataya açık bir durumdur.

İş akışınızı geliştirirken, her iterasyonu manuel olarak test etmek hızla o kadar karmaşık hale gelir ki yavaş ve hataya açık olur. Dahası, bir hata bulursanız hatanın pipeline'ınızda tam olarak nereden kaynaklandığını belirlemek çok zor olacaktır. Test etme işte burada devreye girer.

Test etme, pipeline'ınızın her bir parçasının beklendiği gibi çalıştığını sistematik olarak kontrol etmenizi sağlar. İyi yazılmış testlerin bir geliştirici için faydaları çok büyüktür:

- **Güven**: Testler tüm pipeline'ı kapsadığı için, bir şeyi değiştirmenin başka bir şeyi etkilemediğinden emin olabilirsiniz
- **Güvenilirlik**: Birden fazla geliştirici pipeline üzerinde çalıştığında, diğer geliştiricilerin pipeline'ı ve her bileşeni bozmadığını bilirler.
- **Şeffaflık**: Testler, pipeline'ın nerede başarısız olduğunu gösterir ve sorunu takip etmeyi kolaylaştırır. Ayrıca bir tür dokümantasyon işlevi görürler ve bir sürecin veya iş akışının nasıl çalıştırılacağını gösterirler.
- **Hız**: Testler otomatik olduğu için çok hızlı ve tekrarlı bir şekilde çalıştırılabilir. Yeni hatalar ekleme korkusu daha az olarak hızlı iterasyon yapabilirsiniz.

Yazabileceğimiz birçok farklı test türü vardır:

1. **Modül seviyesi testler**: Bireysel süreçler için
2. **İş akışı seviyesi testler**: Tek bir iş akışı için
3. **Pipeline seviyesi testler**: Pipeline'ın bir bütün olarak testi için
4. **Performans testleri**: Pipeline'ın hız ve verimliliği için
5. **Stres testleri**: Aşırı koşullar altında pipeline'ın performansını değerlendirerek limitlerini belirleme

Bireysel süreçleri test etmek, diğer dillerdeki birim testlere benzerdir. İş akışını veya tüm pipeline'ı test etmek, diğer dillerde entegrasyon testleri olarak adlandırılan, bileşenlerin etkileşimlerini test ettiğimiz teste benzerdir.

[**nf-test**](https://www.nf-test.com/) modül, iş akışı ve pipeline seviyesi test yazmanıza olanak tanıyan bir araçtır. Kısacası, pipeline'ın her bir parçasının _izole bir şekilde_ beklendiği gibi çalışıp çalışmadığını sistematik olarak kontrol etmenizi sağlar.

### Öğrenme hedefleri

Bu yan görevde, pipeline için bir iş akışı seviyesi test ve çağırdığı üç süreç için modül seviyesi testler yazmak üzere nf-test'i kullanmayı öğreneceksiniz.

Bu yan görevin sonunda, aşağıdaki teknikleri etkin bir şekilde kullanabileceksiniz:

- Projenizde nf-test'i başlatma
- Modül seviyesi ve iş akışı seviyesi testler oluşturma
- Yaygın onaylama türlerini ekleme
- Snapshot'ların ne zaman kullanılacağını vs. içerik onaylamalarını anlama
- Tüm bir proje için testleri çalıştırma

Bu beceriler, pipeline projelerinizde kapsamlı bir test stratejisi uygulamanıza yardımcı olarak onların daha sağlam ve sürdürülebilir olmasını sağlayacaktır.

### Ön koşullar

Bu yan görevi üstlenmeden önce:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmalısınız.
- Temel Nextflow kavramlarını ve mekanizmalarını (süreçler, kanallar, operatörler, dosyalarla çalışma, meta veri) kullanma konusunda rahat olmalısınız

---

## 0. Başlayın

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md)'nda açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine gidin

Bu eğitimin dosyalarının bulunduğu dizine geçelim.

```bash
cd side-quests/nf-test
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri gözden geçirin

Ana iş akışı dosyasını ve pipeline'a girdi içeren `greetings.csv` adlı bir CSV dosyası bulacaksınız.

```console title="Dizin içeriği"
.
├── greetings.csv
└── main.nf
```

Dosyaların ayrıntılı açıklaması için [Hello Nextflow'dan ısınma](../hello_nextflow/00_orientation.md) bölümüne bakın.

Test edeceğimiz iş akışı, [Hello Workflow](../hello_nextflow/03_hello_workflow.md)'da oluşturulan Hello iş akışının bir alt kümesidir.

??? example "Hello Nextflow iş akışı ne yapar?"

    [Hello Nextflow](../hello_nextflow/index.md) eğitimini yapmadıysanız, bu basit iş akışının ne yaptığına dair kısa bir genel bakış:

    İş akışı, selamlamalar içeren bir CSV dosyasını alır, bunlar üzerinde ardışık dört dönüştürme adımı çalıştırır ve selamlamaları söyleyen eğlenceli bir karakterin ASCII resmini içeren tek bir metin dosyası çıktılar.

    Dört adım, ayrı modül dosyalarında saklanan Nextflow süreçleri (`sayHello`, `convertToUpper`, `collectGreetings` ve `cowpy`) olarak uygulanır.

    1. **`sayHello`:** Her selamlamayı kendi çıktı dosyasına yazar (örn., "Hello-output.txt")
    2. **`convertToUpper`:** Her selamlamayı büyük harfe dönüştürür (örn., "HELLO")
    3. **`collectGreetings`:** Tüm büyük harf selamlamalarını tek bir toplu dosyada toplar
    4. **`cowpy`:** `cowpy` aracını kullanarak ASCII sanatı oluşturur

    Sonuçlar `results/` adlı bir dizine yayınlanır ve pipeline'ın nihai çıktısı (varsayılan parametrelerle çalıştırıldığında) büyük harfe dönüştürülmüş selamlamaları söyleyen bir karakterin ASCII sanatını içeren düz bir metin dosyasıdır.

    Bu yan görevde, yalnızca ilk iki süreci içeren Hello iş akışının ara bir formunu kullanıyoruz.

Üzerinde çalışacağımız alt küme iki süreçten oluşur: `sayHello` ve `convertToUpper`.
Tam iş akışı kodunu aşağıda görebilirsiniz.

??? example "İş akışı kodu"

    ```groovy title="main.nf"
    /*
    * Pipeline parametreleri
    */
    params.input_file = "greetings.csv"

    /*
    * 'Hello World!' yazdırmak için echo kullan
    */
    process sayHello {

        publishDir 'results', mode: 'copy'

        input:
            val greeting

        output:
            path "${greeting}-output.txt"

        script:
        """
        echo '$greeting' > '$greeting-output.txt'
        """
    }

    /*
    * Selamlamayı büyük harfe dönüştürmek için metin değiştirme aracı kullan
    */
    process convertToUpper {

        publishDir 'results', mode: 'copy'

        input:
            path input_file

        output:
            path "UPPER-${input_file}"

        script:
        """
        cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
        """
    }

    workflow {

        // CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // bir selamlama yayınla
        sayHello(greeting_ch)

        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)
    }
    ```

#### İş akışını çalıştırın

İş akışının beklendiği gibi çalıştığından emin olmak için çalıştıralım.

```bash
nextflow run main.nf
```

```console title="İş akışını çalıştırma sonucu"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

TEBRİKLER! Az önce bir test çalıştırdınız!

"Dur, ne? Sadece iş akışını çalıştırdım ve işe yaradı! Bu nasıl bir test?"

İyi soru!

Az önce ne olduğunu inceleyelim.

İş akışını varsayılan parametrelerle çalıştırdınız, çalıştığını onayladınız ve sonuçlardan memnunsunuz. Bu, test etmenin özüdür. Hello Nextflow eğitim kursunu çalıştıysanız, her bölüme her zaman kullandığımız iş akışını çalıştırarak başladığımızı fark etmişsinizdir, her şeyin doğru şekilde kurulduğunu onaylamak için.

Yazılım testi esasen bu süreci bizim için yapar.

#### Görevi gözden geçirin

Zorluk, herhangi bir değişiklik yapılması durumunda her parçanın beklendiği gibi çalışmaya devam ettiğini doğrulamayı kolaylaştırmak için bu iş akışına nf-test kullanarak standartlaştırılmış testler eklemektir.

#### Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun hedefini ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışır durumda
- [ ] Çalışma dizinini uygun şekilde ayarladım
- [ ] İş akışını başarıyla çalıştırdım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. `nf-test`'i başlatın

`nf-test` paketi, projemiz için test geliştirmeye başlamamız için birkaç şeyi ayarlayan bir başlatma komutu sağlar.

```bash
nf-test init
```

Bu, aşağıdaki çıktıyı üretmelidir:

```bash
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

Ayrıca bir yapılandırma dosyası taslağı içeren bir `tests` dizini oluşturur.

### 1.1. Bir nf-test oluşturun

`nf-test`, nf-test dosyalarını oluşturmak için bir dizi araçla birlikte gelir ve işin çoğunu bizim için yapar. Bunlar `generate` alt komutu altında gelir. Pipeline için bir test oluşturalım:

```bash
nf-test generate pipeline main.nf
```

```console title="Çıktı"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

Bu, `tests` dizini içinde bir `main.nf.test` dosyası oluşturacaktır. Bu bizim pipeline seviyesi test dosyamızdır. `tree tests/` komutunu çalıştırırsanız şöyle bir şey görmelisiniz:

```console title="Test dizini içeriği"
tests/
├── main.nf.test
└── nextflow.config
```

`main.nf.test` dosyası bizim pipeline seviyesi test dosyamızdır. Açıp içeriğine bir göz atalım.

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

Test dosyasının yapısını anlamak için bir saniye ayıralım.

`nextflow_pipeline` bloğu, tüm pipeline seviyesi testler için giriş noktasıdır. Aşağıdakileri içerir:

- `name`: Testin adı.
- `script`: Pipeline betiğinin yolu.

`test` bloğu gerçek testtir. Aşağıdakileri içerir:

- `when`: Testin çalıştırılması gereken koşullar. Bu, pipeline'ı çalıştırmak için kullanılacak parametreleri içerir.
- `then`: Yapılması gereken onaylamalar. Bu, pipeline'ın beklenen sonuçlarını içerir.

Düz Türkçe ile, testin mantığı şöyle okunur:
"**When** bu _parametreler_ bu _pipeline_'a sağlandığında, **then** bu sonuçları görmeyi bekliyoruz."

Bu işlevsel bir test değildir, bir sonraki bölümde bunu nasıl işlevsel bir teste dönüştüreceğimizi göstereceğiz.

### Test Adları Hakkında Bir Not

Yukarıdaki örnekte, yalnızca pipeline'ın başarılı bir şekilde çalışıp çalışmadığını kontrol eden temel bir test için uygun olan varsayılan "Should run without failures" adını kullandık. Ancak, daha spesifik test durumları ekledikçe, gerçekte neyi test ettiğimizi gösteren daha açıklayıcı adlar kullanmalıyız. Örneğin:

- "Should convert input to uppercase" - belirli işlevselliği test ederken
- "Should handle empty input gracefully" - uç durumları test ederken
- "Should respect max memory parameter" - kaynak kısıtlamalarını test ederken
- "Should create expected output files" - dosya oluşturmayı test ederken

İyi test adları şunları içermelidir:

1. Beklenen davranışın ne olduğunu açık hale getirmek için "Should" ile başlamalı
2. Test edilen belirli işlevselliği veya senaryoyu tanımlamalı
3. Test başarısız olursa, hangi işlevselliğin bozuk olduğunu bilecek kadar açık olmalı

Daha sonra daha fazla onaylama ve spesifik test durumu ekledikçe, her testin neyi doğruladığını açık hale getirmek için bu daha açıklayıcı adları kullanacağız.

### 1.2. Testi çalıştırın

Ne olduğunu görmek için testi çalıştıralım.

```bash
nf-test test tests/main.nf.test
```

```console title="nf-test pipeline başarısız"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures' FAILED (4.652s)

  Assertion failed:

  assert workflow.success
         |        |
         workflow false

  Nextflow stdout:

  ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv

   -- Check '/workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/meta/nextflow.log' file for details
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.679s (1 failed)
```

Test başarısız oluyor! Ne oldu?

1. nf-test, `when` bloğundaki ayarları kullanarak pipeline'ı olduğu gibi çalıştırmaya çalıştı:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test, pipeline'ın durumunu kontrol etti ve `when` bloğuyla karşılaştırdı:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

nf-test'in pipeline'ın başarısız olduğunu nasıl bildirdiğine ve Nextflow'dan hata mesajını sağladığına dikkat edin:

```console title="Hata"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

Peki sorun neydi? Pipeline'ın proje dizininde bir greetings.csv dosyası olduğunu hatırlayın. nf-test pipeline'ı çalıştırdığında bu dosyayı arayacaktır, ancak bulamaz. Dosya orada, ne oluyor? Yola bakarsak testin `./nf-test/tests/longHashString/` yolunda gerçekleştiğini görebiliriz. Tıpkı Nextflow gibi, nf-test de her şeyi izole tutmak için her test için yeni bir dizin oluşturur. Veri dosyası orada bulunmadığı için orijinal testteki dosyanın yolunu düzeltmeliyiz.

Test dosyasına geri dönelim ve `when` bloğundaki dosyanın yolunu değiştirelim.

Testte pipeline'ın köküne nasıl işaret edeceğimizi merak ediyor olabilirsiniz. Bu yaygın bir durum olduğu için, nf-test hayatımızı kolaylaştırmak için bir dizi global değişken sağlar. Tam listeyi [burada](https://www.nf-test.com/docs/testcases/global_variables/) bulabilirsiniz, ancak şimdilik pipeline projesinin kökü anlamına gelen `projectDir` değişkenini kullanacağız.

_Önce:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

_Sonra:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
when {
    params {
        input_file = "${projectDir}/greetings.csv"
    }
}
```

İşe yarayıp yaramadığını görmek için testi tekrar çalıştıralım.

```bash title="nf-test pipeline başarılı"
nf-test test tests/main.nf.test
```

```console title="Pipeline başarılı"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

Başarılı! Pipeline başarıyla çalışıyor ve test geçiyor. İstediğiniz kadar çalıştırın ve her zaman aynı sonucu alacaksınız!

Varsayılan olarak, Nextflow çıktısı gizlidir, ancak nf-test'in kesinlikle iş akışını çalıştırdığına kendinizi ikna etmek için `--verbose` bayrağını kullanabilirsiniz:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Pipeline tüm süreçleri çalıştırıyor"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures'
    > Nextflow 24.10.4 is available - Please consider updating your version to it
    > N E X T F L O W  ~  version 24.10.0
    > Launching `/workspaces/training/side-quests/nf-test/main.nf` [zen_ampere] DSL2 - revision: bbf79d5c31
    > [2b/61e453] Submitted process > sayHello (2)
    > [31/4e1606] Submitted process > sayHello (1)
    > [bb/5209ee] Submitted process > sayHello (3)
    > [83/83db6f] Submitted process > convertToUpper (2)
    > [9b/3428b1] Submitted process > convertToUpper (1)
    > [ca/0ba51b] Submitted process > convertToUpper (3)
    PASSED (5.206s)


SUCCESS: Executed 1 tests in 5.239s
```

### 1.3. Onaylamalar ekleyin

Basit bir kontrol, pipeline'ımızın beklediğimiz tüm süreçleri çalıştırdığından ve hiçbirini sessizce atlamadığından emin olmaktır. Pipeline'ımızın 3 selamlamanın her biri için bir `sayHello` ve bir `convertToUpper` olmak üzere 6 süreç çalıştırdığını hatırlayın.

Testımıze pipeline'ın beklenen sayıda süreci çalıştırdığını kontrol etmek için bir onaylama ekleyelim. Ayrıca test adımızı test ettiğimiz şeyi daha iyi yansıtacak şekilde güncelleyelim.

**Önce:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1"
    test("Should run without failures") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
        }

    }
```

**Sonra:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1 11"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

Test adı artık gerçekte neyi doğruladığımızı daha iyi yansıtıyor - sadece pipeline'ın başarısız olmadan çalıştığını değil, aynı zamanda beklenen sayıda süreci çalıştırdığını.

İşe yarayıp yaramadığını görmek için testi tekrar çalıştıralım.

```bash title="nf-test pipeline başarılı"
nf-test test tests/main.nf.test
```

```console title="Pipeline onaylamalarla başarılı"
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

Başarılı! Pipeline başarıyla çalışıyor ve test geçiyor. Artık genel durumun yanı sıra pipeline'ın ayrıntılarını da test etmeye başladık.

### 1.4. Çıktıyı test edin

Testımıze çıktı dosyasının oluşturulduğunu kontrol etmek için bir onaylama ekleyelim. Sonuçları yorumlamayı kolaylaştırmak için bilgilendirici bir adla ayrı bir test olarak ekleyeceğiz.

**Önce:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

**Sonra:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14-33"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }

    test("Should produce correct output files") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert file("$launchDir/results/Bonjour-output.txt").exists()
            assert file("$launchDir/results/Hello-output.txt").exists()
            assert file("$launchDir/results/Holà-output.txt").exists()
            assert file("$launchDir/results/UPPER-Bonjour-output.txt").exists()
            assert file("$launchDir/results/UPPER-Hello-output.txt").exists()
            assert file("$launchDir/results/UPPER-Holà-output.txt").exists()
        }

    }
```

İşe yarayıp yaramadığını görmek için testi tekrar çalıştırın.

```bash title="nf-test pipeline başarılı"
nf-test test tests/main.nf.test
```

```console title="Pipeline dosya onaylamalarıyla başarılı"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

Başarılı! Testler geçer çünkü pipeline başarıyla tamamlandı, doğru sayıda süreç çalıştı ve çıktı dosyaları oluşturuldu. Bu aynı zamanda testleriniz için bilgilendirici adlar sağlamanın ne kadar yararlı olduğunu göstermelidir.

Bu sadece yüzey, pipeline'ın ayrıntılarını kontrol etmek için onaylamalar yazmaya devam edebiliriz, ancak şimdilik pipeline'ın içini test etmeye geçelim.

### Çıkarım

Bir pipeline için nf-test yazmayı biliyorsunuz.

### Sırada ne var?

Bir Nextflow sürecini test etmeyi öğrenin.

---

## 2. Bir Nextflow sürecini test edin

Pipeline'ın her parçası için test yazmak zorunda değiliz, ancak ne kadar çok testımız olursa pipeline hakkında o kadar kapsamlı olabiliriz ve beklendiği gibi çalıştığından o kadar emin olabiliriz. Bu bölümde pipeline'daki her iki süreci de bireysel birimler olarak test edeceğiz.

### 2.1. `sayHello` sürecini test edin

`sayHello` süreci ile başlayalım.

Süreç için testler oluşturmak üzere `nf-test generate` komutunu tekrar kullanalım.

```bash
nf-test generate process main.nf
```

```console title="Çıktı"
> nf-test generate process main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.sayhello.nf.test
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

Şimdilik `main.sayhello.nf.test` dosyasındaki `sayhello` sürecine odaklanalım.

Dosyayı açıp içeriğine bir göz atalım.

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Daha önce olduğu gibi, test ayrıntılarıyla başlıyoruz, ardından `when` ve `then` blokları geliyor. Ancak, sürece girdileri tanımlamamıza izin veren ek bir `process` bloğumuz da var.

İşe yarayıp yaramadığını görmek için testi çalıştıralım.

```bash title="nf-test pipeline başarılı"
nf-test test tests/main.sayhello.nf.test
```

```console title="Süreç testi başarısız"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [1eaad118] 'Should run without failures' FAILED (4.876s)

  Assertion failed:

  assert process.success
         |       |
         |       false
         sayHello

  Nextflow stdout:

  Process `sayHello` declares 1 input but was called with 0 arguments
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.884s (1 failed)
```

Test başarısız oluyor çünkü `sayHello` süreci 1 girdi bildiriyor ancak 0 argümanla çağrıldı. Sürece bir girdi ekleyerek bunu düzeltelim. [Hello Workflow](../hello_nextflow/03_hello_workflow.md)'dan (ve yukarıdaki ısınma bölümünden) `sayHello` sürecimizin sağlamamız gereken tek bir değer girdisi aldığını hatırlayın. Ayrıca test adını test ettiğimiz şeyi daha iyi yansıtacak şekilde düzeltmeliyiz.

**Önce:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Sonra:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

İşe yarayıp yaramadığını görmek için testi tekrar çalıştıralım.

```console title="nf-test pipeline başarılı"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.604s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.611s
```

Başarılı! Test geçer çünkü `sayHello` süreci başarıyla çalıştı ve çıktı oluşturuldu.

### 2.2. Test tarafından oluşturulan snapshot'a göz atın

`tests/main.sayhello.nf.test` dosyasına bakarsak, onaylama bloğunda bir `snapshot()` metodu kullandığını görebiliriz:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

Bu, nf-test'e `sayHello` sürecinin çıktısının bir snapshot'ını oluşturmasını söylüyor. Snapshot dosyasının içeriğine bir göz atalım.

```console title="Snapshot dosyası içeriği"
code tests/main.sayhello.nf.test.snap
```

Burada yazdırmayacağız, ancak süreç ve süreç çıktılarının ayrıntılarını içeren bir JSON dosyası görmelisiniz. Özellikle, şuna benzeyen bir satır görebiliriz:

```json title="Snapshot dosyası içeriği"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

Bu, açıkça test ettiğimiz `sayHello` süreci tarafından oluşturulan çıktıları temsil eder. Testi yeniden çalıştırırsak, program yeni çıktının başlangıçta kaydedilen çıktıyla eşleşip eşleşmediğini kontrol edecektir. Bu, süreç çıktılarının değişmediğini test etmenin hızlı, basit bir yoludur, bu nedenle nf-test bunu varsayılan olarak sağlar.

!!!warning

    Bu, orijinal çalıştırmada kaydettiğimiz çıktının doğru olduğundan emin olmamız gerektiği anlamına gelir!

Gelecekteki geliştirme sırasında, kodda çıktının farklı olmasına neden olan bir şey değişirse, test başarısız olacak ve değişikliğin beklenip beklenmediğini belirlememiz gerekecektir.

- Kodda bir şeyin bozulduğu ortaya çıkarsa, testin geçeceği beklentisiyle düzeltmemiz gerekecektir.
- Beklenen bir değişiklikse (örn., araç iyileştirildi ve sonuçlar daha iyi) o zaman yeni çıktıyı eşleştirilecek referans olarak kabul etmek için snapshot'ı güncellememiz gerekecektir. nf-test bu amaç için bir `--update-snapshot` parametresine sahiptir.

Testi tekrar çalıştırabiliriz ve testin geçmesi gerektiğini görebiliriz:

```console title="snapshot ile nf-test süreci başarılı"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

Başarılı! Test geçer çünkü `sayHello` süreci başarıyla çalıştı ve çıktı snapshot ile eşleşti.

### 2.3. Snapshot'lara Alternatif: Doğrudan İçerik Onaylamaları

Snapshot'lar çıktıdaki herhangi bir değişikliği yakalamak için harika olsa da, bazen tüm dosyanın eşleşmesi konusunda bu kadar katı olmadan belirli içeriği doğrulamak isteyebilirsiniz. Örneğin:

- Çıktının bazı kısımları değişebilir (zaman damgaları, rastgele kimlikler, vb.) ancak belirli anahtar içerik mevcut olmalıdır
- Çıktıda belirli desenler veya değerleri kontrol etmek istediğinizde
- Testi başarıyı neyin oluşturduğu konusunda daha açık hale getirmek istediğinizde

Testımızı belirli içeriği kontrol edecek şekilde nasıl değiştirebileceğimiz aşağıda:

**Önce:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Sonra:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 16 17"
    test("Should run without failures and contain expected greeting") {

        when {
            params {
                // parametreleri burada tanımlayın
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('hello')
            assert !path(process.out[0][0]).readLines().contains('HELLO')
        }

    }
```

nf-test'in süreç çıktılarını liste listeleri olarak gördüğüne dikkat edin, bu nedenle `process.out[0][0]` bu süreçten ilk kanal öğesinin (veya 'emisyonunun') ilk kısmını getiriyor.

Bu yaklaşım:

- Çıktıda tam olarak ne beklediğimizi açık hale getirir
- Çıktıdaki alakasız değişikliklere karşı daha dayanıklıdır
- Testler başarısız olduğunda daha iyi hata mesajları sağlar
- Daha karmaşık doğrulamalara izin verir (regex desenleri, sayısal karşılaştırmalar, vb.)

İşe yarayıp yaramadığını görmek için testi çalıştıralım.

```bash title="nf-test pipeline başarılı"
nf-test test tests/main.sayhello.nf.test
```

```console title="Süreç testi başarısız"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (7.196s)


SUCCESS: Executed 1 tests in 7.208s
```

### 2.4. `convertToUpper` sürecini test edin

`tests/main.converttoupper.nf.test` dosyasını açıp içeriğine bir göz atalım:

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Bu, `sayHello` sürecine benzer bir testtir, ancak `convertToUpper` sürecini test ediyor. Bunun başarısız olacağını biliyoruz çünkü tıpkı `sayHello` gibi, `convertToUpper` süreci tek bir yol girdisi alır, ancak biz bir tane belirtmedik.

Şimdi convertToUpper sürecine büyük harfe dönüştürmek istediğimiz bazı metinler içeren tek bir girdi dosyası sağlamamız gerekiyor. Bunu yapmanın birçok yolu var:

- Test için özel bir dosya oluşturabiliriz
- Mevcut data/greetings.csv dosyasını yeniden kullanabiliriz
- Test içinde anında oluşturabiliriz

Şimdilik, pipeline seviyesi testinde kullandığımız örneği kullanarak mevcut data/greetings.csv dosyasını yeniden kullanalım. Daha önce olduğu gibi, testi test ettiğimiz şeyi daha iyi yansıtacak şekilde adlandırabiliriz, ancak bu sefer içeriği belirli dizeler için kontrol etmek yerine (diğer süreçte yaptığımız gibi) 'snapshot' ile bırakalım.

**Önce:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Sonra:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "${projectDir}/greetings.csv"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Ve testi çalıştırın!

```bash title="nf-test pipeline başarılı"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test süreci convertToUpper başarılı"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.755s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.764s
```

Not: `tests/main.converttoupper.nf.test.snap` adresinde `convertToUpper` süreci için bir snapshot dosyası oluşturduk. Testi tekrar çalıştırırsak, nf-test'in tekrar geçtiğini görmeliyiz.

```bash title="nf-test süreci convertToUpper başarılı"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test süreci convertToUpper başarılı"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s
```

### Çıkarım

Bir Nextflow süreci için test yazmayı ve çalıştırmayı biliyorsunuz.

### Sırada ne var?

Her şey için aynı anda testleri nasıl çalıştıracağınızı öğrenin!

## 3. Tüm depo için testleri çalıştırın

Her bileşen üzerinde nf-test çalıştırmak iyidir, ancak zahmetli ve hataya açıktır. Her şeyi aynı anda test edemez miyiz?

Evet yapabiliriz!

Tüm repo üzerinde nf-test çalıştıralım.

### 3.1. Tüm repo üzerinde nf-test çalıştırın

`nf-test test` komutunu çalıştırarak tüm repo üzerinde nf-test çalıştırabiliriz.

```bash
nf-test test .
```

Not: Her testi içerecek şekilde mevcut dizinimizden her şeyi çalıştırmak için sadece `.` kullanıyoruz!

```console title="nf-test repo başarılı"
> nf-test test .

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [3d26d9af] 'Should run without failures and produce correct output' PASSED (4.155s)

Test Workflow main.nf

  Test [f183df37] 'Should run successfully with correct number of processes' PASSED (3.33s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (3.102s)

Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (2.614s)


SUCCESS: Executed 4 tests in 13.481s
```

Buna bakın! Tek bir komutla her süreç için 1 ve tüm pipeline için 2 olmak üzere 4 test çalıştırdık. Büyük bir kod tabanında bunun ne kadar güçlü olduğunu hayal edin!

---

## Özet

Bu yan görevde, bireysel süreçler için testler oluşturmanın yanı sıra tüm pipeline için uçtan uca testler oluşturmak üzere nf-test'in özelliklerinden yararlanmayı öğrendiniz.
Artık çıktı doğrulama için iki ana yaklaşımın, snapshot'ların ve doğrudan içerik onaylamalarının ve her birinin ne zaman kullanılacağının farkındasınız.
Ayrıca testleri tek tek veya tüm bir proje için nasıl çalıştıracağınızı biliyorsunuz.

Bu teknikleri kendi çalışmanızda uygulamak şunları sağlamanıza olanak tanıyacaktır:

- Kodunuzun beklendiği gibi çalışması
- Değişikliklerin mevcut işlevselliği bozmaması
- Diğer geliştiricilerin güvenle katkıda bulunabilmesi
- Sorunların hızlı bir şekilde tespit edilip düzeltilebilmesi
- Çıktı içeriğinin beklentileri karşılaması

### Anahtar desenler

1. Pipeline seviyesi testler:
   - Temel başarı testi
   - Süreç sayısı doğrulama
   - Çıktı dosyası varlığı kontrolleri
2. Süreç seviyesi testler
3. Çıktı doğrulama için iki yaklaşım:
   - Tam çıktı doğrulama için snapshot'lar kullanma
   - Belirli içerik kontrolleri için doğrudan içerik onaylamaları kullanma
4. Tek bir komutla bir depodaki tüm testleri çalıştırma

### Ek kaynaklar

Daha gelişmiş test özellikleri ve en iyi uygulamalar için [nf-test dokümantasyonuna](https://www.nf-test.com/) göz atın. Şunları yapmak isteyebilirsiniz:

- Testlerinize daha kapsamlı onaylamalar ekleyin
- Uç durumlar ve hata koşulları için testler yazın
- Testleri otomatik olarak çalıştırmak için sürekli entegrasyon kurun
- İş akışı ve modül testleri gibi diğer test türleri hakkında bilgi edinin
- Daha gelişmiş içerik doğrulama tekniklerini keşfedin

**Unutmayın:** Testler, kodunuzun nasıl davranması gerektiğinin canlı dokümantasyonudur. Ne kadar çok test yazarsanız ve onaylamalarınız ne kadar spesifik olursa, pipeline'ınızın güvenilirliğinden o kadar emin olabilirsiniz.

---

## Sırada ne var?

[Yan Görevler menüsüne](./index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ altındaki düğmeye tıklayın.
