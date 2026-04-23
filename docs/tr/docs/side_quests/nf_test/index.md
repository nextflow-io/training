# nf-test ile Test Etme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

İş akışınızın her parçasının beklenen şekilde çalıştığını sistematik olarak doğrulayabilmek; tekrarlanabilirlik ve uzun vadeli bakım açısından kritik öneme sahiptir ve geliştirme sürecinde büyük bir yardım sağlayabilir.

Testin neden bu kadar önemli olduğunu bir dakika düşünelim. Bir iş akışı geliştiriyorsanız, yapacağınız ilk şeylerden biri geçerli olduğunu bildiğiniz ve bir sonuç üretmesi gereken bazı test verileri bulmaktır. İlk süreci pipeline'a ekler ve çalışması için girdilerinize bağlarsınız. Ardından her şeyin çalışıp çalışmadığını kontrol etmek için test verileri üzerinde çalıştırırsınız. Bunun çalıştığını varsayarak bir sonraki sürece geçer ve test verilerini tekrar çalıştırırsınız. Memnun olduğunuz bir pipeline elde edene kadar bu süreci tekrarlarsınız.

Ardından, `--skip_process` gibi basit bir doğru/yanlış parametresi ekleyebilirsiniz. Artık pipeline'ı iki kez çalıştırmanız gerekir; beklenen şekilde çalıştığından emin olmak için her parametreyle bir kez. Ama bir dakika, `--skip_process` parametresinin süreci gerçekten atlayıp atlamadığını nasıl kontrol ederiz? Çıktılara bakmamız veya log dosyalarını kontrol etmemiz gerekir! Bu hem zahmetli hem de hataya açıktır.

Pipeline'ınızı geliştirdikçe, her iterasyonu manuel olarak test etmek hem yavaş hem de hataya açık hale gelecek kadar karmaşık bir hal alacaktır. Üstelik bir hata bulursanız, hatanın pipeline'ınızın tam olarak nereden kaynaklandığını tespit etmek çok zor olacaktır. İşte burada test devreye girer.

Test, pipeline'ınızın her parçasının beklenen şekilde çalışıp çalışmadığını sistematik olarak kontrol etmenizi sağlar. İyi yazılmış testlerin bir geliştirici için faydaları çok büyüktür:

- **Güven**: Testler tüm pipeline'ı kapsadığından, bir şeyi değiştirmenin başka bir şeyi etkilemediğinden emin olabilirsiniz.
- **Güvenilirlik**: Birden fazla geliştirici pipeline üzerinde çalıştığında, diğer geliştiricilerin pipeline'ı ve her bileşeni bozmadığını bilirler.
- **Şeffaflık**: Testler, pipeline'ın nerede başarısız olduğunu gösterir ve sorunu takip etmeyi kolaylaştırır. Aynı zamanda bir süreç veya iş akışının nasıl çalıştırılacağını gösteren bir dokümantasyon işlevi görürler.
- **Hız**: Testler otomatik olduğundan, çok hızlı ve tekrar tekrar çalıştırılabilirler. Yeni hatalar ekleme korkusu olmadan hızlıca iterasyon yapabilirsiniz.

Yazabileceğimiz pek çok farklı test türü vardır:

1. **Modül düzeyinde testler**: Bireysel süreçler için
2. **İş akışı düzeyinde testler**: Tek bir iş akışı için
3. **Pipeline düzeyinde testler**: Bir bütün olarak pipeline için
4. **Performans testleri**: Pipeline'ın hızı ve verimliliği için
5. **Stres testleri**: Sınırlarını belirlemek için pipeline'ın aşırı koşullar altındaki performansını değerlendirme

Bireysel süreçleri test etmek, diğer dillerdeki birim testlerine benzerdir. İş akışını veya tüm pipeline'ı test etmek ise diğer dillerde entegrasyon testleri olarak adlandırılan şeye benzer; burada bileşenlerin etkileşimlerini test ederiz.

[**nf-test**](https://www.nf-test.com/), modül, iş akışı ve pipeline düzeyinde testler yazmanıza olanak tanıyan bir araçtır. Kısacası, pipeline'ın her bir parçasının _izole olarak_ beklenen şekilde çalışıp çalışmadığını sistematik olarak kontrol etmenizi sağlar.

### Öğrenme hedefleri

Bu yan görevde, pipeline için iş akışı düzeyinde bir test ve çağırdığı üç süreç için modül düzeyinde testler yazmak amacıyla nf-test kullanmayı öğreneceksiniz.

Bu yan görevin sonunda aşağıdaki teknikleri etkin bir şekilde kullanabileceksiniz:

- Projenizde nf-test'i başlatma
- Modül düzeyinde ve iş akışı düzeyinde testler oluşturma
- Yaygın assertion türleri ekleme
- Snapshot'ların ne zaman, içerik assertion'larının ne zaman kullanılacağını anlama
- Tüm proje için testleri çalıştırma

Bu beceriler, pipeline projelerinizde kapsamlı bir test stratejisi uygulamanıza yardımcı olacak ve bunların daha dayanıklı ve sürdürülebilir olmasını sağlayacaktır.

### Ön koşullar

Bu yan göreve başlamadan önce şunları yapmanız gerekir:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmanız.
- Temel Nextflow kavramları ve mekanizmalarını (süreçler, kanallar, operatörler, dosyalarla çalışma, meta veri) rahatça kullanabilmeniz.

---

## 0. Başlarken

#### Eğitim kod alanını açın

Henüz yapmadıysanız, eğitim ortamını [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı şekilde açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Bu eğitim için dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/nf-test
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri inceleyin

Bir ana iş akışı dosyası ve pipeline'ın girdisini içeren `greetings.csv` adlı bir CSV dosyası bulacaksınız.

```console title="Directory contents"
.
├── greetings.csv
└── main.nf
```

Dosyaların ayrıntılı açıklaması için [Hello Nextflow'daki ısınma bölümüne](../hello_nextflow/00_orientation.md) bakın.

Test edeceğimiz iş akışı, [Hello Workflow](../hello_nextflow/03_hello_workflow.md) bölümünde oluşturulan Hello iş akışının bir alt kümesidir.

??? example "Hello Nextflow iş akışı ne yapar?"

    [Hello Nextflow](../hello_nextflow/index.md) eğitimini yapmadıysanız, bu basit iş akışının ne yaptığına dair kısa bir genel bakış:

    İş akışı, selamlamaları içeren bir CSV dosyası alır, bunlar üzerinde dört ardışık dönüşüm adımı çalıştırır ve eğlenceli bir karakterin selamlamaları söylediğini gösteren ASCII resmi içeren tek bir metin dosyası çıktısı verir.

    Dört adım, ayrı modül dosyalarında saklanan Nextflow süreçleri (`sayHello`, `convertToUpper`, `collectGreetings` ve `cowpy`) olarak uygulanmıştır.

    1. **`sayHello`:** Her selamlamayı kendi çıktı dosyasına yazar (örn. "Hello-output.txt")
    2. **`convertToUpper`:** Her selamlamayı büyük harfe dönüştürür (örn. "HELLO")
    3. **`collectGreetings`:** Tüm büyük harfli selamlamaları tek bir toplu dosyada toplar
    4. **`cowpy`:** `cowpy` aracını kullanarak ASCII sanatı oluşturur

    Sonuçlar `results/` adlı bir dizine yayımlanır ve pipeline'ın nihai çıktısı (varsayılan parametrelerle çalıştırıldığında), büyük harfli selamlamaları söyleyen bir karakterin ASCII sanatını içeren düz bir metin dosyasıdır.

    Bu yan görevde, yalnızca ilk iki süreci içeren Hello iş akışının ara bir formunu kullanıyoruz. <!-- TODO: change this to use the full finished workflow as suggested in https://github.com/nextflow-io/training/issues/735 -->

Üzerinde çalışacağımız alt küme iki süreçten oluşmaktadır: `sayHello` ve `convertToUpper`.
Tam iş akışı kodunu aşağıda görebilirsiniz.

??? example "İş akışı kodu"

    ```groovy title="main.nf"
    /*
    * Pipeline parametreleri
    */
    params.input_file = "greetings.csv"

    /*
    * Standart çıktıya 'Hello World!' yazdırmak için echo kullan
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

İş akışının beklenen şekilde çalıştığından emin olmak için çalıştıralım.

```bash
nextflow run main.nf
```

```console title="Result of running the workflow"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

TEBRİKLER! Az önce bir test çalıştırdınız!

"Dur bir dakika, ben sadece iş akışını çalıştırdım ve çalıştı! Bu nasıl bir test oluyor?"

Güzel soru!

Az önce ne olduğunu inceleyelim.

İş akışını varsayılan parametrelerle çalıştırdınız, çalıştığını doğruladınız ve sonuçlardan memnunsunuz. İşte testin özü budur. Hello Nextflow eğitim kursunu tamamladıysanız, her bölüme başlamadan önce her şeyin doğru kurulduğunu doğrulamak için başlangıç noktası olarak kullandığımız iş akışını her zaman çalıştırdığımızı fark etmişsinizdir.

Test yazılımı esasen bu süreci bizim için otomatikleştirir.

#### Görevi inceleyin

Göreviniz, yapılacak olası değişiklikler durumunda her parçanın beklenen şekilde çalışmaya devam ettiğini kolayca doğrulayabilmek için nf-test kullanarak bu iş akışına standartlaştırılmış testler eklemektir.

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Hazırlık kontrol listesi

Başlamaya hazır mısınız?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Kod alanım çalışıyor
- [ ] Çalışma dizinimi uygun şekilde ayarladım
- [ ] İş akışını başarıyla çalıştırdım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. `nf-test`'i Başlatma

`nf-test` paketi, projemiz için test geliştirmeye başlamamız amacıyla birkaç şeyi ayarlayan bir başlatma komutu sağlar.

```bash
nf-test init
```

Bu komut aşağıdaki çıktıyı üretmelidir:

```bash
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

Ayrıca bir yapılandırma dosyası taslağı içeren bir `tests` dizini oluşturur.

### 1.1. Bir nf-test oluşturma

`nf-test`, nf-test dosyaları oluşturmak için bir araç seti ile birlikte gelir ve işin büyük bölümünü bizim yerimize yapar. Bu araçlar `generate` alt komutu altında yer alır. Pipeline için bir test oluşturalım:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

Bu komut, `tests` dizini içinde bir `main.nf.test` dosyası oluşturacaktır. Bu, pipeline düzeyindeki test dosyamızdır. `tree tests/` komutunu çalıştırırsanız şuna benzer bir şey görmelisiniz:

```console title="Test directory contents"
tests/
├── main.nf.test
└── nextflow.config
```

`main.nf.test` dosyası, pipeline düzeyindeki test dosyamızdır. Açıp içeriğine bakalım.

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

    test("Should run without failures") {

        when {
            params {
                // parametreleri burada tanımlayın. Örnek:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

Test dosyasının yapısını anlamak için bir dakika ayıralım.

`nextflow_pipeline` bloğu, tüm pipeline düzeyindeki testler için giriş noktasıdır. Şunları içerir:

- `name`: Testin adı.
- `script`: Pipeline betiğinin yolu.

`test` bloğu, gerçek testin kendisidir. Şunları içerir:

- `when`: Testin çalıştırılması gereken koşullar. Bu, pipeline'ı çalıştırmak için kullanılacak parametreleri içerir.
- `then`: Yapılması gereken assertion'lar. Bu, pipeline'ın beklenen sonuçlarını içerir.

Testin mantığını sade bir dille şöyle okuyabiliriz:
"Bu _parametreler_ bu _pipeline_'a sağlandığında, bu sonuçları görmeyi bekliyoruz."

Bu işlevsel bir test değildir; bir sonraki bölümde bunu nasıl işlevsel hale getireceğimizi göstereceğiz.

### Test Adları Hakkında Bir Not

Yukarıdaki örnekte, yalnızca pipeline'ın başarıyla çalışıp çalışmadığını kontrol eden temel bir test için uygun olan varsayılan "Should run without failures" adını kullandık. Ancak daha spesifik test senaryoları ekledikçe, gerçekte neyi test ettiğimizi belirten daha açıklayıcı adlar kullanmalıyız. Örneğin:

- "Should convert input to uppercase" - belirli bir işlevselliği test ederken
- "Should handle empty input gracefully" - uç durumları test ederken
- "Should respect max memory parameter" - kaynak kısıtlamalarını test ederken
- "Should create expected output files" - dosya oluşturmayı test ederken

İyi test adları şunları yapmalıdır:

1. Beklenen davranışın ne olduğunu netleştirmek için "Should" ile başlamalıdır
2. Test edilen belirli işlevselliği veya senaryoyu açıklamalıdır
3. Test başarısız olduğunda hangi işlevselliğin bozulduğunu anlayacak kadar açık olmalıdır

İlerleyen bölümlerde daha fazla assertion ve spesifik test senaryoları ekledikçe, her testin neyi doğruladığını netleştirmek için bu daha açıklayıcı adları kullanacağız.

### 1.2. Testi çalıştırma

Testi çalıştırıp ne olduğunu görelim.

```bash
nf-test test tests/main.nf.test
```

```console title="nf-test pipeline fail"
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

Test başarısız oldu! Ne oldu?

1. nf-test, `when` bloğundaki ayarları kullanarak pipeline'ı olduğu gibi çalıştırmaya çalıştı:

```groovy title="tests/main.nf.test"
when {
    params {
        // parametreleri burada tanımlayın. Örnek:
        // outdir = "tests/results"
    }
}
```

2. nf-test, pipeline'ın durumunu kontrol etti ve `then` bloğuyla karşılaştırdı:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

nf-test'in pipeline'ın başarısız olduğunu bildirdiğine ve Nextflow'dan gelen hata mesajını sağladığına dikkat edin:

```console title="Error"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

Peki sorun neydi? Pipeline'ın proje dizininde bir `greetings.csv` dosyası olduğunu hatırlayın. nf-test pipeline'ı çalıştırdığında bu dosyayı arayacak, ancak bulamayacaktır. Dosya orada, peki ne oluyor? Yola bakarsak, testin `./nf-test/tests/uzunHashDizisi/` yolunda gerçekleştiğini görebiliriz. Nextflow gibi, nf-test de her şeyi izole tutmak için her test için yeni bir dizin oluşturur. Veri dosyası orada bulunmadığından, orijinal testteki dosya yolunu düzeltmemiz gerekir.

Test dosyasına geri dönelim ve `when` bloğundaki dosya yolunu değiştirelim.

Testte pipeline'ın köküne nasıl işaret edeceğimizi merak ediyor olabilirsiniz. Bu yaygın bir durum olduğundan, nf-test işimizi kolaylaştırmak için kullanabileceğimiz bir dizi global değişken sunar. Tam listeyi [buradan](https://www.nf-test.com/docs/testcases/global_variables/) bulabilirsiniz; ancak şimdilik pipeline projesinin kökü anlamına gelen `projectDir` değişkenini kullanacağız.

=== "Sonra"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
    when {
        params {
            input_file = "${projectDir}/greetings.csv"
        }
    }
    ```

=== "Önce"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
    when {
        params {
            // parametreleri burada tanımlayın. Örnek:
            // outdir = "tests/results"
        }
    }
    ```

Çalışıp çalışmadığını görmek için testi tekrar çalıştıralım.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

Başarılı! Pipeline başarıyla çalışıyor ve test geçiyor. İstediğiniz kadar çalıştırın, her seferinde aynı sonucu alacaksınız!

Varsayılan olarak Nextflow çıktısı gizlenir, ancak nf-test'in iş akışını gerçekten çalıştırdığına kendinizi ikna etmek için `--verbose` bayrağını kullanabilirsiniz:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Pipeline runs all processes"
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

### 1.3. Assertion ekleme

Basit bir kontrol, pipeline'ımızın beklediğimiz tüm süreçleri çalıştırıp çalıştırmadığını ve hiçbirini sessizce atlayıp atlamadığını doğrulamaktır. Pipeline'ımızın 6 süreç çalıştırdığını hatırlayın; 3 selamlama için birer tane `sayHello` ve birer tane `convertToUpper`.

Pipeline'ın beklenen sayıda süreci çalıştırıp çalıştırmadığını kontrol etmek için testimize bir assertion ekleyelim. Ayrıca test adını neyi test ettiğimizi daha iyi yansıtacak şekilde güncelleyeceğiz.

=== "Sonra"

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

=== "Önce"

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

Test adı artık gerçekte neyi doğruladığımızı daha iyi yansıtıyor; yalnızca pipeline'ın başarısız olmadan çalışıp çalışmadığını değil, beklenen sayıda süreci çalıştırıp çalıştırmadığını da kontrol ediyoruz.

Çalışıp çalışmadığını görmek için testi tekrar çalıştıralım.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with assertions"
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

Başarılı! Pipeline başarıyla çalışıyor ve test geçiyor. Artık pipeline'ın genel durumunun yanı sıra ayrıntılarını da test etmeye başladık.

### 1.4. Çıktıyı test etme

Çıktı dosyasının oluşturulup oluşturulmadığını kontrol etmek için testimize bir assertion ekleyelim. Sonuçların yorumlanmasını kolaylaştırmak için bunu bilgilendirici bir adla ayrı bir test olarak ekleyeceğiz.

=== "Sonra"

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

=== "Önce"

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

Çalışıp çalışmadığını görmek için testi tekrar çalıştırın.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with file assertions"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

Başarılı! Pipeline başarıyla tamamlandığı, doğru sayıda süreç çalıştığı ve çıktı dosyaları oluşturulduğu için testler geçiyor. Bu aynı zamanda testleriniz için bu bilgilendirici adları sağlamanın ne kadar yararlı olduğunu da gösteriyor.

Bu yalnızca yüzeyin bir kısmıdır; pipeline'ın ayrıntılarını kontrol etmek için assertion yazmaya devam edebiliriz, ancak şimdilik pipeline'ın iç kısımlarını test etmeye geçelim.

### Özetle

Bir pipeline için nf-test nasıl yazılır öğrendiniz.

### Sırada ne var?

Bir Nextflow sürecini nasıl test edeceğinizi öğrenin.

---

## 2. Bir Nextflow Sürecini Test Etme

Pipeline'ın her parçası için test yazmak zorunda değiliz, ancak ne kadar çok testimiz olursa pipeline hakkında o kadar kapsamlı olabilir ve beklenen şekilde çalıştığından o kadar emin olabiliriz. Bu bölümde, pipeline'daki her iki süreci bireysel birimler olarak test edeceğiz.

### 2.1. `sayHello` sürecini test etme

`sayHello` süreciyle başlayalım.

Süreç için testler oluşturmak amacıyla `nf-test generate` komutunu tekrar kullanalım.

```bash
nf-test generate process main.nf
```

```console title="Output"
> nf-test generate process main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.sayhello.nf.test
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

Şimdilik `main.sayhello.nf.test` dosyasındaki `sayhello` sürecine odaklanalım.

Dosyayı açıp içeriğine bakalım.

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

    test("Should run without failures") {

        when {
            params {
                // parametreleri burada tanımlayın. Örnek:
                // outdir = "tests/results"
            }
            process {
                """
                // sürecin girdilerini burada tanımlayın. Örnek:
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

Daha önce olduğu gibi, test ayrıntılarıyla başlıyoruz, ardından `when` ve `then` blokları geliyor. Ancak, sürece girdileri tanımlamamıza olanak tanıyan ek bir `process` bloğumuz da var.

Çalışıp çalışmadığını görmek için testi çalıştıralım.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
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

`sayHello` süreci 1 girdi bildirdiği ancak 0 argümanla çağrıldığı için test başarısız oldu. Sürece bir girdi ekleyerek bunu düzeltelim. [Hello Workflow](../hello_nextflow/03_hello_workflow.md) bölümünden (ve yukarıdaki ısınma bölümünden) hatırlayacağınız üzere, `sayHello` sürecimiz sağlamamız gereken tek bir değer girdisi alır. Ayrıca test adını neyi test ettiğimizi daha iyi yansıtacak şekilde düzeltmeliyiz.

=== "Sonra"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10"
        test("Should run without failures and produce correct output") {

            when {
                params {
                    // parametreleri burada tanımlayın. Örnek:
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

=== "Önce"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
        test("Should run without failures") {

            when {
                params {
                    // parametreleri burada tanımlayın. Örnek:
                    // outdir = "tests/results"
                }
                process {
                    """
                    // sürecin girdilerini burada tanımlayın. Örnek:
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

Çalışıp çalışmadığını görmek için testi tekrar çalıştıralım.

```console title="nf-test pipeline pass"
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

Başarılı! `sayHello` süreci başarıyla çalıştığı ve çıktı oluşturulduğu için test geçiyor.

### 2.2. Test tarafından oluşturulan snapshot'ı inceleme

`tests/main.sayhello.nf.test` dosyasına bakarsak, assertion bloğunda bir `snapshot()` metodu kullandığını görebiliriz:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

Bu, nf-test'e `sayHello` sürecinin çıktısının bir snapshot'ını oluşturmasını söyler. Snapshot dosyasının içeriğine bakalım.

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

Burada yazdırmayacağız, ancak süreç ve süreç çıktılarının ayrıntılarını içeren bir JSON dosyası görmelisiniz. Özellikle şuna benzer bir satır görebiliriz:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

Bu, açıkça test ettiğimiz `sayHello` süreci tarafından oluşturulan çıktıları temsil eder. Testi yeniden çalıştırırsak, program yeni çıktının başlangıçta kaydedilen çıktıyla eşleşip eşleşmediğini kontrol edecektir. Bu, süreç çıktılarının değişmediğini test etmenin hızlı ve basit bir yoludur; bu nedenle nf-test bunu varsayılan olarak sunar.

!!!warning "Uyarı"

    Bu, orijinal çalıştırmada kaydettiğimiz çıktının doğru olduğundan emin olmamız gerektiği anlamına gelir!

Gelecekteki geliştirme sürecinde, kodda çıktının farklı olmasına neden olan bir şey değişirse, test başarısız olacak ve değişikliğin beklenen bir değişiklik olup olmadığını belirlememiz gerekecektir.

- Kodda bir şeyin bozulduğu ortaya çıkarsa, düzeltmemiz gerekecektir; düzeltilen kodun testi geçmesi beklenir.
- Beklenen bir değişiklikse (örn. araç iyileştirildi ve sonuçlar daha iyi) snapshot'ı yeni çıktıyı eşleştirilecek referans olarak kabul edecek şekilde güncellememiz gerekecektir. nf-test bu amaç için `--update-snapshot` parametresine sahiptir.

Testi tekrar çalıştırabiliriz ve testin geçmesi gerekir:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

Başarılı! `sayHello` süreci başarıyla çalıştığı ve çıktı snapshot ile eşleştiği için test geçiyor.

### 2.3. Snapshot'lara Alternatif: Doğrudan İçerik Assertion'ları

Snapshot'lar çıktıdaki herhangi bir değişikliği yakalamak için harika olsa da, bazen tüm dosyanın eşleşmesi konusunda bu kadar katı olmadan belirli içerikleri doğrulamak istersiniz. Örneğin:

- Çıktının bazı kısımları değişebilir (zaman damgaları, rastgele kimlikler vb.) ancak belirli anahtar içerikler mevcut olmalıdır
- Çıktıdaki belirli kalıpları veya değerleri kontrol etmek istediğinizde
- Testin başarıyı neyin oluşturduğu konusunda daha açık olmasını istediğinizde

Testimizi belirli içerikleri kontrol edecek şekilde nasıl değiştirebileceğimiz aşağıda gösterilmiştir:

=== "Sonra"

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

=== "Önce"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
        test("Should run without failures and produce correct output") {

            when {
                params {
                    // parametreleri burada tanımlayın. Örnek:
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

nf-test'in süreç çıktılarını liste listesi olarak gördüğüne dikkat edin; bu nedenle `process.out[0][0]`, bu süreçten gelen ilk kanal öğesinin (veya 'yayının') ilk bölümünü alır.

Bu yaklaşım:

- Çıktıda tam olarak ne beklediğimizi açıkça ortaya koyar
- Çıktıdaki alakasız değişikliklere karşı daha dayanıklıdır
- Testler başarısız olduğunda daha iyi hata mesajları sağlar
- Daha karmaşık doğrulamalara (regex kalıpları, sayısal karşılaştırmalar vb.) olanak tanır

Çalışıp çalışmadığını görmek için testi çalıştıralım.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (7.196s)


SUCCESS: Executed 1 tests in 7.208s
```

### 2.4. `convertToUpper` sürecini test etme

`tests/main.converttoupper.nf.test` dosyasını açıp içeriğine bakalım:

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

    test("Should run without failures") {

        when {
            params {
                // parametreleri burada tanımlayın. Örnek:
                // outdir = "tests/results"
            }
            process {
                """
                // sürecin girdilerini burada tanımlayın. Örnek:
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

Bu, `sayHello` sürecine benzer bir test, ancak `convertToUpper` sürecini test ediyor. `sayHello`'da olduğu gibi, `convertToUpper` süreci de tek bir path girdisi alır ancak biz henüz bir tane belirtmedik; bu nedenle bu testin de başarısız olacağını biliyoruz.

Şimdi `convertToUpper` sürecine büyük harfe dönüştürmek istediğimiz bir metin içeren tek bir girdi dosyası sağlamamız gerekiyor. Bunu yapmanın pek çok yolu vardır:

- Test için özel bir dosya oluşturabiliriz
- Mevcut `data/greetings.csv` dosyasını yeniden kullanabiliriz
- Test içinde anında oluşturabiliriz

Şimdilik, pipeline düzeyindeki testte kullandığımız örneği kullanarak mevcut `data/greetings.csv` dosyasını yeniden kullanalım. Daha önce olduğu gibi, neyi test ettiğimizi daha iyi yansıtmak için testi adlandırabiliriz; ancak bu sefer (diğer süreçte yaptığımız gibi belirli dizeler kontrol etmek yerine) içeriği 'snapshot' olarak bırakalım.

=== "Sonra"

    ```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10"
        test("Should run without failures and produce correct output") {

            when {
                params {
                    // parametreleri burada tanımlayın. Örnek:
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

=== "Önce"

    ```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
        test("Should run without failures") {

            when {
                params {
                    // parametreleri burada tanımlayın. Örnek:
                    // outdir = "tests/results"
                }
                process {
                    """
                    // sürecin girdilerini burada tanımlayın. Örnek:
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

Ve testi çalıştırın!

```bash title="nf-test pipeline pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
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

`convertToUpper` süreci için `tests/main.converttoupper.nf.test.snap` konumunda bir snapshot dosyası oluşturduğumuza dikkat edin. Testi tekrar çalıştırırsak, nf-test'in yeniden geçtiğini görmeliyiz.

```bash title="nf-test process convertToUpper pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s
```

### Özetle

Bir Nextflow süreci için nasıl test yazılacağını ve bunların nasıl çalıştırılacağını öğrendiniz.

### Sırada ne var?

Her şeyi aynı anda nasıl test edeceğinizi öğrenin!

## 3. Tüm Depo için Testleri Çalıştırma

Her bileşen üzerinde nf-test çalıştırmak uygun olsa da zahmetli ve hataya açıktır. Her şeyi aynı anda test edemez miyiz?

Evet, edebiliriz!

nf-test'i tüm depo üzerinde çalıştıralım.

### 3.1. nf-test'i tüm depo üzerinde çalıştırma

`nf-test test` komutunu çalıştırarak nf-test'i tüm depo üzerinde çalıştırabiliriz.

```bash
nf-test test .
```

Mevcut dizinimizdeki her şeyi çalıştırmak için yalnızca `.` kullandığımıza dikkat edin. Bu, her testi dahil edecektir!

```console title="nf-test repo pass"
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

Buna bakın! Tek bir komutla her süreç için 1 ve tüm pipeline için 2 olmak üzere toplam 4 test çalıştırdık. Büyük bir kod tabanında bunun ne kadar güçlü olduğunu hayal edin!

---

## Özet

Bu yan görevde, bireysel süreçler için testler oluşturup çalıştırmak ve tüm pipeline için uçtan uca testler yapmak amacıyla nf-test'in özelliklerinden nasıl yararlanacağınızı öğrendiniz.
Artık çıktı doğrulamanın iki ana yaklaşımından, snapshot'lardan ve doğrudan içerik assertion'larından, ve her birinin ne zaman kullanılacağından haberdarsınız.
Ayrıca testleri tek tek veya tüm proje için nasıl çalıştıracağınızı da biliyorsunuz.

Bu teknikleri kendi çalışmalarınızda uygulamak şunları sağlamanıza olanak tanıyacaktır:

- Kodunuzun beklenen şekilde çalışması
- Değişikliklerin mevcut işlevselliği bozmaması
- Diğer geliştiricilerin güvenle katkıda bulunabilmesi
- Sorunların hızla tespit edilip düzeltilebilmesi
- Çıktı içeriğinin beklentilerle eşleşmesi

### Temel kalıplar

<!-- TODO: Can we add snippets of code below to illustrate? -->

1. Pipeline düzeyinde testler:
   - Temel başarı testi
   - Süreç sayısı doğrulaması
   - Çıktı dosyası varlık kontrolleri
2. Süreç düzeyinde testler
3. Çıktı doğrulamanın iki yaklaşımı:
   - Tam çıktı doğrulaması için snapshot kullanımı
   - Belirli içerik kontrolleri için doğrudan içerik assertion'ları kullanımı
4. Tek bir komutla bir depodaki tüm testleri çalıştırma

### Ek kaynaklar

Daha gelişmiş test özellikleri ve en iyi uygulamalar için [nf-test belgelerine](https://www.nf-test.com/) göz atın. Şunları yapmak isteyebilirsiniz:

- Testlerinize daha kapsamlı assertion'lar ekleme
- Uç durumlar ve hata koşulları için testler yazma
- Testleri otomatik olarak çalıştırmak için sürekli entegrasyon kurma
- İş akışı ve modül testleri gibi diğer test türleri hakkında bilgi edinme
- Daha gelişmiş içerik doğrulama tekniklerini keşfetme

**Unutmayın:** Testler, kodunuzun nasıl davranması gerektiğinin yaşayan bir dokümantasyonudur. Ne kadar çok test yazarsanız ve assertion'larınız ne kadar spesifik olursa, pipeline'ınızın güvenilirliğinden o kadar emin olabilirsiniz.

---

## Sırada ne var?

[Yan Görevler menüsüne](../) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
