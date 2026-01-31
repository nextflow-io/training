# Bölüm 4: Test ekleme

Bu kursun ilk bölümünde, tamamen doğrusal olan ve her örneğin verisini diğerlerinden bağımsız olarak işleyen bir varyant çağırma pipeline'ı oluşturdunuz.

İkinci bölümde, GATK ile ortak varyant çağırma işlemini gerçekleştirmek için kanalları ve kanal operatörlerini nasıl kullanacağınızı gösterdik.

Üçüncü bölümde, pipeline'ı modülerleştirdik.

Eğitimin bu bölümünde, Nextflow ile iyi entegre olan ve pipeline'ınıza hem modül düzeyinde hem de workflow düzeyinde testler eklemeyi kolaylaştıran bir test framework'ü olan [**nf-test**](https://www.nf-test.com/) kullanımını göstereceğiz. Eğitimin bu bölümünü takip edebilmek için Bölüm 1, Bölüm 2 ve Bölüm 3'ü tamamlamış olmanız ve ayrıca nf-test'in temellerini ve test yapmanın neden önemli olduğunu kapsayan [nf-test yan görevi](../../side_quests/nf-test.md)'ni tamamlamış olmanız gerekmektedir.

---

## 0. Isınma

!!! note

    Doğru çalışma dizininde olduğunuzdan emin olun:
    `cd /workspaces/training/nf4-science/genomics`

Bu eğitim kursunun önceki bölümlerini tamamladıysanız, uygun modül dizin yapısına sahip genomik pipeline'ın çalışan bir versiyonuna sahip olmalısınız.

??? abstract "Dizin içeriği"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf
    ```

Bu modül dizini, ihtiyacınız olması durumunda `solutions` dizininde bulunabilir.

Bölüm 3'teki ile aynı workflow ile başlayacağız, bunu sizin için `genomics-4.nf` dosyasında sağladık. [nf-test yan görevi](../../side_quests/nf-test.md) için olduğu gibi, bu pipeline'daki üç process'e birkaç farklı türde test ekleyeceğiz ve ayrıca bir workflow düzeyinde test ekleyeceğiz.

### 0.1. Workflow'un çalıştığını kontrol edin

Test eklemeye başlamadan önce, workflow'un beklendiği gibi çalıştığından emin olun.

```bash
nextflow run genomics-4.nf -resume
```

Bu eğitim kursunu baştan beri takip ediyorsanız, artık çok tanıdık gelmelidir.

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-4.nf` [gloomy_poincare] DSL2 - revision: 43203316e0

    executor >  local (7)
    [18/89dfa4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [30/b2522b] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    [a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Daha önce olduğu gibi, şimdi proje dizininizin içinde bir `work` dizini ve bir `results_genomics` dizini olacak. Aslında bu sonuçları daha sonra testlerimizde kullanacağız. Ancak şu andan itibaren pipeline'ı test etmek için `nf-test` paketini kullanacağız.

### 0.2. `nf-test`'i başlatın

[nf-test yan görevi](../../side_quests/nf-test.md)'nde olduğu gibi, `nf-test` paketini başlatmamız gerekiyor.

```bash
nf-test init
```

??? success "Komut çıktısı"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

    Project configured. Configuration is stored in nf-test.config
    ```

??? abstract "nf-test.config içeriği"

    ```groovy title="nf-test.config"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Ayrıca bir yapılandırma dosyası taslağı içeren bir `tests` dizini oluşturur.

### Çıkarım

Artık genomik pipeline'ımız için testler yazmaya hazırız.

### Sırada ne var?

Process çağrılarının başarılı olduğunu ve doğru çıktılar ürettiğini değerlendiren temel testler yazın.

---

## 1. Başarı ve eşleşen çıktılar için bir process'i test edin

BAM dosyaları için verimli rastgele erişimi etkinleştirmek üzere indeks dosyaları oluşturan `SAMTOOLS_INDEX` process'ini test ederek başlayacağız. Bu iyi bir ilk test durumudur çünkü:

1. Tek, iyi tanımlanmış bir girdiye (bir BAM dosyası) sahiptir
2. Tahmin edilebilir bir çıktı (bir BAI indeks dosyası) üretir
3. Çıktı, aynı girdiler için özdeş olmalıdır

### 1.1. Bir test dosyası taslağı oluşturun

İlk olarak, bir test dosyası taslağı oluşturun:

```bash
nf-test generate process modules/samtools/index/main.nf
```

??? success "Komut çıktısı"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/samtools/index/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/samtools/index/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Bu, `main.nf` ile aynı dizinde bir dosya oluşturur.
Dosya gezgininde dizine gidebilir ve dosyayı açabilirsiniz, aşağıdaki kodu içermelidir:

```groovy title="tests/modules/samtools/index/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"

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

Başlangıç onaylamaları [nf-test yan görevi](../../side_quests/nf-test.md)'nden tanıdık olmalıdır:

- `assert process.success` process'in başarılı bir şekilde çalışmasını ve herhangi bir hata olmadan tamamlanmasını beklediğimizi belirtir.
- `snapshot(process.out).match()` çalıştırma sonucunun önceki bir çalıştırmada (varsa) elde edilen sonuçla özdeş olmasını beklediğimizi belirtir.
  Bunu daha sonra daha ayrıntılı olarak tartışacağız.

Bunu başlangıç noktası olarak kullanarak, samtools index process için doğru test girdilerini ve varsa parametreleri eklememiz gerekiyor.

### 1.2. Test dosyasını taşıyın ve script yolunu güncelleyin

Testi doldurma işine başlamadan önce, dosyayı kesin konumuna taşımamız gerekiyor. Her modül için bir dizin eklememizin nedenlerinden biri, artık her modülün `main.nf` dosyasıyla birlikte konumlanmış bir `tests` dizininde testleri gönderebilmemizdir. O dizini oluşturun ve test dosyasını oraya taşıyın.

```bash
mkdir -p modules/samtools/index/tests
mv tests/modules/samtools/index/main.nf.test modules/samtools/index/tests/
```

Şimdi test dosyasının `script` bölümünü göreceli bir yola basitleştirebiliriz:

=== "Sonra"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "../main.nf"
    process "SAMTOOLS_INDEX"
    ```

=== "Önce"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"
    ```

Bu, teste modülün `main.nf` dosyasını tam yolu belirtmek zorunda kalmadan nerede bulacağını söyler.

### 1.3. SAMTOOLS_INDEX için test girdileri sağlayın

Taslak dosya, `samtools index` girdisine uygun gerçek bir test girdisiyle değiştirmemiz gereken bir yer tutucu içerir. Uygun girdi, `data/bam` dizininde mevcut olan bir BAM dosyasıdır.

=== "Sonra"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        input[0] = file("${projectDir}/data/bam/reads_son.bam")
        """
    }
    ```

=== "Önce"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        // define inputs of the process here. Example:
        // input[0] = file("test-file.txt")
        """
    }
    ```

### 1.4. Testi işlevselliğe göre adlandırın

Daha önce öğrendiğimiz gibi, testi test bağlamında anlam ifade eden bir şeyle yeniden adlandırmak iyi bir uygulamadır.

=== "Sonra"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should index reads_son.bam correctly") {
    ```

    Bu rastgele bir string alır, bu yüzden istediğimiz her şeyi koyabiliriz.
    Burada dosya adına ve formatına atıfta bulunmayı seçiyoruz.

=== "Önce"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should run without failures") {
    ```

### 1.5. Testi çalıştırın ve çıktıyı inceleyin

Testi çalıştırın:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Komut çıktısı"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.717s)
      Snapshots:
        1 created [Should index reads_son.bam correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 7.727s
    ```

Daha önce öğrendiğimiz gibi, bu process'in başarısı hakkındaki temel onaylamayı doğruladı ve process'in çıktısına dayalı bir anlık görüntü dosyası oluşturdu. Anlık görüntü dosyasının içeriğini `tests/modules/samtools/index/tests/main.nf.test.snap` dosyasında görebiliriz:

```json title="modules/samtools/index/tests/main.nf.test.snap" linenums="1"
{
  "Should index reads_son.bam correctly": {
    "content": [
      {
        "0": [
          [
            "reads_son.bam:md5,af5956d9388ba017944bef276b71d809",
            "reads_son.bam.bai:md5,a2ca7b84998218ee77eff14af8eb8ca2"
          ]
        ]
      }
    ],
    "meta": {
      "nf-test": "0.9.3",
      "nextflow": "25.10.2"
    },
    "timestamp": "2026-01-27T15:09:48.394063389"
  }
}
```

Ayrıca testi tekrar çalıştırabilir ve çıktının anlık görüntüyle özdeş olduğu için geçtiğini görebiliriz:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Komut çıktısı"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.938s)


    SUCCESS: Executed 1 tests in 7.987s
    ```

### 1.6. SAMTOOLS_INDEX'e daha fazla test ekleyin

Bazen çeşitli potansiyel sorunlar için test yaptığımızdan emin olmak için bir dizi farklı girdi dosyasını test etmek yararlıdır. Test verilerimizden üçlüdeki anne ve babanın BAM dosyaları için testler ekleyin.

```groovy
    test("Should index reads_mother.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_mother.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

    test("Should index reads_father.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_father.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Ardından testi tekrar çalıştırabilirsiniz:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Komut çıktısı"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.185s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (6.576s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (6.31s)
      Snapshots:
        2 created [Should index reads_father.bam correctly, Should index reads_mother.bam correctly]


    Snapshot Summary:
      2 created

    SUCCESS: Executed 3 tests in 20.117s
    ```

`--update-snapshot` parametresinin etkisine atıfta bulunan uyarıya dikkat edin.

!!! note

    Burada daha önce pipeline'ın bilimsel çıktılarını göstermek için kullandığımız test verilerini kullanıyoruz.
    Bu testleri bir üretim ortamında çalıştırmayı planlıyor olsaydık, test amaçları için daha küçük girdiler üretmiş olurduk.

    Genel olarak, process işlevselliğini değerlendirmek için gerekli ve yeterli olan en küçük veri parçalarını kullanarak birim testleri mümkün olduğunca hafif tutmak önemlidir, aksi takdirde toplam çalışma süresi oldukça ciddi bir şekilde artabilir.
    Düzenli olarak çalıştırılması çok uzun süren bir test paketi, uygunluk adına atlanması muhtemel bir test paketidir.

### Çıkarım

Bir genomik process için ilk modül testinizi yazdınız, `SAMTOOLS_INDEX`'in farklı BAM dosyaları için doğru şekilde indeks dosyaları oluşturduğunu doğruladınız. Test paketi şunları sağlar:

1. Process başarıyla çalışır
2. İndeks dosyaları oluşturulur
3. Çıktılar çalıştırmalar arasında tutarlıdır
4. Process tüm örnek BAM dosyaları için çalışır

### Sırada ne var?

Genomik workflow'umuzdaki diğer process'ler için testlerin nasıl yazılacağını öğrenin, zincirleme process'leri ele almak için kurulum yöntemini kullanın. Ayrıca çıktıların, özellikle VCF dosyalarımızın, beklenen varyant çağrılarını içerip içermediğini değerlendireceğiz.

---

## 2. Zincirleme bir process'e testler ekleyin ve içerik için test edin

`GATK_HAPLOTYPECALLER`'ı test etmek için, process'e girdi olarak `SAMTOOLS_INDEX` çıktısını sağlamamız gerekiyor. Bunu `SAMTOOLS_INDEX`'i çalıştırarak, çıktılarını alarak ve bunları workflow için test verileriyle saklayarak yapabiliriz. Bu aslında cilalı bir pipeline için önerilen yaklaşımdır, ancak nf-test, `setup` yöntemini kullanarak alternatif bir yaklaşım sağlar.

Setup yöntemiyle, test kurulumunun bir parçası olarak `SAMTOOLS_INDEX` process'ini tetikleyebilir ve ardından çıktısını `GATK_HAPLOTYPECALLER` için bir girdi olarak kullanabiliriz. Bunun bir maliyeti var: `GATK_HAPLOTYPECALLER` için testi her çalıştırdığımızda `SAMTOOLS_INDEX` process'ini çalıştırmak zorunda kalacağız. Ancak belki hala workflow'u geliştiriyoruz ve daha sonra değiştirmek zorunda kalabileceğimiz test verilerini önceden oluşturmak istemiyoruz. `SAMTOOLS_INDEX` process'i de çok hızlıdır, bu nedenle belki de çıktılarını önceden oluşturmanın ve saklamanın faydaları ihmal edilebilir düzeydedir. Setup yöntemi şu şekilde çalışır.

### 2.1. Test dosyasını oluşturun ve yerleştirin

Daha önce olduğu gibi, önce dosya taslağını oluştururuz:

```bash
nf-test generate process modules/gatk/haplotypecaller/main.nf
```

??? success "Komut çıktısı"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/haplotypecaller/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/haplotypecaller/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Bu, aşağıdaki test taslağını üretir:

```groovy title="tests/modules/gatk/haplotypecaller/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_HAPLOTYPECALLER"
    script "modules/gatk/haplotypecaller/main.nf"
    process "GATK_HAPLOTYPECALLER"

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

### 2.2. Test dosyasını taşıyın ve script yolunu güncelleyin

Modülün `main.nf` dosyasıyla birlikte konumlanmış test dosyası için bir dizin oluştururuz:

```bash
mkdir -p modules/gatk/haplotypecaller/tests
```

Ve test taslak dosyasını oraya taşırız:

```bash
mv tests/modules/gatk/haplotypecaller/main.nf.test modules/gatk/haplotypecaller/tests/
```

Son olarak, script yolunu güncellemeyi unutmayın:

=== "Sonra"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "../main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

=== "Önce"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "modules/gatk/haplotypecaller/main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

### 2.3. Setup yöntemini kullanarak girdiler sağlayın

`when` bloğundan önce bir `setup` bloğu ekliyoruz, burada orijinal girdi dosyalarımızdan birinde `SAMTOOLS_INDEX` process'inin bir çalıştırmasını tetikleyebiliriz. Ayrıca, daha önce olduğu gibi test adını anlamlı bir şeyle değiştirmeyi unutmayın.

=== "Sonra"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7" hl_lines="1-12"
        test("Should call son's haplotype correctly") {

            setup {
                run("SAMTOOLS_INDEX") {
                    script "../../../samtools/index/main.nf"
                    process {
                        """
                        input[0] =  file("${projectDir}/data/bam/reads_son.bam")
                        """
                    }
                }
            }

            when {
    ```

=== "Önce"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7"  hl_lines="1"
    test("Should run without failures") {

        when {
    ```

Daha sonra test girdilerini belirttiğimiz `when` bloğunda o process'in çıktısına başvurabiliriz:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="20"
        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }
```

Bu değişikliği yapın ve testi tekrar çalıştırın:

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Komut çıktısı"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Snapshots:
        1 created [Should call son's haplotype correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 40.555s
    ```

Ayrıca daha önce olduğu gibi bir anlık görüntü dosyası üretir.

### 2.4. Tekrar çalıştırın ve hatayı gözlemleyin

İlginç bir şekilde, aynı komutu tekrar çalıştırırsanız, bu sefer test başarısız olacaktır.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? failure "Komut çıktısı"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' FAILED (40.123s)

      java.lang.RuntimeException: Different Snapshot:
      [                                                                                           [
          {                                                                                           {
              "0": [                                                                                      "0": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ],                                                                                          ],
              "1": [                                                                                      "1": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "idx": [                                                                                    "idx": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "vcf": [                                                                                    "vcf": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ]                                                                                           ]
          }                                                                                           }
      ]                                                                                           ]

      Nextflow stdout:

      Nextflow stderr:


        Obsolete snapshots can only be checked if all tests of a file are executed successful.


    FAILURE: Executed 1 tests in 40.156s (1 failed)
    ```

Hata mesajı size iki çalıştırma için anlık görüntüler arasında farklılıklar olduğunu söyler; özellikle, VCF dosyaları için md5sum değerleri farklıdır.

Neden? Uzun lafın kısası, HaplotypeCaller aracı VCF başlığına her seferinde farklı olan bir zaman damgası ekler (tanım gereği).
Sonuç olarak, varyant çağrılarının kendileri açısından özdeş içeriğe sahip olsalar bile, dosyaların özdeş md5sum'lara sahip olmasını bekleyemeyiz.

Bununla nasıl başa çıkıyoruz?

### 2.5. Belirli bir varyantı kontrol etmek için bir içerik onaylama yöntemi kullanın

Sorunu çözmenin bir yolu [farklı türde bir onaylama](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions) kullanmaktır.
Bu durumda, özdeşlik onaylamak yerine belirli içerik için kontrol yapacağız.
Daha kesin olarak, aracın VCF dosyasının satırlarını okumasını ve belirli satırların varlığını kontrol etmesini sağlayacağız.

Pratikte, `then` bloğundaki ikinci onaylamayı aşağıdaki gibi değiştiriyoruz:

=== "Sonra"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3 4"
            then {
                assert process.success
                assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_son')
                assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3282	GT:DP:GQ:MIN_DP:PL	0/0:25:72:24:0,72,719')
            }
    ```

=== "Önce"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3"
    then {
        assert process.success
        assert snapshot(process.out).match()
    }
    ```

Burada VCF çıktı dosyasının tam içeriğini okuyoruz ve bir içerik eşleşmesi arıyoruz, bunu küçük bir test dosyasında yapmak sorun değil, ancak bunu daha büyük bir dosyada yapmak istemezsiniz.
Bunun yerine belirli satırları okumayı seçebilirsiniz.

Bu yaklaşım, test için 'sinyal' olarak ne kullanmak istediğimizi daha dikkatli seçmeyi gerektirir.
İyi tarafı, bir analiz aracının daha fazla geliştirme geçirirken 'zor' özellikleri (nadir varyantlar gibi) tutarlı bir şekilde tanımlayıp tanımlayamadığını büyük hassasiyetle test etmek için kullanılabilir.

### 2.6. Tekrar çalıştırın ve başarıyı gözlemleyin

Testi bu şekilde değiştirdiğimizde, testi birden çok kez çalıştırabiliriz ve tutarlı bir şekilde geçecektir.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Komut çıktısı"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)


    SUCCESS: Executed 1 tests in 40.555s
    ```

### 2.7. Daha fazla test ekleyin

Anne ve baba örnekleri için benzer testler ekleyin:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="43"
    test("Should call mother's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_mother.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3278	GT:DP:GQ:MIN_DP:PL	0/0:38:99:37:0,102,1530')
        }
    }

    test("Should call father's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_father.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3281	GT:DP:GQ:MIN_DP:PL	0/0:44:99:42:0,120,1800')
        }
    }
```

### 2.8. Test komutunu çalıştırın

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Komut çıktısı"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (41.47s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (45.556s)


    SUCCESS: Executed 3 tests in 127.586s
    ```

Bu, pipeline'daki bu ikinci adım için temel test planını tamamlar. Üçüncü ve son modül düzeyinde teste geçiyoruz!

### Çıkarım

Şunları yapmayı öğrendiniz:

1. Diğer process'lerin çıktılarına bağlı olan process'leri test etmek
2. VCF çıktı dosyalarında belirli genomik varyantları doğrulamak
3. Belirli içeriği kontrol ederek deterministik olmayan çıktıları ele almak
4. Birden çok örnek üzerinde varyant çağırmayı test etmek

### Sırada ne var?

Ortak genotipleme adımı için önceden oluşturulmuş test verilerini kullanan testlerin nasıl yazılacağını öğrenin.

---

## 3. Önceden oluşturulmuş test verilerini kullanın

Ortak genotipleme adımı için, farklı bir yaklaşım kullanacağız - önceden oluşturulmuş test verilerini kullanma. Bu genellikle şunlar için tercih edilir:

1. Birden çok bağımlılığı olan karmaşık process'ler
2. Çalışması uzun süren process'ler
3. Kararlı, üretim pipeline'ının parçası olan process'ler

### 3.1. Test verisi oluşturun

Bu bölümün başında oluşturduğumuz sonuçları inceleyin:

```bash
tree results_genomics/
```

```console title="Sonuçlar dizini içeriği"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam -> /workspaces/training/nf4-science/genomics/work/42/a3bf19dbfaf1f3672b16a5d5e6a8be/reads_father.bam
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/cf/289c2d264f496d60a69e3e9ba6463e/reads_father.bam.bai
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/work/af/f31a6ade82cc0cf853c4f61c8bc473/reads_mother.bam
    ├── reads_mother.bam.bai -> /workspaces/training/nf4-science/genomics/work/18/89dfa40a3def17e45421e54431a126/reads_mother.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/work/9f/9615dd553d6f13d8bec4f006ac395f/reads_son.bam
    └── reads_son.bam.bai -> /workspaces/training/nf4-science/genomics/work/4d/cb384a97db5687cc9daab002017c7c/reads_son.bam.bai

2 directories, 14 files
```

Ortak genotipleme adımı, indekslerle birlikte haplotip çağırma adımları tarafından üretilen VCF dosyalarına girdi olarak ihtiyaç duyar. O halde sahip olduğumuz sonuçları `jointgenotyping` modülünün test dizinine kopyalayalım.

```bash
mkdir -p modules/gatk/jointgenotyping/tests/inputs/
cp results_genomics/gvcf/*.g.vcf results_genomics/gvcf/*.g.vcf.idx modules/gatk/jointgenotyping/tests/inputs/
```

Şimdi bu dosyaları ortak genotipleme adımı için yazacağımız testin girdileri olarak kullanabiliriz.

### 3.2. Test dosyası taslağını oluşturun

Daha önce olduğu gibi, önce dosya taslağını oluştururuz:

```bash
nf-test generate process modules/gatk/jointgenotyping/main.nf
```

??? success "Komut çıktısı"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/jointgenotyping/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/jointgenotyping/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Bu, aşağıdaki test taslağını üretir:

```groovy title="tests/modules/gatk/jointgenotyping/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"

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

### 3.3. Test dosyasını taşıyın ve script yolunu güncelleyin

Bu sefer zaten modülün `main.nf` dosyasıyla birlikte konumlanmış testler için bir dizinimiz var, bu yüzden test taslak dosyasını oraya taşıyabiliriz:

```bash
mv tests/modules/gatk/jointgenotyping/main.nf.test modules/gatk/jointgenotyping/tests/
```

Ve script yolunu güncellemeyi unutmayın:

=== "Sonra"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "../main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

=== "Önce"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

### 3.4. Girdileri sağlayın

Process girdi tanımlarına göre girdileri doldurun ve testi buna göre yeniden adlandırın:

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="7"
    test("Should call trio's joint genotype correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf")
                ]
                input[1] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf.idx")
                ]
                input[2] = file("${projectDir}/data/ref/intervals.bed")
                input[3] = "family_trio"
                input[4] = file("${projectDir}/data/ref/ref.fasta")
                input[5] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[6] = file("${projectDir}/data/ref/ref.dict")
                """
            }
        }
```

### 3.5. İçerik onaylamalarını kullanın

Ortak genotipleme adımının çıktısı başka bir VCF dosyasıdır, bu yüzden yine bir içerik onaylaması kullanacağız.

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="25"
    then {
        assert process.success
        assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son')
        assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0')
    }
```

Çıktı dosyasındaki belirli bir varyantın içeriğini kontrol ederek, bu test şunları doğrular:

1. Ortak genotipleme process'i başarıyla çalışır
2. Çıktı VCF'si üç örneği de doğru sırayla içerir
3. Belirli bir varyant aşağıdaki bilgilerle doğru şekilde çağrılır:
   - Her örnek için doğru genotipler (baba için 0/1, anne ve oğul için 1/1)
   - Doğru okuma derinlikleri ve genotip kaliteleri
   - Alel frekansı (AF=0.833) gibi popülasyon düzeyinde istatistikler

Tüm dosyayı anlık görüntülemedik, ancak belirli bir varyantı kontrol ederek, ortak genotipleme process'inin beklendiği gibi çalıştığından emin olabiliriz.

### 3.6. Testi çalıştırın

```bash
nf-test test modules/gatk/jointgenotyping/tests/main.nf.test
```

??? success "Komut çıktısı"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (53.827s)


    SUCCESS: Executed 1 tests in 53.837s
    ```

Test geçer ve ortak genotipleme process'imizin şunları doğru yaptığını doğrular:

1. Bireysel örnek VCF'lerini birleştirir
2. Ortak varyant çağırma gerçekleştirir
3. Çalıştırmalar arasında tutarlı genotip çağrılarıyla çok örnekli bir VCF üretir

### Çıkarım

Şunları yapmayı biliyorsunuz:

- Testler için girdi olarak önceden oluşturulmuş sonuçları kullanmak
- Önceden oluşturulmuş test verilerini kullanarak testler yazmak

### Sırada ne var?

Tüm varyant çağırma pipeline'ının uçtan uca çalıştığını doğrulamak için bir workflow düzeyinde test ekleyin.

---

## 4. Workflow düzeyinde test ekleyin

Şimdi BAM dosyalarından ortak genotiplere kadar tam varyant çağırma pipeline'ını test edeceğiz. Bu şunları doğrular:

1. Tüm process'ler doğru şekilde birlikte çalışır
2. Veriler adımlar arasında düzgün akar
3. Son varyant çağrıları tutarlıdır

### 4.1. Workflow testini oluşturun

Tam pipeline için bir test dosyası oluşturun:

```bash
nf-test generate pipeline genomics-4.nf
```

??? success "Komut çıktısı"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/genomics-4.nf'
    Wrote pipeline test file '/workspaces/training/nf4-science/genomics/tests/genomics-4.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Bu, temel bir test taslağı oluşturur:

```groovy title="tests/genomics-4.nf.test" linenums="1"
nextflow_pipeline {

    name "Test Workflow genomics-4.nf"
    script "genomics-4.nf"

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

Sadece adı anlamlı bir şeye düzeltin (bunun neden yararlı olduğunu kısa süre sonra göreceksiniz).

=== "Sonra"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run the pipeline without failures") {
    ```

=== "Önce"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run without failures") {
    ```

!!! note

    Bu durumda test dosyası `nf-test`'in oluşturduğu yerde kalabilir.

### 4.2. Girdi parametrelerini belirtin

Hala girdileri belirtmemiz gerekiyor, bu modül düzeyinde testlere göre workflow düzeyinde biraz farklı yapılır.
Bir profil belirlemek de dahil olmak üzere bunu yapmanın birkaç yolu vardır.
Ancak, daha basit bir yol, `nf-test init`'in `tests` dizininde başlangıçta oluşturduğu `nextflow.config` dosyasında bir `params {}` bloğu kurmaktır.

```groovy title="tests/nextflow.config" linenums="1"
/*
========================================================================================
    Nextflow config file for running tests
========================================================================================
*/

// Output directory for workflow outputs
outputDir = 'results_genomics'

/*
 * Pipeline parametreleri
 */

params {
    // Birincil girdi (satır başına bir tane olmak üzere girdi dosyalarının dosyası)
    reads_bam = "${projectDir}/data/sample_bams.txt"

    // Yardımcı dosyalar
    reference = "${projectDir}/data/ref/ref.fasta"
    reference_index = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict = "${projectDir}/data/ref/ref.dict"
    intervals = "${projectDir}/data/ref/intervals.bed"

    // Son çıktı dosyası için temel ad
    cohort_name = "family_trio"
}
```

Testi çalıştırdığımızda, `nf-test` bu yapılandırma dosyasını alacak ve girdileri buna göre çekecektir.

### 4.3. Workflow testini çalıştırın

```bash
nf-test test tests/genomics-4.nf.test
```

??? success "Komut çıktısı"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (171.019s)


    SUCCESS: Executed 1 tests in 171.056s
    ```

Test geçer ve tam varyant çağırma pipeline'ımızın şunları yaptığını onaylar:

1. Tüm örnekleri başarıyla işler
2. Tüm adımları doğru şekilde zincirleme yapar

### 4.4. TÜM testleri çalıştırın

nf-test'in kolunda bir numara daha var. Tüm testleri aynı anda çalıştırabiliriz! `nf-test.config` dosyasını değiştirin, böylece nf-test her dizinde nf-test dosyalarını arar. Bunu `testsDir` parametresini değiştirerek yapabilirsiniz:

=== "Sonra"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "."
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

=== "Önce"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Şimdi, sadece nf-test'i çalıştırabiliriz ve depomuzda _her bir testi_ çalıştıracaktır:

```bash
nf-test test
```

??? success "Komut çıktısı"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (39.947s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (43.17s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (44.244s)

    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (61.129s)

    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (8.671s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (8.518s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (5.378s)

    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (169.714s)


    SUCCESS: Executed 8 tests in 380.801s
    ```

1 komutta 8 test! Çok sayıda test yapılandırmak için uzun zaman harcadık, ancak onları çalıştırmaya gelince çok hızlı ve kolaydı. Testleri bir kez yazmak için zaman harcadığımızı, böylece onları birçok kez çalıştırmaktan zaman kazanabileceğimizi görebilirsiniz.

Ayrıca, bunu otomatikleştirebiliriz! Siz veya bir meslektaşınız yeni kod eklemeye çalıştığında testlerin her seferinde çalıştığını hayal edin. Pipeline'larımızın yüksek bir standart korumasını sağlamamızın yolu budur.

## Çıkarım

Artık nf-test kullanarak genomik pipeline'ınız için çeşitli testlerin nasıl yazılacağını ve çalıştırılacağını biliyorsunuz. Bu test framework'ü, varyant çağırma workflow'unuzun farklı ortamlarda ve kod değişiklikleri yaptığınızda tutarlı, güvenilir sonuçlar üretmesini sağlamaya yardımcı olur.

Şu kritik bileşenleri test etmeyi öğrendiniz:

- Varyant çağırma için BAM dosyalarını hazırlayan `SAMTOOLS_INDEX` process'i
- Bireysel örneklerde varyantları tanımlayan `GATK_HAPLOTYPECALLER` process'i
- Bir kohort genelinde varyant çağrılarını birleştiren `GATK_JOINTGENOTYPING` process'i

Ayrıca genomik verilere özgü farklı test stratejileri uyguladınız:

- Zaman damgaları gibi deterministik olmayan öğelere rağmen VCF dosyalarının beklenen varyant çağrılarını içerdiğini doğrulamak
- İlgili örnekler arasında uygun varyant tanımlamasını sağlamak için bir aile üçlüsü veri setiyle test etmek
- Çıktı dosyalarınızda belirli genomik koordinatları ve varyant bilgilerini kontrol etmek

Bu test becerileri,
