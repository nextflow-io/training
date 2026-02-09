# Bölüm 3: Bir nf-core modülü kullanın

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core eğitim kursunun bu üçüncü bölümünde, mevcut bir nf-core modülünü nasıl bulacağınızı, kuracağınızı ve pipeline'ınızda nasıl kullanacağınızı gösteriyoruz.

nf-core ile çalışmanın büyük avantajlarından biri, [nf-core/modules](https://github.com/nf-core/modules) deposundan önceden oluşturulmuş, test edilmiş modüllerden yararlanabilme yeteneğidir.
Her süreci sıfırdan yazmak yerine, en iyi uygulamaları takip eden topluluk tarafından sürdürülen modülleri kurabilir ve kullanabilirsiniz.

Bunun nasıl çalıştığını göstermek için, `core-hello` pipeline'ındaki özel `collectGreetings` modülünü nf-core/modules'den `cat/cat` modülü ile değiştireceğiz.

??? info "Bu bölümden nasıl başlanır"

    Kursun bu bölümü, [Bölüm 2: Hello'yu nf-core için yeniden yazın](./02_rewrite_hello.md) bölümünü tamamladığınızı ve çalışan bir `core-hello` pipeline'ınız olduğunu varsayar.

    Bölüm 2'yi tamamlamadıysanız veya bu bölüm için yeni başlamak istiyorsanız, başlangıç noktanız olarak `core-hello-part2` çözümünü kullanabilirsiniz.
    `hello-nf-core/` dizini içinden bu komutu çalıştırın:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    Bu size modül eklemeye hazır, tamamen işlevsel bir nf-core pipeline'ı verir.
    Aşağıdaki komutu çalıştırarak başarıyla çalıştığını test edebilirsiniz:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Uygun bir nf-core modülü bulun ve kurun

İlk olarak, mevcut bir nf-core modülünü nasıl bulacağımızı ve pipeline'ımıza nasıl kuracağımızı öğrenelim.

Birden fazla selamlama dosyasını tek bir dosyada birleştirmek için Unix `cat` komutunu kullanan `collectGreetings` sürecini değiştirmeyi hedefliyoruz.
Dosyaları birleştirmek çok yaygın bir işlemdir, bu nedenle nf-core'da bu amaç için tasarlanmış bir modülün zaten var olabileceğini düşünmek mantıklıdır.

Hadi başlayalım.

### 1.1. nf-core web sitesinde mevcut modüllere göz atın

nf-core projesi, [https://nf-co.re/modules](https://nf-co.re/modules) adresinde merkezi bir modül kataloğu tutar.

Web tarayıcınızda modüller sayfasına gidin ve 'concatenate' aramak için arama çubuğunu kullanın.

![modül arama sonuçları](./img/module-search-results.png)

Gördüğünüz gibi, birçoğu çok özel dosya türlerini birleştirmek için tasarlanmış modüller olmak üzere oldukça fazla sonuç var.
Bunlar arasında, genel amaçlı olan `cat_cat` adlı bir modül görmelisiniz.

!!! note "Modül adlandırma kuralı"

    Alt çizgi (`_`) karakteri, modül adlarında eğik çizgi (`/`) karakterinin yerine kullanılır.

    nf-core modülleri, bir araç birden fazla komut sağladığında `yazılım/komut` adlandırma kuralını takip eder; örneğin `samtools/view` (samtools paketi, view komutu) veya `gatk/haplotypecaller` (GATK paketi, HaplotypeCaller komutu).
    Yalnızca bir ana komut sağlayan araçlar için modüller `fastqc` veya `multiqc` gibi tek seviye kullanır.

Modül belgelerini görüntülemek için `cat_cat` modül kutusuna tıklayın.

Modül sayfası şunları gösterir:

- Kısa bir açıklama: "Sıkıştırılmış veya sıkıştırılmamış dosyaların birleştirilmesi için bir modül"
- Kurulum komutu: `nf-core modules install cat/cat`
- Girdi ve çıktı kanal yapısı
- Mevcut parametreler

### 1.2. Komut satırından mevcut modülleri listeleyin

Alternatif olarak, nf-core araçlarını kullanarak doğrudan komut satırından modülleri de arayabilirsiniz.

```bash
nf-core modules list remote
```

Bu, nf-core/modules deposundaki tüm mevcut modüllerin bir listesini görüntüler, ancak aramak istediğiniz modülün adını zaten bilmiyorsanız biraz daha az kullanışlıdır.
Ancak, biliyorsanız, belirli modülleri bulmak için listeyi `grep`'e yönlendirebilirsiniz:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "Komut çıktısı"

    ```console
    │ cat/cat
    ```

Sadece `grep` yaklaşımının yalnızca arama terimini adında içeren sonuçları çekeceğini unutmayın; bu `cat_cat` için işe yaramaz.

### 1.3. Modül hakkında ayrıntılı bilgi alın

Komut satırından belirli bir modül hakkında ayrıntılı bilgi görmek için `info` komutunu kullanın:

```bash
nf-core modules info cat/cat
```

Bu, modülün girdileri, çıktıları ve temel kullanım bilgileri dahil olmak üzere modül hakkındaki belgeleri görüntüler.

??? success "Komut çıktısı"

    ```console

                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    ╭─ Module: cat/cat  ─────────────────────────────────────────────────╮
    │ 🌐 Repository: https://github.com/nf-core/modules.git              │
    │ 🔧 Tools: cat                                                      │
    │ 📖 Description: A module for concatenation of gzipped or           │
    │ uncompressed files                                                 │
    ╰────────────────────────────────────────────────────────────────────╯
                      ╷                                          ╷
    📥 Inputs        │Description                               │Pattern
    ╺━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
    input[0]         │                                          │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      meta  (map)     │Groovy Map containing sample information  │
                      │e.g. [ id:'test', single_end:false ]      │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      files_in  (file)│List of compressed / uncompressed files   │      *
                      ╵                                          ╵
                          ╷                                 ╷
    📥 Outputs           │Description                      │     Pattern
    ╺━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
    file_out             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      meta  (map)         │Groovy Map containing sample     │
                          │information                      │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      ${prefix}  (file)   │Concatenated file. Will be       │ ${file_out}
                          │gzipped if file_out ends with    │
                          │".gz"                            │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
    versions             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      versions.yml  (file)│File containing software versions│versions.yml
                          ╵                                 ╵

    💻  Installation command: nf-core modules install cat/cat

    ```

Bu, web sitesinde bulabileceğiniz bilgilerle tamamen aynıdır.

### 1.4. cat/cat modülünü kurun

Artık istediğimiz modülü bulduğumuza göre, onu pipeline'ımızın kaynak koduna eklememiz gerekiyor.

İyi haber şu ki, nf-core projesi bu kısmı kolaylaştırmak için bazı araçlar içeriyor.
Özellikle, `nf-core modules install` komutu, kodu almayı ve tek bir adımda projenizde kullanılabilir hale getirmeyi otomatikleştirmeyi mümkün kılar.

Pipeline dizininize gidin ve kurulum komutunu çalıştırın:

```bash
cd core-hello
nf-core modules install cat/cat
```

Araç önce bir depo türü belirtmenizi isteyebilir.
(İstemezse, "Son olarak, araç modülü kurmaya devam edecektir." bölümüne atlayın.)

??? success "Komut çıktısı"

    ```console

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.4.1 - https://nf-co.re


    WARNING  'repository_type' not defined in .nf-core.yml
    ? Is this repository a pipeline or a modules repository? (Use arrow keys)
    » Pipeline
      Modules repository
    ```

Eğer öyleyse, varsayılan yanıtı (`Pipeline`) kabul etmek için enter tuşuna basın ve devam edin.

Araç daha sonra gelecekte bu istemden kaçınmak için projenizin yapılandırmasını değiştirmeyi teklif edecektir.

??? success "Komut çıktısı"

    ```console
        INFO     To avoid this prompt in the future, add the 'repository_type' key to your .nf-core.yml file.
        ? Would you like me to add this config now? [y/n] (y):
    ```

Bu kullanışlı araçtan yararlanmakta fayda var!
Varsayılan yanıtı (evet) kabul etmek için enter tuşuna basın.

Son olarak, araç modülü kurmaya devam edecektir.

??? success "Komut çıktısı"

    ```console
    INFO Config added to '.nf-core.yml'
    INFO Reinstalling modules found in 'modules.json' but missing from directory:
    INFO Installing 'cat/cat'
    INFO Use the following statement to include this module:

        include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
    ```

Komut otomatik olarak:

- Modül dosyalarını `modules/nf-core/cat/cat/` dizinine indirir
- Kurulu modülü izlemek için `modules.json` dosyasını günceller
- İş akışınızda kullanmanız için doğru `include` ifadesini size sağlar

!!! tip

    Modül kurulum komutunu çalıştırmadan önce mevcut çalışma dizininizin pipeline projenizin kök dizini olduğundan her zaman emin olun.

Modülün doğru şekilde kurulduğunu kontrol edelim:

```bash
tree -L 4 modules
```

??? abstract "Dizin içeriği"

    ```console
    modules
    ├── local
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nf-core
        └── cat
            └── cat
                ├── environment.yml
                ├── main.nf
                ├── meta.yml
                └── tests

    5 directories, 7 files
    ```

Ayrıca nf-core yardımcı programından yerel olarak kurulu modülleri listelemesini isteyerek kurulumu doğrulayabilirsiniz:

```bash
nf-core modules list local
```

??? success "Komut çıktısı"

    ```console
    INFO     Repository type: pipeline
    INFO     Modules installed in '.':

    ┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
    ┃ Module Name ┃ Repository      ┃ Version SHA ┃ Message                                ┃ Date       ┃
    ┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
    │ cat/cat     │ nf-core/modules │ 41dfa3f     │ update meta.yml of all modules (#8747) │ 2025-07-07 │
    └─────────────┴─────────────────┴─────────────┴────────────────────────────────────────┴────────────┘
    ```

Bu, `cat/cat` modülünün artık projenizin kaynak kodunun bir parçası olduğunu doğrular.

Ancak, yeni modülü gerçekten kullanmak için onu pipeline'ımıza aktarmamız gerekiyor.

### 1.5. Modül içe aktarmalarını güncelleyin

`core-hello/workflows/hello.nf` iş akışının içe aktarma bölümünde `collectGreetings` modülü için `include` ifadesini `CAT_CAT` için olanla değiştirelim.

Hatırlatma olarak, modül kurulum aracı bize kullanmamız gereken tam ifadeyi verdi:

```groovy title="Kurulum komutu tarafından üretilen içe aktarma ifadesi"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

nf-core kuralının, modülleri içe aktarırken modül adları için büyük harf kullanmak olduğunu unutmayın.

[core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) dosyasını açın ve aşağıdaki değişikliği yapın:

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

nf-core modülü için yolun yerel modüllerden nasıl farklı olduğuna dikkat edin:

- **nf-core modülü**: `'../modules/nf-core/cat/cat/main'` (`main.nf` dosyasına referans verir)
- **Yerel modül**: `'../modules/local/collectGreetings.nf'` (tek dosya referansı)

Modül artık iş akışı için kullanılabilir, bu nedenle tek yapmamız gereken `collectGreetings` çağrısını `CAT_CAT` kullanacak şekilde değiştirmek. Değil mi?

O kadar hızlı değil.

Bu noktada, koda dalıp düzenlemeye başlamak cazip gelebilir, ancak yeni modülün ne beklediğini ve ne ürettiğini dikkatlice incelemek için bir an durmaya değer.

Bunu ayrı bir bölüm olarak ele alacağız çünkü henüz ele almadığımız yeni bir mekanizma içeriyor: metadata map'ler.

!!! note

    İsteğe bağlı olarak `collectGreetings.nf` dosyasını silebilirsiniz:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    Ancak, yerel ve nf-core modülleri arasındaki farkları anlamak için bir referans olarak tutmak isteyebilirsiniz.

### Özet

Bir nf-core modülünü nasıl bulacağınızı ve projenizde kullanılabilir hale getireceğinizi biliyorsunuz.

### Sırada ne var?

Yeni bir modülün ne gerektirdiğini değerlendirin ve onu bir pipeline'a entegre etmek için gerekli önemli değişiklikleri belirleyin.

---

## 2. Yeni modülün gereksinimlerini değerlendirin

Özellikle, modülün **arayüzünü**, yani girdi ve çıktı tanımlarını incelememiz ve değiştirmeye çalıştığımız modülün arayüzüyle karşılaştırmamız gerekiyor.
Bu, yeni modülü doğrudan yerine koyabileceğimizi mi yoksa kablolamada bazı uyarlamalar yapmamız gerekip gerekmediğini belirlememizi sağlayacaktır.

İdeal olarak bu, modülü kurmadan _önce_ yapmanız gereken bir şeydir, ama hey, geç olması hiç olmamasından iyidir.
(Değerinde, artık istemediğinize karar verdiğiniz modüllerden kurtulmak için bir `uninstall` komutu var.)

!!! note

    CAT_CAT süreci, farklı sıkıştırma türleri, dosya uzantıları vb. ile ilgili oldukça akıllıca bir işleme içerir; bunlar burada size göstermeye çalıştığımız şeyle doğrudan ilgili değildir, bu nedenle çoğunu görmezden geleceğiz ve yalnızca önemli olan kısımlara odaklanacağız.

### 2.1. İki modülün arayüzlerini karşılaştırın

Hatırlatma olarak, `collectGreetings` modülümüzün arayüzü şöyle görünüyor:

```groovy title="modules/local/collectGreetings.nf (alıntı)" linenums="1" hl_lines="6-7 10"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
```

`collectGreetings` modülü iki girdi alır:

- `input_files` işlenecek bir veya daha fazla girdi dosyası içerir;
- `batch_name`, çıktı dosyasına çalıştırmaya özgü bir ad atamak için kullandığımız bir değerdir; bu bir metadata biçimidir.

Tamamlandığında, `collectGreetings` `outfile` etiketi ile yayınlanan tek bir dosya yolu çıktısı verir.

Karşılaştırmalı olarak, `cat/cat` modülünün arayüzü daha karmaşıktır:

```groovy title="modules/nf-core/cat/cat/main.nf (alıntı)" linenums="1" hl_lines="11 14"
process CAT_CAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"               , emit: versions
```

CAT_CAT modülü tek bir girdi alır, ancak bu girdi iki şey içeren bir demettir:

- `meta`, metadata içeren bir yapıdır, buna metamap denir;
- `files_in`, işlenecek bir veya daha fazla girdi dosyası içerir; `collectGreetings`'in `input_files` değerine eşdeğerdir.

Tamamlandığında, CAT_CAT çıktılarını iki bölümde sunar:

- Metamap ve birleştirilmiş çıktı dosyasını içeren başka bir demet, `file_out` etiketi ile yayınlanır;
- Kullanılan yazılım sürümü hakkında bilgi yakalayan bir `versions.yml` dosyası, `versions` etiketi ile yayınlanır.

Ayrıca varsayılan olarak, çıktı dosyasının metadata'nın bir parçası olan bir tanımlayıcıya göre adlandırılacağını unutmayın (kod burada gösterilmemiştir).

Sadece koda bakarak bunların hepsini takip etmek çok fazla gibi görünebilir, bu nedenle her şeyin nasıl bir araya geldiğini görselleştirmenize yardımcı olacak bir diyagram:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/module_comparison.svg"
</figure>

İki modülün içerik açısından benzer girdi gereksinimlerine sahip olduğunu (bir dizi girdi dosyası artı bazı metadata) ancak bu içeriğin nasıl paketlendiği konusunda çok farklı beklentilere sahip olduğunu görebilirsiniz.
Versions dosyasını şimdilik görmezden gelecek olursak, ana çıktıları da eşdeğerdir (birleştirilmiş bir dosya), ancak CAT_CAT ayrıca çıktı dosyasıyla birlikte metamap'i de yayınlar.

Paketleme farklılıklarıyla başa çıkmak oldukça kolay olacak, birazdan göreceğiniz gibi.
Ancak, metamap kısmını anlamak için size bazı ek bağlam sunmamız gerekiyor.

### 2.2. Metamap'leri anlamak

Size CAT_CAT modülünün girdisinin bir parçası olarak bir metadata map beklediğini söyledik.
Bunun ne olduğuna daha yakından bakmak için birkaç dakika ayıralım.

**Metadata map**, genellikle kısaca **metamap** olarak adlandırılır; veri birimleri hakkında bilgi içeren Groovy tarzı bir map'tir.
Nextflow pipeline'ları bağlamında, veri birimleri istediğiniz herhangi bir şey olabilir: bireysel örnekler, örnek grupları veya tüm veri kümeleri.

Kurala göre, bir nf-core metamap `meta` olarak adlandırılır ve çıktıları adlandırmak ve veri birimlerini izlemek için kullanılan gerekli `id` alanını içerir.

Örneğin, tipik bir metadata map şöyle görünebilir:

```groovy title="Örnek seviyesinde metamap örneği"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

Veya metadata'nın grup seviyesinde eklendiği bir durumda:

```groovy title="Grup seviyesinde metamap örneği"
[id: 'batch1', date: '25.10.01']
```

Şimdi bunu, girdi dosyalarının bir metamap ile bir demete paketlenmesini bekleyen ve metamap'i de çıktı demetinin bir parçası olarak çıktılayan `CAT_CAT` sürecinin bağlamına koyalım.

```groovy title="modules/nf-core/cat/cat/main.nf (alıntı)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

Sonuç olarak, her veri birimi ilgili metadata ekli olarak pipeline boyunca ilerler.
Sonraki süreçler de bu metadata'ya kolayca erişebilir.

`CAT_CAT` tarafından çıktılanan dosyanın metadata'nın bir parçası olan bir tanımlayıcıya göre adlandırılacağını size söylediğimizi hatırlıyor musunuz?
İlgili kod budur:

```groovy title="modules/nf-core/cat/cat/main.nf (alıntı)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

Bu kabaca şu şekilde çevrilir: eğer harici görev parametre sistemi (`task.ext`) aracılığıyla bir `prefix` sağlanırsa, çıktı dosyasını adlandırmak için onu kullan; aksi takdirde metamap'teki `id` alanına karşılık gelen `${meta.id}` kullanarak bir tane oluştur.

Bu modüle gelen girdi kanalının şöyle içeriklerle geldiğini hayal edebilirsiniz:

```groovy title="Örnek girdi kanal içeriği"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

Ardından çıktı kanal içeriği şöyle çıkar:

```groovy title="Örnek çıktı kanal içeriği"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

Daha önce belirtildiği gibi, `tuple val(meta), path(files_in)` girdi kurulumu tüm nf-core modüllerinde kullanılan standart bir kalıptır.

Umarız bunun ne kadar yararlı olabileceğini görmeye başlayabilirsiniz.
Yalnızca metadata'ya göre çıktıları adlandırmanıza izin vermekle kalmaz, aynı zamanda farklı parametre değerleri uygulamak gibi şeyler de yapabilirsiniz ve belirli operatörlerle birlikte, pipeline boyunca akan verileri gruplayabilir, sıralayabilir veya filtreleyebilirsiniz.

!!! note "Metadata hakkında daha fazla bilgi edinin"

    Nextflow iş akışlarında metadata ile çalışmaya kapsamlı bir giriş için, samplesheet'lerden metadata okumayı ve işlemeyi özelleştirmek için nasıl kullanılacağını içeren [İş akışlarında metadata](../side_quests/metadata) yan görevine bakın.

### 2.3. Yapılacak değişiklikleri özetleyin

İncelediğimiz şeylere dayanarak, `cat/cat` modülünü kullanmak için pipeline'ımızda yapmamız gereken başlıca değişiklikler şunlardır:

- Grup adını içeren bir metamap oluşturun;
- Metamap'i birleştirilecek girdi dosyaları kümesiyle (`convertToUpper`'dan çıkan) bir demete paketleyin;
- Çağrıyı `collectGreetings()`'den `CAT_CAT`'e değiştirin;
- `CAT_CAT` süreci tarafından üretilen demetten çıktı dosyasını `cowpy`'ye geçirmeden önce çıkarın.

Bu işe yaramalı! Artık bir planımız olduğuna göre, dalmaya hazırız.

### Özet

Yeni bir modülün girdi ve çıktı arayüzünü gereksinimlerini belirlemek için nasıl değerlendireceğinizi biliyorsunuz ve metamap'lerin nf-core pipeline'ları tarafından metadata'yı pipeline boyunca akan verilerle yakından ilişkili tutmak için nasıl kullanıldığını öğrendiniz.

### Sırada ne var?

Yeni modülü bir iş akışına entegre edin.

---

## 3. CAT_CAT'i `hello.nf` iş akışına entegre edin

Artık metamap'ler hakkında her şeyi bildiğinize göre (veya en azından bu kursun amaçları için yeterince), yukarıda özetlediğimiz değişiklikleri gerçekten uygulamanın zamanı geldi.

Netlik adına, bunu parçalara ayıracağız ve her adımı ayrı ayrı ele alacağız.

!!! note

    Aşağıda gösterilen tüm değişiklikler `core-hello/workflows/hello.nf` iş akışı dosyasındaki `main` bloğundaki iş akışı mantığına yapılır.

### 3.1. Bir metadata map oluşturun

İlk olarak, nf-core modüllerinin metamap'in en azından bir `id` alanı gerektirdiğini aklımızda tutarak `CAT_CAT` için bir metadata map oluşturmamız gerekiyor.

Başka bir metadata'ya ihtiyacımız olmadığı için, basit tutabilir ve şöyle bir şey kullanabiliriz:

```groovy title="Sözdizimi örneği"
def cat_meta = [id: 'test']
```

Ancak `id` değerini sabit kodlamak istemiyoruz; `params.batch` parametresinin değerini kullanmak istiyoruz.
Yani kod şöyle olur:

```groovy title="Sözdizimi örneği"
def cat_meta = [id: params.batch]
```

Evet, temel bir metamap oluşturmak kelimenin tam anlamıyla bu kadar basit.

Bu satırları `convertToUpper` çağrısından sonra ekleyelim, `collectGreetings` çağrısını kaldırarak:

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // grup adıyla ID olarak metadata map oluştur
        def cat_meta = [ id: params.batch ]

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Bu, `id`'nin grup adımıza ayarlandığı basit bir metadata map oluşturur (test profilini kullanırken `test` olacaktır).

### 3.2. Metadata demetleri içeren bir kanal oluşturun

Ardından, dosya kanalını metadata ve dosyaları içeren demet kanalına dönüştürün:

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // grup adıyla ID olarak metadata map oluştur
        def cat_meta = [ id: params.batch ]

        // demet formatında metadata ve dosyalarla bir kanal oluştur
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // grup adıyla ID olarak metadata map oluştur
        def cat_meta = [ id: params.batch ]

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Eklediğimiz satır iki şey başarır:

- `.collect()` `convertToUpper` çıktısından tüm dosyaları tek bir listeye toplar
- `.map { files -> tuple(cat_meta, files) }` `CAT_CAT`'in beklediği formatta `[metadata, files]` demeti oluşturur

`CAT_CAT` için girdi demetini kurmak için yapmamız gereken tek şey bu.

### 3.3. CAT_CAT modülünü çağırın

Şimdi yeni oluşturulan kanal üzerinde `CAT_CAT`'i çağırın:

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // grup adıyla ID olarak metadata map oluştur
        def cat_meta = [ id: params.batch ]

        // demet formatında metadata ve dosyalarla bir kanal oluştur
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // nf-core cat/cat modülünü kullanarak dosyaları birleştir
        CAT_CAT(ch_for_cat)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // grup adıyla ID olarak metadata map oluştur
        def cat_meta = [ id: params.batch ]

        // demet formatında metadata ve dosyalarla bir kanal oluştur
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Bu, bu değiştirmenin en zor kısmını tamamlar, ancak henüz tam olarak bitirmedik: birleştirilmiş çıktıyı `cowpy` sürecine nasıl geçireceğimizi hala güncellememiz gerekiyor.

### 3.4. `cowpy` için demetten çıktı dosyasını çıkarın

Daha önce, `collectGreetings` süreci doğrudan `cowpy`'ye geçirebileceğimiz bir dosya üretiyordu.
Ancak, `CAT_CAT` süreci çıktı dosyasına ek olarak metamap'i de içeren bir demet üretir.

`cowpy` henüz metadata demetlerini kabul etmediğinden (bunu kursun bir sonraki bölümünde düzelteceğiz), `CAT_CAT` tarafından üretilen demetten çıktı dosyasını `cowpy`'ye vermeden önce çıkarmamız gerekiyor:

=== "Sonra"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // grup adıyla ID olarak metadata map oluştur
        def cat_meta = [ id: params.batch ]

        // demet formatında metadata ve dosyalarla bir kanal oluştur
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenate the greetings
        CAT_CAT(ch_for_cat)

        // cowpy henüz metadata kullanmadığı için demetten dosyayı çıkar
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generate ASCII art of the greetings with cowpy
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Önce"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // grup adıyla ID olarak metadata map oluştur
        def cat_meta = [ id: params.batch ]

        // demet formatında metadata ve dosyalarla bir kanal oluştur
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenate the greetings
        CAT_CAT(ch_for_cat)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

`.map{ meta, file -> file }` işlemi, `CAT_CAT` tarafından üretilen `[metadata, file]` demetinden dosyayı yeni bir kanal olan `ch_for_cowpy`'ye çıkarır.

Ardından, son satırda `collectGreetings.out.outfile` yerine `ch_for_cowpy`'yi `cowpy`'ye geçirmek yeterlidir.

!!! note

    Kursun bir sonraki bölümünde, `cowpy`'yi doğrudan metadata demetleriyle çalışacak şekilde güncelleyeceğiz, bu nedenle bu çıkarma adımı artık gerekli olmayacaktır.

### 3.5. İş akışını test edin

Yeni entegre edilen `cat/cat` modülüyle iş akışının çalıştığını test edelim:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

Bu makul bir hızda çalışmalıdır.

??? success "Komut çıktısı"

    ```console
    N E X T F L O W ~ version 25.04.3

        Launching `./main.nf` [evil_pike] DSL2 - revision: b9e9b3b8de

        Input/output options
          input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
          outdir                    : core-hello-results

        Institutional config options
          config_profile_name       : Test profile
          config_profile_description: Minimal test dataset to check pipeline function

        Generic options
          validate_params           : false
          trace_report_suffix       : 2025-10-30_18-50-58

        Core Nextflow options
          runName                   : evil_pike
          containerEngine           : docker
          launchDir                 : /workspaces/training/hello-nf-core/core-hello
          workDir                   : /workspaces/training/hello-nf-core/core-hello/work
          projectDir                : /workspaces/training/hello-nf-core/core-hello
          userName                  : root
          profile                   : test,docker
          configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

        !! Only displaying parameters that differ from the pipeline defaults !!
        ------------------------------------------------------
        executor >  local (8)
        [b3/f005fd] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
        [08/f923d0] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
        [34/3729a9] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
        [24/df918a] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
        -[core/hello] Pipeline completed successfully-
    ```

`CAT_CAT`'in artık `collectGreetings` yerine süreç yürütme listesinde göründüğüne dikkat edin.

Ve bu kadar! Artık pipeline'daki bu adım için özel prototip düzeyindeki kod yerine sağlam, topluluk tarafından düzenlenmiş bir modül kullanıyoruz.

### Özet

Artık şunları nasıl yapacağınızı biliyorsunuz:

- nf-core modüllerini bulun ve kurun
- Bir nf-core modülünün gereksinimlerini değerlendirin
- Bir nf-core modülüyle kullanmak için basit bir metadata map oluşturun
- Bir nf-core modülünü iş akışınıza entegre edin

### Sırada ne var?

Yerel modüllerinizi nf-core kurallarını takip edecek şekilde nasıl uyarlayacağınızı öğrenin.
Ayrıca size nf-core araçlarını kullanarak bir şablondan yeni nf-core modülleri nasıl oluşturacağınızı göstereceğiz.
