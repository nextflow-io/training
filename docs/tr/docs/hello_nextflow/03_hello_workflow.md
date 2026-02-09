# Bölüm 3: Hello Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [Tüm oynatma listesini](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) Nextflow YouTube kanalında görün.

:green_book: Video metni [burada](./transcripts/03_hello_workflow.md) mevcuttur.
///

Gerçek dünyadaki iş akışlarının çoğu birden fazla adım içerir.
Bu eğitim modülünde, süreçleri çok adımlı bir iş akışında nasıl birbirine bağlayacağınızı öğreneceksiniz.

Bu size aşağıdakileri başarmanın Nextflow yolunu öğretecektir:

1. Verilerin bir süreçten diğerine akmasını sağlamak
2. Birden fazla süreç çağrısından gelen çıktıları tek bir süreç çağrısında toplamak
3. Bir sürece ek parametreler iletmek
4. Bir süreçten çıkan birden fazla çıktıyı işlemek

Bunu göstermek için, Bölüm 1 ve 2'deki alana özgü olmayan Hello World örneği üzerine inşa etmeye devam edeceğiz.
Bu sefer, gerçek iş akışlarını nasıl oluşturduklarını daha iyi yansıtmak için iş akışımızda aşağıdaki değişiklikleri yapacağız:

1. Selamlamayı büyük harfe dönüştüren ikinci bir adım ekleyin.
2. Tüm dönüştürülmüş selamlamaları toplayan ve bunları tek bir dosyaya yazan üçüncü bir adım ekleyin.
3. Son çıktı dosyasını adlandırmak için bir parametre ekleyin ve bunu toplama adımına ikincil bir girdi olarak iletin.
4. Toplama adımının işlenen hakkında basit bir istatistik de raporlamasını sağlayın.

??? info "Bu bölümden nasıl başlanır"

    Kursun bu bölümü, [Hello Nextflow](./index.md) kursunun Bölüm 1-2'sini tamamladığınızı varsayar, ancak bu bölümlerde ele alınan temel bilgilerde rahat iseniz, özel bir şey yapmadan buradan başlayabilirsiniz.

---

## 0. Isınma: `hello-workflow.nf` dosyasını çalıştırın

Başlangıç noktası olarak `hello-workflow.nf` iş akışı betiğini kullanacağız.
Bu eğitim kursunun Bölüm 2'sini tamamlayarak üretilen betiğe eşdeğerdir, ancak `view()` ifadelerini kaldırdık ve çıktı hedefini değiştirdik:

```groovy title="hello-workflow.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
```

Bu diyagram, iş akışının mevcut işleyişini özetlemektedir.
Tanıdık görünmeli, ancak şimdi sürecin çıktılarının girdiler gibi bir kanalda paketlendiğini açıkça gösteriyoruz.
Bu çıktı kanalını bir dakika içinde iyi bir şekilde kullanacağız.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-workflow-channels.svg"
</figure>

Her şeyin çalıştığından emin olmak için, herhangi bir değişiklik yapmadan önce betiği bir kez çalıştırın:

```bash
nextflow run hello-workflow.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [admiring_lamarr] DSL2 - revision: 4d4053520d

    executor >  local (3)
    [b1/5826b5] process > sayHello (2) [100%] 3 of 3 ✔
    ```

Daha önce olduğu gibi, çıktı dosyalarını `output` bloğunda belirtilen konumda bulacaksınız.
Bu bölüm için, `results/hello_workflow/` altındadır.

??? abstract "Dizin içeriği"

    ```console
    results/hello_workflow
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    └── Holà-output.txt
    ```

Bu sizin için çalıştıysa, çok adımlı bir iş akışını nasıl bir araya getireceğinizi öğrenmeye hazırsınız.

---

## 1. İş akışına ikinci bir adım ekleyin

Her selamlamayı büyük harfe dönüştürmek için bir adım ekleyeceğiz.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

Bunun için üç şey yapmamız gerekiyor:

- Büyük harf dönüşümünü yapmak için kullanacağımız komutu tanımlayın.
- Büyük harf dönüştürme komutunu sarmalayan yeni bir süreç yazın.
- Yeni süreci workflow bloğunda çağırın ve `sayHello()` sürecinin çıktısını girdi olarak alacak şekilde ayarlayın.

### 1.1. Büyük harf dönüştürme komutunu tanımlayın ve terminalde test edin

Selamlamaları büyük harfe dönüştürmek için, 'text replacement' (metin değiştirme) anlamına gelen `tr` adlı klasik bir UNIX aracını aşağıdaki sözdizimi ile kullanacağız:

```bash title="Sözdizimi"
tr '[a-z]' '[A-Z]'
```

Bu, aksanlı harfleri hesaba katmayan çok basit bir metin değiştirme tek satırıdır, bu nedenle örneğin 'Holà' 'HOLà' olacaktır, ancak Nextflow kavramlarını göstermek için yeterince iyi bir iş çıkaracaktır ve önemli olan budur.

Test etmek için, `echo 'Hello World'` komutunu çalıştırabilir ve çıktısını `tr` komutuna yönlendirebiliriz:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

Çıktı, `Hello World` dizesinin büyük harf versiyonunu içeren `UPPER-output.txt` adlı bir metin dosyasıdır.

??? abstract "Dosya içeriği"

    ```console title="UPPER-output.txt"
    HELLO WORLD
    ```

Temel olarak iş akışımızla yapmaya çalışacağımız şey budur.

### 1.2. Büyük harf dönüştürme adımını Nextflow süreci olarak yazın

Yeni sürecimizi ilkine göre modelleyebiliriz, çünkü aynı bileşenlerin hepsini kullanmak istiyoruz.

Aşağıdaki süreç tanımını, ilkinin hemen altına iş akışı betiğine ekleyin:

```groovy title="hello-workflow.nf" linenums="20"
/*
 * Selamlamayı büyük harfe dönüştürmek için bir metin değiştirme aracı kullanın
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

Bunda, ilk sürecin çıktısı için yaptığımıza benzer şekilde, ikinci çıktı dosya adını girdi dosya adına göre oluşturuyoruz.

### 1.3. Workflow bloğuna yeni sürece bir çağrı ekleyin

Şimdi Nextflow'a az önce tanımladığımız süreci gerçekten çağırmasını söylememiz gerekiyor.

Workflow bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="44" hl_lines="10-11"
    workflow {

        main:
        // CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // bir selamlama yayınla
        sayHello(greeting_ch)
        // selamlamayı büyük harfe dönüştür
        convertToUpper()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="44"
    workflow {

        main:
        // CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Bu henüz işlevsel değil çünkü `convertToUpper()` sürecine ne girdi yapılması gerektiğini belirtmedik.

### 1.4. İlk sürecin çıktısını ikinci sürece iletin

Şimdi `sayHello()` sürecinin çıktısının `convertToUpper()` sürecine akmasını sağlamamız gerekiyor.

Uygun bir şekilde, Nextflow bir sürecin çıktısını otomatik olarak bir kanala paketler, ısınma bölümündeki diyagramda gösterildiği gibi.
Bir sürecin çıktı kanalına `<process>.out` olarak başvurabiliriz.

Dolayısıyla `sayHello` sürecinin çıktısı, doğrudan `convertToUpper()` içine takabileceğimiz `sayHello.out` adlı bir kanaldır.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep-connector.svg"
</figure>

Workflow bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // selamlamayı büyük harfe dönüştür
        convertToUpper()
    ```

Bunun gibi basit bir durum için (bir çıktıdan bir girdiye), iki süreci bağlamak için yapmamız gereken tek şey bu!

### 1.5. İş akışı çıktı yayınlamasını ayarlayın

Son olarak, ikinci süreçten gelen sonuçları da yayınlamak için iş akışı çıktılarını güncelleyelim.

#### 1.5.1. `workflow` bloğunun `publish:` bölümünü güncelleyin

Workflow bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="56" hl_lines="3"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
    }
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="56"
        publish:
        first_output = sayHello.out
    }
    ```

Mantık daha öncekiyle aynıdır.

#### 1.5.2. `output` bloğunu güncelleyin

`output` bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="61" hl_lines="6-9"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="61"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

Bir kez daha, mantık daha öncekiyle aynıdır.

Bu size çıktı ayarlarını her bir çıktı için çok ayrıntılı bir düzeyde kontrol edebileceğinizi gösterir.
Ne olduğunu görmek için süreçlerden biri için yolları veya yayınlama modunu değiştirmeyi deneyin.

Tabii ki, bu burada bazı bilgileri tekrarladığımız anlamına gelir, bu da tüm çıktılar için konumu aynı şekilde güncellemek istesek rahatsız edici olabilir.
Kursun ilerleyen bölümlerinde, bu ayarları yapılandırılmış bir şekilde birden fazla çıktı için nasıl yapılandıracağınızı öğreneceksiniz.

### 1.6. İş akışını `-resume` ile çalıştırın

Bunu `-resume` bayrağını kullanarak test edelim, çünkü iş akışının ilk adımını zaten başarıyla çalıştırdık.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [high_cantor] DSL2 - revision: d746983511

    executor >  local (3)
    [ab/816321] process > sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [e0/ecf81b] process > convertToUpper (3) [100%] 3 of 3 ✔
    ```

Konsol çıktısında şimdi az önce eklediğimiz yeni sürece karşılık gelen ekstra bir satır var.

Çıktıları `output` bloğunda ayarlandığı gibi `results/hello_workflow` dizininde bulacaksınız.

??? abstract "Dizin içeriği"

    ```console
    results/hello_workflow/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Bu kullanışlı! Ancak ikinci sürece yapılan çağrılardan birinin çalışma dizinine bakmaya değer.

??? abstract "Dizin içeriği"

    ```console
    work/e0/ecf81b4cacc648b9b994218d5b29d7/
    ├── Holà-output.txt -> /workspaces/training/hello-nextflow/work/ab/81632178cd37e9e815959278808819/Holà-output.txt
    └── UPPER-Holà-output.txt
    ```

İki `*-output` dosyası olduğuna dikkat edin: ilk sürecin çıktısı ve ikinci sürecin çıktısı.

İlk sürecin çıktısı orada çünkü Nextflow, yürütme için gereken her şeyi aynı alt dizinde bulundurmak için orada **sahneledi**.

Ancak, aslında ilk süreç çağrısının alt dizinindeki orijinal dosyaya işaret eden sembolik bir bağlantıdır.
Varsayılan olarak, burada yaptığımız gibi tek bir makinede çalışırken, Nextflow girdi ve ara dosyaları sahnelemek için kopyalar yerine sembolik bağlantılar kullanır.

Şimdi, devam etmeden önce, tek yaptığımız şeyin `sayHello` çıktısını `convertToUpper` girdisine bağlamak olduğunu ve iki sürecin seri olarak çalıştırılabildiğini düşünün.
Nextflow, bireysel girdi ve çıktı dosyalarını işleme ve bunları iki komut arasında iletme konusundaki zor işi bizim için yaptı.

Bu, Nextflow kanallarının bu kadar güçlü olmasının nedenlerinden biridir: iş akışı adımlarını birbirine bağlamada yer alan meşakkatli işleri hallederler.

### Özet

Bir adımın çıktısını bir sonraki adıma girdi olarak sağlayarak süreçleri nasıl zincirleme yapacağınızı biliyorsunuz.

### Sırada ne var?

Toplu süreç çağrılarından gelen çıktıları nasıl toplayacağınızı ve bunları tek bir sürece nasıl besleyeceğinizi öğrenin.

---

## 2. Tüm selamlamaları toplamak için üçüncü bir adım ekleyin

Burada yaptığımız gibi, bir kanaldaki öğelerin her birine bir dönüşüm uygulamak için bir süreç kullandığımızda, bazen o sürecin çıktı kanalından öğeleri toplamak ve bunları bir tür analiz veya toplama gerçekleştiren başka bir sürece beslemek isteriz.

Bunu göstermek için, `convertToUpper` süreci tarafından üretilen tüm büyük harf selamlamaları toplayan ve bunları tek bir dosyaya yazan boru hattımıza yeni bir adım ekleyeceğiz.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

Sürprizi bozmak istemem ama bu çok kullanışlı bir operatör içerecek.

### 2.1. Toplama komutunu tanımlayın ve terminalde test edin

İş akışımıza eklemek istediğimiz toplama adımı, birden fazla büyük harf selamlamayı tek bir dosyada birleştirmek için `cat` komutunu kullanacaktır.

Daha önce yaptığımız gibi, beklendiği gibi çalıştığını doğrulamak için komutu terminalde tek başına çalıştıralım.

Terminalinizde aşağıdakileri çalıştırın:

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

Çıktı, orijinal selamlamaların büyük harf versiyonlarını içeren `COLLECTED-output.txt` adlı bir metin dosyasıdır.

??? abstract "Dosya içeriği"

    ```console title="COLLECTED-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

İş akışımızla elde etmek istediğimiz sonuç budur.

### 2.2. Toplama adımını yapmak için yeni bir süreç oluşturun

Yeni bir süreç oluşturalım ve ona `collectGreetings()` diyelim.
Daha önce gördüklerimize dayanarak yazmaya başlayabiliriz.

#### 2.2.1. Sürecin 'açık' kısımlarını yazın

Aşağıdaki süreç tanımını iş akışı betiğine ekleyin:

```groovy title="hello-workflow.nf" linenums="37"
/*
 * Büyük harf selamlamaları tek bir çıktı dosyasında toplayın
 */
process collectGreetings {

    input:
    ???

    output:
    path "COLLECTED-output.txt"

    script:
    """
    cat ??? > 'COLLECTED-output.txt'
    """
}
```

Bu, şimdiye kadar öğrendiklerinize dayanarak güvenle yazabileceğimiz şeydir.
Ancak bu işlevsel değil!
Girdi tanımını ve betik komutunun ilk yarısını dışarıda bırakıyor çünkü bunu nasıl yazacağımızı bulmamız gerekiyor.

#### 2.2.2. `collectGreetings()` için girdileri tanımlayın

`convertToUpper()` sürecinden gelen tüm çağrılardan selamlamaları toplamamız gerekiyor.
Önceki adımdan ne alabileceğimizi biliyoruz?

`convertToUpper()` tarafından çıktılanan kanal, büyük harf selamlamaları içeren bireysel dosyaların yollarını içerecektir.
Bu bir girdi yuvasına denk gelir; basitlik için buna `input_files` diyelim.

Süreç bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          path input_files
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          ???
    ```

Bunun birden fazla dosya içermesini beklememize rağmen `path` önekini kullandığımıza dikkat edin.

#### 2.2.3. Birleştirme komutunu oluşturun

Burada işler biraz karmaşık olabilir, çünkü keyfi sayıda girdi dosyasını işleyebilmemiz gerekiyor.
Özellikle, komutu önceden yazamayız, bu nedenle Nextflow'a çalışma zamanında sürece akan girdilere göre nasıl oluşturacağını söylememiz gerekiyor.

Başka bir deyişle, `[file1.txt, file2.txt, file3.txt]` öğesini içeren bir girdi kanalımız varsa, Nextflow'un bunu `cat file1.txt file2.txt file3.txt` haline getirmesi gerekir.

Neyse ki, betik komutunda basitçe `cat ${input_files}` yazarsak Nextflow bunu bizim için yapmaktan mutluluk duyar.

Süreç bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="54" hl_lines="3"
        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="54"
        script:
        """
        cat ??? > 'COLLECTED-output.txt'
        """
    ```

Teoride bu, keyfi sayıda girdi dosyasını işlemelidir.

!!! tip

    Bazı komut satırı araçları, her girdi dosyası için bir argüman (örneğin `-input`) sağlamayı gerektirir.
    Bu durumda, komutu oluşturmak için biraz ekstra iş yapmamız gerekir.
    Bunun bir örneğini [Nextflow for Genomics](../../nf4_science/genomics/) eğitim kursunda görebilirsiniz.

### 2.3. Toplama adımını iş akışına ekleyin

Şimdi sadece büyük harf dönüştürme adımının çıktısında toplama sürecini çağırmamız gerekiyor.
Bu da `convertToUpper.out` adlı bir kanaldır.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-connector.svg"
</figure>

#### 2.3.1. Süreç çağrılarını bağlayın

Workflow bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)

        // tüm selamlamaları bir dosyada topla
        collectGreetings(convertToUpper.out)
    }
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="75"
        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)
    }
    ```

Bu, `convertToUpper()` çıktısını `collectGreetings()` girdisine bağlar.

#### 2.3.2. İş akışını `-resume` ile çalıştırın

Hadi deneyelim.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Komut çıktısı"

    ```console hl_lines="8"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [mad_gilbert] DSL2 - revision: 6acfd5e28d

    executor >  local (3)
    [79/33b2f0] sayHello (2)         | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (3)   | 3 of 3, cached: 3 ✔
    [47/50fe4a] collectGreetings (1) | 3 of 3 ✔
    ```

Üçüncü adım dahil olmak üzere başarıyla çalışır.

Ancak, son satırdaki `collectGreetings()` için çağrı sayısına bakın.
Sadece bir tane bekliyorduk, ancak üç tane var.

Şimdi son çıktı dosyasının içeriğine bir göz atın.

??? abstract "Dosya içeriği"

    ```console title="results/COLLECTED-output.txt"
    Holà
    ```

Ah hayır. Toplama adımı her selamlama üzerinde ayrı ayrı çalıştırıldı, bu istediğimiz şey DEĞİLDİ.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-no-operator.svg"
</figure>

Nextflow'a üçüncü adımın `convertToUpper()` tarafından çıktılanan kanaldaki tüm öğeler üzerinde çalışmasını istediğimizi açıkça söylemek için bir şeyler yapmamız gerekiyor.

### 2.4. Selamlamaları tek bir girdide toplamak için bir operatör kullanın

Evet, bir kez daha sorunumuzun cevabı bir operatördür.

Özellikle, uygun şekilde adlandırılmış [`collect()`](https://nextflow.io/docs/latest/reference/operator.html#collect) operatörünü kullanacağız.

#### 2.4.1. `collect()` operatörünü ekleyin

Bu sefer biraz farklı görünecek çünkü operatörü bir kanal fabrikası bağlamında eklemiyoruz; onu bir çıktı kanalına ekliyoruz.

`convertToUpper.out` alıyoruz ve `collect()` operatörünü ekliyoruz, bu bize `convertToUpper.out.collect()` veriyor.
Bunu doğrudan `collectGreetings()` süreç çağrısına takabiliriz.

Workflow bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // tüm selamlamaları bir dosyada topla
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // tüm selamlamaları bir dosyada topla
        collectGreetings(convertToUpper.out)
    }
    ```

#### 2.4.2. Bazı `view()` ifadeleri ekleyin

Kanal içeriklerinin öncesi ve sonrası durumlarını görselleştirmek için birkaç `view()` ifadesi de ekleyelim.

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // tüm selamlamaları bir dosyada topla
        collectGreetings(convertToUpper.out.collect())

        // isteğe bağlı view ifadeleri
        convertToUpper.out.view { contents -> "Collect öncesi: $contents" }
        convertToUpper.out.collect().view { contents -> "Collect sonrası: $contents" }
    }
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="73"
        // tüm selamlamaları bir dosyada topla
        collectGreetings(convertToUpper.out.collect())
    }
    ```

`view()` ifadeleri istediğiniz yere gidebilir; okunabilirlik için çağrıdan hemen sonra koyduk.

#### 2.4.3. İş akışını tekrar `-resume` ile çalıştırın

Hadi deneyelim:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    Collect öncesi: /workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt
    Collect öncesi: /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt
    Collect öncesi: /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt
    Collect sonrası: [/workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt, /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt, /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt]
    ```

Başarıyla çalışır, ancak log çıktısı bundan biraz daha dağınık görünebilir (okunabilirlik için temizledik).

Bu sefer üçüncü adım sadece bir kez çağrıldı!
`view()` ifadelerinin çıktısına bakıldığında, aşağıdakileri görüyoruz:

- Her selamlama için bir tane olmak üzere üç `Collect öncesi:` ifadesi: o noktada dosya yolları kanaldaki bireysel öğelerdir.
- Tek bir `Collect sonrası:` ifadesi: üç dosya yolu artık tek bir öğede paketlenmiştir.

Bunu aşağıdaki diyagramla özetleyebiliriz:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-WITH-operator.svg"
</figure>

Son olarak, her şeyin doğru çalıştığından emin olmak için çıktı dosyasının içeriğine bir göz atabilirsiniz.

??? abstract "Dosya içeriği"

    ```console title="results/COLLECTED-output.txt"
    BONJOUR
    HELLO
    HOLà
    ```

Bu sefer son çıktı dosyasında üç selamlamanın hepsi var. Başarı!

!!! note

    Bunu `-resume` olmadan birkaç kez çalıştırırsanız, selamlamaların sırasının bir çalıştırmadan diğerine değiştiğini göreceksiniz.
    Bu size, öğelerin süreç çağrıları boyunca aktığı sıranın tutarlı olmasının garanti edilmediğini gösterir.

#### 2.4.4. Okunabilirlik için `view()` ifadelerini kaldırın

Bir sonraki bölüme geçmeden önce, konsol çıktısını karmaşıklaştırmamak için `view()` ifadelerini silmenizi öneririz.

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="73"
        // tüm selamlamaları bir dosyada topla
        collectGreetings(convertToUpper.out.collect())
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // tüm selamlamaları bir dosyada topla
        collectGreetings(convertToUpper.out.collect())

        // isteğe bağlı view ifadeleri
        convertToUpper.out.view { contents -> "Collect öncesi: $contents" }
        convertToUpper.out.collect().view { contents -> "Collect sonrası: $contents" }
    ```

Bu temelde 2.4.2 noktasından ters işlemdir.

### Özet

Bir grup süreç çağrısından gelen çıktıları nasıl toplayacağınızı ve bunları ortak bir analiz veya toplama adımına nasıl besleyeceğinizi biliyorsunuz.

Özetlemek gerekirse, şimdiye kadar oluşturduğunuz şey budur:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

### Sırada ne var?

Bir sürece birden fazla girdi nasıl iletileceğini öğrenin.

---

## 3. Bir sürece ek parametreler iletin

Sonraki selamlama gruplarını işlerken son sonuçların üzerine yazmadan son çıktı dosyasını belirli bir şey olarak adlandırabilmek istiyoruz.

Bu amaçla, iş akışında aşağıdaki iyileştirmeleri yapacağız:

- Toplayıcı süreci, çıktı dosyası için kullanıcı tanımlı bir ad (`batch_name`) kabul edecek şekilde değiştirin
- İş akışına bir komut satırı parametresi (`--batch`) ekleyin ve bunu toplayıcı sürece iletin

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-batch.svg"
</figure>

### 3.1. Toplayıcı süreci değiştirin

Ek girdiyi bildirmemiz ve bunu çıktı dosya adına entegre etmemiz gerekecek.

#### 3.1.1. Ek girdiyi bildirin

İyi haber: süreç tanımında istediğimiz kadar girdi değişkeni bildirebiliriz.
Buna `batch_name` diyelim.

Süreç bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="3"
        input:
        path input_files
        val batch_name
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="42"
        input:
        path input_files
    ```

Süreçlerinizi istediğiniz kadar girdi bekleyecek şekilde ayarlayabilirsiniz.
Şu anda bunların hepsi gerekli girdiler olarak ayarlanmıştır; iş akışının çalışması için bir değer sağlamanız GEREKİR.

Gerekli ve isteğe bağlı girdileri nasıl yöneteceğinizi Nextflow yolculuğunuzda daha sonra öğreneceksiniz.

#### 3.1.2. Çıktı dosya adında `batch_name` değişkenini kullanın

Değişkeni, daha önce dinamik dosya adları oluşturduğumuz şekilde çıktı dosya adına ekleyebiliriz.

Süreç bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-${batch_name}-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

Bu, süreci iş akışının son çıktısı için belirli bir dosya adı oluşturmak üzere `batch_name` değerini kullanacak şekilde ayarlar.

### 3.2. Bir `batch` komut satırı parametresi ekleyin

Şimdi `batch_name` için değeri sağlamanın ve bunu süreç çağrısına beslemenin bir yoluna ihtiyacımız var.

#### 3.2.1. Parametreyi ayarlamak için `params` kullanın

CLI parametrelerini bildirmek için `params` sistemini nasıl kullanacağınızı zaten biliyorsunuz.
Bunu bir `batch` parametresi bildirmek için kullanalım (varsayılan bir değerle çünkü tembeliz).

Pipeline parametreleri bölümünde, aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="6"
    /*
     * Pipeline parametreleri
     */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="55"
    /*
     * Pipeline parametreleri
     */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

`--input` için gösterdiğimiz gibi, komut satırında `--batch` ile bir değer belirterek bu varsayılan değeri geçersiz kılabilirsiniz.

#### 3.2.2. `batch` parametresini sürece iletin

Parametrenin değerini sürece sağlamak için, onu süreç çağrısına eklememiz gerekiyor.

Workflow bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // tüm selamlamaları bir dosyada topla
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // tüm selamlamaları bir dosyada topla
        collectGreetings(convertToUpper.out.collect())
    ```

Bir sürece birden fazla girdi sağlamak için, bunları çağrı parantezlerinde virgülle ayrılmış olarak listelediğinizi görüyorsunuz.

!!! warning

    Girdileri sürece, sürecin girdi tanım bloğunda listelendikleri TAMAMEN AYNI SIRADA sağlamanız GEREKİR.

### 3.3. İş akışını çalıştırın

Bunu komut satırında bir grup adıyla çalıştırmayı deneyelim.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [confident_rutherford] DSL2 - revision: bc58af409c

    executor >  local (1)
    [79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [b5/f19efe] collectGreetings   | 1 of 1 ✔
    ```

Başarıyla çalışır ve istenen çıktıyı üretir:

??? abstract "Dosya içeriği"

    ```console title="results/COLLECTED-trio-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Şimdi, parametreyi uygun şekilde belirttiğimiz sürece, diğer girdi grupları üzerindeki sonraki çalıştırmalar önceki sonuçları bozmayacaktır.

### Özet

Bir sürece birden fazla girdi nasıl iletileceğini biliyorsunuz.

### Sırada ne var?

Birden fazla çıktı nasıl yayınlanacağını ve bunların uygun şekilde nasıl işleneceğini öğrenin.

---

## 4. Toplayıcı adıma bir çıktı ekleyin

Şimdiye kadar her biri yalnızca bir çıktı üreten süreçler kullanıyorduk.
İlgili çıktılarına, hem bir çıktıyı bir sonraki sürece iletme bağlamında (örn. `convertToUpper(sayHello.out)`) hem de `publish:` bölümü bağlamında (örn. `first_output = sayHello.out`) kullandığımız `<process>.out` sözdizimini kullanarak çok rahat bir şekilde erişebildik.

Bir süreç birden fazla ürettiğinde ne olur?
Birden fazla çıktıyı nasıl işleriz?
Belirli bir çıktıyı seçip kullanabilir miyiz?

Hepsi mükemmel sorular ve kısa cevap evet yapabiliriz!

Birden fazla çıktı ayrı kanallara paketlenecektir.
Ya bu çıktı kanallarına isim vermeyi seçebiliriz, bu da daha sonra onlara ayrı ayrı başvurmayı kolaylaştırır, ya da onlara indeks ile başvurabiliriz.

Gösteri amaçlı olarak, belirli bir girdi grubu için toplanan selamlama sayısını saymak ve bunu bir dosyada raporlamak istediğimizi varsayalım.

### 4.1. Selamlama sayısını saymak ve çıktılamak için süreci değiştirin

Bu, süreç tanımında iki temel değişiklik gerektirecektir: selamlamaları saymanın ve bir rapor dosyası yazmanın bir yoluna ihtiyacımız var, sonra bu rapor dosyasını sürecin `output` bloğuna eklememiz gerekiyor.

#### 4.1.1. Toplanan selamlama sayısını sayın

Uygun bir şekilde, Nextflow süreç tanımının `script:` bloğuna keyfi kod eklememize izin verir, bu da bunun gibi şeyler yapmak için gerçekten kullanışlıdır.

Bu, `input_files` dizisindeki dosya sayısını almak için Nextflow'un yerleşik `size()` fonksiyonunu kullanabileceğimiz ve sonucu bir `echo` komutuyla dosyaya yazabileceğimiz anlamına gelir.

`collectGreetings` süreç bloğunda, aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'Bu grupta ${count_greetings} selamlama vardı.' > '${batch_name}-report.txt'
        """
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="55"
        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

`count_greetings` değişkeni çalışma zamanında hesaplanacaktır.

#### 4.1.2. Rapor dosyasını yayınlayın ve çıktıları adlandırın

Prensipte tek yapmamız gereken rapor dosyasını `output:` bloğuna eklemektir.

Ancak, bunu yaparken, çıktı bildirimlerimize bazı `emit:` etiketleri de ekleyeceğiz. Bunlar, konumsal indeksler kullanmak zorunda kalmak yerine çıktıları ada göre seçmemizi sağlayacaktır.

Süreç bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 3"
        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="46"
        output:
        path "COLLECTED-${batch_name}-output.txt"
    ```

`emit:` etiketleri isteğe bağlıdır ve sadece çıktılardan birine etiket ekleyebilirdik.
Ancak atasözünde dendiği gibi, neden ikisi birden olmasın?

!!! tip

    Bir sürecin çıktılarını `emit:` kullanarak adlandırmazsanız, yine de ilgili (sıfır tabanlı) indekslerini kullanarak onlara ayrı ayrı erişebilirsiniz.
    Örneğin, ilk çıktıyı almak için `<process>.out[0]`, ikinci çıktıyı almak için `<process>.out[1]` vb. kullanırsınız.

    Çıktıları adlandırmayı tercih ediyoruz çünkü aksi takdirde, özellikle süreç çok sayıda çıktı ürettiğinde, yanlış indeksi hata ile almak çok kolaydır.

### 4.2. İş akışı çıktılarını güncelleyin

Artık `collectGreetings` sürecinden iki çıktı geldiğine göre, `collectGreetings.out` çıktısı iki kanal içerir:

- `collectGreetings.out.outfile` son çıktı dosyasını içerir
- `collectGreetings.out.report` rapor dosyasını içerir

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-report.svg"
</figure>

İş akışı çıktılarını buna göre güncellememiz gerekiyor.

#### 4.2.1. `publish:` bölümünü güncelleyin

`workflow` bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4 5"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out
    ```

Gördüğünüz gibi, belirli süreç çıktılarına başvurmak artık önemsizdir.
Bölüm 5'te (Konteynerler) boru hattımıza bir adım daha eklediğimizde, `collectGreetings.out.outfile` dosyasına kolayca başvurabilecek ve onu yeni sürece verebileceğiz (spoiler: yeni süreç `cowpy` olarak adlandırılır).

Ancak şimdilik, iş akışı düzeyindeki çıktıları güncellemeyi bitirelim.

#### 4.2.2. `output` bloğunu güncelleyin

`output` bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-workflow.nf" linenums="86" hl_lines="14-17"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
        batch_report {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Önce"

    ```groovy title="hello-workflow.nf" linenums="80"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

`collected` çıktı tanımını güncellememize gerek yok çünkü bu ad değişmedi.
Sadece yeni çıktıyı eklememiz gerekiyor.

### 4.3. İş akışını çalıştırın

Bunu mevcut selamlama grubuyla çalıştırmayı deneyelim.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [ecstatic_wilson] DSL2 - revision: c80285f8c8

    executor >  local (1)
    [c5/4c6ca9] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [0e/6cbc59] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [02/61ead2] collectGreetings   [100%] 1 of 1 ✔
    ```

`results/hello_workflow/` dizinine bakarsanız, yeni rapor dosyasını, `trio-report.txt` dosyasını bulacaksınız.
İş akışının işlenen selamlama sayısını doğru bir şekilde raporladığını doğrulamak için açın.

??? abstract "Dosya içeriği"

    ```txt title="trio-report.txt"
    Bu grupta 3 selamlama vardı.
    ```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-4-way.svg"
</figure>

CSV'ye daha fazla selamlama eklemekten ve ne olduğunu test etmekten çekinmeyin.

### Özet

Bir sürecin birden fazla adlandırılmış çıktı nasıl yayınlayacağını ve bunların iş akışı düzeyinde uygun şekilde nasıl işleneceğini biliyorsunuz.

Daha genel olarak, süreçleri yaygın şekillerde birbirine bağlamada yer alan temel ilkeleri anlıyorsunuz.

### Sırada ne var?

Ekstra uzun bir mola verin, bunu hak ettiniz.

Hazır olduğunuzda, daha iyi sürdürülebilirlik ve kod verimliliği için kodunuzu nasıl modülerleştireceğinizi öğrenmek için [**Bölüm 4: Hello Modules**](./04_hello_modules.md) bölümüne geçin.

---

## Quiz

<quiz>
Workflow bloğunda bir sürecin çıktısına nasıl erişirsiniz?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

Daha fazla bilgi: [1.4. İlk sürecin çıktısını ikinci sürece iletin](#14-pass-the-output-of-the-first-process-to-the-second-process)
</quiz>

<quiz>
Nextflow'da süreç yürütme sırasını ne belirler?
- [ ] Süreçlerin workflow bloğunda yazıldığı sıra
- [ ] Süreç adına göre alfabetik sıra
- [x] Süreçler arasındaki veri bağımlılıkları
- [ ] Paralel yürütme için rastgele sıra

Daha fazla bilgi: [1.4. İlk sürecin çıktısını ikinci sürece iletin](#14-pass-the-output-of-the-first-process-to-the-second-process)
</quiz>

<quiz>
Tüm çıktıları alt süreç için tek bir listeye toplamak için `???` yerine hangi operatör gelmelidir?

```groovy hl_lines="4"
workflow {
    greetings_ch = Channel.of('Hello', 'Bonjour', 'Hola')
    SAYHELLO(greetings_ch)
    GATHER_ALL(SAYHELLO.out.???)
}
```

- [ ] `flatten()`
- [x] `collect()`
- [ ] `mix()`
- [ ] `join()`

Daha fazla bilgi: [2.4. Selamlamaları tek bir girdide toplamak için bir operatör kullanın](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
`collect()` operatörünü ne zaman kullanmalısınız?
- [ ] Öğeleri paralel olarak işlemek istediğinizde
- [ ] Kanal içeriklerini filtrelemeniz gerektiğinde
- [x] Alt sürecin üst süreçten gelen tüm öğelere ihtiyacı olduğunda
- [ ] Verileri birden fazla sürece bölmek istediğinizde

Daha fazla bilgi: [2.4. Selamlamaları tek bir girdide toplamak için bir operatör kullanın](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
Bir süreçten adlandırılmış bir çıktıya nasıl erişirsiniz?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

Daha fazla bilgi: [4.1.2. Rapor dosyasını yayınlayın ve çıktıları adlandırın](#412-emit-the-report-file-and-name-outputs)
</quiz>

<quiz>
Bir süreçte bir çıktıyı adlandırmak için doğru sözdizimi nedir?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

Daha fazla bilgi: [4.1.2. Rapor dosyasını yayınlayın ve çıktıları adlandırın](#412-emit-the-report-file-and-name-outputs)
</quiz>

<quiz>
Bir sürece birden fazla girdi sağlarken, ne doğru olmalıdır?
- [ ] Tüm girdiler aynı türde olmalıdır
- [ ] Girdiler alfabetik sırayla sağlanmalıdır
- [x] Girdilerin sırası, girdi bloğunda tanımlanan sırayla eşleşmelidir
- [ ] Aynı anda yalnızca iki girdi sağlanabilir

Daha fazla bilgi: [3. Bir sürece birden fazla girdi iletin](#3-pass-more-than-one-input-to-a-process)
</quiz>
