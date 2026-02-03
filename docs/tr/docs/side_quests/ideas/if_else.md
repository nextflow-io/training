# Bölüm 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. İneğe ünlü bilim insanlarından alıntılar yaptırın

Bu bölüm, şimdiye kadar öğrendiklerinizi pratik yapmak için bazı ek alıştırmalar içermektedir.
Bu alıştırmaları yapmak, eğitimin sonraki bölümlerini anlamak için _gerekli değildir_, ancak ineğe ünlü bilim insanlarından alıntılar yaptırmayı çözerek öğrendiklerinizi pekiştirmenin eğlenceli bir yolunu sunar.

```console title="cowsay-output-Grace-Hopper.txt"
  _________________________________________________
 /                                                 \
| Humans are allergic to change. They love to       |
| say, 'We've always done it this way.' I try to fi |
| ght that. That's why I have a clock on my wall th |
| at runs counter-clockwise.                        |
| -Grace Hopper                                     |
 \                                                 /
  =================================================
                                                 \
                                                  \
                                                    ^__^
                                                    (oo)\_______
                                                    (__)\       )\/\
                                                        ||----w |
                                                        ||     ||
```

### 1.1. `hello-containers.nf` betiğini getQuote süreci kullanacak şekilde değiştirin

`containers/data/pioneers.csv` dosyasında bilgisayar ve biyoloji öncülerinin bir listesi bulunmaktadır.
Üst düzeyde, bu alıştırmayı tamamlamak için şunları yapmanız gerekecektir:

- Varsayılan `params.input_file` parametresini `pioneers.csv` dosyasına işaret edecek şekilde değiştirin.
- Her girdi için bir alıntı almak üzere `quote` konteynerini kullanan bir `getQuote` süreci oluşturun.
- `getQuote` sürecinin çıktısını, alıntıyı görüntülemek için `cowsay` sürecine bağlayın.

`quote` konteyner imajı için, önceki ek alıştırmada kendiniz oluşturduğunuz imajı ya da Seqera Containers'dan aldığınız imajı kullanabilirsiniz.

!!! Hint "İpucu"

    getQuote sürecinizin `script` bloğu için iyi bir seçenek şu olabilir:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Bu alıştırmanın çözümünü `containers/solutions/hello-containers-4.1.nf` dosyasında bulabilirsiniz.

### 1.2. Nextflow pipeline'ınızı `quote` ve `sayHello` modlarında çalışacak şekilde değiştirin.

Pipeline'ınıza hem `quote` hem de `sayHello` için tasarlanmış girdileri kabul edebilmesi için dallanma mantığı ekleyin.
İşte bir Nextflow workflow'unda `if` ifadesinin nasıl kullanılacağına dair bir örnek:

```groovy title="hello-containers.nf"
workflow {
    if (params.quote) {
        ...
    }
    else {
        ...
    }
    cowSay(text_ch)
}
```

!!! Hint "İpucu"

    Bir sürecin çıktı kanalına bir isim atamak için `new_ch = processName.out` kullanabilirsiniz.

Bu alıştırmanın çözümünü `containers/solutions/hello-containers-4.2.nf` dosyasında bulabilirsiniz.

### Çıkarımlar

Nextflow'da süreçleri çalıştırmak için konteynerlerin nasıl kullanılacağını ve pipeline'larınıza bazı dallanma mantıkları oluşturmayı öğrendiniz!

### Sırada ne var?

Kendinizi kutlayın, bir mola verin ve biraz su için!

Hazır olduğunuzda, şimdiye kadar öğrendiklerinizi daha gerçekçi bir veri analizi kullanım durumuna nasıl uygulayacağınızı öğrenmek için bu eğitim serisinin 3. Bölümüne geçin.
