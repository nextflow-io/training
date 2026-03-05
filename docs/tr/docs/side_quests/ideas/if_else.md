# Bölüm 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. İneği ünlü bilim insanlarından alıntı yaptırın

Bu bölüm, şimdiye kadar öğrendiklerinizi pratik yapmak için bazı ek alıştırmalar içermektedir.
Bu alıştırmaları yapmak, eğitimin sonraki bölümlerini anlamak için _gerekli değildir_, ancak ineği ünlü bilim insanlarından alıntı yaptırmayı öğrenerek öğrendiklerinizi pekiştirmenin eğlenceli bir yolunu sunar.

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
- Her girdi için bir alıntı getirmek üzere `quote` konteynerini kullanan bir `getQuote` süreci oluşturun.
- `getQuote` sürecinin çıktısını, alıntıyı görüntülemek için `cowsay` sürecine bağlayın.

`quote` konteyner imajı için, önceki ek alıştırmada kendiniz oluşturduğunuz imajı veya Seqera Containers'dan aldığınız imajı kullanabilirsiniz.

!!! Hint "İpucu"

    getQuote sürecinizin `script` bloğu için iyi bir seçim şu olabilir:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Bu alıştırmanın çözümünü `containers/solutions/hello-containers-4.1.nf` dosyasında bulabilirsiniz.

### 1.2. Nextflow pipeline'ınızı `quote` ve `sayHello` modlarında çalışabilecek şekilde değiştirin.

Pipeline'ınıza hem `quote` hem de `sayHello` için tasarlanmış girdileri kabul edebilmesi için dallanma mantığı ekleyin.
İşte bir Nextflow iş akışında `if` ifadesinin nasıl kullanılacağına dair bir örnek:

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

    Bir sürecin çıktı kanalına isim atamak için `new_ch = processName.out` kullanabilirsiniz.

Bu alıştırmanın çözümünü `containers/solutions/hello-containers-4.2.nf` dosyasında bulabilirsiniz.

### Özet

Nextflow'da süreçleri çalıştırmak için konteynerlerin nasıl kullanılacağını ve pipeline'larınıza bazı dallanma mantığının nasıl ekleneceğini öğrendiniz!

### Sırada ne var?

Kutlama yapın, bir esneme molası verin ve biraz su için!

Hazır olduğunuzda, şimdiye kadar öğrendiklerinizi daha gerçekçi bir veri analizi kullanım senaryosuna nasıl uygulayacağınızı öğrenmek için bu eğitim serisinin Bölüm 3'üne geçin.
