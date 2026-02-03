# Bölüm 1: Daha Fazla Konteyner

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[YAPILACAK]

---

## 1. Konteyner imajlarını nasıl bulabilir veya oluşturabilirsiniz

Bazı yazılım geliştiricileri, yazılımları için Docker Hub gibi konteyner kayıt defterlerinde mevcut olan konteyner imajları sağlar, ancak birçoğu sağlamaz.
Bu isteğe bağlı bölümde, Nextflow pipeline'larınızda kullanmak istediğiniz araçlar için konteyner imajı elde etmenin iki yolunu göstereceğiz: Seqera Containers'ı kullanmak ve konteyner imajını kendiniz oluşturmak.

Bu bölümün sonundaki alıştırmada kullanılacak olan `quote` pip paketi için bir konteyner imajı elde edecek/oluşturacaksınız.

### 1.1. Seqera Containers'dan konteyner imajı alın

Seqera Containers, pip ve conda (bioconda dahil) ile yüklenebilen araçlar için konteyner imajları oluşturan ücretsiz bir hizmettir.
[Seqera Containers](https://www.seqera.io/containers/) sitesine gidin ve `quote` pip paketini arayın.

![Seqera Containers](img/seqera-containers-1.png)

`quote` pip paketi için bir konteyner imajı talep etmek üzere "+Add" ve ardından "Get Container" düğmesine tıklayın.

![Seqera Containers](img/seqera-containers-2.png)

Paketin bu sürümü için bir topluluk konteyneri ilk kez oluşturuluyorsa, tamamlanması birkaç dakika sürebilir.
Sizin için oluşturulan konteyner imajının URI'sini (örneğin `community.wave.seqera.io/library/pip_quote:ae07804021465ee9`) kopyalamak için tıklayın.

Artık konteyner imajını kullanarak `quote` komutunu çalıştırabilir ve Grace Hopper'dan rastgele bir söz alabilirsiniz.

```bash
docker run --rm community.wave.seqera.io/library/pip_quote:ae07804021465ee9 quote "Grace Hopper"
```

Çıktı:

```console title="Çıktı"
Humans are allergic to change. They love to say, 'We've always done it
this way.' I try to fight that. That's why I have a clock on my wall
that runs counter-clockwise.
```

### 1.2. Konteyner imajını kendiniz oluşturun

Seqera Containers web sitesindeki bazı derleme ayrıntılarını kullanarak `quote` pip paketi için konteyner imajını kendimiz oluşturalım.
Seqera Containers web sitesine dönün ve "Build Details" düğmesine tıklayın.

Bakacağımız ilk öğe, konteyner imajını oluşturmak için gereken tüm komutları içeren bir tür betik dosyası olan `Dockerfile`'dır.
Her bir bölümün ne yaptığını anlamanıza yardımcı olmak için aşağıdaki Dockerfile'a bazı açıklayıcı yorumlar ekledik.

```Dockerfile title="Dockerfile"
# Micromamba temel docker imajından başla
FROM mambaorg/micromamba:1.5.10-noble
# conda.yml dosyasını konteynere kopyala
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Nextflow'un kullanması için çeşitli yardımcı programları ve conda.yml dosyasındaki paketleri yükle
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
# Konteyneri root kullanıcısı olarak çalıştır
USER root
# PATH ortam değişkenini micromamba kurulum dizinini içerecek şekilde ayarla
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
```

Bakacağımız ikinci öğe, konteyner imajına yüklenmesi gereken paketlerin listesini içeren `conda.yml` dosyasıdır.

```conda.yml title="conda.yml"
channels:
- conda-forge
- bioconda
dependencies:
- pip
- pip:
  - quote==3.0.0 #
```

Bu dosyaların içeriğini `containers/build` dizininde bulunan taslakların içine kopyalayın, ardından konteyner imajını kendiniz oluşturmak için aşağıdaki komutu çalıştırın.

!!! note "Not"

    Konteyner imajını `quote` adı ve `latest` etiketi ile etiketlemek için `-t quote:latest` bayrağını kullanıyoruz.
    Bu sistemi çalıştırırken konteyner imajına atıfta bulunmak için bu etiketi kullanabileceğiz.

```bash
docker build -t quote:latest containers/build
```

Oluşturma işlemi tamamlandıktan sonra, az önce oluşturduğunuz konteyner imajını çalıştırabilirsiniz.

```bash
docker run --rm quote:latest quote "Margaret Oakley Dayhoff"
```

### Önemli Noktalar

Nextflow pipeline'larınızda kullanmak istediğiniz bir araç için konteyner imajı almanın iki farklı yolunu öğrendiniz: Seqera Containers'ı kullanmak ve konteyner imajını kendiniz oluşturmak.

### Sırada ne var?

Bu eğitim serisinin [sonraki bölümüne](./04_hello_genomics.md) devam etmek için ihtiyacınız olan her şeye sahipsiniz.
Ayrıca `quote` konteynerini kullanarak bilgisayar/biyoloji öncüleri hakkında alıntılar almak ve bunları `cowsay` konteynerini kullanarak görüntülemek için isteğe bağlı bir alıştırmayla devam edebilirsiniz.

---

## 2. İneği ünlü bilim insanlarından alıntı yaptırın

Bu bölüm, şimdiye kadar öğrendiklerinizi pratik yapmak için bazı ek alıştırmalar içerir.
Bu alıştırmaları yapmak, eğitimin sonraki bölümlerini anlamak için _gerekli değildir_, ancak ineği ünlü bilim insanlarından alıntı yaptırmayı nasıl yapacağınızı çözerek öğrendiklerinizi pekiştirmenin eğlenceli bir yolunu sağlar.

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

### 2.1. getQuote process kullanmak için `hello-containers.nf` betiğini değiştirin

`containers/data/pioneers.csv` dosyasında bilgisayar ve biyoloji öncülerinin bir listesi var.
Üst düzeyde, bu alıştırmayı tamamlamak için şunları yapmanız gerekecek:

- Varsayılan `params.input_file` parametresini `pioneers.csv` dosyasını işaret edecek şekilde değiştirin.
- Her girdi için bir alıntı almak üzere `quote` konteynerini kullanan bir `getQuote` process oluşturun.
- Alıntıyı görüntülemek için `getQuote` process'inin çıktısını `cowsay` process'ine bağlayın.

`quote` konteyner imajı için, önceki ek alıştırmada kendiniz oluşturduğunuz imajı veya Seqera Containers'dan aldığınız imajı kullanabilirsiniz.

!!! tip "İpucu"

    getQuote process'inizin `script` bloğu için iyi bir seçenek şu olabilir:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Bu alıştırmanın çözümünü `containers/solutions/hello-containers-4.1.nf` dosyasında bulabilirsiniz.

### 2.2. Nextflow pipeline'ınızı `quote` ve `sayHello` modlarında çalışabilecek şekilde değiştirin

Pipeline'ınıza hem `quote` hem de `sayHello` için tasarlanmış girdileri kabul etmesine izin verecek bazı dallanma mantığı ekleyin.
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

!!! tip "İpucu"

    Bir process'in çıktı kanalına bir ad atamak için `new_ch = processName.out` kullanabilirsiniz.

Bu alıştırmanın çözümünü `containers/solutions/hello-containers-4.2.nf` dosyasında bulabilirsiniz.

### Önemli Noktalar

Nextflow'da process'leri çalıştırmak için konteynerleri nasıl kullanacağınızı ve pipeline'larınıza bazı dallanma mantığını nasıl ekleyeceğinizi biliyorsunuz!

### Sırada ne var?

Kutlama yapın, bir germe molası verin ve biraz su için!

Hazır olduğunuzda, şimdiye kadar öğrendiklerinizi daha gerçekçi bir veri analizi kullanım durumuna nasıl uygulayacağınızı öğrenmek için bu eğitim serisinin Bölüm 3'üne geçin.
