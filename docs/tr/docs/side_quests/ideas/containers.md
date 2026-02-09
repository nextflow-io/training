# Bölüm 1: Daha Fazla Konteyner

[YAPILACAK]

---

## 1. Konteyner imajlarını nasıl bulabilir veya oluşturabilirsiniz

Bazı yazılım geliştiricileri, yazılımları için Docker Hub gibi konteyner kayıt defterlerinde mevcut olan konteyner imajları sağlar, ancak birçoğu sağlamaz.
Bu isteğe bağlı bölümde, Nextflow pipeline'larınızda kullanmak istediğiniz araçlar için konteyner imajı almanın iki yolunu göstereceğiz: Seqera Containers kullanarak ve konteyner imajını kendiniz oluşturarak.

Bu bölümün sonundaki alıştırmada kullanılacak olan `quote` pip paketi için bir konteyner imajı alacak/oluşturacaksınız.

### 1.1. Seqera Containers'dan bir konteyner imajı alın

Seqera Containers, pip ve conda (bioconda dahil) ile kurulabilen araçlar için konteyner imajları oluşturan ücretsiz bir hizmettir.
[Seqera Containers](https://www.seqera.io/containers/) adresine gidin ve `quote` pip paketini arayın.

![Seqera Containers](img/seqera-containers-1.png)

`quote` pip paketi için bir konteyner imajı talep etmek üzere "+Add" ve ardından "Get Container" seçeneklerine tıklayın.

![Seqera Containers](img/seqera-containers-2.png)

Paketin bu sürümü için ilk kez bir topluluk konteyneri oluşturuluyorsa, tamamlanması birkaç dakika sürebilir.
Sizin için oluşturulan konteyner imajının URI'sini (örn. `community.wave.seqera.io/library/pip_quote:ae07804021465ee9`) kopyalamak için tıklayın.

Artık konteyner imajını kullanarak `quote` komutunu çalıştırabilir ve Grace Hopper'dan rastgele bir söz alabilirsiniz.

```bash
docker run --rm community.wave.seqera.io/library/pip_quote:ae07804021465ee9 quote "Grace Hopper"
```

Çıktı:

```console title="Output"
Humans are allergic to change. They love to say, 'We've always done it
this way.' I try to fight that. That's why I have a clock on my wall
that runs counter-clockwise.
```

### 1.2. Konteyner imajını kendiniz oluşturun

`quote` pip paketi için konteyner imajını kendimiz oluşturmak üzere Seqera Containers web sitesinden bazı oluşturma ayrıntılarını kullanalım.
Seqera Containers web sitesine dönün ve "Build Details" düğmesine tıklayın.

Bakacağımız ilk öğe, konteyner imajını oluşturmak için gereken tüm komutları içeren bir tür betik dosyası olan `Dockerfile`'dır.
Her bir bölümün ne yaptığını anlamanıza yardımcı olmak için aşağıdaki Dockerfile'a bazı açıklayıcı yorumlar ekledik.

```Dockerfile title="Dockerfile"
# micromamba temel docker imajından başla
FROM mambaorg/micromamba:1.5.10-noble
# conda.yml dosyasını konteynere kopyala
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Nextflow'un kullanması için çeşitli yardımcı programları ve conda.yml dosyasındaki paketleri kur
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

Bakacağımız ikinci öğe, konteyner imajına kurulması gereken paketlerin listesini içeren `conda.yml` dosyasıdır.

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
    Bu sisteme çalıştırırken konteyner imajına başvurmak için bu etiketi kullanabileceğiz.

```bash
docker build -t quote:latest containers/build
```

Oluşturma işlemi tamamlandıktan sonra, az önce oluşturduğunuz konteyner imajını çalıştırabilirsiniz.

```bash
docker run --rm quote:latest quote "Margaret Oakley Dayhoff"
```

### Özet

Nextflow pipeline'larınızda kullanmak istediğiniz bir araç için konteyner imajı almanın iki farklı yolunu öğrendiniz: Seqera Containers kullanarak ve konteyner imajını kendiniz oluşturarak.

### Sırada ne var?

Bu eğitim serisinin [bir sonraki bölümüne](./04_hello_genomics.md) devam etmek için ihtiyacınız olan her şeye sahipsiniz.
Ayrıca `quote` konteynerini kullanarak bilgisayar/biyoloji öncüleri hakkında alıntılar almak ve bunları `cowsay` konteyneri kullanarak çıktı olarak vermek için isteğe bağlı bir alıştırma ile devam edebilirsiniz.

---

## 2. İneği ünlü bilim insanlarından alıntı yaptırın

Bu bölüm, şimdiye kadar öğrendiklerinizi pratik yapmak için bazı ek alıştırmalar içerir.
Bu alıştırmaları yapmak, eğitimin sonraki bölümlerini anlamak için _gerekli değildir_, ancak ineği ünlü bilim insanlarından alıntı yaptırmayı nasıl yapacağınızı çözerek öğrendiklerinizi pekiştirmenin eğlenceli bir yolunu sunar.

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

### 2.1. `hello-containers.nf` betiğini bir getQuote süreci kullanacak şekilde değiştirin

`containers/data/pioneers.csv` dosyasında bilgisayar ve biyoloji öncülerinin bir listesi var.
Üst düzeyde, bu alıştırmayı tamamlamak için şunları yapmanız gerekecek:

- Varsayılan `params.input_file` değerini `pioneers.csv` dosyasına işaret edecek şekilde değiştirin.
- Her girdi için bir alıntı almak üzere `quote` konteynerini kullanan bir `getQuote` süreci oluşturun.
- Alıntıyı görüntülemek için `getQuote` sürecinin çıktısını `cowsay` sürecine bağlayın.

`quote` konteyner imajı için, önceki ek alıştırmada kendiniz oluşturduğunuz imajı veya Seqera Containers'dan aldığınız imajı kullanabilirsiniz.

!!! tip "İpucu"

    getQuote sürecinizin `script` bloğu için iyi bir seçim şu olabilir:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Bu alıştırmanın bir çözümünü `containers/solutions/hello-containers-4.1.nf` dosyasında bulabilirsiniz.

### 2.2. Nextflow pipeline'ınızı `quote` ve `sayHello` modlarında çalıştırılabilmesine izin verecek şekilde değiştirin.

Pipeline'ınıza hem `quote` hem de `sayHello` için tasarlanmış girdileri kabul etmesine izin vermek üzere bazı dallanma mantığı ekleyin.
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

!!! tip "İpucu"

    Bir sürecin çıktı kanalına bir ad atamak için `new_ch = processName.out` kullanabilirsiniz.

Bu alıştırmanın bir çözümünü `containers/solutions/hello-containers-4.2.nf` dosyasında bulabilirsiniz.

### Özet

Nextflow'da süreçleri çalıştırmak için konteynerleri nasıl kullanacağınızı ve pipeline'larınıza bazı dallanma mantığını nasıl ekleyeceğinizi biliyorsunuz!

### Sırada ne var?

Kutlama yapın, bir esneme molası verin ve biraz su için!

Hazır olduğunuzda, şimdiye kadar öğrendiklerinizi daha gerçekçi bir veri analizi kullanım durumuna nasıl uygulayacağınızı öğrenmek için bu eğitim serisinin Bölüm 3'üne geçin.
