# Bölüm 2: nf-core/molkart'ı Çalıştırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bölüm 1'de, Nextflow çalıştırmasının temellerini anlamak için basit bir Hello World iş akışı çalıştırdık.
Şimdi gerçek dünyadan bir biyogörüntüleme pipeline'ı çalıştıracağız: **nf-core/molkart**.

Bu pipeline, Resolve Bioscience'dan gelen Molecular Cartography uzamsal transkriptomik verilerini işler.
Ancak burada öğreneceğiniz Nextflow kalıpları, herhangi bir nf-core pipeline'ı veya üretim iş akışı için geçerlidir.

## 1. nf-core pipeline'larını anlama

Pipeline'ı çalıştırmadan önce, nf-core'un ne olduğunu ve iş akışlarını çalıştırmak için neden önemli olduğunu anlayalım.

### 1.1. nf-core nedir?

[nf-core](https://nf-co.re/), topluluk tarafından yönlendirilen yüksek kaliteli Nextflow pipeline'ları koleksiyonudur.
Tüm nf-core pipeline'ları aynı yapı ve kurallara uyar; bu da birini çalıştırmayı öğrendiğinizde herhangi birini çalıştırabileceğiniz anlamına gelir.

nf-core pipeline'larının temel özellikleri:

- **Standartlaştırılmış yapı**: Tüm pipeline'lar tutarlı parametre adlarına ve kullanım kalıplarına sahiptir
- **Yerleşik test verileri**: Her pipeline, hızlı doğrulama için test profilleri içerir
- **Kapsamlı dokümantasyon**: Detaylı kullanım talimatları ve parametre açıklamaları
- **Kalite kontrol**: MultiQC kullanarak otomatik kalite kontrol raporları
- **Konteyner desteği**: Tekrarlanabilirlik için önceden hazırlanmış konteynerler

!!! tip "nf-core hakkında daha fazla bilgi edinmek ister misiniz?"

    nf-core pipeline geliştirmeye derinlemesine bir giriş için [Build with nf-core](../../hello_nf-core/index.md) eğitim kursuna göz atın.
    Sıfırdan nf-core pipeline'ları oluşturmayı ve özelleştirmeyi kapsar.

### 1.2. molkart pipeline'ı

![nf-core/molkart pipeline'ı](img/molkart.png)

[nf-core/molkart](https://nf-co.re/molkart) pipeline'ı, uzamsal transkriptomik görüntüleme verilerini birkaç aşamadan geçirerek işler:

1. **Görüntü ön işleme**: Izgara deseni doldurma ve isteğe bağlı kontrast iyileştirme
2. **Hücre segmentasyonu**: Çoklu algoritma seçenekleri (Cellpose, Mesmer, ilastik, Stardist)
3. **Nokta ataması**: Transkript noktalarını segmente edilmiş hücrelere atama
4. **Kalite kontrol**: Kapsamlı kalite kontrol raporları oluşturma

Temel çıktılar şunlardır:

- Hücre-transkript sayım tabloları
- Segmentasyon maskeleri
- MultiQC kalite kontrol raporu

---

## 2. Test verileri ile molkart çalıştırma

Başlamadan önce, kodunu inceleyebilmemiz için molkart deposunu yerel olarak klonlayalım:

```bash
cd /workspaces/training/nf4-science/imaging
git clone --branch 1.2.0 --depth 1 https://github.com/nf-core/molkart
```

Bu, tam pipeline kaynak kodunu içeren bir `molkart/` dizini oluşturur.

!!! note "Neden yerel olarak klonluyoruz?"

    Genellikle, nf-core pipeline'larını doğrudan GitHub'dan `nextflow run nf-core/molkart -r 1.2.0` kullanarak çalıştırırsınız.
    Nextflow, istenen pipeline sürümünü sizin için otomatik olarak `$HOME/.nextflow/assets/nf-core/molkart` dizinine indirir ve oradan çalıştırır.
    Ancak bu eğitim için, kodu daha kolay inceleyebilmemiz amacıyla pipeline'ı farklı bir yerel dizine klonluyoruz.

### 2.1. Konteyner gereksinimlerini anlama

Tam pipeline'ı çalıştırmadan önce, konteynerların neden nf-core pipeline'ları için gerekli olduğunu öğrenelim.

Pipeline'ın parametrelerini bir parametre dosyası kullanarak sağlayacağız.
Parametre dosyası, her parametreyi ve değerini listeleyen bir YAML dosyasıdır; bu sayede yazılan değerler (tam sayılar gibi) bozulmadan korunur ve komut satırı kısa tutulur.

Çalışma dizininde zaten bir `params.yaml` dosyası bulunmaktadır:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Bu parametreleri açıklayalım:

- `input`: Örnek meta verilerini içeren samplesheet'in yolu
- `mindagap_tilesize`, `mindagap_boxsize`, `mindagap_loopnum`: Izgara deseni doldurma için parametreler
- `clahe_pyramid_tile`: Kontrast iyileştirme için çekirdek boyutu
- `segmentation_method`: Hücre segmentasyonu için hangi algoritma(lar)ın kullanılacağı
- `outdir`: Sonuçların nereye kaydedileceği

Bu parametreleri kullanarak pipeline'ı çalıştırmayı deneyelim:

```bash
nextflow run ./molkart -params-file params.yaml
```

!!! Warning "Bu komut başarısız olacak - bu kasıtlı!"

    Neden gerekli olduklarını göstermek için kasıtlı olarak bunu konteynerler olmadan çalıştırıyoruz.

Birkaç dakika sonra, şuna benzer bir hata göreceksiniz:

??? failure "Komut çıktısı"

    ```console
    ERROR ~ Error executing process > 'NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)'

    Caused by:
      Process `NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)` terminated with an error exit status (127)

    Command executed:

      duplicate_finder.py \
          spots.txt \
          90

    Command exit status:
      127

    Command error:
      .command.sh: line 3: duplicate_finder.py: command not found
    ```

**Ne oluyor?**

`command not found` hatası (çıkış durumu 127), Nextflow'un `duplicate_finder.py` çalıştırmayı denediği ancak sisteminizde bulamadığı anlamına gelir.
Bunun nedeni:

1. Pipeline'ın özel biyoinformatik yazılımlarının kurulu olmasını beklemesi
2. Bu araçların (örneğin `duplicate_finder.py`, `apply_clahe.dask.py`, vb.) standart Linux dağıtımlarının bir parçası olmaması
3. Konteynerler olmadan, Nextflow'un komutları doğrudan yerel makinenizde çalıştırmaya çalışması

**Bu araçların nereden gelmesi gerekiyor?**

Yazılım gereksinimlerini nasıl bildirdiğini görmek için işlem modüllerinden birini inceleyelim.

CLAHE ön işleme modülünü açın:

```bash
code molkart/modules/local/clahe/main.nf
```

5. satıra bakın - şunu göreceksiniz:

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

Bu satır Nextflow'a şunu söyler: "Bu işlemi çalıştırmak için, gerekli tüm yazılımları içeren `ghcr.io/schapirolabor/molkart-local:v0.0.4` Docker imajını kullan."

Her işlem, hangi konteyner imajının gerekli araçları sağladığını bildirir.
Ancak Nextflow bu konteynerları yalnızca ona söylerseniz kullanır!

**Çözüm: Yapılandırmada Docker'ı etkinleştirin**

### 2.2. Docker'ı yapılandırma ve pipeline'ı başlatma

Docker'ı etkinleştirmek için `nextflow.config` dosyasında `docker.enabled` değerini `false`'dan `true`'ya değiştirmemiz gerekir.

Yapılandırma dosyasını açın:

```bash
code nextflow.config
```

`docker.enabled = false` ifadesini `docker.enabled = true` olarak değiştirin:

```groovy
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}
```

Şimdi pipeline'ı tekrar çalıştırın; bu sefer daha sonra karşılaştırabilmek için üç segmentasyon yönteminin tamamını çalıştırın:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "mesmer,cellpose,stardist"
```

Bu sefer Nextflow:

1. Yapılandırmadan `docker.enabled = true` ayarını okuyacak
2. Gerekli Docker imajlarını çekecek (yalnızca ilk seferde)
3. Her işlemi belirtilen konteyner içinde çalıştıracak
4. Tüm araçlar konteynerler içinde mevcut olduğu için başarıyla çalışacak

!!! Tip "Konteynerlar neden önemlidir"

    Çoğu nf-core pipeline'ı konteynerizasyonu (Docker, Singularity, Podman, vb.) **gerektirir** çünkü:

    - Standart ortamlarda bulunmayan özel biyoinformatik yazılımları kullanırlar
    - Konteynerler tekrarlanabilirliği sağlar - aynı yazılım sürümleri her yerde çalışır
    - Düzinelerce aracı ve bağımlılıklarını manuel olarak kurmanız gerekmez

    Nextflow'da konteynerler hakkında daha fazla bilgi için, Hello Nextflow eğitiminden [Hello Containers](../../hello_nextflow/05_hello_containers.md) bölümüne bakın.

### 2.3. Çalıştırmayı izleme

Pipeline çalışırken, buna benzer bir çıktı göreceksiniz:

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `./molkart/main.nf` [exotic_banach] revision: e4152308ec

    WARN: Unrecognized config option 'validation.help.enabled'
    WARN: Unrecognized config option 'validation.defaultIgnoreParams'
    WARN: Unrecognized config option 'validation.monochromeLogs'

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/molkart 1.2.0
    ------------------------------------------------------
    Segmentation methods and options
      segmentation_method: mesmer,cellpose,stardist

    Image preprocessing
      mindagap_boxsize   : 7
      mindagap_loopnum   : 100
      clahe_kernel       : 25
      mindagap_tilesize  : 90
      clahe_pyramid_tile : 368

    Input/output options
      input              : data/samplesheet.csv
      outdir             : results

    Generic options
      trace_report_suffix: 2026-06-23_16-25-30

    Core Nextflow options
      runName            : exotic_banach
      containerEngine    : docker
      launchDir          : /workspaces/training/nf4-science/imaging
      workDir            : /workspaces/training/nf4-science/imaging/work
      projectDir         : /workspaces/training/nf4-science/imaging/molkart
      userName           : root
      profile            : standard
      configFiles        : /workspaces/training/nf4-science/imaging/molkart/nextflow.config, /workspaces/training/nf4-science/imaging/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.10650748

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/molkart/blob/master/CITATIONS.md

    executor >  local (22)
    [b4/e57ff1] NFC…NDAGAP_MINDAGAP (mem_only) | 2 of 2 ✔
    [2e/cf8910] NFC…T:MOLKART:CLAHE (mem_only) | 2 of 2 ✔
    [94/b1ef70] NFC…RT:CREATE_STACK (mem_only) | 1 of 1 ✔
    [58/4b6426] NFC…DUPLICATEFINDER (mem_only) | 1 of 1 ✔
    [c2/ef7002] NFC…DEEPCELL_MESMER (mem_only) | 1 of 1 ✔
    [4f/ddd328] NFC…OLKART:STARDIST (mem_only) | 1 of 1 ✔
    [8c/dd0362] NFC…OLKART:CELLPOSE (mem_only) | 1 of 1 ✔
    [08/2ddcb5] NFC…KART:MASKFILTER (mem_only) | 3 of 3 ✔
    [56/e30de0] NFC…LKART:SPOT2CELL (mem_only) | 3 of 3 ✔
    [b0/5aa635] NFC…:CREATE_ANNDATA (mem_only) | 3 of 3 ✔
    [5e/39d5d0] NFC…LKART:MOLKARTQC (mem_only) | 3 of 3 ✔
    [8e/8bd365] NFCORE_MOLKART:MOLKART:MULTIQC | 1 of 1 ✔
    -[nf-core/molkart] Pipeline completed successfully-
    Completed at: 23-Jun-2026 16:31:40
    Duration    : 6m 9s
    CPU hours   : (a few seconds)
    Succeeded   : 22
    ```

Bu çıktının, pipeline'ın takip ettiği nf-core kuralları nedeniyle Hello World örneğimizden daha detaylı olduğuna dikkat edin:

- Pipeline, sürümünü ve logosunu gösterir
- Yapılandırma parametreleri görüntülenir
- Birden fazla işlem paralel olarak çalışır (birden fazla işlem satırıyla gösterilir)
- İşlem adları tam modül yolunu içerir (örneğin, `NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP`)

### 2.4. İşlem çalıştırmasını anlama

`executor > local (22)` satırı size şunları söyler:

- **executor**: Hangi hesaplama ortamının kullanıldığı (`local` = makineniz)
- **(22)**: Başlatılan toplam görev sayısı

Her işlem satırı şunları gösterir:

- **Hash** (`[b4/e57ff1]`): Çalışma dizini tanımlayıcısı (önceki gibi)
- **İşlem adı**: Tam modül yolu ve işlem adı
- **Girdi tanımlayıcısı**: Parantez içinde örnek adı
- **İlerleme**: Görev sayısı ve tamamlanma durumu (örneğin, `1 of 1 ✔`)

### Özetle

Bir nf-core pipeline'ını test verileri ile nasıl başlatacağınızı ve çalıştırma çıktısını nasıl yorumlayacağınızı biliyorsunuz.

### Sırada ne var?

Sonuçları nerede bulacağınızı ve nasıl yorumlayacağınızı öğrenin.

---

## 3. Çıktıları bulma ve inceleme

Pipeline başarıyla tamamlandığında, bir tamamlanma mesajı ve çalıştırma özeti göreceksiniz.

### 3.1. Sonuçlar dizinini bulma

Varsayılan olarak, nf-core pipeline'ları çıktıları `outdir` parametresiyle belirtilen bir dizine yazar; bunu `results/` olarak ayarladık.

İçeriği listeleyin:

```bash
tree results/
```

Birkaç alt dizin görmelisiniz:

```console title="results/"
results/
├── anndata/
├── clahe/
├── mindagap/
├── molkartqc/
├── multiqc/
├── pipeline_info/
├── segmentation/
├── spot2cell/
└── stack/
```

Her alt dizin, pipeline'ın belirli bir aşamasından çıktılar içerir:

- **mindagap/**: MindaGap ön işleme adımından ızgara doldurulmuş görüntüler
- **clahe/**: CLAHE ön işlemeden kontrast artırılmış görüntüler
- **stack/**: Segmentasyon için oluşturulan çok kanallı görüntü yığınları
- **segmentation/**: Farklı algoritmalardan segmentasyon sonuçları (cellpose/, mesmer/, stardist/, filtered_masks/)
- **spot2cell/**: Hücre-transkript sayım tabloları
- **anndata/**: Hücre-transkript matrisleri ve uzamsal koordinatları içeren AnnData nesneleri
- **molkartqc/**: Nokta ataması için kalite kontrol metrikleri
- **multiqc/**: Kapsamlı kalite kontrol raporu
- **pipeline_info/**: Çalıştırma raporları ve loglar

### 3.2. MultiQC raporunu inceleme

MultiQC raporu, tüm pipeline adımlarından kalite metriklerini toplayan kapsamlı bir HTML dosyasıdır.

Raporu dosya tarayıcısında açın ve ardından doğrudan VS Code'da işlenmiş olarak görmek için "Show Preview" düğmesine tıklayın.

Rapor şunları içerir:

- Tüm örnekler için genel istatistikler
- Ön işleme metrikleri
- Segmentasyon kalite metrikleri
- Tespit edilen hücre ve nokta sayısı

!!! Tip "İpucu"

    MultiQC raporları genellikle tüm nf-core pipeline'larına dahil edilir.
    Her zaman pipeline çalıştırması ve veri kalitesi hakkında üst düzey bir genel bakış sağlarlar.

### 3.3. Hücre-transkript tablolarını inceleme

En önemli bilimsel çıktı, hücre-transkript sayım tablosudur.
Bu, her hücrede her transkriptten kaç tanesinin tespit edildiğini söyler.

spot2cell dizinine gidin:

```bash
ls results/spot2cell/
```

Şunlar gibi dosyalar bulacaksınız:

- `cellxgene_mem_only_cellpose.csv`: Cellpose segmentasyonu kullanarak hücre-transkript tablosu
- `cellxgene_mem_only_mesmer.csv`: Mesmer segmentasyonu kullanarak hücre-transkript tablosu
- `cellxgene_mem_only_stardist.csv`: Stardist segmentasyonu kullanarak hücre-transkript tablosu

Bu test veri setinde yalnızca 1 örnek çalıştırdık, ancak gerçek bir deneyde bu tabloları her örnek için elde ederdik.
Nextflow'un birden fazla segmentasyon yöntemini paralel olarak işleyerek sonuçları karşılaştırmayı nasıl kolaylaştırdığına dikkat edin.

### 3.4. Çalıştırma raporlarını görüntüleme

Nextflow otomatik olarak birkaç çalıştırma raporu oluşturur.

pipeline_info dizinini kontrol edin:

```bash
ls results/pipeline_info/
```

Anahtar dosyalar:

- **execution_report.html**: Zaman çizelgesi ve kaynak kullanım görselleştirmesi
- **execution_timeline.html**: İşlem çalıştırmasının Gantt grafiği
- **execution_trace.txt**: Detaylı görev çalıştırma metrikleri
- **pipeline_dag.html**: İş akışı yapısını gösteren yönlendirilmiş döngüsüz grafik

Kaynak kullanımını görmek için çalıştırma raporunu açın:

```bash
code results/pipeline_info/execution_report.html
```

Bu şunları gösterir:

- Her işlemin ne kadar sürdüğü
- CPU ve bellek kullanımı
- Hangi görevlerin önbellekte olduğu ve hangilerinin çalıştırıldığı

!!! Tip "İpucu"

    Bu raporlar, kaynak tahsisini optimize etmek ve performans sorunlarını gidermek için son derece kullanışlıdır.

### Özetle

Pipeline çıktılarını nasıl bulacağınızı, kalite kontrol raporlarını nasıl inceleyeceğinizi ve çalıştırma metriklerine nasıl erişeceğinizi biliyorsunuz.

### Sırada ne var?

Çalışma dizini ve Nextflow'un ara dosyaları nasıl yönettiği hakkında bilgi edinin.

---

## 4. Çalışma dizinini keşfetme

Hello World örneğimizde olduğu gibi, tüm gerçek çalışma `work/` dizininde gerçekleşir.

### 4.1. Çalışma dizini yapısını anlama

Çalışma dizini, çalıştırılan her görev için bir alt dizin içerir.
22 göreve sahip bu pipeline çalıştırması için 22 çalışma alt dizini olacaktır.

Çalışma dizinini listeleyin:

```bash
ls -d work/*/*/ | head -5
```

Bu, ilk 5 görev dizinini gösterir.

### 4.2. Bir görev dizinini inceleme

Konsol çıktısından segmentasyon işlem hash'lerinden birini seçin (örneğin, `[3m/4n5o6p]`) ve içine bakın:

```bash
ls -la work/3m/4n5o6p*/
```

Şunları göreceksiniz:

- **.command.\*** dosyaları: Nextflow çalıştırma betikleri ve loglar (önceki gibi)
- **Hazırlanmış girdi dosyaları**: Gerçek girdi dosyalarına sembolik bağlantılar
- **Çıktı dosyaları**: Segmentasyon maskeleri, ara sonuçlar, vb.

Hello World'den temel fark:

- Gerçek pipeline'lar büyük girdi dosyalarını hazırlar (görüntüler, referans verileri)
- Çıktı dosyaları oldukça büyük olabilir (segmentasyon maskeleri, işlenmiş görüntüler)
- Görev başına birden fazla girdi ve çıktı dosyası

!!! Tip "İpucu"

    Bir süreç başarısız olursa, çalışma dizinine gidebilir, hata mesajları için `.command.err` dosyasını inceleyebilir ve hatta sorunu ayıklamak için `.command.sh` dosyasını manuel olarak yeniden çalıştırabilirsiniz.

### 4.3. Çalışma dizini temizliği

Çalışma dizini, birden fazla pipeline çalıştırması üzerinden oldukça büyüyebilir.
Bölüm 1'de öğrendiğimiz gibi, eski çalıştırmalardan çalışma dizinlerini kaldırmak için `nextflow clean` kullanabilirsiniz.

Ancak büyük ara dosyalara sahip nf-core pipeline'ları için düzenli olarak temizlik yapmak özellikle önemlidir.

### Özetle

nf-core pipeline'larının çalışma dizinlerini nasıl düzenlediğini ve hata ayıklama için bireysel görevleri nasıl inceleyeceğinizi anlıyorsunuz.

### Sırada ne var?

Nextflow önbelleği ve başarısız pipeline çalıştırmalarını nasıl devam ettireceğinizi öğrenin.

---

## 5. Bir pipeline çalıştırmasını devam ettirme

Nextflow'un en güçlü özelliklerinden biri, bir pipeline'ı hata noktasından devam ettirme yeteneğidir.

### 5.1. Önbellek mekanizması

Bir pipeline'ı `-resume` ile çalıştırdığınızda, Nextflow:

1. Her görev için önbelleği kontrol eder
2. Girdiler, kod ve parametreler özdeşse, önbellekteki sonucu yeniden kullanır
3. Yalnızca değişen veya başarısız olan görevleri yeniden çalıştırır

Bu, çalıştırmanın geç aşamalarında hataların oluşabileceği uzun süreli pipeline'lar için gereklidir.

### 5.2. molkart ile resume deneme

Aynı komutu tekrar çalıştırın, ancak `-resume` ekleyin:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "mesmer,cellpose,stardist" -resume
```

Şuna benzer bir çıktı görmelisiniz:

```console
executor >  local (1)
[43/e03702] NFC…NDAGAP_MINDAGAP (mem_only) | 2 of 2, cached: 2 ✔
[2e/cf8910] NFC…T:MOLKART:CLAHE (mem_only) | 2 of 2, cached: 2 ✔
[94/b1ef70] NFC…RT:CREATE_STACK (mem_only) | 1 of 1, cached: 1 ✔
[58/4b6426] NFC…DUPLICATEFINDER (mem_only) | 1 of 1, cached: 1 ✔
[c2/ef7002] NFC…DEEPCELL_MESMER (mem_only) | 1 of 1, cached: 1 ✔
[4f/ddd328] NFC…OLKART:STARDIST (mem_only) | 1 of 1, cached: 1 ✔
[8c/dd0362] NFC…OLKART:CELLPOSE (mem_only) | 1 of 1, cached: 1 ✔
[08/2ddcb5] NFC…KART:MASKFILTER (mem_only) | 3 of 3, cached: 3 ✔
[56/e30de0] NFC…LKART:SPOT2CELL (mem_only) | 3 of 3, cached: 3 ✔
[b0/5aa635] NFC…:CREATE_ANNDATA (mem_only) | 3 of 3, cached: 3 ✔
[5e/39d5d0] NFC…LKART:MOLKARTQC (mem_only) | 3 of 3, cached: 3 ✔
[73/239f45] NFCORE_MOLKART:MOLKART:MULTIQC | 1 of 1 ✔
-[nf-core/molkart] Pipeline completed successfully-
```

Her ön işleme ve segmentasyon sürecindeki `cached: N` ek açıklamasına dikkat edin; bu görevler yeniden çalıştırılmak yerine önbellekten kullanıldı.

### 5.3. Resume ne zaman yararlıdır

Resume özellikle şu durumlarda değerlidir:

- Bir pipeline kaynak sınırları nedeniyle başarısız olduğunda (yetersiz bellek, zaman sınırı aşımı)
- Üst akış adımlarını yeniden çalıştırmadan alt akış işlemlerini değiştirmeniz gerektiğinde
- Veri indirme sırasında ağ bağlantınız düştüğünde
- Hesaplamayı yeniden yapmadan ek çıktılar eklemek istediğinizde

!!! Warning "Uyarı"

    Resume yalnızca girdi verilerini, pipeline kodunu veya parametreleri değiştirmemişseniz çalışır.
    Bunlardan herhangi birini değiştirirseniz, Nextflow etkilenen görevleri doğru bir şekilde yeniden çalıştıracaktır.

### Özetle

Başarılı görevleri tekrarlamadan pipeline'ları verimli bir şekilde yeniden çalıştırmak için `-resume` kullanmayı biliyorsunuz.

### Sırada ne var?

Artık nf-core/molkart'ı test verileri ile çalıştırabildiğinize göre, kendi veri setleriniz için nasıl yapılandıracağınızı öğrenmeye hazırsınız.
