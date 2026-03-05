# Bölüm 3: Girdileri düzenleme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bölüm 2'de molkart'ı komut satırında birden fazla parametre ile çalıştırdık.
Şimdi girdileri yönetmek için iki daha iyi yaklaşım öğreneceğiz: **parametre dosyaları** ve **örnek tabloları**.

## 1. Parametre dosyalarını kullanma

### 1.1. Uzun komut satırlarıyla ilgili sorun

Bölüm 2'deki komutumuzu hatırlayın:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results
```

Bu işe yarar, ancak tekrar üretmek, paylaşmak veya değiştirmek zordur.
Aynı analizi önümüzdeki ay tekrar çalıştırmanız gerekirse ne olur?
Bir iş arkadaşınız tam ayarlarınızı kullanmak isterse ne olur?

### 1.2. Çözüm: Parametre dosyası kullanın

`params.yaml` adında bir dosya oluşturun:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Artık komutunuz şu hale gelir:

```bash
nextflow run ./molkart -params-file params.yaml -resume
```

Bu kadar! Parametre dosyası tam yapılandırmanızı belgelendirir ve tekrar çalıştırmayı veya paylaşmayı kolaylaştırır.

### 1.3. Parametreleri geçersiz kılma

Komut satırından belirli parametreleri geçersiz kılmaya devam edebilirsiniz:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

Yukarıdaki satır, `segmentation_method`'u `stardist` olarak değiştirir ve `--outdir` adını `params.yaml` dosyasındaki parametreler yerine `stardist_results` olarak ayarlar.
Ayrıca, `-resume` bayrağının önceki çalıştırmadaki ön işleme sonuçlarını yeniden kullanmamıza izin vererek zaman tasarrufu sağladığını görebilirsiniz.
Bu kalıbı, ardışık düzenin farklı varyasyonlarını hızlıca test etmek için kullanabilirsiniz.

### Çıkarım

Parametre dosyaları analizlerinizi tekrar üretilebilir ve paylaşımı kolay hale getirir.
Gerçek analiz çalışmaları için bunları kullanın.

### Sırada ne var?

Örnek tablolarının birden fazla örnek hakkındaki bilgileri nasıl düzenlediğini öğrenin.

---

## 2. Örnek tablosu kalıbı

### 2.1. Örnek tablosu nedir?

Örnek tablosu, girdi örneklerinizi açıklayan bir CSV dosyasıdır.
Her satır bir örnektir ve sütunlar o örnek için dosyaları ve meta verileri belirtir.

Kullandığımız örnek tablosuna bakalım:

```bash
cat data/samplesheet.csv
```

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
```

Sütunlar şunlardır:

- `sample`: Benzersiz örnek tanımlayıcısı
- `nuclear_image`: Çekirdek boyama görüntüsü (TIFF)
- `spot_table`: Transkript noktaları (TXT)
- `membrane_image`: Membran boyama görüntüsü (TIFF, isteğe bağlı)

### 2.2. Dosya yolları

Örnek tabloları birden fazla yol türünü kabul eder:

- **URL'ler**: Nextflow otomatik olarak indirir (yukarıda gösterildiği gibi)
- **Yerel yollar**: `data/nuclear.tiff` veya `/absolute/path/to/nuclear.tiff`
- **Bulut depolama**: `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

Aynı örnek tablosunda yol türlerini karıştırabilirsiniz.

### 2.3. Kendi örnek tablonuzu oluşturma

Önce test veri dosyalarını yerel olarak indirelim:

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

Şimdi örnek tablosunu bu yerel dosyalara referans verecek şekilde değiştirelim:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! warning "Uyarı"

    Örnek tablosundaki yolların, örnek tablosunun bulunduğu yere değil, Nextflow'u **çalıştırdığınız** yere göre göreceli olduğuna dikkat edin.

Son olarak, yerel dosya yollarına sahip örnek tablosuyla nf-core/molkart'ı bir kez daha çalıştıralım:

`nextflow run ./molkart -params-file params.yaml -resume`

Gördüğünüz gibi, Nextflow bu çalıştırmayı dosyalar Github'dan indirildiğinde olduğu gibi benzer şekilde yürütür. Bu, Nextflow'un harika özelliklerinden biridir, verinin nerede bulunduğuna bakılmaksızın verileri sizin için uygun şekilde hazırlar.

### Çıkarım

Örnek tabloları, dosya yollarıyla birlikte meta verilerinizi açıkça tanımlamanıza olanak tanıyan bir şekilde çok örnekli veri setlerini düzenler.
Çoğu nf-core ardışık düzeni bu kalıbı kullanır.

### Sırada ne var?

Artık girdileri ele aldığımıza göre, Nextflow ardışık düzenlerini farklı hesaplama ortamları için nasıl yapılandıracağımızı keşfedelim.
