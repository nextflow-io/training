# Bölüm 4: Yapılandırma

Bölüm 1-3'te Nextflow'u nasıl çalıştıracağımızı, bir nf-core pipeline'ını nasıl çalıştıracağımızı ve parametre dosyaları ile samplesheet'ler kullanarak girdileri nasıl yöneteceğimizi öğrendik.
Şimdi **yapılandırma dosyaları** ve **profiller** kullanarak pipeline'ları farklı hesaplama ortamları için nasıl yapılandıracağımızı keşfedeceğiz.

## Öğrenme hedefleri

Bu bölümün sonunda şunları yapabileceksiniz:

- Nextflow'un yapılandırmayı birden fazla kaynaktan nasıl çözdüğünü anlama
- Konteynerler ve test için nf-core yerleşik profillerini kullanma
- Farklı hesaplama ortamları için özel profiller oluşturma
- Süreç etiketlerini kullanarak kaynak taleplerini özelleştirme
- Kısıtlı ortamlarda kaynak limitlerini yönetme
- `nextflow config` ile çözümlenmiş yapılandırmayı inceleme

---

## 1. Nextflow yapılandırmasını anlama

### 1.1. Yapılandırma dosyası nedir?

Nextflow, **iş akışı mantığını** (ne yapılacağı) **yürütme ayarlarından** (nasıl ve nerede yapılacağı) ayırmak için yapılandırma dosyalarını kullanır.

Yapılandırma dosyaları şunları kontrol eder:

- Konteyner motorları (Docker, Singularity, Conda)
- Hesaplama kaynakları (CPU'lar, bellek, zaman)
- Yürütme platformları (yerel, HPC, bulut)
- Pipeline parametreleri

### 1.2. Yapılandırma önceliği

Nextflow yapılandırmayı birden fazla kaynaktan yükler ve sonraki kaynaklar önceki kaynakları geçersiz kılar:

1. **Pipeline yapılandırması**: Pipeline deposundaki `nextflow.config`
2. **Dizin yapılandırması**: Mevcut çalışma dizininizdeki `nextflow.config`
3. **Kullanıcı yapılandırması**: `~/.nextflow/config`
4. **Komut satırı**: Doğrudan iletilen parametreler ve seçenekler

Bu katmanlı yaklaşım, varsayılanları pipeline'da tutmanıza, kullanıcıya özel ayarlarla geçersiz kılmanıza ve komut satırında hızlı ayarlamalar yapmanıza olanak tanır.

### 1.3. Mevcut yapılandırmamız

Kullandığımız yapılandırmaya bakalım:

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

Bölüm 2'deki `docker.enabled = true` satırını yorum satırı haline getirelim veya geri değiştirelim ve molkart'ta bunun yerine bir profil kullanarak aynı sonucu nasıl elde edebileceğimizi bulalım.

---

## 2. Profilleri kullanma

### 2.1. Profiller nedir?

Profiller, `nextflow run` komutu aracılığıyla `-profile` bayrağı ile etkinleştirilebilen adlandırılmış yapılandırma kümeleridir.
Yapılandırma dosyalarını düzenlemeden farklı hesaplama senaryoları arasında geçiş yapmayı kolaylaştırırlar.

Tüm nf-core pipeline'ları kullanabileceğimiz bir dizi varsayılan profille birlikte gelir.

### 2.2. Yerleşik profilleri inceleme

Pipeline kod tabanıyla ilişkili `molkart/nextflow.config` dosyasında bunları inceleyelim:

```bash
code molkart/nextflow.config
```

`profiles` bloğunu arayın:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
    docker {
        docker.enabled          = true
        singularity.enabled     = false
        conda.enabled           = false
    }
    singularity {
        singularity.enabled     = true
        singularity.autoMounts  = true
        docker.enabled          = false
    }
    conda {
        conda.enabled           = true
        docker.enabled          = false
        conda.channels          = ['conda-forge', 'bioconda']
    }
}
```

Yaygın konteyner profilleri:

- `docker`: Docker konteynerlerini kullan (yerel geliştirme için en yaygın)
- `singularity`: Singularity/Apptainer kullan (HPC'de yaygın)
- `conda`: Conda ortamlarını kullan
- `apptainer`: Apptainer konteynerlerini kullan

### 2.3. nextflow.config yerine profillerle yeniden çalıştırma

Artık yerel `nextflow.config` dosyamızdaki docker yapılandırmasını devre dışı bıraktığımıza ve profilleri anladığımıza göre, `-profile` bayrağını kullanarak pipeline'ı yeniden çalıştıralım.

Daha önce Bölüm 3'te özel parametrelerimizle bir `params.yaml` dosyası oluşturduk.
Şimdi bunu yerleşik Docker profiliyle birleştirebiliriz:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

Her bayrağın ne yaptığını inceleyelim:

- `-profile docker`: molkart'ın `nextflow.config` dosyasından Docker profilini etkinleştirir, bu da `docker.enabled = true` ayarını yapar
- `-params-file params.yaml`: Tüm pipeline parametrelerini YAML dosyamızdan yükler
- `-resume`: Önceki çalıştırmalardan önbelleğe alınmış sonuçları yeniden kullanır

`-resume` kullandığımız için Nextflow, son çalıştırmadan bu yana bir şeyin değişip değişmediğini kontrol edecektir.
Parametreler, girdiler ve kod aynıysa, tüm görevler önbellekten alınacak ve pipeline neredeyse anında tamamlanacaktır.

```console title="Output (excerpt)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

Tüm süreçlerin `cached: 2` veya `cached: 1` gösterdiğine dikkat edin - hiçbir şey yeniden yürütülmedi!

### 2.4. Test profilleri

Test profilleri, pipeline'ın çalıştığını doğrulamanıza olanak tanıyan varsayılan girdi parametrelerini ve veri dosyalarını belirtmenin hızlı yollarını sağlar.
nf-core pipeline'ları her zaman en az iki test profili içerecektir:

- `test`: Hızlı test için küçük veri seti ve hızlı parametreler
- `test_full`: Daha büyük verilerle daha kapsamlı test

`includeConfig` yönergesi kullanılarak dahil edilen molkart'taki `test` profiline daha yakından bakalım:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

Bu, pipeline'ı `-profile test` ile her çalıştırdığımızda Nextflow'un yapılandırmayı `conf/test.config` dosyasından yükleyeceği anlamına gelir.

```groovy title="molkart/conf/test.config (excerpt)"
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    input = 'https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv'
    mindagap_tilesize = 90
    mindagap_boxsize = 7
    mindagap_loopnum = 100
    clahe_pyramid_tile = 368
    segmentation_method = "mesmer,cellpose,stardist"
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

Bu profilin daha önce `params.yaml` dosyamızda kullandığımız parametrelerin aynısını içerdiğine dikkat edin.

Birden fazla profili virgülle ayırarak etkinleştirebilirsiniz.
Params dosyamıza ihtiyaç duymadan pipeline'ımızı test etmek için bunu kullanalım:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

Bu şunları birleştirir:

- `docker`: Docker konteynerlerini etkinleştir
- `test`: Test veri setini ve parametrelerini kullan

Profiller soldan sağa uygulanır, bu nedenle aynı değerleri ayarlarlarsa sonraki profiller önceki profilleri geçersiz kılar.

### Özet

nf-core pipeline'ları konteynerler, test ve özel ortamlar için yerleşik profillerle birlikte gelir.
İhtiyacınız olan yapılandırmayı oluşturmak için birden fazla profili birleştirebilirsiniz.

### Sırada ne var?

Farklı hesaplama ortamları için kendi özel profillerinizi nasıl oluşturacağınızı öğrenin.

---

## 3. Özel profiller oluşturma

### 3.1. Yerel geliştirme ve HPC'de yürütme arasında geçiş için profiller oluşturma

İki senaryo için özel profiller oluşturalım:

1. Docker ile yerel geliştirme
2. Slurm zamanlayıcı ve Singularity ile üniversite HPC'si

`nextflow.config` dosyanıza aşağıdakileri ekleyin:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'standard_queue'
        singularity.cacheDir = '/shared/containers'
    }
}
```

Artık ortamlar arasında kolayca geçiş yapabilirsiniz:

```bash
# Yerel geliştirme için
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# HPC için (mevcut olduğunda)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! Note "Not"

    Slurm zamanlayıcısına erişimimiz olmadığı için bu eğitim ortamında HPC profilini test edemeyiz.
    Ancak bu, gerçek dünya kullanımı için nasıl yapılandıracağınızı gösterir.

### 3.2. Yapılandırmayı incelemek için `nextflow config` kullanma

`nextflow config` komutu, pipeline'ı çalıştırmadan tamamen çözümlenmiş yapılandırmayı gösterir.

Varsayılan yapılandırmayı görüntüleyin:

```bash
nextflow config ./molkart
```

Belirli bir profille yapılandırmayı görüntüleyin:

```bash
nextflow config -profile local_dev ./molkart
```

Bu şunlar için son derece kullanışlıdır:

- Yapılandırma sorunlarını giderme
- Hangi değerlerin gerçekten kullanılacağını anlama
- Birden fazla profilin nasıl etkileşime girdiğini kontrol etme

### Özet

Özel profiller, tek bir komut satırı bayrağıyla farklı hesaplama ortamları arasında geçiş yapmanıza olanak tanır.
Çalıştırmadan önce çözümlenmiş yapılandırmayı incelemek için `nextflow config` kullanın.

### Sırada ne var?

nf-core'un süreç etiketi sistemini kullanarak bireysel süreçler için kaynak taleplerini nasıl özelleştireceğinizi öğrenin.

---

## 4. Kaynak taleplerini özelleştirme

### 4.1. nf-core pipeline'larında süreç etiketlerini anlama

Basitlik için nf-core pipeline'ları, tüm pipeline'larda kaynak tahsisini standartlaştırmak için [**süreç etiketlerini**](https://www.nextflow.io/docs/latest/reference/process.html#process-label) kullanır.
Her süreç, sırasıyla düşük, orta veya yüksek hesaplama kaynağı gereksinimlerini tanımlamak için `process_low`, `process_medium` veya `process_high` gibi bir etiketle işaretlenir.
Bu etiketler, pipeline'ın `conf/` dizininde bulunan yapılandırma dosyalarından birinde belirli kaynak taleplerine dönüştürülür.

```groovy title="molkart/conf/base.config (excerpt)"
process {
    cpus   = { 1      * task.attempt }
    memory = { 6.GB   * task.attempt }
    time   = { 4.h    * task.attempt }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_single {
        cpus   = { 1                   }
        memory = { 6.GB * task.attempt }
        time   = { 4.h  * task.attempt }
    }
    withLabel:process_low {
        cpus   = { 2     * task.attempt }
        memory = { 12.GB * task.attempt }
        time   = { 4.h   * task.attempt }
    }
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 8.h   * task.attempt }
    }
    withLabel:process_high {
        cpus   = { 12    * task.attempt }
        memory = { 72.GB * task.attempt }
        time   = { 16.h  * task.attempt }
    }
}
```

`task.attempt` çarpanına dikkat edin - bu, pipeline `process.maxRetries > 1` ile ayarlanmışsa, sonraki görev yeniden denemelerinin daha fazla kaynak talep etmesine olanak tanır.

### 4.2. Belirli süreçler için kaynakları geçersiz kılma

İnce ayarlı kontrol için, bireysel süreçleri isme göre hedefleyin:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

Bu pipeline'ı yukarıdaki geçersiz kılma ile çalıştırmayı denersek, `CELLPOSE` süreci etiketi tarafından tanımlanan varsayılan yerine 16 CPU ve 32 GB bellek talep edecektir.
Bu, mevcut ortamımızda o kadar RAM'imiz olmadığı için pipeline'ın başarısız olmasına neden olacaktır.
Bir sonraki bölümde bu tür başarısızlıkları nasıl önleyeceğimizi öğreneceğiz.

!!! Tip "İpucu"

    Süreç adlarını bulmak için pipeline yürütme çıktısına bakın veya `.nextflow.log` dosyasını kontrol edin.
    Süreç adları `WORKFLOW:SUBWORKFLOW:PROCESS` desenini takip eder.

### Özet

nf-core pipeline'ları kaynak tahsisini standartlaştırmak için süreç etiketlerini kullanır.
Kaynakları etikete göre (birden fazla süreci etkiler) veya isme göre (belirli bir süreci etkiler) geçersiz kılabilirsiniz.

### Sırada ne var?

GitHub Codespaces gibi kısıtlı ortamlarda kaynak limitlerini nasıl yöneteceğinizi öğrenin.

---

## 5. Kısıtlı ortamlarda kaynakları yönetme

### 5.1. Kaynak limitleri sorunu

Molkart'ı 16 CPU ve 32 GB bellek talep eden bir süreçle çalıştırmayı denersek (bölüm 4.2'de gösterildiği gibi), mevcut ortamımızda o kadar kaynağımız olmadığı için başarısız olurdu.
Daha büyük düğümlere sahip bir küme ortamında, bu tür talepler zamanlayıcıya gönderilirdi.

GitHub Codespaces gibi kısıtlı ortamlarda, limitler olmadan Nextflow, mevcut kaynakları aşan süreçleri çalıştırmayı reddederdi.

### 5.2. Kaynak limitlerini ayarlama

`resourceLimits` yönergesi, kaynak taleplerini belirtilen değerlerde sınırlar:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

Bu Nextflow'a şunu söyler: "Herhangi bir süreç 2 CPU'dan veya 7 GB bellekten fazlasını talep ederse, bunun yerine bu limitlere sınırla."

### 5.3. Özel profillere kaynak limitleri ekleme

Özel profillerinizi uygun limitler içerecek şekilde güncelleyin:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
        process.resourceLimits = [
            cpus: 2,
            memory: 7.GB
        ]
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'batch'
        process.resourceLimits = [
            cpus: 32,
            memory: 128.GB,
            time: 24.h
        ]
    }
}
```

!!! Warning "Uyarı"

    Kaynak limitlerini çok düşük ayarlamak süreçlerin başarısız olmasına veya yavaş çalışmasına neden olabilir.
    Pipeline'ın daha az bellek yoğun algoritmalar kullanması veya verileri daha küçük parçalar halinde işlemesi gerekebilir.

### Özet

Süreç kaynak taleplerini sınırlayarak kaynak kısıtlı ortamlarda pipeline'ları çalıştırmak için `resourceLimits` kullanın.
Farklı profiller, ortamlarına uygun farklı limitlere sahip olabilir.

### Sırada ne var?

Temel Nextflow for Bioimaging eğitimini tamamladınız!

---

## Sonuç

Artık Nextflow pipeline'larını farklı hesaplama ortamları için nasıl yapılandıracağınızı anlıyorsunuz.

Öğrendiğiniz temel beceriler:

- **Yapılandırma önceliği**: Nextflow'un ayarları birden fazla kaynaktan nasıl çözdüğü
- **nf-core profilleri**: Konteynerler, test ve yardımcı programlar için yerleşik profilleri kullanma
- **Özel profiller**: Farklı ortamlar için kendi profillerinizi oluşturma
- **Süreç etiketleri**: Etikete göre kaynak taleplerini anlama ve geçersiz kılma
- **Kaynak limitleri**: `resourceLimits` ile kısıtlı ortamları yönetme
- **Yapılandırma incelemesi**: Ayarları hata ayıklamak ve doğrulamak için `nextflow config` kullanma

Bu yapılandırma becerileri herhangi bir Nextflow pipeline'ına aktarılabilir ve iş akışlarını yerel makineler, HPC kümeleri ve bulut platformlarında verimli bir şekilde çalıştırmanıza yardımcı olacaktır.

### Sırada ne var?

Nextflow for Bioimaging kursunu tamamladığınız için tebrikler!

Sonraki adımlar:

- Geri bildirim sağlamak için kurs anketini doldurun
- İş akışları geliştirme hakkında daha fazla bilgi edinmek için [Hello Nextflow](../hello_nextflow/index.md) bölümüne göz atın
- nf-core araçlarına daha derinlemesine dalmak için [Hello nf-core](../hello_nf-core/index.md) bölümünü keşfedin
- [Eğitim koleksiyonlarındaki](../training_collections/index.md) diğer kurslara göz atın
