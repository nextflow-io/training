```markdown
# Bölüm 1: Hesaplama Ortamınıza Uyarlama

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[Nextflow Run](../nextflow_run/index.md) kursunda bir pipeline'ın girdilerini, parametrelerini ve çıktılarını yapılandırdınız.
Bu kurs, tablonun diğer yarısını ele alır: iş akışı kodunu değiştirmeden, bir pipeline'ın yürütülmesini çalıştığı hesaplama ortamına uyarlamayı öğreneceksiniz.

!!! example "Senaryo"

    Pipeline'ınızı Docker kullanarak dizüstü bilgisayarınızda geliştirip test ettiniz.
    Şimdi onu başkalarına devretmeniz gerekiyor: bir iş arkadaşınızın yalnızca Conda kurulu, kurumunuzun HPC kümesi ise işlerin kendi zamanlayıcısı ve kaynak sınırları aracılığıyla gönderilmesini bekliyor.
    Bunların hiçbiri pipeline'ın kendisini yeniden yazmayı gerektirmemelidir.

Aynı pipeline kodu tüm bu ortamlarda çalışabilir; çünkü bunların hiçbiri iş akışına sabit olarak kodlanmamıştır.
Yazılım paketleme, yürütme platformu ve kaynak tahsisi; kodun üzerine katmanlanan yapılandırma aracılığıyla kontrol edilir. Bu kurs tam olarak bunu ele alır: aynı pipeline'ı kod değil, yapılandırma değiştirerek yeni bir ortama nasıl uyarlayacağınızı öğreneceksiniz.

---

## 1. Yazılım Paketleme Teknolojisi Seçimi

[Nextflow Run](../nextflow_run/index.md) kursunda, `nextflow.config` dosyasında Docker'a alternatif olarak önceden ayarlanmış bir `conda` profili gördünüz.
Burada aynı geçişi kendiniz oluşturacak ve bir sürecin Conda ile gerçekten kullanılabilir hale gelmesi için neler gerektiğini göreceksiniz.

### 1.1. Docker'ı Devre Dışı Bırakma ve Conda'yı Etkinleştirme

`docker.enabled` değerini `false` olarak değiştirin ve Conda'yı etkinleştiren bir yönerge ekleyin.

=== "Sonra"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Bu ayar, Nextflow'un Conda paketi belirtilmiş her süreç için Conda ortamları oluşturmasına ve kullanmasına olanak tanır.
`cowpy` sürecinin henüz bir Conda paketi yok; şimdi bunu tamamen yapılandırma üzerinden ekleyelim.

### 1.2. Yapılandırma Aracılığıyla Conda Paketi Ekleme

`conda` yönergesi, `modules/cowpy.nf` dosyasında `container` yönergesinin zaten bulunduğu gibi, doğrudan süreç tanımına eklenebilir. Ancak bunu yapmak zorunda değilsiniz: `withName` kullanarak bu yönergeyi yapılandırma dosyasından, yalnızca `cowpy` süreciyle sınırlı olacak şekilde ayarlayabilirsiniz.

=== "Sonra"

    ```groovy title="nextflow.config" linenums="6" hl_lines="3-5"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
    }
    ```

Bu ayar, pipeline kodunda zaten bulunan `container` yönergesinin yerini almaz; o koda hiç dokunmadan yanına bir alternatif ekler.

!!! tip "İpucu"

    [Seqera Containers](https://seqera.io/containers/) araması, bir araç için Conda paketi URI'sini bulmak amacıyla kullanışlı bir yoldur; bu araçtan konteyner oluşturmayı planlamasanız bile.

### 1.3. Conda Kullanımını Doğrulamak İçin İş Akışını Çalıştırma

```bash
nextflow run main.nf --batch conda
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_pike] revision: c3c85dec78

    executor >  local (8)
    [6d/d48030] sayHello (2)       | 3 of 3 ✔
    [f5/7a9d76] convertToUpper (1) | 3 of 3 ✔
    [1c/79b693] collectGreetings   | 1 of 1 ✔
    Creating env using conda: conda-forge::cowpy==1.1.5 [cache /workspaces/training/execution-config/work/conda/env-898314d566668b6587ad714ae06b8520]
    [bb/64b67c] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-conda-output.txt

      batch_report: full_pipeline/conda-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-conda-output.txt
    ```

Bu işlem, arka planda farklı mekanizmalar kullanılsa da Docker ile çalıştırmakla aynı çıktıyı üretir: Nextflow bir konteyner imajı çekmek yerine Conda paketini alır ve ondan bir ortam oluşturur.

!!! info "Bilgi"

    Yeni bir Conda ortamı oluşturmak, ilk seferinde bir konteyner çekmekten biraz daha uzun sürebilir; ancak burada kullanılan paket küçük olduğundan işlem hızlı tamamlanacaktır.

Şimdi bu kursun geri kalanı için Docker'a geri dönün.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

??? tip "Docker ve Conda'yı Birlikte Kullanma"

    Bu ayarlar süreç bazında kapsamlandırıldığından, Docker ve Conda'yı karıştırabilirsiniz: bazı süreçler Docker, diğerleri Conda kullanır; bu tercih her araç için neyin mevcut olduğuna bağlıdır.
    Aynı süreç için hem `container` yönergesi (pipeline kodunda) hem de `conda` yönergesi (burada, yapılandırmadan) ayarlanmışsa ve her iki paketleme sistemi de etkinleştirilmişse, Nextflow konteynerlere öncelik verir.

### Özetle

Bir sürecin hangi yazılım paketleme teknolojisini kullanacağını nasıl yapılandıracağınızı ve Docker ile Conda arasında nasıl geçiş yapacağınızı öğrendiniz.

### Sırada ne var?

Nextflow'un görevlerinizi gerçekten çalıştırmak için kullandığı yürütme platformunu nasıl değiştireceğinizi öğrenin.

---

## 2. Yürütme Platformu Seçimi

Şimdiye kadar çalıştırdığınız her pipeline, yerel yürütücüyü kullandı: her görev, Nextflow'un kendisiyle aynı makinede çalışır.
Nextflow mevcut CPU ve bellek kaynaklarını kontrol eder; yeterli kaynak serbest kalana kadar görevleri bekletir.

Yerel yürütücü kullanışlıdır; ancak tek bir makinenin ötesine ölçeklenemez.
Nextflow, HPC zamanlayıcıları (Slurm, LSF, SGE, PBS ve diğerleri) ile bulut platformları (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes ve daha fazlası) dahil olmak üzere [pek çok farklı yürütme arka ucunu](https://nextflow.io/docs/latest/executor.html) destekler.

### 2.1. Farklı Bir Arka Uç Hedefleme

Yürütücü, `executor` adlı bir süreç yönergesiyle ayarlanır.
Varsayılan değeri `local` olduğundan, aşağıdaki yapılandırma örtük olarak geçerlidir:

```groovy title="Yerleşik yapılandırma"
process {
    executor = 'local'
}
```

Farklı bir arka uç hedeflemek için yönergeyi istediğiniz yürütücüye ayarlayın.

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Uyarı"

    Eğitim ortamı bir HPC kümesine bağlı değildir; dolayısıyla bunu burada çalıştırmanız mümkün değildir.

### 2.2. Arka Uca Özgü Sözdizimi Soyutlanmıştır

Çoğu HPC platformu, iş gönderimlerinde CPU, bellek ve kuyruk adı gibi kaynak taleplerinin kendi sözdizimiyle belirtilmesini gerektirir.
`my-science-work` adlı bir kuyrukta 8 CPU ve 4 GB RAM için yapılan aynı talep, zamanlayıcıya göre tamamen farklı görünür.

??? abstract "Örnekler"

    ```bash title="SLURM yapılandırması / sbatch ile gönderim"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="PBS yapılandırması / qsub ile gönderim"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=8
    #PBS -l mem=4gb
    ```

    ```bash title="SGE yapılandırması / qsub ile gönderim"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=8
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Nextflow tüm bunları soyutlar: `cpus`, `memory` ve `queue` gibi standartlaştırılmış özellikleri bir kez belirtirsiniz (tam liste için [süreç yönergelerine](https://nextflow.io/docs/latest/reference/process.html#process-directives) bakın); Nextflow bunları çalışma zamanında uygun arka uca özgü betiklere dönüştürür.

### 2.3. Nextflow'un Gerçekte Ne Çalıştırdığını İnceleme

Bu dönüşüm yalnızca bir yapılandırma dosyası kolaylığı değildir: yerel yürütücüyle bile şu an inceleyebileceğiniz somut bir şeyle desteklenir.
[Nextflow Run, bölüm 1.3](../nextflow_run/01_run_nextflow.md#13-explore-the-work-directory)'te `work/` dizini altındaki bir görev dizininin içine baktınız ve Nextflow'un çalıştırdığı tam komutu içeren `.command.sh` dosyasını buldunuz.
Aynı dizinde henüz incelemediğiniz bir dosya daha bulunur: `.command.run`.

```bash
cat work/0a/0df4a1*/.command.run
```

??? success "Komut çıktısı (alıntı)"

    ```console
    #!/bin/bash
    ### ---
    ### name: 'convertToUpper (3)'
    ### container: 'null'
    ### outputs:
    ### - 'UPPER-Bonjour-output.txt'
    ### ...
    set -e
    set -u
    ...
    nxf_launch() {
        /bin/bash -ue /workspaces/training/nextflow-run/work/0a/0df4a1028c2001758b1841cff92fc7/.command.sh
    }
    ...
    ```

`.command.run`, Nextflow'un yürütme için gerçekten ilettiği betiktir.
`.command.sh` dosyasını, onu gerçekten çalıştırmak için gereken her şeyle sarar: ortam kurulumu, girdi/çıktı hazırlama ve sonucu Nextflow'a geri bildirme.
`local` yürütücüyle Nextflow bu betiği aynı makinede doğrudan çalıştırır.

Farklı bir `executor` ayarladığınızda değişen tam olarak budur.
Slurm veya PBS gibi bir HPC zamanlayıcısı için Nextflow aynı türde sarmalayıcı betiği oluşturur; [2.2](#22-backend-specific-syntax-is-abstracted-away) bölümünde gördüğünüz zamanlayıcıya özgü başlığı (`cpus`, `memory` ve `queue` ayarlarınızdan dönüştürülmüş olarak) ekler ve sonucu o zamanlayıcının kendi gönderim komutuna, örneğin Slurm için `sbatch`'e iletir.
Bundan sonra Nextflow, yerel bir süreci doğrudan izlemek yerine iş durumu için zamanlayıcıyı sorgular.
Bulut toplu iş arka uçları, gönderim komutu yerine API çağrılarıyla yönetildiğinden biraz farklı çalışır; ancak temel fikir aynıdır: aynı görev betiği çalışır, yalnızca nasıl başlatıldığı ve izlendiği değişir.

### Özetle

Farklı hesaplama altyapısını hedeflemek için yürütücüyü nasıl değiştireceğinizi, Nextflow'un arka uca özgü gönderim sözdizimini nasıl soyutladığını ve bir görev farklı bir arka uçta çalıştığında arka planda gerçekte neler olduğunu öğrendiniz.

### Sırada ne var?

[Bölüm 2](./02_resources_and_retries.md)'ye geçin; burada hesaplama kaynaklarını nasıl profilleyeceğinizi ve tahsis edeceğinizi, görev başarısızlıklarını ise yeniden denemelerle nasıl ele alacağınızı öğreneceksiniz.

---

## Özet

Bu bölümde şunları öğrendiniz:

- Docker ve Conda arasında yazılım paketleme teknolojisini değiştirme
- Süreç tanımına `conda` yönergesi ekleme
- `executor` yönergesiyle yürütme platformunu değiştirme
- Nextflow'un bir görev için gerçekte ne oluşturduğunu ve çalıştırdığını, bunun yürütücüler arasında nasıl değiştiğini inceleme
```
