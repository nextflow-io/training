# Bölüm 2: Pipeline Yürütmesini Yapılandırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[Bölüm 1](./01_run_demo.md)'de nf-core/demo pipeline'ını test profili kullanarak buldunuz ve çalıştırdınız.
Şimdi pipeline yürütmesinin nasıl yapılandırılacağına bakıyoruz: parametre ayarlama, doğrulamayı anlama ve kaynak tahsisi ile araç argümanlarını özelleştirme.

[Hello Config](../hello_nextflow/06_hello_config.md)'da açıklandığı gibi, pipeline kodunu değiştirmeden pipeline'ımızın hangi veriler üzerinde ve nasıl çalışacağını değiştirebilmek istiyoruz.
Bu amaçla Nextflow, pipeline yapılandırmasını kontrol etmenin birden fazla yolunu destekler; bu durum başlangıçta biraz bunaltıcı gelebilir.

nf-core projesi, yapılandırma öğelerini düzenlemek için kurallar belirler ve üst düzeyde iki tür yapılandırmayı birbirinden ayırır: **pipeline parametreleri** ve dar anlamda **yapılandırma**.

- **Pipeline parametreleri** (`params` sistemi aracılığıyla ayarlanır) genellikle girdi dosyaları, araç davranış bayrakları ve analiz parametrelerini içerir.
- Dar anlamda **yapılandırma**, pipeline'ın nasıl çalıştırıldığına ilişkin lojistiği ifade eder; yani yürütücü, hesaplama kaynağı tahsisleri ve benzerleri.

<figure class="excalidraw">
    --8<-- "docs/en/docs/nfcore_use/img/params_vs_config.excalidraw.svg"
</figure>

Önce pipeline parametrelerini ele alarak başlayalım, ardından dar anlamda yapılandırmaya bakacağız.

---

## 1. Pipeline Parametreleri

Tüm nf-core pipeline'ları için, `--help` bayrağını kullanarak doğrudan komut satırından pipeline parametrelerinin tam listesini alabilirsiniz; bu bayrağın kendisi de bir pipeline parametresidir.

### 1.1. `--help` ile Parametre Listesini Alma

Demo pipeline için yardım komutunu çalıştırın:

```bash
nextflow run nf-core/demo --help
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [adoring_meucci] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>


    Input/output options
      --input                       [string] Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string] The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
      --email                       [string] Email address for completion summary.
      --multiqc_title               [string] MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string] Name of iGenomes reference.
      --fasta                       [string] Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean] Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]          Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string] Display the help message.
      --help_full                   [boolean]         Display the full detailed help message.
      --show_hidden                 [boolean]         Display hidden parameters in the help message (only works when --help or --help_full are provided).
     !! Hiding 19 param(s), use the `--showHidden` parameter to show them !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

Gördüğünüz gibi çıktı, parametreleri kategorilere göre gruplandırır (Girdi/çıktı seçenekleri, Referans genom seçenekleri vb.) ve her biri için tür ile açıklama bilgisi sunar.

Bu kategorilendirme, aşağıda daha ayrıntılı ele alınan bir şema dosyası tarafından belirlenir.
Yalın Nextflow pipeline'larında `--help`, yalnızca geliştirici tarafından manuel olarak uygulanmışsa çalışır.

!!! tip "İpucu"

    Varsayılan olarak gizlenen `--publish_dir_mode` veya `--monochrome_logs` gibi ek parametreleri görmek için `--help --show_hidden` kullanın.

### 1.2. Parametre Değerlerini Ayarlama

[Hello Config](../hello_nextflow/06_hello_config.md)'da ele alındığı gibi, parametre değerlerini komut satırında `--param_name` ile ayarlayabilir ya da bir dizi parametreyi YAML dosyasında toplayıp `-params-file` ile iletebilirsiniz.
Her iki yaklaşım da nf-core pipeline'larında aynı şekilde çalışır.

Örneğin, kırpma adımını atlamak için `skip_trim` boolean parametresini `true` olarak ayarlamak istiyoruz.
Çalışma dizininizde bu değer önceden ayarlanmış `my_params.yml` adlı bir params dosyası bulunmaktadır:

```yaml title="my_params.yml"
skip_trim: true
```

`-params-file` ile iletin:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-notrim -params-file my_params.yml
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [focused_heisenberg] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results-notrim

    Process skipping options
      skip_trim                 : true

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_22-08-47

    Core Nextflow options
      revision                  : master
      runName                   : focused_heisenberg
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               | 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

`SEQTK_TRIM` süreci artık çıktıda görünmüyor.

!!! warning "Uyarı: Parametre girdileri hakkında önemli kısıtlamalar"

    **Komut satırında boolean parametreleri ayarlama**

    Nextflow 26.04 sürümünden itibaren, komut satırında sağlanan tüm değerler string olarak yazılır.
    `skip_trim` gibi bir boolean parametre için, bunu yalın bir bayrak olarak (`--skip_trim`) veya `--skip_trim true` şeklinde geçirmek, **string** `"true"` olarak değerlendirilir ve şema doğrulamasında başarısız olur:

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    Bir boolean parametreyi gerçek `true`/`false` değerine ayarlamak için yukarıda gösterildiği gibi `-params-file` kullanın ya da bir yapılandırma dosyasında ayarlayın.
    String, integer ve dosya yolu parametreleri bu durumdan etkilenmez ve doğrudan komut satırında ayarlanabilir.
    Bu kurs, boolean parametreler için bu kalıbı her yerde kullanır.

    **Özel yapılandırma dosyaları kullanma**

    `-c` ile iletilen özel bir yapılandırma dosyasında pipeline parametrelerini ayarlamak teknik olarak mümkün olsa da, Nextflow'un yapılandırma öncelik kurallarına bağlı olarak bu, pipeline'ın kendi `nextflow.config` dosyasında önceden ayarlanmış varsayılanları geçersiz kılmayabilir.
    Komut satırında `--param_name` veya `-params-file` kullanmak daha güvenilirdir; çünkü bunlar her zaman öncelik taşır.

    Pratik bir kural olarak: `--help` çıktısında görünüyorsa, bir yapılandırma dosyası yerine komut satırı veya params dosyası aracılığıyla ayarlayın.

### 1.3. Parametre Doğrulama

İlginç bir bilgi: `--help` komutu tüm nf-core pipeline'larında çalışır; çünkü nf-core projesi, geliştiricilerin tüm pipeline parametrelerini bir JSON şema dosyasında (`nextflow_schema.json`) resmi olarak tanımlamasını zorunlu kılar.
Bu şema, her parametrenin türünü, açıklamasını, varsayılan değerini ve grubunu kaydeder.

`--help` çıktısını desteklemenin yanı sıra, şema dosyası başlatma sırasında otomatik doğrulamayı da mümkün kılar.
Bu, Nextflow'un ilettiğiniz her parametrenin var olup olmadığını ve uygun bir değer verilip verilmediğini (uygun türde, izin verilen değer aralığında vb.) kontrol edebildiği anlamına gelir.

Bunu [girdi doğrulama bölümünde](../nfcore_build/04_input_validation.md) daha ayrıntılı ele alıyoruz; ancak demo pipeline'a geçersiz parametre girdisi vererek bunu şimdiden uygulamada görebilirsiniz.

#### 1.3.1. Tanınmayan Parametreler

Var olmayan bir parametre geçirmeyi deneyin:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --foobar "invalid"
```

Konsol çıktısı bir uyarı içerir:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

Pipeline çalışmaya devam eder; ancak uyarı, `--foobar`'ın tanınan bir parametre olmadığını hemen bildirir.
Bu, `--outdir` yerine `--outDir` kullanmak gibi kritik olmayan yazım hatalarına dikkatinizi çekmek için tasarlanmıştır; bu sayede zaman ve hesaplama kaynağı israfından kaçınabilirsiniz.

#### 1.3.2. Geçersiz Parametre Değerleri

Doğrulama, parametre **değerlerini** de kontrol eder.
`--skip_trim` parametresi bir boolean bayraktır; bu nedenle string değer geçirmek pipeline'ın hemen başarısız olmasına neden olur:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

Pipeline, herhangi bir süreç çalışmadan önce durur; böylece başarısız veya hatalı bir yürütmeden kaçınmış olursunuz.
[1.2](#12-set-parameter-values) bölümünde belirtildiği gibi, komut satırı değerleri string olarak yazıldığından boolean parametreler, komut satırında geçirilmek yerine params dosyasında gerçek `true`/`false` değerine ayarlanmalıdır.

### 1.4. Girdi Doğrulama

Aynı doğrulama mantığı, girdi dosyalarının geçerliliğini kontrol etmek için de kullanılabilir.
Örneğin, bir pipeline ana veri girdisi olarak bir samplesheet bekliyorsa (ki bu durum pek çok nf-core pipeline'ında geçerlidir), geliştirici girdi dosyasının nasıl yapılandırılması gerektiğini açıklayan bir girdi şeması (parametre şemasından ayrı) sağlayabilir.

Ardından Nextflow, çalışma zamanında sağlanan girdi dosyasının geçerli olup olmadığını kontrol edebilir.

Bunu da [girdi doğrulama bölümünde](../nfcore_build/04_input_validation.md) daha ayrıntılı ele alıyoruz; ancak demo pipeline'a geçersiz bir girdi samplesheet'i vererek bunu şimdiden uygulamada görebilirsiniz.

`nf-core/demo` pipeline'ı, `sample`, `fastq_1` ve `fastq_2` sütunlarına sahip bir CSV dosyası bekler.
Bu, beklenen yapıyı, sütun türlerini ve kısıtlamaları belirten bir şema dosyasında (`assets/schema_input.json`) tanımlanmıştır.

??? abstract "Girdiler için şema dosyası"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

Şema, `sample` ve `fastq_1`'in zorunlu olduğunu, `fastq_2`'nin ise isteğe bağlı olduğunu belirtir (hem çift uçlu hem de tek uçlu veriyi destekler).
Dosya yolları, varlık ve uzantı kalıbı açısından doğrulanır.

Bunu göstermek için çalışma dizininizde `malformed_samplesheet.csv` adlı hatalı biçimlendirilmiş bir samplesheet sağlıyoruz:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

Bu samplesheet'te zorunlu `fastq_1` sütunu eksik ve `fastq_2`'de var olmayan bir dosya yolu bulunuyor.

Demo pipeline'ı `malformed_samplesheet.csv` dosyasını girdi olarak kullanarak çalıştırın:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

Gördüğünüz gibi pipeline hemen başarısız olur ve **tüm** doğrulama hatalarını bir arada raporlar.
nf-schema ilk hatada durmaz; her sorunu toplar ve birlikte listeler; böylece sorunları tek tek keşfetmek yerine hepsini bir seferde düzeltebilirsiniz.

Her hata, soruna neden olan tam girdiyi ve alanı tanımlar; bu sayede samplesheet'inizi düzelttikten sonra pipeline'ı, Nextflow'un dosya yoluna gerçekten erişmeye çalıştığı ileri bir noktada başarısız olmayacağından emin olarak yeniden başlatabilirsiniz.

Geliştiriciler için tüm bunlar, [nf-core ile Geliştirme'nin 4. Bölümünde](../nfcore_build/04_input_validation.md) daha ayrıntılı ele alınmaktadır.

### Özetle

Bir pipeline'ın parametrelerinin tam listesini `--help` ile nasıl alacağınızı, bunları komut satırı veya params dosyası aracılığıyla nasıl ayarlayacağınızı ve Nextflow'un hem parametre değerlerini hem de girdi dosyalarını pipeline şemalarına göre nasıl doğruladığını öğrendiniz.

### Sırada ne var?

Diğer yapılandırma türü hakkında bilgi edinin: kaynak tahsisi ve araç argümanlarını kapsayan, pipeline'ın nasıl çalıştığına ilişkin yapılandırma.

---

## 2. Yapılandırma

Dar anlamda yapılandırma, pipeline'ın **nasıl** çalıştığını kontrol eder: kaynak tahsisi, araca özgü argümanlar, görevlerin nerede yürütüleceği ve hangi yazılım paketleme sisteminin kullanılacağı.

nf-core pipeline'ları, `nextflow.config` ve `conf/` dizininde varsayılan yapılandırmayı içerir.
Herhangi bir şeyi geçersiz kılmadan önce, varsayılanların nerede bulunduğunu bilmek faydalıdır.

### 2.1. Yapılandırma Dosyalarını İnceleme

[Bölüm 1](./01_run_demo.md)'de pipeline kaynak kodunun `$NXF_HOME/assets` altında bulunduğunu gördünüz.
[Bölüm 1](./01_run_demo.md)'de oluşturduğunuz `pipelines` sembolik bağlantısını kullanarak, mevcut yapılandırma dosyalarını listelemek için şunu çalıştırın:

```bash
ls pipelines/nf-core/demo/conf/
```

```console
base.config
containers_conda_lock_files_amd64.config
containers_conda_lock_files_arm64.config
containers_docker_amd64.config
containers_docker_arm64.config
containers_singularity_https_amd64.config
containers_singularity_https_arm64.config
containers_singularity_oras_amd64.config
containers_singularity_oras_arm64.config
igenomes.config
igenomes_ignored.config
modules.config
test.config
test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/nfcore_use/img/nfcore_config_files.excalidraw.svg"
</figure>

En önemli yapılandırma dosyaları şunlardır:

- **`conf/base.config`**: Süreçlere CPU, bellek ve süre atayan kaynak etiketlerini (`process_low`, `process_medium`, `process_high`) tanımlar. Bir sürecin beklenenden fazla kaynak kullandığını gördüğünüzde, bu varsayılanlar buradan gelir.
- **`conf/modules.config`**: Süreç başına araç argümanlarını (`ext.args`) ve çıktı yayımlama ayarlarını (`publishDir`) belirler. Her aracın varsayılan olarak hangi argümanları aldığını görmek için bu dosyayı açın.
- **`conf/test.config`**: [Bölüm 1](./01_run_demo.md)'de kullandığınız test profili; `resourceLimits` aracılığıyla kaynakları sınırlar ve bir test samplesheet'i ayarlar. `-profile test` ile etkinleştirilir.
  Tam boyutlu bir test veri kümesiyle çalıştırmak için, kıyaslama açısından kullanışlı olan `conf/test_full.config` da mevcuttur.

Merkezi `nextflow.config`, yukarıdakilerin tümünü yükler ve her şey için uygun varsayılan değerleri ayarlar.

Bu dosyalarda belirtilen ayarlardan herhangi birini değiştirmek istiyorsanız, bu dosyaların hiçbirini doğrudan değiştirmeyin.
Bunun yerine kendi yapılandırma dosyanızı oluşturun ve `-c` ile iletin.
Belirttiğiniz değerler, diğer dosyalarda ayarlanan varsayılan değerleri geçersiz kılar.

Bunu pratikte deneyelim.

### 2.2. Süreç Kaynaklarını ve Araç Argümanlarını Özelleştirme

nf-core modülleri iki yaygın yapılandırma geçersiz kılma türünü destekler: **kaynak tahsisi** (CPU, bellek, süre) ve `ext.args` aracılığıyla **araç argümanları**.

Pek çok komut satırı aracının, pipeline parametresi olarak sunulacak kadar yaygın kullanılmayan argümanları vardır.
`ext.args` kuralı, bu argümanları altta yatan araca pipeline parametresi yerine bir yapılandırma dosyası aracılığıyla iletmenizi sağlar.

Çalışma dizininizde sağlanan `custom.config` dosyası her iki geçersiz kılmayı da göstermektedir:

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

İlk blok, `FASTQC` kaynak tahsisini geçersiz kılar.
Varsayılan olarak `FASTQC`, `base.config`'deki `process_medium` etiketini kullanır ve 6 CPU ile 36 GB bellek tahsis eder; burada ise 2 CPU ve 4 GB ile sınırlandırıyoruz.

İkinci blok, `ext.args` aracılığıyla `SEQTK_TRIM`'e ekstra bir argüman iletir.
`-b 5` bayrağı, `seqtk trimfq`'ya kalite kırpmanın yanı sıra her okumanın başından 5 baz kırpmasını söyler.

Pipeline'ı bu yapılandırmayla çalıştırın:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-custom -c custom.config
```

??? success "Komut çıktısı"

    ```console
    executor >  local (8)
    [95/b32876] NFCORE_DEMO:DEMO:FASTQC (SAMPLE1_PE)     | 3 of 3 ✔
    [17/428668] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE1_PE) | 3 of 3 ✔
    [cf/85991a] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [3c/94a7a0] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

`-c` bayrağı, yapılandırmanızı pipeline'ın yerleşik yapılandırmasının üzerine ekler.

`ext.args` geçersiz kılmasının etkili olduğunu doğrulamak için, çalıştırma çıktısından `SEQTK_TRIM` work dizini hash'ini bulun (örn. `work/17/428668...`) ve içindeki `.command.sh` dosyasını kontrol edin:

```bash
cat work/17/428668/.command.sh
```

??? success "Komut çıktısı"

    ```console
    #!/usr/bin/env bash -e -u -o pipefail
    printf "%s\n" sample1_R1.fastq.gz sample1_R2.fastq.gz | while read f;
    do
        seqtk \
            trimfq \
            -b 5 \
            $f \
            | gzip --no-name > SAMPLE1_PE_$(basename $f)
    done
    ...
    ```

`seqtk trimfq` komutunda `-b 5`'i görmelisiniz.

`ext.args` hakkında bilmeniz gereken önemli bir nokta: bir modülün zaten ayarlanmış varsayılan bir değeri varsa, değeriniz ona eklenmek yerine onu **tamamen değiştirir**.
Örneğin, `FASTQC`'nin `conf/modules.config`'de varsayılan olarak `ext.args = '--quiet'` ayarı bulunur:

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args   = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}",
        ]
    }
```

`FASTQC` için `ext.args = '--kmers 8'` ayarlarsanız, `--quiet` bayrağı artık uygulanmaz.
Her ikisini de korumak için `ext.args = '--quiet --kmers 8'` olarak ayarlayın.

`ext.args`'ı geçersiz kılmadan önce her zaman bir modülün varsayılan yapılandırmasını kontrol etmelisiniz.

### Özetle

nf-core pipeline yapılandırma varsayılanlarının nerede bulunduğunu ve özel bir yapılandırma dosyasıyla kaynak tahsislerini ve araç argümanlarını nasıl geçersiz kılacağınızı öğrendiniz.

### Sırada ne var?

Öğrendiklerinizi gerçek bir üretim pipeline'ına uygulayacağınız [Bölüm 3](./03_run_production_pipeline.md)'e geçin.

---

## Özet

Bu bölümde şunları öğrendiniz:

- Yardım alma, parametre ayarlama ve parametre ile girdi doğrulamasını anlama
- Yapılandırma dosyaları aracılığıyla kaynak tahsisini ve araç argümanlarını özelleştirme
