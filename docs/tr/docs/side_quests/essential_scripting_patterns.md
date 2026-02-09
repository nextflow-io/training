# Temel Nextflow Kodlama Kalıpları

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow, Java Sanal Makinesi üzerinde çalışan bir programlama dilidir. Nextflow, [Groovy](http://groovy-lang.org/) üzerine kurulu olsa ve sözdiziminin çoğunu onunla paylaşsa da, Nextflow "uzantılı Groovy"den fazlasıdır -- tam olarak belirtilmiş bir [sözdizimi](https://nextflow.io/docs/latest/reference/syntax.html) ve [standart kütüphane](https://nextflow.io/docs/latest/reference/stdlib.html) ile bağımsız bir dildir.

Değişkenler, map'ler ve listeler için temel sözdiziminin ötesine geçmeden çok sayıda Nextflow kodu yazabilirsiniz. Çoğu Nextflow eğitimi iş akışı düzenlemesine (kanallar, process'ler ve veri akışı) odaklanır ve sadece bununla bile şaşırtıcı derecede ileri gidebilirsiniz.

Ancak, veri manipüle etmeniz, karmaşık dosya adlarını ayrıştırmanız, koşullu mantık uygulamanız veya sağlam üretim iş akışları oluşturmanız gerektiğinde, kodunuzun iki ayrı yönünü düşünmek yardımcı olur: **dataflow** (kanallar, operatörler, process'ler ve workflow'lar) ve **scripting** (closure'lar, fonksiyonlar ve process script'leri içindeki kod). Bu ayrım bir nebze keyfi olsa da -hepsi Nextflow kodu- ne zaman pipeline'ınızı düzenleyip ne zaman veriyi manipüle ettiğinizi anlamak için yararlı bir zihinsel model sağlar. Her ikisinde de ustalaşmak, açık ve bakımı yapılabilir iş akışları yazma yeteneğinizi önemli ölçüde artırır.

### Öğrenme hedefleri

Bu yan görev, sizi temel kavramlardan üretime hazır kalıplara kadar uygulamalı bir yolculuğa çıkarır.
Basit bir CSV okuyan iş akışını gerçekçi zorluklarla adım adım geliştirerek karmaşık bir biyoinformatik pipeline'a dönüştüreceğiz:

- **Sınırları anlamak:** Dataflow işlemleri ile scripting arasında ayrım yapmak ve bunların birlikte nasıl çalıştığını anlamak
- **Veri manipülasyonu:** Güçlü operatörler kullanarak map'leri ve koleksiyonları çıkartmak, dönüştürmek ve alt kümelemek
- **String işleme:** Regex kalıpları ile karmaşık dosya adlandırma şemalarını ayrıştırmak ve değişken enterpolasyonunda ustalaşmak
- **Yeniden kullanılabilir fonksiyonlar:** Daha temiz, daha bakımı kolay iş akışları için karmaşık mantığı adlandırılmış fonksiyonlara çıkartmak
- **Dinamik mantık:** Farklı girdi türlerine uyum sağlayan process'ler oluşturmak ve dinamik kaynak tahsisi için closure'lar kullanmak

### Gereksinimleri

Bu öğretici için, şunları bildiğinizi varsayıyoruz:

- Temel Nextflow iş akışları oluşturma (**hello-nextflow** modülümüzü tamamladıysanız buna hazırsınız)
- Process'ler, iş akışları ve kanallar gibi Nextflow temellerinin temel anlayışı
- Temel komut satırı (bash/shell) kavramları

## 1. Başlarken: Basit Analizden Karmaşık Pipeline'lara

Hadi bir CSV dosyasından veri yükleyip çözümleyen ve basit, düz metin bir rapor oluşturan basit bir workflow ile başlayalım. Yolculuğumuz boyunca, bu iş akışını adım adım daha karmaşık biyoinformatik pipeline'a dönüştüreceğiz - gerçek dünya zorluklarıyla her adımda karşılaşacağız.

### 1.1. Temel CSV İşleme

Önce, aşağıdaki dosya yapısı ile bir proje oluşturun:

```bash
mkdir -p essential_patterns/{data,modules}
cd essential_patterns
```

İlk CSV veri dosyasını aşağıdaki içerikle oluşturalım:

```bash
cat > data/samples.csv << 'EOF'
id,condition,read1,read2
sample1,treated,data/sample1_1.fq.gz,data/sample1_2.fq.gz
sample2,control,data/sample2_1.fq.gz,data/sample2_2.fq.gz
sample3,treated,data/sample3_1.fq.gz,data/sample3_2.fq.gz
EOF
```

Şimdi, basit bir Nextflow workflow dosyası oluşturalım - bu dosya CSV'yi okur, basit bir rapor oluşturur ve ek metadata (bu durumda `condition`) toplar:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
    // CSV dosyasını oku ve parse et
    Channel
        .fromPath('data/samples.csv')
        .splitCsv(header: true)
        .map { row ->
            // Meta nesnesini ve girdi dosya yollarını yarat
            meta = [id: row.id, condition: row.condition]
            reads = [row.read1, row.read2]
            tuple(meta, reads)
        }
        .set { ch_samples }

    // Her numune için basit bir rapor üret
    ch_samples.map { meta, reads ->
        log.info "İşleniyor: ${meta.id} (${meta.condition})"
        "Sample ${meta.id} (${meta.condition}) - Okuma dosyaları: ${reads.join(', ')}"
    }
    .collectFile(name: 'report.txt', newLine: true)
    .view { "Rapor dosyası oluşturuldu: $it" }
}
```

Şimdi, bu iş akışını çalıştıralım:

```bash
nextflow run main.nf
```

### 1.2. İş akışını bir Process'e Genişletme

Şimdi, işlerimizi daha kapsamlı ve genişletilebilir kılmak için süreci bir process'e genişletelim. Bu, kodunuzun organizasyonunu geliştirir ve ek process'leri entegre etmeyi kolaylaştırır:

```groovy title="modules/report.nf" linenums="1"
process GENERATE_REPORT {
    publishDir "results", mode: 'copy'

    input:
    tuple val(meta), val(reads)

    output:
    path "${meta.id}_report.txt"

    exec:
    def report = """
    Sample Report
    ============
    ID: ${meta.id}
    Condition: ${meta.condition}
    Read files: ${reads.join(', ')}
    """

    file("${meta.id}_report.txt").text = report
}
```

Şimdi, ana iş akışımızı güncelleyelim:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENERATE_REPORT } from './modules/report'

workflow {
    // CSV dosyasını oku ve parse et
    Channel
        .fromPath('data/samples.csv')
        .splitCsv(header: true)
        .map { row ->
            // Meta nesnesini ve girdi dosya yollarını yarat
            meta = [id: row.id, condition: row.condition]
            reads = [row.read1, row.read2]
            tuple(meta, reads)
        }
        .set { ch_samples }

    // Rapor oluşturmak için process'i kullan
    GENERATE_REPORT(ch_samples)
}
```

Çalıştırın ve process kullanarak oluşturulan raporları görün:

```bash
nextflow run main.nf
```

### 1.3. Veri Doğrulama ve Hata İşleme Ekleme

Gerçek dünyada, veriler genellikle hatalıdır. Eksik değerleri ve geçersiz girişleri işlemeniz gerekir. Şimdi, veri doğrulama ekleyelim:

```groovy title="main.nf" linenums="1" hl_lines="14-22"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENERATE_REPORT } from './modules/report'

workflow {
    // CSV dosyasını oku ve parse et
    Channel
        .fromPath('data/samples.csv')
        .splitCsv(header: true)
        .map { row ->
            // Gerekli alanların mevcut olup olmadığını kontrol et
            if (!row.id || !row.read1 || !row.read2) {
                log.warn "Eksik alanlar ile satır atlandı: $row"
                return null  // Filter operatörü için null dön
            }

            // Meta nesnesini ve girdi dosya yollarını yarat
            meta = [id: row.id, condition: row.condition ?: 'unknown']
            reads = [row.read1, row.read2]
            tuple(meta, reads)
        }
        .filter { it != null }  // Eksik verilere sahip satırları filtrele
        .set { ch_samples }

    // Rapor oluşturmak için process'i kullan
    GENERATE_REPORT(ch_samples)
}
```

Bu kod artık eksik gerekli alanların kontrolünü sağlar ve `condition` belirtilmediği takdirde bir varsayılan değer atar. `?:` operatörü (Elvis operatörü olarak bilinir) burada kullanılır; bunu daha sonra bu öğreticide daha ayrıntılı bir şekilde inceleyeceğiz.

### 1.4. Gerçek Dünya Dosyalarını Ele Almak

Şimdi, veri dosyalarınızın gerçekten var olduğundan emin olmak için bir adım daha ileri gidelim. `File` nesneleri üzerinde `exists()` methodunu kullanacağız ve eğer dosyalar bulunamazsa bir hata günlüğü oluşturacağız:

```groovy title="main.nf" linenums="1" hl_lines="14-22 25-34"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENERATE_REPORT } from './modules/report'

workflow {
    // CSV dosyasını oku ve parse et
    Channel
        .fromPath('data/samples.csv')
        .splitCsv(header: true)
        .map { row ->
            // Gerekli alanların mevcut olup olmadığını kontrol et
            if (!row.id || !row.read1 || !row.read2) {
                log.warn "Eksik alanlar ile satır atlandı: $row"
                return null  // Filter operatörü için null dön
            }

            // Meta nesnesini ve girdi dosya yollarını yarat
            meta = [id: row.id, condition: row.condition ?: 'unknown']
            reads = [row.read1, row.read2]
            tuple(meta, reads)
        }
        .filter { it != null }  // Eksik verilere sahip satırları filtrele
        .map { meta, reads ->
            // Okuma dosyalarının mevcut olduğunu doğrula
            def read1_file = file(reads[0])
            def read2_file = file(reads[1])

            if (!read1_file.exists() || !read2_file.exists()) {
                log.warn "Bir veya daha fazla okuma dosyası bulunamadı: ${reads.join(', ')}"
                return null
            }

            tuple(meta, [read1_file, read2_file])
        }
        .filter { it != null }  // Kayıp dosyalara sahip satırları filtrele
        .set { ch_samples }

    // Rapor oluşturmak için process'i kullan
    GENERATE_REPORT(ch_samples)
}
```

!!! note "Channel Operatörleri Nasıl Çalışır"

    Yukarıdaki örnekte, `.map{}` ve `.filter{}` operatörlerini kullanarak bir kanal üzerinde dönüşümler gerçekleştirdik.
    Bu, herhangi bir Groovy koleksiyonuyla yapabileceğiniz şeylere benzer, ancak **paralel olarak** ve **tembel (lazy) değerlendirme** ile gerçekleştirilir.

    `.map{}` operatörü, her öğeyi dönüştürür (kanaldaki her mesaj için closure çağırılır). Her öğe için, `.filter{}` bir boolean değeri geliştirdiğinden dahil edilip edilmeyeceğini belirler.

### 1.5. Dosya Adı Şablonlarıyla Çalışma

Biyoinformatik iş akışlarında, genellikle veri dosyalarının belirli bir adlandırma kalıbını takip etmesini beklersiniz. Şimdi, dosya adlarından metadata çıkaran bir örnek ekleyelim:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENERATE_REPORT } from './modules/report'

// Dosya adından metadata çıkaran fonksiyon
def extractMetaFromFilename(filename) {
    // Örnek formatı: sample_condition_rep.fq.gz
    def matcher = filename =~ /(.+)_(.+)_rep(\d+).*\.fq\.gz$/

    if (matcher.matches()) {
        def sample = matcher[0][1]
        def condition = matcher[0][2]
        def replicate = matcher[0][3]
        return [sample: sample, condition: condition, replicate: replicate]
    }

    return null
}

workflow {
    // Veri klasöründen dosya yolları için bir kanal oluştur
    Channel
        .fromFilePairs('data/*_{1,2}.fq.gz', checkIfExists: true)
        .map { id, files ->
            def meta = extractMetaFromFilename(id)

            if (!meta) {
                log.warn "Dosya adından meta çıkarılamadı: $id"
                return null
            }

            tuple(meta, files)
        }
        .filter { it != null }
        .set { ch_samples }

    // Rapor oluşturmak için process'i kullan
    GENERATE_REPORT(ch_samples)
}
```

Bu örnek, düzenli ifade eşleşmesini kullanarak dosya adından bilgi çıkarır. Biyoinformatik iş akışlarında sıklıkla karşılaşacağınız yaygın bir pattern'dir.

### 1.6. Zorluklarla Yüzleşmek

Buraya kadar, veri okuma, process oluşturma, hata işleme ve daha fazlasını içeren bir iş akışı geliştirdik. Ancak, daha güçlü ve karmaşık iş akışları için daha fazla pattern'e ihtiyacımız var:

1. **Veri Yapılarını İşlemek**: Bir `meta` map'ini, liste içindeki listelerle nasıl işlersiniz?
2. **Karmaşık Dosya Adı Desenleri**: Karmaşık regex pattern'lerini nasıl kullanırsınız?
3. **Yeniden Kullanılabilir Fonksiyonlar**: Kod duplikasyonunu nasıl önlersiniz?
4. **Koşullu İş Akışı Mantığı**: Koşullara dayalı olarak farklı process'ler nasıl çalıştırırsınız?
5. **Dinamik Kaynak Tahsisi**: Her process için nasıl dinamik kaynaklar ayırırsınız?

Bu zorlukların her birini, bu dersin geri kalanında ele alacağız.

## 2. Karmaşık Veri Yapıları ve Dil Özellikleri

Gerçek dünya iş akışlarında, karmaşık veri yapılarını manipüle etmek ve dilsel yeteneklerden yararlanmak için ek becerilere ihtiyacımız var.

### 2.1. Map ve Koleksiyonları İşleme

Bir genomik varyantları işleme örneği ile map ve koleksiyon işlemeyi gösterelim:

```groovy title="modules/variant_caller.nf" linenums="1"
process CALL_VARIANTS {
    input:
    tuple val(meta), path(bam), path(bai)
    path(ref_fasta)

    output:
    tuple val(meta), path("${meta.id}.vcf")

    script:
    """
    # Chromosome ayarlarına dayalı komutları üret
    gatk HaplotypeCaller \
        -R $ref_fasta \
        -I $bam \
        --output ${meta.id}.vcf \
        ${meta.intervals ? "--intervals ${meta.intervals}" : ""}
    """
}
```

Şimdi, ana iş akışında bazı karmaşık map manipülasyonları kullanacağız:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CALL_VARIANTS } from './modules/variant_caller'

workflow {
    // Örnek metadata'sını içeren bir kanal oluştur
    Channel
        .fromPath('data/samples.csv')
        .splitCsv(header: true)
        .map { row ->
            // Varsayılan metadata ile başla
            def meta = [
                id: row.id,
                intervals: null,
                is_tumor: false,
                options: [:]
            ]

            // Koşul bilgisi varsa ekle
            if (row.condition) {
                meta.condition = row.condition

                // Tumor numuneleri için ek ayarlar ekle
                if (row.condition == 'tumor') {
                    meta.is_tumor = true
                    meta.options.put('contamination', 0.02)
                    meta.options.put('tumor_lod', 3.5)
                }
            }

            // İnterval bilgisi varsa ekle
            if (row.intervals) {
                meta.intervals = row.intervals
            }

            tuple(meta, file(row.bam), file("${row.bam}.bai"))
        }
        .set { ch_bam_files }

    // Referans genomu
    ch_ref = Channel.fromPath('data/reference.fasta')

    // Varyant çağrıcı process'i çalıştır
    CALL_VARIANTS(ch_bam_files, ch_ref)
}
```

!!! note "Map Manipülasyonu"

    Önceki örnekte, başlangıç değerleriyle bir map oluşturduk, ardından koşullu ifadelerle özellikleri ekledik veya değiştirdik. Bu esnek, karmaşık metadata işleme mekanizması sağlar.

    `[:]` sözdizimi boş bir map oluşturur ve `put()` methodu ona yeni değerler ekler. Sık kullanılan alternatif, sadece `map.key = value` sözdizimini kullanmaktır.

### 2.2. Process Script Bloklarında Scripting

Process'lerin gerçek dünya uyumluluğunu arttırmak için, karmaşık komutlar oluşturmak üzere script bloklarında scripting'i kullanabilirsiniz:

```groovy title="modules/genomics_db.nf" linenums="1"
process CREATE_GENOMICS_DB {
    input:
    tuple val(cohort_name), val(sample_list), path(gvcfs), path(indices)

    output:
    tuple val(cohort_name), path("${cohort_name}_gdb")

    script:
    // Örnek adları ve gVCF dosyaları içeren bir map oluştur
    def sample_map_items = []
    sample_list.eachWithIndex { sample_name, idx ->
        sample_map_items << "${sample_name}=${gvcfs[idx]}"
    }

    // -V argümanlarını yapılandır
    def sample_map = sample_map_items.join(" -V ")

    """
    gatk GenomicsDBImport \
        -V ${sample_map} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
}
```

Process script bloklarında scripting kullanmanın bu kalıpları son derece güçlüdür ve birçok senaryoda uygulanabilir - değişken girdi türlerini işlemekten dosya koleksiyonlarından karmaşık komut satırı argümanları oluşturmaya kadar, process'lerinizi gerçek dünyadaki verilerin çeşitli gereksinimlerine gerçekten uyumlu hale getirir.

### 2.3. Değişken Enterpolasyonu: Nextflow ve Shell Değişkenleri

Process script'lerinde Nextflow değişkenleri, shell değişkenleri ve komut ikamelerini karıştırırız, her biri farklı enterpolasyon sözdizimi ile. Yanlış sözdizimi kullanmak hatalara neden olur. Bunları, bir işlem raporu oluşturan bir process ile keşfedelim.

`modules/generate_report.nf` modül dosyasına bir göz atalım:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    # İşlem raporu oluşturuyoruz

    # İş akışı kimliğini elde et (Bash değişkeni)
    workflow_id=\$(date +%Y%m%d_%H%M%S)

    # İşlem ID'sini Bash'te tanımla
    process_id="PROCESS_\$RANDOM"

    # Bir rapor oluştur - hem Nextflow hem de Bash değişkenlerini kullanarak
    cat <<EOF > ${meta.id}_report.txt
    ===============================
    İŞLEM RAPORU
    ===============================
    İş akışı: \$workflow_id
    Process: \$process_id

    Numune ID: ${meta.id}
    Durum: ${meta.condition ?: 'Bilinmiyor'}

    Dosyalar işlendi:
    - ${reads[0]}
    - ${reads[1]}

    Dil İstatistikleri:
    Nextflow: ${nextflow.version}
    Java: \$(java -version 2>&1 | grep version | awk '{print \$3}' | tr -d '"')
    ===============================
    EOF
    """
}
```

Bu kod değişken enterpolasyonu için üç farklı sözdizimi kullanır:

1. **Nextflow değişkenleri** (`${meta.id}`, `${nextflow.version}`, vb.) - bunlar Nextflow tarafından değerlendirilir ve script bloğuna gönderilmeden önce değerleriyle değiştirilir.

2. **Kaçan shell değişkenleri** (`\$workflow_id`, `\$process_id`, vb.) - bunlar shell tarafından değerlendirilir, `\` ön eki Nextflow'un bunları değerlendirmesini önler.

3. **Kaçan komut ikamesi** (`\$(komut)`) - bunlar shell tarafından yürütülür, `\` ön eki Nextflow'un bunları çalıştırmasını önler.

Doğru enterpolasyon seçiminde başarısızlık, hatalarla sonuçlanabilir.

### 2.4. Güvenli Gezinme ve Elvis Operatörü

Nextflow, nesneleri ve map'leri gezinmek için Groovy'den iki güçlü operatör devralmıştır: Güvenli Gezinme Operatörü (`?.`) ve Elvis Operatörü (`?:`). Bunları kompleks genomik metadata'yı ele almak için kullanalım:

```groovy title="modules/genome_metrics.nf" linenums="1"
process GENOME_METRICS {
    input:
    tuple val(meta), path(vcf)

    output:
    path "${meta.id}_metrics.txt"

    script:
    // Güvenli gezinme - alt özelliklere bakmadan önce null kontrolü yapar
    def ploidy = meta.genome?.ploidy ?: 2
    def build = meta.genome?.build ?: 'hg38'
    def species = meta.genome?.species ?: 'human'

    // Metrikleri hesaplarken güvenli gezinme ve varsayılan değerler kullan
    """
    # Varyant metriklerini hesapla
    bcftools stats ${vcf} > stats.txt

    # Özel metrik sonuçlarını hazırla
    echo "Numune: ${meta.id}" > ${meta.id}_metrics.txt
    echo "Referans Yapı: ${build}" >> ${meta.id}_metrics.txt
    echo "Tür: ${species}" >> ${meta.id}_metrics.txt
    echo "Ploidi: ${ploidy}" >> ${meta.id}_metrics.txt

    # Filtreleme için minimum kaliteyi belirle
    min_qual=${meta.filters?.min_qual ?: 30}
    echo "Kaliteli varyantlar: \$(grep -v ^# ${vcf} | awk '\$6>='\${min_qual} | wc -l)" >> ${meta.id}_metrics.txt
    """
}
```

Bu örnekte:

- `?.` operatörü, `meta.genome` null ise `meta.genome?.ploidy`'nin de null döndürmesini sağlar, NullPointerException'ları önler
- `?:` operatörü, sol taraftaki değer null veya false ise sağ taraftaki değeri döndürür, pratik varsayılan değerlere izin verir

Ana iş akışında nasıl kullanıldığına bakalım:

```groovy title="main.nf" linenums="1"
workflow {
    Channel
        .fromPath('data/samples.csv')
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.id,
                condition: row.condition
            ]

            // Genom bilgisi tanımlıysa ekle
            if (row.build) {
                meta.genome = [
                    build: row.build,
                    species: row.species ?: (row.build.startsWith('hg') ? 'human' : 'unknown')
                ]

                // Ploidi belirtilmişse ekle
                if (row.ploidy) {
                    meta.genome.ploidy = row.ploidy.toInteger()
                }
            }

            // Filtreler belirtilmişse ekle
            if (row.min_qual) {
                meta.filters = [min_qual: row.min_qual.toInteger()]
            }

            tuple(meta, file(row.vcf))
        }
        .set { ch_vcf_files }

    GENOME_METRICS(ch_vcf_files)
}
```

!!! tip "Güvenli Gezinme ve Elvis Operatörlerinin Gücü"

    Bu operatörler birleştiğinde, eklenmemiş olabilecek özellikler için uzun koşullu ifadeleri önemli ölçüde kısaltır. Değişken ve koşullu olarak mevcut olan özelliklere sahip karmaşık veri yapılarına sahip olduğunuzda özellikle değerlidir.

### 2.5. Regex Kalıpları

Düzenli ifadeler, karmaşık dosya adlandırma şemalarını analiz etmede güçlü araçlardır. Bir hücre hattı numunelerini içeren bir örnek üzerinden nasıl kullanıldığını görelim:

```groovy title="modules/cell_line.nf" linenums="1"
process CELL_LINE_REPORT {
    input:
    path(fastq_files)

    output:
    path "cell_lines_summary.txt"

    exec:
    // Hücre hattı metriklerini içerecek boş bir map oluştur
    def cell_lines = [:]

    // Her fastq dosyası için
    fastq_files.each { fastq ->
        // Dosya adından hücre hattı bilgisini çıkar
        // Beklenen format: CELL-LINE_CONDITION_REP#_R#.fastq.gz
        // Örnek: MCF7_TREATED_REP1_R1.fastq.gz
        def matcher = fastq.getName() =~ /^([^_]+)_([^_]+)_REP(\d+)_R(\d)\.fastq\.gz$/

        if (matcher.matches()) {
            def cell_line = matcher.group(1)  // İlk grup hücre hattı adı
            def condition = matcher.group(2)  // İkinci grup koşul
            def replicate = matcher.group(3)  // Üçüncü grup replikattır
            def read_num = matcher.group(4)   // Dördüncü grup okuma numarasıdır

            // Cell_lines map'inde hücre hattı için bir giriş oluştur veya güncelle
            if (!cell_lines.containsKey(cell_line)) {
                cell_lines[cell_line] = [conditions: [:], total_files: 0]
            }

            // Koşul için bir giriş oluştur veya güncelle
            if (!cell_lines[cell_line].conditions.containsKey(condition)) {
                cell_lines[cell_line].conditions[condition] = [:]
            }

            // Replikat için bir giriş oluştur veya güncelle
            if (!cell_lines[cell_line].conditions[condition].containsKey(replicate)) {
                cell_lines[cell_line].conditions[condition][replicate] = []
            }

            // Dosyayı ekle ve sayaçları güncelle
            cell_lines[cell_line].conditions[condition][replicate] << fastq
            cell_lines[cell_line].total_files++
        }
    }

    // Özeti dosyaya yaz
    def summary = file("cell_lines_summary.txt")
    summary.text = ""

    cell_lines.each { cell_line, data ->
        summary.text += "Hücre Hattı: ${cell_line}\n"
        summary.text += "Toplam Dosyalar: ${data.total_files}\n"

        data.conditions.each { condition, replicates ->
            summary.text += "  Koşul: ${condition}\n"

            replicates.each { replicate, files ->
                summary.text += "    Replikat ${replicate}: ${files.size()} dosya(lar)\n"
            }
        }
        summary.text += "\n"
    }
}
```

Bu process, dosya adından meta bilgileri çıkarmak ve organize etmek için düzenli ifade eşleşmesini kullanır. Groovy'nin pattern eşleşme sözdizimi `=~` ile düzenli ifade oluşturur, `.matches()` tüm desenin eşleşip eşleşmediğini kontrol eder ve `.group(#)` methodu eşleşen grupları çıkarır.

Şimdi, bu process'i çalıştıracak bir iş akışı oluşturalım:

```groovy title="main.nf" linenums="1"
workflow {
    // FASTQ dosyaları için kanal oluştur
    Channel
        .fromPath('data/*.fastq.gz')
        .collect()  // Tüm dosyaları tek bir liste olarak topla
        .set { ch_fastq_files }

    CELL_LINE_REPORT(ch_fastq_files)
}
```

!!! note "Regex Eşleştirmesinde Parantezler"

    Düzenli ifadelerde, parantezler yakalama gruplarını tanımlar. Örneğin, `([^_]+)` bir yakalama grubudur ve underscore olmayan karakterlerin bir veya daha fazla tekrarıyla eşleşir. Bu grupların değerlerini daha sonra `matcher.group(1)`, `matcher.group(2)`, vb. ile çıkarabiliriz.

    `matcher.group(0)` her zaman **tüm** eşleşmeyi içerir, `matcher.group(1)` ilk yakalama grubunu içerir, vb.

### 2.6. Çok Satırlı Dizeleri İşleme

Nextflow'da, heredoc sözdizimini kullanarak çok satırlı komutlar veya yapılandırma dosyaları oluşturabilirsiniz:

```groovy title="modules/config_generator.nf" linenums="1"
process GENERATE_CONFIG {
    input:
    tuple val(meta), path(input_files)

    output:
    tuple val(meta), path("${meta.id}_config.yaml")

    exec:
    // Dinamik YAML yapılandırma dosyası oluştur
    def tool_options = ""

    // Metadata'daki ek seçenekleri ekle
    meta.options?.each { key, value ->
        if (value instanceof Boolean) {
            tool_options += "${key}: ${value ? 'true' : 'false'}\n"
        } else if (value instanceof List) {
            tool_options += "${key}: [${value.join(', ')}]\n"
        } else {
            tool_options += "${key}: ${value}\n"
        }
    }

    // Girdileri için bir bölüm ekle
    def inputs_section = ""
    input_files.each { file ->
        inputs_section += "  - ${file}\n"
    }

    // Dinamik yapılandırmayı doldur
    def config_content = """
    # ${meta.id} için otomatik oluşturulmuş yapılandırma
    # Tarih: ${new Date().format('yyyy-MM-dd')}

    sample_id: ${meta.id}
    analysis_type: ${meta.type ?: 'default'}

    # Varsayılan araç seçenekleri
    options:
      threads: 4
      memory: 8G
      ${tool_options}

    # Girdi dosyaları
    inputs:
    ${inputs_section}
    """

    // Yapılandırma dosyasını yaz
    file("${meta.id}_config.yaml").text = config_content
}
```

Çok satırlı dizeler, değişken enterpolasyonu, koşullu içerikleri ve metadata'ya göre dinamik yapılandırmayı mümkün kılar. Bu pattern, özellikle çeşitli özelliklere dayalı olarak karmaşık araç yapılandırmaları oluşturmak için kullanışlıdır.

## 3. Yeniden Kullanılabilir Fonksiyonlar ve Closure'lar

Daha temiz, daha bakımı kolay kod için, sık kullanılan işlemleri yeniden kullanılabilir fonksiyonlara veya closure'lara çıkarmak en iyi uygulamadır.

### 3.1. Dosya Adı Ayrıştırma için Yeniden Kullanılabilir Fonksiyonlar

Önceden karmaşık bir dosya adı ayrıştırıcı oluşturmak, iş akışınızı basitleştirir:

```groovy title="lib/functions.nf" linenums="1"
// Ortak işlevleri içeren bir kütüphane tanımlayın

// FASTQ dosya adlarını ayrıştır ve bir metadata map'ı döndür
def parseFastqName(filename) {
    // Önce uzantı olmadan temel adı al
    def basename = filename.getName().replaceAll(/\.fastq\.gz$|\.fq\.gz$/, "")

    // Farklı adlandırma formatları için çeşitli kalıpları deneyebilirsiniz

    // Format 1: SAMPLENAME_CONDITION_REP#_R#
    def pattern1 = /^(.+)_(.+)_REP(\d+)_R(\d)$/
    def matcher1 = basename =~ pattern1

    if (matcher1.matches()) {
        return [
            id: matcher1[0][1],
            condition: matcher1[0][2],
            replicate: matcher1[0][3].toInteger(),
            read_number: matcher1[0][4].toInteger(),
            format: 'standard'
        ]
    }

    // Format 2: SAMPLENAME.CONDITION.R#
    def pattern2 = /^(.+)\.(.+)\.R(\d)$/
    def matcher2 = basename =~ pattern2

    if (matcher2.matches()) {
        return [
            id: matcher2[0][1],
            condition: matcher2[0][2],
            read_number: matcher2[0][3].toInteger(),
            format: 'alternate'
        ]
    }

    // Diğer formatları burada ekleyin...

    // Eğer hiçbiri eşleşmezse, minimum ayrıştırma yap
    log.warn "Dosya adı için tanınan şablon yok: ${filename}"

    // R1/R2 için son deneyin
    def readPattern = /_R(\d)$/
    def readMatcher = basename =~ readPattern

    def readNum = readMatcher.find() ? readMatcher[0][1].toInteger() : null

    return [
        id: basename.replaceAll(/_R\d$/, ""),
        read_number: readNum,
        format: 'unknown'
    ]
}

// Referans dosyalarını doğrulamak için bir yardımcı fonksiyon
def validateReference(ref_file, genome_build) {
    if (!ref_file.exists()) {
        exit 1, "Hata: Referans dosyası bulunamadı: ${ref_file}"
    }

    // Genome build'e göre doğrulama kontrollerini buraya ekleyin
    if (genome_build.startsWith('hg')) {
        // İnsan genomu kontrolleri
        // ...
    } else if (genome_build.startsWith('mm')) {
        // Fare genomu kontrolleri
        // ...
    }

    return true
}

// İki dosya listesini metadata'ya göre birleştir
def mergeLists(list1, list2, comparator) {
    def result = []

    list1.each { item1 ->
        def match = list2.find { item2 -> comparator(item1, item2) }
        if (match) {
            result << [item1, match]
        }
    }

    return result
}
```

İş akışında bu fonksiyonları nasıl kullanacağımızı görelim:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Özel işlevleri içe aktar
include { parseFastqName; validateReference } from './lib/functions'

workflow {
    // FASTQ dosyaları için kanal oluştur
    Channel
        .fromPath('data/*.fastq.gz')
        .map { file ->
            def meta = parseFastqName(file)
            tuple(meta, file)
        }
        .branch {
            read1: it[0].read_number == 1
            read2: it[0].read_number == 2
        }
        .set { ch_reads }

    // R1 ve R2'yi numune ID'lerine göre eşleştir
    ch_reads.read1
        .map { meta, file -> tuple(meta.id, meta, file) }
        .set { ch_read1 }

    ch_reads.read2
        .map { meta, file -> tuple(meta.id, meta, file) }
        .set { ch_read2 }

    ch_read1
        .join(ch_read2, by: 0) // ID'ye göre eşleştir
        .map { id, meta1, file1, meta2, file2 ->
            // Read1 meta'sını temel al, ancak her ikisindeki faydalı bilgileri birleştir
            def meta = meta1 + [
                read_files: [file1, file2],
                read_count: 2
            ]

            tuple(meta, [file1, file2])
        }
        .set { ch_read_pairs }

    // Şimdi doğrulanmış referans dosyalarıyla eşleştir
    ch_ref = Channel.fromPath('data/reference.fasta')

    // Verilen genomik yapı için referans geçerliliğini doğrula
    validateReference(file('data/reference.fasta'), 'hg38')
}
```

Bu örnekte, `parseFastqName` ve `validateReference` işlevlerini kullanarak iş akışımızın daha temiz ve modüler olmasını sağladık.

### 3.2. Dinamik Kaynak Ayarlama için Closure'lar

Yeniden kullanılabilir closure'lar, kaynak gereksinimlerinin proses bazında dinamik olarak yapılandırılması için oldukça kullanışlıdır:

```groovy title="conf/resources.config" linenums="1"
// Kaynak yapılandırması için closure'lar tanımla
def standard_resources(task) {
    // Hafıza ve CPU'yu görev adına ve girdi boyutuna göre tahsis et
    def base_memory = 8.GB
    def base_cpus = 4

    // Çok fazla veri ile çalıştığını bildiğimiz belirli görevler için hafıza artır
    if (task.name == 'GENOME_METRICS' || task.name == 'CALL_VARIANTS') {
        base_memory = 16.GB
    }

    // Çoklu-thread özellikleri olan görevler için CPU'ları artır
    if (task.name == 'CALL_VARIANTS' || task.name == 'MAPPING') {
        base_cpus = 8
    }

    // Girdi dosya boyutuna göre belleği ölçeklendir (varsa)
    if (task.attempt > 1) {
        base_memory = base_memory * task.attempt
    }

    return [
        memory: base_memory,
        cpus: base_cpus,
        time: 8.hour
    ]
}

// Küçük görevler için kaynakları optimize et
def light_resources(task) {
    def memory = 2.GB

    if (task.attempt > 1) {
        memory = memory * task.attempt
    }

    return [
        memory: memory,
        cpus: 2,
        time: 1.hour
    ]
}

// HPC-specific görevler için özel kaynaklar
def hpc_resources(task) {
    // HPC-specific ayarları dahil et
    def resources = standard_resources(task)
    resources.queue = 'high-memory'
    resources.clusterOptions = '--account=genomics'
    return resources
}

// SLURM'a özel ayarlar
def slurm_resources(task) {
    def resources = standard_resources(task)
    resources.clusterOptions = "--account=genomics --partition=${task.name == 'CALL_VARIANTS' ? 'highmem' : 'normal'}"
    return resources
}
```

Bunu `nextflow.config` dosyasında nasıl kullanacağınız:

```groovy title="nextflow.config" linenums="1"
// Kaynak closure'larını içe aktar
includeConfig 'conf/resources.config'

process {
    // Varsayılan olarak tüm process'lere standart kaynaklar uygula
    executor = 'local'

    withName: 'GENERATE_REPORT|GENOME_METRICS' {
        // Basit görevler için daha düşük kaynaklar
        memory = { light_resources(task).memory }
        cpus = { light_resources(task).cpus }
        time = { light_resources(task).time }
    }

    withName: 'CALL_VARIANTS|CREATE_GENOMICS_DB' {
        // Standart kaynakları kullan, görev adı ve girişimler dahil
        memory = { standard_resources(task).memory }
        cpus = { standard_resources(task).cpus }
        time = { standard_resources(task).time }
    }
}

// Ortama özgü profiller
profiles {
    hpc {
        process {
            executor = 'sge'

            // HPC-specific kaynakları kullan
            memory = { hpc_resources(task).memory }
            cpus = { hpc_resources(task).cpus }
            time = { hpc_resources(task).time }
            queue = { hpc_resources(task).queue }
            clusterOptions = { hpc_resources(task).clusterOptions }
        }
    }

    slurm {
        process {
            executor = 'slurm'

            // SLURM-specific kaynakları kullan
            memory = { slurm_resources(task).memory }
            cpus = { slurm_resources(task).cpus }
            time = { slurm_resources(task).time }
            clusterOptions = { slurm_resources(task).clusterOptions }
        }
    }
}
```

Bu örnekteki closure'lar, görev adı ve girişim sayısı gibi bağlamsal bilgilere dayalı olarak kaynakların dinamik olarak yapılandırılmasını sağlar. Bu yaklaşım, veri boyutuna ve görev karmaşıklığına göre ölçeklenebilen, yüksek oranda yapılandırılabilir iş akışları sağlar.

### 3.3. Fonksiyonlarla Koşullu Mantık

Karmaşık, koşullu mantığı fonksiyonlara çıkarmak, iş akışı okunabilirliğini geliştirebilir:

```groovy title="modules/alignment_workflow.nf" linenums="1"
include { MAPPING as MAPPING_WITH_ALT } from './mapping'
include { MAPPING as MAPPING_STANDARD } from './mapping'
include { QC_CHECK } from './qc'

// Mantık için yardımcı işlev
def shouldUseAltPipeline(meta) {
    // İnsan genomu ve tumor numunesi için alt pipeline'ı kullan
    if (meta.species == 'human' && meta.condition == 'tumor') {
        return true
    }

    // Alternatif referans bulunan metadatada özellikle belirtildiyse
    if (meta.options?.use_alt_pipeline) {
        return true
    }

    // 36bp'den kısa okuma uzunluğu için özel hizalama kullan
    if (meta.read_length && meta.read_length < 36) {
        return true
    }

    // Varsayılan olarak standart pipeline
    return false
}

workflow ALIGNMENT_WORKFLOW {
    take:
    reads_ch     // meta ve okuma dosyalarıyla tuple kanalı

    main:
    // Koşullu mantığı uygulamak için kanalı ikiye ayır
    reads_ch.branch {
        alt_pipeline: { meta, reads -> shouldUseAltPipeline(meta) }
        standard_pipeline: { meta, reads -> !shouldUseAltPipeline(meta) }
    }.set { split_reads }

    // Her pipeline için ayrı modül çalıştır
    MAPPING_WITH_ALT(split_reads.alt_pipeline)
    MAPPING_STANDARD(split_reads.standard_pipeline)

    // Her numuneyi QC kontrolünden geçir
    QC_CHECK(MAPPING_WITH_ALT.out.mix(MAPPING_STANDARD.out))

    emit:
    bams = QC_CHECK.out // QC kontrolünden geçen BAM dosyalarını çıkart
}
```

Ana iş akışında bu modülü nasıl kullanacağımızı görelim:

```groovy title="main.nf" linenums="1" hl_lines="18-20 31"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Modulü içe aktar
include { ALIGNMENT_WORKFLOW } from './modules/alignment_workflow'
include { parseFastqName } from './lib/functions'

workflow {
    // FASTQ dosyaları için kanal oluştur
    Channel
        .fromFilePairs('data/*_{1,2}.fastq.gz')
        .map { id, files ->
            // Dosya adından bilgileri çıkar
            def meta = parseFastqName(files[0])

            // Gerekirse ek metadata ekle
            if (id.contains('tumor')) {
                meta.condition = 'tumor'
                meta.options = [use_alt_pipeline: true]
            }

            // Okuma uzunluğunu belirle (gerçekte, bunun için BBTools gibi bir araç kullanmak isteyebilirsiniz)
            meta.read_length = 75 // örnek değer

            tuple(meta, files)
        }
        .set { ch_reads }

    // Hizalama iş akışını çalıştır
    ALIGNMENT_WORKFLOW(ch_reads)
}
```

Bu yaklaşım, karmaşık koşullu mantığı tek bir yardımcı fonksiyona (`shouldUseAltPipeline`) çıkarak modülariteyi artırır ve okunabilirliği geliştirir. Bu yapılandırmada, koşullu hizalama mantığını değiştirmek istediğimizde, sadece bir yerde değişiklik yapmamız yeterlidir.

## 4. Gerçek Dünya Karmaşıklıkları Ele Alma

Son bölümümüzde, gerçek dünya pipeline'larında karşılaşabileceğiniz karmaşıklıklardan bazılarını ele alacağız.

### 4.1. Ortama göre Adapte Olan Workflow'lar

Üretime hazır workflow'lar, farklı çalışma ortamlarına ve kullanıcı gereksinimlerine uyum sağlayabilmelidir:

```groovy title="modules/environment_detector.nf" linenums="1"
// Mevcut ortamı algıla ve ilgili yapılandırmaları ayarla
def detectEnvironment() {
    def env = [:]

    // İşletim sistemi ve kaynak kullanılabilirliğini algıla
    env.os = System.properties['os.name'].toLowerCase()
    env.is_linux = env.os.contains('linux')
    env.is_mac = env.os.contains('mac')
    env.is_windows = env.os.contains('windows')

    // Kullanılabilir bellek ve CPU'yu algıla
    env.max_memory = Runtime.getRuntime().maxMemory()
    env.available_cpus = Runtime.getRuntime().availableProcessors()

    // Docker kullanılabilirliğini kontrol et
    def docker_process = "docker info".execute()
    docker_process.waitFor()
    env.has_docker = docker_process.exitValue() == 0

    // Singularity kullanılabilirliğini kontrol et
    def singularity_process = "singularity --version".execute()
    singularity_process.waitFor()
    env.has_singularity = singularity_process.exitValue() == 0

    // GPU kullanılabilirliğini kontrol et (Linux'ta)
    if (env.is_linux) {
        def nvidia_process = "nvidia-smi".execute()
        nvidia_process.waitFor()
        env.has_gpu = nvidia_process.exitValue() == 0
    } else {
        env.has_gpu = false
    }

    return env
}

// Bir ortam algılayıcı çalıştır
def environment = detectEnvironment()

// Ortama dayalı olarak bir konteyner seçin
def selectContainer() {
    if (params.container) {
        return params.container
    }

    // Docker kullanılabilirse, varsayılan docker imajını kullan
    if (environment.has_docker) {
        return "docker://quay.io/biocontainers/biopython:1.78"
    }

    // Singularity kullanılabilirse, varsayılan singularity imajını kullan
    if (environment.has_singularity) {
        return "https://depot.galaxyproject.org/singularity/biopython:1.78--py39h38f01e4_0"
    }

    // Konteyner yok, uyarı
    log.warn "Konteynerler kullanılamıyor - yerel kuruluma güveniliyor"
    return null
}

// GPU gerektiren görevler için kullanılabilirlik kontrol et
def hasGpuSupport() {
    return environment.has_gpu && (environment.has_docker || environment.has_singularity)
}

// Ortama dayalı olarak maksimum kaynakları hesapla
def getMaxResources() {
    // Toplam CPU'nun %75'i ve toplam belleğin %80'ini kullan
    def max_cpus = Math.max(1, (environment.available_cpus * 0.75).intValue())
    def max_memory = (environment.max_memory * 0.8).intValue()

    return [cpus: max_cpus, memory: max_memory]
}
```

Şimdi, bu ortam algılamayı ana nextflow.config dosyasında nasıl kullanabiliriz:

```groovy title="nextflow.config" linenums="1"
// Ortam algılama fonksiyonlarını içe aktar
includeConfig 'modules/environment_detector.nf'

// Ortama göre konteyner seç
process.container = selectContainer()

// Kaynak sınırlarını hesapla
def max_resources = getMaxResources()

// GPU desteği varsa GPU kullanan profili etkinleştir
if (hasGpuSupport()) {
    includeConfig 'conf/gpu.config'
}

// Maksimum kaynak tahsisi
executor {
    $local {
        cpus = max_resources.cpus
        memory = max_resources.memory
    }
}

// İşletim sistemine dayalı profile'ları etkinleştir
profiles {
    standard {
        process.executor = 'local'
    }

    mac {
        process.executor = 'local'
        docker.runOptions = '-e TMPDIR=/tmp'
    }

    linux {
        process.executor = 'local'
    }

    // Ortama göre uygun profile'ı otomatik seç
    if (environment.is_mac) {
        includeConfig 'conf/mac.config'
    } else if (environment.is_linux) {
        includeConfig 'conf/linux.config'
    }
}
```

Bu örnekte, çalışma zamanında ortama göre adapte olabilen akıllı bir iş akışı yapılandırması oluşturduk. Kullanılabilir kaynakları, işletim sistemini ve konteyner kullanılabilirliğini otomatik olarak algılar ve en uygun ayarları seçer.

### 4.2. İç İçe Fonksiyonlar

Birbiriyle ilişkili olan fonksiyonları gruplamak, karmaşık scripti daha temiz hale getirebilir:

```groovy title="lib/utils.nf" linenums="1"
// Biyoinformatik araçlarıyla çalışmak için yardımcı fonksiyonlar
class BioTools {
    // FASTQ kalitesini kontrol et
    static Map checkFastqQuality(path) {
        // Karşılaşılacak kalite değerlerini belirleyin
        def detectQualityEncoding = { file ->
            // Phred+33 veya Phred+64 kodlamasını tespit et
            // Gerçekte, bu BBTools, FastQC vb. kullanır
            return 'Phred+33' // Örnek değer
        }

        // Kısa okuma sekansı sayısı
        def countShortReads = { file, min_length=30 ->
            // Minimum uzunluğun altındaki okumaları say
            // Gerçekte, awk/BBTools/seqtk gibi araçlar kullanın
            return 0 // Örnek değer
        }

        // N içeren okumaların sayısı
        def countNContent = { file ->
            // N içeren okumaları say
            return 0 // Örnek değer
        }

        return [
            encoding: detectQualityEncoding(path),
            short_reads: countShortReads(path),
            n_content: countNContent(path)
        ]
    }

    // BAM kalitesini kontrol et
    static Map checkBamQuality(path) {
        // Eşleme kalitesini kontrol et
        def checkMappingQuality = { file ->
            // Ortalama eşleme kalitesi (MAPQ) değeri
            return 30 // Örnek değer
        }

        // Kapsam derinliğini kontrol et
        def checkCoverage = { file ->
            // Ortalama kapsam derinliği
            return 30 // Örnek değer
        }

        // İç içe fonksiyonları kullanarak sonuçları döndür
        return [
            mapq: checkMappingQuality(path),
            coverage: checkCoverage(path)
        ]
    }
}

// Genomik metadata işleme fonksiyonları
class GenomicMeta {
    // Genom yapılandırması oluştur
    static Map createGenomeConfig(build) {
        // Genom yapısını varsayılan değerlerle başlat
        def config = [
            build: build,
            species: 'unknown'
        ]

        // Genom yapı tipini algıla
        def detectSpecies = { genome_build ->
            if (genome_build.startsWith('hg')) return 'human'
            if (genome_build.startsWith('mm')) return 'mouse'
            if (genome_build.startsWith('dm')) return 'fruitfly'
            if (genome_build.startsWith('TAIR')) return 'arabidopsis'
            return 'unknown'
        }

        // Kromozom uzunluklarını ata
        def getChromosomeSizes = { genome_build, species ->
            // Gerçekte, bu genome fasta'dan veya bir veritabanından alınır
            return [
                'chr1': 248956422,
                'chr2': 242193529
                // ... diğer kromozomlar
            ]
        }

        // Metadatayı tamamla
        config.species = detectSpecies(build)
        config.chromosomes = getChromosomeSizes(build, config.species)

        return config
    }
}
```

Bu fonksiyonları ana iş akışınızda nasıl kullanabileceğiniz:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Kütüphane fonksiyonlarını içe aktar
// import ilkel fonksiyon olmadığı için def anahtar sözcüğü kullanmıyoruz
import './lib/utils' as utils

// Utility sınıflarını adlandırma alanlarında açıkça kullan
def biotools = utils.BioTools
def genomic = utils.GenomicMeta

workflow {
    // Okuma dosyaları için kanal oluştur
    Channel
        .fromPath('data/*.fastq.gz')
        .map { file ->
            def qc_results = biotools.checkFastqQuality(file)

            // Çok düşük kaliteli dosyaları filtreleme
            if (qc_results.short_reads > 1000 || qc_results.n_content > 500) {
                log.warn "Düşük kaliteli dosya atlandı: $file"
                return null
            }

            tuple(file, qc_results)
        }
        .filter { it != null }
        .set { ch_filtered_reads }

    // Genom yapılandırması oluştur
    def genome_config = genomic.createGenomeConfig('hg38')
    log.info "Genom yapılandırması: $genome_config"
}
```

İç içe fonksiyon yaklaşımı, ilişkili işlemleri düzenli paketlerde gruplandırmanıza olanak tanır. Aynı zamanda, kodu temiz tutmaya ve işlevselliği bağlantılı birimlere ayırmaya yardımcı olur.

### 4.3. Etkin Dosya Ayrıştırma için Regex Desenleri

Genom ayrıştırma için daha gelişmiş bir regex deseni kullanalım:

```groovy title="modules/genome_parser.nf" linenums="1"
// Çeşitli genomik dosya adı formatlarını ayrıştırmak için işlev
def parseGenomicFilename(filename) {
    // FASTQ deseni
    // Format: SAMPLE_CONDITION_REPn_Rn.fastq.gz
    // Örnek: GM12878_CONTROL_REP3_R1.fastq.gz
    def fastq_pattern = ~/^([^_]+)_([^_]+)(?:_REP(\d+))?(?:_R([12]))?.*\.(?:fastq|fq)\.gz$/

    // BAM deseni
    // Format: SAMPLE.CONDITION.sorted.bam
    // Örnek: GM12878.CONTROL.sorted.bam
    def bam_pattern = ~/^([^\.]+)\.([^\.]+)(?:\.REP(\d+))?.*\.bam$/

    // VCF deseni
    // Format: SAMPLE-CONDITION-platform.vcf.gz
    // Örnek: GM12878-CONTROL-illumina.vcf.gz
    def vcf_pattern = ~/^([^\-]+)\-([^\-]+)(?:\-([^\-]+))?.*\.vcf(?:\.gz)?$/

    // Dosya temel adını al (yol yok)
    def basename = filename.getName()

    // Her deseni dene
    def matcher = fastq_pattern.matcher(basename)
    if (matcher.matches()) {
        return [
            id: matcher[0][1],
            condition: matcher[0][2],
            replicate: matcher[0][3] ? matcher[0][3].toInteger() : 1,
            read_number: matcher[0][4] ? matcher[0][4].toInteger() : null,
            file_type: 'fastq'
        ]
    }

    matcher = bam_pattern.matcher(basename)
    if (matcher.matches()) {
        return [
            id: matcher[0][1],
            condition: matcher[0][2],
            replicate: matcher[0][3] ? matcher[0][3].toInteger() : 1,
            file_type: 'bam'
        ]
    }

    matcher = vcf_pattern.matcher(basename)
    if (matcher.matches()) {
        return [
            id: matcher[0][1],
            condition: matcher[0][2],
            platform: matcher[0][3],
            file_type: 'vcf'
        ]
    }

    // Eğer hiçbir kalıp eşleşmezse, en azından dosya tipi ve adı döndür
    def extension = basename.lastIndexOf('.') > 0 ?
        basename.substring(basename.lastIndexOf('.') + 1) : 'unknown'

    log.warn "Bilinen şablonla eşleşmeyen dosya adı: $basename"
    return [
        id: basename.replaceAll(/\.[^\.]+$/, ""),
        file_type: extension
    ]
}
```

Bu desenleri bir iş akışında nasıl kullanacağınıza bir örnek:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Genom ayrıştırma fonksiyonunu içe aktar
include { parseGenomicFilename } from './modules/genome_parser'

// Farklı dosya türlerine göre işleme
workflow {
    // Tüm dosya türleri için kanal oluştur
    Channel
        .fromPath('data/*')
        .map { file ->
            def meta = parseGenomicFilename(file)
            tuple(meta, file)
        }
        .branch {
            fastq: it[0].file_type == 'fastq'
            bam: it[0].file_type == 'bam'
            vcf: it[0].file_type == 'vcf'
            other: true
        }
        .set { ch_files }

    // Her tür için loglama yap
    ch_files.fastq.view { meta, file -> "FASTQ: ${meta.id} (${meta.condition}), R${meta.read_number}" }
    ch_files.bam.view { meta, file -> "BAM: ${meta.id} (${meta.condition})" }
    ch_files.vcf.view { meta, file -> "VCF: ${meta.id} (${meta.condition})" }
    ch_files.other.view { meta, file -> "BAŞKA: ${file}" }
}
```

Bu yaklaşım, farklı dosya adı formatlarını ayrıştırmak ve bunları aynı iş akışına koyabilmek için çok daha sağlam bir çözüm sunar.

### 4.4. Özet: Pratik Nextflow Patternleri

Tebrikler! Şimdi, gerçek dünya Nextflow iş akışlarında sıklıkla kullanılan bir dizi önemli kodlama patternini ele aldınız:

1. **Veri Manipülasyonu**

   - Map ve koleksiyon işleme
   - Metadata organizasyonu
   - Dosya yollarını tuple ve metadata ile eşleştirme

2. **String İşleme**

   - Değişken enterpolasyonu (Nextflow, shell ve komut ikamesi)
   - Düzenli ifadeler (regex) ile karmaşık şekilde adlandırılmış dosyaları ayrıştırma
   - Çok satırlı dizeleri ve heredoc sözdizimini kullanma

3. **Kod Organizasyonu**

   - Yeniden kullanılabilir fonksiyonlar oluşturma
   - Sınıf benzeri yapılarla ilgili fonksiyonları gruplama
   - İş akışı organizasyonu için closure'ları kullanma

4. **Koşullu ve Dinamik Davranış**
   - Güvenli gezinme ve Elvis operatörleri
   - Ortama göre adapte olan iş akışları
   - Metadata'ya göre dinamik kaynak tahsisi

Bu konseptleri ustalaştığınızda, basit bir CSV işleme iş akışından kapsamlı bir genomik pipeline'a geçişinizi gördünüz. Bu kalıpları kendi iş akışlarınızda kullanarak, daha sağlam, okunaklı ve bakımı yapılabilir kod yazabilirsiniz.

### Özet

Bu yan görevde, Nextflow'un altındaki dilin temellerini derinlemesine inceledik. Şimdi, veri yapılarını çıkarabilir, dönüştürebilir ve manipüle edebilir, karmaşık dosya adı şemalarını işleyebilir, koşullu mantığı uygulayabilir, kodu yeniden kullanılabilir fonksiyonlara ve closure'lara çıkarabilir ve pipeline'larınızı gerçek dünya koşullarına adapte edebilirsiniz.

Bu becerileri uygulamak, herhangi bir Nextflow workflow'unu daha sağlam, daha iyi organize edilmiş ve bakımı daha kolay hale getirecektir. Scripting yetenekleriniz dataflow organizasyonuyla el ele gittiğinde, neredeyse her türlü veri işleme gereksinimini karşılayabilen, gerçekten güçlü iş akışları oluşturabilirsiniz.

## Bir sonraki adımlar

- **[nf-core](https://nf-co.re/)** ve Nextflow'un gerçek dünya uygulamalarını keşfetmeye devam edin.
- **[Standart kütüphane](https://nextflow.io/docs/latest/script.html)** de kullanılabilecek başka fonksiyonlara ve özelliklere göz atın.
- Gerçek iş akışları oluşturmak için **[DSL2 referans](https://nextflow.io/docs/latest/dsl2.html)** belgesini daha derinlemesine inceleyin.
- Diğer yan görevlere göz atın; **[Splitting & Grouping](./splitting_and_grouping.md)** ve **[Metadata](./metadata.md)** genellikle bu temel scripting kavramlarıyla iyi bir şekilde birlikte kullanılır.
