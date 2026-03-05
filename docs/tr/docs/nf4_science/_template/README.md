# nf4_science Kurs Şablonu

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu, `nf4_science/` altında yeni bir alana özgü kurs oluşturmak için genel bir şablondur.
Genomik ve RNAseq kurslarında bulunan yaygın kalıplara dayanmaktadır.

## Nasıl kullanılır

1. Kursunuzu oluşturmak için bu dizini ve ilgili betik dizinini kopyalayın:

   ```bash
   # Dokümanlar
   cp -r docs/en/docs/nf4_science/_template docs/en/docs/nf4_science/{your_domain}

   # Betikler
   cp -r nf4-science/_template nf4-science/{your_domain}
   ```

2. İş akışı betiğini yeniden adlandırın:

   ```bash
   mv nf4-science/{your_domain}/_template.nf nf4-science/{your_domain}/{your_domain}.nf
   ```

3. Modül dosyalarını araçlarınıza uyacak şekilde yeniden adlandırın:

   ```bash
   mv nf4-science/{your_domain}/modules/_tool_a.nf nf4-science/{your_domain}/modules/{tool_a}.nf
   mv nf4-science/{your_domain}/modules/_tool_b.nf nf4-science/{your_domain}/modules/{tool_b}.nf
   ```

4. Dokümanlarda ve betiklerde tüm `{PLACEHOLDER}` değerlerini arayın ve içeriğinizle değiştirin.

5. `mkdocs.yml` nav'a ekleyin (aşağıdaki kod parçasına bakın).

6. `data/` dizinini test veri kümeleriyle doldurun.

7. `solutions/` dizinini her bölüm için çalışan kodla oluşturun.

## mkdocs.yml nav kod parçası

```yaml
- { Your Domain }:
    - nf4_science/{your_domain}/index.md
    - nf4_science/{your_domain}/00_orientation.md
    - nf4_science/{your_domain}/01_method.md
    - nf4_science/{your_domain}/02_single_sample.md
    - nf4_science/{your_domain}/03_multi_sample.md
    - nf4_science/{your_domain}/survey.md
    - nf4_science/{your_domain}/next_steps.md
```

## Yer tutucu referansı

### Kurs düzeyinde yer tutucular

| Yer Tutucu                   | Açıklama                          | Örnek (Genomik)                            |
| ---------------------------- | --------------------------------- | ------------------------------------------ |
| `{DOMAIN}`                   | Alan adı (başlık harfleriyle)     | Genomics                                   |
| `{DOMAIN_DIR}`               | Dizin adı (küçük harflerle)       | genomics                                   |
| `{METHOD}`                   | Analiz yöntemi adı                | variant calling                            |
| `{METHOD_SHORT_DESCRIPTION}` | Tek satırlık yöntem açıklaması    | variant calling with GATK                  |
| `{ACCESSORY_FILES}`          | Kullanılan yardımcı dosya türleri | index files and reference genome resources |

### Araç yer tutucuları

| Yer Tutucu                            | Açıklama                         | Örnek (Genomik)                                     |
| ------------------------------------- | -------------------------------- | --------------------------------------------------- |
| `{TOOL_A}` / `{TOOL_B}`               | Araç görünen adları              | Samtools / GATK                                     |
| `{TOOL_A_URL}` / `{TOOL_B_URL}`       | Araç dokümantasyon URL'leri      | https://www.htslib.org/                             |
| `{TOOL_A_CONTAINER_URI}`              | Tam konteyner URI'si             | community.wave.seqera.io/library/samtools:1.20--... |
| `{TOOL_A_MODULE}` / `{TOOL_B_MODULE}` | Modül dosya adları (.nf olmadan) | samtools_index / gatk_haplotypecaller               |
| `{TOOL_A_PROCESS_NAME}`               | BÜYÜK HARF süreç adı             | SAMTOOLS_INDEX                                      |
| `{TOOL_A_COMMAND}`                    | Çalıştırılacak kabuk komutu      | samtools index '<input_bam>'                        |

### Girdi/Çıktı yer tutucuları

| Yer Tutucu             | Açıklama                       | Örnek (Genomik)      |
| ---------------------- | ------------------------------ | -------------------- |
| `{PRIMARY_INPUT_TYPE}` | Birincil girdi dosya türü      | BAM file             |
| `{PRIMARY_PARAM_NAME}` | Nextflow parametre adı         | reads_bam            |
| `{TEST_INPUT_PATH}`    | data/ içindeki test girdi yolu | bam/reads_mother.bam |
| `{FINAL_OUTPUT_TYPE}`  | Son pipeline çıktı türü        | joint-called VCFs    |

### İçerik yer tutucuları

| Yer Tutucu              | Açıklama                               |
| ----------------------- | -------------------------------------- |
| `{PART2_SUMMARY}`       | Bölüm 2 kapsamının tek satırlık özeti  |
| `{AGGREGATION_SUMMARY}` | Toplama adımının tek satırlık özeti    |
| `{PROCESS_LIST}`        | Virgülle ayrılmış süreç adları listesi |
| `{TYPEFORM_ID}`         | Anket için Typeform embed ID'si        |

## Pedagojik yapı

Şablon üç bölümlük bir yapıyı takip eder:

1. **Bölüm 1 (01_method.md)**: Metodolojiyi anlamak için Docker konteynerlerinde manuel test
2. **Bölüm 2 (02_single_sample.md)**: Komutları Nextflow'a sarma; tek örnek, ardından toplu işlem
3. **Bölüm 3 (03_multi_sample.md)**: Kanal operatörleri kullanarak çok örnekli toplama ekleme

### Temel kurallar

- Kod değişikliklerini göstermek için `hl_lines` ile **Önce/Sonra sekmeli kod blokları** kullanın
- Daraltılabilir beklenen çıktı için `??? success "Komut çıktısı"` kullanın
- Daraltılabilir dizin ağaçları için `??? abstract "Dizin içeriği"` kullanın
- Her ana bölümü **Özet** ve **Sırada ne var?** alt bölümleriyle bitirin
- Ana numaralı bölümleri ayırmak için `---` yatay çizgiler kullanın
- Öğrencilerin doldurması için iskelet dosyaları (iş akışı + modüller) sağlayın
- Çözümleri bölümlere göre düzenleyin (`solutions/part2/`, `solutions/part3/`)

## Dizin yapısı

```
docs/en/docs/nf4_science/{domain}/
├── index.md                    # Course overview with frontmatter
├── 00_orientation.md           # Environment setup
├── 01_method.md                # Manual testing in containers
├── 02_single_sample.md         # Single-sample Nextflow implementation
├── 03_multi_sample.md          # Multi-sample aggregation
├── survey.md                   # Typeform feedback survey
├── next_steps.md               # Course summary and suggestions
└── img/                        # Diagrams (.excalidraw.svg, .png)

nf4-science/{domain}/
├── {domain}.nf                 # Skeleton workflow file
├── nextflow.config             # Minimal config (docker.enabled = true)
├── data/                       # Test datasets and resources
│   └── samplesheet.csv         # Sample metadata
├── modules/                    # Skeleton module files
│   ├── {tool_a}.nf
│   └── {tool_b}.nf
└── solutions/
    ├── part2/                  # Complete Part 2 solution
    │   ├── {domain}-2.nf
    │   ├── nextflow.config
    │   └── modules/
    └── part3/                  # Complete Part 3 solution
        ├── {domain}-3.nf
        ├── nextflow.config
        └── modules/
```
