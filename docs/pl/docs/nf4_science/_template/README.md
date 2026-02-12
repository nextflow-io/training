# Szablon kursu nf4_science

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

To jest ogólny szablon do tworzenia nowego kursu specyficznego dla danej dziedziny w ramach `nf4_science/`.
Opiera się na wspólnych wzorcach znalezionych w kursach Genomics i RNAseq.

## Jak używać

1. Skopiuj ten katalog i odpowiadający mu katalog skryptów, aby utworzyć swój kurs:

   ```bash
   # Dokumentacja
   cp -r docs/en/docs/nf4_science/_template docs/en/docs/nf4_science/{your_domain}

   # Skrypty
   cp -r nf4-science/_template nf4-science/{your_domain}
   ```

2. Zmień nazwę skryptu workflow'a:

   ```bash
   mv nf4-science/{your_domain}/_template.nf nf4-science/{your_domain}/{your_domain}.nf
   ```

3. Zmień nazwy plików modułów, aby pasowały do Twoich narzędzi:

   ```bash
   mv nf4-science/{your_domain}/modules/_tool_a.nf nf4-science/{your_domain}/modules/{tool_a}.nf
   mv nf4-science/{your_domain}/modules/_tool_b.nf nf4-science/{your_domain}/modules/{tool_b}.nf
   ```

4. Wyszukaj wszystkie wartości `{PLACEHOLDER}` w dokumentacji i skryptach, a następnie zastąp je swoją treścią.

5. Dodaj do nawigacji w `mkdocs.yml` (zobacz fragment poniżej).

6. Wypełnij katalog `data/` zestawami danych testowych.

7. Uzupełnij katalog `solutions/` działającym kodem dla każdej części.

## Fragment nawigacji mkdocs.yml

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

## Wykaz symboli zastępczych

### Symbole zastępcze na poziomie kursu

| Symbol zastępczy             | Opis                                  | Przykład (Genomics)                        |
| ---------------------------- | ------------------------------------- | ------------------------------------------ |
| `{DOMAIN}`                   | Nazwa dziedziny (z wielkimi literami) | Genomics                                   |
| `{DOMAIN_DIR}`               | Nazwa katalogu (małe litery)          | genomics                                   |
| `{METHOD}`                   | Nazwa metody analizy                  | variant calling                            |
| `{METHOD_SHORT_DESCRIPTION}` | Jednoliniowy opis metody              | variant calling with GATK                  |
| `{ACCESSORY_FILES}`          | Typy używanych plików pomocniczych    | index files and reference genome resources |

### Symbole zastępcze narzędzi

| Symbol zastępczy                      | Opis                                | Przykład (Genomics)                                 |
| ------------------------------------- | ----------------------------------- | --------------------------------------------------- |
| `{TOOL_A}` / `{TOOL_B}`               | Nazwy wyświetlane narzędzi          | Samtools / GATK                                     |
| `{TOOL_A_URL}` / `{TOOL_B_URL}`       | Adresy URL dokumentacji narzędzi    | https://www.htslib.org/                             |
| `{TOOL_A_CONTAINER_URI}`              | Pełny URI kontenera                 | community.wave.seqera.io/library/samtools:1.20--... |
| `{TOOL_A_MODULE}` / `{TOOL_B_MODULE}` | Nazwy plików modułów (bez .nf)      | samtools_index / gatk_haplotypecaller               |
| `{TOOL_A_PROCESS_NAME}`               | Nazwa procesu WIELKIMI LITERAMI     | SAMTOOLS_INDEX                                      |
| `{TOOL_A_COMMAND}`                    | Polecenie powłoki do uruchomienia   | samtools index '<input_bam>'                        |

### Symbole zastępcze wejścia/wyjścia

| Symbol zastępczy       | Opis                                  | Przykład (Genomics)  |
| ---------------------- | ------------------------------------- | -------------------- |
| `{PRIMARY_INPUT_TYPE}` | Główny typ pliku wejściowego          | BAM file             |
| `{PRIMARY_PARAM_NAME}` | Nazwa parametru Nextflow'a            | reads_bam            |
| `{TEST_INPUT_PATH}`    | Ścieżka do danych testowych w data/   | bam/reads_mother.bam |
| `{FINAL_OUTPUT_TYPE}`  | Końcowy typ wyjścia pipeline'u        | joint-called VCFs    |

### Symbole zastępcze treści

| Symbol zastępczy        | Opis                                                |
| ----------------------- | --------------------------------------------------- |
| `{PART2_SUMMARY}`       | Jednoliniowe podsumowanie zakresu Części 2          |
| `{AGGREGATION_SUMMARY}` | Jednoliniowe podsumowanie kroku agregacji           |
| `{PROCESS_LIST}`        | Lista nazw procesów oddzielonych przecinkami        |
| `{TYPEFORM_ID}`         | Identyfikator osadzenia Typeform dla ankiety        |

## Struktura pedagogiczna

Szablon podąża za trójczęściowym łukiem:

1. **Część 1 (01_method.md)**: Ręczne testowanie w kontenerach Docker w celu zrozumienia metodologii
2. **Część 2 (02_single_sample.md)**: Opakowanie poleceń w Nextflow'a; pojedyncza próbka, następnie przetwarzanie wsadowe
3. **Część 3 (03_multi_sample.md)**: Dodanie agregacji wielu próbek przy użyciu operatorów kanałów

### Kluczowe konwencje

- Używaj **bloków kodu z zakładkami Przed/Po** z `hl_lines`, aby pokazać zmiany w kodzie
- Używaj `??? success "Command output"` dla zwijanego oczekiwanego wyjścia
- Używaj `??? abstract "Directory contents"` dla zwijanego drzewa katalogów
- Kończ każdą główną sekcję podsekcjami **Podsumowanie** i **Co dalej?**
- Używaj poziomych linii `---` do oddzielania głównych numerowanych sekcji
- Dostarczaj pliki szkieletowe (workflow + moduły) do wypełnienia przez uczestników
- Organizuj rozwiązania według części (`solutions/part2/`, `solutions/part3/`)

## Struktura katalogów

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
