# Template del Corso nf4_science

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Questo è un template generico per creare un nuovo corso specifico per dominio sotto `nf4_science/`.
È basato sui pattern comuni trovati nei corsi di Genomics e RNAseq.

## Come utilizzarlo

1. Copiate questa directory e la directory degli script corrispondente per creare il vostro corso:

   ```bash
   # Docs
   cp -r docs/en/docs/nf4_science/_template docs/en/docs/nf4_science/{your_domain}

   # Scripts
   cp -r nf4-science/_template nf4-science/{your_domain}
   ```

2. Rinominate lo script del flusso di lavoro:

   ```bash
   mv nf4-science/{your_domain}/_template.nf nf4-science/{your_domain}/{your_domain}.nf
   ```

3. Rinominate i file dei moduli per corrispondere ai vostri strumenti:

   ```bash
   mv nf4-science/{your_domain}/modules/_tool_a.nf nf4-science/{your_domain}/modules/{tool_a}.nf
   mv nf4-science/{your_domain}/modules/_tool_b.nf nf4-science/{your_domain}/modules/{tool_b}.nf
   ```

4. Cercate tutti i valori `{PLACEHOLDER}` nei documenti e negli script e sostituiteli con i vostri contenuti.

5. Aggiungete al nav di `mkdocs.yml` (vedete lo snippet qui sotto).

6. Popolate la directory `data/` con dataset di test.

7. Costruite la directory `solutions/` con codice funzionante per ogni parte.

## Snippet nav per mkdocs.yml

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

## Riferimento dei placeholder

### Placeholder a livello di corso

| Placeholder                  | Descrizione                                | Esempio (Genomics)                         |
| ---------------------------- | ------------------------------------------ | ------------------------------------------ |
| `{DOMAIN}`                   | Nome del dominio (title case)              | Genomics                                   |
| `{DOMAIN_DIR}`               | Nome della directory (minuscolo)           | genomics                                   |
| `{METHOD}`                   | Nome del metodo di analisi                 | variant calling                            |
| `{METHOD_SHORT_DESCRIPTION}` | Descrizione del metodo in una riga         | variant calling with GATK                  |
| `{ACCESSORY_FILES}`          | Tipi di file accessori utilizzati          | index files and reference genome resources |

### Placeholder degli strumenti

| Placeholder                           | Descrizione                       | Esempio (Genomics)                                  |
| ------------------------------------- | --------------------------------- | --------------------------------------------------- |
| `{TOOL_A}` / `{TOOL_B}`               | Nomi visualizzati degli strumenti | Samtools / GATK                                     |
| `{TOOL_A_URL}` / `{TOOL_B_URL}`       | URL della documentazione          | https://www.htslib.org/                             |
| `{TOOL_A_CONTAINER_URI}`              | URI completo del container        | community.wave.seqera.io/library/samtools:1.20--... |
| `{TOOL_A_MODULE}` / `{TOOL_B_MODULE}` | Nomi dei file modulo (senza .nf)  | samtools_index / gatk_haplotypecaller               |
| `{TOOL_A_PROCESS_NAME}`               | Nome del processo in MAIUSCOLO    | SAMTOOLS_INDEX                                      |
| `{TOOL_A_COMMAND}`                    | Il comando shell da eseguire      | samtools index '<input_bam>'                        |

### Placeholder di input/output

| Placeholder            | Descrizione                       | Esempio (Genomics)   |
| ---------------------- | --------------------------------- | -------------------- |
| `{PRIMARY_INPUT_TYPE}` | Tipo di file di input primario    | BAM file             |
| `{PRIMARY_PARAM_NAME}` | Nome del parametro Nextflow       | reads_bam            |
| `{TEST_INPUT_PATH}`    | Percorso all'input di test in data/ | bam/reads_mother.bam |
| `{FINAL_OUTPUT_TYPE}`  | Tipo di output finale della pipeline | joint-called VCFs    |

### Placeholder dei contenuti

| Placeholder             | Descrizione                                      |
| ----------------------- | ------------------------------------------------ |
| `{PART2_SUMMARY}`       | Riepilogo in una riga dell'ambito della Parte 2  |
| `{AGGREGATION_SUMMARY}` | Riepilogo in una riga del passaggio di aggregazione |
| `{PROCESS_LIST}`        | Lista separata da virgole dei nomi dei processi  |
| `{TYPEFORM_ID}`         | ID di embed Typeform per il sondaggio            |

## Struttura pedagogica

Il template segue un arco in tre parti:

1. **Parte 1 (01_method.md)**: Test manuale nei container Docker per comprendere la metodologia
2. **Parte 2 (02_single_sample.md)**: Incapsulare i comandi in Nextflow; singolo campione, poi batch
3. **Parte 3 (03_multi_sample.md)**: Aggiungere aggregazione multi-campione usando operatori di canali

### Convenzioni chiave

- Usate **blocchi di codice con tab Before/After** con `hl_lines` per mostrare le modifiche al codice
- Usate `??? success "Command output"` per output atteso comprimibile
- Usate `??? abstract "Directory contents"` per alberi di directory comprimibili
- Terminate ogni sezione principale con sottosezioni **Takeaway** e **Cosa c'è dopo?**
- Usate linee orizzontali `---` per separare le sezioni numerate principali
- Fornite file scheletro (flusso di lavoro + moduli) da completare per gli studenti
- Organizzate le soluzioni per parte (`solutions/part2/`, `solutions/part3/`)

## Struttura delle directory

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
