# Plantilla de curs nf4_science

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Aquesta és una plantilla genèrica per crear un nou curs específic de domini sota `nf4_science/`.
Es basa en els patrons comuns trobats als cursos de Genòmica i RNAseq.

## Com utilitzar-la

1. Copia aquest directori i el directori de scripts corresponent per crear el teu curs:

   ```bash
   # Documentació
   cp -r docs/en/docs/nf4_science/_template docs/en/docs/nf4_science/{your_domain}

   # Scripts
   cp -r nf4-science/_template nf4-science/{your_domain}
   ```

2. Reanomena el script del workflow:

   ```bash
   mv nf4-science/{your_domain}/_template.nf nf4-science/{your_domain}/{your_domain}.nf
   ```

3. Reanomena els fitxers de mòdul perquè coincideixin amb les teves eines:

   ```bash
   mv nf4-science/{your_domain}/modules/_tool_a.nf nf4-science/{your_domain}/modules/{tool_a}.nf
   mv nf4-science/{your_domain}/modules/_tool_b.nf nf4-science/{your_domain}/modules/{tool_b}.nf
   ```

4. Cerca tots els valors `{PLACEHOLDER}` a la documentació i els scripts i substitueix-los pel teu contingut.

5. Afegeix-ho a la navegació de `mkdocs.yml` (consulta el fragment a continuació).

6. Omple el directori `data/` amb conjunts de dades de prova.

7. Construeix el directori `solutions/` amb codi funcional per a cada part.

## Fragment de navegació de mkdocs.yml

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

## Referència de marcadors de posició

### Marcadors de posició a nivell de curs

| Marcador de posició          | Descripció                             | Exemple (Genòmica)                         |
| ---------------------------- | -------------------------------------- | ------------------------------------------ |
| `{DOMAIN}`                   | Nom del domini (majúscules inicials)   | Genomics                                   |
| `{DOMAIN_DIR}`               | Nom del directori (minúscules)         | genomics                                   |
| `{METHOD}`                   | Nom del mètode d'anàlisi               | variant calling                            |
| `{METHOD_SHORT_DESCRIPTION}` | Descripció d'una línia del mètode      | variant calling with GATK                  |
| `{ACCESSORY_FILES}`          | Tipus de fitxers accessoris utilitzats | index files and reference genome resources |

### Marcadors de posició d'eines

| Marcador de posició                   | Descripció                          | Exemple (Genòmica)                                  |
| ------------------------------------- | ----------------------------------- | --------------------------------------------------- |
| `{TOOL_A}` / `{TOOL_B}`               | Noms de visualització de les eines  | Samtools / GATK                                     |
| `{TOOL_A_URL}` / `{TOOL_B_URL}`       | URLs de documentació de les eines   | https://www.htslib.org/                             |
| `{TOOL_A_CONTAINER_URI}`              | URI complet del contenidor          | community.wave.seqera.io/library/samtools:1.20--... |
| `{TOOL_A_MODULE}` / `{TOOL_B_MODULE}` | Noms de fitxer de mòdul (sense .nf) | samtools_index / gatk_haplotypecaller               |
| `{TOOL_A_PROCESS_NAME}`               | Nom del procés en MAJÚSCULES        | SAMTOOLS_INDEX                                      |
| `{TOOL_A_COMMAND}`                    | La comanda shell a executar         | samtools index '<input_bam>'                        |

### Marcadors de posició d'entrada/sortida

| Marcador de posició    | Descripció                          | Exemple (Genòmica)   |
| ---------------------- | ----------------------------------- | -------------------- |
| `{PRIMARY_INPUT_TYPE}` | Tipus de fitxer d'entrada principal | BAM file             |
| `{PRIMARY_PARAM_NAME}` | Nom del paràmetre Nextflow          | reads_bam            |
| `{TEST_INPUT_PATH}`    | Ruta a l'entrada de prova a data/   | bam/reads_mother.bam |
| `{FINAL_OUTPUT_TYPE}`  | Tipus de sortida final del pipeline | joint-called VCFs    |

### Marcadors de posició de contingut

| Marcador de posició     | Descripció                                     |
| ----------------------- | ---------------------------------------------- |
| `{PART2_SUMMARY}`       | Resum d'una línia de l'abast de la Part 2      |
| `{AGGREGATION_SUMMARY}` | Resum d'una línia del pas d'agregació          |
| `{PROCESS_LIST}`        | Llista separada per comes de noms de processos |
| `{TYPEFORM_ID}`         | ID d'incrustació de Typeform per a l'enquesta  |

## Estructura pedagògica

La plantilla segueix un arc de tres parts:

1. **Part 1 (01_method.md)**: Proves manual en contenidors Docker per entendre la metodologia
2. **Part 2 (02_single_sample.md)**: Embolcallar comandes en Nextflow; mostra única, després per lots
3. **Part 3 (03_multi_sample.md)**: Afegir agregació multi-mostra utilitzant operadors de canal

### Convencions clau

- Utilitza **blocs de codi amb pestanyes Abans/Després** amb `hl_lines` per mostrar canvis de codi
- Utilitza `??? success "Command output"` per a sortida esperada desplegable
- Utilitza `??? abstract "Directory contents"` per a arbres de directoris desplegables
- Finalitza cada secció principal amb subseccions **Conclusió** i **Què segueix?**
- Utilitza regles horitzontals `---` per separar seccions principals numerades
- Proporciona fitxers esquelet (workflow + mòduls) perquè els estudiants els omplin
- Organitza les solucions per part (`solutions/part2/`, `solutions/part3/`)

## Estructura de directoris

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
