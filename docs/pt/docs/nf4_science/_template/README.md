# Modelo de Curso nf4_science

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Este é um modelo genérico para criar um novo curso específico de domínio em `nf4_science/`.
Ele é baseado nos padrões comuns encontrados nos cursos de Genômica e RNAseq.

## Como usar

1. Copie este diretório e o diretório de scripts correspondente para criar seu curso:

   ```bash
   # Documentação
   cp -r docs/en/docs/nf4_science/_template docs/en/docs/nf4_science/{your_domain}

   # Scripts
   cp -r nf4-science/_template nf4-science/{your_domain}
   ```

2. Renomeie o script do fluxo de trabalho:

   ```bash
   mv nf4-science/{your_domain}/_template.nf nf4-science/{your_domain}/{your_domain}.nf
   ```

3. Renomeie os arquivos de módulo para corresponder às suas ferramentas:

   ```bash
   mv nf4-science/{your_domain}/modules/_tool_a.nf nf4-science/{your_domain}/modules/{tool_a}.nf
   mv nf4-science/{your_domain}/modules/_tool_b.nf nf4-science/{your_domain}/modules/{tool_b}.nf
   ```

4. Procure por todos os valores `{PLACEHOLDER}` na documentação e scripts e substitua-os pelo seu conteúdo.

5. Adicione ao nav do `mkdocs.yml` (veja o trecho abaixo).

6. Preencha o diretório `data/` com conjuntos de dados de teste.

7. Construa o diretório `solutions/` com código funcional para cada parte.

## Trecho do nav do mkdocs.yml

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

## Referência de placeholders

### Placeholders de nível de curso

| Placeholder                  | Descrição                                | Exemplo (Genômica)                                      |
| ---------------------------- | ---------------------------------------- | ------------------------------------------------------- |
| `{DOMAIN}`                   | Nome do domínio (título em maiúsculas)   | Genomics                                                |
| `{DOMAIN_DIR}`               | Nome do diretório (minúsculas)           | genomics                                                |
| `{METHOD}`                   | Nome do método de análise                | variant calling                                         |
| `{METHOD_SHORT_DESCRIPTION}` | Descrição do método em uma linha         | variant calling with GATK                               |
| `{ACCESSORY_FILES}`          | Tipos de arquivos acessórios utilizados  | index files and reference genome resources              |

### Placeholders de ferramentas

| Placeholder                           | Descrição                          | Exemplo (Genômica)                                  |
| ------------------------------------- | ---------------------------------- | --------------------------------------------------- |
| `{TOOL_A}` / `{TOOL_B}`               | Nomes de exibição das ferramentas  | Samtools / GATK                                     |
| `{TOOL_A_URL}` / `{TOOL_B_URL}`       | URLs da documentação das ferramentas | https://www.htslib.org/                             |
| `{TOOL_A_CONTAINER_URI}`              | URI completo do contêiner          | community.wave.seqera.io/library/samtools:1.20--... |
| `{TOOL_A_MODULE}` / `{TOOL_B_MODULE}` | Nomes dos arquivos de módulo (sem .nf) | samtools_index / gatk_haplotypecaller               |
| `{TOOL_A_PROCESS_NAME}`               | Nome do processo em MAIÚSCULAS     | SAMTOOLS_INDEX                                      |
| `{TOOL_A_COMMAND}`                    | O comando shell a ser executado    | samtools index '<input_bam>'                        |

### Placeholders de entrada/saída

| Placeholder            | Descrição                           | Exemplo (Genômica)   |
| ---------------------- | ----------------------------------- | -------------------- |
| `{PRIMARY_INPUT_TYPE}` | Tipo de arquivo de entrada primário | BAM file             |
| `{PRIMARY_PARAM_NAME}` | Nome do parâmetro Nextflow          | reads_bam            |
| `{TEST_INPUT_PATH}`    | Caminho para entrada de teste em data/ | bam/reads_mother.bam |
| `{FINAL_OUTPUT_TYPE}`  | Tipo de saída final do pipeline     | joint-called VCFs    |

### Placeholders de conteúdo

| Placeholder             | Descrição                                      |
| ----------------------- | ---------------------------------------------- |
| `{PART2_SUMMARY}`       | Resumo de uma linha do escopo da Parte 2       |
| `{AGGREGATION_SUMMARY}` | Resumo de uma linha da etapa de agregação      |
| `{PROCESS_LIST}`        | Lista separada por vírgulas de nomes de processos |
| `{TYPEFORM_ID}`         | ID de incorporação do Typeform para a pesquisa |

## Estrutura pedagógica

O modelo segue um arco de três partes:

1. **Parte 1 (01_method.md)**: Teste manual em contêineres Docker para entender a metodologia
2. **Parte 2 (02_single_sample.md)**: Encapsular comandos em Nextflow; amostra única, depois em lote
3. **Parte 3 (03_multi_sample.md)**: Adicionar agregação de múltiplas amostras usando operadores de canal

### Convenções principais

- Use **blocos de código com abas Antes/Depois** com `hl_lines` para mostrar mudanças no código
- Use `??? success "Command output"` para saída esperada recolhível
- Use `??? abstract "Directory contents"` para árvores de diretório recolhíveis
- Termine cada seção principal com subseções **Conclusão** e **O que vem a seguir?**
- Use réguas horizontais `---` para separar seções numeradas principais
- Forneça arquivos esqueleto (fluxo de trabalho + módulos) para os alunos preencherem
- Organize as soluções por parte (`solutions/part2/`, `solutions/part3/`)

## Estrutura de diretórios

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
