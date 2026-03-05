---
title: RNAseq
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Escrever um fluxo de trabalho linear para aplicar métodos básicos de processamento e QC de RNAseq
    - Manipular arquivos específicos do domínio, como FASTQ e recursos de genoma de referência, de forma apropriada
    - Manipular dados de sequenciamento single-end e paired-end
    - Aproveitar o paradigma de dataflow do Nextflow para paralelizar o processamento de RNAseq por amostra
    - Agregar relatórios de QC através de múltiplas etapas e amostras usando operadores de canal relevantes
  audience_prerequisites:
    - "**Público:** Este curso é destinado a pesquisadores em transcriptômica e áreas relacionadas que desejam desenvolver ou personalizar pipelines de análise de dados."
    - "**Habilidades:** Assume-se alguma familiaridade com a linha de comando, conceitos básicos de script e formatos de arquivo comuns de RNAseq."
    - "**Pré-requisitos:** Conceitos e ferramentas fundamentais do Nextflow abordados no [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow para RNAseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Um curso prático aplicando Nextflow a um caso de uso real de transcriptômica: processamento de RNAseq bulk com Trim Galore, HISAT2 e FastQC.**

Este curso se baseia no treinamento para iniciantes [Hello Nextflow](../../hello_nextflow/) e demonstra como usar Nextflow no contexto específico de análise de RNAseq bulk.
Você implementará um pipeline de processamento que remove sequências adaptadoras, alinha as reads a um genoma de referência e realiza controle de qualidade (QC) em várias etapas.

<!-- additional_information -->

## Visão geral do curso

Este curso é prático, com exercícios orientados a objetivos estruturados para introduzir informações gradualmente.

Você começará executando as ferramentas de processamento manualmente no terminal para entender a metodologia, depois construirá progressivamente um pipeline Nextflow que automatiza e escala a análise.

### Plano de aulas

Dividimos isso em três partes que focam em aspectos específicos da aplicação do Nextflow a um caso de uso de RNAseq.

| Capítulo do curso                                                       | Resumo                                                                                                                         | Duração estimada |
| ----------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ | ---------------- |
| [Parte 1: Visão geral do método](./01_method.md)                        | Compreendendo a metodologia de processamento de RNAseq e executando as ferramentas manualmente                                 | 30 min           |
| [Parte 2: Implementação de amostra única](./02_single-sample.md)        | Construindo um pipeline que remove adaptadores, alinha e faz QC de uma única amostra, depois escalando para múltiplas amostras | 60 min           |
| [Parte 3: Implementação multi-amostra paired-end](./03_multi-sample.md) | Estendendo o pipeline para manipular dados paired-end e agregar relatórios de QC através das amostras                          | 45 min           |

Ao final deste curso, você será capaz de aplicar conceitos e ferramentas fundamentais do Nextflow a um caso de uso típico de RNAseq.

Pronto para fazer o curso?

[Começar :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
