---
title: Nextflow para Genômica
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Escrever um fluxo de trabalho linear para aplicar chamada de variantes a uma única amostra
    - Manipular arquivos acessórios, como arquivos de índice e recursos de genoma de referência, adequadamente
    - Aproveitar o paradigma de dataflow do Nextflow para paralelizar a chamada de variantes por amostra
    - Implementar chamada conjunta multi-amostra usando operadores de canal relevantes
  audience_prerequisites:
    - "**Público:** Este curso é destinado a pesquisadores em genômica e áreas relacionadas que desejam desenvolver ou personalizar pipelines de análise de dados."
    - "**Habilidades:** Assume-se alguma familiaridade com a linha de comando, conceitos básicos de script e formatos de arquivo comuns de genômica."
    - "**Pré-requisitos:** Conceitos e ferramentas fundamentais do Nextflow cobertos em [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow para Genômica

**Um curso prático aplicando Nextflow a um caso de uso real de genômica: chamada de variantes com GATK.**

Este curso se baseia no treinamento para iniciantes [Hello Nextflow](../../hello_nextflow/) e demonstra como usar o Nextflow no contexto específico do domínio de genômica.
Você implementará um pipeline de chamada de variantes com [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), um pacote de software amplamente utilizado para analisar dados de sequenciamento de alto rendimento.

<!-- additional_information -->

## Visão geral do curso

Este curso é prático, com exercícios orientados a objetivos estruturados para introduzir informações gradualmente.

Você começará executando as ferramentas de chamada de variantes manualmente no terminal para entender a metodologia, depois construirá progressivamente um pipeline Nextflow que automatiza e escala a análise.

### Plano de aulas

Dividimos isso em três partes que focam em aspectos específicos da aplicação do Nextflow a um caso de uso de genômica.

| Capítulo do curso                                                               | Resumo                                                                                                      | Duração estimada |
| ------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- | ---------------- |
| [Parte 1: Visão geral do método](./01_method.md)                                | Compreendendo a metodologia de chamada de variantes e executando as ferramentas manualmente                 | 30 min           |
| [Parte 2: Chamada de variantes por amostra](./02_per_sample_variant_calling.md) | Construindo um pipeline que indexa arquivos BAM e chama variantes, depois escalando para múltiplas amostras | 60 min           |
| [Parte 3: Chamada conjunta em uma coorte](./03_joint_calling.md)                | Adicionando genotipagem conjunta multi-amostra usando operadores de canal para agregar saídas por amostra   | 45 min           |

Ao final deste curso, você será capaz de aplicar conceitos e ferramentas fundamentais do Nextflow a um caso de uso típico de genômica.

Pronto para fazer o curso?

[Começar :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
