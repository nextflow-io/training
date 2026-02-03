---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Executar e gerenciar a execução de fluxos de trabalho Nextflow
    - Encontrar e interpretar saídas (resultados) e arquivos de log
    - Reconhecer os componentes principais do Nextflow em um fluxo de trabalho simples de múltiplas etapas
    - Configurar a execução de pipelines para rodar em plataformas computacionais comuns, incluindo HPC e nuvem
    - Resumir as melhores práticas para reprodutibilidade, portabilidade e reutilização de código que tornam os pipelines FAIR, incluindo modularidade de código e contêineres de software
  audience_prerequisites:
    - "**Público:** Este curso é projetado para alunos que são completamente novos no Nextflow e querem executar pipelines existentes."
    - "**Habilidades:** Alguma familiaridade com a linha de comando, conceitos básicos de script e formatos de arquivo comuns é assumida."
    - "**Domínio:** Os exercícios são todos independentes de domínio, então nenhum conhecimento científico prévio é necessário."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run é uma introdução prática para executar análises de dados reproduzíveis e escaláveis.**

Trabalhando através de exemplos práticos e exercícios guiados, você aprenderá os fundamentos do uso do Nextflow, incluindo como executar pipelines, gerenciar arquivos e dependências de software, paralelizar a execução sem esforço e executar fluxos de trabalho em diferentes ambientes computacionais.

Você levará consigo as habilidades e confiança para começar a executar fluxos de trabalho com Nextflow.

<!-- additional_information -->

## Visão geral do curso

### O que você fará

Este curso é prático, com exercícios orientados a objetivos estruturados para introduzir informações gradualmente.

Você executará várias versões de um pipeline Nextflow que processa entradas de texto.
Você começará com uma versão simples que consiste em uma única etapa e eventualmente progredirá para uma versão de múltiplas etapas que recebe um arquivo CSV de entradas de texto tabulares, executa algumas etapas de transformação e produz um único arquivo de texto contendo uma imagem ASCII de um personagem dizendo o texto transformado.

Este curso foca na execução de pipelines (nomeado após o comando principal `nextflow run`).
Se você está procurando uma introdução ao desenvolvimento de pipelines Nextflow, veja [Hello Nextflow](../hello_nextflow/index.md).

### Plano de aula

Dividimos isso em três partes que focarão em aspectos específicos da execução e gerenciamento de pipelines escritos em Nextflow.

| Capítulo do curso                                        | Resumo                                                                                                                             | Duração estimada |
| -------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- | ---------------- |
| [Parte 1: Operações básicas de execução](./01_basics.md) | Executar e gerenciar a execução de um fluxo de trabalho simples                                                                    | 30 mins          |
| [Parte 2: Executar pipelines reais](./02_pipeline.md)    | Processar entradas complexas, executar fluxos de trabalho de múltiplas etapas, usar contêineres e paralelizar execução sem esforço | 60 mins          |
| [Parte 3: Configuração de execução](./03_config.md)      | Personalizar o comportamento do pipeline e otimizar o uso em diferentes ambientes computacionais                                   | 60 mins          |

Ao final deste curso, você estará bem preparado para enfrentar os próximos passos em sua jornada para executar fluxos de trabalho reproduzíveis para suas necessidades de computação científica.

Pronto para fazer o curso?

[Começar a aprender :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
