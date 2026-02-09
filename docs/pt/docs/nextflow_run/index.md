---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Iniciar e gerenciar a execução de fluxos de trabalho Nextflow
    - Encontrar e interpretar saídas (resultados) e arquivos de log
    - Reconhecer componentes principais do Nextflow em um fluxo de trabalho simples de múltiplas etapas
    - Configurar a execução de pipelines para rodar em plataformas de computação comuns, incluindo HPC e nuvem
    - Resumir as melhores práticas para reprodutibilidade, portabilidade e reutilização de código que tornam pipelines FAIR, incluindo modularidade de código e contêineres de software
  audience_prerequisites:
    - "**Público:** Este curso é projetado para alunos que são completamente novos no Nextflow e desejam executar pipelines existentes."
    - "**Habilidades:** Assume-se alguma familiaridade com a linha de comando, conceitos básicos de script e formatos de arquivo comuns."
    - "**Domínio:** Os exercícios são todos agnósticos de domínio, portanto nenhum conhecimento científico prévio é necessário."
---

# Nextflow Run

**Nextflow Run é uma introdução prática à execução de fluxos de trabalho de análise de dados reprodutíveis e escaláveis.**

Trabalhando através de exemplos práticos e exercícios guiados, você aprenderá os fundamentos do uso do Nextflow, incluindo como executar pipelines, gerenciar arquivos e dependências de software, paralelizar a execução sem esforço e executar fluxos de trabalho em diferentes ambientes computacionais.

Você adquirirá as habilidades e a confiança para começar a executar fluxos de trabalho com Nextflow.

<!-- additional_information -->

## Visão geral do curso

### O que você fará

Este curso é prático, com exercícios orientados a objetivos estruturados para introduzir informações gradualmente.

Você executará várias versões de um pipeline Nextflow que processa entradas de texto.
Você começará com uma versão simples que consiste em uma única etapa e, eventualmente, progredirá para uma versão de múltiplas etapas que recebe um arquivo CSV de entradas de texto tabulares, executa algumas etapas de transformação e gera um único arquivo de texto contendo uma imagem ASCII de um personagem dizendo o texto transformado.

Este curso foca na execução de pipelines (nomeado após o comando principal `nextflow run`).
Se você está procurando uma introdução ao desenvolvimento de pipelines Nextflow, veja [Hello Nextflow](../hello_nextflow/index.md).

### Plano de aula

Dividimos isso em três partes que focarão em aspectos específicos da execução e gerenciamento de pipelines escritos em Nextflow.

| Capítulo do curso                                  | Resumo                                                                                                                                  | Duração estimada |
| -------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- | ---------------- |
| [Parte 1: Operações básicas](./01_basics.md)       | Iniciando e gerenciando a execução de um fluxo de trabalho simples                                                                      | 30 min           |
| [Parte 2: Executar pipelines reais](./02_pipeline.md) | Processando entradas complexas, executando fluxos de trabalho de múltiplas etapas, usando contêineres e paralelizando a execução sem esforço | 60 min           |
| [Parte 3: Configuração de execução](./03_config.md)   | Personalizando o comportamento do pipeline e otimizando o uso em diferentes ambientes computacionais                                    | 60 min           |

Ao final deste curso, você estará bem preparado para enfrentar os próximos passos em sua jornada para executar fluxos de trabalho reprodutíveis para suas necessidades de computação científica.

Pronto para fazer o curso?

[Começar a aprender :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
