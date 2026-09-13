---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Iniciar e gerenciar pipelines Nextflow pela linha de comando
    - Entender como canais e operadores permitem fluxos de trabalho eficientes com múltiplas entradas e múltiplas etapas
    - Usar contêineres para gerenciar dependências de software e garantir reprodutibilidade
    - Configurar a execução de pipelines e saídas
    - Gerar relatórios de execução, inspecionar o histórico de execuções anteriores e limpar diretórios de trabalho antigos
    - Executar pipelines diretamente de repositórios remotos como o GitHub
  audience_prerequisites:
    - "**Público:** Este curso é projetado para alunos que são completamente novos no Nextflow e querem executar pipelines existentes."
    - "**Habilidades:** Alguma familiaridade com a linha de comando, conceitos básicos de script e formatos de arquivo comuns é assumida."
    - "**Domínio:** Os exercícios são todos independentes de domínio, então nenhum conhecimento científico prévio é necessário."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run é uma introdução prática para executar análises de dados reproduzíveis e escaláveis.**

Trabalhando através de uma série de exercícios orientados a objetivos, você aprenderá os fundamentos para iniciar e gerenciar pipelines Nextflow, entenderá como canais e operadores permitem o processamento paralelo de múltiplas entradas e usará contêineres para gerenciar dependências de software.

Você levará consigo as habilidades e confiança para começar a executar fluxos de trabalho com Nextflow.

<!-- additional_information -->

## Visão geral do curso

Este curso é prático, com exercícios orientados a objetivos estruturados para introduzir informações gradualmente.

Você executará várias versões de um pipeline Nextflow que processa entradas de texto, começando com uma versão simples de uma única etapa e progredindo para uma versão de múltiplas etapas que recebe um arquivo CSV de entradas, executa algumas etapas de transformação e produz um único arquivo de texto contendo arte ASCII gerada por uma ferramenta em contêiner.

Este curso foca na execução de pipelines (nomeado após o comando principal `nextflow run`).
Se você está procurando uma introdução ao desenvolvimento de pipelines Nextflow, veja [Hello Nextflow](../hello_nextflow/index.md).

!!! note "Nota"

    Procurando a versão anterior deste curso? Ela foi substituída pela versão nesta página, mas ainda pode ser acessada na [versão 3.6.1](https://training.nextflow.io/3.6.1/nextflow_run/) do site de treinamento.

### Plano de aula

| Capítulo do curso                                                               | Resumo                                                                                                                | Duração estimada |
| ------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- | ---------------- |
| [Parte 1: Executar Nextflow](./01_run_nextflow.md)                              | Iniciar e gerenciar pipelines Nextflow, e entender os mecanismos essenciais do fluxo de trabalho                      | 25 mins          |
| [Parte 2: Configurar o pipeline](./02_configure_pipeline.md)                    | Configurar a execução do pipeline e as saídas usando `nextflow.config`                                                | 20 mins          |
| [Parte 3: Gerenciar execuções de fluxos de trabalho](./03_manage_executions.md) | Gerar relatórios de execução, inspecionar o histórico de execuções anteriores e limpar diretórios de trabalho antigos | 10 mins          |
| [Parte 4: Executar pipelines remotos](./04_remote_repositories.md)              | Executar um pipeline diretamente do GitHub e fixá-lo em uma revisão específica                                        | 10 mins          |

Ao final deste curso, você estará bem preparado para enfrentar os próximos passos em sua jornada para executar fluxos de trabalho reproduzíveis para suas necessidades de computação científica.

Pronto para fazer o curso?

[Começar a aprender :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
