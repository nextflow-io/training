---
title: Configuração de Execução
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Alternar a tecnologia de empacotamento de software entre Docker e Conda
    - Selecionar uma plataforma de execução e entender como o Nextflow adapta a execução de tarefas a ela
    - Controlar alocações de recursos computacionais e reexecutar automaticamente tarefas que falham
    - Definir e combinar perfis para alternar entre configurações predefinidas
  audience_prerequisites:
    - "**Público:** Este curso é destinado a pessoas que já sabem como executar pipelines Nextflow localmente e desejam configurar a execução com mais profundidade."
    - "**Habilidades:** É assumida alguma familiaridade com a linha de comando."
    - "**Cursos:** É necessário ter concluído o [Nextflow Run](../nextflow_run/index.md) ou estar confortável em executar um pipeline local com `nextflow run`."
---

# Configuração de Execução

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


**Configure Execution é uma introdução prática à adaptação da execução de pipelines Nextflow a diferentes ambientes computacionais.**

Trabalhando em exercícios orientados a objetivos, você aprenderá como alternar a tecnologia de empacotamento de software, selecionar uma plataforma de execução, controlar alocações de recursos computacionais e reexecuções, e agrupar configurações em perfis intercambiáveis.

Você sairá com as habilidades e a confiança para configurar a execução de pipelines Nextflow como um profissional.

<!-- additional_information -->

## Visão geral do curso

Este curso é prático e se baseia nas habilidades abordadas no [Nextflow Run](../nextflow_run/index.md).

Você utilizará o mesmo pipeline de múltiplas etapas daquele curso e adaptará progressivamente sua configuração para diferentes ambientes computacionais, agrupando tudo em perfis que podem ser alternados em tempo de execução.

### Plano de aulas

| Capítulo do curso                                                                   | Resumo                                                                                     | Duração estimada |
| ----------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ | ---------------- |
| [Parte 1: Adapte ao seu ambiente computacional](./01_packaging_and_execution.md)    | Alternar a tecnologia de empacotamento de software e selecionar uma plataforma de execução | 20 min           |
| [Parte 2: Gerencie recursos computacionais e falhas](./02_resources_and_retries.md) | Controlar alocações de recursos e reexecutar automaticamente tarefas que falham            | 15 min           |
| [Parte 3: Use perfis para alternar configurações](./03_profiles.md)                 | Definir e combinar perfis, e inspecionar a configuração completamente resolvida            | 15 min           |

Ao final deste curso, você estará confortável para configurar pipelines Nextflow para uma variedade de ambientes computacionais e alternar entre eles com o mínimo de esforço.

Pronto para começar o curso?

[Começar a aprender :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
