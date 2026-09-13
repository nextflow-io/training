---
title: Use nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Encontrar, recuperar e executar pipelines da comunidade nf-core
    - Configurar a execução de pipelines usando parâmetros e arquivos de configuração
    - Entender como os pipelines nf-core validam parâmetros e dados de entrada
    - Executar um pipeline em escala de produção (nf-core/rnaseq) e substituir suas alocações de recursos padrão
  audience_prerequisites:
    - "**Público:** Este curso é destinado a pessoas que já sabem executar pipelines Nextflow locais, são novas no nf-core e querem executar pipelines da comunidade existentes."
    - "**Habilidades:** É esperada alguma familiaridade com a linha de comando, conceitos básicos de scripting e formatos de arquivo comuns."
    - "**Cursos:** É necessário ter concluído o [Nextflow Run](../nextflow_run/index.md) ou estar confortável em executar um pipeline local com `nextflow run`."
    - "**Área:** Os exercícios utilizam pipelines de bioinformática, mas não é necessário conhecimento científico prévio na área."
---

# Use nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Use nf-core é uma introdução prática para encontrar, executar e configurar pipelines da comunidade nf-core.**

Trabalhando com exemplos práticos e exercícios guiados, você aprenderá a encontrar e recuperar pipelines nf-core, executá-los usando seus perfis de teste integrados e personalizar sua execução por meio de parâmetros e arquivos de configuração.

Você sairá com as habilidades e a confiança necessárias para começar a executar pipelines nf-core em suas próprias análises.

<!-- additional_information -->

## Visão geral do curso

Este curso é prático, com exercícios orientados a objetivos estruturados para introduzir as informações gradualmente.

Você começará com o `nf-core/demo`, um pipeline mínimo mantido pelo projeto nf-core para fins de treinamento, e depois aplicará o que aprendeu ao `nf-core/rnaseq`, um pipeline de produção amplamente utilizado para análise de sequenciamento de RNA em bulk.

Este curso foca na execução de pipelines.
Se você está procurando uma introdução ao desenvolvimento de pipelines compatíveis com nf-core, consulte [Build with nf-core](../nfcore_build/index.md).

### Plano de aulas

| Capítulo do curso                                                    | Resumo                                                                                                          | Duração estimada |
| -------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | ---------------- |
| [Parte 1: Executar um pipeline de demonstração](./01_run_demo.md)    | Encontrar e recuperar um pipeline nf-core e executá-lo usando seu perfil de teste                               | 20 min           |
| [Parte 2: Configurar a execução do pipeline](./02_configure_execution.md) | Definir parâmetros, entender a validação e personalizar a alocação de recursos e os argumentos das ferramentas | 20 min           |
| [Parte 3: Executar um pipeline de produção](./03_run_production_pipeline.md) | Baixar e executar o nf-core/rnaseq e substituir suas alocações de recursos padrão                         | 20 min           |

Ao final deste curso, você será capaz de aproveitar a riqueza de pipelines da comunidade oferecidos pelo projeto nf-core.

Pronto para começar o curso?

[Começar a aprender :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
