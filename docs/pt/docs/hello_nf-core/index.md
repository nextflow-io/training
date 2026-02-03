---
title: Hello nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Recuperar, executar e gerenciar a execução de pipelines nf-core
    - Descrever a estrutura do código e a organização do projeto dos pipelines nf-core
    - Criar um pipeline básico compatível com nf-core a partir de um template
    - Atualizar um fluxo de trabalho Nextflow simples para se adequar aos padrões nf-core
    - Adicionar módulos nf-core a um pipeline compatível com nf-core
    - Contribuir com seus próprios módulos para o nf-core
    - Validar entradas e parâmetros usando ferramentas nf-core
  audience_prerequisites:
    - "**Público:** Este curso é projetado para alunos que já estão familiarizados com Nextflow básico e querem aprender a usar recursos e melhores práticas do nf-core."
    - "**Habilidades:** Familiaridade com a linha de comando, conceitos básicos de script e formatos de arquivo comuns é assumida."
    - "**Cursos:** Deve ter completado o curso [Hello Nextflow](../hello_nextflow/index.md) ou equivalente."
    - "**Domínio:** Os exercícios são todos agnósticos de domínio, portanto nenhum conhecimento científico prévio é necessário."
---

# Hello nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello nf-core é uma introdução prática ao uso de recursos e melhores práticas do nf-core.**

![nf-core logo](./img/nf-core-logo.png#only-light)
![nf-core logo](./img/nf-core-logo-darkbg.png#only-dark)

Através de exemplos práticos e exercícios guiados, você aprenderá a usar e desenvolver módulos e pipelines compatíveis com nf-core, e a utilizar as ferramentas nf-core de forma eficaz.

Você adquirirá as habilidades e a confiança para começar a desenvolver pipelines de acordo com as melhores práticas do nf-core.

<!-- additional_information -->

## Visão geral do curso

Este curso é projetado para ser prático, com exercícios orientados a objetivos estruturados para introduzir informações gradualmente.

Você será apresentado ao [**nf-core**](https://nf-co.re/), um esforço comunitário para desenvolver e manter um conjunto curado de pipelines científicos construídos usando Nextflow, bem como ferramentas e diretrizes relevantes que promovem desenvolvimento aberto, testes e revisão por pares ([Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x), [Genome Biol 26, 228 (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9)).

Os pipelines desenvolvidos pela comunidade nf-core são projetados para serem modulares, escaláveis e portáteis, permitindo que pesquisadores os adaptem e executem facilmente usando seus próprios dados e recursos computacionais.
As diretrizes de melhores práticas aplicadas pelo projeto garantem ainda que os pipelines sejam robustos, bem documentados e validados contra conjuntos de dados do mundo real.
Isso ajuda a aumentar a confiabilidade e reprodutibilidade das análises científicas e, em última análise, permite que pesquisadores acelerem suas descobertas científicas.

Não cobriremos tudo o que há para saber sobre pipelines nf-core neste curso, porque o nf-core abrange muitos recursos e convenções desenvolvidos pela comunidade ao longo dos anos.
Em vez disso, focaremos nos conceitos essenciais que ajudarão você a começar e entender como o nf-core funciona.

### Plano de aulas

Dividimos isso em cinco partes, cada uma focando em aspectos específicos do uso de recursos nf-core.

| Capítulo do curso                                                  | Resumo                                                                                                                                                                         | Duração estimada |
| ------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------- |
| [Parte 1: Execute um pipeline demo](./01_run_demo.md)              | Execute um pipeline nf-core existente e examine sua estrutura de código para ter uma noção do que torna esses pipelines diferentes de fluxos de trabalho Nextflow básicos      | 30 minutos       |
| [Parte 2: Reescreva Hello para nf-core](./02_rewrite_hello.md)     | Adapte um fluxo de trabalho existente ao modelo de estrutura nf-core, começando pelo fluxo de trabalho simples produzido no curso [Hello Nextflow](../hello_nextflow/index.md) | 60 minutos       |
| [Parte 3: Use um módulo nf-core](./03_use_module.md)               | Explore a biblioteca de módulos da comunidade e aprenda a integrar módulos pré-construídos e testados que encapsulam ferramentas de bioinformática comuns                      | 30 minutos       |
| [Parte 4: Crie um módulo nf-core](./04_make_module.md)             | Crie seu próprio módulo no estilo nf-core usando a estrutura específica, convenções de nomenclatura e requisitos de metadados estabelecidos pelo nf-core                       | 30 minutos       |
| [Parte 5: Adicione validação de entrada](./05_input_validation.md) | Implemente validação de entrada tanto para parâmetros de linha de comando quanto para arquivos de dados de entrada usando nf-schema                                            | 30 minutos       |

Ao final deste curso, você será capaz de aproveitar a enorme riqueza de recursos oferecidos pelo projeto nf-core.

Pronto para fazer o curso?

[Começar a aprender :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
