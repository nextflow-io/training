---
title: Nextflow para Genômica
hide:
  - toc
---

# Nextflow para Genômica

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Este curso de treinamento é destinado a pesquisadores em genômica e áreas relacionadas que estão interessados em desenvolver ou personalizar pipelines de análise de dados.
Ele se baseia no treinamento para iniciantes [Hello Nextflow](../../hello_nextflow/) e demonstra como usar o Nextflow no contexto específico do domínio de genômica.

Especificamente, este curso demonstra como implementar um pipeline simples de chamada de variantes com [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), um pacote de software amplamente utilizado para analisar dados de sequenciamento de alto rendimento.

Vamos começar! Clique no botão "Open in GitHub Codespaces" abaixo para iniciar o ambiente de treinamento (preferencialmente em uma aba separada) e continue lendo enquanto ele carrega.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Objetivos de aprendizagem

Ao trabalhar neste curso, você aprenderá como aplicar conceitos e ferramentas fundamentais do Nextflow a um caso de uso típico de genômica.

Ao final deste workshop, você será capaz de:

- Escrever um fluxo de trabalho linear para aplicar chamada de variantes a uma única amostra
- Manipular arquivos acessórios, como arquivos de índice e recursos de genoma de referência, adequadamente
- Aproveitar o paradigma de dataflow do Nextflow para paralelizar a chamada de variantes por amostra
- Implementar chamada de variantes multi-amostra usando operadores de canal relevantes
- Implementar testes de pipeline por etapa e de ponta a ponta que lidam com idiossincrasias específicas de genômica adequadamente

<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Pré-requisitos

O curso assume alguma familiaridade mínima com o seguinte:

- Ferramentas e formatos de arquivo comumente usados neste domínio científico
- Experiência com a linha de comando
- Conceitos e ferramentas fundamentais do Nextflow cobertos no treinamento para iniciantes [Hello Nextflow](../../hello_nextflow/)

Para requisitos técnicos e configuração de ambiente, consulte o mini-curso [Configuração de Ambiente](../../envsetup/).
