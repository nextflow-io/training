---
title: Nextflow para RNAseq
hide:
  - toc
---

# Nextflow para RNAseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Este curso de treinamento é destinado a pesquisadores em transcriptômica e áreas relacionadas que estão interessados em desenvolver ou personalizar pipelines de análise de dados.
Ele se baseia no treinamento para iniciantes [Hello Nextflow](../../hello_nextflow/) e demonstra como usar Nextflow no contexto específico de análise de RNAseq bulk.

Especificamente, este curso demonstra como implementar um pipeline simples de processamento de RNAseq bulk para remover sequências adaptadoras, alinhar as reads a um genoma de referência e realizar controle de qualidade (QC) em várias etapas.

Vamos começar! Clique no botão "Open in GitHub Codespaces" abaixo para iniciar o ambiente de treinamento (preferencialmente em uma aba separada), depois continue lendo enquanto ele carrega.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Objetivos de aprendizagem

Ao trabalhar neste curso, você aprenderá como aplicar conceitos e ferramentas fundamentais do Nextflow a um caso de uso típico de RNAseq.

Ao final deste workshop, você será capaz de:

- Escrever um fluxo de trabalho linear para aplicar métodos básicos de processamento e QC de RNAseq
- Manipular arquivos específicos do domínio, como FASTQ e recursos de genoma de referência, de forma apropriada
- Manipular dados de sequenciamento single-end e paired-end
- Aproveitar o paradigma de dataflow do Nextflow para paralelizar o processamento de RNAseq por amostra
- Agregar relatórios de QC através de múltiplas etapas e amostras usando operadores de canal relevantes

<!-- TODO
- Configurar a execução do pipeline e gerenciar e otimizar alocações de recursos
- Implementar testes por etapa e de ponta a ponta do pipeline que lidam apropriadamente com as idiossincrasias específicas de RNAseq
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Pré-requisitos

O curso assume alguma familiaridade mínima com o seguinte:

- Ferramentas e formatos de arquivo comumente usados neste domínio científico
- Experiência com a linha de comando
- Conceitos e ferramentas fundamentais do Nextflow abordados no treinamento para iniciantes [Hello Nextflow](../../hello_nextflow/)

Para requisitos técnicos e configuração do ambiente, consulte o mini-curso [Configuração do Ambiente](../../envsetup/).
