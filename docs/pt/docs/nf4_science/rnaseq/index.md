---
title: Nextflow para RNAseq
hide:
  - toc
---

# Nextflow para RNAseq

Este curso de treinamento é destinado a pesquisadores em transcriptômica e áreas relacionadas que estão interessados em desenvolver ou personalizar pipelines de análise de dados.
Ele se baseia no treinamento para iniciantes [Hello Nextflow](../../hello_nextflow/) e demonstra como usar o Nextflow no contexto específico de análise de RNAseq em massa (bulk).

Especificamente, este curso demonstra como implementar um pipeline simples de processamento de RNAseq em massa para aparar sequências adaptadoras, alinhar as leituras a um genoma de referência e realizar controle de qualidade (QC) em várias etapas.

Vamos começar! Clique no botão "Open in GitHub Codespaces" abaixo para iniciar o ambiente de treinamento (de preferência em uma aba separada), e continue lendo enquanto ele carrega.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Objetivos de aprendizagem

Ao trabalhar neste curso, você aprenderá como aplicar conceitos e ferramentas fundamentais do Nextflow a um caso de uso típico de RNAseq.

Ao final deste workshop, você será capaz de:

- Escrever um fluxo de trabalho linear para aplicar métodos básicos de processamento e QC de RNAseq
- Manipular arquivos específicos do domínio, como FASTQ e recursos de genoma de referência, de forma apropriada
- Manipular dados de sequenciamento single-end e paired-end
- Aproveitar o paradigma de fluxo de dados do Nextflow para paralelizar o processamento de RNAseq por amostra
- Agregar relatórios de QC em múltiplas etapas e amostras usando operadores de canal relevantes

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Pré-requisitos

O curso assume alguma familiaridade mínima com o seguinte:

- Ferramentas e formatos de arquivo comumente usados neste domínio científico
- Experiência com a linha de comando
- Conceitos e ferramentas fundamentais do Nextflow abordados no treinamento para iniciantes [Hello Nextflow](../../hello_nextflow/).

Para requisitos técnicos e configuração do ambiente, consulte o mini-curso [Configuração do Ambiente](../../envsetup/).
